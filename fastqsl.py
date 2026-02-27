import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os, subprocess, pickle

def fastqsl(Bx, By=None, Bz=None, *, xa=None, ya=None, za=None, spherical=False, \
            xperiod=False, yperiod=False, zperiod=False, \
            xreg=None, yreg=None, zreg=None, csFlag=False, factor=4, delta=None, \
            lon_delta=None, lat_delta=None, r_delta=None, arc_delta=None, seed=None, \
            RK4Flag=False, step=1.0, tol=1.0e-4, maxsteps=None, \
            scottFlag=False, inclineFlag=False, r_local=0., nthreads=0, silent=False, \
            B_out=False, CurlB_out=False, length_out=False, twist_out=False, \
            rF_out=False, targetB_out=False, targetCurlB_out=False, \
            path_out=False, loopB_out=False, loopCurlB_out=False, \
            odir=None, fname=None, save_file=False, preview=False, tmp_dir=None, keep_tmp=False):
# ------------------------------------------------------------
    # check input
    B3flag= By is not None and Bz is not None
    
    Bx_ndim=np.array(Bx).ndim
    sBx    =np.array(Bx).shape

    if B3flag:
        if Bx_ndim != 3 or np.array(By).ndim != 3 or np.array(Bz).ndim != 3: 
            raise Exception('Bx, By and Bz must be 3D arrays!')
        if sBx != np.array(By).shape or sBx != np.array(Bz).shape:
            raise Exception('Bx, By and Bz must have the same dimensions!')
    else: 
        # the dimensions of Bx are (nz,ny,nx,3)
        if Bx_ndim != 4 or sBx[3] != 3: raise Exception('Something is wrong with the magnetic field')

    nz, ny, nx = sBx[0:3]

    if nx < 2 or ny < 2 or nz < 2: raise Exception('The thickness of Bx should not be smaller than 2')
    elif not spherical:
        if nx == 2: xperiod=True
        if ny == 2: yperiod=True
        if nz == 2: zperiod=True

    stretchFlag = (xa is not None) and (ya is not None) and (za is not None)
    if stretchFlag:
        if (nx != len(xa) or ny != len(ya) or nz != len(za)):
            raise Exception('The size of xa, ya and za must be consistant with the dimensions of the magnetic field')
    elif spherical: raise Exception('xa, ya and za should be specified in spherical coordinates')
# ------------------------------------------------------------
    #  understand the output grid

    # provide qsl.seed even sflag is False
    launch_out = seed is True

    sflag = seed is not None and not launch_out

    if sflag:
        if isinstance(seed, str):
            if seed == 'original':
                seed = np.zeros((nz, ny, nx, 3), dtype='f4')
                if stretchFlag:
                    for i in range(nx): seed[:, :, i, 0] = xa[i]
                    for i in range(ny): seed[:, i, :, 1] = ya[i]
                    for i in range(nz): seed[i, :, :, 2] = za[i]
                else:
                    for i in range(nx): seed[:, :, i, 0] = i
                    for i in range(ny): seed[:, i, :, 1] = i
                    for i in range(nz): seed[i, :, :, 2] = i
            elif seed == 'original_bottom':
                seed = np.zeros((ny, nx, 3), dtype='f4')
                if stretchFlag:
                    for i in range(nx): seed[:, i, 0] = xa[i]
                    for i in range(ny): seed[i, :, 1] = ya[i]
                    seed[:, :, 2] = za[0]
                else:
                    for i in range(nx): seed[:, i, 0] = i
                    for i in range(ny): seed[i, :, 1] = i
            else: raise ValueError('Something is wrong with seed')
        else: seed=np.array(seed, dtype='f4')

        out_dim= seed.ndim-1
        if seed.shape[out_dim] != 3: raise ValueError('Something is wrong with seed')
        nq1=seed.shape[out_dim-1] if out_dim >= 1 else 1
        nq2=seed.shape[out_dim-2] if out_dim >= 2 else 1
        nq3=seed.shape[out_dim-3] if out_dim == 3 else 1

        bflag, cflag, vflag = False, False, False
    else :
        # if spherical and xa(nx-1) is very close to 2*!pi but not equivalent, 
        # xa(nx-1) will be reset to 2*!pi.
        # then if ~preset_xreg, xreg[1] will be reset to 2*!pi
        preset_xreg= xreg is not None

        # if (.not. preset_yreg .and. (pole_south .and. pole_north) in fastqsl.x, 
        # yreg will be reset to [-!pi/2., !pi/2.]
        preset_yreg= yreg is not None

        if stretchFlag:
            if xreg is None: xreg = [xa[0], xa[-1]]
            if yreg is None: yreg = [ya[0], ya[-1]]
            if zreg is None: zreg = [za[0], za[0]]
        else :
            if xreg is None: xreg = [0, nx-1]
            if yreg is None: yreg = [0, ny-1]
            if zreg is None: zreg = [0, 0]
        
        if len(xreg) != 2 or len(yreg) != 2 or len(zreg) != 2 : 
            raise Exception('xreg, yreg and zreg must be 2-element arrays!')
        if sum([xreg[0] == xreg[1], yreg[0] == yreg[1], zreg[0] == zreg[1]]) >= 2 \
        or (csFlag and (zreg[0] == zreg[1])):
            raise Exception('Something is wrong with the cross section')
        if (spherical and csFlag):
            p0= [np.cos(yreg[0])*np.cos(xreg[0]), np.cos(yreg[0])*np.sin(xreg[0]), np.sin(yreg[0])]
            p1= [np.cos(yreg[1])*np.cos(xreg[1]), np.cos(yreg[1])*np.sin(xreg[1]), np.sin(yreg[1])]
            # if p0 \cdot p1 eq 1 or -1, there are many great circles can pass p0 and p1
            if np.dot(p0,p1) > 0.999: raise Exception('The great circle can not be clearly defined')
        zmin = za[0] if stretchFlag else 0.0
        bflag = zreg[0] == zmin and zreg[1] == zmin
        vflag = (xreg[1] != xreg[0]) and (yreg[1] != yreg[0]) and (zreg[1] != zreg[0]) and (not csFlag)
        cflag = not (vflag or bflag)
        out_dim = 3 if vflag else 2

	    # grid spacing for output
        if stretchFlag:
            if spherical:
                if delta is None:
                    if lon_delta is None: lon_delta = (xa[nx-1]-xa[0])/((nx-1.0)*factor)
                    if lat_delta is None: lat_delta = (ya[ny-1]-ya[0])/((ny-1.0)*factor)
                    if   r_delta is None:   r_delta = (za[nz-1]-za[0])/((nz-1.0)*factor)
                else:
                    if lon_delta is None: lon_delta = delta/za[0]
                    if lat_delta is None: lat_delta = delta/za[0]
                    if   r_delta is None:   r_delta = delta
                if     arc_delta is None: arc_delta = min([lon_delta, lat_delta])
            elif           delta is None:     delta = (xa[nx-1]-xa[0])/((nx-1.0)*factor)
        elif               delta is None:     delta = 1.0/factor
    # end if sflag
# ------------------------------------------------------------
    if maxsteps is None: maxsteps=int(4*(nx+ny+nz)/step) if RK4Flag else 4*(nx+ny+nz)

    scottFlag       = scottFlag and (maxsteps != 0)
    twist_out       = twist_out and (maxsteps != 0)
    length_out      = length_out and (maxsteps != 0)

    rF_out          = rF_out and (maxsteps != 0)
    targetB_out     = targetB_out and (maxsteps != 0)
    targetCurlB_out = targetCurlB_out and (maxsteps != 0)

    if (out_dim == 3 and path_out) :
        print('invoking path_out for a 3D ouput grid is not allowed, path_out is changed to False')
    
    path_out        = path_out and (maxsteps != 0) and out_dim < 3
    loopB_out       = loopB_out and path_out
    loopCurlB_out   = loopCurlB_out and path_out

    verbose         = not silent
# ------------------------------------------------------------
    # the directory for output
    # https://learn.microsoft.com/en-us/dotnet/standard/io/file-path-formats#canonicalize-separators
    # Actually Windows also accept '/', set as below just for a better readability in Windows
    cdir = os.getcwd()+os.sep 
    if preview or save_file:
        if odir is not None:
            if odir[-1] != os.sep: odir=odir+os.sep
        else: odir= cdir+'fastqsl'+os.sep
        os.makedirs(odir, exist_ok=True)
# ------------------------------------------------------------   
    # the temporary directory for the data transmission between fastqsl.x and fastqsl.py 
    if tmp_dir is not None: 
        if tmp_dir[-1] != os.sep : tmp_dir=tmp_dir+os.sep
    else: tmp_dir= cdir+'tmpFastQSL'+os.sep
    old_tmp_dir=os.path.exists(tmp_dir)
    if old_tmp_dir:
        for file in os.listdir(tmp_dir): 
            if file[-4:] == '.bin': os.remove(tmp_dir+file)
    else: os.makedirs(tmp_dir, exist_ok=True)
# ------------------------------------------------------------
    # transmit data to fastqsl.x
    with open(tmp_dir+'head.bin','wb') as file:
        file.write(np.array([step, tol, r_local], dtype='f4'))
        file.write(np.array([maxsteps, RK4Flag, inclineFlag, \
        launch_out, B_out, CurlB_out, length_out, twist_out, \
		rF_out, targetB_out, targetCurlB_out, \
		path_out, loopB_out, loopCurlB_out, \
		sflag, bflag, cflag, vflag, nthreads, scottFlag, verbose], dtype='i4'))

    with open(tmp_dir+'bfield.bin','wb') as file:
        file.write(np.array([nx, ny, nz, B3flag, spherical, preview, xperiod, yperiod, zperiod], dtype='i4'))
        file.write(np.array(Bx, dtype='f4', order='C'))
        if B3flag: file.write(np.array([By,Bz], dtype='f4', order='C'))

    if stretchFlag:
        with open(tmp_dir+'axis.bin','wb') as file:
            file.write(np.array(xa, dtype='f4'))
            file.write(np.array(ya, dtype='f4'))
            file.write(np.array(za, dtype='f4'))
    if sflag:
        with open(tmp_dir+'dim_seed.bin','wb') as file: file.write(np.array([nq1,nq2,nq3], dtype='i4'))
        with open(tmp_dir+'seed.bin','wb') as file: file.write(np.array(seed, dtype='f4', order='C'))
    else:
        deltas=[arc_delta, lon_delta, lat_delta, r_delta] if spherical else [delta,delta,delta,delta]
        with open(tmp_dir+'head_region.bin','wb') as file:
            file.write(np.array([xreg,yreg,zreg,deltas[0:2],deltas[2:4]], dtype='f4'))
            file.write(np.array([csFlag, preset_xreg, preset_yreg], dtype='i4'))
# ------------------------------------------------------------
    # computed by fastqsl.x
    os.chdir(tmp_dir)

    # please specify the path
    # the following r can avoid the potential problem of '\n' from os.sep ='\' in Windows
    subprocess.run(r'/path/of/fastqsl.x', shell=True)
    os.chdir(cdir)
# ################################### retrieving results ######################################
# make the dictionary qsl
    qsl={}

    if maxsteps != 0:
        if RK4Flag: qsl.update({'step':np.array(step, dtype="f4")})
        else: qsl.update({'tol':np.array(tol, dtype="f4")})

    # the output grid   
    if not sflag:
        with open(tmp_dir+'tail_region.bin','rb') as file:
            nq1, nq2, nq3, normal_index = np.fromfile(file, dtype="i4", count=4) 
            dummy = np.fromfile(file, dtype="f4")
            xreg=dummy[0:2]
            yreg=dummy[2:4]
            zreg=dummy[4:6]

        if spherical:
            qsl.update({'lon_reg':xreg, 'lat_reg':yreg, 'r_reg':zreg})
            if normal_index == -1: qsl.update({'arc_delta':np.array(arc_delta,'f4')})
            if normal_index in [0, 2]: qsl.update({'lat_delta':np.array(lat_delta,'f4')})
            if normal_index in [1, 2]: qsl.update({'lon_delta':np.array(lon_delta,'f4')})
            if normal_index != 2 or vflag: qsl.update({'r_delta':np.array(r_delta,'f4')})
        else: qsl.update({'xreg':xreg, 'yreg':yreg, 'zreg':zreg, 'delta':np.array(delta,'f4')})
    # end if not sflag

    if os.path.exists(tmp_dir+'q_local.bin'): qsl.update({'r_local':np.array(r_local, dtype="f4")})

    if   out_dim== 0:
        dim =(1)
        dim3=(3)
    elif out_dim== 1: 
        dim =(nq1)
        dim3=(nq1,3)
    elif out_dim== 2: 
        dim =(nq2,nq1)
        dim3=(nq2,nq1,3)
    elif out_dim== 3: 
        dim =(nq3,nq2,nq1)
        dim3=(nq3,nq2,nq1,3)

    qsl_data0=[ \
    ['axis1',      'f4', (nq1, 2)], \
    ['seed',       'f4', dim3], \
    ['q',          'f4', dim], \
    ['q_perp',     'f4', dim], \
    ['q_local',    'f4', dim], \
    ['rboundary',  'i1', dim], \
    ['sign2d',     'i2', (nq2, nq1)], \
    ['length',     'f4', dim], \
    ['twist',      'f4', dim], \
    ['B',          'f4', dim3], \
    ['CurlB',      'f4', dim3], \
    ['rFs',        'f4', dim3], \
    ['rFe',        'f4', dim3], \
    ['Bs',         'f4', dim3], \
    ['Be',         'f4', dim3], \
    ['CurlBs',     'f4', dim3], \
    ['CurlBe',     'f4', dim3], \
    ['path',       'f4', (-1,3)], \
    ['loopB',      'f4', (-1,3)], \
    ['loopCurlB',  'f4', (-1,3)], \
    ['index_seed', 'i4', dim]]

    qsl_data=[str3 for str3 in qsl_data0 if os.path.exists(tmp_dir+str3[0]+'.bin')]    
    if path_out:
        with open(tmp_dir+'indexes.bin','rb') as file: indexes=np.fromfile(file, dtype='i8')

    for str3 in qsl_data:
        name=str3[0]
        with open(tmp_dir+name+'.bin','rb') as file: 
            dummy=np.fromfile(file, dtype=str3[1]).reshape(str3[2])
        if name in ['path', 'loopB','loopCurlB']:
            if   out_dim <=1: qsl.update({name:[dummy[indexes[i]:indexes[i+1],:] \
                                      for i in range(nq1)]})
            elif out_dim ==2: qsl.update({name:[[dummy[indexes[i+j*nq1]:indexes[i+j*nq1+1],:] \
                                      for i in range(nq1)] for j in range(nq2)]})
            elif out_dim ==3: qsl.update({name:[[[dummy[indexes[i+j*nq1+k*nq1*nq2]:indexes[i+j*nq1+k*nq1*nq2+1],:] \
                                      for i in range(nq1)] for j in range(nq2)] for k in range(nq3)]})
        else: qsl.update({name:dummy})
# ------------------------------------------------------------
# the name of .pkl file
    if fname is None and (preview or save_file):

        if sflag: head_str='seed'
        if bflag: head_str='bottom'
        if cflag: head_str='cs'
        if vflag: head_str='volume'

        if not (spherical or sflag):
            decimal3=round(1000.*(delta-np.floor(delta)))
            for i in range(1,5): 
                if (decimal3 % 10**i) != 0: break
            if i >= 4: delta_str='_delta'+f'{round(delta):d}'
            else: delta_str='_delta' + ('{:0.'+f'{4-i:d}'+'f}').format(float(delta))
        else: delta_str=''

        if cflag and not csFlag:
            if   normal_index== 0:
                cut_coordinate=xreg[0]
                cut_str0='_lon' if spherical else '_x'
            elif normal_index== 1:
                cut_coordinate=yreg[0]
                cut_str0='_lat' if spherical else '_y'
            elif normal_index== 2:
                cut_coordinate=zreg[0]
                cut_str0='_r'   if spherical else '_z'

            decimal3=round(1000.*(cut_coordinate-np.floor(cut_coordinate)))
            for i in range(1,5): 
                if (decimal3 % 10**i) != 0: break
            if i >= 4: cut_str=cut_str0+ f'{round(cut_coordinate):d}'
            else:      cut_str=cut_str0+ ('{:0.'+f'{4-i:d}'+'f}').format(float(cut_coordinate))
        else: cut_str=''

        fname = head_str + delta_str + cut_str
        
    if save_file:
        with open(odir+fname+'.pkl', 'wb') as file: pickle.dump(qsl, file)
# ------------------------------------------------------------
    if preview:
        if verbose: print('\n'+'These images are produced:')

        if stretchFlag:
            with open(tmp_dir+'magnetogram.bin','rb') as file:
                nx_mag, ny_mag = np.fromfile(file, dtype='i4', count=2)
                mag_delta  = np.fromfile(file, dtype='f4', count=1)[0]
                magnetogram= np.fromfile(file, dtype='f4').reshape(ny_mag, nx_mag)
                
            extent=[xa[0], xa[0]+mag_delta*nx_mag, ya[0], ya[0]+mag_delta*ny_mag]
        else: 
            magnetogram=Bz[0,:,:].copy() if B3flag else Bx[0,:,:,2].copy()
            nx_mag, ny_mag= nx, ny
            extent=[0, nx_mag-1, 0, ny_mag-1]  

        # mark the area for calculation or plot seed and their field lines on the magnetogram
        plt.figure(figsize=(nx_mag*0.02, ny_mag*0.02), dpi=100)
        if spherical: two_pi=2*np.pi
        if sflag:
            if out_dim <=1:

                if path_out:
                    for i in range(nq1): 
                        x_path=qsl['path'][i][:,0].copy()
                        if spherical: x_path=np.mod(x_path-xa[0], two_pi)+xa[0]
                        plt.plot(x_path, qsl['path'][i][:,1], '.g', markersize=1)

                if out_dim ==0:
                     x_seed=np.array(seed[0],'f4')
                     y_seed=np.array(seed[1],'f4')
                else: 
                    x_seed=np.array(seed[:,0],'f4')
                    y_seed=np.array(seed[:,1],'f4')

                if spherical: x_seed=np.mod(x_seed-xa[0], two_pi)+xa[0]
                plt.plot(x_seed,y_seed, 'xr', markersize=5)
            elif out_dim ==2:
                x_margin=np.concatenate((seed[:,0,0],seed[:,-1,0],seed[0,:,0],seed[-1,:,0]))
                y_margin=np.concatenate((seed[:,0,1],seed[:,-1,1],seed[0,:,1],seed[-1,:,1]))
            elif out_dim ==3:
                x_margin=np.concatenate((seed[0,:,0,0],seed[0,:,-1,0],seed[0,0,:,0],seed[0,-1,:,0]))
                y_margin=np.concatenate((seed[0,:,0,1],seed[0,:,-1,1],seed[0,0,:,1],seed[0,-1,:,1]))
        else:
            if csFlag:
                if spherical:
                    x_margin=qsl['axis1'][:,0]
                    y_margin=qsl['axis1'][:,1]
                else:
                    x_margin=xreg
                    y_margin=yreg
            else:
                x_margin=[xreg[0],xreg[1],xreg[1],xreg[0],xreg[0]]
                y_margin=[yreg[0],yreg[0],yreg[1],yreg[1],yreg[0]]
        # end if sflag

        if out_dim >=2: 
            style='-r' if len(x_margin) <=5 else '.r'

            plt.plot(x_margin,y_margin, style, markersize=1)

            # some regions may stand across at longitude = 0 or 2*!pi
            if spherical: 
                plt.plot(np.array(x_margin)+two_pi, y_margin, style, markersize=1)
                plt.plot(np.array(x_margin)-two_pi, y_margin, style, markersize=1)

        magnetogram[np.isnan(magnetogram)]=0.
        magnetogram[np.isinf(magnetogram)]=0.  
        scale_top=  np.max(np.abs(magnetogram))/4.0 
        if scale_top > 1000.: scale_top=1000.

        # scale_top=100.
        plt.imshow(magnetogram,extent=extent,origin='lower',vmin=-scale_top,vmax=scale_top,cmap='gray')
        # plt.colorbar()
        plt.axis('off')
        plt.savefig(odir+fname+'_magnetogram.png', bbox_inches='tight', pad_inches=0)
        plt.close()
        if verbose: print(odir+fname+'_magnetogram.png')
        # ------------------------------------------------------------
        # preview q/q_perp, length, twist

        # if vflag and the bottom plane is included, 'sign2d.bin' also can be found
        plot_bottom=os.path.exists(tmp_dir+'sign2d.bin')

        if maxsteps != 0 and (out_dim ==2 or plot_bottom):
            rb_tmp=qsl['rboundary'] if out_dim ==2 else qsl['rboundary'][0,:,:]
            rbs = rb_tmp // 10
            rbe = np.mod(rb_tmp,10)
            cmap_rboundary = ListedColormap([(0,0,0), (0.5,0.5,0.5), (1,1,1), \
            (0.75,0,0), (1,0.5,0), (0,0.5,0), (0,0.5,1), (0,0,0.5), (0.5,0,0.5)], N=9)

            if plot_bottom:
                rb_target=np.zeros((nq2,nq1),'i1')
                rb_target[qsl['sign2d'] == 1]=rbe[qsl['sign2d'] == 1]
                rb_target[qsl['sign2d'] ==-1]=rbs[qsl['sign2d'] ==-1]

                # for the case of qsl['sign2d'] ==0
                rb_target[rbs == rbe] = rbs[rbs == rbe]

                plt.imsave(odir+fname+'_rb_target.png', rb_target, \
                           vmin=-0.5, vmax=8.5, origin='lower', cmap=cmap_rboundary)
                if verbose: print(odir+fname+'_rb_target.png')
            else:
                plt.imsave(odir+fname+'_rbs.png', rbs, \
                           vmin=-0.5,vmax=8.5, origin='lower', cmap=cmap_rboundary)
                plt.imsave(odir+fname+'_rbe.png', rbe, \
                           vmin=-0.5,vmax=8.5, origin='lower', cmap=cmap_rboundary)
                if verbose: 
                    print(odir+fname+'_rbs.png')
                    print(odir+fname+'_rbe.png')
                
                closed=np.full((nq2,nq1),1.)
                closed[rb_tmp == 11]=0.5
                plt.imsave(odir+fname+'_mark_closed.png', closed, \
                           vmin=0.,vmax=1., origin='lower', cmap='gray')
                print(odir+fname+'_mark_closed.png')
            # end if plot_bottom

            # q/q_perp/q_local
            q_strs=[str3[0] for str3 in qsl_data if str3[0] in ['q', 'q_perp','q_local' ]]
            
            for q_str in q_strs:
                q_tmp = qsl[q_str].copy() if out_dim ==2 else qsl[q_str][0,:,:].copy()

                # 1. for white color
                q_tmp[np.isnan(q_tmp)]=1.
                q_tmp[np.isinf(q_tmp)]=1.

                if plot_bottom:
                    plt.imsave(odir+fname+'_slog'+q_str+'.png', \
                    np.log10(q_tmp)*qsl['sign2d'], vmin=-5., vmax=5., origin='lower', cmap='bwr')
                    if verbose: print(odir+fname+'_slog'+q_str+'.png')

                    if targetB_out and q_str=='q': q_tmp1=q_tmp.copy()

                    # make white color for open field line
                    # q_tmp[rb_tmp != 11]=1. 
                    # plt.imsave(odir+fname+'_slog'+q_str+'.png', \
                    # np.log10(q_tmp)*qsl['sign2d'], vmin=-5., vmax=5., origin='lower', cmap='bwr')
                    # if verbose: print(odir+fname+'_slog'+q_str+'_orig.png')              
                else:
                    plt.imsave(odir+fname+ '_log'+q_str+'.png', \
                    np.log10(q_tmp), vmin=1., vmax=5., origin='lower', cmap='gray')
                    if verbose: print(odir+fname+'_log'+q_str+'.png')

                    # make white color for open field line
                    # q_tmp[rb_tmp != 11]=1. 
                    # plt.imsave(odir+fname+ '_log'+q_str+'.png', \
                    # np.log10(q_tmp), vmin=1., vmax=5., origin='lower', cmap='gray')
                    # if verbose: print(odir+fname+'_log'+q_str+'_orig.png')
            # end q_str in q_strs:

            if targetB_out and plot_bottom and (q_tmp1 is not None):

                # In Titov (2007), q = N^2 / Delta, and Delta = Bnr
                Bnr=np.zeros((nq2,nq1),'f4')

                Bs=qsl['Bs'] if out_dim ==2 else qsl['Bs'][0,:,:,:]
                Be=qsl['Be'] if out_dim ==2 else qsl['Be'][0,:,:,:]
                
                for j in range(nq2):
                    for i in range(nq1): 
                        if   qsl['sign2d'][j,i] ==  1 and rbe[j,i] in range(1,7):
                            Bn_target=Be[j,i,(6-rbe[j,i])//2]
                            if Bn_target != 0.: Bnr[j,i]=np.abs(Bs[j,i,2]/Bn_target)
                        elif qsl['sign2d'][j,i] == -1 and rbs[j,i] in range(1,7):
                            Bn_target=Bs[j,i,(6-rbs[j,i])//2]
                            if Bn_target != 0.: Bnr[j,i]=np.abs(Be[j,i,2]/Bn_target)
                    
                N=np.sqrt(q_tmp1 * Bnr)  # Priest and Demoulin (1995)
            
                # 1. for black color
                N[N == 0.    ]=1.
                N[np.isinf(N)]=1.
                N[np.isnan(N)]=1.
                plt.imsave(odir+fname+'_N.png', np.log10(N), \
                           vmin=0.5, vmax=3., origin='lower', cmap='gray')
                if verbose: print(odir+fname+'_N.png')
                    
                # 1. for white color, should be done after N is calculated
                Bnr[Bnr == 0.    ]=1.
                Bnr[np.isinf(Bnr)]=1.
                Bnr[np.isnan(Bnr)]=1.

                plt.imsave(odir+fname+'_lg_Bnr.png', np.log10(Bnr), \
                           vmin=-2., vmax=2., origin='lower', cmap='bwr')
                if verbose: print(odir+fname+'_lg_Bnr.png')
            # end if targetB_out

            if length_out and out_dim == 2:
                length_top=2.*(nz-1) if not stretchFlag else 2.*(za[nz-1]-za[0])
                length_tmp=qsl['length'].copy()
                length_tmp[np.isnan(length_tmp)]=0. # 0. for black color
                length_tmp[np.isinf(length_tmp)]=0.
                plt.imsave(odir+fname+'_length.png', length_tmp, \
                    vmin=0., vmax=length_top, origin='lower', cmap='gray')
                print(odir+fname+'_length.png')

            if twist_out and out_dim == 2:
                twist_tmp=qsl['twist'].copy() #*(-1.)
                twist_tmp[np.isnan(twist_tmp)]=0. # 0. for white color
                twist_tmp[np.isinf(twist_tmp)]=0.
                plt.imsave(odir+fname+'_twist.png', twist_tmp, \
                    vmin=-2., vmax=2, origin='lower', cmap='bwr')
                print(odir+fname+'_twist.png')

        # end if maxsteps != 0 and (out_dim ==2 or plot_bottom)
    # end if preview
# ------------------------------------------------------------
    if not keep_tmp:
        for str3 in qsl_data: os.remove(tmp_dir+str3[0]+'.bin')
        if path_out: os.remove(tmp_dir+'indexes.bin')
        if not sflag: os.remove(tmp_dir+'tail_region.bin')
        if preview and stretchFlag: os.remove(tmp_dir+'magnetogram.bin')
        if not old_tmp_dir: os.rmdir(tmp_dir)

    if verbose:
        print('\n'+'Elements in qsl:')
        for key in qsl.keys():
            if   isinstance(qsl[key], np.ndarray):
                content = qsl[key] if qsl[key].size <=3 else qsl[key].shape
                print('{0:<20}{1:<10}'.format("qsl['"+key+"']", qsl[key].dtype.name), content)
            elif isinstance(qsl[key], list): print('{0:<19}'.format("qsl['"+key+"']"), "list")
        if save_file:
            print('Try:')
            print('with open("'+odir+fname+'.pkl"'+', "rb") as file: qsl = pickle.load(file)')
    return qsl
