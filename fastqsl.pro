PRO fastqsl, Bx, By, Bz, xa=xa, ya=ya, za=za, spherical=spherical,                $
	        xperiod=xperiod, yperiod=yperiod, zperiod=zperiod,                    $
            xreg=xreg, yreg=yreg, zreg=zreg, csFlag=csFlag,                       $
            factor=factor, delta=delta, lon_delta=lon_delta, lat_delta=lat_delta, $
            r_delta=r_delta, arc_delta=arc_delta, seed=seed,                      $
            RK4Flag=RK4Flag, step=step, tol=tol, maxsteps=maxsteps,               $
			scottFlag=scottFlag, inclineFlag=inclineFlag, r_local=r_local,        $
			silent=silent, nthreads=nthreads, B_out=B_out, CurlB_out=CurlB_out,   $
            length_out=length_out, twist_out=twist_out, rF_out=rF_out,            $
            targetB_out=targetB_out, targetCurlB_out=targetCurlB_out,             $
            path_out=path_out, loopB_out=loopB_out, loopCurlB_out=loopCurlB_out,  $
            odir=odir, fname=fname, save_file=save_file, compress=compress,       $
			preview=preview, tmp_dir=tmp_dir, keep_tmp=keep_tmp, qsl=qsl
;------------------------------------------------------------
; check input
B3flag = N_PARAMS() eq 3
sbx=size(Bx)
if B3flag then begin
	sby=size(By) & sbz=size(Bz)
	if sbx[0] ne 3 or sby[0] ne 3 or sbz[0] ne 3 then message, 'Bx, By and Bz must be 3D arrays!'
	if sbx[1] ne sby[1] or sbx[1] ne sbz[1] or $
	   sbx[2] ne sby[2] or sbx[2] ne sbz[2] or $
	   sbx[3] ne sby[3] or sbx[3] ne sbz[3] then message, 'Bx, By and Bz must have the same dimensions!'
	nx=sbz[1] & ny=sbz[2] & nz=sbz[3]
endif else begin 
    ; the dimensions of Bx are (3,nx,ny,nz)
	if sbx[0] ne 4 or sbx[1] ne 3 then message, 'Something is wrong with the magnetic field'
	nx=sbx[2] & ny=sbx[3] & nz=sbx[4]
endelse

stretchFlag= keyword_set(xa) and keyword_set(ya) and keyword_set(za)
spherical= keyword_set(spherical)

if stretchFlag then begin
	if (nx ne n_elements(xa)) or (ny ne n_elements(ya)) or (nz ne n_elements(za)) then $
	message, 'the size of xa, ya and za must be consistant with the dimensions of the magnetic field'
endif else if spherical then message, 'xa, ya and za should be specified in spherical coordinates'

if nx lt 2 or ny lt 2 or nz lt 2 then begin
	message, 'The thickness of Bx should not be smaller than 2!'
endif else begin
	if spherical then begin
		xperiod=0
		yperiod=0
		zperiod=0
	endif else begin
		if nx eq 2 then xperiod=1
		if ny eq 2 then yperiod=1
		if nz eq 2 then zperiod=1
		xperiod=keyword_set(xperiod)
		yperiod=keyword_set(yperiod)
		zperiod=keyword_set(zperiod)
	endelse
endelse
;------------------------------------------------------------
; understand the output grid
csFlag=keyword_set(csFlag)

;  provide qsl.seed even sflag eq 0
launch_out = keyword_set(seed) and n_elements(seed) eq 1 and (size(seed,/tname) ne 'STRING') 

sflag= keyword_set(seed) and ~launch_out

if sflag then begin
	
	if (size(seed,/tname) eq 'STRING') then begin
		if seed eq 'original' then begin
			seed=fltarr(3, nx, ny, nz)
			if stretchFlag then begin
				for i=0, nx-1 do seed[0,i,*,*]=xa[i]
				for i=0, ny-1 do seed[1,*,i,*]=ya[i]
				for i=0, nz-1 do seed[2,*,*,i]=za[i]
			endif else begin
				for i=0, nx-1 do seed[0,i,*,*]=i
				for i=0, ny-1 do seed[1,*,i,*]=i
				for i=0, nz-1 do seed[2,*,*,i]=i
			endelse
		endif else if seed eq 'original_bottom' then begin
			seed=fltarr(3, nx, ny)
			if stretchFlag then begin
				for i=0, nx-1 do seed[0,i,*]=xa[i]
				for i=0, ny-1 do seed[1,*,i]=ya[i]
				seed[2,*,*]=za[0]
			endif else begin
				for i=0, nx-1 do seed[0,i,*]=i
				for i=0, ny-1 do seed[1,*,i]=i
			endelse
		endif else message, 'Something is wrong with seed'
	endif

	sz_seed=size(seed)
	if (sz_seed[1] ne 3) then message, 'Something is wrong with seed'

	out_dim=sz_seed[0]-1
	if out_dim ge 1 then nq1=sz_seed[2] else nq1=1
	if out_dim ge 2 then nq2=sz_seed[3] else nq2=1
	if out_dim eq 3 then nq3=sz_seed[4] else nq3=1
	
	bflag=0B
	vflag=0B
	cflag=0B
endif else begin

	; if spherical and xa[nx-1] is very close to 2*!pi but not equivalent, 
	; xa[nx-1] will be reset to 2*!pi.
	; then if ~preset_xreg, xreg[1] will be reset to 2*!pi
	preset_xreg= keyword_set(xreg)

	; if (.not. preset_yreg .and. (pole_south .and. pole_north) in fastqsl.x, 
	; yreg will be reset to [-!pi/2., !pi/2.]
	preset_yreg= keyword_set(yreg)
	
	if (stretchFlag) then begin
		if ~keyword_set(xreg) then xreg=[xa[0], xa[nx-1]]
		if ~keyword_set(yreg) then yreg=[ya[0], ya[ny-1]]
		if ~keyword_set(zreg) then zreg=[za[0], za[0]]
	endif else begin
		if ~keyword_set(xreg) then xreg=[0, nx-1]
		if ~keyword_set(yreg) then yreg=[0, ny-1]
		if ~keyword_set(zreg) then zreg=[0,    0]
	endelse

	if n_elements(xreg) ne 2 or n_elements(yreg) ne 2 or n_elements(zreg) ne 2 then $
	message, 'xreg, yreg and zreg must be 2-element arrays!'

	if (total([xreg[0] eq xreg[1], yreg[0] eq yreg[1], zreg[0] eq zreg[1]], /int) ge 2)  or $
	(csFlag and (zreg[0] eq zreg[1])) then message, 'Something is wrong with the cross section'

	if (spherical and csFlag) then begin
		p0= [cos(yreg[0])*[cos(xreg[0]), sin(xreg[0])], sin(yreg[0])]
		p1= [cos(yreg[1])*[cos(xreg[1]), sin(xreg[1])], sin(yreg[1])]
		; if p0 \cdot p1 eq 1 or -1, there are many great circles can pass p0 and p1
		if (abs(TRANSPOSE(p0)#p1) gt 0.999) then message, 'The great circle can not be clearly defined'
	endif
	
	if stretchFlag then zmin=za[0] else zmin=0.
	bflag = zreg[0] eq zmin and zreg[1] eq zmin
	vFlag = xreg[1] ne xreg[0] and $
			yreg[1] ne yreg[0] and $
			zreg[1] ne zreg[0] and ~csflag
	cFlag = ~(vflag or bflag)
	if vflag then out_dim=3 else out_dim=2

	; grid spacing for output
	if ~keyword_set(factor) then factor=4
	if stretchFlag then begin
		if spherical then begin
			if keyword_set(delta) then begin
				if ~keyword_set(lon_delta) then lon_delta = delta/za[0]
				if ~keyword_set(lat_delta) then lat_delta = delta/za[0]
				if ~keyword_set(  r_delta) then   r_delta = delta
			endif else begin
				if ~keyword_set(lon_delta) then lon_delta = (xa[nx-1]-xa[0])/((nx-1.0)*factor) 
				if ~keyword_set(lat_delta) then lat_delta = (ya[ny-1]-ya[0])/((ny-1.0)*factor) 
				if ~keyword_set(  r_delta) then   r_delta = (za[nz-1]-za[0])/((nz-1.0)*factor)
			endelse 
			if     ~keyword_set(arc_delta) then arc_delta = min([lon_delta, lat_delta])
		endif else if  ~keyword_set(delta) then     delta = (xa[nx-1]-xa[0])/((nx-1.0)*factor)
	endif else if      ~keyword_set(delta) then     delta = 1.0/factor
endelse
;------------------------------------------------------------
if ~keyword_set(nthreads) then nthreads = 0L
; it is surprised that !CPU.HW_NCPU is 1 in GDL on a MacOS, let's correct nthreads in fastqsl.x
; if nthreads is 0, nthreads will be reset to OMP_GET_NUM_PROCS()-2
; if nthreads gt OMP_GET_NUM_PROCS(), nthreads will be reset to OMP_GET_NUM_PROCS()

RK4Flag = keyword_set(RK4Flag)
if ~keyword_set(step) then step= 1.
if ~keyword_set(tol)  then tol = 10.0^(-4.)
if (N_ELEMENTS(maxsteps) eq 0) then $
if RK4flag then maxsteps=long(4L*(nx+ny+nz)/step) else maxsteps=4L*(nx+ny+nz)
inclineFlag=keyword_set(inclineFlag)
if ~keyword_set(r_local) then r_local=0.

B_out           = keyword_set(B_out)
CurlB_out       = keyword_set(CurlB_out)

; if ~B_out and ~CurlB_out and ~launch_out and maxsteps eq 0 then begin
; 	print, 'nothing to do, maybe invoke B_out or CurlB_out, or set seed=1, or set non-zero maxsteps'
; 	return
; endif

scottFlag       = keyword_set(scottFlag)   and (maxsteps ne 0)
twist_out       = keyword_set(twist_out)   and (maxsteps ne 0)
length_out      = keyword_set(length_out)  and (maxsteps ne 0)

rF_out          = keyword_set(rF_out)      and (maxsteps ne 0)
targetB_out     = keyword_set(targetB_out) and (maxsteps ne 0)
targetCurlB_out = keyword_set(targetCurlB_out) and (maxsteps ne 0)

if (out_dim eq 3 and keyword_set(path_out)) then print, $
'invoking path_out for a 3D ouput grid is not allowed, path_out is changed to 0'

path_out        = keyword_set(path_out) and (maxsteps ne 0) and out_dim lt 3
loopB_out       = keyword_set(loopB_out) and path_out
loopCurlB_out   = keyword_set(loopCurlB_out) and path_out

preview         = keyword_set(preview) ; or n_elements(preview) eq 0
save_file       = keyword_set(save_file)
verbose         =~keyword_set(silent)
;------------------------------------------------------------
; the directory for output

; https://learn.microsoft.com/en-us/dotnet/standard/io/file-path-formats#canonicalize-separators
; Actually Windows also accept '/', set as below just for a better readability in Windows
os_sep=PATH_SEP()

cd, current = cdir
IF STRMID(cdir, STRLEN(cdir)-1) NE os_sep THEN cdir=cdir+os_sep

; if preview or save_file then begin
	if keyword_set(odir) then begin
		IF STRMID(odir, STRLEN(odir)-1) NE os_sep THEN odir=odir+os_sep
	endif else odir= cdir+'fastqsl'+os_sep
	if ~file_test(odir) then file_mkdir, odir
; endif
;------------------------------------------------------------
; the temporary directory for the data transmission between fastqsl.x and fastqsl.pro
if keyword_set(tmp_dir) then begin
	IF STRMID(tmp_dir, STRLEN(tmp_dir)-1) NE os_sep THEN tmp_dir=tmp_dir+os_sep
endif else tmp_dir= cdir+'tmpFastQSL'+os_sep

old_tmp_dir=file_test(tmp_dir)
if old_tmp_dir then begin
	dummy=file_search(tmp_dir, '*.bin', count=nf)
	if nf gt 0 then file_delete, dummy
endif else file_mkdir, tmp_dir
;------------------------------------------------------------
;  transmit data to fastqsl.x
get_lun, unit

openw,  unit, tmp_dir+'head.bin'
writeu, unit, float([step, tol, r_local]), $
              long([maxsteps, RK4Flag, inclineFlag, $
              launch_out, B_out, CurlB_out, length_out, twist_out, $
		      rF_out, targetB_out, targetCurlB_out, $
			  path_out, loopB_out, loopCurlB_out, $
			  sflag, bflag, cflag, vflag, nthreads, scottFlag, verbose])
close,  unit

openw,  unit, tmp_dir+'bfield.bin'
writeu, unit, long([nx, ny, nz, B3flag, spherical, $
preview, xperiod, yperiod, zperiod])
writeu, unit, float(Bx)
if B3flag then writeu, unit, float(By), float(Bz)
close,  unit

if (stretchFlag) then begin
	openw,  unit, tmp_dir+'axis.bin'
	writeu, unit, float(xa), float(ya), float(za)
	close,  unit
endif

if sflag then begin
	openw,  unit, tmp_dir+'dim_seed.bin'
	writeu, unit, long([nq1,nq2,nq3])
	close,  unit
	openw,  unit, tmp_dir+'seed.bin'
	writeu, unit, float(seed)
	close,  unit
endif else begin
	if spherical then deltas=[arc_delta, lon_delta, lat_delta, r_delta] $
	                 else deltas=fltarr(4)+delta

	openw,  unit, tmp_dir+'head_region.bin'
	writeu, unit, float([xreg, yreg, zreg, deltas]), long([csflag, preset_xreg, preset_yreg])
	close,  unit
endelse
;------------------------------------------------------------
; computed by fastqsl.x
cd, tmp_dir
; please specify the path
spawn, '~/Desktop/QSLS/update/fastqsl.x'
; spawn, '/data/QSLS/update/fastqsl.x'
; spawn, '/path/of/fastqsl.x'
cd, cdir
; ################################### retrieving results ######################################
; make the structure QSL
if (out_dim eq 0) then begin
	array_float='0.'
	array_byte='0B'
	array_long='0L'
	array_ptr='ptr_new()'
endif else begin
	; if out_dim eq 2, fltarr(nq1, nq2, nq3) will become fltarr(nq1, nq2) by IDL
	array_float='fltarr(nq1, nq2, nq3)'
	array_byte='bytarr(nq1, nq2, nq3)'
	array_long='lonarr(nq1, nq2, nq3)'
	array_ptr='ptrarr(nq1, nq2, nq3)'
endelse
array_vfloat='fltarr(3, nq1, nq2, nq3)'

qsl_data0=[ $
['axis1', 'fltarr(2, nq1)'], $
['seed', array_vfloat], $
['q', array_float], $
['q_perp', array_float], $
['q_local', array_float], $
['rboundary', array_byte], $
['sign2d', 'intarr(nq1, nq2)'], $
['length', array_float], $
['twist', array_float], $
['B', array_vfloat], $
['CurlB', array_vfloat], $
['rFs', array_vfloat], $
['rFe', array_vfloat], $
['Bs', array_vfloat], $
['Be', array_vfloat], $
['CurlBs', array_vfloat], $
['CurlBe', array_vfloat], $o
['path', array_ptr], $
['loopB', array_ptr], $
['loopCurlB', array_ptr], $
['index_seed', array_long]]

if ~sflag then begin
	nq1=0L & nq2=0L & nq3=0L & normal_index=0L
	xreg=[0.,0.] & yreg=[0.,0.] & zreg=[0.,0.]
	openr,  unit, tmp_dir+'tail_region.bin'
	readu,  unit, nq1, nq2, nq3, normal_index, xreg, yreg, zreg
	close,  unit
endif

openw,  unit, tmp_dir+'qsl_structure.pro'
printf, unit, 'pro qsl_structure, QSL, nq1, nq2, nq3, xreg, yreg, zreg, $'
printf, unit, 'delta, arc_delta, lon_delta, lat_delta, r_delta, step, tol, r_local'

printf, unit, 'get_lun, unit'

if path_out then begin
	printf, unit, 'n_loops=long64(nq1)*nq2*nq3'
	printf, unit, 'indexes=lon64arr(n_loops+1)'
	printf, unit, 'openr, unit, "'+tmp_dir+'indexes.bin"'
	printf, unit, 'readu, unit, indexes'
	printf, unit, 'close, unit'
	printf, unit, 'dummy=fltarr(3,indexes[n_loops])'
endif

n_data=0
for i=0, n_elements(qsl_data0)/2-1 do begin
	name=qsl_data0[0,i]
	file=tmp_dir+name+'.bin'
	if file_test(file) then begin
		n_data=n_data+1
		if n_data eq 1 then qsl_data=name else qsl_data=[qsl_data,name]
		printf, unit, name+'='+qsl_data0[1,i]
		printf, unit, 'openr, unit, "'+file+'"'
		if (qsl_data0[1,i] eq array_ptr) then begin
			printf, unit, 'readu, unit, dummy'
			printf, unit, 'for i=0L, n_loops-1 do $'
			printf, unit, name+'[i]=ptr_new(dummy[*,indexes[i]:indexes[i+1]-1])'		
		endif else printf, unit, 'readu, unit, '+name
		printf, unit, 'close, unit'
	endif
endfor

printf, unit, 'free_lun, unit, /force'

strs='QSL={ $'

if maxsteps ne 0 then $
if rk4flag then strs=[strs, 'step:float(step), $'] else strs=[strs, 'tol:float(tol), $'] 

if ~sflag then begin
	if spherical then begin
		strs=[strs, 'lon_reg:xreg, $','lat_reg:yreg, $','r_reg:zreg, $']
		if normal_index eq -1 then strs=[strs, 'arc_delta:float(arc_delta), $']
		if normal_index eq 0 or normal_index eq 2 then strs=[strs, 'lat_delta:float(lat_delta), $']
		if normal_index eq 1 or normal_index eq 2 then strs=[strs, 'lon_delta:float(lon_delta), $']
		if (normal_index ne 2) or vflag then strs=[strs, 'r_delta:float(r_delta), $']
	endif else strs=[strs, 'xreg:xreg, $','yreg:yreg, $','zreg:zreg, $','delta:float(delta), $']
endif

if file_test(tmp_dir+'q_local.bin') then strs=[strs,'r_local:float(r_local), $']

for i=0, n_data-1 do strs=[strs, qsl_data[i]+':'+qsl_data[i]+', $']

n_strs=n_elements(strs)
strs[n_strs-1]=(strsplit(strs[n_strs-1],',',/extract))[0]+'}'
for i=0, n_strs-1 do printf, unit, strs[i]

printf, unit, 'end'
close,  unit

cd, tmp_dir
resolve_routine, 'qsl_structure'
cd, cdir

qsl_structure, QSL, nq1, nq2, nq3, xreg, yreg, zreg, $
delta, arc_delta, lon_delta, lat_delta, r_delta, step, tol, r_local
;------------------------------------------------------------	
; the name of .sav file
if ~keyword_set(fname) then begin

    if sflag then head_str='seed'
    if bflag then head_str='bottom'
    if cflag then head_str='cs'
    if vflag then head_str='volume'
	
	if ~(spherical or sflag) then begin
		decimal3=round(1000.*(delta-floor(delta)))
		for i = 1, 4 do if (decimal3 mod 10^i) ne 0 then break
		; string(0.99,'(i0)')='0'; string(round(0.99), '(i0)')='1'
		if i ge 4 then delta_str='_delta'+string(round(delta), '(i0)') $
		          else delta_str='_delta'+string(delta, '(f0.'+string(4-i)+')')
	endif else delta_str=''

	if cFlag and ~csflag then begin

		case normal_index of
			0: begin
				cut_coordinate=xreg[0]
				if spherical then cut_str0='_lon' else cut_str0='_x'
			end
			1: begin
				cut_coordinate=yreg[0]
				if spherical then cut_str0='_lat' else cut_str0='_y'
			end
			2: begin
				cut_coordinate=zreg[0]
				if spherical then cut_str0='_r'   else cut_str0='_z'
			end
		endcase 

		decimal3=round(1000.*(cut_coordinate -floor(cut_coordinate)))
		for i = 1, 4 do if (decimal3 mod 10^i) ne 0 then break

		if i ge 4 then cut_str=cut_str0+string(round(cut_coordinate), '(i0)') $
		          else cut_str=cut_str0+string(cut_coordinate, '(f0.'+string(4-i)+')')

	endif else cut_str=''

	fname = head_str + delta_str + cut_str
endif

if save_file then save, filename=odir+fname+'.sav', QSL, compress=compress
;------------------------------------------------------------
if preview then begin
if verbose then begin
	print, ''
	print, 'These images are produced:'
endif

; mark the area for calculation or plot seed and their field lines on the magnetogram
if stretchFlag then begin
	nx_mag=0L & ny_mag=0L & mag_delta=0.0
	openr, unit, tmp_dir+'magnetogram.bin'
	readu, unit, nx_mag, ny_mag, mag_delta
	magnetogram=fltarr(nx_mag, ny_mag)
	readu, unit, magnetogram
	close, unit
endif else begin
	if B3flag then magnetogram=Bz[*,*,0] else magnetogram=reform(Bx[2,*,*,0])
	nx_mag=nx
	ny_mag=ny
endelse

cur_device=!D.name
SET_PLOT, 'Z'
DEVICE, SET_RESOLUTION=[nx_mag, ny_mag]

abnormal = WHERE(~FINITE(magnetogram))
if (abnormal[0] ne -1) then magnetogram[abnormal]=0.
scale_top= max(abs(magnetogram))/4.0 < 1000.0
tv, bytscl(temporary(magnetogram), min=-scale_top, max=scale_top, top=253B)

if spherical then two_pi=2*!pi
if sflag then begin
	if (out_dim le 1) then begin
	
		x_seed=reform(seed[0,*])
		y_seed=reform(seed[1,*])
		if stretchFlag then begin
			x_seed=x_seed-xa[0]
			if spherical then x_seed= ((x_seed mod two_pi)+two_pi) mod two_pi
			x_seed=x_seed/mag_delta
			y_seed=(y_seed-ya[0])/mag_delta
		endif
		plots, x_seed, y_seed, /dev, PSym=2, color=254B
		
		if path_out then begin
			for i=0, nq1-1 do begin
				x_path=reform((*QSL.path[i])[0,*])
				y_path=reform((*QSL.path[i])[1,*])
				if stretchFlag then begin
					x_path=x_path-xa[0]
					; remap x_path to be in [0, 2*!pi)
					if spherical then x_path= ((x_path mod two_pi)+two_pi) mod two_pi
					x_path=x_path/mag_delta
					y_path=(y_path-ya[0])/mag_delta
				endif
				plots, x_path, y_path, /dev, PSym=3, color=255B
			endfor
		endif

	endif else begin
		; seed[0,*,0]=seed[0,*,0,0], even when out_dim is 3
		x_margin=[reform(seed[0,*,0]    ), $
			      reform(seed[0,0,*]    ), $
			      reform(seed[0,*,nq2-1]), $
		          reform(seed[0,nq1-1,*])]
				   
		y_margin=[reform(seed[1,*,0]    ), $
		          reform(seed[1,0,*]    ), $
		          reform(seed[1,*,nq2-1]), $
		          reform(seed[1,nq1-1,*])]
	endelse
endif else begin
	if CSflag then begin
		if spherical then begin
			x_margin=reform(QSL.axis1(0,*))
			y_margin=reform(QSL.axis1(1,*))
		endif else begin
			x_margin=xreg
			y_margin=yreg
		endelse
	endif else begin
		x_margin=[xreg[0],xreg[1],xreg[1],xreg[0],xreg[0]]
		y_margin=[yreg[0],yreg[0],yreg[1],yreg[1],yreg[0]]
	endelse
endelse

if (out_dim ge 2) then begin
	if (spherical and CSflag) then PSym=3 else PSym=0

	if stretchFlag then begin

		; some regions may stand across at longitude = 0 or 2*!pi
		if spherical then begin
			x_margin1=(x_margin-two_pi-xa[0])/mag_delta
			x_margin2=(x_margin+two_pi-xa[0])/mag_delta
		endif

		x_margin=(x_margin-xa[0])/mag_delta
		y_margin=(y_margin-ya[0])/mag_delta
	endif
		
	; prevent missing margin line from tiny numerical errors
	if spherical then begin
		abnormal=where(y_margin lt 0. and y_margin ge -1.) 
		if (abnormal[0] ne -1) then y_margin[abnormal]=0.
		abnormal=where(y_margin gt ny_mag-1. and y_margin le ny_mag+2.) 
		if (abnormal[0] ne -1) then y_margin[abnormal]=ny_mag-1.
	endif else begin
		abnormal=where(x_margin gt nx_mag-1. and x_margin le nx_mag) 
		if (abnormal[0] ne -1) then x_margin[abnormal]=nx_mag-1.
		abnormal=where(y_margin gt ny_mag-1. and y_margin le ny_mag) 
		if (abnormal[0] ne -1) then y_margin[abnormal]=ny_mag-1.
	endelse
		
	plots, x_margin , y_margin, /dev, PSym=PSym, color=254B
	if (spherical) then begin
		plots, x_margin1, y_margin, /dev, PSym=PSym, color=254B
		plots, x_margin2, y_margin, /dev, PSym=PSym, color=254B
	endif
endif

r_mag=[bindgen(254)+1B, 255B,   0B]
g_mag=[bindgen(254)+1B,   0B, 255B]
b_mag=[bindgen(254)+1B,   0B,   0B]
write_png, odir+fname+'_magnetogram.png', TVRD(), r_mag, g_mag, b_mag
set_plot, cur_device

if verbose then print, odir+fname+'_magnetogram.png'
;------------------------------------------------------------
; preview q/q_perp, length, twist

; if vflag and the bottom plane is included, 'sign2d.bin' also can be found
plot_bottom=file_test(tmp_dir+'sign2d.bin')

if (maxsteps ne 0 and (out_dim eq 2 or plot_bottom)) then begin

	r_doppler= [bindgen(127)*2B, REPLICATE(255B, 129)]
	b_doppler= REVERSE(r_doppler)
	g_doppler= [r_doppler[0:127], b_doppler[128:255]]
	
	r_boundary=[0B, 128B, 255B, 192B, 255B,   0B,   0B,   0B, 128B, REPLICATE(255B, 247)]
	g_boundary=[0B, 128B, 255B,   0B, 128B, 128B, 128B,   0B,   0B, REPLICATE(255B, 247)]
	b_boundary=[0B, 128B, 255B,   0B,   0B,   0B, 255B, 128B, 128B, REPLICATE(255B, 247)]

	rb_tmp= qsl.rboundary[*,*,0] ; This approach works whether out_dim eq 2 or 3
	rbs=rb_tmp/10B
	rbe=rb_tmp mod 10B
	open_mark=where(rb_tmp ne 11)

	if plot_bottom then begin
		rb_target=bytarr(nq1,nq2)
		dummy=where(qsl.sign2d eq 1)
		if (dummy[0] ne -1) then rb_target[dummy]=rbe[dummy]
		dummy=where(qsl.sign2d eq -1)
		if (dummy[0] ne -1) then rb_target[dummy]=rbs[dummy]
		dummy=where(qsl.sign2d eq 0)
		if (dummy[0] ne -1) then begin
			dummy2=where(rbs[dummy] eq rbe[dummy])
			if (dummy2[0] ne -1) then rb_target[dummy[dummy2]]=rbs[dummy[dummy2]]
		endif

		write_png, odir+fname+'_rb_target.png', rb_target, r_boundary, g_boundary, b_boundary
		if verbose then print, odir+fname+'_rb_target.png'

	endif else begin
		write_png, odir+fname+'_rbs.png', rbs, r_boundary, g_boundary, b_boundary
		write_png, odir+fname+'_rbe.png', rbe, r_boundary, g_boundary, b_boundary

		if verbose then print, odir+fname+'_rbs.png'
		if verbose then print, odir+fname+'_rbe.png'
		
		;  gray for closed field lines
		closed=bytarr(nq1, nq2)+128B              

		; white for opened field lines
		if (open_mark[0] ne -1) then closed[open_mark]=255B 
		
		write_png, odir+fname+'_mark_closed.png', closed
		if verbose then print, odir+fname+'_mark_closed.png'
	endelse

	; q/q_perp
	for i=0, n_data-1 do begin
		if qsl_data[i] eq 'q' or qsl_data[i] eq 'q_perp' or qsl_data[i] eq 'q_local' then begin
			q_str=qsl_data[i]
			if q_str eq 'q'       then q_tmp=qsl.q[*,*,0]
			if q_str eq 'q_perp'  then q_tmp=qsl.q_perp[*,*,0]
			if q_str eq 'q_local' then q_tmp=qsl.q_local[*,*,0]

			abnormal=WHERE(~FINITE(q_tmp))
			if (abnormal[0] ne -1) then q_tmp[abnormal]=1. ; 1. for white color

			if plot_bottom then begin
				im=bytscl(alog10(q_tmp)*qsl.sign2d, min=-5, max=5)
				write_png, odir+fname+'_slog'+q_str+'.png', im, r_doppler, g_doppler, b_doppler
				if verbose then print, odir+fname+'_slog'+q_str+'.png'

				if targetB_out and q_str eq 'q' then q_tmp1=q_tmp

				; if (open_mark[0] ne -1) then im[open_mark]=127B
				; write_png, odir+fname+'_slog'+q_str+'_orig.png', im, r_doppler, g_doppler, b_doppler
				; if verbose then print, odir+fname+'_slog'+q_str+'_orig.png'
			endif else begin
				im=bytscl(alog10(q_tmp), min=1, max=5)
				write_png, odir+fname+'_log'+q_str+'.png', im
				if verbose then print, odir+fname+'_log'+q_str+'.png'

				; if (open_mark[0] ne -1) then im[open_mark]=0B
				; write_png, odir+fname+'_log'+q_str+'.png', im
				; if verbose then print, odir+fname+'_log'+q_str+'_orig.png'
			endelse
		endif
	endfor

	if targetB_out and plot_bottom and (n_elements(q_tmp1) ne 0) then begin

		; In Titov (2007), q = N^2 / Delta, and Delta = Bnr
		Bnr=fltarr(nq1,nq2)

		; just for saving time in the following do loops
		bs=qsl.bs[*,*,*,0]
		be=qsl.be[*,*,*,0]
		sign2d=qsl.sign2d

		for j=0, nq2-1 do begin
		for i=0, nq1-1 do begin
			if            sign2d[i,j] eq  1 and rbe[i,j] ge 1 and rbe[i,j] le 6 then begin
				Bn_target=be[(6-rbe[i,j])/2,i,j]
				if Bn_target ne 0. then Bnr[i,j]=abs(bs[2,i,j]/Bn_target)
			endif else if sign2d[i,j] eq -1 and rbs[i,j] ge 1 and rbs[i,j] le 6 then begin
				Bn_target=bs[(6-rbs[i,j])/2,i,j]
				if Bn_target ne 0. then Bnr[i,j]=abs(be[2,i,j]/Bn_target)
			endif
		endfor
		endfor

		N=sqrt(q_tmp1 * Bnr) ; Priest and Demoulin (1995)

		abnormal=where(N eq 0.0 or ~FINITE(N))
		if (abnormal[0] ne -1) then N[abnormal]=1.

		im=bytscl(alog10(N), min=0.5, max=3, /nan)
		write_png, odir+fname+'_N.png', im
		if verbose then print, odir+fname+'_N.png'

		; 1. for white color, should be done after N is calculated
		abnormal=where(Bnr eq 0.0 or ~FINITE(Bnr))
		if (abnormal[0] ne -1) then Bnr[abnormal]=1. 

		im=bytscl(alog10(Bnr), min=-2, max=2)
		write_png, odir+fname+'_lg_Bnr.png', im, r_doppler, g_doppler, b_doppler
		if verbose then print, odir+fname+'_lg_Bnr.png'
	endif

	if length_out and out_dim eq 2 then begin
		if stretchFlag then length_top=2.*(za[nz-1]-za[0]) else length_top=2.*(nz-1)
		im=bytscl(qsl.length, min=0, max=length_top, /nan)
		write_png, odir+fname+'_length.png', im
		if verbose then print, odir+fname+'_length.png'
	endif

	if twist_out and out_dim eq 2 then begin
		twist_tmp=qsl.twist
		abnormal=WHERE(~FINITE(twist_tmp))
		if (abnormal[0] ne -1) then twist_tmp[abnormal]=0. ; 0. for white color
		im=bytscl(twist_tmp, min=-2, max=2)
		write_png, odir+fname+'_twist.png', im, r_doppler, g_doppler, b_doppler
		if verbose then print, odir+fname+'_twist.png'
	endif
endif
endif ; preview
;------------------------------------------------------------
free_lun, unit, /force

if ~keyword_set(keep_tmp) then begin
	if old_tmp_dir then begin
		file_delete, tmp_dir+'qsl_structure.pro'
		if n_data gt 0 then file_delete, tmp_dir+qsl_data+'.bin'
		if path_out then file_delete, tmp_dir+'indexes.bin'
		if ~sflag then file_delete, tmp_dir+'tail_region.bin'
		if preview and stretchFlag then file_delete, tmp_dir+'magnetogram.bin'
	endif else file_delete, tmp_dir, /recursive
endif

if (verbose) then begin
	print, ''
	print, 'Elements in qsl:'
	names=tag_names(qsl)
	for i=0, n_tags(qsl)-1 do begin
		tname_i=size(qsl.(i),/tname)
		type_tail=string('', '(A'+string(8-strlen(tname_i),'(A)')+')')
		n_i=n_elements(qsl.(i))
		if n_i gt 3 or (tname_i eq 'POINTER' and n_i gt 1) then begin
			sz_i=size(qsl.(i))
			content='  Array['
			for j=1, sz_i[0]-1 do content=content+string(sz_i[j], '(i0)')+', '
			content=content+string(sz_i(sz_i[0]), '(i0)')+']'
		endif else content=qsl.(i)
		name_tail=string('', '(A'+string(15-strlen(names[i]),'(A)')+')')
		print, 'qsl.'+strlowcase(names[i])+name_tail+tname_i+type_tail, content
	endfor
	if save_file then begin
		print, "Try:"
		print, "restore, '"+odir+fname+".sav'"
	endif
endif

END
