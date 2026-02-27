import numpy as np
import pickle, os
from fastqsl import fastqsl 
from convert_coordinate import convert_coordinate

def bfield_charge4_cartersian(bfile, ixh=45, nz=61, xmax=2.5):
    nx=ixh*2+1
    xa=np.linspace(-xmax, xmax, nx)
    B_delta=xmax/ixh
    ya=xa.copy()
    za=np.linspace(0., (nz-1)*B_delta, nz)

    Bx=np.zeros((nz, nx, nx),'f4')
    By=np.zeros((nz, nx, nx),'f4')
    Bz=np.zeros((nz, nx, nx),'f4')

    ri=np.array([[-1.5,0.0,-0.5], \
                 [-0.5,0.0,-0.5], \
                 [ 0.5,0.0,-0.5], \
                 [ 1.5,0.0,-0.5]],'f4')
    qi=np.array([1.,-1.,1.,-1.],'f4')

    for k in range(nz):
        for j in range(ixh+1):
            for i in range(ixh+1):
                rp=np.array([xa[i], ya[j], za[k]],'f4')
                Bp=np.zeros((3),'f4')
                for s in range(4):
                    dr=rp-ri[s,:]
                    Bp=Bp+qi[s]*dr/(np.sum(dr**2))**1.5
                Bx[k,j,i]=Bp[0]
                By[k,j,i]=Bp[1]
                Bz[k,j,i]=Bp[2]
    Bx[:,0:ixh+1,ixh+1:nx]= np.flip(Bx[:,0:ixh+1,0:ixh],2)
    By[:,0:ixh+1,ixh+1:nx]=-np.flip(By[:,0:ixh+1,0:ixh],2)
    Bz[:,0:ixh+1,ixh+1:nx]=-np.flip(Bz[:,0:ixh+1,0:ixh],2)

    Bx[:,ixh+1:nx,:]= np.flip(Bx[:,0:ixh,:],1)
    By[:,ixh+1:nx,:]=-np.flip(By[:,0:ixh,:],1)
    Bz[:,ixh+1:nx,:]= np.flip(Bz[:,0:ixh,:],1)

    Bvec=np.zeros((nz, nx, nx,3),'f4')
    Bvec[:,:,:,0]=Bx
    Bvec[:,:,:,1]=By
    Bvec[:,:,:,2]=Bz

    with open(bfile, 'wb') as file: pickle.dump((Bvec, Bx, By, Bz, xa, ya, za), file)
# ------------------------------------------------------------
def bfield_charge4_spherical(bfile, n_lon=181, n_lat=91, n_r=31):

    B_lon=np.zeros((n_r, n_lat, n_lon),'f4')
    B_lat=np.zeros((n_r, n_lat, n_lon),'f4')
    B_r  =np.zeros((n_r, n_lat, n_lon),'f4')

    ri=np.array([[0.7*np.pi, 0.0, 0.8], \
                 [0.9*np.pi, 0.1, 0.8], \
                 [1.1*np.pi,-0.1, 0.8], \
                 [1.3*np.pi, 0.1, 0.8]],'f4')

    ri_car=ri.copy()
    for s in range(4): ri_car[s,:]= ri[s,2]* \
    np.array([np.cos(ri[s,1])*np.cos(ri[s,0]), np.cos(ri[s,1])*np.sin(ri[s,0]), np.sin(ri[s,1])],'f4')

    qi=np.array([1.,-1.,1.,-1.],'f4')

    lon_rad= np.linspace(0., 2.*np.pi, n_lon)
    lat_rad= np.linspace(-np.pi/2., np.pi/2., n_lat)
    radius = np.linspace(1., 2.5, n_r)

    cos_lon=np.cos(lon_rad)
    sin_lon=np.sin(lon_rad)
    cos_lat=np.cos(lat_rad)
    sin_lat=np.sin(lat_rad)
    
    for i in range(n_lon):
        e_lon=np.array([-sin_lon[i], cos_lon[i], 0.],'f4')
        for j in range(n_lat):
            e_lat=np.array([-sin_lat[j]*cos_lon[i], -sin_lat[j]*sin_lon[i], cos_lat[j]],'f4')
            e_r  =np.array([ cos_lat[j]*cos_lon[i],  cos_lat[j]*sin_lon[i], sin_lat[j]],'f4')
            for k in range(n_r):
                rp_car=radius[k]*e_r
                Bp_car=np.zeros((3),'f4')
                for s in range(4):
                    dr=rp_car-ri_car[s,:]
                    Bp_car=Bp_car+qi[s]*dr/(np.sum(dr**2)**1.5)
            
                B_lon[k,j,i]=np.dot(e_lon,Bp_car)
                B_lat[k,j,i]=np.dot(e_lat,Bp_car)
                B_r  [k,j,i]=np.dot(e_r  ,Bp_car)
    with open(bfile, 'wb') as file: pickle.dump((B_lon, B_lat, B_r, lon_rad, lat_rad, radius), file)
# ------------------------------------------------------------
# this way gives the identical field as bfield_charge4_spherical, but is slower
def bfield2_charge4_spherical(bfile, n_lon=181, n_lat=91, n_r=31):

    B_lon=np.zeros((n_r, n_lat, n_lon),'f4')
    B_lat=np.zeros((n_r, n_lat, n_lon),'f4')
    B_r  =np.zeros((n_r, n_lat, n_lon),'f4')

    ri=np.array([[0.7*np.pi, 0.0, 0.8], \
                 [0.9*np.pi, 0.1, 0.8], \
                 [1.1*np.pi,-0.1, 0.8], \
                 [1.3*np.pi, 0.1, 0.8]],'f4')

    ri_car= convert_coordinate(ri, mode='lon_lat_r_to_xyz')

    qi=np.array([1.,-1.,1.,-1.],'f4')

    lon_rad= np.linspace(0., 2.*np.pi, n_lon)
    lat_rad= np.linspace(-np.pi/2., np.pi/2., n_lat)
    radius = np.linspace(1., 2.5, n_r)

    cos_lon=np.cos(lon_rad)
    sin_lon=np.sin(lon_rad)
    cos_lat=np.cos(lat_rad)
    sin_lat=np.sin(lat_rad)

    seed_original=np.zeros((n_r, n_lat, n_lon, 3),'f4')
    for i in range(n_lon): seed_original[:,:,i,0]=lon_rad[i]
    for i in range(n_lat): seed_original[:,i,:,1]=lat_rad[i]
    for i in range(n_r):   seed_original[i,:,:,2]= radius[i]
    seed_car= convert_coordinate(seed_original, mode='lon_lat_r_to_xyz')

    B_car=np.zeros((n_r, n_lat, n_lon, 3),'f4')
    for k in range(n_r):
        for j in range(n_lat):
            for i in range(n_lon):
                for s in range(4):
                    dr=seed_car[k,j,i,:]-ri_car[s,:]
                    B_car[k,j,i,:]=B_car[k,j,i,:]+qi[s]*dr/(np.sum(dr**2)**1.5)

    (seed_original, Bvec)= convert_coordinate(seed_car, B_car, mode='xyz_to_lon_lat_r')
    B_lon=Bvec[:,:,:,0]
    B_lat=Bvec[:,:,:,1]
    B_r  =Bvec[:,:,:,2]
    with open(bfile, 'wb') as file: pickle.dump((B_lon, B_lat, B_r, lon_rad, lat_rad, radius), file)
# ------------------------------------------------------------
# Examples for cartersian grid
bfile='charge4_cartersian.pkl'
cdir = os.getcwd()+os.sep 
if bfile not in os.listdir(cdir): bfield_charge4_cartersian(bfile)
with open(bfile, 'rb') as file: (Bvec, Bx, By, Bz, xa, ya, za) = pickle.load(file)

nx=len(xa)
ixh=(nx-1)/2
kend=len(za)-1

# images of Figure 4 in Zhang, P., Chen, J.*, Liu, R. and Wang, C., 2022, ApJ, 937, 26
fastqsl(Bx, By, Bz, fname='method1_z0', preview=True)
fastqsl(Bx, By, Bz, fname='method2_z0', preview=True, scottFlag=True, RK4Flag=True)
fastqsl(Bvec, xreg=[0,nx-1], yreg=[ixh,ixh], zreg=[0,kend/2], \
        fname='method1_y0', preview=True)
fastqsl(Bx, By, Bz, xreg=[0,nx-1], yreg=[ixh,ixh], zreg=[0,kend/2], \
        fname='method2_y0', preview=True, scottFlag=True)

# An example of calculating in a cross section which is tilted to x-axis and y-axis, 
# and with stretched (actually uniformed) grids
fastqsl(Bx, By, Bz, xa=xa, ya=ya, za=za, delta=0.01, \
xreg=[-2,0], yreg=[1,0], zreg=[0,2], csFlag=True, \
fname='tilted_cs', RK4Flag=True, step=2.0, odir= 'fastqsl/', twist_out=True, preview=True)

# An example of calculating in a box volume
fastqsl(Bx, By, Bz, xreg=[ixh/2,ixh], yreg=[ixh/4,ixh], zreg=[kend/4, kend/2], \
delta=0.8, tol=1.0e-3, odir='fastqsl', nthreads=12, preview=True)

# only exporting CurlB
qsl=fastqsl(Bx, By, Bz, CurlB_out=True, maxsteps=0, seed='original', odir= 'fastqsl', save_file=True)

# Figure 4 of Chen (2026), see fname+'_logq_local.png'
fastqsl(Bx, By, Bz, xa=xa, ya=ya, za=za, \
xreg=[-1.7, 1.7], yreg=[0,0], zreg=[0,1.3], preview=True, r_local=2, fname='r_local2')
fastqsl(Bx, By, Bz, xa=xa, ya=ya, za=za, \
xreg=[-1.7, 1.7], yreg=[0,0], zreg=[0,1.3], preview=True, r_local=1, fname='r_local1')
fastqsl(Bx, By, Bz, xa=xa, ya=ya, za=za, \
xreg=[-1.7, 1.7], yreg=[0,0], zreg=[0,1.3], preview=True, r_local=0.2, fname='r_local0.2')
fastqsl(Bx, By, Bz, xa=xa, ya=ya, za=za, \
xreg=[-1.7, 1.7], yreg=[0,0], zreg=[0,1.3], preview=True, r_local=0.05, fname='r_local0.05')
# ------------------------------------------------------------
# Examples for spherical grid
bfile='charge4_spherical.pkl'
cdir = os.getcwd()+os.sep 
if bfile not in os.listdir(cdir): bfield_charge4_spherical(bfile)
with open(bfile, 'rb') as file: (b_lon, b_lat, b_r, lon_rad, lat_rad, radius) = pickle.load(file)

# Q at r=1
fastqsl(b_lon, b_lat, b_r, xa=lon_rad, ya=lat_rad, za=radius, spherical=True, \
factor=4, preview=True, \
xreg=[0.,2*np.pi], yreg=[-np.pi/2, np.pi/2])


# trace field lines from two points
qsl=fastqsl(b_lon, b_lat, b_r, xa=lon_rad, ya=lat_rad, za=radius, spherical=True, \
fname='spherical_seed_path', preview=True, \
seed=[[np.pi*0.85, 0.1, 1.], [np.pi*1.1, -0.2, 1.2]], path_out=True, loopB_out=True)

# transfrom *qsl['path'][1] in (longitude, latitude, radius() to path1_car in (x, y, z), and
# transfrom *qsl['loopB'][1] in (B_lon, B_lat, B_r) to loopB1_car in (B_x, B_y, B_z)
path1_car, loopB1_car= convert_coordinate(qsl['path'][1], qsl['loopB'][1], mode='lon_lat_r_to_xyz')