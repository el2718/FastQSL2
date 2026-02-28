pro bfield_charge4_cartersian, bfile, ixh=ixh, nz=nz, xmax=xmax

if ~keyword_set(ixh)  then ixh =45
if ~keyword_set(nz)   then nz  =61
if ~keyword_set(xmax) then xmax=2.25

nx=ixh*2+1
B_delta=xmax/ixh

Bx=fltarr(nx,nx,nz)
By=fltarr(nx,nx,nz)
Bz=fltarr(nx,nx,nz)

ri=[[-1.5,0.0,-0.5], $
    [ 1.5,0.0,-0.5], $
    [ 0.5,0.0,-0.5], $
    [-0.5,0.0,-0.5]]

qi=[1.,-1.,1.,-1.]

xa=findgen(nx)*B_delta-ixh*B_delta
ya=xa
za=findgen(nz)*B_delta


for k=0, nz-1 do begin
for j=0,ixh do begin
for i=0,ixh do begin

	rp=[xa[i], ya[j], za[k]]
	Bp=fltarr(3)
	for s=0,3 do begin
		dr=rp-ri[*,s]
		Bp=Bp+qi[s]*dr/(total(dr^2.))^1.5
	endfor

	Bx[i,j,k]=Bp[0]
	By[i,j,k]=Bp[1]
	Bz[i,j,k]=Bp[2]
endfor
endfor
endfor

Bx[ixh+1:nx-1,0:ixh,*]= REVERSE(Bx[0:ixh-1,0:ixh,*],1)
By[ixh+1:nx-1,0:ixh,*]=-REVERSE(By[0:ixh-1,0:ixh,*],1)
Bz[ixh+1:nx-1,0:ixh,*]=-REVERSE(Bz[0:ixh-1,0:ixh,*],1)

Bx[*,ixh+1:nx-1,*]= REVERSE(Bx[*,0:ixh-1,*],2)
By[*,ixh+1:nx-1,*]=-REVERSE(By[*,0:ixh-1,*],2)
Bz[*,ixh+1:nx-1,*]= REVERSE(Bz[*,0:ixh-1,*],2)

Bvec=fltarr(3,nx,nx,nz)
Bvec[0,*,*,*]=Bx
Bvec[1,*,*,*]=By
Bvec[2,*,*,*]=Bz

save, filename=bfile, Bvec, Bx, By, Bz, xa, ya, za
end
;------------------------------------------------------------
pro bfield_charge4_spherical, bfile, n_lon=n_lon, n_lat=n_lat, n_r=n_r

if ~keyword_set(n_lon) then n_lon=181
if ~keyword_set(n_lat) then n_lat=91
if ~keyword_set(n_r)   then n_r=31

B_lon=fltarr(n_lon, n_lat, n_r)
B_lat=fltarr(n_lon, n_lat, n_r)
B_r  =fltarr(n_lon, n_lat, n_r)

ri=[[0.7*!pi, 0.0, 0.8], $
    [0.9*!pi, 0.1, 0.8], $
    [1.1*!pi,-0.1, 0.8], $
    [1.3*!pi, 0.1, 0.8]]   

ri_car=fltarr(3, 4)
for s=0, 3 do ri_car[*,s]= ri[2,s]* $
[cos(ri[1,s])*[cos(ri[0,s]), sin(ri[0,s])], sin(ri[1,s])]

qi=[1.,-1.,1.,-1.]

lon_rad=findgen(n_lon)/(n_lon-1.)*2.*!pi
lat_rad=findgen(n_lat)/(n_lat-1.)*!pi-!pi/2.
radius=findgen(n_r)/(n_r-1.)*1.5+1.

cos_lon=cos(lon_rad)
sin_lon=sin(lon_rad)
cos_lat=cos(lat_rad)
sin_lat=sin(lat_rad)

for i=0, n_lon-1 do begin
    e_lon=[-sin_lon[i], cos_lon[i], 0.]
for j=0, n_lat-1 do begin
    e_lat=[-sin_lat[j]*cos_lon[i], -sin_lat[j]*sin_lon[i], cos_lat[j]]
    e_r  =[ cos_lat[j]*cos_lon[i],  cos_lat[j]*sin_lon[i], sin_lat[j]]
for k=0, n_r-1 do begin
    rp_car=radius[k]*e_r
	Bp_car=fltarr(3)
	for s=0, 3 do begin
		dr=rp_car-ri_car[*,s]
		Bp_car=Bp_car+qi[s]*dr/(total(dr^2.)^1.5)
	endfor
    
	B_lon[i,j,k]=total(e_lon*Bp_car)
	B_lat[i,j,k]=total(e_lat*Bp_car)
	B_r  [i,j,k]=total(e_r  *Bp_car)
endfor
endfor
endfor

save, filename=bfile, B_lon, B_lat, B_r, lon_rad, lat_rad, radius

end
;------------------------------------------------------------
; this way gives the identical field as bfield_charge4_spherical, but is slower
pro bfield2_charge4_spherical, bfile, n_lon=n_lon, n_lat=n_lat, n_r=n_r

if ~keyword_set(n_lon) then n_lon=181
if ~keyword_set(n_lat) then n_lat=91
if ~keyword_set(n_r)   then n_r=31

ri=[[0.7*!pi, 0.0, 0.8], $
    [0.9*!pi, 0.1, 0.8], $
    [1.1*!pi,-0.1, 0.8], $
    [1.3*!pi, 0.1, 0.8]]
   
ri_car= convert_coordinate(ri, mode='lon_lat_r_to_xyz')
qi= [1.,-1.,1.,-1.]

lon_rad=findgen(n_lon)/(n_lon-1.)*2.*!pi
lat_rad=findgen(n_lat)/(n_lat-1.)*!pi-!pi/2.
radius =findgen(n_r)  /(n_r-1.)  *1.5+1.

seed_original=fltarr(3, n_lon, n_lat, n_r)
for i=0, n_lon-1 do seed_original[0,i,*,*]=lon_rad[i]
for i=0, n_lat-1 do seed_original[1,*,i,*]=lat_rad[i]
for i=0, n_r  -1 do seed_original[2,*,*,i]= radius[i]
seed_car= convert_coordinate(seed_original, mode='lon_lat_r_to_xyz')

B_car= fltarr(3, n_lon, n_lat, n_r)

for k=0, n_r  -1 do begin
for j=0, n_lat-1 do begin
for i=0, n_lon-1 do begin
	for s=0, 3 do begin
		dr=seed_car[*,i,j,k]-ri_car[*,s]
		B_car[*,i,j,k]=B_car[*,i,j,k]+qi[s]*dr/(total(dr^2)^1.5)
	endfor
endfor
endfor
endfor

seed_original= convert_coordinate(seed_car, B_car, v1=Bvec, mode='xyz_to_lon_lat_r')
B_lon=reform(Bvec[0,*,*,*])
B_lat=reform(Bvec[1,*,*,*])
B_r  =reform(Bvec[2,*,*,*])
save, filename=bfile, B_lon, B_lat, B_r, lon_rad, lat_rad, radius

end
;------------------------------------------------------------
; Examples for cartersian grid

;create a quadrupole field file
bfile='charge4_cartersian.sav'
if ~file_test(bfile) then bfield_charge4_cartersian, bfile
restore, bfile

nx=n_elements(xa)
ixh=(nx-1)/2
kend=n_elements(za)-1

;images of Figure 4 in Zhang, P., Chen, J.*, Liu, R. and Wang, C., 2022, ApJ, 937, 26
fastqsl, Bx, By, Bz, fname='method1_z0', /preview
fastqsl, Bx, By, Bz, fname='method2_z0', /scottFlag, /rk4, /preview
fastqsl, Bvec, xreg=[0,nx-1], yreg=[ixh,ixh], zreg=[0,kend/2], fname='method1_y0', /preview
fastqsl, Bvec, xreg=[0,nx-1], yreg=[ixh,ixh], zreg=[0,kend/2], fname='method2_y0', /scott, /preview

;An example of calculating in a cross section which is tilted to x-axis and y-axis, 
;and with stretched (actually uniformed) grids
fastqsl, Bx, By, Bz, xa=xa, ya=ya, za=za, delta=0.01, $
xreg=[-2,0], yreg=[1,0], zreg=[0,2], /csflag, $
fname='tilted_cs', /rk4, step=2.0, odir= 'fastqsl/', /twist, /preview

;An example of calculating in a box volume
fastqsl, Bx, By, Bz, xreg=[ixh/2,ixh], yreg=[ixh/4,ixh], zreg=[kend/4, kend/2], $
delta=0.8, tol=1.0e-3, odir= 'fastqsl', nthreads=12, /preview

; only exporting curlB
fastqsl, Bx, By, Bz, /curlB_out, maxsteps=0, seed='original', odir= 'fastqsl', qsl=qsl, /save_file

; Figure 4 of Chen (2026), see fname+'_logq_local.png'
fastqsl, Bx, By, Bz, xa=xa, ya=ya, za=za, $
xreg=[-1.7, 1.7], yreg=[0,0], zreg=[0,1.3], /preview, r_local=2, fname='r_local2'
fastqsl, Bx, By, Bz, xa=xa, ya=ya, za=za, $
xreg=[-1.7, 1.7], yreg=[0,0], zreg=[0,1.3], /preview, r_local=1, fname='r_local1'
fastqsl, Bx, By, Bz, xa=xa, ya=ya, za=za, $
xreg=[-1.7, 1.7], yreg=[0,0], zreg=[0,1.3], /preview, r_local=0.2, fname='r_local0.2'
fastqsl, Bx, By, Bz, xa=xa, ya=ya, za=za, $
xreg=[-1.7, 1.7], yreg=[0,0], zreg=[0,1.3], /preview, r_local=0.05, fname='r_local0.05'
;------------------------------------------------------------
; Examples for spherical grid
bfile='charge4_spherical.sav'
if ~file_test(bfile) then bfield_charge4_spherical, bfile
restore, bfile

; Q at r=1
fastqsl, b_lon, b_lat, b_r, xa=lon_rad, ya=lat_rad, za=radius, /spherical, $
factor=4, /preview, $
xreg=[0.,2*!Pi], yreg=[-!pi/2, !pi/2]


; trace field lines from two points
fastqsl, b_lon, b_lat, b_r, xa=lon_rad, ya=lat_rad, za=radius, /spherical, $
fname='spherical_seed_path', /preview, $
seed=[[!pi*0.85, 0.1, 1.], [!pi*1.1, -0.2, 1.2]], /path, /loopB, $ 
qsl=qsl

; transfrom *qsl.path[1] in (longitude, latitude, radius() to path1_car in (x, y, z), and
; transfrom *qsl.loopB[1] in (B_lon, B_lat, B_r) to loopB1_car in (B_x, B_y, B_z)
path1_car= convert_coordinate(*qsl.path[1], *qsl.loopB[1], mode='lon_lat_r_to_xyz', v1out=loopB1_car)

; wsa_par, b_lon, b_lat, b_r, lon_rad, lat_rad, radius, fs=fs, theta_b=theta_b, /preview, nth=12

end
