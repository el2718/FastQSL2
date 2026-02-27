module fields
implicit none
logical:: curlB_field_Flag, dbdc_field_Flag, stretchFlag, &
spherical, periodFlag(0:2), period_lon, south_pole, north_pole
integer:: binary_index_top, pend(0:2), dend(0:2)
integer, allocatable:: binary_values(:)
real:: pi, half_pi, two_pi, NaN, lat_pole, lat_pole2, pmin(0:2), pmax(0:2), period(0:2)
real, allocatable, target:: Bfield(:, :, :, :), CurlB_field(:, :, :, :), dbdc_field(:, :, :, :, :)
real, allocatable:: cos_lat_tmp(:)

type axis_stretch
	logical:: uni_Flag
	real:: da_uni
	real, allocatable :: pa(:), da(:), coef_diff(:,:)
	integer:: binary_index_start, index_try_start
endtype axis_stretch
type(axis_stretch), target :: axis(0:2)

type pole_field
	real, allocatable:: bfield(:,:,:,:), curlB_field(:,:,:,:), dbdc_field(:,:,:,:,:)
    real:: darc, over2darc, pmin(0:1), pmax(0:1), origin(0:1)
    integer:: aend
endtype pole_field
type(pole_field), target:: south, north

! interface can't use this type if site_info is defined in module trace
type site_info
	real:: v(0:8), dvds(0:8), B(0:2), CurlB(0:2), alpha, ds_factor, &
	v_yin(0:8), dvds_yin(0:8), B_yin(0:2), CurlB_yin(0:2)
	logical:: alphaFlag, yinFlag, scottFlag, scottLaunch
endtype site_info

procedure(), pointer:: round_weight

contains

subroutine readB
implicit none
logical:: preview, southFlag, B3flag
integer:: i, j, k, s, nx, ny, nz, nx_mag, ny_mag, aend1, round(0:1,0:2), j1, j2
real:: weight(0:1,0:1,0:1), mag_delta, clat_pole, dlast, dperiod, &
vp_yin(0:2), vp(0:2), bp_yin(0:2), bp(0:2), e_yin(0:2, 0:1), e_yang(0:2, 0:1), fj2
real, allocatable:: field_tmp(:, :, :, :), magnetogram(:, :), lon_tmp(:), lat_tmp(:)
real, pointer:: ax_tmp(:)
type(pole_field), pointer:: pole
!------------------------------------------------------------
! read Bx, By, Bz
open(1, file='bfield.bin', access='stream', status='old')
read(1) nx, ny, nz, B3flag, spherical, preview, periodFlag

pend = [nx, ny, nz] - 1
dend = pend - 1

if (B3flag) then
	allocate(field_tmp(0:pend(0), 0:pend(1), 0:pend(2), 0:2))
else
	allocate(field_tmp(0:2, 0:pend(0), 0:pend(1), 0:pend(2)))
endif
read(1) field_tmp
close(1, status='delete')
!------------------------------------------------------------
inquire(file='axis.bin', exist=stretchFlag)

if (stretchFlag) then

	binary_index_top=0

	open(1, file='axis.bin', access='stream', status='old')
	do s=0, 2
		allocate(axis(s)%pa(0:pend(s)))
		read(1) axis(s)%pa
	enddo
	close(1, status='delete')

	if (spherical) then ! period check
		dperiod= (axis(0)%pa(pend(0))-axis(0)%pa(0))-two_pi
		dlast = axis(0)%pa(pend(0))-axis(0)%pa(pend(0)-1)
		periodFlag(0)= dperiod .ge. -dlast*1.01 .and. dperiod .le. dlast*0.01
		periodFlag(1:2)=.false.

		if (periodFlag(0)) then
			period(0)=two_pi
			if (dperiod .lt. -dlast*0.01) then
				allocate(lon_tmp(0:pend(0)))
				lon_tmp=axis(0)%pa
				deallocate(axis(0)%pa)
				pend(0)=pend(0)+1
				dend(0)=pend(0)-1
				allocate(axis(0)%pa(0:pend(0)))
				axis(0)%pa(0:pend(0)-1)=lon_tmp
				deallocate(lon_tmp)				
			endif
			axis(0)%pa(pend(0))= axis(0)%pa(0) + period(0)
		endif
	endif
!------------------------------------------------------------
	do s=0, 2
		pmin(s)=axis(s)%pa(0)
		pmax(s)=axis(s)%pa(pend(s))
		allocate(axis(s)%da(0:dend(s)))
		axis(s)%da=axis(s)%pa(1:pend(s))-axis(s)%pa(0:dend(s))
		call diff_coefficent(s)

		! axis(s)%da(dend(s)) do not affect the value of round(0, i)
		axis(s)%uni_Flag= minval(abs(axis(s)%da(0:dend(s)-1))) &
			    .ge. 0.99*maxval(abs(axis(s)%da(0:dend(s)-1))) .or. (pend(s) .eq. 1)

		if (axis(s)%uni_Flag) then
			if (pend(s) .eq. 1) then
				axis(s)%da_uni=axis(s)%pa(1)-axis(s)%pa(0)
			else
				axis(s)%da_uni=(axis(s)%pa(dend(s))-axis(s)%pa(0))/dend(s)
			endif
		else
			axis(s)%binary_index_start=floor(dlog10(dble(dend(s)))/dlog10(2.0D0))
			if (axis(s)%binary_index_start .gt. binary_index_top) &
			binary_index_top=axis(s)%binary_index_start
		endif
	enddo
! stop
	if (binary_index_top .ne. 0) then
		allocate(binary_values(0:binary_index_top))
		forall(i=0:binary_index_top) binary_values(i)=2**i
		forall(s=0:2, .not. axis(s)%uni_Flag) &
		axis(s)%index_try_start=binary_values(axis(s)%binary_index_start)
	endif

	round_weight => round_weight_stretch
else
	pmin= 0.
	pmax= pend
	round_weight => round_weight_uniform
endif

period_lon = periodFlag(0) .and. spherical

if (.not. spherical) forall(s=0:2, periodFlag(s)) period(s)=pmax(i)-pmin(i)
!------------------------------------------------------------
allocate(Bfield(0:2, 0:pend(0), 0:pend(1), 0:pend(2)))

if (B3flag) then
! switch the indexes order for a better efficiency in subroutine interpolate
!$OMP PARALLEL DO PRIVATE(s), schedule(DYNAMIC)
	do s=0, 2
		Bfield(s,0:nx-1,0:ny-1,0:nz-1)=field_tmp(:,:,:,s)
	enddo
!$OMP END PARALLEL DO
else
	Bfield(:,0:nx-1,0:ny-1,0:nz-1)=field_tmp
endif
deallocate(field_tmp)
if (periodFlag(0)) Bfield(:,pend(0),:,:)=Bfield(:,0,:,:)
if (periodFlag(1)) Bfield(:,:,pend(1),:)=Bfield(:,:,0,:)
if (periodFlag(2)) Bfield(:,:,:,pend(2))=Bfield(:,:,:,0)
!------------------------------------------------------------
if (curlB_field_Flag) then

	allocate(curlB_field(0:2, 0:pend(0), 0:pend(1), 0:pend(2)))
	if (spherical) then
		allocate(cos_lat_tmp(0:pend(1)))
		cos_lat_tmp=cos(axis(1)%pa)
	endif
	!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC)
	do k=0, pend(2)
	do j=0, pend(1)
	do i=0, pend(0)
		call curlB_grid(i, j, k, curlB_field(:,i,j,k))
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
	if (spherical) deallocate(cos_lat_tmp)
	
endif
!------------------------------------------------------------
if (dbdc_field_Flag) then
	allocate(dbdc_field(0:2, 0:2, 0:pend(0), 0:pend(1), 0:pend(2)))
	!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC)
	do k=0, pend(2)
	do j=0, pend(1)
	do i=0, pend(0)
		call dbdc_grid(i, j, k, dbdc_field(:,:,i,j,k))
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
endif
!------------------------------------------------------------
! interpolate magnetogram on uniformed grids from the stretched input
if (stretchFlag .and. preview) then
	mag_delta=minval([axis(0)%da(0:dend(0)-1), axis(1)%da])

	nx_mag=nint((pmax(0)-pmin(0))/mag_delta)+1
	ny_mag=nint((pmax(1)-pmin(1))/mag_delta)+1
	
	allocate(magnetogram(0:nx_mag-1, 0:ny_mag-1))
	do j=0, ny_mag-1
	do i=0, nx_mag-1
		call round_weight(pmin + mag_delta * [i, j, 0], round, weight)
		magnetogram(i, j)=sum(weight(:, :, 0)*Bfield(2, round(:,0), round(:,1), 0))
	enddo
	enddo

	open(1, file='magnetogram.bin', access='stream', status='replace')
	write(1) nx_mag, ny_mag, mag_delta, magnetogram
	close(1)
	deallocate(magnetogram)
endif
!------------------------------------------------------------
south_pole = .false.
north_pole = .false.
if (.not. period_lon) return

! print*, period_lon
! return

lat_pole = 3./8. *Pi
! lat_pole = 31./64. *Pi

clat_pole= half_pi-lat_pole
lat_pole2= half_pi+0.01

do s=0, 1
	southFlag = s .eq. 0
	if (southFlag) then
		pole => south
		dlast= axis(1)%pa(1)-axis(1)%pa(0)
		south_pole = abs(pmin(1) + half_pi) .lt. 1.01*dlast
		if (.not. south_pole) cycle
		j1=1
	else
		pole => north
		dlast= axis(1)%pa(pend(1))-axis(1)%pa(pend(1)-1)
		north_pole = abs(pmax(1) - half_pi) .lt. 1.01*dlast
		if (.not. north_pole) cycle
		j1=pend(1)-1
	endif

	pole%darc= minval(axis(1)%da)
	pole%over2darc=0.5/pole%darc

	pole%aend=ceiling(clat_pole/pole%darc)+1
	aend1=pole%aend+1

	allocate(pole%Bfield(0:2, -aend1:aend1, -aend1:aend1, 0:pend(2)))
	allocate(lon_tmp(-aend1:aend1))
	allocate(lat_tmp(-aend1:aend1))

	forall(i= - aend1:aend1) lat_tmp(i)= i * pole%darc
	pole%origin(1)=0.
	pole%origin(0)=pi*(1.5-s)
	lon_tmp= lat_tmp+pole%origin(0)

	pole%pmax= pole%origin + pole%aend * pole%darc
	pole%pmin= pole%origin - pole%aend * pole%darc

	do k= 0, pend(2)
	do j= - aend1, aend1
	do i= - aend1, aend1
		if (abs(lat_tmp(j)) .le. half_pi-dlast) &
		call initializeB_yin([lon_tmp(i), lat_tmp(j), axis(2)%pa(k)], pole%Bfield(:,i,j,k))
	enddo
	enddo
	enddo

	do k= 0, pend(2)
		pole%Bfield(:,0,0,k)=0.
		do i= 0, pend(0)
			vp=[axis(0)%pa(i), axis(1)%pa(j1), axis(2)%pa(k)]
			bp=Bfield(:,i,j1,k)
			call vp_yinyang(vp_yin, vp, .false., e_yin, e_yang)
			forall(i=0:1) bp_yin(i)=bp(0)*dot_product(e_yang(:, 0), e_yin(:, i))+&
                                    bp(1)*dot_product(e_yang(:, 1), e_yin(:, i))     
			              bp_yin(2)=bp(2)
			pole%Bfield(:,0,0,k)= pole%Bfield(:,0,0,k)+ bp_yin
		enddo
		pole%Bfield(:,0,0,k)= pole%Bfield(:,0,0,k)/(pend(0)+1)
	enddo

	do k= 0, pend(2)
	do i= - aend1, aend1
		do j= 1, aend1
			if (lat_tmp(j) .le. half_pi-dlast) then
				j2=j
				fj2=float(j)
				exit
			endif
		enddo
		do j= -j2+1, -1
			pole%Bfield(:,i,j,k)= pole%Bfield(:,0,0,k)* (-j)/fj2+ pole%Bfield(:,0,-j2,k)* (j2+j)/fj2
		enddo
		do j= 1, j2-1
			pole%Bfield(:,i,j,k)= pole%Bfield(:,0,0,k)*   j /fj2+ pole%Bfield(:,0, j2,k)* (j2-j)/fj2
		enddo
	enddo
	enddo


	if (curlB_field_Flag) then
		allocate(pole%curlB_field(0:2, -pole%aend:pole%aend, -pole%aend:pole%aend, 0:pend(2)))
		allocate(cos_lat_tmp(-aend1:aend1))
		cos_lat_tmp=cos(lat_tmp)
		!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC)
		do k= 0, pend(2)
		do j= - pole%aend, pole%aend
		do i= - pole%aend, pole%aend
			call curlB_grid_pole(i, j, k, pole%curlB_field(:,i,j,k), southFlag)
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
		deallocate(cos_lat_tmp)
	endif

	if (dbdc_field_Flag) then
		allocate(pole%dbdc_field(0:2, 0:2, -pole%aend:pole%aend, -pole%aend:pole%aend, 0:pend(2)))
		!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC)
		do k= 0, pend(2)
		do j= - pole%aend, pole%aend
		do i= - pole%aend, pole%aend
			call dbdc_grid_pole(i, j, k, pole%dbdc_field(:,:,i,j,k), southFlag)
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
	endif
	deallocate(lon_tmp, lat_tmp)
enddo

! print*,south_pole,north_pole
END subroutine readB


subroutine diff_grid(i, j, k, s, index_diff, coef_diff, B3)
integer:: i, j, k, r, s, index_diff(0:2)
real:: coef_diff(0:2), B3(0:2,0:2)
!------------------------------------------------------------
if (s .eq. 0) r=i
if (s .eq. 1) r=j
if (s .eq. 2) r=k

if (periodFlag(s) .and. (r .eq. 0 .or. r .eq. pend(s))) then
	index_diff=[pend(s)-1, 0, 1]
	if (.not. stretchFlag) coef_diff=[-0.5, 0., 0.5]
else if (r .eq. 0) then
	index_diff=[0, 1, 2]
	if (.not. stretchFlag) coef_diff=[-1.5, 2.,-0.5]
else if (r .eq. pend(s)) then
	index_diff=[-2, -1, 0] + r
	if (.not. stretchFlag) coef_diff=[ 0.5,-2., 1.5]
else
	index_diff=[-1,  0, 1] + r
	if (.not. stretchFlag) coef_diff=[-0.5, 0., 0.5]
endif

if (stretchFlag) coef_diff = axis(s)%coef_diff(:, r)
	
if (s .eq. 0) B3=Bfield(:, index_diff, j, k)
if (s .eq. 1) B3=Bfield(:, i, index_diff, k)
if (s .eq. 2) B3=Bfield(:, i, j, index_diff)

end subroutine diff_grid

!dbdc, b:\vec(B)/B, c: coordinates
subroutine dbdc_grid(i, j, k, dbdc)
implicit none
integer:: i, j, k, r, s, t, index_diff(0:2)
real:: dbdc(0:2, 0:2), coef_diff(0:2), B3(0:2,0:2)
!------------------------------------------------------------
dbdc=0.
do s=0, 2
	call diff_grid(i, j, k, s, index_diff, coef_diff, B3)
	do t=0, 2
		if (coef_diff(t) .ne. 0.) &
		dbdc(s,:)=dbdc(s,:)+coef_diff(t)*normalize(B3(:,t))
	enddo
enddo
END subroutine dbdc_grid


subroutine curlB_grid(i, j, k, CurlBp)
implicit none
integer:: i, j, k, r, s, t, index_diff(0:2), index2(0:1)
real:: gradBp(0:2,0:2), CurlBp(0:2), coef_diff(0:2), B3(0:2,0:2), hs(0:2,0:2)
!------------------------------------------------------------
gradBp=0.

do s=0, 2
	call diff_grid(i, j, k, s, index_diff, coef_diff, B3)

	if (spherical) then
		if (s .eq. 0) hs(1:2,:) = 1./(axis(2)%pa(k)*cos_lat_tmp(j))
		if (s .eq. 1) then
			hs(0,:) = cos_lat_tmp(index_diff)/(axis(2)%pa(k)*cos_lat_tmp(j))
			hs(2,:) = 1./axis(2)%pa(k)
		endif
		if (s .eq. 2) then
			hs(0,:) = axis(2)%pa(index_diff)/axis(2)%pa(k)
			hs(2,:) = hs(0,:)
		endif
	endif

	index2=mod(s+[1,2],3)

	do t=0, 2
		if (coef_diff(t) .ne. 0.) then
			if (spherical) then
				gradBp(s,index2)=gradBp(s,index2)+coef_diff(t)*(B3(index2,t)*hs(index2,t))
			else
				gradBp(s,index2)=gradBp(s,index2)+coef_diff(t)*B3(index2,t)
			endif
		endif
	enddo	
enddo

curlBp=[gradBp(1,2)-gradBp(2,1), &
        gradBp(2,0)-gradBp(0,2), &
        gradBp(0,1)-gradBp(1,0)]

END subroutine curlB_grid


subroutine diff_coefficent(s)
implicit none
integer:: i, i_end, s
real, pointer:: da(:), coef(:,:)
!------------------------------------------------------------
i_end=pend(s)
allocate(axis(s)%coef_diff(0:2,0:i_end))
da   => axis(s)%da
coef => axis(s)%coef_diff

forall(i=1:i_end-1)
	coef(0,i)=-da(i)/(da(i-1)*(da(i)+da(i-1)))
	coef(2,i)= da(i-1)/(da(i)*(da(i)+da(i-1)))
end forall

if (periodFlag(s)) then
	coef(0,0)=-da(0)/(da(i_end-1)*(da(0)+da(i_end-1)))
	coef(2,0)= da(i_end-1)/(da(0)*(da(0)+da(i_end-1)))
	coef(0:2:2, i_end)=coef(0:2:2, 0)
else 
	coef(0,0)=-(2.0*da(0)+da(1))/(da(0)*(da(0)+da(1)))
	coef(2,0)=-da(0)/(da(1)*(da(0)+da(1)))

	coef(0,i_end)=  da(i_end-1)/(da(i_end-2)*(da(i_end-1)+da(i_end-2)))
	coef(2,i_end)= (2.0*da(i_end-1)+da(i_end-2))/(da(i_end-1)*(da(i_end-1)+da(i_end-2)))
endif

! this is a mathematical identity
coef(1, :) = - coef(0, :) - coef(2, :)

END subroutine diff_coefficent


function xy2lon(xy)
implicit none
real:: xy(0:1), xy2lon, cos_lon
!------------------------------------------------------------
if (all(xy .eq. 0.)) then 
	xy2lon=0.
else
	cos_lon=xy(0)/norm2s(xy)
	if (.not. (cos_lon .lt.  1.)) then
		xy2lon = 0.
	else if   (cos_lon .le. -1.)  then
		xy2lon = pi
	else if   (xy(1) .ge. 0.)     then
		xy2lon=          acos(cos_lon)
	else
		xy2lon= two_pi - acos(cos_lon)
	endif
endif
end function xy2lon


function vp_car2spherical(vp_car)
implicit none
real:: vp_car2spherical(0:2), vp_car(0:2)
!------------------------------------------------------------
vp_car2spherical(2)=norm2s(vp_car(0:2))
vp_car2spherical(1)=asin(vp_car(2)/vp_car2spherical(2))
vp_car2spherical(0)=xy2lon(vp_car(0:1))
end function vp_car2spherical


function vp_spherical2car(vp)
implicit none
real:: vp_spherical2car(0:2), vp(0:2), cos_p(0:1), sin_p(0:1)
!------------------------------------------------------------
cos_p=cos(vp(0:1))
sin_p=sin(vp(0:1))
vp_spherical2car = vp(2)* [cos_p(1)*[cos_p(0), sin_p(0)], sin_p(1)]
end function vp_spherical2car


function normalize(vector)
implicit none
real::vector(0:2), normalize(0:2)
normalize=vector/norm2s(vector)
end function normalize


! This approach is slightly faster than norm2()
function norm2s(vector)
implicit none
real::vector(:), norm2s !, vector(0:2)
norm2s=sqrt(sum(vector**2))
end function norm2s


function distance(vp, vp1)
implicit none
real:: distance, vp(0:2), vp1(0:2)
!------------------------------------------------------------
if (spherical) then
	distance= norm2s(vp_spherical2car(vp1)-vp_spherical2car(vp))
else
	distance= norm2s(vp1-vp)
endif
end function distance


function inside(vp)
implicit none
logical:: inside
real:: vp(0:2)
!------------------------------------------------------------
inside = &
(periodFlag(0) .or. (vp(0)>=pmin(0) .and. vp(0)<=pmax(0))) .and. &
((south_pole .or. periodFlag(1) .or.vp(1)>=pmin(1))        .and. &
 (north_pole .or. periodFlag(1) .or.vp(1)<=pmax(1)))       .and. &
(periodFlag(2) .or. (vp(2)>=pmin(2) .and. vp(2)<=pmax(2)))
end function inside


function inside_yin(vp_yin)
implicit none
type(pole_field), pointer::pole
real:: vp_yin(0:2)
logical:: inside_yin
!------------------------------------------------------------
if (vp_yin(0) .gt. pi) then 
	if (south_pole) then
 		pole => south
	else
		inside_yin= .false.
		return
	endif
else
	if (north_pole) then
		pole => north
	else
		inside_yin= .false.
		return
	endif
endif

inside_yin= all(pole%pmin < vp_yin(0:1) .and. pole%pmax > vp_yin(0:1)) 

end function inside_yin


function lie_boundary(vp)
implicit none
logical:: lie_boundary
real:: vp(0:2)
!---------------------------------------------------------------------------
lie_boundary = &
(.not. periodFlag(0) .and. vp(0)==pmin(0) .or. vp(0)==pmax(0))    .or. &
(.not. south_pole .and. .not. periodFlag(1) .and. vp(1)==pmin(1)) .or. &
(.not. north_pole .and. .not. periodFlag(1) .and. vp(1)==pmax(1)) .or. &
(.not. periodFlag(2) .and. vp(2)==pmin(2) .or. vp(2)==pmax(2))
end function lie_boundary


subroutine round_weight_uniform(vp, round, weight)
implicit none
real:: w(0:1,0:2), weight(0:1,0:1,0:1), vp(0:2), vpi
integer:: round(0:1,0:2), i, j, k
!------------------------------------------------------------
do i=0, 2
	vpi=vp(i)
	if (periodFlag(i)) vpi = modulo(vpi, period(i))

	if (.not. (vpi .gt. 0.0)) then
	! compared with vpi .le. 0.0, this way can prevent the crash from vpi=NaN (by B=0)
		round(0,i)=0
		w(1,i)=0.0
	else if (vpi .ge. pmax(i)) then
		round(0,i)=dend(i)
		w(1,i)=1.0
	else
		round(0,i)=floor(vpi)
		w(1,i)=vp(i)-round(0,i)
	endif
enddo

round(1,:)=round(0,:)+1
w(0,:)=1.0-w(1,:)
forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*w(k,2)

end subroutine round_weight_uniform


subroutine round_weight_stretch(vp, round, weight)
implicit none
real:: w(0:1,0:2), vp(0:2), vpi, weight(0:1,0:1,0:1)
integer:: i, j, k, round(0:1, 0:2), binary_index, index_try
!------------------------------------------------------------
do i=0, 2
	vpi=vp(i)
	if (periodFlag(i)) vpi = modulo(vpi-pmin(i), period(i)) + pmin(i)

	if (.not. (vpi .gt. pmin(i))) then
		round(0, i)=0
		w(1,i)=0.0
	else if   (vpi .ge. pmax(i))  then
		round(0, i)=dend(i)
		w(1,i)=1.0
	else
!------------------------------------------------------------
		if (axis(i)%uni_Flag) then
			round(0, i)=floor((vpi-pmin(i))/axis(i)%da_uni)
			if (round(0, i) .gt. dend(i)) round(0, i)=dend(i)
		else
!------------------------------------------------------------
			! this way is slower than binary (tree) search
			! round(0, i) = count(vpi .ge. axis(i)%pa(1:dend(i)))
!------------------------------------------------------------
			! binary (tree) search
			binary_index = axis(i)%binary_index_start
			index_try = axis(i)%index_try_start

			do while(binary_index .ge. 1)
				binary_index = binary_index-1
				if (vpi .ge. axis(i)%pa(index_try)) then
					if (index_try + binary_values(binary_index) .le. dend(i)) &
					index_try = index_try + binary_values(binary_index)
				else
					index_try = index_try - binary_values(binary_index)
				endif
			enddo
			
			if (vpi .ge. axis(i)%pa(index_try)) then
				round(0,i) = index_try
			else
				round(0,i) = index_try-1
			endif
	
		endif
!------------------------------------------------------------
		w(1,i) = (vpi-axis(i)%pa(round(0,i))) / axis(i)%da(round(0,i))
	endif
enddo

round(1,:)=round(0,:)+1
w(0,:)= 1.0 - w(1,:)
forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*w(k,2)

end subroutine round_weight_stretch


subroutine dbdc_grid_pole(i, j, k, dbdc, southflag)
implicit none
type(pole_field), pointer::pole
integer:: i, j, k, t, k_diff(0:2)
real:: dbdc(0:2, 0:2)
logical:: southflag
!------------------------------------------------------------
if (southflag) then
	pole => south
else
	pole => north
endif

if (k .eq. 0) then
	k_diff=[0,1,2]
else if (k .eq. pend(2)) then
	k_diff=[-2,-1,0]+k
else
	k_diff=[-1,0,1]+k
endif

dbdc=0.0

dbdc(0,:)= (normalize(pole%Bfield(:, i+1, j, k))-normalize(pole%Bfield(:, i-1, j, k)))*pole%over2darc
dbdc(1,:)= (normalize(pole%Bfield(:, i, j+1, k))-normalize(pole%Bfield(:, i, j-1, k)))*pole%over2darc

do t=0, 2
	if (axis(2)%coef_diff(t, k) .ne. 0.) &
	dbdc(2,:)=dbdc(2,:) + axis(2)%coef_diff(t, k)* normalize(pole%Bfield(:, i, j, k_diff(t)))
enddo

END subroutine dbdc_grid_pole


subroutine curlB_grid_pole(i, j, k, CurlBp, southflag)
implicit none
integer:: i, j, k, s, k_diff(0:2)
real:: gradBp(0:2,0:2), CurlBp(0:2), hs(0:2,0:2,0:2)
! hs(ib, coef_diff, diff)
logical:: southflag
type(pole_field), pointer::pole
!------------------------------------------------------------
if (k .eq. 0) then
	k_diff=[0,1,2]
else if (k .eq. pend(2)) then
	k_diff=[-2,-1,0]+k
else
	k_diff=[-1,0,1]+k
endif

if (southflag) then 
	pole => south
else
	pole => north
endif

hs(1:2,0:2:2,0)=1./(axis(2)%pa(k)*cos_lat_tmp(j))
hs(0,0:2:2,1)=cos_lat_tmp(j-1:j+1:2)*hs(1,0,0)
hs(2,0:2:2,1)=1./axis(2)%pa(k)

hs(:,0,0:1)=-hs(:,0,0:1) ! for central difference

hs(0,:,2)=axis(2)%pa(k_diff)/axis(2)%pa(k)
hs(1,:,2)=hs(0,:,2)

! if spherical, gradBp is not identical to \nabla \vec{B}, just for convenience of coding
forall(s=1:2)   gradBp(0,s)=sum(hs(s,0:2:2,0)*pole%Bfield(s, i-1:i+1:2, j, k))*pole%over2darc
forall(s=0:2:2) gradBp(1,s)=sum(hs(s,0:2:2,1)*pole%Bfield(s, i, j-1:j+1:2, k))*pole%over2darc
forall(s=0:1)   gradBp(2,s)=sum(axis(2)%coef_diff(:, k)*hs(s,:,2)*pole%Bfield(s, i, j, k_diff))

curlBp=[gradBp(1,2)-gradBp(2,1), &
        gradBp(2,0)-gradBp(0,2), &
        gradBp(0,1)-gradBp(1,0)]
END subroutine curlB_grid_pole


subroutine vp_yinyang(vp_yin, vp, toyang, e_yin, e_yang)
implicit none
real:: vp(0:2), vp_yin(0:2), cos_yin(0:1), sin_yin(0:1), cos_yang(0:1), sin_yang(0:1)
real, optional:: e_yin(0:2, 0:1), e_yang(0:2, 0:1)
logical:: toyang
!------------------------------------------------------------
if (toyang) then
	cos_yin=cos(vp_yin(0:1))
	sin_yin=sin(vp_yin(0:1))

	vp(0)=xy2lon([sin_yin(1),cos_yin(0)*cos_yin(1)])
	vp(1)=asin(cos_yin(1)*sin_yin(0))
	vp(2)=vp_yin(2)

	cos_yang=cos(vp(0:1))
	sin_yang=sin(vp(0:1))
else
	cos_yang=cos(vp(0:1))
	sin_yang=sin(vp(0:1))

	vp_yin(0)=xy2lon([cos_yang(1)*sin_yang(0), sin_yang(1)])
	vp_yin(1)=asin(cos_yang(1)*cos_yang(0))
	vp_yin(2)=vp(2)

	cos_yin=cos(vp_yin(0:1))
	sin_yin=sin(vp_yin(0:1))
endif

if (present(e_yin)) then
	e_yin(:, 0)=[0., -sin_yin(0), cos_yin(0)]
	e_yin(:, 1)=[cos_yin(1), -sin_yin(1)*[cos_yin(0), sin_yin(0)]]

	e_yang(:, 0)=[-sin_yang(0), cos_yang(0), 0.]
	e_yang(:, 1)=[-sin_yang(1)*[cos_yang(0), sin_yang(0)], cos_yang(1)]
endif
end subroutine vp_yinyang


subroutine initializeB_yin(vp_yin, bp_yin)
implicit none
integer:: i, round(0:1,0:2)
real:: vp(0:2), vp_yin(0:2), bp(0:2), bp_yin(0:2), weight(0:1,0:1,0:1), &
e_yang(0:2, 0:1), e_yin(0:2,0:1)
!------------------------------------------------------------
call vp_yinyang(vp_yin, vp, .true., e_yin, e_yang)
call round_weight(vp, round, weight)
forall(i=0:2) bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
forall(i=0:1) bp_yin(i)=bp(0)*dot_product(e_yang(:, 0), e_yin(:, i))+&
                        bp(1)*dot_product(e_yang(:, 1), e_yin(:, i))     
bp_yin(2)=bp(2)

end subroutine initializeB_yin


subroutine round_weight_pole(vp, round, weight, southflag)
implicit none
real:: w(0:1,0:2), vp(0:2), weight(0:1,0:1,0:1), p_lonlat(0:1)
integer:: i, j, k, round(0:1, 0:2), binary_index, index_try
logical:: southflag
!------------------------------------------------------------
if (southflag) then
	p_lonlat=(vp(0:1)-south%origin)/south%darc
else
	p_lonlat=(vp(0:1)-north%origin)/north%darc
endif
round(0, 0:1)=floor(p_lonlat)
w(1, 0:1)=p_lonlat-round(0, 0:1)

i=2
	if (.not. (vp(i) .gt. pmin(i))) then
		round(0, i)=0
		w(1,i)=0.0
	else if   (vp(i) .ge. pmax(i))  then
		round(0, i)=dend(i)
		w(1,i)=1.0
	else
!------------------------------------------------------------
		if (axis(i)%uni_Flag) then
			round(0, i)=floor((vp(i)-pmin(i))/axis(i)%da_uni)
			if (round(0, i) .gt. dend(i)) round(0, i)=dend(i)
		else
!------------------------------------------------------------
			! binary (tree) search
			binary_index = axis(i)%binary_index_start
			index_try = axis(i)%index_try_start

			do while(binary_index .ge. 1)
				binary_index = binary_index-1
				if (vp(i) .ge. axis(i)%pa(index_try)) then
					if (index_try + binary_values(binary_index) .le. dend(i)) &
					index_try = index_try + binary_values(binary_index)
				else
					index_try = index_try - binary_values(binary_index)
				endif
			enddo
			
			if (vp(i) .ge. axis(i)%pa(index_try)) then
				round(0,i) = index_try
			else
				round(0,i) = index_try-1
			endif

		endif
!------------------------------------------------------------
		w(1,i) = (vp(i)-axis(i)%pa(round(0,i))) / axis(i)%da(round(0,i))
	endif

round(1,:)=round(0,:)+1
w(0,:)= 1.0 - w(1,:)
forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*w(k,2)

end subroutine round_weight_pole

end module fields
