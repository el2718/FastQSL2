module fields
implicit none
logical:: CurlBvec_Flag, dbdc_field_Flag, stretchFlag, &
spherical, periodFlag(0:2), period_lon, south_pole, north_pole, A_input, &
keep_tmp, magnetogram_out ! thest two are inputed in fastqsl.f90
integer:: binary_index_top, pend(0:2), dend(0:2)
integer, allocatable:: binary_values(:)
real:: pi, half_pi, two_pi, NaN, lat_pole, lat_pole2, pmin(0:2), pmax(0:2), period(0:2)
real, allocatable, target:: Bvec(:,:,:,:), CurlBvec(:,:,:,:), &
dbdc_field(:,:,:,:,:), Avec(:,:,:,:)
real, allocatable:: cos_lat_tmp(:)

type axis_stretch
	logical:: uni_Flag
	real:: da_uni
	real, allocatable :: pa(:), da(:), coef_diff(:,:)
	integer:: binary_index_start, index_try_start
endtype axis_stretch
type(axis_stretch), target :: axis(0:2)

type pole_field
	real, allocatable:: Bvec(:,:,:,:), CurlBvec(:,:,:,:), dbdc_field(:,:,:,:,:), Avec(:,:,:,:)
    real:: darc, over2darc, pmin(0:1), pmax(0:1), origin(0:1)
    integer:: aend
endtype pole_field
type(pole_field), target:: south, north

! interface can't use this type if site_info is defined in module trace
type site_info
	real:: v(0:8), dvds(0:8), B(0:2), CurlB(0:2), A(0:2), ds_factor, &
	v_yin(0:8), dvds_yin(0:8), B_yin(0:2), CurlB_yin(0:2), A_yin(0:2), private(0:9)
	logical:: CurlBFlag, Aflag, yinFlag, scottFlag, scottLaunch
endtype site_info

procedure(), pointer:: round_weight

contains

subroutine readB
implicit none
logical:: southFlag, CurlB_input
integer:: i, j, k, s, t, nx, ny, nz, nx_mag, ny_mag, aend1, round(0:1,0:2), j1, j2
real:: weight(0:1,0:1,0:1), mag_delta, clat_pole, dlast, dperiod, &
e_yin(0:2, 0:1), e_yang(0:2, 0:1), fj2, vp_yin(0:2), vp(0:2), &
bp_yin(0:2), bp(0:2), ap_yin(0:2), ap(0:2), curlbp_yin(0:2), curlbp(0:2)
real, allocatable:: field_tmp(:,:,:,:), magnetogram(:,:), lon_tmp(:), lat_tmp(:)
real, pointer:: ax_tmp(:)
type(pole_field), pointer:: pole
!------------------------------------------------------------
! read Bx, By, Bz
open(1, file='field.bin', access='stream', status='old')
read(1) nx, ny, nz, stretchFlag, spherical, periodFlag, CurlB_input, A_input

pend = [nx, ny, nz] - 1
dend = pend - 1
!------------------------------------------------------------
if (stretchFlag) then

	binary_index_top=0
	do s=0, 2
		allocate(axis(s)%pa(0:pend(s)))
		read(1) axis(s)%pa
	enddo

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
allocate(field_tmp(0:2, 0:nx-1, 0:ny-1, 0:nz-1))
read(1) field_tmp
!------------------------------------------------------------
allocate(Bvec(0:2, 0:pend(0), 0:pend(1), 0:pend(2)))
Bvec(:, 0:nx-1, 0:ny-1, 0:nz-1)=field_tmp
if (periodFlag(0)) Bvec(:,pend(0),:,:)=Bvec(:,0,:,:)
if (periodFlag(1)) Bvec(:,:,pend(1),:)=Bvec(:,:,0,:)
if (periodFlag(2)) Bvec(:,:,:,pend(2))=Bvec(:,:,:,0)
!------------------------------------------------------------
if (CurlB_input) then
	read(1) field_tmp
	allocate(CurlBvec(0:2, 0:pend(0), 0:pend(1), 0:pend(2)))
	CurlBvec(:, 0:nx-1, 0:ny-1, 0:nz-1)=field_tmp
	if (periodFlag(0)) CurlBvec(:,pend(0),:,:)=CurlBvec(:,0,:,:)
	if (periodFlag(1)) CurlBvec(:,:,pend(1),:)=CurlBvec(:,:,0,:)
	if (periodFlag(2)) CurlBvec(:,:,:,pend(2))=CurlBvec(:,:,:,0)
endif
if (A_input) then
	read(1) field_tmp
	allocate(Avec(0:2, 0:pend(0), 0:pend(1), 0:pend(2)))
	Avec(:, 0:nx-1, 0:ny-1, 0:nz-1)=field_tmp
	if (periodFlag(0)) Avec(:,pend(0),:,:)=Avec(:,0,:,:)
	if (periodFlag(1)) Avec(:,:,pend(1),:)=Avec(:,:,0,:)
	if (periodFlag(2)) Avec(:,:,:,pend(2))=Avec(:,:,:,0)
endif
!------------------------------------------------------------
deallocate(field_tmp)
if (keep_tmp) then
	close(1, status='keep')
else
	close(1, status='delete')
endif
!------------------------------------------------------------
if (CurlBvec_Flag .and. .not. CurlB_input) then
	allocate(CurlBvec(0:2, 0:pend(0), 0:pend(1), 0:pend(2)))
	if (spherical) then
		allocate(cos_lat_tmp(0:pend(1)))
		cos_lat_tmp=cos(axis(1)%pa)
	endif
	!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC)
	do k=0, pend(2)
	do j=0, pend(1)
	do i=0, pend(0)
		call CurlB_grid(i, j, k, CurlBvec(:,i,j,k))
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
if (magnetogram_out) then

	if (stretchFlag) then
		mag_delta=minval([axis(0)%da(0:dend(0)-1), axis(1)%da])

		nx_mag=nint((pmax(0)-pmin(0))/mag_delta)+1
		ny_mag=nint((pmax(1)-pmin(1))/mag_delta)+1
		
		allocate(magnetogram(0:nx_mag-1, 0:ny_mag-1))
		do j=0, ny_mag-1
		do i=0, nx_mag-1
			call round_weight(pmin + mag_delta * [i, j, 0], round, weight)
			magnetogram(i, j)=sum(weight(:,:, 0)*Bvec(2, round(:,0), round(:,1), 0))
		enddo
		enddo
	else
		nx_mag=nx
		ny_mag=ny
		mag_delta=1.
		allocate(magnetogram(0:nx_mag-1, 0:ny_mag-1))
		magnetogram=Bvec(2, 0:nx_mag-1, 0:ny_mag-1, 0)
	endif

	open(1, file='magnetogram.bin', access='stream', status='replace')
	write(1) nx_mag, ny_mag, mag_delta, magnetogram
	close(1)
	deallocate(magnetogram)
endif
!------------------------------------------------------------
south_pole = .false.
north_pole = .false.
if (.not. period_lon) return

lat_pole = 3./8. *Pi

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

	allocate(pole%Bvec(0:2, -aend1:aend1, -aend1:aend1, 0:pend(2)))
	if (CurlB_input) allocate(pole%CurlBvec(0:2, -aend1:aend1, -aend1:aend1, 0:pend(2)))
	if (A_input) allocate(pole%Avec(0:2, -aend1:aend1, -aend1:aend1, 0:pend(2)))
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
		if (abs(lat_tmp(j)) .le. half_pi-dlast) then
			vp_yin=[lon_tmp(i), lat_tmp(j), axis(2)%pa(k)]
			call vp_yinyang(vp_yin, vp, .true., e_yin, e_yang)
			call round_weight(vp, round, weight)
			forall(t=0:2) bp(t)=sum(weight*Bvec(t, round(:,0), round(:,1), round(:,2)))
			forall(t=0:1) bp_yin(t)=bp(0)*dot_product(e_yang(:, 0), e_yin(:, t))+&
									bp(1)*dot_product(e_yang(:, 1), e_yin(:, t))     
			bp_yin(2)=bp(2)
			pole%Bvec(:,i,j,k)=bp_yin
			if (CurlB_input) then
				forall(t=0:2) CurlBp(t)=sum(weight*CurlBvec(t, round(:,0), round(:,1), round(:,2)))
				forall(t=0:1) CurlBp_yin(t)=CurlBp(0)*dot_product(e_yang(:, 0), e_yin(:, t))+&
										    CurlBp(1)*dot_product(e_yang(:, 1), e_yin(:, t))     
				CurlBp_yin(2)=CurlBp(2)
				pole%CurlBvec(:,i,j,k)=CurlBp_yin
			endif
			if (A_input) then
				forall(t=0:2) ap(t)=sum(weight*Avec(t, round(:,0), round(:,1), round(:,2)))
				forall(t=0:1) ap_yin(t)=ap(0)*dot_product(e_yang(:, 0), e_yin(:, t))+&
										ap(1)*dot_product(e_yang(:, 1), e_yin(:, t))     
				ap_yin(2)=ap(2)
				pole%Avec(:,i,j,k)=ap_yin
			endif
		endif
	enddo
	enddo
	enddo

	do k= 0, pend(2)
		pole%Bvec(:,0,0,k)=0.
		if (CurlB_input) pole%CurlBvec(:,0,0,k)=0.
		if (A_input) pole%Avec(:,0,0,k)=0.
		do i= 0, pend(0)
			vp=[axis(0)%pa(i), axis(1)%pa(j1), axis(2)%pa(k)]
			call vp_yinyang(vp_yin, vp, .false., e_yin, e_yang)
			bp=Bvec(:,i,j1,k)
			forall(t=0:1) bp_yin(t)=bp(0)*dot_product(e_yang(:, 0), e_yin(:, t))+&
			                        bp(1)*dot_product(e_yang(:, 1), e_yin(:, t))     
			              bp_yin(2)=bp(2)
			pole%Bvec(:,0,0,k)= pole%Bvec(:,0,0,k)+ bp_yin
			if (CurlB_input) then 
				CurlBp=CurlBvec(:,i,j1,k)
				forall(t=0:1) CurlBp_yin(t)=CurlBp(0)*dot_product(e_yang(:, 0), e_yin(:, t))+&
										    CurlBp(1)*dot_product(e_yang(:, 1), e_yin(:, t))     
							  CurlBp_yin(2)=CurlBp(2)
				pole%CurlBvec(:,0,0,k)= pole%CurlBvec(:,0,0,k)+ CurlBp_yin
			endif
			if (A_input) then 
				ap=Avec(:,i,j1,k)
				forall(t=0:1) ap_yin(t)=ap(0)*dot_product(e_yang(:, 0), e_yin(:, t))+&
										ap(1)*dot_product(e_yang(:, 1), e_yin(:, t))     
							  ap_yin(2)=ap(2)
				pole%Avec(:,0,0,k)= pole%Avec(:,0,0,k)+ ap_yin
			endif
		enddo
		pole%Bvec(:,0,0,k)= pole%Bvec(:,0,0,k)/(pend(0)+1)
		if (CurlB_input) pole%CurlBvec(:,0,0,k)= pole%CurlBvec(:,0,0,k)/(pend(0)+1)
		if (A_input) pole%Avec(:,0,0,k)= pole%Avec(:,0,0,k)/(pend(0)+1)
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
			pole%Bvec(:,i,j,k)= pole%Bvec(:,0,0,k)* (-j)/fj2+ pole%Bvec(:,0,-j2,k)* (j2+j)/fj2
			if (CurlB_input) pole%CurlBvec(:,i,j,k)= pole%CurlBvec(:,0,0,k)* (-j)/fj2+ pole%CurlBvec(:,0,-j2,k)* (j2+j)/fj2
			if (A_input) pole%Avec(:,i,j,k)= pole%Avec(:,0,0,k)* (-j)/fj2+ pole%Avec(:,0,-j2,k)* (j2+j)/fj2
		enddo
		do j= 1, j2-1
			pole%Bvec(:,i,j,k)= pole%Bvec(:,0,0,k)*   j /fj2+ pole%Bvec(:,0, j2,k)* (j2-j)/fj2
			if (CurlB_input) pole%CurlBvec(:,i,j,k)= pole%CurlBvec(:,0,0,k)*   j /fj2+ pole%CurlBvec(:,0, j2,k)* (j2-j)/fj2
			if (A_input) pole%Avec(:,i,j,k)= pole%Avec(:,0,0,k)*   j /fj2+ pole%Avec(:,0, j2,k)* (j2-j)/fj2
		enddo
	enddo
	enddo

	if (CurlBvec_Flag .and. .not. CurlB_input) then
		allocate(pole%CurlBvec(0:2, -pole%aend:pole%aend, -pole%aend:pole%aend, 0:pend(2)))
		allocate(cos_lat_tmp(-aend1:aend1))
		cos_lat_tmp=cos(lat_tmp)
		!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC)
		do k= 0, pend(2)
		do j= - pole%aend, pole%aend
		do i= - pole%aend, pole%aend
			call CurlB_grid_pole(i, j, k, pole%CurlBvec(:,i,j,k), southFlag)
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
	
if (s .eq. 0) B3=Bvec(:, index_diff, j, k)
if (s .eq. 1) B3=Bvec(:, i, index_diff, k)
if (s .eq. 2) B3=Bvec(:, i, j, index_diff)

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


subroutine CurlB_grid(i, j, k, CurlBp)
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

CurlBp=[gradBp(1,2)-gradBp(2,1), &
        gradBp(2,0)-gradBp(0,2), &
        gradBp(0,1)-gradBp(1,0)]

END subroutine CurlB_grid


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
coef(1,:) = - coef(0,:) - coef(2,:)

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
vp_car2spherical(2)=norm2s(vp_car)
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

dbdc(0,:)= (normalize(pole%Bvec(:, i+1, j, k))-normalize(pole%Bvec(:, i-1, j, k)))*pole%over2darc
dbdc(1,:)= (normalize(pole%Bvec(:, i, j+1, k))-normalize(pole%Bvec(:, i, j-1, k)))*pole%over2darc

do t=0, 2
	if (axis(2)%coef_diff(t, k) .ne. 0.) &
	dbdc(2,:)=dbdc(2,:) + axis(2)%coef_diff(t, k)* normalize(pole%Bvec(:, i, j, k_diff(t)))
enddo

END subroutine dbdc_grid_pole


subroutine CurlB_grid_pole(i, j, k, CurlBp, southflag)
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
forall(s=1:2)   gradBp(0,s)=sum(hs(s,0:2:2,0)*pole%Bvec(s, i-1:i+1:2, j, k))*pole%over2darc
forall(s=0:2:2) gradBp(1,s)=sum(hs(s,0:2:2,1)*pole%Bvec(s, i, j-1:j+1:2, k))*pole%over2darc
forall(s=0:1)   gradBp(2,s)=sum(axis(2)%coef_diff(:, k)*hs(s,:,2)*pole%Bvec(s, i, j, k_diff))

CurlBp=[gradBp(1,2)-gradBp(2,1), &
        gradBp(2,0)-gradBp(0,2), &
        gradBp(0,1)-gradBp(1,0)]
END subroutine CurlB_grid_pole


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


module trace
use fields
implicit none
integer:: maxsteps, maxsteps_foot, normal_index
logical:: RK4flag, diff_flag, inclineFlag, traceflag, int_private_out(0:9), privateFlag, &
B_out, CurlB_out, rf_out, targetB_out, targetCurlB_out, loopCurlB_out
real:: a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, &
b1, b3, b4, b5, b6, ce1, ce3, ce4, ce5, ce6, &                   ! for RKF45
step, min_step, max_step, min_step_foot, tol, min_incline, r_local, r_local_square
! r: remote; b:boundary, F:foot; s/e: start/end, where B vector points to inside/outside
! rbs/rbe: boundary mark for start/end point of the field line, see subroutine trim_size
! rFs/rFe: coordinates of the start/end point of the field line
! bs/bp/be: B at rFs/vp/rFe
! CurlBs/CurlBp/CurlBe: CurlB at rFs/vp/rFe
! bnp, normal component of B to the surface
! tangent, the angle between B and surface .lt. arctan(min_incline)
! q/q_perp: value of q/q_perp given by Scott (2017)
! label: label of exporting fieldlines
type line_info
	real:: rFs(0:2), rFe(0:2), ev3(0:2), bp(0:2), bnp, &
	bs(0:2), be(0:2), CurlBp(0:2), CurlBs(0:2), CurlBe(0:2), &
	q, q_perp, q_local, brn_s, brn_e, rFs_yin(0:2), rFe_yin(0:2), int_private(0:9)
	logical:: tangent, get, scottFlag, path_out, loopB_out, loopCurlB_out, &
	q_local_Flag, local_s_flag, local_e_flag, s_yinflag, e_yinflag
	integer:: rbs, rbe, its, ite
	real, allocatable:: path(:,:), loopB(:,:), loopCurlB(:,:)
endtype line_info

interface
	subroutine interpolate_site(site, rk_first)
	use fields
	type(site_info), target :: site
	logical, optional:: rk_first
	end subroutine interpolate_site
end interface

procedure(interpolate_site), pointer:: interpolate

contains

include 'privates.f90'

subroutine cal_yinyang(site, toyang, dvdsflag)
implicit none
type(site_info), target :: site
real, target:: e_yang(0:2, 0:1), e_yin(0:2, 0:1)
real, pointer:: vector1(:), vector2(:), e1(:,:), e2(:,:), dvds1(:), dvds2(:)
real:: matrix0(0:1, 0:1), matrix1(0:1, 0:1), matrix2(0:1, 0:1), &
matrix3(0:1, 0:1), matrix4(0:1, 0:1)
integer:: i, j
logical:: toyang
logical, optional:: dvdsflag
!------------------------------------------------------------
call vp_yinyang(site%v_yin(0:2), site%v(0:2), toyang, e_yin, e_yang)

if (site%scottFlag .or. present(dvdsflag)) then

if (toyang) then
	vector1 =>site%v_yin
	vector2 =>site%v
	e1      =>e_yin
	e2      =>e_yang
else
	vector1 =>site%v
	vector2 =>site%v_yin
	e1      =>e_yang
	e2      =>e_yin
endif

forall(i=0:1, j=0:1) matrix0(i,j)=dot_product(e1(:, i), e2(:, j))

if (site%scottFlag) then
	vector2(5:8:3)=vector1(5:8:3)
	forall(i=0:1, j=1:2) vector2(i+j*3)=dot_product(vector1(j*3:j*3+1), matrix0(:, i))
endif

if (present(dvdsflag)) then
	if (dvdsflag) then
		if (toyang) then
			dvds1 => site%dvds_yin
			dvds2 => site%dvds
		else
			dvds1 => site%dvds
			dvds2 => site%dvds_yin
		endif

		!\partial (vector2(0:1))/\partial (vector1(0:1))
		matrix1(0,:)=matrix0(0,:)*cos(vector1(1))
		matrix1(:,0)=matrix1(:,0)/cos(vector2(1))
		
		forall(i=0:1) dvds2(i)=dot_product(dvds1(0:1), matrix1(:, i))
		dvds2(2)=dvds1(2)

		if (site%scottFlag) then
			matrix2(:,0)= matrix0(:,1)
			matrix2(:,1)=-matrix0(:,0)

			matrix3(0,:)= matrix0(1,:)
			matrix3(1,:)=-matrix0(0,:)

			!d martrix0/ds
			matrix4= dvds2(0)*sin(vector2(1))*matrix2 + dvds1(0)*sin(vector1(1))*matrix3

			forall(i=0:1, j=1:2) & 
			dvds2(i+j*3)=dot_product(  dvds1(j*3:j*3+1), matrix0(:,i))+&
						 dot_product(vector1(j*3:j*3+1), matrix4(:,i))
			dvds2(5:8:3)=dvds1(5:8:3)
		endif
	else ! if dvdsflag is .false., toyang is .true. already
		forall(i=0:1) site%B(i)=dot_product(site%B_yin(0:1), matrix0(:, i))
		site%B(2)=site%B_yin(2)
		if (site%CurlBFlag) then
			forall(i=0:1) site%CurlB(i)=dot_product(site%CurlB_yin(0:1), matrix0(:, i))
			site%CurlB(2)=site%CurlB_yin(2)
		endif
		if (site%Aflag) then
			forall(i=0:1) site%A(i)=dot_product(site%A_yin(0:1), matrix0(:, i))
			site%A(2)=site%A_yin(2)
		endif
	endif
endif
endif
end subroutine cal_yinyang


subroutine RK4(dt, site, site1)
implicit none
real:: dt, ds, k2(0:8), k3(0:8)
real, pointer:: vector(:), vector1(:), dvds1(:), k1(:)
type(site_info), target:: site, site1
!------------------------------------------------------------
! the unit of ds is same as the physical unit (e.g. the unit of xreg, yreg, zreg)
! the unit of dt is the scale of the cell (a self-adaptive fashion that varying from cell to cell)
ds=dt/site%ds_factor

if (site%yinflag) then
	vector  => site %v_yin
	k1      => site %dvds_yin
	vector1 => site1%v_yin
	dvds1   => site1%dvds_yin
else
	vector  => site %v
	k1      => site %dvds
	vector1 => site1%v
	dvds1   => site1%dvds
endif

site1%yinflag   =site%yinflag

vector1 = vector + ds*  1./3.*k1
call interpolate(site1)
k2 = dvds1

vector1 = vector + ds*(-1./3.*k1+k2)
call interpolate(site1)
k3 = dvds1

vector1 = vector + ds*(k1-k2+k3)
call interpolate(site1)

vector1 = vector + ds*0.125*(k1+3.*(k2+k3)+dvds1)

if (site%yinflag) call cal_yinyang(site1, .true.)

end subroutine RK4


subroutine RKF45(dt, site, site1, tol_this, dt_executed)
implicit none
logical:: repeat_flag
integer:: rb
real:: k2(0:8), k3(0:8), k4(0:8), k5(0:8), dvp(0:2), &
dt, dt_executed, ds, ds0, trim_factor, error, tol_this, tol_ds, incline
real, pointer:: vector(:), vector1(:), dvds1(:), k1(:)
type(site_info), target:: site, site1
!------------------------------------------------------------
tol_ds=tol_this/site%ds_factor

if (site%yinflag) then
	vector  => site %v_yin
	k1      => site %dvds_yin
	vector1 => site1%v_yin
	dvds1   => site1%dvds_yin
else
	vector  => site %v
	k1      => site %dvds
	vector1 => site1%v
	dvds1   => site1%dvds
endif

site1%yinflag   =site%yinflag

repeat_flag=.true.

do while (repeat_flag)
	ds=dt/site%ds_factor
	
	vector1 = vector + ds* a21*k1
	call interpolate(site1)
	k2 = dvds1

	vector1 = vector + ds*(a31*k1+ a32*k2)
	call interpolate(site1)
	k3 = dvds1

	vector1 = vector + ds*(a41*k1+ a42*k2+ a43*k3)
	call interpolate(site1)
	k4 = dvds1

	vector1 = vector + ds*(a51*k1+ a52*k2+ a53*k3+ a54*k4)
	call interpolate(site1)
	k5 = dvds1

	vector1 = vector + ds*(a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5)
	call interpolate(site1)

	vector1 = vector + ds*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*dvds1)

	dvp = ds*(ce1*k1(0:2)+ce3*k3(0:2)+ce4*k4(0:2)+ce5*k5(0:2)+ce6*dvds1(0:2))

	if (spherical) then
		! error = distance(vp1-dvp, vp1) is mathematically correct. 
		! While vp1 and dvp are float arrays, and dvp is always very small, 
		! results distance(vp1-dvp, vp1) eq 0. for most cases. 
		error = norm2s([dvp(0)*vector(2)*cos(vector(1)), dvp(1)*vector(2), dvp(2)])
	else
		error = norm2s(dvp)
	endif
	
	dt_executed = dt
	
	if (site%yinFlag) call cal_yinyang(site1, .true.)
!------------------------------------------------------------
	if (inside(site1%v(0:2))) then
		repeat_flag = (error .gt. tol_ds) .and. (abs(dt_executed) .gt. min_step)
		if (error .gt. 0.) then
			dt= dt_executed * ((tol_ds/error)**0.2)*0.9
			if (.not. (abs(dt) .ge. min_step)) dt=sign(min_step, dt_executed)
			if (abs(dt) .gt. max_step)         dt=sign(max_step, dt_executed)
		else
			dt= sign(max_step, dt_executed)
		endif
	else
		repeat_flag = abs(dt_executed) .gt. min_step
		if (repeat_flag) then
			call trim_size(site, site1, rb, ds0, trim_factor, incline)
	
			! then if a next do loop exist, repeat_flag mostly will be .false. in that loop
			! because dt will \approx 0.5*sign(min_step, dt) in that loop
			dt = sign(abs(dt_executed*trim_factor)-0.5*min_step, dt_executed)

			if (.not. (abs(dt) .ge. min_step) .or. (incline .eq. min_incline)) dt = sign(min_step, dt_executed)
		endif
	endif
enddo

end subroutine RKF45


subroutine correct_foot(launch_step, dt_executed, site0, site1, rb, identical)
implicit none
logical:: launch_step, identical, flagRK4B, repeat_flag
integer:: rb, rb_index, it
real:: dt, dt_executed, ds0, trim_factor, incline, &
vector1_orig(0:8), v1_yang_orig(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8)
real, pointer:: dvds(:), dvds0(:), dvds1(:), vector(:), vector0(:), vector1(:)
type(site_info), target:: site0, site, site1
!------------------------------------------------------------
call trim_size(site0, site1, rb, ds0, trim_factor, incline)

! if just take one step to be outside in trace_bline, retrace by rk4 with min_step_foot here.
! this will improve the quality of Q around PIL
identical = lie_boundary(site0%v(0:2)) .and. (.not. launch_step)
if (rb .eq. 0 .or. rb .eq. 8 .or. identical) return
!------------------------------------------------------------
if (site0%yinflag) then
	vector0 => site0%v_yin
	dvds0   => site0%dvds_yin
	vector  => site %v_yin
	dvds    => site %dvds_yin
	vector1 => site1%v_yin
	dvds1   => site1%dvds_yin
	v1_yang_orig=site1%v
else
	vector0 => site0%v
	dvds0   => site0%dvds
	vector  => site %v
	dvds    => site %dvds
	vector1 => site1%v
	dvds1   => site1%dvds
endif
vector1_orig=vector1
!------------------------------------------------------------
! the output vector1 is expected to be inside, 
! the residual distance to the boundary \approx 0.5*min_step_foot
dt=sign(abs(dt_executed*trim_factor)-0.5*min_step_foot, dt_executed)
if (abs(dt) .ge. min_step_foot) call RK4(dt, site0, site1)
!------------------------------------------------------------
site%yinFlag  =site0%yinFlag
site%scottFlag=site0%scottFlag
site%ds_factor=site0%ds_factor
if (inside(site1%v(0:2))) then
	if (lie_boundary(site1%v(0:2))) return
	vector = vector1
	call interpolate(site)
else
	vector = vector0
	dvds   = dvds0
endif

dt=sign(min_step_foot, dt_executed)
it=0
repeat_flag=.true.
do while (repeat_flag)
	if (it .ne. 0) then
		vector = vector1
		call interpolate(site)
	endif
	call RK4(dt, site, site1)
	it=it+1
	repeat_flag = inside(site1%v(0:2)) .and. (it .lt. maxsteps_foot)
end do
!------------------------------------------------------------
if (inside(site1%v(0:2))) then
	if (lie_boundary(site1%v(0:2))) return
	vector  = vector0
	dvds    = dvds0
	vector1 = vector1_orig
	if (site0%yinFlag) site1%v = v1_yang_orig
else
	if (site0%yinFlag) call cal_yinyang(site, .true.)
	if (lie_boundary(site%v(0:2))) then
		vector1 = vector
		if (site0%yinFlag) site1%v = site%v
		identical = all(vector0(0:2) .eq. vector(0:2))
		return
	endif
	dt_executed=dt
	call trim_size(site, site1, rb, ds0, trim_factor, incline)
endif

if (rb .ge. 7) return
!------------------------------------------------------------
flagRK4B = incline .gt. min_incline
rb_index=(6-rb)/2
if (flagRK4B) then
	! d vp/d r_i=Bp/B_i

	k1 = dvds/dvds(rb_index)

	vector1 = vector + ds0*  1./3.*k1
	call interpolate(site1)
	k2 = dvds1/dvds1(rb_index)

	vector1 = vector + ds0*(-1./3.*k1+k2)
	call interpolate(site1)
	k3 = dvds1/dvds1(rb_index)

	vector1 = vector + ds0*(k1-k2+k3)
	call interpolate(site1)
	k4 = dvds1/dvds1(rb_index)

	flagRK4B= all(abs([k2(0:2),k3(0:2),k4(0:2)]) .lt. 20.)
endif

if (flagRK4B) then
	vector1 = vector + ds0*0.125*(k1+3.*(k2+k3)+k4)
else
	vector1 = vector*(1.-trim_factor) + vector1*trim_factor
endif

if (site0%yinFlag) call cal_yinyang(site1, .true.)

! avoid the potential error caused by machine precision
if (mod(rb, 2) .eq. 1) then
	site1%v(rb_index)=pmin(rb_index)
else
	site1%v(rb_index)=pmax(rb_index)
endif

end subroutine correct_foot


subroutine trim_size(site, site1, rb, ds0, trim_factor, incline)
implicit none
logical:: boundary_mark(1:6)
integer:: i, rb, rb_index
real:: vp_mid(0:2), vp(0:2), vp1(0:2), ds0, trim_factor, incline
type(site_info), target:: site, site1
!------------------------------------------------------------
vp =site %v(0:2)
vp1=site1%v(0:2)

boundary_mark(1)= (.not. periodFlag(2)) .and. .not. (vp1(2) .ge. pmin(2))
boundary_mark(2)= (.not. periodFlag(2)) .and. .not. (vp1(2) .le. pmax(2))
boundary_mark(3)= (.not. periodFlag(1)) .and. .not. (vp1(1) .ge. pmin(1)) .and. .not. south_pole
boundary_mark(4)= (.not. periodFlag(1)) .and. .not. (vp1(1) .le. pmax(1)) .and. .not. north_pole
boundary_mark(5)= (.not. periodFlag(0)) .and. .not. (vp1(0) .ge. pmin(0))
boundary_mark(6)= (.not. periodFlag(0)) .and. .not. (vp1(0) .le. pmax(0))
!------------------------------------------------------------
incline = 1. ! give a return if rb .eq. 0, 7, 8

if      (count(boundary_mark) .eq. 0) then
	rb=0  !inside
else if (count(boundary_mark) .eq. 1) then
	! rb=findloc(boundary_mark, .true. , 1)
	! findloc is supported from gfortran 9.0, while many servers use gfortran 4.9
	do i = 1, 6
		if (boundary_mark(i)) then
			rb=i
			exit
		endif
	enddo
	
	rb_index=(6-rb)/2

	if (mod(rb, 2) .eq. 1) then
		ds0= pmin(rb_index) - vp(rb_index)
	else
		ds0= pmax(rb_index) - vp(rb_index)
	endif

	trim_factor = ds0 / (vp1(rb_index) - vp(rb_index))

	if (spherical) then
		if (site%yinflag) then
			! rb_index can only be 2 in this case
			incline= abs(site%b_yin(2)/norm2s(site%b_yin))
		else
			incline= abs(site%b(rb_index)/norm2s(site%b))
		endif
	else
		incline=abs(site%dvds(rb_index))
	endif

	if (incline .lt. min_incline) incline = min_incline

else if (.not. (all(boundary_mark(1:2)) .or. all(boundary_mark(3:4)) .or. all(boundary_mark(5:6)))) then
	rb=7  ! end at edges or corners

	vp_mid=(vp1+vp)*0.5
	trim_factor=0.5
	do while (.not. inside(vp_mid) .and. distance(vp_mid, vp) .gt. min_step_foot)
		trim_factor=trim_factor*0.5
		vp_mid=vp1*trim_factor + (1.-trim_factor)*vp
	enddo
else ! NaN found in vp1, caused by terminated at where B is 0/NaN
	rb=8
	trim_factor=0. ! to make dt = sign(min_step, dt_executed) in subroutine rkf45
endif

end subroutine trim_size


function normalize_cross_product(v1, v2)
implicit none
real:: v1(0:2), v2(0:2), normalize_cross_product(0:2)
!------------------------------------------------------------
normalize_cross_product= &
normalize([v1(1)*v2(2)-v1(2)*v2(1), &
           v1(2)*v2(0)-v1(0)*v2(2), &
	       v1(0)*v2(1)-v1(1)*v2(0)])
end function normalize_cross_product


subroutine set_u_v(bp, u0, v0)
implicit none
integer:: maxdim, index1
real:: bp(0:2), u0(0:2), v0(0:2)
!------------------------------------------------------------
maxdim=maxloc(abs(bp), dim=1)-1
index1=mod(maxdim+1,3)

u0(maxdim)= bp(index1)
u0(index1)=-bp(maxdim)
u0(mod(maxdim+2,3))= 0.
u0=normalize(u0)

v0 = normalize_cross_product(bp, u0)

end subroutine set_u_v


subroutine interpolate_spherical(site, rk_first)
implicit none
type(site_info), target:: site
type(pole_field), pointer:: pole
logical, optional:: rk_first
logical:: southflag, yinflag
integer:: round(0:1,0:2), i, j, k
real:: weight(0:1,0:1,0:1), r, sin_lat, cos_lat, &
dbdc_cell(0:2,0:2,0:1,0:1,0:1), dbdcp(0:2,0:2), da(0:2), Ap(0:2)
real, pointer:: bp(:), vector(:), dvds(:), CurlBp(:)
!------------------------------------------------------------
if (present(rk_first)) then
	yinflag = (south_pole .and. (-site%v(1) .gt. lat_pole) .and. (-site%v(1) .le. lat_pole2)) &
	    .or.  (north_pole .and. ( site%v(1) .gt. lat_pole) .and. ( site%v(1) .le. lat_pole2))

	! a RK step or correct_foot will given site%v finally,
	! do not need to process the case of (.not. yinflag .and. site%yinflag) 
	if (yinflag .and. .not. site%yinflag) call cal_yinyang(site, .false.)

	! site%yinFlag is inherited from the last RK step, for the next step, it should be updated
	site%yinFlag = yinflag
else
	if (site%yinflag) then
		if (inside_yin(site%v_yin(0:2))) then
			yinflag= .true.
		else
			call cal_yinyang(site, .true.)
			yinflag= .false.
		endif
	else
		yinflag= .false.
	endif
endif
!------------------------------------------------------------
if (yinflag) then
	southflag = site%v_yin(0) .gt. pi
	if (southflag) then 
		pole => south
	else
		pole => north
	endif
	vector => site%v_yin
	dvds   => site%dvds_yin
	bp     => site%b_yin
	call round_weight_pole(site%v_yin(0:2), round, weight, southflag)
	forall(i=0:2) Bp(i)=sum(weight* pole%Bvec(i, round(:,0), round(:,1), round(:,2)))
else
	vector => site%v
	dvds   => site%dvds
	bp     => site%b
	call round_weight(site%v(0:2), round, weight)
	forall(i=0:2) Bp(i)=sum(weight* Bvec(i, round(:,0), round(:,1), round(:,2)))
endif
!------------------------------------------------------------
r=vector(2)
cos_lat=cos(vector(1))
dvds(0:2)= normalize(bp)/[r*cos_lat, r, 1.]
!------------------------------------------------------------
if (present(rk_first)) then
	if (site%CurlBFlag) then
		if (yinflag) then
			forall(i=0:2) site%CurlB_yin(i)=sum(weight*pole%CurlBvec(i, round(:,0), round(:,1), round(:,2)))
		else
			forall(i=0:2) site%CurlB(i)=sum(weight*CurlBvec(i, round(:,0), round(:,1), round(:,2)))
		endif
	endif
	if(site%AFlag) then
		if (yinflag) then
			forall(i=0:2) site%A_yin(i)=sum(weight*pole%Avec(i, round(:,0), round(:,1), round(:,2)))
		else
			forall(i=0:2) site%A(i)=sum(weight*Avec(i, round(:,0), round(:,1), round(:,2)))
		endif
	endif
	
	if (site%scottLaunch) call set_u_v(bp, vector(3:5), vector(6:8))

	! provide site%B/CurlB/A
	if (site%yinFlag) call cal_yinyang(site, .true., .false.)

	if (rk_first) return ! interpolate_foot is true

	if (site%yinFlag) then
		da(0:1)=pole%darc
		da(2)=axis(2)%da(round(0,2))	
	else
		forall(i=0:2) da(i)=axis(i)%da(round(0,i))
	endif
	site%ds_factor=norm2s(dvds(0:2)/da)
		
endif
!------------------------------------------------------------
if (site%scottFlag) then
	if (dbdc_field_Flag) then
		if (yinflag) then
			forall(i=0:2, j=0:2) &
			dbdcp(i,j)=sum(weight*pole%dbdc_field(i, j, round(:,0), round(:,1), round(:,2)))
		else
			forall(i=0:2, j=0:2) &
			dbdcp(i,j)=sum(weight*dbdc_field(i, j, round(:,0), round(:,1), round(:,2)))
		endif
	else
	!this output is identical as the upper, it don't require dbdc_field while takes more time
		do k=0,1
		do j=0,1
		do i=0,1
			if (weight(i,j,k) .ne. 0.0) then
				if (yinflag) then
					call dbdc_grid_pole(round(i,0), round(j,1), round(k,2), dbdc_cell(:,:,i,j,k), southflag)
				else
					call dbdc_grid(round(i,0), round(j,1), round(k,2), dbdc_cell(:,:,i,j,k))
				endif
			else
				!avoid NaN
				dbdc_cell(:,:,i,j,k)=0.0
			endif
		enddo
		enddo
		enddo
		forall(i=0:2, j=0:2) dbdcp(i,j)=sum(weight*dbdc_cell(i,j,:,:,:))
	endif
	
	sin_lat=sin(vector(1))

	dbdcp(0,  :)=dbdcp(0,  :) /(r*cos_lat)
	dbdcp(1,  :)=dbdcp(1,  :) / r
	dbdcp(:,  0)=dbdcp(:,  0) + [-dvds(1)*sin_lat/cos_lat + dvds(2)/r, dvds(0)*sin_lat, -dvds(0)*cos_lat]
	dbdcp(1:2,1)=dbdcp(1:2,1) + [dvds(2)/r, - dvds(1)]

	forall(i=0:2, j=1:2) dvds(3*j+i)= dot_product(vector(3*j:3*j+2), dbdcp(:,i))
endif
!------------------------------------------------------------
! yinflag determines the real field used for interpolation in this subroutine, 
! site%yinflag determines the vector/dvds used for RK step,
! if yinflag do not provide the corresponding vector/dvds for RK step, a conversion should be executed

! if present(rk_first) is .true., site%yinflag and yinFlag are identical
! if present(rk_first) is .false. and site%yinflag is .false., yinFlag is .false. already
if (.not. yinFlag .and. site%yinflag) call cal_yinyang(site, .false., .true.)

end subroutine interpolate_spherical


subroutine interpolate_cartesian(site, rk_first)
implicit none
type(site_info), target:: site
logical, optional:: rk_first
integer:: round(0:1,0:2), i, j, k
real:: weight(0:1,0:1,0:1), dbdc_cell(0:2,0:2,0:1,0:1,0:1), dbdcp(0:2,0:2), da(0:2), Ap(0:2)
!------------------------------------------------------------
call round_weight(site%v(0:2), round, weight)
forall(i=0:2) site%B(i)=sum(weight* Bvec(i, round(:,0), round(:,1), round(:,2)))
site%dvds(0:2)=normalize(site%b)
!------------------------------------------------------------
if (present(rk_first)) then
	if (CurlBvec_Flag) &
	forall(i=0:2) site%CurlB(i)=sum(weight*CurlBvec(i, round(:,0), round(:,1), round(:,2)))
	if (A_input) &
	forall(i=0:2) site%A(i)=sum(weight*Avec(i, round(:,0), round(:,1), round(:,2)))

	if (rk_first) return ! interpolate_foot is true

	if (stretchFlag) then
		forall(i=0:2) da(i)=axis(i)%da(round(0,i))
		site%ds_factor=norm2s(site%dvds(0:2)/da)
	else
		site%ds_factor=1.0
	endif

	if (site%scottLaunch) call set_u_v(site%b, site%v(3:5), site%v(6:8))
endif
!--------------------------------------------------------------------------
if (site%scottFlag) then
	if (dbdc_field_Flag) then
		forall(i=0:2, j=0:2) dbdcp(i,j)=sum(weight*dbdc_field(i, j, round(:,0), round(:,1), round(:,2)))
	else
	!this output is identical as the upper, it don't require dbdc_field while takes more time
		do k=0,1
		do j=0,1
		do i=0,1
			if (weight(i,j,k) .ne. 0.0) then
				call dbdc_grid(round(i,0), round(j,1), round(k,2), dbdc_cell(:,:,i,j,k))
			else
				!avoid NaN
				dbdc_cell(:,:,i,j,k)=0.0
			endif
		enddo
		enddo
		enddo
		forall(i=0:2, j=0:2) dbdcp(i,j)=sum(weight*dbdc_cell(i,j,:,:,:))
	endif

	forall(i=0:2, j=1:2) site%dvds(3*j+i)= dot_product(site%v(3*j:3*j+2), dbdcp(:,i))
endif

end subroutine interpolate_cartesian


subroutine locate_path_r(vp, site0, site_r, sign_dt, vr2vp)
implicit none
logical:: repeat_flag
integer:: sign_dt, it
real:: vp(0:2), vp_car(0:2), vr2vp(0:2), distance0, dt, a, b, c, f, r0(0:2), r1(0:2)
type(site_info), target:: site0, site_r, site_a, site_b
type(site_info), pointer:: site, site1, site_tmp
!------------------------------------------------------------
site_a=site0
site_b=site_r

site_a%Aflag=.false.
site_a%CurlBFlag=.false.
site_b%Aflag=.false.
site_b%CurlBFlag=.false.

site => site_a
site1=> site_b
!------------------------------------------------------------
if (spherical) then
	vp_car=vp_spherical2car(vp)
	r0=vp_spherical2car(site %v(0:2))-vp_car
	r1=vp_spherical2car(site1%v(0:2))-vp_car
else
	r0=site %v(0:2)-vp
	r1=site1%v(0:2)-vp
endif

 a=   dot_product(r0-r1, r0-r1)
 b=2.*dot_product(r0-r1, r1)
 c=   dot_product(r1, r1) - r_local_square
 f=(-b-sqrt(b**2.-4.*a*c))/(2.*a)
! this is unnecessary
! if (f .gt. 1. .or. f .lt. 0.) f=(-b+sqrt(b**2.-4.*a*c))/(2.*a)

dt= sign_dt*abs((1.-f)*distance(site%v(0:2), site1%v(0:2))*site%ds_factor-0.5*min_step)
call RK4(dt, site, site1)

if (distance(vp, site1%v(0:2)) .lt. r_local) then
	site_tmp => site
	site     => site1
	site1    => site_tmp
	call interpolate(site, .false.)
endif
!------------------------------------------------------------
repeat_flag=.true.
dt= sign_dt * min_step
it=0
do while (repeat_flag)
	call RK4(dt, site, site1)
	it= it+1
	distance0=distance(vp, site1%v(0:2))
	repeat_flag= distance0 .lt. r_local .and. it .le. 1000 .and. inside(site1%v(0:2))
	if (repeat_flag) then
		site_tmp => site
		site     => site1
		site1    => site_tmp
		call interpolate(site, .false.)
	endif
enddo

if (distance0 .lt. r_local) then
	site => site1
	site1=> site_r
endif
!------------------------------------------------------------
if (spherical) then
	r0=vp_spherical2car(site %v(0:2))-vp_car
	r1=vp_spherical2car(site1%v(0:2))-vp_car
else
	r0=site %v(0:2)-vp
	r1=site1%v(0:2)-vp
endif

 a=   dot_product(r0-r1, r0-r1)
 b=2.*dot_product(r0-r1, r1)
 c=   dot_product(r1, r1) - r_local_square
 f=(-b-sqrt(b**2.-4.*a*c))/(2.*a)

dt= sign_dt*(1.-f)*distance(site%v(0:2), site1%v(0:2))*site%ds_factor
call RK4(dt, site, site_r)

if (spherical) then
	vr2vp= vp_spherical2car(site_r%v(0:2))-vp_car
else
	vr2vp= site_r%v(0:2)-vp
endif

end subroutine locate_path_r


subroutine trace_bline(vp, info)
! vp: vector position of the launch point
implicit none
logical:: intFlag, CurlBFlag, AFlag, identical, repeat_flag, interpolate_foot, exist_vr
integer:: i, sign_down, sign_up, sign_forward, it, sign_dt, rb, e_index, s_index
real:: vp(0:2), dt, dt_executed, step_this, tol_this, dL, int2private(0:9), &
Bn_s, Bn_e, us(0:2), ue(0:2), vs(0:2), ve(0:2), us1(0:2), ue1(0:2), vs1(0:2), ve1(0:2), &
bs2(0:2), be2(0:2), us2(0:2), ue2(0:2), vs2(0:2), ve2(0:2), &
b_car(0:2), cos_p(0:1), sin_p(0:1), incline, vr2vp(0:2), vr(0:2), brn
real, pointer:: bp(:), bs(:), be(:)
type(line_info), target:: info
type(site_info), target:: site_a, site_b, site_p, site_s, site_e, site_r
type(site_info), pointer:: site, site1, site_tmp
!------------------------------------------------------------
intFlag = info%get .and. privateFlag
if (intFlag) int2private=0.
CurlBFlag = info%get .and. CurlBvec_Flag
AFlag = info%get .and. A_input
if (info%q_local_Flag) then
	info%local_s_flag=.false.
	info%local_e_flag=.false.
endif
!------------------------------------------------------------
site_p%v(0:2)=vp
site_p%scottFlag=info%scottFlag
site_p%scottLaunch=info%scottFlag
site_p%AFlag = AFlag
site_p%CurlBFlag= CurlBFlag
site_p%yinFlag=.false. ! this would be changed by the following command
call interpolate(site_p, .false.)
if (intFlag) call privates(site_p)
bp => info%bp
bs => info%bs
be => info%be
bp =  site_p%b
if (site_p%CurlBFlag) info%CurlBp=site_p%CurlB
!------------------------------------------------------------
if (.not. inside(vp)) then
	info%rbs=7; info%rFs=NaN
	info%rbe=7; info%rFe=NaN
	
	if (info%get) then
		info%int_private=NaN
		bs =NaN
		be =NaN
	endif

	if (info%scottFlag) then
		info%q = NaN
		if (.not. diff_flag) info%q_perp = NaN
	endif

	if (info%path_out) then
		info%its=0
		info%ite=0
		info%path(:,0) = vp
		if (info%loopB_out) info%loopB(:,0) = site_p%B
		if (info%loopCurlB_out) info%loopCurlB(:,0) = site_p%CurlB
	endif
	return
endif
!------------------------------------------------------------
if (info%get .and. diff_flag) then
	if (Normal_index .eq. -1) then
		if (spherical) then
			cos_p=cos(vp(0:1)); sin_p=sin(vp(0:1))
			
			b_car(0)=-bp(0)*sin_p(0) + (-bp(1)*sin_p(1) + bp(2)*cos_p(1))*cos_p(0)
			b_car(1)= bp(0)*cos_p(0) + (-bp(1)*sin_p(1) + bp(2)*cos_p(1))*sin_p(0)
			b_car(2)=                    bp(1)*cos_p(1) + bp(2)*sin_p(1)
		else
			b_car   = bp
		endif
		info%bnp=dot_product(b_car, info%ev3)
	else
		info%bnp=bp(normal_index)
	endif

	incline=abs(info%bnp/norm2s(bp))
	
	info%tangent = .not. (incline .ge. min_incline)

	if (info%tangent) then
		if (lie_boundary(vp)) then
		! In most cases, lie_boundary(vp) is .true. makes key_diff become .true.
		! the rest cases will trace with Scott (2017)
			incline=min_incline
			info%tangent=.false.
		else
		! In most cases, info%tangent is .true. makes key_trace4 become .true., 
		! the rest cases will trace with Scott (2017)
			incline=1.
		endif
	endif
else
	incline=1.
endif

if (.not. inclineFlag) incline=1.

! Equations (20) (21) in Zhang (2022)
if (RK4flag) then
	step_this= step * incline
	if (step_this .lt. min_step) step_this= min_step
else
	step_this= min_step
	tol_this = tol * incline **1.5
endif
!------------------------------------------------------------
! check if only one direction should be traced
if (.not. traceflag) then
	! info%rbs/info%rbe could be changed after the following check
	info%rbs=0
	info%rbe=0
endif

sign_down = -1
sign_up   =  1

do i=0, 2
	if (periodFlag(i)) cycle
		
	if      (vp(i) .eq. pmin(i)) then
		if (south_pole .and. i .eq. 1) cycle
		rb=5-2*i
		sign_forward= 1
	else if (vp(i) .eq. pmax(i)) then
		if (north_pole .and. i .eq. 1) cycle
		rb=6-2*i
		sign_forward=-1
	else
		cycle
	endif
	
	if      (bp(i)*sign_forward .gt. 0.) then
		sign_down= 1
		info%rbs=rb
		site_s=site_p
		if (info%path_out) info%its=0
	else if (bp(i)*sign_forward .lt. 0.) then
		sign_up  =-1
		info%rbe=rb
		site_e=site_p
		if (info%path_out) info%ite=0
	endif
enddo
!------------------------------------------------------------
if (.not. traceflag) then
	info%rFs=vp
	info%rFe=vp
	if (info%get) then
		bs=bp
		be=bp
	endif
	if (info%scottFlag) then
		info%q = NaN
		if (.not. diff_flag) info%q_perp = NaN
	endif
	return
endif 
!------------------------------------------------------------
if (info%path_out) then
	info%path(:,0)=vp
	if (info%loopB_out) info%loopB(:,0)=site_p%B
	if (info%loopCurlB_out) info%loopCurlB(:,0)=site_p%CurlB
endif
!------------------------------------------------------------
site_a%scottLaunch=.false.
site_a%scottFlag=info%scottFlag
site_a%AFlag=AFlag
site_a%CurlBFlag= CurlBFlag
site_b%scottLaunch=.false.
site_b%scottFlag=info%scottFlag
site_b%AFlag=AFlag
site_b%CurlBFlag= CurlBFlag

do sign_dt = sign_down, sign_up, 2

	it= 0
	dt= step_this*sign_dt
	site  => site_p
	site1 => site_a
	repeat_flag=.true.
	interpolate_foot= .false. 
	exist_vr=.false.
	
	do while (repeat_flag)
		if (RK4flag) then
			call RK4  (dt, site, site1)
		else
			call RKF45(dt, site, site1, tol_this, dt_executed)
		endif

		repeat_flag = inside(site1%v(0:2)) .and. (abs(it+sign_dt) .lt. maxsteps)

		if (.not. repeat_flag) then
			if (RK4flag) dt_executed=dt
			call correct_foot(it .eq. 0, dt_executed, site, site1, rb, identical)
			if (identical) site1 => site
			
			! if rb .eq. 8, NaN found in site1%v(0:2), then site1 should not be interpolated
			! if key_trace4, other 4 lines don't need this interpolation at foots
			interpolate_foot= rb .ne. 8 .and. .not. identical .and. (info%get .or. info%scottFlag)
		endif

		if (repeat_flag .or. interpolate_foot) then

			it=it+sign_dt

			call interpolate(site1, interpolate_foot)

			if (info%path_out) then
				info%path(:,it)=site1%v(0:2)
				if (info%loopB_out) info%loopB(:,it)=site1%B
				if (info%loopCurlB_out) info%loopCurlB(:,it)=site1%CurlB
			endif

			if (intFlag) then
				dL = distance(site%v(0:2), site1%v(0:2))
				! dL = dt_executed/site%ds_factor, this is actually better than the upper,
				! while is is more difficult to explain to the public

				call privates(site1)
				int2private = int2private + (site%private+site1%private)*dL
			endif

			if (info%q_local_Flag .and. (.not. exist_vr)) then

				exist_vr= distance(site1%v(0:2), vp) .ge. r_local
				
				if (exist_vr) then
					site_r=site1
					
					call locate_path_r(vp, site, site_r, sign_dt, vr2vp)
					
					call interpolate(site_r, .true.)
					if (spherical) then
						cos_p=cos(site_r%v(0:1))
						sin_p=sin(site_r%v(0:1))
						b_car(0)=-site_r%B(0)*sin_p(0) + (-site_r%B(1)*sin_p(1) + site_r%B(2)*cos_p(1))*cos_p(0)
						b_car(1)= site_r%B(0)*cos_p(0) + (-site_r%B(1)*sin_p(1) + site_r%B(2)*cos_p(1))*sin_p(0)
						b_car(2)=                          site_r%B(1)*cos_p(1) + site_r%B(2)*sin_p(1)
					else
						b_car = site_r%B
					endif
					brn=dot_product(b_car, normalize(vr2vp))
					if (sign_dt .eq. -1) then
						info%brn_s=brn
						info%local_s_flag=.true.
						if (info%scottFlag) then
							bs2=site_r%B
							us2=site_r%v(3:5)-dot_product(site_r%v(3:5),bs2)/dot_product(bs2,bs2)*bs2
							vs2=site_r%v(6:8)-dot_product(site_r%v(6:8),bs2)/dot_product(bs2,bs2)*bs2
						endif
					else
						info%brn_e=brn
						info%local_e_flag=.true.
						if (info%scottFlag) then
							be2=site_r%B
							ue2=site_r%v(3:5)-dot_product(site_r%v(3:5),be2)/dot_product(be2,be2)*be2
							ve2=site_r%v(6:8)-dot_product(site_r%v(6:8),be2)/dot_product(be2,be2)*be2
						endif
					endif
				endif !exist_vr
			endif
		endif

		if (repeat_flag) then
			if (abs(it) .eq. 1) then
				site  => site_a
				site1 => site_b
			else
				! switch site, site1
				site_tmp => site
				site     => site1
				site1    => site_tmp
			endif
		else
			if (sign_dt .eq. -1) then
				info%its=it
				site_s=site1
				info%rbs=rb
			else
				info%ite=it
				site_e=site1
				info%rbe=rb
			endif
		endif
	enddo
enddo
!------------------------------------------------------------
info%rFs=site_s%v(0:2)
info%rFe=site_e%v(0:2)
bs=site_s%b
be=site_e%b

info%s_yinflag=site_s%yinflag
if (site_s%yinflag) info%rFs_yin=site_s%v_yin(0:2)
info%e_yinflag=site_e%yinflag
if (site_e%yinflag) info%rFe_yin=site_e%v_yin(0:2)

if (info%get) then
	if (intFlag) info%int_private= int2private * 0.5
	if (targetCurlB_out) then
		info%CurlBs=site_s%CurlB
		info%CurlBe=site_e%CurlB
	endif
endif
!------------------------------------------------------------
! Scott_2017_ApJ_848_117
if (.not. info%scottFlag) return

us=site_s%v(3:5)
vs=site_s%v(6:8)
ue=site_e%v(3:5)
ve=site_e%v(6:8)

if (all([[info%rbs, info%rbe] .ge. 1, [info%rbs, info%rbe] .le. 6])) then
	! if site_s%yinflag/site_e%yinflag, s_index/e_index can only be 2
	s_index=(6-info%rbs)/2; Bn_s=bs(s_index)
	e_index=(6-info%rbe)/2; Bn_e=be(e_index)
	us1=us-us(s_index)/Bn_s*bs
	vs1=vs-vs(s_index)/Bn_s*bs
	ue1=ue-ue(e_index)/Bn_e*be
	ve1=ve-ve(e_index)/Bn_e*be

	info%q= abs(dot_product(ue1,ue1)*dot_product(vs1,vs1)  &
		   +    dot_product(us1,us1)*dot_product(ve1,ve1)  &
		   -2.0*dot_product(ue1,ve1)*dot_product(us1,vs1)) &
		     / (dot_product(bp, bp) / abs(Bn_s*Bn_e))
else
	info%q = NaN
endif
!------------------------------------------------------------
if (info%q_local_Flag) then
	if (info%local_s_flag) then
		Bn_s=norm2s(bs2)
	else
		us2=us1
		vs2=vs1
	endif

	if (info%local_e_flag) then
		Bn_e=norm2s(be2)
	else
		ue2=ue1
		ve2=ve1
	endif
	info%q_local=abs(dot_product(ue2,ue2)*dot_product(vs2,vs2)   &
				+    dot_product(us2,us2)*dot_product(ve2,ve2)   &
				-2.0*dot_product(ue2,ve2)*dot_product(us2,vs2))  &
				  / (dot_product(bp, bp) / abs(Bn_s*Bn_e))
endif
!------------------------------------------------------------
if (.not. diff_flag) then

	us1=us-dot_product(us,bs)/dot_product(bs,bs)*bs
	vs1=vs-dot_product(vs,bs)/dot_product(bs,bs)*bs
	ue1=ue-dot_product(ue,be)/dot_product(be,be)*be
	ve1=ve-dot_product(ve,be)/dot_product(be,be)*be

	info%q_perp=abs(dot_product(ue1,ue1)*dot_product(vs1,vs1)   &
				+   dot_product(us1,us1)*dot_product(ve1,ve1)   &
			   -2.0*dot_product(ue1,ve1)*dot_product(us1,vs1))  &
				 / (dot_product(bp, bp) / (norm2s(bs)*norm2s(be)))
endif

end subroutine trace_bline

end module trace


subroutine show_time(percent)
implicit none
real:: percent
integer:: times(8)
!------------------------------------------------------------
call date_and_time(VALUES=times)
print '(9X,F6.2,1X,"%",7X,I2.2,":",I2.2,":",I2.2)', percent, times(5:7)
end subroutine show_time

module compute
use trace
implicit none
integer:: iend, jend, nq1, nq2, nq3, ijend(1:2)
integer(1), target, allocatable:: rbs(:,:), rbe(:,:)
integer(1), allocatable:: rboundary(:,:)
integer(2), allocatable:: sign2d(:,:)! IDL do not have the type of 8-bit signed integer
real:: deltas(-1:2), delta_i, delta_j, xreg(0:1), yreg(0:1), zreg(0:1), ev3(0:2)
real, allocatable:: q(:,:), q_perp(:,:), &
seed(:,:,:), b_layer(:,:,:), CurlB_layer(:,:,:), bnp2d(:,:), &
bs_layer(:,:,:), be_layer(:,:,:), CurlBs_layer(:,:,:), CurlBe_layer(:,:,:), &
q_local(:,:), brn_s(:,:), brn_e(:,:)
real, allocatable, target:: rFs(:,:,:), rFe(:,:,:), rFs_yin(:,:,:), rFe_yin(:,:,:)
logical:: vflag, bflag, cflag, sFlag, scottFlag, diff_seed, pole_j0, pole_jend, &
targetB_flag, sign2dFlag, allocate_path, path_out, loopB_out, q_local_Flag
logical, allocatable:: tangent(:,:), local_s_flag(:,:), local_e_flag(:,:)
logical, allocatable, target:: s_yinFlag(:,:), e_yinFlag(:,:)
type line
	real, allocatable:: path(:,:), loopB(:,:), loopCurlB(:,:)
endtype line
type(line), allocatable :: lines(:)

type int_line
	real, allocatable:: data(:,:)
endtype int_line
type(int_line) :: int_private(0:9)

integer, allocatable:: index_seed(:), loop_size(:)

contains

subroutine q_bridge(i, j)
implicit none
integer:: i, j, ip, id, its, ite
logical:: pole_j
type(line_info):: info
real:: weight(0:1,0:1,0:1)
integer:: round(0:1, 0:2)
!------------------------------------------------------------
pole_j= (pole_j0 .and. j .eq. 0) .or. (pole_jend .and. j .eq. jend)
if (pole_j .and. i .ne. 0) return
!------------------------------------------------------------
info%path_out= allocate_path
info%loopB_out= loopB_out
info%loopCurlB_out= loopCurlB_out
info%q_local_Flag= q_local_Flag

if (info%path_out) then
	id=i+j*nq1
	allocate(info%path(0:2,-maxsteps:maxsteps))
	if (info%loopB_out) allocate(info%loopB(0:2,-maxsteps:maxsteps))
	if (info%loopCurlB_out) allocate(info%loopCurlB(0:2,-maxsteps:maxsteps))
endif

info%scottFlag = scottFlag .or. pole_j
info%get = .true.
if (diff_flag) then
	if (diff_seed) then
		info%ev3 = seed_ev3(i,j)
	else
		info%ev3 = ev3
	endif
endif

call trace_bline(seed(:,i,j), info)

if (info%path_out) then
	its=info%its
	ite=info%ite
	index_seed(id)= -its
	loop_size(id)= ite-its+1
	allocate(lines(id)%path(0:2, 0:ite-its))
	lines(id)%path = info%path(:, its:ite)
	deallocate(info%path)
	if (info%loopB_out) then 
		allocate(lines(id)%loopB(0:2, 0:ite-its))
		lines(id)%loopB = info%loopB(:, its:ite)
		deallocate(info%loopB)
	endif
	if (info%loopCurlB_out) then
		allocate(lines(id)%loopCurlB(0:2, 0:ite-its))
		lines(id)%loopCurlB = info%loopCurlB(:, its:ite)
		deallocate(info%loopCurlB)
	endif
endif
!------------------------------------------------------------
rboundary(i, j)= 10*info%rbs + info%rbe
rFs(:, i, j)=info%rFs
rFe(:, i, j)=info%rFe

b_layer(:, i, j)=info%bp

if (targetB_flag) then
	bs_layer(:, i, j)=info%bs
	be_layer(:, i, j)=info%be
endif

if (CurlB_out) CurlB_layer(:, i, j)=info%CurlBp
if (privateFlag) then
	do ip =0, 9 
		if (int_private_out(ip)) &
		int_private(ip)%data(i, j)=info%int_private(ip)
	enddo
endif
if (targetCurlB_out) then 
	CurlBs_layer(:, i, j)=info%CurlBs
	CurlBe_layer(:, i, j)=info%CurlBe
endif

if (scottFlag) then
	q(i,j)=info%q
	q_perp(i,j)=info%q_perp
	if (q_local_flag) q_local(i,j)=info%q_local
else if (diff_flag) then
	if (pole_j) then
		q(i,j)=info%q
		if (q_local_flag) q_local(i,j)=info%q_local
	endif
	rbs(i, j)=info%rbs
	rbe(i, j)=info%rbe
	tangent(i,j)=info%tangent
	bnp2d(i, j)=info%bnp

	if (q_local_flag) then
		local_s_flag(i,j)=info%local_s_flag
		local_e_flag(i,j)=info%local_e_flag
		brn_s(i,j)=info%brn_s
		brn_e(i,j)=info%brn_e
	endif

	s_yinFlag(i,j)=info%s_yinFlag
	e_yinFlag(i,j)=info%e_yinFlag

	if (south_pole .or. north_pole) then
		if (info%s_yinFlag) rFs_yin(:,i,j)=info%rFs_yin
		if (info%e_yinFlag) rFe_yin(:,i,j)=info%rFe_yin
	endif
Endif

if (sign2dFlag) then
	if      (info%bp(2) .gt. 0.0) then
		sign2d(i,j)= 1
	else if (info%bp(2) .lt. 0.0) then
		sign2d(i,j)=-1
	else ! info%bp(2) = 0 or NaN
		sign2d(i,j)= 0
	endif
Endif

END subroutine q_bridge


subroutine q_diff(i,j)
! method 3 of Pariat (2012), some problematic sites are filled with Scott (2017)
implicit none
logical:: key_rb16, key_nB, key_trace4, key_diff, local_trace4, local_diff, local_launch, exist_vr
logical, pointer:: yinFlag(:,:)
integer:: i, j, k, sign_dt, it, it_end, s_index, e_index, &
i_dim, i_diff, id4(1:4), ij(1:2), diff_index(0:2,1:2)
integer, pointer:: dims(:)
integer, target :: dims_s(1:2), dims_e(1:2)
real:: delta_diff(1:2), bn_square, Bn_s, Bn_e, Dmatrix(1:2, 1:2),  &
rF4_s(1:2, 1:4), rF4_e(1:2, 1:4), u0(0:2), v0(0:2), &
g_e(1:2), g_s(1:2), g_e0(0:2), g_s0(0:2), gh, cos_p(0:1), sin_p(0:1), &
cos_tmp, sin_tmp, rF4_s_local(0:2, 1:4), rF4_e_local(0:2, 1:4), vr2vp(0:2), &
rf3(1:2, 0:2), seed3(0:2, 0:2), coef(0:2,1:2), d0, d1, arrow_seed(0:2, 1:2), vp_yin(0:2)
real, pointer:: bp_car(:), vp_car(:), vp4_car(:,:), rF_tmp(:,:,:), diff(:,:), rF_yin(:,:,:)
real, target :: diff_s(1:2, 1:2), diff_e(1:2, 1:2), bp(0:2), bp_car_tmp(0:2), &
vp(0:2), vp_car_tmp(0:2), vp4(0:2, 1:4), vp4_car_tmp(0:2, 1:4)
type(line_info):: info
type(site_info), target :: site, site_r
!------------------------------------------------------------
if (.not. inside(seed(:,i,j))) then
	q(i,j)=NaN
	if (q_local_Flag) q_local(i,j)=NaN
	return
endif 

if ((pole_j0 .and. j .eq. 0) .or. (pole_jend .and. j .eq. jend)) return
!------------------------------------------------------------
if (tangent(i,j)) then
	key_nB=.true.
else
	ij=[i,j]
	do i_diff=1, 2
		if (ij(i_diff) .eq. 0) then
			diff_index(:,i_diff)=ij(i_diff) + [0,1,2]
			if (.not. diff_seed) coef(:,i_diff) = [-3., 4., -1.]
		else if(ij(i_diff) .eq. ijend(i_diff)) then
			diff_index(:,i_diff)=ij(i_diff) + [-2,-1,0]
			if (.not. diff_seed) coef(:,i_diff) = [1., -4., 3.]
		else
			diff_index(:,i_diff)=ij(i_diff) + [-1,0,1]
			if (.not. diff_seed) coef(:,i_diff) = [-1., 0., 1.]
		endif
	enddo

	! for a tangent field line, step and tol is not adjusted by incline
	! then the accuracy of rs/re is not good enough for finite difference 
	key_nB = inclineFlag .and. any([tangent(diff_index(:, 1), j), tangent(i, diff_index(:, 2))])

	if (diff_seed .and. .not. key_nB) then
		do i_diff = 1, 2
			if (i_diff .eq. 1) then
				seed3=seed(:, diff_index(:,1), j)
			else
				seed3=seed(:, i, diff_index(:,2))
			endif
			d0=distance(seed3(:,0), seed3(:,1))
			d1=distance(seed3(:,1), seed3(:,2))

			if (d0 .le. 0. .or. d1 .eq. 0.) then
				key_nB=.true.
				exit
			endif

			arrow_seed(:, i_diff)=seed3(:,2)-seed3(:,0)
			delta_diff(i_diff) = (d0+d1)/2.

			if (diff_index(0, i_diff) .eq. 0) then
				coef(0,i_diff)=-(2.+ d1/d0)
				coef(2,i_diff)=-d0/d1
			else if (diff_index(2, i_diff) .eq. ijend(i_diff)) then
				coef(0,i_diff)= d1/d0
				coef(2,i_diff)= 2.+ d0/d1
			else
				coef(0,i_diff)=-d1/d0
				coef(2,i_diff)= d0/d1
			endif
			coef(1,i_diff) = - coef(0,i_diff) - coef(2,i_diff)
		enddo

		cos_tmp = cos2vector(arrow_seed(:,1), arrow_seed(:,2), gh) ! gh=sin_tmp

		if (delta_diff(1)* 100. < delta_diff(2) .or. &
		    delta_diff(2)* 100. < delta_diff(1) .or. gh .le. 0.01) key_nB=.true.

	endif
endif
!------------------------------------------------------------
key_rb16 = all([[rbs(i,j), rbe(i,j)] .ge. 1, [rbs(i,j), rbe(i,j)] .le. 6])

if (key_rb16) then
	s_index=(6-rbs(i,j))/2
	e_index=(6-rbe(i,j))/2
	dims_s=mod(s_index+[1,2], 3)
	dims_e=mod(e_index+[1,2], 3)
endif
!------------------------------------------------------------
! these values will be checked for a few times later, and then could be .false.
key_trace4 =       key_nB .and. key_rb16
key_diff   = .not. key_nB .and. key_rb16

if (q_local_flag) then
	local_trace4 =     key_nB .and. (local_s_flag(i,j) .or. local_e_flag(i,j))
	local_diff = .not. key_nB .and. (local_s_flag(i,j) .or. local_e_flag(i,j)) &
    .and. ((i .ne. 0) .and. (i .ne. iend) .and. (j .ne. 0) .and. (j .ne. jend))
	site  %scottFlag=.false.
	site  %scottLaunch=.false.
	site_r%scottFlag=.false.
	site_r%scottLaunch=.false.
else
	local_trace4 = .false.
	local_diff   = .false.
endif

local_launch = q_local_flag .and. &
((i .ne. 0) .and. (i .ne. iend) .and. (j .ne. 0) .and. (j .ne. jend))
!------------------------------------------------------------
vp=seed(:, i, j)

! cos(), sin() are pretty time consuming
if (spherical) cos_p=cos(vp(0:1))
!------------------------------------------------------------
if (key_trace4 .or. local_trace4) then
! use the plane perpendicular to the field line to calculate Q

	bp= b_layer(:, i, j)
	bn_square= dot_product(bp, bp)
	gh=1.

	if (spherical) then
		bp_car  => bp_car_tmp
		vp_car  => vp_car_tmp
		vp4_car => vp4_car_tmp

		sin_p=sin(vp(0:1))

		bp_car(0)=-bp(0)*sin_p(0) + (-bp(1)*sin_p(1) + bp(2)*cos_p(1))*cos_p(0)
		bp_car(1)= bp(0)*cos_p(0) + (-bp(1)*sin_p(1) + bp(2)*cos_p(1))*sin_p(0)
		bp_car(2)=                    bp(1)*cos_p(1) + bp(2)*sin_p(1)
		
		vp_car = vp(2)* [cos_p(1)*[cos_p(0), sin_p(0)], sin_p(1)]
	else
		bp_car  => bp
		vp_car  => vp
		vp4_car => vp4
	endif

	if (diff_seed) then
		if (j .eq. 0) then 
			delta_diff=distance(vp, seed(:,i,j+1))
		else
			delta_diff=distance(vp, seed(:,i,j-1))
		endif
	else
		if (spherical) then
			select case(normal_index)
			case(-1)
				delta_diff=[vp(2)*deltas(-1), deltas(2)]
			case(0)
				delta_diff=[vp(2)*deltas( 1), deltas(2)]
			case(1)
				delta_diff=[vp(2)*cos_p(1)*deltas(0), deltas(2)]
			case(2)
				delta_diff=[vp(2)*cos_p(1)*deltas(0), vp(2)*deltas(1)]
			end select
		else
			delta_diff=[delta_i, delta_j]
		endif
	endif

	! this choice is still unclear for all cases
	! choose it just from practice 
	delta_diff=delta_diff/2.

	call set_u_v(bp_car, u0, v0)
	vp4_car(:,1)= vp_car+delta_diff(1)*u0
	vp4_car(:,2)= vp_car+delta_diff(2)*v0
	vp4_car(:,3)= vp_car-delta_diff(1)*u0
	vp4_car(:,4)= vp_car-delta_diff(2)*v0

	if (any(delta_diff .gt. r_local)) then
		local_launch=.false.
		local_trace4=.false.
	endif
	if (spherical) then
		do k=1, 4
			vp4(:,k)=vp_car2spherical(vp4_car(:,k))
		enddo
	endif

	do k=1, 4
		if (.not. inside(vp4(:,k))) then
			key_trace4=.false.
			local_trace4=.false.
			exit
		endif
	enddo

	if (key_trace4 .or. local_trace4) then
		info%scottFlag= .false.
		info%q_local_Flag= .false.
		info%get= .false.
		info%path_out= local_trace4
		if (info%path_out) then
			allocate(info%path(0:2,-maxsteps:maxsteps))
			info%loopB_out= .false.
			info%loopCurlB_out= .false.
		endif

		do k= 1, 4
			call trace_bline(vp4(:,k), info)

			! key_trace4 shoud be check for 4 times
			if ((info%rbs .ne. rbs(i,j)) .or. (info%rbe .ne. rbe(i,j))) key_trace4=.false.
			if (.not. key_trace4 .and. .not. local_trace4) exit

			if (key_trace4) then
				
				if (s_yinFlag(i,j)) then 
					if (info%s_yinFlag) then
						rF4_s(:,k)=info%rFs_yin(0:1)
					else
						call vp_yinyang(vp_yin, info%rFs, .false.)
						rF4_s(:,k)=vp_yin(0:1)
					endif
				else
					rF4_s(:,k)=info%rFs(dims_s)
				endif

				if (e_yinFlag(i,j)) then 
					if (info%e_yinFlag) then
						rF4_e(:,k)=info%rFe_yin(0:1)
					else
						call vp_yinyang(vp_yin, info%rFe, .false.)
						rF4_e(:,k)=vp_yin(0:1)
					endif
				else
					rF4_e(:,k)=info%rFe(dims_e)
				endif
			endif
			
			if (local_trace4) then
				do sign_dt = -1, 1, 2
					if (sign_dt .eq. -1) then
						if (.not. local_s_flag(i,j)) cycle
						it_end= info%its
					else
						if (.not. local_e_flag(i,j)) cycle
						it_end= info%ite
					endif

					do it = sign_dt, it_end, sign_dt

						exist_vr= distance(info%path(:, it), vp) .ge. r_local
						if (exist_vr) then
							site  %v(0:2)= info%path(:, it-sign_dt)
							site_r%v(0:2)= info%path(:, it)
							site%yinFlag=.false.
							call interpolate(site, .false.)
							call locate_path_r(vp, site, site_r, sign_dt, vr2vp)

							if (sign_dt .eq. -1) then
								rF4_s_local(:,k)= vr2vp
							else
								rF4_e_local(:,k)= vr2vp
							endif
							exit
						endif
						
						if (it .eq. it_end) then
							if (sign_dt .eq. -1) then
								local_s_flag(i,j)=.false.
							else
								local_e_flag(i,j)=.false.
							endif
						endif
					enddo ! it =0, it_end, sign_dt
				enddo ! sign_dt =-1, 1, 2
			endif ! local_trace4
		enddo ! k= 1, 4

		if (info%path_out) deallocate(info%path)

		if (key_trace4) then
			diff_s = rF4_s(:, 1:2)-rF4_s(:, 3:4)
			diff_e = rF4_e(:, 1:2)-rF4_e(:, 3:4)

			if (period_lon) then
			do i_dim= 1, 2
				if (dims_s(i_dim) .eq. 0 .and. .not. s_yinFlag(i,j)) &
				diff_s(i_dim, :)=modulo(diff_s(i_dim, :)+pi, two_pi)-pi
				if (dims_e(i_dim) .eq. 0 .and. .not. e_yinFlag(i,j)) &
				diff_e(i_dim, :)=modulo(diff_e(i_dim, :)+pi, two_pi)-pi
			enddo
			endif
		endif

	endif
endif
!------------------------------------------------------------
if (key_diff .or. local_diff) then

	bn_square= bnp2d(i,j)**2.
	if (.not. diff_seed) Then
		if (spherical) then
			select case(Normal_index)
			case(-1:0)
				gh= vp(2)**2.
			case(1)
				gh=(vp(2)*cos_p(1))**2.
			case(2)
				gh= vp(2)**4. * cos_p(1)**2.
			end select
		else
			gh=1.
		endif
		delta_diff = [delta_i, delta_j]
	endif
!------------------------------------------------------------
	if (any([rbs(diff_index(:,1),j),rbs(i,diff_index(:,2))] .ne. rbs(i,j)) .or. &
		any([rbe(diff_index(:,1),j),rbe(i,diff_index(:,2))] .ne. rbe(i,j))) key_diff = .false.

	if (key_diff) then
		do i_diff = 1, 2
			do sign_dt = -1, 1, 2
				if (sign_dt .eq. -1) then
					rF_tmp => rFs
					diff   => diff_s
					dims   => dims_s
					yinFlag=> s_yinFlag
					rF_yin => rFs_yin
				else
					rF_tmp => rFe
					diff   => diff_e
					dims   => dims_e
					yinFlag=> e_yinFlag
					rF_yin => rFe_yin
				endif
				
				if (yinFlag(i,j)) then
					if (i_diff .eq. 1) then
						do k=0, 2
							if (yinFlag(diff_index(k, 1), j)) then
								rF3(:, k)=rf_yin(0:1, diff_index(k, 1), j)
							else
								call vp_yinyang(vp_yin, rF_tmp(:, diff_index(k,1), j), .false.)
								rF3(:,k)=vp_yin(0:1)
							endif
						enddo
					else
						do k=0, 2
							if (yinFlag(i, diff_index(k, 2))) then
								rF3(:, k)=rf_yin(0:1, i, diff_index(k, 2))
							else
								call vp_yinyang(vp_yin, rF_tmp(:, i, diff_index(k, 2)), .false.)
								rF3(:,k)=vp_yin(0:1)
							endif
						enddo
					endif
				else
					if (i_diff .eq. 1) then
						rF3=rF_tmp(dims, diff_index(:,1), j)
					else
						rF3=rF_tmp(dims, i, diff_index(:,2))
					endif
				endif

				! the differences of longitute could be about 2 pi, then remap it to [-pi, pi)
				if (period_lon .and. .not. yinFlag(i,j)) then
				do i_dim= 1, 2
					if (dims(i_dim) .eq. 0) &
					rf3(i_dim, 1:2)=rf3(i_dim, 0)-pi + modulo(rf3(i_dim, 1:2)-(rf3(i_dim, 0)-pi), two_pi)
				enddo
				endif

				forall(i_dim=1:2) diff(i_dim, i_diff)= dot_product(coef(:, i_diff), rf3(i_dim, :))

			enddo
		enddo
	endif
!------------------------------------------------------------
	if (local_diff) then
		id4= i + j*nq1 + [-1, -nq1, 1, nq1]
		do k = 1, 4
			if (distance(lines(id4(k))%path(:, index_seed(id4(k))), vp) .gt. r_local) then
				 local_diff  = .false.
				 local_launch= .false.
				 exit
			endif
			do sign_dt = -1, 1, 2
				if (sign_dt .eq. -1) then
					if (.not. local_s_flag(i,j)) cycle
					it_end= 0
				else
					if (.not. local_e_flag(i,j)) cycle
					it_end= loop_size(id4(k))-1
				endif

				do it = index_seed(id4(k))+sign_dt, it_end, sign_dt

					exist_vr= distance(lines(id4(k))%path(:, it), vp) .ge. r_local
					if (exist_vr) then
						site  %v(0:2)= lines(id4(k))%path(:, it-sign_dt)
						site_r%v(0:2)= lines(id4(k))%path(:, it)
						site%yinFlag=.false.
						call interpolate(site, .false.)
						call locate_path_r(vp, site, site_r, sign_dt, vr2vp)
						
						if (sign_dt .eq. -1) then
							rF4_s_local(:,k)= vr2vp
						else
							rF4_e_local(:,k)= vr2vp
						endif
						exit
					endif
					
					if (it .eq. it_end) then
						if (sign_dt .eq. -1) then
							local_s_flag(i,j)=.false.
						else
							local_e_flag(i,j)=.false.
						endif
					endif
				enddo ! it =0, it_end, sign_dt
			enddo ! sign_dt =-1, 1, 2
		enddo ! k = 1, 4
	endif ! local_diff
endif
!------------------------------------------------------------
if (key_trace4 .or. key_diff) then

	if (spherical) then
		if (s_yinFlag(i,j)) then
			g_s =[(rFs(2,i,j)*cos(rFs_yin(1,i,j)))**2., rFs(2,i,j)**2]
		else
			g_s0=[(rFs(2,i,j)*cos(rFs(1,i,j)))**2., rFs(2,i,j)**2., 1.]
			g_s =g_s0(dims_s)
		endif
		if (e_yinFlag(i,j)) then
			g_e =[(rFe(2,i,j)*cos(rFe_yin(1,i,j)))**2., rFe(2,i,j)**2]
		else
			g_e0=[(rFe(2,i,j)*cos(rFe(1,i,j)))**2., rFe(2,i,j)**2., 1.]
			g_e =g_e0(dims_e)
		endif
	else
		g_s=1.; g_e=1.
	endif

	Bn_s=bs_layer(s_index,i,j)
	Bn_e=be_layer(e_index,i,j)

	Dmatrix(1,1)=  diff_e(1, 1)*diff_s(2, 2) - diff_e(1, 2)*diff_s(2, 1)
	Dmatrix(1,2)= -diff_e(1, 1)*diff_s(1, 2) + diff_e(1, 2)*diff_s(1, 1)
	Dmatrix(2,1)=  diff_e(2, 1)*diff_s(2, 2) - diff_e(2, 2)*diff_s(2, 1)
	Dmatrix(2,2)= -diff_e(2, 1)*diff_s(1, 2) + diff_e(2, 2)*diff_s(1, 1)

	q(i,j)=abs((g_e(1)*g_s(2)*Dmatrix(1,1)**2 +  &
		        g_e(1)*g_s(1)*Dmatrix(1,2)**2 +  &
		        g_e(2)*g_s(2)*Dmatrix(2,1)**2 +  &
		        g_e(2)*g_s(1)*Dmatrix(2,2)**2) * &
		        Bn_s*Bn_e/(bn_square* gh* (4.*delta_diff(1)*delta_diff(2))**2.))

else if (key_rb16) then
	info%scottFlag=.true.
	info%get=.false.
	info%path_out=.false.
	info%q_local_Flag=q_local_flag
	call trace_bline(vp, info)
	q(i,j)=info%q
else
	q(i,j)=NaN
endif

if (q_local_flag) then

	if ((local_diff .or. local_trace4) .and. local_s_flag(i,j)) then
		g_s = r_local_square

		diff_s(1, 1)= acos(cos2vector(rF4_s_local(:, 1), rF4_s_local(:, 3)))
		diff_s(2, 1)= 0.
		cos_tmp = cos2vector(rF4_s_local(:,3)-rF4_s_local(:,1), &
		                     rF4_s_local(:,4)-rF4_s_local(:,2), sin_tmp)
		diff_s(:, 2)= acos(cos2vector(rF4_s_local(:, 2), rF4_s_local(:, 4)))*[cos_tmp, sin_tmp]
		Bn_s=brn_s(i,j)
	else if (.not. ((key_trace4 .or. key_diff) .and. local_launch)) then
		diff_s=NaN
	endif

	! if (.not. ((local_diff .or. local_trace4) .and. local_s_flag(i,j))
	!     .and. ((key_trace4 .or. key_diff) .and. local_launch) )
	! diff_s is calculated from the boundary mapping

	if ((local_diff .or. local_trace4) .and. local_e_flag(i,j)) then
		g_e = r_local_square

		diff_e(1, 1)= acos(cos2vector(rF4_e_local(:, 1), rF4_e_local(:, 3)))
		diff_e(2, 1)= 0.
		cos_tmp = cos2vector(rF4_e_local(:,3)-rF4_e_local(:,1), &
			                 rF4_e_local(:,4)-rF4_e_local(:,2), sin_tmp)
		diff_e(:, 2)= acos(cos2vector(rF4_e_local(:, 2), rF4_e_local(:, 4)))*[cos_tmp, sin_tmp]
		Bn_e=brn_e(i,j)
	else if (.not. ((key_trace4 .or. key_diff) .and. local_launch)) then
		diff_e=NaN
	endif

	Dmatrix(1,1)=  diff_e(1, 1)*diff_s(2, 2) - diff_e(1, 2)*diff_s(2, 1)
	Dmatrix(1,2)= -diff_e(1, 1)*diff_s(1, 2) + diff_e(1, 2)*diff_s(1, 1)
	Dmatrix(2,1)=  diff_e(2, 1)*diff_s(2, 2) - diff_e(2, 2)*diff_s(2, 1)
	Dmatrix(2,2)= -diff_e(2, 1)*diff_s(1, 2) + diff_e(2, 2)*diff_s(1, 1)

	q_local(i,j)=abs((g_e(1)*g_s(2)*Dmatrix(1,1)**2 +  &
		              g_e(1)*g_s(1)*Dmatrix(1,2)**2 +  &
		              g_e(2)*g_s(2)*Dmatrix(2,1)**2 +  &
		              g_e(2)*g_s(1)*Dmatrix(2,2)**2) * &
				      Bn_s*Bn_e/(bn_square* gh* (4.*delta_diff(1)*delta_diff(2))**2.))

	! if q_local(i,j) is NaN
	if (.not. (q_local(i,j) .ge. 0.) .and. local_launch) q_local(i,j)=info%q_local
	
endif

end subroutine q_diff


subroutine compute_layer
implicit none
integer:: i, j, ip, label, label0, loop_end
!------------------------------------------------------------
if (sign2dFlag) allocate(sign2d(0:iend, 0:jend))
!------------------------------------------------------------
! if sflag, jend can be 0 in some case, so 'DO i= 0, iend' should be outside
!$OMP PARALLEL DO PRIVATE(i,j), schedule(DYNAMIC)
DO i= 0, iend
DO j= 0, jend
	call q_bridge(i, j)
enddo
enddo
!$OMP END PARALLEL DO
!------------------------------------------------------------
if (pole_j0 .or. pole_jend) then
do j= 0, jend, jend
	if (j .eq.    0 .and. .not. pole_j0  ) cycle
	if (j .eq. jend .and. .not. pole_jend) cycle
	if (allocate_path) then
		label0= j*nq1
		index_seed(label0+1:label0+iend)= index_seed(label0)
		loop_size(label0+1:label0+iend)= loop_size(label0)
		loop_end= loop_size(label0)-1
	endif

	DO i= 1, iend
		rFs(:, i, j)=rFs(:, 0, j)
		rFe(:, i, j)=rFe(:, 0, j)
		rboundary(i, j)=rboundary(0, j)
		b_layer(:, i, j)=b_layer(:, 0, j)
		if (privateFlag) then
			do ip =0, 9 
				if (int_private_out(ip)) &
				int_private(ip)%data(i, j)=int_private(ip)%data(0, j)
			enddo
		endif
		if (CurlB_out) CurlB_layer(:, i, j)=CurlB_layer(:, 0, j)
		
		if (targetB_flag) then 
			bs_layer(:, i, j)=bs_layer(:, 0, j)
			be_layer(:, i, j)=be_layer(:, 0, j)
		endif 

		q(i,j)=q(0, j)
		if (q_local_Flag) q_local(i,j)=q_local(0, j)
		if (scottFlag) then	
			q_perp(i,j)=q_perp(0, j)
		else if (diff_flag) then
			rbs(i, j)=rbs(0, j)
			rbe(i, j)=rbe(0, j)
			tangent(i, j)=tangent(0, j)
			bnp2d(i, j)= bnp2d(0, j)
			if (south_pole .or. north_pole) then
				rFs_yin(:, i, j)=rFs_yin(:, 0, j)
				rFe_yin(:, i, j)=rFe_yin(:, 0, j)
				s_yinFlag(i,j)=s_yinFlag(0,j)
				e_yinFlag(i,j)=e_yinFlag(0,j)
			endif
		endif

		if (targetCurlB_out) then 
			CurlBs_layer(:, i, j)=CurlBs_layer(:, 0, j)
			CurlBe_layer(:, i, j)=CurlBe_layer(:, 0, j)
		endif

		if (sign2dFlag) sign2d(i, j)= sign2d(0, j)

		if (allocate_path) then
			label= i+label0
			allocate(lines(label)%path(0:2, 0:loop_end))
			lines(label)%path = lines(label0)%path
			if (loopB_out) then
				allocate(lines(label)%loopB(0:2, 0:loop_end))
				lines(label)%loopB = lines(label0)%loopB 
			endif
			if (loopCurlB_out) then
				allocate(lines(label)%loopCurlB(0:2, 0:loop_end))
				lines(label)%loopB = lines(label0)%loopB
			endif
		endif
	enddo
enddo
endif
!------------------------------------------------------------
if (sign2dFlag) then
	open(1, file='sign2d.bin', access='stream', status='replace')
	write(1) sign2d
	close(1)
	deallocate(sign2d)
endif
!------------------------------------------------------------
if (diff_flag) then
!$OMP PARALLEL DO collapse(2) PRIVATE(i,j), schedule(DYNAMIC)
DO j= 0, jend
DO i= 0, iend
	call q_diff(i, j)
enddo
enddo
!$OMP END PARALLEL DO
endif
!------------------------------------------------------------
where(q .lt. 2.) q=2.
if (scottFlag) where(q_perp .lt. 2.) q_perp=2.
if (q_local_Flag) where(q_local .lt. 2.) q_local=2.

end subroutine compute_layer


function seed_ev3(i, j)
implicit none
integer:: i, j
real:: arrow_i(0:2), arrow_j(0:2), p0(0:2), p1(0:2), seed_ev3(0:2)
!------------------------------------------------------------
if (i .eq. 0) then
	if (spherical) then
		p0 = vp_spherical2car(seed(:,i,j))
	else
		p0 = seed(:,i,j)
	endif
else
	if (spherical) then
		p0 = vp_spherical2car(seed(:,i-1,j))
	else
		p0 = seed(:,i-1,j)
	endif
endif
if (i .eq. iend) then
	if (spherical) then
		p1 = vp_spherical2car(seed(:,i,j))
	else
		p1 = seed(:,i,j)
	endif
else
	if (spherical) then
		p1 = vp_spherical2car(seed(:,i+1,j))
	else
		p1 = seed(:,i+1,j)
	endif
endif

arrow_i=p1-p0

if (j .eq. 0) then
	if (spherical) then
		p0 = vp_spherical2car(seed(:,i,j))
	else
		p0 = seed(:,i,j)
	endif
else
	if (spherical) then
		p0 = vp_spherical2car(seed(:,i,j-1))
	else
		p0 = seed(:,i,j-1)
	endif
endif
if (j .eq. jend) then
	if (spherical) then
		p1 = vp_spherical2car(seed(:,i,j))
	else
		p1 = seed(:,i,j)
	endif
else
	if (spherical) then
		p1 = vp_spherical2car(seed(:,i,j+1))
	else
		p1 = seed(:,i,j+1)
	endif
endif
arrow_j=p1-p0

seed_ev3 = normalize_cross_product(arrow_i, arrow_j)

end function seed_ev3


function cos2vector(vector, vector1, sin2vector)
implicit none
real:: vector(0:2), vector1(0:2), cos2vector
real, optional:: sin2vector
!------------------------------------------------------------
cos2vector = dot_product(vector, vector1)/(norm2s(vector)*norm2s(vector1))
if (cos2vector .le. -1.) cos2vector=-1.
if (cos2vector .ge.  1.) cos2vector= 1.
if (present(sin2vector)) sin2vector= sqrt(1.-cos2vector**2.)
end function cos2vector


subroutine initialize_region
implicit none
integer:: i, j, nqx, nqy, nqz
real:: ev1(0:2), ev2(0:2), i_shift, j_shift, &
point0(0:2), point1(0:2), point2(0:2), point_orig(0:2), &
p0(0:2), p1(0:2), p2_spherical(0:2), arc, sin_arc, max_da(0:2)
real, allocatable:: axis1(:,:)
logical:: csFlag, preset_xreg, preset_yreg
!------------------------------------------------------------
! receive information from fastqsl.pro/fastqsl.py
open(1, file='head_region.bin', access='stream', status='old')
read(1) xreg, yreg, zreg, deltas, csFlag, preset_xreg, preset_yreg
close(1, status='delete')
!------------------------------------------------------------
if (period_lon .and. .not. preset_xreg) xreg(1)= two_pi

if (csflag) then
	normal_index = -1

	delta_i=deltas(-1)
	delta_j=deltas( 2)

	if ((-xreg(0) .eq. xreg(1)) .and. (-yreg(0) .eq. yreg(1))) then
		if (spherical) then
			p1= [cos(yreg(1))*[cos(xreg(1)), sin(xreg(1))], sin(yreg(1))]
			arc= acos(cos2vector([1., 0., 0.], p1, sin_arc))
			nq1= int(arc/delta_i)*2+1
		else
			nq1= int(norm2s([xreg(1), yreg(1)])/delta_i)*2+1
		endif	
	else
		if (spherical) then
			p0= [cos(yreg(0))*[cos(xreg(0)), sin(xreg(0))], sin(yreg(0))]
			p1= [cos(yreg(1))*[cos(xreg(1)), sin(xreg(1))], sin(yreg(1))]
			arc= acos(cos2vector(p0, p1, sin_arc))
			nq1= int(arc/delta_i)+1
		else
			nq1= int(norm2s([xreg(1)-xreg(0), yreg(1)-yreg(0)])/delta_i)+1
		endif
	endif
	
	nq2=int(abs(zreg(1)-zreg(0))/delta_j)+1
else
	if (-xreg(0) .eq. xreg(1)) then
		! make the symmetry around x=0
		nqx= int(xreg(1)/deltas(0))*2+1
	else
		nqx= int((xreg(1)-xreg(0))/deltas(0))+1
	endif

	if (spherical) then
		if (south_pole .and. .not. preset_yreg) yreg(0)= -half_pi
		if (north_pole .and. .not. preset_yreg) yreg(1)=  half_pi
	endif
	if (-yreg(0) .eq. yreg(1)) then
		! make the symmetry around y=0
		nqy= int(yreg(1)/deltas(1))*2+1
	else
		nqy= int((yreg(1)-yreg(0))/deltas(1))+1
	endif
	
	nqz= int((zreg(1)-zreg(0))/deltas(2))+1
	
	if (xreg(1) .eq. xreg(0)) then
		normal_index=0; nq1=nqy; nq2=nqz; delta_i=deltas(1); delta_j=deltas(2)
	else if (yreg(1) .eq. yreg(0)) then
		normal_index=1; nq1=nqx; nq2=nqz; delta_i=deltas(0); delta_j=deltas(2)
	else 
		normal_index=2; nq1=nqx; nq2=nqy; delta_i=deltas(0); delta_j=deltas(1)
	endif
endif

if (vflag) then
	nq3=nqz
else
	nq3=1
endif

iend=nq1-1
jend=nq2-1
!------------------------------------------------------------
! 2D grids for calculation
! ev1, ev2: two cartesian elementary vectors of the cross section
! ev3: the cartesian normal direction of the cross section

allocate(seed(0:2, 0:iend, 0:jend))

if (csFlag .and. spherical) then
	if ((-xreg(0) .eq. xreg(1)) .and. (-yreg(0) .eq. yreg(1))) then
		seed(0:1, iend/2, 0)=0.
		do i=iend/2+1, iend
			p2_spherical = vp_car2spherical([1.,0.,0.]*sin(arc-delta_i*(i-iend/2))+ &
			                                        p1*sin(delta_i*(i-iend/2))/sin_arc)
			seed(0:1, i, 0)      =  p2_spherical(0:1)
			seed(0:1, iend-i, 0) = -p2_spherical(0:1)
		enddo
	else
		do i=0, iend
			p2_spherical = vp_car2spherical((p0*sin(arc-delta_i*i)+p1*sin(delta_i*i))/sin_arc)
			seed(0:1, i, 0) = p2_spherical(0:1)
		enddo
	endif

	forall(j=1:jend) seed(0:1, :, j) = seed(0:1, :, 0)
	forall(j=0:jend) seed(  2, :, j) = zreg(0) + delta_j*j

	ev3=normalize_cross_product(p0, p1)
else
	point0=[xreg(0),yreg(0),zreg(0)]

	! default values, can be adjusted if a sysmetry is found
	point_orig=point0
	i_shift=0
	j_shift=0

	select case(Normal_index)
	case(-1)
		point1=[xreg(1),yreg(1),zreg(0)]
		point2=[xreg(0),yreg(0),zreg(1)]

		ev1=normalize(point1-point0)
		ev2=normalize(point2-point0)
		
		ev3=normalize_cross_product(ev1, ev2)
		if (xreg(0) .eq. -xreg(1) .and. yreg(0) .eq. -yreg(1)) then
			i_shift=iend/2
			point_orig(0:1)=0.
		endif
	case(0)
		ev1=[0.,1.,0.]; ev2=[0.,0.,1.]
		if (yreg(0) .eq. -yreg(1)) then
			i_shift=iend/2
			point_orig(1)=0.
		endif
	case(1)
		ev1=[1.,0.,0.]; ev2=[0.,0.,1.]
		if (xreg(0) .eq. -xreg(1)) then
			i_shift=iend/2
			point_orig(0)=0.
		endif
	case(2)
		ev1=[1.,0.,0.]; ev2=[0.,1.,0.]
		if (xreg(0) .eq. -xreg(1)) then
			i_shift=iend/2
			point_orig(0)=0.
		endif
		if (yreg(0) .eq. -yreg(1)) then
			j_shift=jend/2
			point_orig(1)=0.
		endif
	end select
	forall(i=0:iend, j=0:jend) seed(:,i,j) = point_orig + (i-i_shift)*delta_i*ev1 + (j-j_shift)*delta_j*ev2
endif

if (csFlag) then
	allocate(axis1(0:1, 0:iend))
	axis1=seed(0:1, :, 0)
	open(1, file='axis1.bin', access='stream', status='replace')
	write(1) axis1
	close(1)
	deallocate(axis1)
endif
!------------------------------------------------------------
! inform fastqsl.pro/fastqsl.py these; 
! xreg, yreg, zreg(1) may be changed due to flooring of giving nq1, nq2, nq3
xreg=[seed(0, 0, 0), seed(0, iend, jend)]
yreg=[seed(1, 0, 0), seed(1, iend, jend)]
if (vflag) then
	zreg(1)=zreg(0)+deltas(2)*(nq3-1)
else
	zreg(1)=seed(2, iend, jend)
endif

open(1, file='tail_region.bin', access='stream', status='replace')
write(1) nq1, nq2, nq3, normal_index, xreg, yreg, zreg
close(1)
!------------------------------------------------------------
pole_j0  = south_pole .and. normal_index .eq. 2 .and. (abs(yreg(0)+half_pi) .lt. deltas(1)/32.)
pole_jend= north_pole .and. normal_index .eq. 2 .and. (abs(half_pi-yreg(1)) .lt. deltas(1)/32.)
!------------------------------------------------------------
if (stretchFlag) then
	forall(i=0:2) max_da(i)=maxval(axis(i)%da)
	if (spherical) then
		select case (normal_index)
		case (-1)
			min_step=minval([deltas(-1)/maxval(max_da(0:1)), deltas(2)/max_da(2)])
		case (0)
			min_step=minval(deltas(1:2)  /max_da(1:2))
		case (1)
			min_step=minval(deltas(0:2:2)/max_da(0:2:2))
		case (2)
			min_step=minval(deltas(0:1)  /max_da(0:1))
		end select
	else
		min_step=delta_i/maxval(max_da)
	endif
else
	min_step=delta_i
endif
min_step=min_step/2.

end subroutine initialize_region

end module compute


program fastqsl
use compute
implicit none
logical:: qflag, verbose, launch_out
real(8):: tcalc, tnow, tend, omp_get_wtime
integer:: i, k, nthreads, i2end, OMP_GET_NUM_PROCS
integer(8), allocatable:: indexes(:)
integer(1):: ip
character(len=1) :: str_aux
!------------------------------------------------------------
open(1, file='head.bin', access='stream', status='old')
read(1) step, tol, r_local, maxsteps, RK4Flag, inclineFlag, &
        launch_out, B_out, CurlB_out, &
		rF_out, targetB_out, targetCurlB_out, &
		path_out, loopB_out, loopCurlB_out, &
		sflag, bflag, cflag, vflag, nthreads, scottFlag, &
		verbose, keep_tmp, magnetogram_out, int_private_out
close(1, status='delete')
!------------------------------------------------------------
! if (verbose) call system_clock(tnow)
if (verbose) tnow=omp_get_wtime()
NaN = transfer(2143289344, 1.0)
pi = 3.141592653589793
half_pi=pi/2.
two_pi =pi*2.
!------------------------------------------------------------
if (sflag) then
	open(1, file='dim_seed.bin', access='stream', status='old')
	read(1) nq1, nq2, nq3
	close(1, status='delete')

	launch_out=.false.

	iend=nq1-1
	jend=nq2-1
	allocate(seed(0:2, 0:iend, 0:jend))
	
	pole_j0  =.false.
	pole_jend=.false.

	! if Q is calculated by Scott (2017), min_step=0.125 is good enough
	! If Q is calculated by Pariat (2012), min_step \times grid spacing of {xyz}a
	! should be smaller than the grid spacing of seed; usually min_step=0.125 is good enough
	min_step=0.125

	Normal_index = -1
endif
!------------------------------------------------------------
! https://www.openmp.org/spec-html/5.0/openmpsu112.html
if (nthreads .gt. OMP_GET_NUM_PROCS()) nthreads=OMP_GET_NUM_PROCS()
if (nthreads .eq. 0) nthreads=OMP_GET_NUM_PROCS()-2
CALL OMP_set_num_threads(nthreads)

privateFlag=any(int_private_out)
traceflag = maxsteps .ne. 0

call set_CurlBvec_Flag

dbdc_field_Flag  = traceflag .and. .not. (sflag .and. scottFlag .and. nq1 .le. 1000)

call readB

if (spherical) then
	interpolate => interpolate_spherical
else
	interpolate => interpolate_cartesian
endif

if (.not. sflag) call initialize_region

diff_flag= iend .ge. 2 .and. jend .ge. 2 .and. .not. scottFlag
diff_seed= diff_flag .and. sflag
!------------------------------------------------------------
! parameters for tracing
inclineFlag = inclineFlag .and. diff_flag
min_incline = 0.05

! max steps for subroutine correct_foot
if (RK4flag) then
	if (step .lt. min_step) min_step= step
	min_step_foot= min_step*0.25
	maxsteps_foot= step/min_step_foot*4
else
	min_step_foot= min_step*0.25
	maxsteps_foot= min_step/min_step_foot*4
	max_step= minval(pend)/4.

	a21=   1./4.
	a31=   3./32.;   a32=    9./32.
	a41=1932./2197.; a42=-7200./2197.; a43=  7296./2197.
	a51= 439./216.;  a52=-8.;          a53=  3680./513.;   a54= -845./4104.
	a61=  -8./27.;   a62= 2.;          a63= -3544./2565.;  a64= 1859./4104.; a65=-11./40.
	b1 =  16./135.;   b3= 6656./12825.; b4= 28561./56430.;  b5=   -9./50.;    b6=  2./55. 
	ce1=   1./360.;  ce3= -128./4275.; ce4= -2197./75240.; ce5=    1./50.;   ce6=  2./55.	
endif
!------------------------------------------------------------
ijend=[iend, jend]
! allocate arrays for module compute
allocate(q(0:iend, 0:jend))
allocate(rboundary(0:iend, 0:jend))
allocate(rFs(0:2, 0:iend, 0:jend))
allocate(rFe(0:2, 0:iend, 0:jend))
allocate(b_layer(0:2, 0:iend, 0:jend))

q_local_Flag= r_local .gt. 0. .and. (diff_flag .or. scottFlag)
if (q_local_Flag) then
	r_local_square=r_local**2.
	allocate(q_local(0:iend, 0:jend))
	if (diff_flag) then 
		allocate(local_s_flag(0:iend, 0:jend))
		allocate(local_e_flag(0:iend, 0:jend))
		allocate(brn_s(0:iend, 0:jend))
		allocate(brn_e(0:iend, 0:jend))
	endif
endif

targetB_flag= diff_flag .or. targetB_out
if (targetB_flag) then 
	allocate(bs_layer(0:2, 0:iend, 0:jend))
	allocate(be_layer(0:2, 0:iend, 0:jend))
endif

if (privateFlag) then
	do ip = 0, 9 
		if (int_private_out(ip)) &
		allocate(int_private(ip)%data(0:iend, 0:jend))
	enddo
endif

if (scottFlag) allocate(q_perp(0:iend, 0:jend))
if (diff_flag) then
	allocate(rbs(0:iend, 0:jend))
	allocate(rbe(0:iend, 0:jend))
	allocate(tangent(0:iend, 0:jend))
	allocate(bnp2d(0:iend, 0:jend))
	allocate(s_yinFlag(0:iend, 0:jend))
	allocate(e_yinFlag(0:iend, 0:jend))
	if (south_pole .or. north_pole) then
		allocate(rFs_yin(0:2, 0:iend, 0:jend))
		allocate(rFe_yin(0:2, 0:iend, 0:jend))
	endif
endif
if (CurlB_out) allocate(CurlB_layer(0:2, 0:iend, 0:jend))
if (targetCurlB_out) then
	allocate(CurlBs_layer(0:2, 0:iend, 0:jend))
	allocate(CurlBe_layer(0:2, 0:iend, 0:jend))
endif

allocate_path = path_out .or. (q_local_Flag .and. diff_flag)
if (allocate_path) then
	i2end=nq1*nq2-1
	allocate(lines(0:i2end))
	allocate(loop_size(0:i2end))
	allocate(index_seed(0:i2end))
	allocate(indexes(0:i2end))
	indexes(0)=0
endif
!------------------------------------------------------------
! announce the start of the computation
if (verbose) then
	print*, '  _____________________________________'
	print*, '        schedule         time'
endif
!------------------------------------------------------------
! Fortran use unit 0 for error, unit 5 for input (keyboard) and unit 6 for output (screen), these units should not be used
! unit 1 is used in compute_layer, if traceflag is .false., scottFlag, int_private_out, rF_out 
! are already set to .false. in fastqsl.pro/fastqsl.py

qflag= traceflag .and. (diff_flag .or. scottFlag)

if (traceflag)       open(2,  file='rboundary.bin', access='stream', status='replace')
if (qflag)           open(3,  file='q.bin',         access='stream', status='replace')
if (scottFlag)       open(4,  file='q_perp.bin',    access='stream', status='replace')
if (rF_out)          open(9,  file='rFs.bin',       access='stream', status='replace')
if (rF_out)          open(10, file='rFe.bin',       access='stream', status='replace')
if (b_out)           open(11, file='B.bin',         access='stream', status='replace')
if (CurlB_out)       open(12, file='CurlB.bin',     access='stream', status='replace')
if (targetB_out)     open(13, file='Bs.bin',        access='stream', status='replace')
if (targetB_out)     open(14, file='Be.bin',        access='stream', status='replace')
if (targetCurlB_out) open(15, file='CurlBs.bin',    access='stream', status='replace')
if (targetCurlB_out) open(16, file='CurlBe.bin',    access='stream', status='replace')

! launch_out, sflag can not both be .true. at any case
if (launch_out)      open(17, file='seed.bin',      access='stream', status='replace')
if (sflag)           open(17, file='seed.bin',      access='stream', status='old')
if (q_local_Flag)    open(18, file='q_local.bin',   access='stream', status='replace')

if (path_out) then
                     open(20, file='index_seed.bin',access='stream', status='replace')
                     open(21, file='indexes.bin',   access='stream', status='replace')
                     open(22, file='path.bin',      access='stream', status='replace')
if (loopB_out)       open(23, file='loopB.bin',     access='stream', status='replace')
if (loopCurlB_out)   open(24, file='loopCurlB.bin', access='stream', status='replace')
endif

do ip=0, 9
	if (int_private_out(ip)) then
		write(str_aux,"(i0)") ip
		open(30+ip, file='int_private'//trim(str_aux)//'.bin', access='stream', status='replace')
	endif
enddo

do k=0, nq3-1
	if (verbose .and. mod(k, nthreads) .eq. 0) call show_time(float(k)/nq3*100.0)

	if (vflag) then
		seed(2, :, :) = zreg(0) + k*deltas(2)
		sign2dFlag= seed(2,0,0) .eq. pmin(2) .and. traceFlag
	else if (sflag) then
		read(17) seed
		sign2dFlag= all(seed(2,:,:) .eq. pmin(2)) .and. traceFlag .and. nq1 .ge. 2 .and. nq2 .ge. 2
	else
		sign2dFlag= bflag .and. traceFlag
	endif

	call compute_layer

	if (traceflag)       write(2)  rboundary
	if (qflag)           write(3)  q
	if (scottFlag)       write(4)  q_perp
	if (rF_out)          write(9)  rFs
	if (rF_out)          write(10) rFe
	if (b_out)           write(11) b_layer
	if (CurlB_out)       write(12) CurlB_layer
	if (targetB_out)     write(13) bs_layer
	if (targetB_out)     write(14) be_layer
	if (targetCurlB_out) write(15) CurlBs_layer
	if (targetCurlB_out) write(16) CurlBe_layer
	if (launch_out)      write(17) seed
    if (q_local_Flag)    write(18) q_local

	if (path_out) then
		write(20) index_seed

		do i=1, i2end
			indexes(i)=indexes(i-1)+loop_size(i-1)
		enddo
		write(21) indexes
		indexes(0)=indexes(i2end)+loop_size(i2end)

		do i=0, i2end
			write(22) lines(i)%path
		enddo
		if (loopB_out) then
			do i=0, i2end
				write(23) lines(i)%loopB
			enddo
		endif
		if (loopCurlB_out) then
			do i=0, i2end
				write(24) lines(i)%loopCurlB
			enddo
		endif
	endif

	if (allocate_path) then
		do i=0, i2end
			deallocate(lines(i)%path)
			if (loopB_out) deallocate(lines(i)%loopB)
			if (loopCurlB_out) deallocate(lines(i)%loopCurlB)
		enddo
	endif
	
	do ip=0, 9
		if (int_private_out(ip)) &
		write(30+ip) int_private(ip)%data
	enddo
	
enddo

if (traceflag)             close(2)
if (qflag)                 close(3)
if (scottFlag)             close(4)
if (rF_out)                close(9)
if (rF_out)                close(10)
if (b_out)                 close(11)
if (CurlB_out)             close(12)
if (targetB_out)           close(13)
if (targetB_out)           close(14)
if (targetCurlB_out)       close(15)
if (targetCurlB_out)       close(16)
if (launch_out .or. sflag) close(17)
if (q_local_Flag)          close(18) 

if (path_out) then
                           close(20)
write(21) indexes(0)
                           close(21)
                           close(22)
if (loopB_out)             close(23)
if (loopCurlB_out)         close(24)
endif

do ip=0, 9
	if (int_private_out(ip)) close(30+ip)
enddo
!------------------------------------------------------------
! house keeping
deallocate(q, rboundary, rFs, rFe, b_layer, seed, Bvec)
if (allocate_path) deallocate(lines, loop_size, indexes, index_seed)
if (CurlB_out) deallocate(CurlB_layer)
if (targetB_flag) deallocate(bs_layer, be_layer)
if (targetCurlB_out) deallocate(CurlBs_layer, CurlBe_layer)
if (scottFlag) deallocate(q_perp)
if (diff_flag) deallocate(rbs, rbe, tangent, bnp2d, s_yinFlag, e_yinFlag)
if (diff_flag .and. (south_pole .or. north_pole)) deallocate(rFs_yin, rFe_yin)

if (q_local_Flag) then
	deallocate(q_local)
	if (diff_flag) deallocate(local_s_flag, local_e_flag, brn_s, brn_e)
endif

if (privateFlag) then
	do ip = 0, 9 
		if (int_private_out(ip)) &
		deallocate(int_private(ip)%data)
	enddo
endif

if ( dbdc_field_Flag) deallocate( dbdc_field)
if (CurlBvec_Flag) deallocate(CurlBvec)
if (A_input)       deallocate(Avec)
if (stretchFlag) then
	do i=0, 2
		deallocate(axis(i)%pa, axis(i)%da, axis(i)%coef_diff)
	enddo
	if (binary_index_top .ne. 0) deallocate(binary_values)
endif
round_weight => null()
interpolate  => null()
!------------------------------------------------------------
if (verbose) then
	! when 100.00% is printed, everything is done in fastqsl.x
	call show_time(100.0)

	tend=omp_get_wtime()
	
	tcalc=tend-tnow !s
	if (tcalc .ge. 3600.) then
		print '(F7.2,   " hours elapsed in fastqsl.x")', tcalc/3600.
	else if (tcalc .ge. 60.) then
		print '(F7.2, " minutes elapsed in fastqsl.x")', tcalc/60.
	else
		print '(F7.2, " seconds elapsed in fastqsl.x")', tcalc
	endif
endif
!------------------------------------------------------------
! In Windows, the pop-up window for fastqsl.exe can not be closed automatically
! call get_environment_variable("HOME", str_aux)
! if (str_aux .ne. "/") &
! call system('taskkill /im fastqsl.exe /f') 

! another way to kill the pop-up window
! call abort
end program fastqsl