module trace
use fields
implicit none
integer:: maxsteps, maxsteps_foot, normal_index
logical:: RK4flag, diff_flag, inclineFlag, traceflag, &
twist_out, length_out, B_out, CurlB_out, rf_out, targetB_out, targetCurlB_out
real:: a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, &
b1, b3, b4, b5, b6, ce1, ce3, ce4, ce5, ce6, &                   ! for RKF45
step, min_step, max_step, min_step_foot, tol, min_incline, over8pi, r_local, r_local_square

! r: remote; b:boundary, F:foot; s/e: start/end, where B vector points to inside/outside
! rbs/rbe: boundary mark for start/end point of the field line, see subroutine trim_size
! rFs/rFe: coordinates of the start/end point of the field line
! length/twist: length/twist of the field line
! bs/bp/be: B at rFs/vp/rFe
! curlBs/curlBp/curlBe: curlB at rFs/vp/rFe
! bnp, normal component of B to the surface
! tangent, the angle between B and surface .lt. arctan(min_incline)
! q/q_perp: value of q/q_perp given by Scott (2017)
! label: label of exporting fieldlines
type line_info
	real:: length, twist, rFs(0:2), rFe(0:2), ev3(0:2), bp(0:2), bnp, &
	bs(0:2), be(0:2), curlBp(0:2), curlBs(0:2), curlBe(0:2), &
	q, q_perp, q_local, brn_s, brn_e, rFs_yin(0:2), rFe_yin(0:2)
	logical:: tangent, get, scottFlag, path_out, loopB_out, loopCurlB_out, &
	q_local_Flag, local_s_flag, local_e_flag, s_yinflag, e_yinflag
	integer:: rbs, rbe, its, ite
	real, allocatable:: path(:, :), loopB(:, :), loopCurlB(:, :)
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

subroutine cal_yinyang(site, toyang, dvdsflag)
implicit none
type(site_info), target :: site
real, target:: e_yang(0:2, 0:1), e_yin(0:2, 0:1)
real, pointer:: vector1(:), vector2(:), e1(:, :), e2(:, :), dvds1(:), dvds2(:)
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

		if (site%alphaFlag) then
			forall(i=0:1) site%CurlB(i)=dot_product(site%CurlB_yin(0:1), matrix0(:, i))
			site%CurlB(2)=site%CurlB_yin(2)
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
dbdc_cell(0:2,0:2,0:1,0:1,0:1), dbdcp(0:2,0:2), da(0:2)
real, pointer:: bp(:), vector(:), dvds(:), curlBp(:), &
b3d(:,:,:,:), curlb3d(:,:,:,:), dbdc3d(:,:,:,:,:)
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
	b3d    => pole%bfield
	vector => site%v_yin
	dvds   => site%dvds_yin
	bp     => site%b_yin
	call round_weight_pole(site%v_yin(0:2), round, weight, southflag)
else
	b3d    => Bfield
	vector => site%v
	dvds   => site%dvds
	bp     => site%b
	call round_weight(site%v(0:2), round, weight)
endif
!------------------------------------------------------------
forall(i=0:2) Bp(i)=sum(weight* b3d(i, round(:,0), round(:,1), round(:,2)))
r=vector(2)
cos_lat=cos(vector(1))
dvds(0:2)= normalize(bp)/[r*cos_lat, r, 1.]
!------------------------------------------------------------
if (present(rk_first)) then
	if(site%alphaFlag) then
		if (yinflag) then
			curlB3d => pole%curlB_field
			curlbp  => site%curlb_yin
		else
			curlB3d => curlB_field
			curlbp  => site%curlb
		endif
		forall(i=0:2) CurlBp(i)=sum(weight*CurlB3d(i, round(:,0), round(:,1), round(:,2)))
		site%alpha=dot_product(CurlBp, bp)/dot_product(bp, bp)
	endif

	if (site%scottLaunch) call set_u_v(bp, vector(3:5), vector(6:8))

	! provide site%B/CurlB
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
			dbdc3d => pole%dbdc_field
		else
			dbdc3d => dbdc_field
		endif

		forall(i=0:2, j=0:2) dbdcp(i,j)=sum(weight*dbdc3d(i, j, round(:,0), round(:,1), round(:,2)))
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
real:: weight(0:1,0:1,0:1), dbdc_cell(0:2,0:2,0:1,0:1,0:1), dbdcp(0:2,0:2), da(0:2)
!------------------------------------------------------------
call round_weight(site%v(0:2), round, weight)
forall(i=0:2) site%B(i)=sum(weight* Bfield(i, round(:,0), round(:,1), round(:,2)))
site%dvds(0:2)=normalize(site%b)
!------------------------------------------------------------
if (present(rk_first)) then
	if(site%alphaFlag) then
		forall(i=0:2) site%CurlB(i)=sum(weight*curlB_field(i, round(:,0), round(:,1), round(:,2)))
		site%alpha=dot_product(site%CurlB, site%b)/dot_product(site%b, site%b)
	endif

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

site_a%alphaFlag=.false.
site_b%alphaFlag=.false.

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
logical:: alphaFlag, sum_line, identical, repeat_flag, interpolate_foot, exist_vr
integer:: i, sign_down, sign_up, sign_forward, it, sign_dt, rb, e_index, s_index
real:: vp(0:2), dt, dt_executed, step_this, tol_this, dL, int2alpha, &
Bn_s, Bn_e, us(0:2), ue(0:2), vs(0:2), ve(0:2), us1(0:2), ue1(0:2), vs1(0:2), ve1(0:2), &
bs2(0:2), be2(0:2), us2(0:2), ue2(0:2), vs2(0:2), ve2(0:2), &
b_car(0:2), cos_p(0:1), sin_p(0:1), incline, vr2vp(0:2), vr(0:2), brn
real, pointer:: bp(:), bs(:), be(:)
type(line_info), target:: info
type(site_info), target:: site_a, site_b, site_p, site_s, site_e, site_r
type(site_info), pointer:: site, site1, site_tmp
!------------------------------------------------------------
alphaFlag = info%get .and. (twist_out .or. info%loopCurlB_out)
sum_line = info%get .and. (length_out .or. alphaFlag)
if (alphaFlag) int2alpha=0.
if (length_out) info%length=0.
if (info%q_local_Flag) then
	info%local_s_flag=.false.
	info%local_e_flag=.false.
endif
!------------------------------------------------------------
site_p%v(0:2)=vp
site_p%scottFlag=info%scottFlag
site_p%scottLaunch=info%scottFlag
site_p%alphaFlag = alphaFlag .or. CurlB_out .or. (lie_boundary(vp) .and. targetCurlB_out)
site_p%yinFlag=.false. ! this would be changed by the following command
call interpolate(site_p, .false.)
bp => info%bp
bs => info%bs
be => info%be
bp =  site_p%b
if (info%get .and. curlb_out) info%curlBp=site_p%CurlB
!------------------------------------------------------------
if (.not. inside(vp)) then
	info%rbs=7; info%rFs=NaN
	info%rbe=7; info%rFe=NaN
	
	if (info%get) then
		info%length=NaN
		info%twist =NaN
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
site_b%scottLaunch=.false.
site_b%scottFlag=info%scottFlag

do sign_dt = sign_down, sign_up, 2

	site_a%alphaFlag=alphaFlag
	site_b%alphaFlag=alphaFlag
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
			if (interpolate_foot) site1%alphaFlag = alphaFlag .or. targetCurlB_out
		endif

		if (repeat_flag .or. interpolate_foot) then

			it=it+sign_dt

			call interpolate(site1, interpolate_foot)

			if (info%path_out) then
				info%path(:,it)=site1%v(0:2)
				if (info%loopB_out) info%loopB(:,it)=site1%B
				if (info%loopCurlB_out) info%loopCurlB(:,it)=site1%CurlB
			endif

			if (sum_line) then
				dL = distance(site%v(0:2), site1%v(0:2))
				if (length_out) info%length = info%length + dL
				if (alphaFlag) int2alpha = int2alpha + (site%alpha+site1%alpha)*dL
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
	if (alphaFlag) info%twist = int2alpha * over8pi
	if (targetCurlB_out) then
		info%curlBs=site_s%curlB
		info%curlBe=site_e%curlB
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
