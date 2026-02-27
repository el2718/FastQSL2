include 'fields.f90'
include 'trace_bline.f90'

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
integer(1), target, allocatable:: rbs(:, :), rbe(:, :)
integer(1), allocatable:: rboundary(:, :)
integer(2), allocatable:: sign2d(:, :)! IDL do not have the type of 8-bit signed integer
real:: deltas(-1:2), delta_i, delta_j, xreg(0:1), yreg(0:1), zreg(0:1), ev3(0:2)
real, allocatable:: q(:, :), q_perp(:, :), length(:, :), twist(:, :), &
seed(:, :, :), b_layer(:, :, :), curlb_layer(:, :, :), bnp2d(:, :), &
bs_layer(:, :, :), be_layer(:, :, :), curlbs_layer(:, :, :), curlbe_layer(:, :, :), &
q_local(:, :), brn_s(:,:), brn_e(:,:)
real, allocatable, target:: rFs(:, :, :), rFe(:, :, :), rFs_yin(:, :, :), rFe_yin(:, :, :)
logical:: vflag, bflag, cflag, sFlag, scottFlag, diff_seed, pole_j0, pole_jend, &
targetB_flag, sign2dFlag, allocate_path, path_out, loopB_out, loopCurlB_out, q_local_Flag
logical, allocatable:: tangent(:, :), local_s_flag(:, :), local_e_flag(:, :)
logical, allocatable, target:: s_yinFlag(:, :), e_yinFlag(:, :)
type line
	real, allocatable:: path(:, :), loopB(:, :), loopCurlB(:, :)
endtype line
type(line), allocatable :: lines(:)

integer, allocatable:: index_seed(:), loop_size(:)

contains

subroutine q_bridge(i, j)
implicit none
integer:: i, j, ip, its, ite
logical:: pole_j
type(line_info):: info
character(len=4) ::i_str, j_str
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
	ip=i+j*nq1
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
	index_seed(ip)= -its
	loop_size(ip)= ite-its+1
	allocate(lines(ip)%path(0:2, 0:ite-its))
	lines(ip)%path = info%path(:, its:ite)
	deallocate(info%path)
	if (info%loopB_out) then 
		allocate(lines(ip)%loopB(0:2, 0:ite-its))
		lines(ip)%loopB = info%loopB(:, its:ite)
		deallocate(info%loopB)
	endif
	if (info%loopCurlB_out) then
		allocate(lines(ip)%loopCurlB(0:2, 0:ite-its))
		lines(ip)%loopCurlB = info%loopCurlB(:, its:ite)
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

if (length_out) length(i, j)=info%length
if (twist_out) twist(i, j)=info%twist
if (curlb_out) curlb_layer(:, i, j)=info%curlBp

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
i_dim, i_diff, ip4(1:4), ij(1:2), diff_index(0:2,1:2)
integer, pointer:: dims(:)
integer, target :: dims_s(1:2), dims_e(1:2)
real:: delta_diff(1:2), bn_square, Bn_s, Bn_e, Dmatrix(1:2, 1:2),  &
rF4_s(1:2, 1:4), rF4_e(1:2, 1:4), u0(0:2), v0(0:2), &
g_e(1:2), g_s(1:2), g_e0(0:2), g_s0(0:2), gh, cos_p(0:1), sin_p(0:1), &
cos_tmp, sin_tmp, rF4_s_local(0:2, 1:4), rF4_e_local(0:2, 1:4), vr2vp(0:2), &
rf3(1:2, 0:2), seed3(0:2, 0:2), coef(0:2,1:2), d0, d1, arrow_seed(0:2, 1:2), vp_yin(0:2)
real, pointer:: bp_car(:), vp_car(:), vp4_car(:,:), rF_tmp(:, :, :), diff(:, :), rF_yin(:, :, :)
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
				coef(:,i_diff) = [-1., 0., 1.]
			endif
			coef(1,i_diff) = - coef(0,i_diff) - coef(2,i_diff)
		enddo

		if (delta_diff(1)* 100. < delta_diff(2) .or. delta_diff(2)* 100. < delta_diff(1)) key_nB=.true.

		cos_tmp = cos2vector(arrow_seed(:,1), arrow_seed(:,2), sin_tmp)
		if (sin_tmp .le. 0.01) then
			key_nB=.true.
		else
			gh=sin_tmp
		endif

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
	site  %alphaFlag=.false.
	site  %scottFlag=.false.
	site  %scottLaunch=.false.
	site_r%alphaFlag=.false.
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
		ip4= i + j*nq1 + [-1, -nq1, 1, nq1]
		do k = 1, 4
			if (distance(lines(ip4(k))%path(:, index_seed(ip4(k))), vp) .gt. r_local) then
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
					it_end= loop_size(ip4(k))-1
				endif

				do it = index_seed(ip4(k))+sign_dt, it_end, sign_dt

					exist_vr= distance(lines(ip4(k))%path(:, it), vp) .ge. r_local
					if (exist_vr) then
						site  %v(0:2)= lines(ip4(k))%path(:, it-sign_dt)
						site_r%v(0:2)= lines(ip4(k))%path(:, it)
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

	q(i,j)=abs((g_e(1)*g_s(2)*Dmatrix(1,1)**2. +  &
		        g_e(1)*g_s(1)*Dmatrix(1,2)**2. +  &
		        g_e(2)*g_s(2)*Dmatrix(2,1)**2. +  &
		        g_e(2)*g_s(1)*Dmatrix(2,2)**2.) * &
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

	q_local(i,j)=abs((g_e(1)*g_s(2)*Dmatrix(1,1)**2. +  &
		              g_e(1)*g_s(1)*Dmatrix(1,2)**2. +  &
		              g_e(2)*g_s(2)*Dmatrix(2,1)**2. +  &
		              g_e(2)*g_s(1)*Dmatrix(2,2)**2.) * &
				      Bn_s*Bn_e/(bn_square* gh* (4.*delta_diff(1)*delta_diff(2))**2.))

	! if q_local(i,j) is NaN
	if (.not. (q_local(i,j) .ge. 0.) .and. local_launch) q_local(i,j)=info%q_local
	
endif

end subroutine q_diff


subroutine compute_layer
implicit none
integer:: i, j, label, label0, loop_end
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
		if (length_out) length(i, j)=length(0, j)
		if (twist_out) twist(i, j)=twist(0, j)
		if (curlb_out) curlb_layer(:, i, j)=curlb_layer(:, 0, j)
		
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
!------------------------------------------------------------
if (sign2dFlag) then
	open(1, file='sign2d.bin', access='stream', status='replace')
	write(1) sign2d
	close(1)
	deallocate(sign2d)
endif
!------------------------------------------------------------
if (diff_flag) then
!$OMP PARALLEL DO PRIVATE(i,j), schedule(DYNAMIC)
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
real:: ev1(0:2), ev2(0:2), point0(0:2), point1(0:2), point2(0:2), &
p0(0:2), p1(0:2), p2_spherical(0:2), arc, sin_arc, max_da(0:2)
real, allocatable:: axis1(:, :)
logical:: csFlag, preset_xreg, preset_yreg
!------------------------------------------------------------
! receive information from fastqsl.pro
open(1, file='head_region.bin', access='stream', status='old')
read(1) xreg, yreg, zreg, deltas, csFlag, preset_xreg, preset_yreg
close(1, status='delete')
!------------------------------------------------------------
if (period_lon .and. .not. preset_xreg) xreg(1)= two_pi

if (csflag) then
	normal_index = -1

	delta_i=deltas(-1)
	delta_j=deltas( 2)

	if (spherical) then
		p0= [cos(yreg(0))*[cos(xreg(0)), sin(xreg(0))], sin(yreg(0))]
		p1= [cos(yreg(1))*[cos(xreg(1)), sin(xreg(1))], sin(yreg(1))]
		arc= acos(cos2vector(p0, p1, sin_arc))
		nq1= int(arc/delta_i)+1
	else
		nq1= int(norm2s([xreg(1)-xreg(0), yreg(1)-yreg(0)])/delta_i)+1
	endif
	
	nq2=int(abs(zreg(1)-zreg(0))/delta_j)+1
else
	nqx= int((xreg(1)-xreg(0))/deltas(0))+1
	nqz= int((zreg(1)-zreg(0))/deltas(2))+1

	if (spherical) then
		if (south_pole .and. .not. preset_yreg) yreg(0)= -half_pi
		if (north_pole .and. .not. preset_yreg) yreg(1)=  half_pi
		if (-yreg(0) .eq. yreg(1)) then
			! make the symmetry around latitude=0
			nqy= int(yreg(1)/deltas(1))*2+1
			yreg(0)=-(nqy-1)/2*deltas(1)
		else
			nqy= int((yreg(1)-yreg(0))/deltas(1))+1
		endif
	else
		nqy= int((yreg(1)-yreg(0))/deltas(1))+1
	endif
	
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
	do i=0, iend
		p2_spherical = vp_car2spherical((p0*sin(arc-delta_i*i)+p1*sin(delta_i*i))/sin_arc)
		seed(0:1, i, 0) = p2_spherical(0:1)
	enddo

	forall(j=1:jend) seed(0:1, :, j) = seed(0:1, :, 0)
	forall(j=0:jend) seed(  2, :, j) = zreg(0) + delta_j*j

	ev3=normalize_cross_product(p0, p1)
else

	point0=[xreg(0),yreg(0),zreg(0)]

	select case(Normal_index)
	case(-1)
		point1=[xreg(1),yreg(1),zreg(0)]
		point2=[xreg(0),yreg(0),zreg(1)]

		ev1=normalize(point1-point0)
		ev2=normalize(point2-point0)
		
		ev3=normalize_cross_product(ev1, ev2)
	case(0)
		ev1=[0.,1.,0.]; ev2=[0.,0.,1.]
	case(1)
		ev1=[1.,0.,0.]; ev2=[0.,0.,1.]
	case(2)
		ev1=[1.,0.,0.]; ev2=[0.,1.,0.]
	end select

	forall(i=0:iend, j=0:jend) seed(:, i, j) = point0 + i*delta_i*ev1 + j*delta_j*ev2
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
! inform fastqsl.pro these; 
! xreg(1), yreg(1), zreg(1) may be changed due to flooring of giving nq1, nq2, nq3

xreg(1)=seed(0, iend, jend)
yreg(1)=seed(1, iend, jend)
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
real:: tcalc
integer:: i, k, nthreads, loop_end, OMP_GET_NUM_PROCS, tnow, tend
integer(8), allocatable:: indexes(:)
!------------------------------------------------------------
open(1, file='head.bin', access='stream', status='old')
read(1) step, tol, r_local, maxsteps, RK4Flag, inclineFlag, &
        launch_out, B_out, CurlB_out, length_out, twist_out, &
		rF_out, targetB_out, targetCurlB_out, &
		path_out, loopB_out, loopCurlB_out, &
		sflag, bflag, cflag, vflag, nthreads, scottFlag, verbose
close(1, status='delete')
!------------------------------------------------------------
if (verbose) call system_clock(tnow)
NaN = transfer(2143289344, 1.0)
pi = 3.141592653589793
over8pi = 1./(8.0*pi)
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

traceflag = maxsteps .ne. 0
curlB_field_Flag = twist_out .or. curlB_out .or. targetCurlB_out .or. loopCurlB_out
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
	if (step .lt. min_step) min_step=step
	min_step_foot=min_step*0.25
	maxsteps_foot=    step/min_step_foot*4
else
	step=min_step
	min_step_foot=min_step*0.25
	maxsteps_foot=min_step/min_step_foot*4
	max_step=minval(pend)/4.

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
if (length_out) allocate(length(0:iend, 0:jend))
if (twist_out)  allocate( twist(0:iend, 0:jend))
if (scottFlag) 	allocate(q_perp(0:iend, 0:jend))
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
if (curlb_out) allocate(curlb_layer(0:2, 0:iend, 0:jend))
if (targetCurlB_out) then
	allocate(curlbs_layer(0:2, 0:iend, 0:jend))
	allocate(curlbe_layer(0:2, 0:iend, 0:jend))
endif

allocate_path = path_out .or. (q_local_Flag .and. diff_flag)
if (allocate_path) then
	loop_end=nq1*nq2-1
	allocate(lines(0:loop_end))
	allocate(loop_size(0:loop_end))
	allocate(index_seed(0:loop_end))
	allocate(indexes(0:loop_end))
	indexes(0)=0
endif
!------------------------------------------------------------
! announce the start of the computation
if (verbose) then
	print*, '  _____________________________________'
	print*, '        schedule         time'
endif
! call system_clock(tnow)
!------------------------------------------------------------
! Fortran use unit 0 for error, unit 5 for input (keyboard) and unit 6 for output (screen), these units should not be used
! unit 1 is used in compute_layer
! if traceflag is .false., scottFlag, twist_out, length_out, rF_out are already set to .false. in fastqsl.pro

qflag= traceflag .and. (diff_flag .or. scottFlag)

if (traceflag)       open(2,  file='rboundary.bin', access='stream', status='replace')
if (qflag)           open(3,  file='q.bin',         access='stream', status='replace')
if (scottFlag)       open(4,  file='q_perp.bin',    access='stream', status='replace')
if (twist_out)       open(7,  file='twist.bin',     access='stream', status='replace')
if (length_out)      open(8,  file='length.bin',    access='stream', status='replace')
if (rF_out)          open(9,  file='rFs.bin',       access='stream', status='replace')
if (rF_out)          open(10, file='rFe.bin',       access='stream', status='replace')
if (b_out)           open(11, file='B.bin',         access='stream', status='replace')
if (curlB_out)       open(12, file='CurlB.bin',     access='stream', status='replace')
if (targetB_out)     open(13, file='Bs.bin',        access='stream', status='replace')
if (targetB_out)     open(14, file='Be.bin',        access='stream', status='replace')
if (targetCurlB_out) open(15, file='CurlBs.bin',    access='stream', status='replace')
if (targetCurlB_out) open(16, file='CurlBe.bin',    access='stream', status='replace')

! launch_out, sflag can not both be .true. at any case
if (launch_out)      open(17, file='seed.bin',      access='stream', status='replace')
if (sflag)           open(17, file='seed.bin',      access='stream', status='old')
if (q_local_Flag)    open(18, file='q_local.bin',   access='stream', status='replace')

if (path_out) then
                     open(19, file='index_seed.bin',access='stream', status='replace')
                     open(20, file='indexes.bin',   access='stream', status='replace')
                     open(21, file='path.bin',      access='stream', status='replace')
if (loopB_out)       open(22, file='loopB.bin',     access='stream', status='replace')
if (loopCurlB_out)   open(23, file='loopCurlB.bin', access='stream', status='replace')
endif

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
	if (twist_out)       write(7)  twist
	if (length_out)      write(8)  length
	if (rF_out)          write(9)  rFs
	if (rF_out)          write(10) rFe
	if (b_out)           write(11) b_layer
	if (curlB_out)       write(12) curlB_layer
	if (targetB_out)     write(13) bs_layer
	if (targetB_out)     write(14) be_layer
	if (targetCurlB_out) write(15) curlBs_layer
	if (targetCurlB_out) write(16) curlBe_layer
	if (launch_out)      write(17) seed
    if (q_local_Flag)    write(18) q_local

	if (path_out) then
		write(19) index_seed

		do i=1, loop_end
			indexes(i)=indexes(i-1)+loop_size(i-1)
		enddo
		write(20) indexes
		indexes(0)=indexes(loop_end)+loop_size(loop_end)

		do i=0, loop_end
			write(21) lines(i)%path
		enddo
		if (loopB_out) then
			do i=0, loop_end
				write(22) lines(i)%loopB
			enddo
		endif
		if (loopCurlB_out) then
			do i=0, loop_end
				write(23) lines(i)%loopCurlB
			enddo
		endif
	endif

	if (allocate_path) then
		do i=0, loop_end
			deallocate(lines(i)%path)
			if (loopB_out) deallocate(lines(i)%loopB)
			if (loopCurlB_out) deallocate(lines(i)%loopCurlB)
		enddo
	endif
enddo

if (traceflag)             close(2)
if (qflag)                 close(3)
if (scottFlag)             close(4)
if (twist_out)             close(7)
if (length_out)            close(8)
if (rF_out)                close(9)
if (rF_out)                close(10)
if (b_out)                 close(11)
if (curlB_out)             close(12)
if (targetB_out)           close(13)
if (targetB_out)           close(14)
if (targetCurlB_out)       close(15)
if (targetCurlB_out)       close(16)
if (launch_out .or. sflag) close(17)
if (q_local_Flag)          close(18)

if (path_out) then
                           close(19)
write(20) indexes(0)
                           close(20)
                           close(21)
if (loopB_out)             close(22)
if (loopCurlB_out)         close(23)
endif
!------------------------------------------------------------
! house keeping
deallocate(q, rboundary, rFs, rFe, b_layer, seed, Bfield)
if (allocate_path) deallocate(lines, loop_size, indexes, index_seed)
if (curlB_out) deallocate(curlb_layer)
if (targetB_flag) deallocate(bs_layer, be_layer)
if (targetCurlB_out) deallocate(CurlBs_layer, CurlBe_layer)
if (scottFlag) deallocate(q_perp)
if (diff_flag) deallocate(rbs, rbe, tangent, bnp2d, s_yinFlag, e_yinFlag)
if (diff_flag .and. (south_pole .or. north_pole)) deallocate(rFs_yin, rFe_yin)

if (q_local_Flag) then
	deallocate(q_local)
	if (diff_flag) deallocate(local_s_flag, local_e_flag, brn_s, brn_e)
endif
if (twist_out) deallocate(twist)
if (length_out) deallocate(length)
if ( dbdc_field_Flag) deallocate( dbdc_field)
if (curlB_field_Flag) deallocate(curlB_field)
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

	call system_clock(tend)
	tcalc=tend-tnow !ms
	if (tcalc .ge. 3.6e6) then
		print '(F7.2, " hours elapsed in fastqsl.x")', tcalc/3.6e6
	else if (tcalc .ge. 6.e4) then
		print '(F7.2, " minutes elapsed in fastqsl.x")', tcalc/6.e4
	else
		print '(F7.2, " seconds elapsed in fastqsl.x")', tcalc/1.e3
	endif
endif
!------------------------------------------------------------
! In Windows, the pop-up window for fastqsl.exe can not be closed automatically
! call system('taskkill /im fastqsl.exe /f') 

! another way to kill the pop-up window
! call abort
end program fastqsl
