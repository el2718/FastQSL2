module share
integer:: mode
real(8) :: pi, two_pi
logical:: r4flag, present1
real(4), allocatable :: dummy(:)

interface
	pure subroutine convert0(coor, matrix)
    implicit none
    real(8), intent(inout):: coor(0:2)
    real(8), intent(out):: matrix(0:2, 0:2)
    end subroutine convert0
end interface
procedure(convert0), pointer:: convert

contains

pure function xy2lon(xy)
implicit none
real(8), intent(in):: xy(0:1)
real(8):: xy2lon, cos_lon
!------------------------------------------------------------
if (all(xy .eq. 0.)) then 
	xy2lon=0.0D0
else
	cos_lon=xy(0)/norm2(xy)
	if (.not. (cos_lon .lt.  1.)) then
		xy2lon = 0.0D0
	else if   (cos_lon .le. -1.)  then
		xy2lon = pi
	else if   (xy(1) .ge. 0.)     then
		xy2lon =          acos(cos_lon)
	else
		xy2lon = two_pi - acos(cos_lon)
	endif
endif
end function xy2lon

! mode:
! 0: 'xyz_to_lon_lat_r'
! 1: 'lon_lat_r_to_xyz'
! 2: 'xyz_to_lon2_lat2_r'
! 3: 'lon2_lat2_r_to_xyz'
! 4: 'lon_lat_r_to_lon2_lat2_r'
! 5: 'lon2_lat2_r_to_lon_lat_r'

pure subroutine convert_0123(coor, matrix)
implicit none
real(8), intent(inout):: coor(0:2)
real(8), intent(out):: matrix(0:2, 0:2)
real(8) :: sin01(0:1), cos01(0:1), coor_in(0:2)
!------------------------------------------------------------
coor_in=coor
if (mode .eq. 0 .or. mode .eq. 2) then
    coor(2)= norm2(coor_in)
    if (mode .eq. 2) coor_in = cshift(coor_in, -1)
    coor(1)= asin(coor_in(2)/coor(2))
    coor(0)= xy2lon(coor_in(0:1))
else if (mode .eq. 1 .or. mode .eq. 3) then
    sin01= sin(coor_in(0:1))
    cos01= cos(coor_in(0:1))
    coor = coor_in(2)*[cos01(1)*[cos01(0), sin01(0)], sin01(1)]
    if (mode .eq. 3) coor= cshift(coor, -1)
endif

if (present1) then
    if (mode .eq. 0 .or. mode .eq. 2) then
        sin01= sin(coor(0:1))
        cos01= cos(coor(0:1))
    endif
    ! stack e_lon, e_lat, e_r for mode 0 'xyz_to_lon_lat_r' 
    matrix=reshape([-sin01(0),           cos01(0),      0.D0, &
                    -sin01(1)*[cos01(0), sin01(0)], cos01(1), &
                     cos01(1)*[cos01(0), sin01(0)], sin01(1)], (/3,3/))
    
    if (mode .eq. 2 .or. mode .eq. 3) matrix = cshift(matrix, -1, 1)
    if (mode .eq. 1 .or. mode .eq. 3) matrix = transpose(matrix)
endif
end subroutine convert_0123


pure subroutine convert_45(coor, matrix)
implicit none
real(8), intent(inout):: coor(0:2)
real(8), intent(out):: matrix(0:2, 0:2)
integer:: i, j 
real(8):: cos_1(0:1), sin_1(0:1), cos_2(0:1), sin_2(0:1), e1(0:2, 0:1), e2(0:2, 0:1)
!------------------------------------------------------------
cos_1=cos(coor(0:1))
sin_1=sin(coor(0:1))

if (mode .eq. 5) then
	coor(0)= xy2lon([sin_1(1), cos_1(0)*cos_1(1)])
	coor(1)= asin(cos_1(1)*sin_1(0))
else
	coor(0)= xy2lon([cos_1(1)*sin_1(0), sin_1(1)])
	coor(1)= asin(cos_1(1)*cos_1(0))
endif

if (present1) then
    cos_2=cos(coor(0:1))
    sin_2=sin(coor(0:1))
    
    e1(:,0)=[-sin_1(0), cos_1(0), 0.D0]
    e1(:,1)=[-sin_1(1)*[cos_1(0), sin_1(0)], cos_1(1)]
    e2(:,0)=[-sin_2(0), cos_2(0), 0.D0]
    e2(:,1)=[-sin_2(1)*[cos_2(0), sin_2(0)], cos_2(1)]

    if (mode .eq. 5) then
        e1 = cshift(e1, -1, 1)
    else
        e2 = cshift(e2, -1, 1)
    endif

    forall(i=0:1, j=0:1) matrix(i, j)=dot_product(e1(:, i), e2(:, j))
    matrix(0:1,2)=0.D0
    matrix(2,0:1)=0.D0
    matrix(2,2)=1.D0
endif
end subroutine convert_45


subroutine io_bin(name, array, readflag)
character(len=*):: name
real(8) :: array(:)
logical :: readflag
!------------------------------------------------------------
if (readflag) then
    open(1, file=name, access='stream', status='old')
    if (r4flag) then
        read(1) dummy
        array = dummy
    else
        read(1) array
    endif
    close(1)
else 
    open(1, file=name, access='stream', status='replace')
    if (r4flag) then
        dummy= array
        write(1) dummy
    else
        write(1) array
    endif
    close(1)
endif
end subroutine io_bin


end module share

program main
use share
implicit none
integer:: nthreads, OMP_GET_NUM_PROCS
integer(8):: k, ndata
logical::  present2, present3, present4
real(8):: matrix(0:2, 0:2)
real(8), allocatable :: coordinate(:), v1(:), v2(:), v3(:), v4(:)
character(len=1) :: str_aux
!------------------------------------------------------------
open(1, file='head.bin', access='stream', status='old')
read(1) mode, nthreads, r4flag, ndata
close(1)
!------------------------------------------------------------
if (r4flag) allocate(dummy(0:ndata-1))
allocate(coordinate(0:ndata-1))
call io_bin('coordinate.bin', coordinate, .true.)

inquire(file='v1.bin', exist=present1)
if (present1) then
    allocate(v1(0:ndata-1))
    call io_bin('v1.bin', v1, .true.)
endif

inquire(file='v2.bin', exist=present2)
if (present2) then
    allocate(v2(0:ndata-1))
    call io_bin('v2.bin', v2, .true.)
endif

inquire(file='v3.bin', exist=present3)
if (present3) then
    allocate(v3(0:ndata-1))
    call io_bin('v3.bin', v3, .true.)
endif

inquire(file='v4.bin', exist=present4)
if (present4) then
    allocate(v4(0:ndata-1))
    call io_bin('v4.bin', v4, .true.)
endif
!------------------------------------------------------------
! https://www.openmp.org/spec-html/5.0/openmpsu112.html
if (nthreads .gt. OMP_GET_NUM_PROCS()) nthreads=OMP_GET_NUM_PROCS()
if (nthreads .eq. 0) nthreads=OMP_GET_NUM_PROCS()-2
CALL OMP_set_num_threads(nthreads)

if (mode .le. 3) then
    convert => convert_0123 
else
    convert => convert_45
endif

two_pi=6.28318530717958647692D0
pi    =3.14159265358979323846D0

!!$OMP PARALLEL DO PRIVATE(k, matrix), schedule(static)
do k=0, ndata/3-1
    call convert(coordinate(k*3:k*3+2), matrix)
    if (present1) then
        v1(k*3:k*3+2)= MATMUL(v1(k*3:k*3+2), matrix)
        if (present2) v2(k*3:k*3+2)= MATMUL(v2(k*3:k*3+2), matrix)
        if (present3) v3(k*3:k*3+2)= MATMUL(v3(k*3:k*3+2), matrix)
        if (present3) v4(k*3:k*3+2)= MATMUL(v4(k*3:k*3+2), matrix)
    endif
enddo
!!$OMP END PARALLEL DO
!------------------------------------------------------------
call io_bin('coordinate_out.bin', coordinate, .false.)
deallocate(coordinate)

if (present1) then
    call io_bin('v1out.bin', v1, .false.)
    deallocate(v1)
endif

if (present2) then
    call io_bin('v2out.bin', v2, .false.)
    deallocate(v2)
endif

if (present3) then
    call io_bin('v3out.bin', v3, .false.)
    deallocate(v3)
endif

if (present4) then
    call io_bin('v4out.bin', v4, .false.)
    deallocate(v4)
endif

if (r4flag) deallocate(dummy)
!------------------------------------------------------------
! If the pop-up window for fastqsl.exe cannot be closed automatically on some Windows systems, please uncomment this line
! call system('taskkill /im convert_coordinate.exe /f')

! another way to kill the pop-up window
! call abort
end program main