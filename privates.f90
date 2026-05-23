subroutine privates(site)
implicit none
type(site_info), target:: site
real:: vp(0:2)            ! vector position, i.e. the coordinates
real, pointer:: Bp(:)     ! magnetic field at vp
real, pointer:: CurlBp(:) ! CurlB field at vp
real, pointer:: Ap(:)     ! the additional field at vp
!------------------------------------------------------------
if (site%yinFlag) then
! if spherical, and |latitude=site%v(1)| .ge. lat_pole
    vp=site%v_yin(0:2)
    Bp => site%B_yin
    CurlBp => site%CurlB_yin
    Ap => site%A_yin
else
    vp=site%v(0:2)
    Bp => site%B
    CurlBp => site%CurlB
    Ap => site%A
endif
!------------------------------------------------------------
! for length
if (int_private_out(0)) &
site%private(0)= 1.

! CurlB \cdot B / (4*\pi*B^2); for twist
if (int_private_out(1)) &
site%private(1)= dot_product(CurlBp, Bp)/(4.*pi*dot_product(Bp, Bp))

! |CurlB|^2
if (int_private_out(2)) &
site%private(2)= dot_product(CurlBp, CurlBp)

! |CurlB|/|B|; norm2s() is faster than norm2()
if (int_private_out(3)) &
site%private(3)= norm2s(CurlBp)/norm2s(Bp)

! A \cdot b
! if (int_private_out(4)) &
! site%private(4) = dot_product(Ap, Bp/norm2s(Bp))
end subroutine privates


subroutine set_CurlBvec_Flag
implicit none
CurlBvec_Flag = CurlB_out .or. targetCurlB_out .or. loopCurlB_out &
.or. int_private_out(1) .or. int_private_out(2) .or. int_private_out(3)
end subroutine set_CurlBvec_Flag