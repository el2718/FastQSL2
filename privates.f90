function privates(Bp, CurlBp, Ap)
implicit none
real:: privates(0:9), Bp(0:2), CurlBp(0:2), Ap(0:2)
!------------------------------------------------------------
! for length
privates(0)= 1.

! CurlB \cdot B / (4*\pi*B^2); for twist
privates(1)= dot_product(CurlBp, Bp)/(4.*pi*dot_product(Bp, Bp))

 ! |CurlB|^2
privates(2)= dot_product(CurlBp, CurlBp)

! |CurlB|/|B|; norm2s() is faster than norm2()
privates(3)= norm2s(CurlBp)/norm2s(Bp)

! A \cdot b
! privates(4) = dot_product(Ap, Bp/norm2s(Bp))
end function privates


subroutine set_CurlBvec_Flag
implicit none
CurlBvec_Flag = CurlB_out .or. targetCurlB_out .or. loopCurlB_out &
.or. int_private_out(1) .or. int_private_out(2) .or. int_private_out(3)
end subroutine set_CurlBvec_Flag