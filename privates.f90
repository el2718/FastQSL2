function privates(Bp, CurlBp, Ap)
implicit none
real:: privates(0:9), Bp(0:2), CurlBp(0:2), Ap(0:2)
!------------------------------------------------------------
! for length
privates(0)= 1.

! CurlB cdot B / (4*\pi*B^2); for twist
privates(1)= sum(CurlBp*Bp)/(4.*pi*sum(Bp*Bp))

 ! |CurlB|^2
privates(2)= sum(CurlBp*CurlBp)

! |CurlB|/|B|; norm2s() is faster than norm2()
! privates(3)= norm2s(CurlBp)/norm2s(Bp) 
privates(3)= norm2s(CurlBp)/norm2s(Ap)

! A \ cdot b
! privates(4) = sum (Ap*Bp/ norm2s (Bp))
end function privates