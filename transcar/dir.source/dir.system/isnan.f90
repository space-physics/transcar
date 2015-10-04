! don't use this unless you have a junky compiler without isnan()

!logical function isnan(x)
!
!real :: x
!	
!!isnan=(x.ne.x)
!isnan=.false.
!if (exponent(x)==129.and.fraction(x)/=0.) then
!  isnan=.true.
!endif
!
!end function isnan


pure logical function isnant(x,nx)
implicit none
!
integer,intent(in) :: nx
real,intent(in) :: x(nx)
!
integer i

isnant=.false.
	
do i=1,nx
  !if (exponent(x(i))==129.and.fraction(x(i))/=0.) then
  if (isnan(x(i))) then
    isnant=.true.
    return
  endif
enddo

end function isnant
	
	

