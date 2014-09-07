logical function isnan(x)

real :: x
	
!isnan=(x.ne.x)
isnan=.false.
if (exponent(x)==129.and.fraction(x)/=0.) then
  isnan=.true.
endif
return
end function isnan
	
logical function isnant(x,nx)
	
integer :: nx
real :: x(nx)
real :: y
logical :: isnan
external isnan
	
isnant=.false.
	
do i=1,nx
  if (exponent(x(i))==129.and.fraction(x(i))/=0.) then
    isnant=.true.
    return
  endif
enddo
return
end function isnant
	
	

