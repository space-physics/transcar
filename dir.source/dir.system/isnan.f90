pure logical function isnant(x,nx)
implicit none

integer,intent(in) :: nx
real,intent(in) :: x(nx)

integer i

isnant=.false.

do i=1,nx
  if (isnan(x(i))) then
    isnant=.true.
    return
  endif
enddo

end function isnant
