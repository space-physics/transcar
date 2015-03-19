
Implicit none

integer :: i,j
real*8,dimension(2346) :: latitude,longitude,utem2
Integer, parameter :: npt=2346,ndeg=8,mdeg=8,rang=(2*mdeg+1)*(ndeg+1),&
			rang1=(ndeg+1)*(mdeg+1)
Integer :: nb_data,ind_data(npt),dpar11
real*8  :: latmin=40.d0,latmax=72.d0

Real*8  :: A(npt*rang),B(rang*npt)

Real*8  :: coef_r(rang),coef_i(rang)


!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!
open(unit=4,file='CoefP1.lis',form='formatted',status='unknown')
do i=1,2346	
	read(4,*) longitude(i),latitude(i),utem2(i)
	read(4,*)
end do
close(4)
dpar11=2346
print* , dpar11
!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!

!! CALCUL LES COEFFICIENTS

call geometrie(dpar11,nb_data,ind_data,longitude,latitude,latmin,&
		latmax,ndeg,mdeg,rang,A,B)

!!!!!!!!!!!!!!!!!!!!!! POTENTIEL

open(20,file='CoefRePot1.lis',form='formatted',status='unknown')

call coefficient(ndeg,mdeg,rang,dpar11,nb_data,ind_data,B,utem2,coef_r,coef_i)

do i=1,rang1
	write(20,*) coef_r(i),coef_i(i)
end do

close(20)

END
