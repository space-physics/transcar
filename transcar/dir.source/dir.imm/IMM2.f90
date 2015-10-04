SUBROUTINE IMM 

! Lire les sorties ELECTRO et les ecrire dans un fichier 
! Transformer cettes sorties au TRANSCAR

Implicit none

integer :: i,ndg,mdg
real*8,parameter :: pi=3.14159265358979d0
integer,parameter :: dpar11=(49*80)
real*8,dimension(dpar11) :: latitude
real*8,dimension(dpar11) :: longitude
character*20 :: nomf
real*8, dimension(dpar11) :: utem2, FEtem2, Etem2, alite2,&
			utem2_fit,FEtem2_fit,Etem2_fit,alite2_fit	

Integer, parameter :: npt=49*80,ndeg=8,mdeg=8,rang=(2*mdeg+1)*(ndeg+1),&
			rang1=(ndeg+1)*(mdeg+1),rangcoef=1000
Integer :: mode,nb_data,nb_data_pot,ind_data(npt),ind_data_pot(npt)

Real*8  :: A(npt*rang),B(rang*npt),Apot(npt*rang),Bpot(rang*npt)
real*8  :: latmin=50.d0,latmax=72.d0,latequi=40.d0, rmin=1E-6, Lmin,Lmax

!Real*8  :: coef_r(rang),coef_i(rang)
Real*8	:: coef(rang)

real*8 :: coefpot(rangcoef),coefflux(rangcoef),coefte(rangcoef),coefalg(rangcoef),&
	      phipot(rangcoef),phicour(rangcoef),phienerg(rangcoef),phifluxF(rangcoef)
					
common/phi/coefpot,coefflux,coefte,coefalg,Lmin,Lmax,ndg,mdg

Lmin=latmin
Lmax=latmax
ndg=ndeg
mdg=mdeg

open(unit=4,file='dir.source/dir.imm/Electro.dat',form='formatted',status='unknown')

read(4,*) (utem2(i),i=1,dpar11)   ! Potentiel sans corotation
read(4,*) (FEtem2(i),i=1,dpar11)  ! Flux d'energie
read(4,*) (Etem2(i),i=1,dpar11)   ! Temperature electronique
read(4,*) (Alite2(i),i=1,dpar11)  ! Courants alignes

close(4)

FEtem2(1:dpar11)=FEtem2(1:dpar11)*1000.
Etem2(1:dpar11) =Etem2(1:dpar11)*1000.
Alite2(1:dpar11)=Alite2(1:dpar11)*1.e10

where(FEtem2(1:dpar11)<1.)
  FEtem2=1.
endwhere
where(Etem2(1:dpar11)<1.)
  Etem2=1.
endwhere
FEtem2=log(FEtem2)
Etem2=log(Etem2)

!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!
! On lit les coordonnees    dans 2   tableaux longitude et   latitude
! Elles sont donnees dans 2 fichiers 'longitude.lis' et 'latitude.lis'

open(unit=4,file='dir.source/dir.imm/latitude.lis',form='formatted',status='unknown')
read(4,*) (latitude(i),i=1,dpar11)
close(4)

open(unit=4,file='dir.source/dir.imm/longitude.lis',form='formatted',status='unknown')
read(4,*) (longitude(i),i=1,dpar11)
do i=1,dpar11
	longitude(i)=longitude(i)*180.d0/pi	
end do
close(4)

open(10,file='dir.output/simu',form='formatted',status='unknown')

do i=1,dpar11
write(10,*)longitude(i),latitude(i),utem2(i),FEtem2(i),Etem2(i),Alite2(i)
enddo
close (10)



!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!

!! CALCUL LES COEFFICIENTS

mode=0
call geometrie(dpar11,nb_data,ind_data,longitude,latitude,latmin,&
		latmax,latequi,mode,ndeg,mdeg,rang,A,B)
mode=1
call geometrie(dpar11,nb_data_pot,ind_data_pot,longitude,latitude,latmin,&
		latmax,latequi,mode,ndeg,mdeg,rang,Apot,Bpot)
!!!!!!!!!!!!!!!!!!!!!! POTENTIEL

open(10,file='dir.output/coefpot.lis',form='formatted',status='unknown')
open(20,file='dir.output/coefflux.lis',form='formatted',status='unknown')
open(30,file='dir.output/coefenergie.lis',form='formatted',status='unknown')
open(40,file='dir.output/coefcour.lis',form='formatted',status='unknown')

call coefficient(ndeg,mdeg,rang,dpar11,nb_data,ind_data_pot,Bpot,utem2,coefpot)
open(4,file='/home/blelly/dir.transcar/dir.output/coefpot3.lis',form='formatted',status='unknown')

do i=1,(ndg+1)*(2*mdg+1)
  read(4,*)coefpot(i)
enddo
close(4)
do i=1,rang
  write(10,*) coefpot(i)
end do
!do i=1,rang1
!	coefpot(i)=dcmplx(coef_r(i),coef_i(i))
!	write(10,*) coef_r(i),coef_i(i)
!end do

close(10)

!!!!!!!!!!!!!!!!! FLUX D'ENERGIE ELECTRONIQUE


call coefficient(ndeg,mdeg,rang,dpar11,nb_data,ind_data,B,FEtem2,coefflux)
do i=1,rang
  write(20,*) coefflux(i)
end do

!do i=1,rang1
!	coefflux(i)=dcmplx(coef_r(i),coef_i(i))
!	write(20,*) coef_r(i),coef_i(i)
!end do
close(20)

!!!!!!!!!!!!!!!!!!! TEMPERATURE ELECTRONIQUE

call coefficient(ndeg,mdeg,rang,dpar11,nb_data,ind_data,B,Etem2,coefte)
do i=1,rang
  write(30,*) coefte(i)
end do

!do i=1,rang1
!	coefte(i)=dcmplx(coef_r(i),coef_i(i))
!	write(30,*) coef_r(i),coef_i(i)
!end do
close(30)
!!!!!!!!!!!!!!!!!!!!!!! COURANTS ALIGNES

call coefficient(ndeg,mdeg,rang,dpar11,nb_data,ind_data,B,Alite2,coefalg)
do i=1,rang
  write(40,*) coefalg(i)
end do

!do i=1,rang1
!	coefalg(i)=dcmplx(coef_r(i),coef_i(i))
!	write(40,*) coef_r(i),coef_i(i)
!end do

close(40)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VERIFIER EN CALCULANT LES FONCTIONS FIT
! Electrofit.dat pour verifier les coefficients obtenus

open(9,file='dir.source/dir.imm/Electrofit.dat',form='formatted',status='unknown')

call projection(nb_data,ind_data,rang,A,B,utem2,utem2_fit)
write(9,*)(utem2_fit(i),i=1,dpar11)


call projection(nb_data,ind_data,rang,A,B,FEtem2,FEtem2_fit)

do i=1,dpar11
	FEtem2_fit(i)=exp(FEtem2_fit(i))/1000.
end do
write(9,*)(FEtem2_fit(i),i=1,dpar11)


call projection(nb_data,ind_data,rang,A,B,Etem2,Etem2_fit)

do i=1,dpar11
	Etem2_fit(i)=exp(Etem2_fit(i))/1000.
end do
write(9,*)(Etem2_fit(i),i=1,dpar11)


call projection(nb_data,ind_data,rang,A,B,Alite2,Alite2_fit)

do i=1,dpar11
	Alite2_fit(i)=Alite2_fit(i)/1E7
end do

write(9,*)(Alite2_fit(i),i=1,dpar11)

close(9)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*, 'Deja fait'
END SUBROUTINE IMM
