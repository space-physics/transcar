!SUBROUTINE IMM

!!! UTILISER IMM POUR SORTIR POTENTIEL, FLUX D'ENERGIE, TEMPERATURE ENERGETIQUE
!!!             COURANTS ALIGNES ET LES ECRIRE DANS "ELECTRO.DAT"

!!!              TRANSFORMER CETTES SORTIES AU TRANSCAR

IMPLICIT NONE

integer :: i,ndg,mdg,ie
real*8,parameter :: pi=3.14159265358979d0
integer,parameter :: dpar11=(49*80)

real*8,dimension(dpar11) :: latitude
real*8,dimension(dpar11) :: longitude


real*4, dimension(dpar11) :: utem2, FEtem2, Etem2, alite2,&
				FEtem3,Etem3,FEion,Fion
real*8,dimension(dpar11) :: u2,FE2,E2,A2,utem2_fit,FEtem2_fit,Etem2_fit,&
				Alite2_fit
				
common/transm/utem2,FEtem3,Etem3,alite2
common/pot/FEtem2,Etem2

Integer, parameter :: npt=49*80,ndeg=8,mdeg=8,rang=(2*mdeg+1)*(ndeg+1),&
			rang1=(ndeg+1)*(mdeg+1),rangcoef=1000
Integer :: nb_data,ind_data(npt)

real*8  :: latmin=40.d0,latmax=72.d0,rmin=1E-6,Lmin,Lmax

Real*8  :: coef_r(rang),coef_i(rang)
Real*8  :: A(rang*npt),B(rang*npt)

complex*16 :: coefpot(rangcoef),coefflux(rangcoef),coefte(rangcoef),coefalg(rangcoef),&
	      phipot(rangcoef),phicour(rangcoef),phienerg(rangcoef),phifluxF(rangcoef)
	      
				
common/phi/coefpot,coefflux,coefte,coefalg,Lmin,Lmax,ndg,mdg

Logical :: continu

Lmin=latmin
Lmax=latmax
ndg=ndeg
mdg=mdeg

!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!
! ON LIT LES COORDONNEES   DANS 2   TABLEAUX LONGITUDE ET LATITUDE
! ELLES SONT DONNEES DANS 2 FICHIERS 'longitude.lis' ET 'latitude.lis'

!open(unit=4,file='dir.source/dir.imm/dir.in_out/latitude.lis',form='formatted',status='unknown')
open(unit=4,file='dir.in_out/latitude.lis',form='formatted',status='unknown')

read(4,*) (latitude(i),i=1,dpar11)
close(4)

!open(unit=4,file='dir.source/dir.imm/dir.in_out/longitude.lis',form='formatted',status='unknown')
open(unit=4,file='dir.in_out/longitude.lis',form='formatted',status='unknown')

read(4,*) (longitude(i),i=1,dpar11)

do i=1,dpar11
	longitude(i)=longitude(i)*180.d0/pi	
end do

close(4)

!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!!

! ON VA ECRIRE LES SORTIES DE IMM DANS LE FICHIER 'Electro.dat'


! Compteur d'iterations
ie=0

continu=.true.

do while(continu)
 
	CALL DRIVE11(ie,continu)

	!open(unit=4,file='dir.source/dir.imm/dir.in_out/Electro.dat',form='formatted',status='unknown')

	!write(4,*) (utem2(i),i=1,dpar11)   ! POTENTIEL SANS CORATATION
	!write(4,*) (FEtem2(i),i=1,dpar11)  ! FLUX D'ENERGIE
	!write(4,*) (Etem2(i),i=1,dpar11)   ! TEMPERATURE ELECTRONIQUE
	!write(4,*) (Alite2(i),i=1,dpar11)  ! Courants alignes

	!close(4)

	u2=dble(utem2)
	FE2=dble(FEtem2)
	E2=dble(Etem2)
	A2=dble(alite2)

	do i=1,dpar11
	
		FE2(i)=FE2(i)*1000.
		E2(i)=E2(i)*1000.
	
		if (FE2(i)< rmin) FE2(i)=1.
		if (E2(i) < rmin) E2(i)=1.
	
		FE2(i)=log(FE2(i))
		E2(i) =log(E2(i))

		A2(i)=A2(i)*1.E7		

	end do


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! CALCUL LES COEFFICIENTS

	call geometrie(dpar11,nb_data,ind_data,longitude,latitude,latmin,&
			latmax,ndeg,mdeg,rang,A,B)



	!!!!!!!!!!!!!!!!!!!!!! POTENTIEL

	call coefficient(ndeg,mdeg,rang,dpar11,nb_data,ind_data,B,u2,coef_r,coef_i)

	do i=1,rang1
		coefpot(i)=(-coef_r(i),coef_i(i))
	end do


	!!!!!!!!!!!!!!!!! FLUX D'ENERGIE ELECTRONIQUE

	call coefficient(ndeg,mdeg,rang,dpar11,nb_data,ind_data,B,FE2,coef_r,coef_i)

	do i=1,rang1
		coefflux(i)=(-coef_r(i),coef_i(i))
	end do



	!!!!!!!!!!!!!!!!!!! TEMPERATURE ELECTRONIQUE

	call coefficient(ndeg,mdeg,rang,dpar11,nb_data,ind_data,B,E2,coef_r,coef_i)

	do i=1,rang1
		coefte(i)=(-coef_r(i),coef_i(i))
	end do


	!!!!!!!!!!!!!!!!!!!!!!! COURANTS ALIGNES

	call coefficient(ndeg,mdeg,rang,dpar11,nb_data,ind_data,B,A2,coef_r,coef_i)

	do i=1,rang1 
		coefalg(i)=(-coef_r(i),coef_i(i))
	end do

end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! VERIFIER EN CALCULANT LES FONCTIONS FIT
	! Electrofit.dat pour verifier les coefficients obtenus

!	open(9,file='dir.source/dir.imm/dir.in_out/Electrofit.dat',form='formatted',status='unknown')
	open(9,file='dir.in_out/Electrofit.dat',form='formatted',status='unknown')

	call projection(nb_data,ind_data,rang,A,B,u2,utem2_fit)
	write(9,*)(utem2_fit(i),i=1,dpar11)


	call projection(nb_data,ind_data,rang,A,B,FE2,FEtem2_fit)
	do i=1,dpar11
		FEtem2_fit(i)=exp(FEtem2_fit(i))/1000.
	end do
	write(9,*)(FEtem2_fit(i),i=1,dpar11)


	call projection(nb_data,ind_data,rang,A,B,E2,Etem2_fit)
	do i=1,dpar11
		Etem2_fit(i)=exp(Etem2_fit(i))/1000.
	end do
	write(9,*)(Etem2_fit(i),i=1,dpar11)


	call projection(nb_data,ind_data,rang,A,B,A2,Alite2_fit)

	do i=1,dpar11
		Alite2_fit(i)=Alite2_fit(i)/1E7
	end do

	write(9,*)(Alite2_fit(i),i=1,dpar11)

	close(9)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


print*, 'Deja fait'

!END SUBROUTINE IMM
END
