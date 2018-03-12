program modelisation
use,intrinsic:: ieee_arithmetic, only: ieee_is_nan
Implicit none


! Parametres de normalisation
Real,parameter		:: N_o = 1.e12,T_o = 1000.
Real,parameter		:: Tmin=50.,Nliminf=1.e0
Real				:: C_o,P_o,Q_o,to,Ro

! Constantes terrestres
Real,parameter		:: go=9.80665,rayon=6378.

Integer,parameter		:: nbpt=500

Integer				:: npoint,Ipos1,Iposn,Iposnp

Integer,parameter 		:: nb_ion=6,ie=0
Integer,parameter		:: indi_O=1,indi_H=2,indi_N=3,indi_N2=4,indi_NO=5,indi_O2=6

Real,dimension (nbpt)		:: N1old,N1new,U1old,U1new,T1pold,T1told,T1old,T1pnew,T1tnew,T1new,q1old,q1new
Real,dimension(nbpt)		:: Po,Pr1,Lo1,D8n1,D7n1
Real,dimension(nbpt)		:: Vel1c,Vel1m,Vel1e,Vel1q
Real,dimension (nbpt)		:: N2old,N2new,U2old,U2new,T2pold,T2told,T2old,T2pnew,T2tnew,T2new,q2old,q2new
Real,dimension(nbpt)		:: Ph,Pr2,Lo2,D8n2,D7n2
Real,dimension(nbpt)		:: Vel2c,Vel2m,Vel2e,Vel2q
Real,dimension (nbpt)		:: N3old,N3new,U3old,U3new,T3pold,T3told,T3old,T3pnew,T3tnew,T3new,q3old,q3new
Real,dimension(nbpt)		:: Pn,Pr3,Lo3,D8n3,D7n3
Real,dimension(nbpt)		:: Vel3c,Vel3m,Vel3e,Vel3q
Real,dimension (nbpt)		:: N4old,N4new,Umold,Umnew,Tmpold,Tmtold,Tmold,Tmpnew,Tmtnew,Tmnew
Real,dimension (nbpt)		:: N5old,N5new,N6old,N6new
Real,dimension(nbpt)		:: Pn2,Pr4,Lo4,D8n4,D7n4,Pr5,Lo5,D8n5,D7n5,Po2,Pr6,Lo6,D8n6,D7n6
Real,dimension(nbpt)		:: Velmc,Velmm,Velme,Velmq

Integer,parameter		:: nb_neutre=6
Integer,parameter		:: indn_H=1,indn_N=2,indn_O=3,indn_N2=4,indn_NO=5,indn_O2=6

Real,dimension(nbpt)		:: Tn,Un,Vn,Wn,Nh,Nn,No,Nn2,Nno,No2

!Constantes de normalisation des equations
Real				:: Co(0:nb_ion),Cij(0:nb_ion,0:nb_ion),Qo(0:nb_ion)
!-------------------------------------------------------------------------------

!tableaux des masses
Real,dimension(0:nb_ion)	:: me_mi,mkg_ion
Real,dimension(nb_neutre)	:: mkg_neutre
!-------------------------------------------------------------------------------

!Coefficients pour le calcul des termes de collisions
Real,dimension(0:nb_ion,0:nb_ion)		:: muij,mirij,Ac,Dco,Dc1,Dc4
Real,dimension(0:nb_ion,0:nb_ion)		:: nuijo
Real,dimension(0:nb_ion,nb_neutre)		:: muin,mirin,mnrin,Am,Dmo,Dm1,Dm4
Real						:: lnCou


Real,dimension(nbpt)		:: alt,alt_geo,div_z,dV2
Real				:: div_tmp,U1r,T1r,Tmr,Teff,Teff2,Teff3,Teff4,reac

!-------------------------------------------------------------------------------

! reactions chimiques

Integer,parameter			:: nb_reac=27
Real,dimension(nb_reac)			:: rk
Real,	parameter,dimension(nb_reac)	:: rko	=(/2.5e-17,		&	!	O+  + H    --> H+  + O
						   1.0e-17,		&	! 	O+  + N2   --> NO+ + N 	energie dependante
                                                   1.0e-18,		&	! 	O+  + NO   --> NO+ + O 	energie dependante
                                                   1.0e-17,		&	! 	O+  + O2   --> O2+ + O 	energie dependante
                                                   2.2e-17,		&	! 	H+  + O    --> O+  + H

                                                   5.e-19,		&	! 	N+  + O    --> O+  + N
                                                   3.6e-18,		&	! 	N+  + H    --> H+  + N
                                                   2.e-17,		&	! 	N+  + NO   --> NO+ + N
                                                   3.7e-17,     	&	! 	N+  + O2   --> O+  + NO
                                                   2.6e-18,		&	! 	N+  + O2   --> NO+ + O
                                                   3.1e-16,		&	! 	N+  + O2   --> O2+ + N
                                                   		
                                                   0.,			&	! 	N2+ + O    --> O+  + N2	energie dependante
                                                   0.,			&	! 	N2+ + O    --> NO+ + N 	energie dependante
                                                   3.3e-16, 		&	! 	N2+ + NO   --> NO+ + N2
                                                   1.5e-14,  		&	! 	N2+ + O2   --> O2+ + N2
                                                   1.66e-12,  		&	! 	N2+ + e-   --> N   + N

                                                   5.36e-11,		&	!	NO+ + e-   --> O   + N

                                                   1.2e-22,		&	! 	O2+ + N    --> NO+ + O
                                                   5.e-22,  		&	! 	O2+ + N2   --> NO+ + NO
                                                   4.5e-16, 		&	! 	O2+ + NO   --> NO+ + O2
                                                   3.69e-12,		&	! 	O2+ + e-   --> O   + O
                                                   	
                                                   6.e-13, 		&	! 	NO  + hnu  --> NO+ + e-
                                                   8.3e-12,            	&	! 	NO  + hnu  --> N   + O
                                                   1.5e-18,		&	! 	NO  + N(4) --> N2  + O2
                                                   7.e-17,             	&	! 	NO  + N(2) --> N2  + O2
                                                                	
                                                   5.3e-18+3.5e-18,	&	! 	N(2)+ O2   --> NO  + O
                                                   1.5e-17            	/)	! 	N(4)+ O2   --> NO  + O
!---------------------------------------------------------------------------------------------------------------------------------

Call coefmass(nb_ion,nb_neutre,masse_ion,masse_neutre,muij,muin,mirij,mirin,mnrin,Ac,Am,Dco,Dc1,Dc4,Dmo,Dm1,Dm4)
Call collifreq(nb_ion,nb_neutre,no,T_o,muij,mirij,muin,mnrin,nuijo,nuino,lnCou)


! ***** Coefficients des reactions resonnantes (Schunk et Nagy 1980 p823) *****
nuino(indi_H,indn_H)  =2.65e-16
nuino(indi_O,indn_O)  =3.67e-17
nuino(indi_O2,indn_O2)=2.59e-17
nuino(indi_N2,indn_N2)=5.14e-17
nuino(indi_H,indn_O)  =6.61e-17
! ***** Coefficients des collisions electron-neutre (Schunk et Nagy 1980 p822) *****
nuino(indi_e,indn_H)  =4.5e-15
nuino(indi_e,indn_O)  =8.9e-17
nuino(indi_e,indn_O2) =1.82e-16
nuino(indi_e,indn_CO2)=3.68e-14
nuino(indi_e,indn_N2) =2.33e-17
nuino(indi_e,indn_O)  =8.9e-17



!************************************************************************
!*               DEFINITION DES PARAMETRES DE REFERENCE			*
!************************************************************************

Co=sqrt(kb*T_o/mkg_ion)        ! Definition des temps, vitesse, longueur, flux de chaleur de reference
qo=no*kb*T_o*Co
to=Co(iref)/go
Ro=Co(iref)*to



do i=1,npoint

  div_z(i)=-3.*alt_geo_1(i)
  T1r=T1new(i)*T_o
  U1r=U1new(i)*C1o
  Tmr=Tmnew(i)*T_o
  div_tmp=1.2e-4*U1r**2

!********************
! H+ + O -> O+ + H  *
!********************

  rk(1)=rko(1)*sqrt(Tn(i)+T1r/massion(1)+div_tmp)

!********************
! O+ + H -> H+ + O  *
!********************

  rk(5)=rko(5)*sqrt(Tn(i)/massion(1)+T1r+div_tmp)

!*********************
! O+ + N2 -> NO+ + N *		St-Maurice J.-P. et P.J. Laneville
!*********************

  Teff=mirin(indi_O,indn_N2)*Tn(i)	&
      +mnrin(indi_O,indn_N2)*T1r	&
      +muin(indi_O,indn_N2)*dV2
  Teff=Teff/1000.
  Teff2=Teff*Teff
  Teff3=Teff*Teff2
  Teff4=Teff3*Teff
  if (Teff <= 3725.) then
    reac  = 1.71676e-1 		&
          - 2.39978e-1*Teff     &
          + 1.48088e-1*Teff2	&
          - 3.43783e-2*Teff3  	&
          + 7.89577e-3*Teff4
  else
    reac  =-1.52489e0		&
          + 2.55704e-1*Teff	&
          + 1.32293e-1*Teff2	&
          - 4.84659e-3*Teff3  	&
          + 5.77477e-5*Teff4
  endif
  rk(2)=rko(2)*reac

!*********************
! O+ + NO -> NO+ + O *		St-Maurice J.-P. et P.J. Laneville
!*********************

  Teff=mirin(indi_O,indn_NO)*Tn(i)	&
      +mnrin(indi_O,indn_NO)*T1r	&
      +muin(indi_O,indn_NO)*dV2
  Teff=Teff/1000.
  Teff2=Teff*Teff
  Teff3=Teff*Teff2
  Teff4=Teff3*Teff
  if (Teff <= 3800.) then
    reac = 6.40408e0 		&
         - 4.46293e-1*Teff     	&
         + 8.50114e-1*Teff2	&
         - 1.15374e-1*Teff3  	&
         + 8.17746e-3*Teff4
  else
    reac =-7.48318e-1		&
         + 7.71673e-1*Teff	&
         + 3.41289e-1*Teff2	&
         - 9.83096e-3*Teff3  	&
         + 9.58846e-5*Teff4
  endif
  rk(3)=rko(3)*reac

!*********************
! O+ + O2 -> O2+ + O *		St-Maurice J.-P. et P.J. Laneville
!*********************

  Teff=mirin(indi_O,indn_O2)*Tn(i)	&
      +mnrin(indi_O,indn_O2)*T1r	&
      +muin(indi_O,indn_O2)*dV2
  Teff=Teff/1000.
  Teff2=Teff*Teff
  Teff3=Teff*Teff2
  if (Teff <= 4800.) then
    Teff4=Teff3*Teff
    reac = 2.78932e0		&
         - 2.30871e0*Teff    	&
         + 9.64093e-1*Teff2	&
         - 1.28611e-1*Teff3  	&
         + 6.26046e-3*Teff4
  else
    reac =-1.74046e0		&
         + 1.00776e0*Teff	&
         - 2.65793e-3*Teff2	&
         - 1.49035e-7*Teff3
  endif
  rk(4)=rko(4)*reac

!**********************
! N2+ + O -> O+ + N2  *
!**********************

  Teff=mirin(indi_N2,indn_O)*Tn(i)	&
      +mnrin(indi_N2,indn_O)*Tmr	&
      +muin(indi_N2,indn_O)*dV2
  if (Teff <= 1500.) then
    rk(12)=3.7131e-17*Teff**.23
  else
    rk(12)=3.7318e-17*Teff**.41
  endif

  rk(12)=rk(12)*No(i)

  rk(1)=rk(1)*No(i)


!********************
! N+ + O -> O+ + N  *
!********************

  rk(6)=rko(6)*No(i)

!**********************
! N2+ + O -> NO+ + N  *
!**********************

  rk(13)=rko(13)*No(i)


  rk(5)=rk(5)*Nh(i)

!********************
! N+ + H -> H+ + N  *
!********************

  rk(7)=rko(7)*Nh(i)


  rk(4)=rk(4)*No2(i)


!**********************
! N+ + O2 -> O+ + NO  *
!**********************

  rk(9)=rko(9)*No2(i)

!**********************
! N+ + O2 -> NO+ + O  *
!**********************

  rk(10)=rko(10)*No2(i)

!**********************
! N+ + O2 -> O2+ + N  *
!**********************

  rk(11)=rko(11)*No2(i)

!************************
! N2+ + O2 -> O2+ + N2  *
!************************

  rk(15)=rko(15)*No2(i)

!**********************
! O2+ + N -> NO+ + O  *
!**********************

  rk(18)=rko(18)*Nn(i)


  rk(2)=rk(2)*Nn2(i)


!************************
! O2+ + N2 -> NO+ + NO  *
!************************

  rk(19)=rko(19)*Nn2(i)


  rk(3)=rk(3)*Nno(i)

!**********************
! N+ + NO -> NO+ + N  *
!**********************

  rk(8)=rko(8)*Nno(i)

!************************
! N2+ + NO -> NO+ + N2  *
!************************

  rk(14)=rko(14)*Nno(i)

!************************
! O2+ + NO -> NO+ + O2  *
!************************

  rk(20)=rko(20)*Nno(i)

!*********************
! N2+ + e- -> N + N  *
!*********************

  rk(16)=rko(16)*Tenew(i)**(-.39)*Nenew(i)

!*********************
! NO+ + e- -> O + N  *
!*********************

  rk(17)=rko(17)*Tenew(i)**(-.85)*Nenew(i)

!*********************
! O2+ + e- -> O + O  *
!*********************

  rk(21)=rko(21)*Tenew(i)**(-.55)*Nenew(i)



  Pr1=Po(i)+rk(5)*N2new(i)+(rk(6)+rk(9))*N3new(i)+rk(12)*N4new(i)
  Lo1=rk(1)+rk(2)+rk(3)+rk(4)

  Pr2=Ph(i)+rk(1)*N2new(i)+rk(7)*N3new(i)
  Lo2=rk(5)

  Pr3=Pn(i)
  Lo3=rk(6)+rk(7)+r(8)+rk(9)+rk(10)+rk(11)

  Pr4=Pn2(i)
  Lo4=rk(12)+rk(13)+rk(14)+rk(15)+rk(16)

  Pr5=(rk(2)+rk(3))*N1new(i)+(rk(8)+rk(10))*N3new(i)		&
     +(rk(13)+rk(14))*N4new(i)+(rk(18)+rk(19)+rk(20))*N6new(i)
  Lo5=rk(17)

  Pr6=Po2(i)+rk(4)*N1new(i)+rk(11)*N3new(i)+rk(15)*N4new(i)
  Lo6=rk(18)+rk(19)+rk(20)+rk(21)

  D8n1(i)=Pr1
  D7n1(i)=-Lo1+div_z(i)*U1new(i)

  D8n2(i)=Pr2
  D7n2(i)=-Lo2+div_z(i)*U2new(i)

  D8n3(i)=Pr3
  D7n3(i)=-Lo3+div_z(i)*U3new(i)


  div_tmp=div_z(i)*Umnew(i)

  D8n4(i)=Pr4
  D7n4(i)=-Lo4+div_tmp

  D8n5(i)=Pr5
  D7n5(i)=-Lo5+div_tmp

  D8n6(i)=Pr6
  D7n6(i)=-Lo6+div_tmp

enddo

!------------------------------------!
! Equation de continuite de l'ion O+ !
!------------------------------------!

call velocity(Vel1c,Ipos1,Iposnp,deltat_2)

call sources(Ipos1,Iposn,deltat_2,8,zero,D8n1,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7n1,0.,0.)

lbc=1.
rbc=(N1old(npoint)/N1old(npoint-1))
call lcpfct(N1old,N1new,Ipos1,Iposn,lbc,0.,0.,N1new(np),.false.,1)


if (any(ieee_is_nan(N1new))) then
  print*,'probleme lors du calcul de N2new dans la boucle 1'
  do i=1,npoint
    print*,alt(i),N1old(i),Vel1c(i),D8n1(i),D7n1(i)
  enddo
  error stop
endif
do i=1,npoint
  N1new(i)=max(N1new(i),r_min)
enddo
N1new(np)=min(1.,N1new(npoint-1)/N1new(npoint-2))*N1new(npoint)

!------------------------------------!
! Equation de continuite de l'ion H+ !
!------------------------------------!

call velocity(Vel2c,Ipos1,Iposnp,deltat_2)

call sources(Ipos1,Iposn,deltat_2,8,zero,D8n2,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7n2,0.,0.)

lbc=1.
rbc=(N2old(npoint)/N2old(npoint-1))
call lcpfct(N2old,N2new,Ipos1,Iposn,lbc,0.,0.,N2new(np),.false.,1)


if (isnant(N2new,npoint)) then
  print*,'probleme lors du calcul de N1new dans la boucle 1'
  do i=1,npoint
    print*,alt(i),N2old(i),Vel2c(i),D8n2(i),D7n2(i)
  enddo
  error stop
endif
do i=1,npoint
  N2new(i)=max(N2new(i),r_min)
enddo
N2new(np)=min(1.,N2new(npoint-1)/N2new(npoint-2))*N2new(npoint)

!------------------------------------!
! Equation de continuite de l'ion N+ !
!------------------------------------!

call velocity(Vel3c,Ipos1,Iposnp,deltat_2)

call sources(Ipos1,Iposn,deltat_2,8,zero,D8n3,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7n3,0.,0.)

lbc=1.
rbc=(N3old(npoint)/N3old(npoint-1))
call lcpfct(N3old,N3new,Ipos1,Iposn,lbc,0.,0.,N3new(np),.false.,1)


if (isnant(N3new,npoint)) then
  print*,'probleme lors du calcul de N3new dans la boucle 1'
  do i=1,npoint
    print*,alt(i),N3old(i),Vel3c(i),D8n3(i),D7n3(i)
  enddo
  error stop
endif
do i=1,npoint
  N3new(i)=max(N3new(i),r_min)
enddo
N3new(np)=min(1.,N3new(npoint-1)/N3new(npoint-2))*N3new(npoint)

!-------------------------------------!
! Equation de continuite de l'ion N2+ !
!-------------------------------------!

call velocity(Velmc,Ipos1,Iposnp,deltat_2)

call sources(Ipos1,Iposn,deltat_2,8,zero,D8n4,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7n4,0.,0.)

lbc=1.
rbc=(N4old(npoint)/N4old(npoint-1))
call lcpfct(N4old,N4new,Ipos1,Iposn,lbc,0.,0.,N4new(np),.false.,1)


if (isnant(N4new,npoint)) then
  print*,'probleme lors du calcul de N4new dans la boucle 1'
  do i=1,npoint
    print*,alt(i),N4old(i),Velmc(i),D8n4(i),D7n4(i)
  enddo
  error stop
endif
do i=1,npoint
  N4new(i)=max(N4new(i),r_min)
enddo
N4new(np)=min(1.,N4new(npoint-1)/N4new(npoint-2))*N4new(npoint)

!-------------------------------------!
! Equation de continuite de l'ion NO+ !
!-------------------------------------!

call sources(Ipos1,Iposn,deltat_2,8,zero,D8n5,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7n5,0.,0.)

lbc=1.
rbc=(N5old(npoint)/N5old(npoint-1))
call lcpfct(N5old,N5new,Ipos1,Iposn,lbc,0.,0.,N5new(np),.false.,1)


if (isnant(N5new,npoint)) then
  print*,'probleme lors du calcul de N5new dans la boucle 1'
  do i=1,npoint
    print*,alt(i),N5old(i),Velmc(i),D8n5(i),D7n5(i)
  enddo
  error stop
endif
do i=1,npoint
  N5new(i)=max(N5new(i),r_min)
enddo
N5new(np)=min(1.,N5new(npoint-1)/N5new(npoint-2))*N5new(npoint)

!-------------------------------------!
! Equation de continuite de l'ion O2+ !
!-------------------------------------!

call sources(Ipos1,Iposn,deltat_2,8,zero,D8n6,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7n6,0.,0.)

lbc=1.
rbc=(N6old(npoint)/N6old(npoint-1))
call lcpfct(N6old,N4new,Ipos1,Iposn,lbc,0.,0.,N6new(np),.false.,1)


if (isnant(N6new,npoint)) then
  print*,'probleme lors du calcul de N6new dans la boucle 1'
  do i=1,npoint
    print*,alt(i),N6old(i),Velmc(i),D8n6(i),D7n6(i)
  enddo
  error stop
  
endif
do i=1,npoint
  N6new(i)=max(N6new(i),r_min)
enddo
N6new(np)=min(1.,N6new(npoint-1)/N6new(npoint-2))*N6new(npoint)


end

