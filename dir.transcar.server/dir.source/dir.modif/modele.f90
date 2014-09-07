subroutine modele(buffer,point_temp,ipos_i,ipos_o)

implicit none

! parametres globaux

type t_position
  Integer			:: z,n1,n2,n3,n4,n5,n6
  Integer			:: u1,u2,u3,u4,u5,u6,ue
  Integer			:: t1p,t1t,t2p,t2t,t3p,t3t,t4p,t4t,t5p,t5t,t6p,t6t,tep,tet
  Integer			:: q1,q2,q3,q4,q5,q6,qe
  Integer			:: po,ph,pn,pn2,po2,heat
  Integer			:: nno,uno,no,nh,nn,nn2,no2,tn,un,vn,wn
  Integer			:: nes,jes,tes,qes
end type t_position





Integer, parameter	:: npt=256,nbuf=10000
Integer,parameter	:: nb_ion=6,nb_spi=nb_ion+1,	indi_O =1,	&	! indice de positionnement pour O+
							indi_H =2,	&	! indice de positionnement pour H+
							indi_N =3,	&	! indice de positionnement pour N+
							indi_N2=4,	&	! indice de positionnement pour N2+
							indi_NO=5,	&	! indice de positionnement pour NO+
							indi_O2=6,	&	! indice de positionnement pour O2+
							indi_e =nb_spi   	! indice de positionnement pour les electrons

Integer,parameter	:: nb_spn=6,			indn_O =1,	&	! indice de positionnement pour O
							indn_H =2,	&	! indice de positionnement pour H
							indn_N =3,	&	! indice de positionnement pour N
							indn_N2=4,	&	! indice de positionnement pour N2
							indn_NO=5,	&	! indice de positionnement pour NO
							indn_O2=6		! indice de positionnement pour O2

							
!variables d'entrees sorties

Real			:: buffer(nbuf)
type(t_position)	:: ipos_i,ipos_o
integer			:: iapproxi,iapproxo,point_temp



!  parametres geometriques

Real	:: z(npt),alt(npt),G(npt),Gr(npt),alt_geo(npt),alt_geo_1(npt)

  ! parametres normalises pour le modele

Real	:: Radn(npt),Radl,Radr,extra(5)

! Variables internes pour la resolution

Real	:: C2a(npt),D2a(npt),C2b(npt),D2b(npt),C2c(npt),D2c(npt),C2d(npt),D2d(npt)
Real	:: y0(npt),y1(npt),chim(nb_ion,nb_ion),diff(nb_ion+nb_spi,nb_ion+nb_spi)

Real	:: dens_i(nb_spi),dens_i_1(nb_spi),temperature(nb_spi),vitesse(nb_spi)
Real	:: Tij(nb_spi,nb_spi),Tij_r(nb_spi,nb_spi),Tij_r_1(nb_spi,nb_spi),Tij_r_3(nb_spi,nb_spi)
Real	:: Tr(nb_spi,nb_spn),Tin(nb_spi,nb_spn),friction(nb_spi)
Real	:: nuij(nb_spi,nb_spi),nuin(nb_spi,nb_spn),coefnu(nb_spi,nb_spn)
Real	:: dens_n(nb_spn)


logical	:: flag_atmos,flag_cinetic,flag_conv
logical	:: isnant

! especes ioniques

Real	:: N1new(npt),U1new(npt),T1pnew(npt),T1tnew(npt),T1new(npt),Q1new(npt)
Real	:: N1old(npt),U1old(npt),T1pold(npt),T1told(npt),T1old(npt),Q1old(npt)
Real	:: Vel1(npt),Vel1_2(npt)

Real	:: N2new(npt),U2new(npt),T2pnew(npt),T2tnew(npt),T2new(npt),Q2new(npt)
Real	:: N2old(npt),U2old(npt),T2pold(npt),T2told(npt),T2old(npt),Q2old(npt)
Real	:: Vel2(npt),Vel2_2(npt)

Real	:: N3new(npt),U3new(npt),T3pnew(npt),T3tnew(npt),T3new(npt),Q3new(npt)
Real	:: N3old(npt),U3old(npt),T3pold(npt),T3told(npt),T3old(npt),Q3old(npt)
Real	:: Vel3(npt),Vel3_2(npt)

Real	:: N4new(npt),U4new(npt),T4pnew(npt),T4tnew(npt),T4new(npt),Q4new(npt)
Real	:: N4old(npt),U4old(npt),T4pold(npt),T4told(npt),T4old(npt),Q4old(npt)
Real	:: Vel4(npt),Vel4_2(npt)

Real	:: N5new(npt),U5new(npt),T5pnew(npt),T5tnew(npt),T5new(npt),Q5new(npt)
Real	:: N5old(npt),U5old(npt),T5pold(npt),T5told(npt),T5old(npt),Q5old(npt)
Real	:: Vel5(npt),Vel5_2(npt)

Real	:: N6new(npt),U6new(npt),T6pnew(npt),T6tnew(npt),T6new(npt),Q6new(npt)
Real	:: N6old(npt),U6old(npt),T6pold(npt),T6told(npt),T6old(npt),Q6old(npt)
Real	:: Vel6(npt),Vel6_2(npt)

Real	:: Nenew(npt),Uenew(npt),Tepnew(npt),Tetnew(npt),Tenew(npt),Qenew(npt)
Real	:: Neold(npt),Ueold(npt),Tepold(npt),Tetold(npt),Teold(npt),Qeold(npt)
Real	:: Vele(npt),Vele_2(npt)

! Atmosphere neutre

Real	:: Nnonew(npt),Nnoold(npt),Unonew(npt),Unoold(npt),Velno(npt),Velno_2(npt)
Real	:: Nh(npt),Nn(npt),No(npt),Nn2(npt),Nno(npt),No2(npt)
Real	:: Ph(npt),Pn(npt),Po(npt),Pn2(npt),         Po2(npt),Heat(npt)
Real	:: Un(npt),Vn(npt),Wn(npt),Tn(npt),vent(2,npt)

! Conditions geophysiques

Real	:: f107(3),Ap(7)

! reactions chimiques

Integer,parameter			:: nb_reac=27
Real					:: rk(nb_reac)
Real,	parameter,dimension(nb_reac)	:: rko	=(/2.5e-17,		&	!	O+  + H    --> H+  + O
						   0.,			&	! 	O+  + N2   --> NO+ + N 	energie dependante
                                                   0.,			&	! 	O+  + NO   --> NO+ + O 	energie dependante
                                                   0.,			&	! 	O+  + O2   --> O2+ + O 	energie dependante
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


interface


  ! routine de resolution numerique

  subroutine convec(iyd,temps,kp,dlongeo,dlatgeo,dlonmag,dlatmag,dlonref,dt,pot,flag_pot)
  integer	:: iyd,kp
  real*8	:: temps,dlongeo,dlatgeo,dlonmag,dlatmag,dlonref,dt,pot
  logical	:: flag_pot

  end subroutine convec
  subroutine  lcpfct(Rhoo,Rhon,I1,In,lbc,vlbc,rbc,vrbc,flag,mode)
    integer	:: I1,In,mode
    logical	:: flag
    real	:: Rhoo(:),Rhon(:),lbc,vlbc,rbc,vrbc
  end subroutine lcpfct

  subroutine makegrid(Rado,Radn,I1,In,alpha)
    integer	:: I1,In
    real	:: Rado(:),Radn(:),alpha
  end subroutine makegrid

  subroutine velocity(Vel,I1,In,dt)
    integer	:: I1,In
    real	:: Vel(:),dt
  end subroutine velocity

  subroutine  sources(I1,In,dt,mode,C2,D2,D2l,D2r)
    integer	:: I1,In,mode
    real	:: C2(:),D2(:),D2l,D2r,dt
  end subroutine sources

  subroutine solve_ODE(n,y0,y1,a,mat,dt)
    integer	:: n
    real	:: dt
    real*8	:: y0(:),y1(:),a(n,n),mat(:,:)
  end subroutine solve_ODE

  ! routines systeme

  subroutine pas_de_temps(iyd,sec,dt,dto)
    integer	:: iyd
    real	:: sec,dt,dto
  end subroutine pas_de_temps

  subroutine jour_mois(iyd,jour,mois)
    integer	:: iyd,jour,mois
  end subroutine jour_mois

  ! entrees-sorties

  subroutine write_out(buffer,longbuf,longrec,nrec)
    integer	:: longbuf,longrec,nrec
    real	:: buffer(:)
  end subroutine write_out

  subroutine lec_indices(an,mois,jour,heure,ap,f107)
    integer	:: an,mois,jour,heure
    real	:: ap(7),f107(3)
  end subroutine lec_indices

  ! conditions geophysiques et atmospheriques

  subroutine msis90(iyd,temps,z,latgeo,longeo,f107,ap,mass,d,t,w)
    integer	:: iyd,mass
    real	:: temps,z,latgeo,longeo,f107(3),ap(7),d(8),t(2),w(2)
  end subroutine msis90

  !  couplage cinetique

  subroutine transelec(npt,iyd,temps,latgeo,longeo,f107,ap,chi,			&
    		   Ne,Te,T1,nb_alt,Nh,No,No2,Nn2,Nn,Tn,indlim,jpreci,		&
                   N_o,T_o,kiappel,z,Heat,Ph,Po,Po2,Pn2,Pn,Nes,Jes,Tes,qes)    	

    integer	:: npt,nb_alt,iyd,kiappel,indlim,jpreci
    real	:: temps,latgeo,longeo,f107(3),ap(7),N_o,T_o
    real	:: Ne(:),Te(:),Nh(:),No(:),No2(:),Nn2(:),Nn(:),Tn(:),z(:)
    real	:: Heat(:),Ph(:),Po(:),Po2(:),Pn2(:),Pn(:),Nes(:),Jes(:),Tes(:),qes(:)
  end subroutine transelec


end interface





do iz=1,nb_alt
  ipos=(iz+1)*ncol
  alt(iz)      = buffer(ipos+ipos_i%z)

  N1new(iz)    = buffer(ipos+ipos_i%n1)/N_o
  N2new(iz)    = buffer(ipos+ipos_i%n2)/N_o
  N3new(iz)    = buffer(ipos+ipos_i%n3)/N_o
  N4new(iz)    = buffer(ipos+ipos_i%n4)/N_o
  N5new(iz)    = buffer(ipos+ipos_i%n5)/N_o
  N6new(iz)    = buffer(ipos+ipos_i%n6)/N_o
  U1new(iz)    = buffer(ipos+ipos_i%u1)/Co(1)
  U2new(iz)    = buffer(ipos+ipos_i%u2)/Co(2)
  U3new(iz)    = buffer(ipos+ipos_i%u3)/Co(3)
  U4new(iz)    = buffer(ipos+ipos_i%u4)/Co(4)
  U5new(iz)    = buffer(ipos+ipos_i%u5)/Co(5)
  U6new(iz)    = buffer(ipos+ipos_i%u6)/Co(6)
  Uenew(iz)    = buffer(ipos+ipos_i%ue)/Co(indi_e)
  T1pnew(iz)   = buffer(ipos+ipos_i%t1p)/T_o
  T1tnew(iz)   = buffer(ipos+ipos_i%t1t)/T_o
  T2pnew(iz)   = buffer(ipos+ipos_i%t2p)/T_o
  T2tnew(iz)   = buffer(ipos+ipos_i%t2t)/T_o
  T3pnew(iz)   = buffer(ipos+ipos_i%t3p)/T_o
  T3tnew(iz)   = buffer(ipos+ipos_i%t3t)/T_o
  T4pnew(iz)   = buffer(ipos+ipos_i%t4p)/T_o
  T4tnew(iz)   = buffer(ipos+ipos_i%t4t)/T_o
  T5pnew(iz)   = buffer(ipos+ipos_i%t5p)/T_o
  T5tnew(iz)   = buffer(ipos+ipos_i%t5t)/T_o
  T6pnew(iz)   = buffer(ipos+ipos_i%t6p)/T_o
  T6tnew(iz)   = buffer(ipos+ipos_i%t6t)/T_o
  Tepnew(iz)   = buffer(ipos+ipos_i%tep)/T_o
  Tetnew(iz)   = buffer(ipos+ipos_i%tet)/T_o
  q1new(iz)    = buffer(ipos+ipos_i%q1)/q_o(1)
  q2new(iz)    = buffer(ipos+ipos_i%q2)/q_o(2)
  q3new(iz)    = buffer(ipos+ipos_i%q3)/q_o(3)
  q4new(iz)    = buffer(ipos+ipos_i%q4)/q_o(4)
  q5new(iz)    = buffer(ipos+ipos_i%q5)/q_o(5)
  q6new(iz)    = buffer(ipos+ipos_i%q6)/q_o(6)
  qenew(iz)    = buffer(ipos+ipos_i%qe)/q_o(indi_e)
  Nnonew(iz)   = buffer(ipos+ipos_i%nno)/N_o
  Unonew(iz)   = buffer(ipos+ipos_i%uno)/Co(indi_NO)
  Po(iz)       = buffer(ipos+ipos_i%po)*to/N_o
  Ph(iz)       = buffer(ipos+ipos_i%ph)*to/N_o
  Pn(iz)       = buffer(ipos+ipos_i%pn)*to/N_o
  Pn2(iz)      = buffer(ipos+ipos_i%pn2)*to/N_o
  Po2(iz)      = buffer(ipos+ipos_i%po2)*to/N_o
  Heat(iz)     = buffer(ipos+ipos_i%heat)

  Nenew(iz)    = N1new(iz)+N2new(iz)+N3new(iz)+N4new(iz)+N5new(iz)+N6new(iz)

  T1new(iz)    = r1_3*T1pnew(iz)+r2_3*T1tnew(iz)
  T2new(iz)    = r1_3*T2pnew(iz)+r2_3*T2tnew(iz)
  T3new(iz)    = r1_3*T3pnew(iz)+r2_3*T3tnew(iz)
  T4new(iz)    = r1_3*T4pnew(iz)+r2_3*T4tnew(iz)
  T5new(iz)    = r1_3*T5pnew(iz)+r2_3*T4tnew(iz)
  T6new(iz)    = r1_3*T6pnew(iz)+r2_3*T4tnew(iz)
  Tenew(iz)    = r1_3*Tepnew(iz)+r2_3*Tetnew(iz)

  alt_geo  (iz) = (Re+1000.*alt(iz))/Ro
  alt_geo_1(iz) = 1./alt_geo(iz)
  G        (iz) = (Re/Ro/alt_geo(iz))**2

enddo



	
Radn (1) = (3.*alt_geo(1)-alt_geo(2))/2.
Radn(np) = (3.*alt_geo(nb_alt)-alt_geo(nb_alt-1))/2.
gr   (1) = go*.5*(3.*G(1)-G(2))*amu/kb*1000.

do i=2,nb_alt
  Radn(iz) = .5*(alt_geo(iz)+alt_geo(iz-1))
  gr  (iz) = go*.5*(G(iz)+G(iz-1))*amu/kb*1000.
enddo


call makegrid(Radn,Radn,Ipos1,Iposnp,alpha)


extra(1) = alt_geo(nb_alt-2)
extra(2) = alt_geo(nb_alt-1)
extra(3) = alt_geo(nb_alt)
extra(4) = (alt_geo(nb_alt-2)+alt_geo(nb_alt-1)+alt_geo(nb_alt))/3.
extra(5) = (alt_geo(nb_alt-2)**2+alt_geo(nb_alt-1)**2+alt_geo(nb_alt)**2)/3.
extra(5) = extra(5)-extra(4)**2



!debut de la boucle temporelle
temps_ecriture=0.
temps_atmos=0.
temps_cinetic=0.

nrec=0

dt_max=0.
do iz=1,nb_alt
  vmax=max(abs(U1new(iz)*Co(1)),abs(U2new(iz)*Co(2)))
  vmax=max(vmax,abs(U3new(iz)*Co(3)))
  vmax=max(vmax,abs(U4new(iz)*Co(4)))
  vmax=max(vmax,abs(U5new(iz)*Co(5)))
  vmax=max(vmax,abs(U6new(iz)*Co(6)))
enddo


do while (temps_TU<=temps_fin.or.iyd_TU/=iyd_fin)	! debut de la boucle temporelle sur la ligne de champ

  !------------------------------------------------------------!	
  ! 	memorisation des resultats de la boucle precedente     !
  !------------------------------------------------------------!

  do iz=1,nb_alt

    N1old(iz) = N1new(iz)
    N2old(iz) = N2new(iz)
    N3old(iz) = N3new(iz)
    N4old(iz) = N4new(iz)
    N6old(iz) = N6new(iz)
    N5old(iz) = N5new(iz)
    Neold(iz) = Nenew(iz)

    Nnoold(iz) = Nnonew(iz)

    U1old(iz) = U1new(iz)
    U2old(iz) = U2new(iz)
    U3old(iz) = U3new(iz)
    U4old(iz) = U4new(iz)
    U5old(iz) = U5new(iz)
    U6old(iz) = U6new(iz)
    Ueold(iz) = Uenew(iz)

    Unoold(iz) = Unonew(iz)

    T1old (iz) = T1new (iz)
    T1pold(iz) = T1pnew(iz)
    T1told(iz) = T1tnew(iz)
    T2old (iz) = T2new (iz)
    T2pold(iz) = T2pnew(iz)
    T2told(iz) = T2tnew(iz)
    T3old (iz) = T3new (iz)
    T3pold(iz) = T3pnew(iz)
    T3told(iz) = T3tnew(iz)
    T4old (iz) = T4new (iz)
    T4pold(iz) = T4pnew(iz)
    T4told(iz) = T4tnew(iz)
    T5old (iz) = T5new (iz)
    T5pold(iz) = T5pnew(iz)
    T5told(iz) = T5tnew(iz)
    T6old (iz) = T6new (iz)
    T6pold(iz) = T6pnew(iz)
    T6told(iz) = T6tnew(iz)
    Teold (iz) = Tenew (iz)
    Tepold(iz) = Tepnew(iz)
    Tetold(iz) = Tetnew(iz)

    q1old(iz) = q1new(iz)
    q2old(iz) = q2new(iz)
    q3old(iz) = q3new(iz)
    q4old(iz) = q4new(iz)
    q5old(iz) = q5new(iz)
    q6old(iz) = q6new(iz)
    qeold(iz) = qenew(iz)

  enddo


  dt_max=max(dt_max,vmax*rln(iz))
  dt_max=5.*dt_max/Ro

  call pas_de_temps(iyd_TU,temps_TU,dt,postint,dto,postinto)
  deltat=dt/to
  deltat_2=.5d0*deltat

  if (temps_ecoule.ge.temps_ecriture) then

    temps_ecriture=temps_ecriture+step_sortie
    nrec=nrec+1

    buffer( 1)=nb_alt
    buffer( 2)=ncol_o
    buffer( 3)=iannee
    buffer( 4)=imois
    buffer( 5)=ijour
    buffer( 6)=iheure
    buffer( 7)=iminute
    buffer( 8)=seconde
    buffer( 9)=intpas
    buffer(10)=longeo
    buffer(11)=latgeo
    buffer(12)=lonmag
    buffer(13)=latmag
    buffer(14)=tmag
    buffer(15)=f107(2)
    buffer(16)=f107(3)
    buffer(17)=ap(2)
    buffer(18)=kptime
    buffer(19)=dTinf
    buffer(20)=dUinf
    buffer(21)=cofo
    buffer(22)=cofh
    buffer(23)=cofn
    buffer(24)=chi0
    buffer(25)=Fe0
    buffer(26)=Ee0
    buffer(27)=Fi0
    buffer(28)=Ei0
    buffer(29)=Bmag
    buffer(30)=dipangle
    buffer(31)=Enord
    buffer(32)=Eest
    buffer(33)=vperpnord
    buffer(34)=vperpest
    buffer(35)=vhorizon
    buffer(36)=vpara
    buffer(37)=iapprox_o
    buffer(38)=ddp
    buffer(39)=Jtop
    buffer(60)=1
    do iz=1,nb_alt
      ipos=(iz+1)*ncol
      buffer(ipos+ipos_o%z)     = alt(iz)
      buffer(ipos+ipos_o%n1)    = N1new(iz)*N_o
      buffer(ipos+ipos_o%n2)    = N2new(iz)*N_o
      buffer(ipos+ipos_o%n3)    = N3new(iz)*N_o
      buffer(ipos+ipos_o%n4)    = N4new(iz)*N_o
      buffer(ipos+ipos_o%n5)    = N5new(iz)*N_o
      buffer(ipos+ipos_o%n6)    = N6new(iz)*N_o
      buffer(ipos+ipos_o%u1)    = U1new(iz)*C_o(1)
      buffer(ipos+ipos_o%u2)    = U2new(iz)*C_o(2)
      buffer(ipos+ipos_o%u3)    = U3new(iz)*C_o(3)
      buffer(ipos+ipos_o%u4)    = U4new(iz)*C_o(4)
      buffer(ipos+ipos_o%u5)    = U5new(iz)*C_o(5)
      buffer(ipos+ipos_o%u6)    = U6new(iz)*C_o(6)
      buffer(ipos+ipos_o%ue)    = Uenew(iz)*C_o(indi_e)
      buffer(ipos+ipos_o%t1p)   = T1pnew(iz)*T_o
      buffer(ipos+ipos_o%t1t)   = T1tnew(iz)*T_o
      buffer(ipos+ipos_o%t2p)   = T2pnew(iz)*T_o
      buffer(ipos+ipos_o%t2t)   = T2tnew(iz)*T_o
      buffer(ipos+ipos_o%t3p)   = T3pnew(iz)*T_o
      buffer(ipos+ipos_o%t3t)   = T3tnew(iz)*T_o
      buffer(ipos+ipos_o%t4p)   = T4pnew(iz)*T_o
      buffer(ipos+ipos_o%t4t)   = T4tnew(iz)*T_o
      buffer(ipos+ipos_o%t5p)   = T5pnew(iz)*T_o
      buffer(ipos+ipos_o%t5t)   = T5tnew(iz)*T_o
      buffer(ipos+ipos_o%t6p)   = T6pnew(iz)*T_o
      buffer(ipos+ipos_o%t6t)   = T6tnew(iz)*T_o
      buffer(ipos+ipos_o%tep)   = Tepnew(iz)*T_o
      buffer(ipos+ipos_o%tet)   = Tetnew(iz)*T_o
      buffer(ipos+ipos_o%q1)    = q1new(iz)*q_o(1)
      buffer(ipos+ipos_o%q2)    = q2new(iz)*q_o(2)
      buffer(ipos+ipos_o%q3)    = q3new(iz)*q_o(3)
      buffer(ipos+ipos_o%q4)    = q4new(iz)*q_o(4)
      buffer(ipos+ipos_o%q5)    = q5new(iz)*q_o(5)
      buffer(ipos+ipos_o%q6)    = q6new(iz)*q_o(6)
      buffer(ipos+ipos_o%qe)    = qenew(iz)*q_o(indi_e)
      buffer(ipos+ipos_o%nno)   = Nnonew(iz)*N_o
      buffer(ipos+ipos_o%uno)   = Unonew(iz)*Co(indi_NO)
      buffer(ipos+ipos_o%po)    = Po(iz)*N_o/to
      buffer(ipos+ipos_o%ph)    = Ph(iz)*N_o/to
      buffer(ipos+ipos_o%pn)    = Pn(iz)*N_o/to
      buffer(ipos+ipos_o%pn2)   = Pn2(iz)*N_o/to
      buffer(ipos+ipos_o%po2)   = Po2(iz)*N_o/to
      buffer(ipos+ipos_o%heat)  = Heat(iz)
      buffer(ipos+ipos_o%no)    = No(iz)
      buffer(ipos+ipos_o%nh)    = Nh(iz)
      buffer(ipos+ipos_o%nn)    = Nn(iz)
      buffer(ipos+ipos_o%nn2)   = Nn2(iz)
      buffer(ipos+ipos_o%no2)   = No2(iz)
      buffer(ipos+ipos_o%tn)    = Tn(iz)
      buffer(ipos+ipos_o%un)    = Un(iz)
      buffer(ipos+ipos_o%vn)    = Vn(iz)
      buffer(ipos+ipos_o%wn)    = Wn(iz)
      buffer(ipos+ipos_o%nes)   = Nes(iz)
      buffer(ipos+ipos_o%jes)   = Jes(iz)
      buffer(ipos+ipos_o%tes)   = Tes(iz)
      buffer(ipos+ipos_o%qes)   = Qes(iz)
    enddo

    call write_out(buffer,longbuf,longrec,nrectemps)

  endif

1000	   format(' temps ecoule:  ',i8,' tube numero: ',i3)
     write(*,1000) int(temps),itube

  temps_ecoule  =temps_ecoule + dt
  temps_TU     = temps_TU     + dt

  if (temps_TU>=86400.d0) then
    iyd_TU=iyd_TU+1
    temps_TU=temps_TU-86400.d0
  endif
  iannee=iyd_TU/1000
  call jour_mois(iyd_TU,ijour,imois)
  iheure=int(temps_TU/3600.)
  iminute=int(temps_TU/60.-iheure*60.)
  seconde=temps_TU-iheure*3600.-iminute*60.
  if (seconde.ge.60.) then
    iminute=iminute+1
    seconde=seconde-60.
  endif
  heure_TU=iheure+iminute/60.
  call lec_indices(iannee,imois,ijour,heure_TU,ap,f107)

  kp=ap2kp(ap(2))
  ikp=kp*3



  if (temps_ecoule.ge.temps_atmos) then
    temps_atmos=temps_atmos+step_atmos
    flag_atmos=.true.
  endif

  if (temps_ecoule.ge.temps_cinetic) then
    temps_cinetic=temps_cinetic+step_cinetic
    flag_cinetic=.true.
  endif


  !--------------------------------------!	
  ! 	Appel de l'atmosphere neutre     !
  !--------------------------------------!

  if (flag_atmos) then
    flag_atmos=.false.
    do i=1,nb_alt
      if (z(iz).le.zcira) then

    ! 	    ... MSIS90 + HWM93
        call msis90(iyd,temps,z(iz),latgeo,longeo,f107,ap,48,d,t,w)
    	
        Nh (iz)=d(7)*cofh
        No (iz)=d(2)*cofo
        No2(iz)=d(4)*cofo2
        Nn2(iz)=d(3)*cofn2
        Nn (iz)=d(8)*cofn
        Tn (iz)=t(2)

        vent(1,iz)=w(1)
        vent(2,iz)=w(2)

      else

        Tn(iz)=Tn(iz-1)
        fact=exp(gr(iz)/Tn(iz)*(z(iz-1)-z(iz)))
        Nh (iz) = Nh (iz-1)*fact
        No (iz) = No (iz-1)*fact**massn(indn_O)
        No2(iz) = No2(iz-1)*fact**massn(indn_O)
        Nn2(iz) = Nn2(iz-1)*fact**massn(indn_O)
        Nn (iz) = Nn (iz-1)*fact**massn(indn_O)
        Tn (iz) = Tn (iz-1)

        vent(1,iz)=vent(1,iz-1)
        vent(2,iz)=vent(2,iz-1)

      endif
    enddo

  !--------------------------------------!	
  ! 	Appel du transport cinetique     !
  !--------------------------------------!

  if (flag_cinetic) then
    flag_cinetic=.false.
    call transelec(npt,iyd,temps,latgeo,longeo,f107,ap,chi,				&
    		   Nenew,Tenew,T1new,nb_alt,Nh,No,No2,Nn2,Nn,Tn,indlim,jpreci,		&
                   N_o,T_o,kiappel,z,Heat,Ph,Po,Po2,Pn2,Pn,Nes,Jes,Tes,qes)
  endif




  !------------------------------------------------!	
  ! 	premiere demi-boucle temporelle a dt/2     !
  !------------------------------------------------!

  if (flag_conv) then
    call convec(iyd_TU,temps_TU,kp,dlongeo,dlatgeo,dlonmag,dlatmag,dlonref,.5*dt,pot,flag_pot)
    longeo=dlongeo
    latgeo=dlatgeo
    vtrans=vpara/Co(indi_H)
  else
    call convec(iyd_TU,temps_TU,kp,dlongeo,dlatgeo,dlonmag,dlatmag,dlonref,dx0,pot,flag_pot)
    longeo=dlongeo
    latgeo=dlatgeo
    vtrans=0.
  endif

  if (vparaB/=0.) vtrans=vparaB/Co(indi_H)

  do i=2,nb_alt
    Velno(iz) = Co(indn_NO)*.5*(Unonew(iz)+Unonew(iz-1))

    Vel1(iz) = Cij(indi_O ,indi_H)*.5*(U1new(iz)+U1new(iz-1))+vtrans
    Vel2(iz) =                     .5*(U2new(iz)+U2new(iz-1))+vtrans
    Vel3(iz) = Cij(indi_N ,indi_H)*.5*(U3new(iz)+U3new(iz-1))+vtrans
    Vel4(iz) = Cij(indi_N2,indi_H)*.5*(U4new(iz)+U4new(iz-1))+vtrans
    Vel5(iz) = Cij(indi_NO,indi_H)*.5*(U5new(iz)+U5new(iz-1))+vtrans
    Vel6(iz) = Cij(indi_O2,indi_H)*.5*(U6new(iz)+U6new(iz-1))+vtrans
    Vele(iz) = Cij(indi_e ,indi_H)*.5*(Uenew(iz)+Uenew(iz-1))+vtrans
  enddo

  Velno (1) = Velno(2)
  Velno(np) = Velno(nb_alt)

  Vel1 (1) = Vel1(2)
  Vel1(np) = Vel1(nb_alt)
  Vel2 (1) = Vel2(2)
  Vel2(np) = Vel2(nb_alt)
  Vel3 (1) = Vel3(2)
  Vel3(np) = Vel3(nb_alt)
  Vel4 (1) = Vel4(2)
  Vel4(np) = Vel4(nb_alt)
  Vel5 (1) = Vel5(2)
  Vel5(np) = Vel5(nb_alt)
  Vel6 (1) = Vel6(2)
  Vel6(np) = Vel6(nb_alt)
  Vele (1) = Vele(2)
  Vele(np) = Vele(nb_alt)

  do iz=1,np
    Velno_2(iz) = .5*Velno(iz)

    Vel1_2(iz) = .5*(Vel1(iz)+vtrans)
    Vel2_2(iz) = .5*(Vel2(iz)+vtrans)
    Vel3_2(iz) = .5*(Vel3(iz)+vtrans)
    Vel4_2(iz) = .5*(Vel4(iz)+vtrans)
    Vel5_2(iz) = .5*(Vel5(iz)+vtrans)
    Vel6_2(iz) = .5*(Vel6(iz)+vtrans)
    Vele_2(iz) = .5*(Vele(iz)+vtrans)
  enddo



  !-------------------------------------!
  !       equations de continuite       !
  !-------------------------------------!

  !	ion O+

  call velocity(Vel1,Ipos1,Iposnp,deltat_2)
  lbc=sqrt(N1old(2)/N1old(3))
  rbc=bclimld(Radn(np),N1old,extra,nb_alt)
  call lcpfct(N1old,N1new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,1)

  if (isnant(N1new,nb_alt)) then
    print*,'probleme lors du calcul de N1new dans la boucle 1 apres traitement du transport'
    goto 246
  endif
	
  !	ion H+

  call velocity(Vel2,Ipos1,Iposnp,deltat_2)
  lbc=sqrt(N2old(2)/N2old(3))
  rbc=bclimld(Radn(np),N2old,extra,nb_alt)
  lbc=1.
  rbc=(N2old(nb_alt)/N2old(nb_alt-1))
  call lcpfct(N2old,N2new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,1)

  if (isnant(N2new,nb_alt)) then
    print*,'probleme lors du calcul de N2new dans la boucle 1 apres traitement du transport'
    goto 246
  endif
	
  !	ion N+

  call velocity(Vel3,Ipos1,Iposnp,deltat_2)
  lbc=sqrt(N3old(2)/N3old(3))
  rbc=bclimld(Radn(np),N3old,extra,nb_alt)
  call lcpfct(N3old,N3new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,1)

  if (isnant(N3new,nb_alt)) then
    print*,'probleme lors du calcul de N3new dans la boucle 1 apres traitement du transport'
    goto 246
  endif
	
  !	ion N2+

  call velocity(Vel4,Ipos1,Iposnp,deltat_2)
  lbc=sqrt(N4old(2)/N4old(3))
  rbc=bclimld(Radn(np),N4old,extra,nb_alt)
  call lcpfct(N4old,N4new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,1)

  if (isnant(N4new,nb_alt)) then
    print*,'probleme lors du calcul de N4new dans la boucle 1 apres traitement du transport'
    goto 246
  endif
	
  !	ion NO+

  call velocity(Vel5,Ipos1,Iposnp,deltat_2)
  lbc=sqrt(N5old(2)/N5old(3))
  rbc=bclimld(Radn(np),N5old,extra,nb_alt)
  call lcpfct(N5old,N5new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,1)

  if (isnant(N5new,nb_alt)) then
    print*,'probleme lors du calcul de N5new dans la boucle 1 apres traitement du transport'
    goto 246
  endif
	
  !	ion O2+

  call velocity(Vel6,Ipos1,Iposnp,deltat_2)
  lbc=sqrt(N6old(2)/N6old(3))
  rbc=bclimld(Radn(np),N6old,extra,nb_alt)
  call lcpfct(N6old,N6new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,1)

  if (isnant(N6new,nb_alt)) then
    print*,'probleme lors du calcul de N6new dans la boucle 1 apres traitement du transport'
    goto 246
  endif



  !-------------------------------------!
  !        equations d'impulsion        !
  !-------------------------------------!

  !	ion O+

  call velocity(Vel1_2,Ipos1,Iposnp,deltat_2)
	
  do iz=1,nb_alt

  !  C2a(iz)=-Cij(indi_O,indi_H)*Tepnew(iz)
  !  D2a(iz)=log(Nenew(iz)*Tepnew(iz))
  !  C2b(iz)=-Cij(indi_O,indi_H)*T1pnew(iz)*N1new(iz)/xn1(iz)
  !  D2b(iz)=log(xn1(iz)*T1pnew(iz))
    C2a(iz)=-Cij(indi_O,indi_H)/Nenew(iz)
    D2a(iz)=Nenew(iz)*Tepnew(iz)
    C2b(iz)=-Cij(indi_O,indi_H)/N1new(iz)
    D2b(iz)=N1new(iz)*T1pnew(iz)

  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  lbc=1.
  rbc=1.
  call lcpfct(U1old,U1new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  if (isnant(U1new,nb_alt)) then
    print*,'probleme lors du calcul de U1new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	ion H+

  call velocity(Vel2_2,Ipos1,Iposnp,deltat_2)
	
  do iz=1,nb_alt

  !  C2a(iz)=-Tepnew(iz)
  !  D2a(iz)=log(Nenew(iz)*Tepnew(iz))
  !  C2b(iz)=-*T2pnew(iz)*N2new(iz)/xn2(iz)
  !  D2b(iz)=log(xn2(iz)*T2pnew(iz))
    C2a(iz)=-1./Nenew(iz)
    D2a(iz)=Nenew(iz)*Tepnew(iz)
    C2b(iz)=-1./N2new(iz)
    D2b(iz)=N2new(iz)*T2pnew(iz)

  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  lbc=1.
  rbc=1.
  call lcpfct(U2old,U2new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  if (isnant(U2new,nb_alt)) then
    print*,'probleme lors du calcul de U2new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	ion N+

  call velocity(Vel3_2,Ipos1,Iposnp,deltat_2)
	
  do iz=1,nb_alt

  !  C2a(iz)=-Cij(indi_N,indi_H)*Tepnew(iz)
  !  D2a(iz)=log(Nenew(iz)*Tepnew(iz))
  !  C2b(iz)=-Cij(indi_N,indi_H)*T3pnew(iz)*N3new(iz)/xn3(iz)
  !  D2b(iz)=log(xn3(iz)*T3pnew(iz))
    C2a(iz)=-Cij(indi_N,indi_H)/Nenew(iz)
    D2a(iz)=Nenew(iz)*Tepnew(iz)
    C2b(iz)=-Cij(indi_N,indi_H)/N3new(iz)
    D2b(iz)=N3new(iz)*T3pnew(iz)

  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  lbc=1.
  rbc=1.
  call lcpfct(U3old,U3new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  if (isnant(U3new,nb_alt)) then
    print*,'probleme lors du calcul de U3new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	ion N2+

  call velocity(Vel4_2,Ipos1,Iposnp,deltat_2)
	
  do iz=1,nb_alt

  !  C2a(iz)=-Cij(indi_N2,indi_H)*Tepnew(iz)
  !  D2a(iz)=log(Nenew(iz)*Tepnew(iz))
  !  C2b(iz)=-Cij(indi_N2,indi_H)*T4pnew(iz)*N4new(iz)/xn4(iz)
  !  D2b(iz)=log(xn4(iz)*T4pnew(iz))
    C2a(iz)=-Cij(indi_N2,indi_H)/Nenew(iz)
    D2a(iz)=Nenew(iz)*Tepnew(iz)
    C2b(iz)=-Cij(indi_N2,indi_H)/N4new(iz)
    D2b(iz)=N4new(iz)*T4pnew(iz)

  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  lbc=1.
  rbc=1.
  call lcpfct(U4old,U4new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  if (isnant(U4new,nb_alt)) then
    print*,'probleme lors du calcul de U4new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	ion NO+

  call velocity(Vel5_2,Ipos1,Iposnp,deltat_2)
	
  do iz=1,nb_alt

  !  C2a(iz)=-Cij(indi_O,indi_H)*Tepnew(iz)
  !  D2a(iz)=log(Nenew(iz)*Tepnew(iz))
  !  C2b(iz)=-Cij(indi_NO,indi_H)*T5pnew(iz)*N5new(iz)/xn5(iz)
  !  D2b(iz)=log(xn5(iz)*T5pnew(iz))
    C2a(iz)=-Cij(indi_NO,indi_H)/Nenew(iz)
    D2a(iz)=Nenew(iz)*Tepnew(iz)
    C2b(iz)=-Cij(indi_NO,indi_H)/N5new(iz)
    D2b(iz)=N5new(iz)*T5pnew(iz)

  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  lbc=1.
  rbc=1.
  call lcpfct(U5old,U5new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  if (isnant(U5new,nb_alt)) then
    print*,'probleme lors du calcul de U5new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	ion O2+

  call velocity(Vel6_2,Ipos1,Iposnp,deltat_2)
	
  do iz=1,nb_alt

  !  C2a(iz)=-Cij(indi_O2,indi_H)*Tepnew(iz)
  !  D2a(iz)=log(Nenew(iz)*Tepnew(iz))
  !  C2b(iz)=-Cij(indi_O2,indi_H)*T6pnew(iz)*N6new(iz)/xn6(iz)
  !  D2b(iz)=log(xn6(iz)*T6pnew(iz))
    C2a(iz)=-Cij(indi_O2,indi_H)/Nenew(iz)
    D2a(iz)=Nenew(iz)*Tepnew(iz)
    C2b(iz)=-Cij(indi_O2,indi_H)/N6new(iz)
    D2b(iz)=N6new(iz)*T6pnew(iz)

  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  lbc=1.
  rbc=1.
  call lcpfct(U6old,U6new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  if (isnant(U6new,nb_alt)) then
    print*,'probleme lors du calcul de U6new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !-------------------------------------!
  !     equations du flux de chaleur    !
  !-------------------------------------!

  !	ion O+
	
  call velocity(Vel1,Ipos1,Iposnp,deltat_2)

  do iz=1,nb_alt

    C2a(iz)=-2.2*Cij(indi_O,indi_H)*q1new(iz)
    D2a(iz)=U1new(iz)
    C2b(iz)=-Cij(indi_O,indi_H)*(r11_18*T1pnew(iz)+r8_9*T1tnew(iz))*N1new(iz)
    D2b(iz)=T1pnew(iz)
    C2c(iz)=-Cij(indi_O,indi_H)*(r17_9*T1pnew(iz)-r8_9*T1tnew(iz))*N1new(iz)
    D2c(iz)=T1tnew(iz)
    C2d(iz)=r4_9*Cij(indi_O,indi_H)*(T1pnew(iz)-T1tnew(iz))**2
    D2d(iz)=N1new(iz)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
  D2cl=(D2c(1)+D2c(2))/2.
  D2cr=ylimd(Radn(np),D2c,extra,nb_alt)
  D2dl=(D2d(1)+D2d(2))/2.
  D2dr=ylimd(Radn(np),D2d,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
  call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),q1old,extra,nb_alt)
  rbc=1.
  call lcpfct(q1old,q1new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  q1new(nb_alt)=max(0.,q1new(nb_alt))

  if (isnant(q1new,nb_alt)) then
    print*,'probleme lors du calcul de q1new dans la boucle 1 apres traitement du transport'
    goto 246
  endif


  !	ion H+
	
  call velocity(Vel2,Ipos1,Iposnp,deltat_2)

  do iz=1,nb_alt

    C2a(iz)=-2.2*q2new(iz)
    D2a(iz)=U2new(iz)
    C2b(iz)=-(r11_18*T2pnew(iz)+r8_9*T2tnew(iz))*N2new(iz)
    D2b(iz)=T1pnew(iz)
    C2c(iz)=-(r17_9*T2pnew(iz)-r8_9*T2tnew(iz))*N2new(iz)
    D2c(iz)=T2tnew(iz)
    C2d(iz)=r4_9*(T2pnew(iz)-T2tnew(iz))**2
    D2d(iz)=N2new(iz)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
  D2cl=(D2c(1)+D2c(2))/2.
  D2cr=ylimd(Radn(np),D2c,extra,nb_alt)
  D2dl=(D2d(1)+D2d(2))/2.
  D2dr=ylimd(Radn(np),D2d,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
  call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),q2old,extra,nb_alt)
  rbc=1.
  call lcpfct(q2old,q2new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  q2new(nb_alt)=max(0.,q2new(nb_alt))

  if (isnant(q2new,nb_alt)) then
    print*,'probleme lors du calcul de q2new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	ion N+
	
  call velocity(Vel3,Ipos1,Iposnp,deltat_2)

  do iz=1,nb_alt

    C2a(iz)=-2.2*Cij(indi_N,indi_H)*q3new(iz)
    D2a(iz)=U3new(iz)
    C2b(iz)=-Cij(indi_N,indi_H)*(r11_18*T3pnew(iz)+r8_9*T3tnew(iz))*N3new(iz)
    D2b(iz)=T3pnew(iz)
    C2c(iz)=-Cij(indi_N,indi_H)*(r17_9*T3pnew(iz)-r8_9*T3tnew(iz))*N3new(iz)
    D2c(iz)=T3tnew(iz)
    C2d(iz)=r4_9*Cij(indi_N,indi_H)*(T3pnew(iz)-T3tnew(iz))**2
    D2d(iz)=N3new(iz)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
  D2cl=(D2c(1)+D2c(2))/2.
  D2cr=ylimd(Radn(np),D2c,extra,nb_alt)
  D2dl=(D2d(1)+D2d(2))/2.
  D2dr=ylimd(Radn(np),D2d,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
  call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),q3old,extra,nb_alt)
  rbc=1.
  call lcpfct(q3old,q3new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  q3new(nb_alt)=max(0.,q3new(nb_alt))

  if (isnant(q3new,nb_alt)) then
    print*,'probleme lors du calcul de q3new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	ion N2+
	
  call velocity(Vel4,Ipos1,Iposnp,deltat_2)

  do iz=1,nb_alt

    C2a(iz)=-2.2*Cij(indi_N,indi_H)*q4new(iz)
    D2a(iz)=U4new(iz)
    C2b(iz)=-Cij(indi_N,indi_H)*(r11_18*T4pnew(iz)+r8_9*T4tnew(iz))*N4new(iz)
    D2b(iz)=T4pnew(iz)
    C2c(iz)=-Cij(indi_N,indi_H)*(r17_9*T4pnew(iz)-r8_9*T4tnew(iz))*N4new(iz)
    D2c(iz)=T4tnew(iz)
    C2d(iz)=r4_9*Cij(indi_N,indi_H)*(T4pnew(iz)-T4tnew(iz))**2
    D2d(iz)=N4new(iz)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
  D2cl=(D2c(1)+D2c(2))/2.
  D2cr=ylimd(Radn(np),D2c,extra,nb_alt)
  D2dl=(D2d(1)+D2d(2))/2.
  D2dr=ylimd(Radn(np),D2d,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
  call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),q4old,extra,nb_alt)
  rbc=1.
  call lcpfct(q4old,q4new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  q4new(nb_alt)=max(0.,q4new(nb_alt))

  if (isnant(q4new,nb_alt)) then
    print*,'probleme lors du calcul de q4new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	ion NO+
	
  call velocity(Vel5,Ipos1,Iposnp,deltat_2)

  do iz=1,nb_alt

    C2a(iz)=-2.2*Cij(indi_N,indi_H)*q5new(iz)
    D2a(iz)=U5new(iz)
    C2b(iz)=-Cij(indi_N,indi_H)*(r11_18*T5pnew(iz)+r8_9*T5tnew(iz))*N5new(iz)
    D2b(iz)=T5pnew(iz)
    C2c(iz)=-Cij(indi_N,indi_H)*(r17_9*T5pnew(iz)-r8_9*T5tnew(iz))*N5new(iz)
    D2c(iz)=T5tnew(iz)
    C2d(iz)=r4_9*Cij(indi_N,indi_H)*(T5pnew(iz)-T5tnew(iz))**2
    D2d(iz)=N5new(iz)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
  D2cl=(D2c(1)+D2c(2))/2.
  D2cr=ylimd(Radn(np),D2c,extra,nb_alt)
  D2dl=(D2d(1)+D2d(2))/2.
  D2dr=ylimd(Radn(np),D2d,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
  call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),q5old,extra,nb_alt)
  rbc=1.
  call lcpfct(q5old,q5new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  q5new(nb_alt)=max(0.,q5new(nb_alt))

  if (isnant(q5new,nb_alt)) then
    print*,'probleme lors du calcul de q5new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	ion O2+
	
  call velocity(Vel6,Ipos1,Iposnp,deltat_2)

  do iz=1,nb_alt

    C2a(iz)=-2.2*Cij(indi_N,indi_H)*q6new(iz)
    D2a(iz)=U6new(iz)
    C2b(iz)=-Cij(indi_N,indi_H)*(r11_18*T6pnew(iz)+r8_9*T6tnew(iz))*N6new(iz)
    D2b(iz)=T6pnew(iz)
    C2c(iz)=-Cij(indi_N,indi_H)*(r17_9*T6pnew(iz)-r8_9*T6tnew(iz))*N6new(iz)
    D2c(iz)=T6tnew(iz)
    C2d(iz)=r4_9*Cij(indi_N,indi_H)*(T6pnew(iz)-T6tnew(iz))**2
    D2d(iz)=N6new(iz)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
  D2cl=(D2c(1)+D2c(2))/2.
  D2cr=ylimd(Radn(np),D2c,extra,nb_alt)
  D2dl=(D2d(1)+D2d(2))/2.
  D2dr=ylimd(Radn(np),D2d,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
  call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),q6old,extra,nb_alt)
  rbc=1.
  call lcpfct(q6old,q6new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  q6new(nb_alt)=max(0.,q6new(nb_alt))

  if (isnant(q6new,nb_alt)) then
    print*,'probleme lors du calcul de q6new dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  !	electrons
	
  call velocity(Vele,Ipos1,Iposnp,deltat_2)

  do iz=1,nb_alt

    C2a(iz)=-2.2*Cij(indi_N,indi_H)*qenew(iz)
    D2a(iz)=Uenew(iz)
    C2b(iz)=-Cij(indi_N,indi_H)*(r11_18*Tepnew(iz)+r8_9*Tetnew(iz))*Nenew(iz)
    D2b(iz)=Tepnew(iz)
    C2c(iz)=-Cij(indi_N,indi_H)*(r17_9*Tepnew(iz)-r8_9*Tetnew(iz))*Nenew(iz)
    D2c(iz)=Tetnew(iz)
    C2d(iz)=r4_9*Cij(indi_N,indi_H)*(Tepnew(iz)-Tetnew(iz))**2
    D2d(iz)=Nenew(iz)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)
  D2cl=(D2c(1)+D2c(2))/2.
  D2cr=ylimd(Radn(np),D2c,extra,nb_alt)
  D2dl=(D2d(1)+D2d(2))/2.
  D2dr=ylimd(Radn(np),D2d,extra,nb_alt)
	
  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
  call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
  call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),qeold,extra,nb_alt)
  rbc=1.
  call lcpfct(qeold,qenew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  qenew(nb_alt)=Qetop

  if (isnant(qenew,nb_alt)) then
    print*,'probleme lors du calcul de qenew dans la boucle 1 apres traitement du transport'
    goto 246
  endif
  	
  !-------------------------------------!
  !      equations de temperatures      !
  !-------------------------------------!


  ! ion O+


  call velocity(Vel1,Ipos1,Iposnp,deltat_2)

  ! temperature parallele

  do i=1,nb_alt
   C2a(i)=-Cij(indi_O,indi_H)*T1pnew(i)
   D2a(i)=U1new(i)
   C2b(i)=-1.2*Cij(indi_O,indi_H)/N1new(i)
   D2b(i)=q1new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T1pold,extra,nb_alt)
  rbc=1.
  call lcpfct(T1pold,T1pnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do i=1,nb_alt
    T1pnew(i)=max(T1pnew(i),T_min)
  enddo

  if (isnant(T1pnew,nb_alt)) then
    print*,'probleme lors du calcul de T1pnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif



  ! temperature perpendiculaire

  do i=1,nb_alt
    C2a(i)=Cij(indi_O,indi_H)*T1tnew(i)
    D2a(i)=U1new(i)
    C2b(i)=-.4*Cij(indi_O,indi_H)/N1new(i)
    D2b(i)=q1new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T1told,extra,nb_alt)
  rbc=1.
  call lcpfct(T1told,T1tnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do  i=1,nb_alt
    T1tnew(i)=max(T1tnew(i),T_min)
    T1new(i)=(T1pnew(i)+2.*T1tnew(i))/3.
  enddo


  if (isnant(T1tnew,nb_alt)) then
    print*,'probleme lors du calcul de T1tnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  ! ion H+


  call velocity(Vel2,Ipos1,Iposnp,deltat_2)

  ! temperature parallele

  do i=1,nb_alt
   C2a(i)=-T2pnew(i)
   D2a(i)=U2new(i)
   C2b(i)=-1.2/N2new(i)
   D2b(i)=q2new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T2pold,extra,nb_alt)
  rbc=1.
  call lcpfct(T2pold,T2pnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do i=1,nb_alt
    T2pnew(i)=max(T2pnew(i),T_min)
  enddo

  if (isnant(T2pnew,nb_alt)) then
    print*,'probleme lors du calcul de T2pnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif



  ! temperature perpendiculaire

  do i=1,nb_alt
    C2a(i)=T2tnew(i)
    D2a(i)=U2new(i)
    C2b(i)=-.4/N2new(i)
    D2b(i)=q2new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T22told,extra,nb_alt)
  rbc=1.
  call lcpfct(T2told,T2tnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do  i=1,nb_alt
    T2tnew(i)=max(T2tnew(i),T_min)
    T2new(i)=(T2pnew(i)+2.*T2tnew(i))/3.
  enddo


  if (isnant(T2tnew,nb_alt)) then
    print*,'probleme lors du calcul de T2tnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  ! ion N+


  call velocity(Vel3,Ipos1,Iposnp,deltat_2)

  ! temperature parallele

  do i=1,nb_alt
   C2a(i)=-Cij(indi_N,indi_H)*T3pnew(i)
   D2a(i)=U3new(i)
   C2b(i)=-1.2*Cij(indi_N,indi_H)/N3new(i)
   D2b(i)=q3new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T3pold,extra,nb_alt)
  rbc=1.
  call lcpfct(T3pold,T3pnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do i=1,nb_alt
    T3pnew(i)=max(T3pnew(i),T_min)
  enddo

  if (isnant(T3pnew,nb_alt)) then
    print*,'probleme lors du calcul de T3pnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif



  ! temperature perpendiculaire

  do i=1,nb_alt
    C2a(i)=Cij(indi_N,indi_H)*T3tnew(i)
    D2a(i)=U3new(i)
    C2b(i)=-.4*Cij(indi_N,indi_H)/N3new(i)
    D2b(i)=q3new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T3told,extra,nb_alt)
  rbc=1.
  call lcpfct(T3told,T3tnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do  i=1,nb_alt
    T3tnew(i)=max(T3tnew(i),T_min)
    T3new(i)=(T3pnew(i)+2.*T3tnew(i))/3.
  enddo


  if (isnant(T3tnew,nb_alt)) then
    print*,'probleme lors du calcul de T3tnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  ! ion N2+


  call velocity(Vel4,Ipos1,Iposnp,deltat_2)

  ! temperature parallele

  do i=1,nb_alt
   C2a(i)=-Cij(indi_N2,indi_H)*T4pnew(i)
   D2a(i)=U4new(i)
   C2b(i)=-1.2*Cij(indi_N2,indi_H)/N4new(i)
   D2b(i)=q4new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T4pold,extra,nb_alt)
  rbc=1.
  call lcpfct(T4pold,T4pnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do i=1,nb_alt
    T4pnew(i)=max(T4pnew(i),T_min)
  enddo

  if (isnant(T4pnew,nb_alt)) then
    print*,'probleme lors du calcul de T4pnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif



  ! temperature perpendiculaire

  do i=1,nb_alt
    C2a(i)=Cij(indi_N2,indi_H)*T4tnew(i)
    D2a(i)=U4new(i)
    C2b(i)=-.4*Cij(indi_N2,indi_H)/N4new(i)
    D2b(i)=q4new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T4told,extra,nb_alt)
  rbc=1.
  call lcpfct(T4told,T4tnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do  i=1,nb_alt
    T4tnew(i)=max(T4tnew(i),T_min)
    T4new(i)=(T4pnew(i)+2.*T4tnew(i))/3.
  enddo


  if (isnant(T4tnew,nb_alt)) then
    print*,'probleme lors du calcul de T4tnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  ! ion NO+


  call velocity(Vel5,Ipos1,Iposnp,deltat_2)

  ! temperature parallele

  do i=1,nb_alt
   C2a(i)=-Cij(indi_NO,indi_H)*T5pnew(i)
   D2a(i)=U5new(i)
   C2b(i)=-1.2*Cij(indi_NO,indi_H)/N5new(i)
   D2b(i)=q5new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T5pold,extra,nb_alt)
  rbc=1.
  call lcpfct(T5pold,T5pnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do i=1,nb_alt
    T5pnew(i)=max(T5pnew(i),T_min)
  enddo

  if (isnant(T5pnew,nb_alt)) then
    print*,'probleme lors du calcul de T5pnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif



  ! temperature perpendiculaire

  do i=1,nb_alt
    C2a(i)=Cij(indi_NO,indi_H)*T5tnew(i)
    D2a(i)=U5new(i)
    C2b(i)=-.4*Cij(indi_NO,indi_H)/N5new(i)
    D2b(i)=q5new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T5told,extra,nb_alt)
  rbc=1.
  call lcpfct(T5told,T5tnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do  i=1,nb_alt
    T5tnew(i)=max(T5tnew(i),T_min)
    T5new(i)=(T5pnew(i)+2.*T5tnew(i))/3.
  enddo


  if (isnant(T5tnew,nb_alt)) then
    print*,'probleme lors du calcul de T5tnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  ! ion O2+


  call velocity(Vel6,Ipos1,Iposnp,deltat_2)

  ! temperature parallele

  do i=1,nb_alt
   C2a(i)=-Cij(indi_O2,indi_H)*T6pnew(i)
   D2a(i)=U6new(i)
   C2b(i)=-1.2*Cij(indi_O2,indi_H)/N6new(i)
   D2b(i)=q6new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T6pold,extra,nb_alt)
  rbc=1.
  call lcpfct(T6pold,T6pnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do i=1,nb_alt
    T6pnew(i)=max(T6pnew(i),T_min)
  enddo

  if (isnant(T6pnew,nb_alt)) then
    print*,'probleme lors du calcul de T6pnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif



  ! temperature perpendiculaire

  do i=1,nb_alt
    C2a(i)=Cij(indi_O2,indi_H)*T6tnew(i)
    D2a(i)=U6new(i)
    C2b(i)=-.4*Cij(indi_O2,indi_H)/N6new(i)
    D2b(i)=q6new(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),T6told,extra,nb_alt)
  rbc=1.
  call lcpfct(T6told,T6tnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do  i=1,nb_alt
    T6tnew(i)=max(T6tnew(i),T_min)
    T6new(i)=(T6pnew(i)+2.*T6tnew(i))/3.
  enddo


  if (isnant(T6tnew,nb_alt)) then
    print*,'probleme lors du calcul de T6tnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif

  ! electrons


  call velocity(Vele,Ipos1,Iposnp,deltat_2)

  ! temperature parallele

  do i=1,nb_alt
   C2a(i)=-Cij(indi_O,indi_H)*Tepnew(i)
   D2a(i)=Uenew(i)
   C2b(i)=-1.2*Cij(indi_O,indi_H)/Nenew(i)
   D2b(i)=qenew(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),Tepold,extra,nb_alt)
  rbc=1.
  call lcpfct(Tepold,Tepnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do i=1,nb_alt
    Tepnew(i)=max(Tepnew(i),T_min)
  enddo

  if (isnant(Tepnew,nb_alt)) then
    print*,'probleme lors du calcul de Tepnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif



  ! temperature perpendiculaire

  do i=1,nb_alt
    C2a(i)=Cij(indi_O,indi_H)*Tetnew(i)
    D2a(i)=Uenew(i)
    C2b(i)=-.4*Cij(indi_O,indi_H)/Nenew(i)
    D2b(i)=qenew(i)
  enddo

  D2al=(D2a(1)+D2a(2))/2.
  D2ar=ylimd(Radn(np),D2a,extra,nb_alt)
  D2bl=(D2b(1)+D2b(2))/2.
  D2br=ylimd(Radn(np),D2b,extra,nb_alt)

  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)

  lbc=1.
  rbc=bclimd(Radn(nb_alt+1),Tetold,extra,nb_alt)
  rbc=1.
  call lcpfct(Tetold,Tetnew,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

  do  i=1,nb_alt
    Tetnew(i)=max(Tetnew(i),T_min)
    Tenew(i)=(Tepnew(i)+2.*Tetnew(i))/3.
  enddo


  if (isnant(Tetnew,nb_alt)) then
    print*,'probleme lors du calcul de Tetnew dans la boucle 1 apres traitement du transport'
    goto 246
  endif


  do iz=1,nb_alt


    !************************************************************************
    !*                          QUELQUES CALCULS UTILES                     *
    !************************************************************************

    dens_n(indi_O )=No(iz)
    dens_n(indi_H )=Nh(iz)
    dens_n(indi_N )=Nn(iz)
    dens_n(indi_N2)=Nn2(iz)
    dens_n(indi_NO)=Nno(iz)
    dens_n(indi_O2)=No2(iz)
    dens_n(indi_O1)=No1(iz)

    dens_i(indi_O )=N1new(iz)
    dens_i(indi_H )=N2new(iz)
    dens_i(indi_N )=N3new(iz)
    dens_i(indi_N2)=N4new(iz)
    dens_i(indi_NO)=N5new(iz)
    dens_i(indi_O2)=N6new(iz)
    dens_i(indi_e )=Nenew(iz)

    do i=1,indi_e
      dens_i_1(i)=1./dens_i(i)
    enddo

    vitesse(indi_O )=U1new(iz)*Co(indi_O )
    vitesse(indi_H )=U2new(iz)*Co(indi_H )
    vitesse(indi_N )=U3new(iz)*Co(indi_N )
    vitesse(indi_N2)=U4new(iz)*Co(indi_N2)
    vitesse(indi_NO)=U5new(iz)*Co(indi_NO)
    vitesse(indi_O2)=U6new(iz)*Co(indi_O2)
    vitesse(indi_e )=Uenew(iz)*Co(indi_e )

    temperature(indi_O )=T1new(iz)*T_o
    temperature(indi_H )=T2new(iz)*T_o
    temperature(indi_N )=T3new(iz)*T_o
    temperature(indi_N2)=T4new(iz)*T_o
    temperature(indi_NO)=T5new(iz)*T_o
    temperature(indi_O2)=T6new(iz)*T_o
    temperature(indi_e )=Tenew(iz)*T_o

    do indi_ion=1,nb_ion
      friction(indi_ion)= (vitesse(indi_ion)-Un(iz))**2+(vperpnord-Vn(iz))**2+(vperpest-Wn(iz))**2
      do j_ion=1,nb_ion
        Tij_r(indi_ion,j_ion)= mirij(indi_ion,j_ion)*temperature(j_ion)+mirij(j_ion,indi_ion)*temperature(indi_ion)
        Tij_r_1(indi_ion,j_ion)=1./Tij_r(indi_ion,j_ion)
        Tij_r_3(indi_ion,j_ion)=Tij_r_1(indi_ion,j_ion)*sqrt(Tij_r_1(indi_ion,j_ion))
        Tij(indi_ion,j_ion)=Tij_r(indi_ion,j_ion)/T_o
      end do

      do j_neutre=1,nb_spn
        Tr(indi_ion,j_neutre)=(temperature(indi_ion)+Tn(iz))/2.
        Tin(indi_ion,j_neutre)=mirin(indi_ion,j_neutre)*Tn(iz)+mnrin(indi_ion,j_neutre)*temperature(indi_ion)
      end do
    end do

    dTen=Tn(iz)-temperature(indi_e)
    if (abs(dTen)<1.e-2) dTen=sign(1.e-2,dTen)
    temperature(indi_e)=Tn(iz)-dTen
    Tn_1=1./Tn(iz)
    Te_1=1./temperature(indi_e)
    Teh=sqrt(temperature(indi_e))
    Teh_3=temperature(indi_e)*Teh



!    nu_omega=(nukN2(i)+nukO2(i)+nukO(i)+nulN2(i)+nulO2(i)+nulO(i)+numN2(i)+numO2(i)+numO(i))/omega(i)*omz
!    coef_conv=1./(1.+nu_omega**2)
!    Tperp(i)=cofterp*coef_conv*Vm_2(i)



    !***********************************************************************
    !                            REACTIONS CHIMIQUES                       *
    !***********************************************************************
	

  !********************
  ! O+ + H -> H+ + O  *
  !********************

    rk(1)=rko(1)*sqrt(Tin(indi_H,indn_O))*dens_n(indn_H)

  !*********************
  ! O+ + N2 -> NO+ + N *		St-Maurice J.-P. et P.J. Laneville
  !*********************

    Teff=Tin(indi_O,indn_N2)
    Teff2=Teff*Teff
    Teff3=Teff*Teff2
    Teff4=Teff3*Teff
    if (Teff <= 3725.) then
      rk(2) =  1.71676e-18 		&
             - 2.39978e-21*Teff     	&
             + 1.48088e-24*Teff2	&
             - 3.43783e-28*Teff3  	&
             + 7.89577e-32*Teff4
    else
      rk(2) =- 1.52489e-17		&
             + 2.55704e-21*Teff		&
             + 1.32293e-24*Teff2	&
             - 4.84659e-29*Teff3  	&
             + 5.77477e-34*Teff4
    endif
    rk(2)=rk(2)*dens_n(ind_N2)


  !*********************
  ! O+ + NO -> NO+ + O *		St-Maurice J.-P. et P.J. Laneville
  !*********************

    Teff=Tin(indi_O,indn_NO)
    Teff2=Teff*Teff
    Teff3=Teff*Teff2
    Teff4=Teff3*Teff
    if (Teff <= 3800.) then
      rk(3) =  6.40408e-18 		&
             - 4.46293e-22*Teff     	&
             + 8.50114e-25*Teff2	&
             - 1.15374e-28*Teff3  	&
             + 8.17746e-33*Teff4
    else
      rk(3) =- 7.48318e-19		&
             + 7.71673e-22*Teff		&
             + 3.41289e-25*Teff2	&
             - 9.83096e-30*Teff3  	&
             + 9.58846e-35*Teff4
    endif
    rk(3)=rk(3)*dens_n(ind_NO)

  !*********************
  ! O+ + O2 -> O2+ + O *		St-Maurice J.-P. et P.J. Laneville
  !*********************

    Teff=Tin(indi_O,indn_O2)
    Teff2=Teff*Teff
    Teff3=Teff*Teff2
    if (Teff <= 4800.) then
      Teff4=Teff3*Teff
      rk(4)  =  2.78932e-17		&
              - 2.30871e-20*Teff     	&
              + 9.64093e-24*Teff2	&
              - 1.28611e-27*Teff3  	&
              + 6.26046e-32*Teff4
    else
      rk(4)  =- 1.74046e-17		&
              + 1.00776e-20*Teff	&
              - 2.65793e-26*Teff2	&
              - 1.49035e-30*Teff3
    endif
    rk(4)=rk(4)*dens_n(ind_O2)

  !********************
  ! H+ + O -> O+ + H  *
  !********************

    rk(5)=rko(5)*dens_n(ind_O)

  !********************
  ! N+ + O -> O+ + N  *
  !********************

    rk(6)=rko(6)*dens_n(ind_O)

  !********************
  ! N+ + H -> H+ + N  *
  !********************

    rk(7)=rko(7)*dens_n(ind_H)

  !**********************
  ! N+ + NO -> NO+ + N  *
  !**********************

    rk(8)=rko(8)*dens_n(ind_NO)

  !**********************
  ! N+ + O2 -> O+ + NO  *
  !**********************

    rk(9)=rko(9)*dens_n(ind_O2)

  !**********************
  ! N+ + O2 -> NO+ + O  *
  !**********************

    rk(10)=rko(10)*dens_n(ind_O2)

  !**********************
  ! N+ + O2 -> O2+ + N  *
  !**********************

    rk(11)=rko(11)*dens_n(ind_O2)

  !**********************
  ! N2+ + O -> O+ + N2  *
  !**********************

    Teff=Tin(indi_N2,indn_O)
!    if (Teff <= 1500.) then par raison de continuite
    if (Teff <= 1480.) then
      rk(12) = 3.71e-17*Teff**(-.23)
    else
      rk(12) = 3.47e-19*Teff**(.41)
    endif
    rk(12)=rk(12)*dens_n(ind_O)

  !**********************
  ! N2+ + O -> NO+ + N  *
  !**********************

  !  if (Teff <= 1500.) then  par raison de continuite
    if (Teff <= 1410.) then
      rk(13) = 1.72e-15*Teff**(-.44)
    else
      rk(13) = 1.66e-17*Teff**(.2)
    endif
    rk(13)=rk(13)*dens_n(ind_O)

  !************************
  ! N2+ + NO -> NO+ + N2  *
  !************************

    rk(14)=rko(14)*dens_n(ind_NO)

  !************************
  ! N2+ + O2 -> O2+ + N2  *
  !************************

    rk(15)=rko(15)*dens_n(ind_O2)

  !*********************
  ! N2+ + e- -> N + N  *
  !*********************

    rk(16)=rko(16)*temperature(indi_e)**(-.39)*dens_i(indi_e)

  !*********************
  ! NO+ + e- -> O + N  *
  !*********************

    rk(17)=rko(17)*temperature(indi_e)**(-.85)*dens_i(indi_e)

  !**********************
  ! O2+ + N -> NO+ + O  *
  !**********************

    rk(18)=rko(18)*dens_n(ind_O2)

  !************************
  ! O2+ + N2 -> NO+ + NO  *
  !************************

    rk(19)=rko(19)*dens_n(ind_O2)

  !************************
  ! O2+ + NO -> NO+ + O2  *
  !************************

    rk(20)=rko(20)*dens_n(indn_NO)

  !*********************
  ! O2+ + e- -> O + O  *
  !*********************

    rk(21)=rko(21)*temperature(indi_e)**(-.85)*dens_i(indi_e)

      	
  ! ion O+

    y0  (indi_O )         = N1old(iz)
    y1  (indi_O )         = (N1new(iz)-N1old(iz))/deltat_2+Po(iz)

    chim(indi_O ,indi_O ) =-rk(1)				&		! O+  + H  -> H+  + O
                           -rk(2)				&		! O+  + N2 -> NO+ + N
                           -rk(3)				&		! O+  + NO -> NO+ + O
                           -rk(4)				&		! O+  + O2 -> O2+ + O
	 	           -3.*Cij(indi_O ,indi_H)*U1new(iz)*alt_geo_1(iz)	! divergence du tube
		
    chim(indi_O ,indi_H ) = rk(5)						! H+  + O  -> O+  + H
		
    chim(indi_O ,indi_N ) = rk(6)				&		! N+  + O  -> O+  + N
                           +rk(9)						! N+  + O2 -> O+  + NO

    chim(indi_O ,indi_N2) = rk(12)						! N2+ + O  -> O+  + N2

    chim(indi_O ,indi_NO) = 0.d0

    chim(indi_O ,indi_O2) = 0.d0

  ! ion H+

    y0  (indi_H )	  = N2old(iz)
    y1  (indi_H )	  = (N2new(iz)-N2old(iz))/deltat_2+Ph(iz)

    chim(indi_H ,indi_O ) = rk(1)						! O+  + H  -> H+  + O

    chim(indi_H ,indi_H ) =-rk(5)				&		! H+  + O  -> O+  + H	
	                   -3.*U2new(iz)*alt_geo_1(iz)				! divergence du tube

    chim(indi_H ,indi_N ) = rk(7)						! N+  + H  -> H+  + N

    chim(indi_H ,indi_N2) = 0.d0

    chim(indi_H ,indi_NO) = 0.d0

    chim(indi_H ,indi_O2) = 0.d0

  ! ion N+

    y0  (indi_N )	  = N3old(iz)
    y1  (indi_N )	  = (N3new(iz)-N3old(iz))/deltat_2+0.21*Pn2(iz)

    chim(indi_N ,indi_O ) = 0.d0
    chim(indi_N ,indi_H ) = 0.d0
    chim(indi_N ,indi_N ) = rk(6)				&		! N+  + O  -> O+  + N
                           -rk(7)				&		! N+  + H  -> H+  + N
	                   -rk(8)				&		! N+  + NO -> NO+ + N
	                   -rk(9)				&		! N+  + O2 -> O+  + NO
			   -rk(10)				&		! N+  + O2 -> NO+ + O
			   -rk(11)				&		! N+  + O2 -> O2+ + N
			   -3.*Cij(indi_N ,indi_H)*U3new(iz)*alt_geo_1(iz)	! divergence du tube

    chim(indi_N ,indi_N2) = 0.d0

    chim(indi_N ,indi_NO) = 0.d0

    chim(indi_N ,indi_O2) = 0.d0


  ! ion N2+

    y0  (indi_N2)	  = N4old(iz)
    y1  (indi_N2)	  = (N4new(iz)-N4old(iz))/deltat_2+0.79*Pn2(iz)

    chim(indi_N2,indi_O ) = 0.d0

    chim(indi_N2,indi_H ) = 0.d0

    chim(indi_N2,indi_N ) = 0.d0

    chim(indi_N2,indi_N2) =-rk(12)				&		! N2+ + O  -> O+  + N2
                           -rk(13)				&		! N2+ + O  -> NO+ + N
                           -rk(14)				&		! N2+ + NO -> NO+ + N2
                           -rk(15)				&		! N2+ + O2 -> O2+ + N2
                           -rk(16)				&		! N2+ + e- -> N   + N
			   -3.*Cij(indi_N2,indi_H)*U3new(iz)*alt_geo_1(iz)	! divergence du tube

    chim(indi_N2,indi_NO) = 0.d0

    chim(indi_N2,indi_O2) = 0.d0

  ! ion NO+

    y0  (indi_NO)	  = N5old(iz)
    y1  (indi_NO)	  = (N5new(iz)-N5old(iz))/deltat_2+rko(22)*dens_n(indn_NO)

    chim(indi_NO,indi_O ) = rk(2)						! O+  + N2 -> NO+ + N

    chim(indi_NO,indi_H ) = 0.d0

    chim(indi_NO,indi_N ) = rk(8)						! N+  + NO -> NO+ + N

    chim(indi_NO,indi_N2) = rk(13)				&		! N2+ + O  -> NO+ + N
                           +rk(14)						! N2+ + NO -> NO+ + N2

    chim(indi_NO,indi_NO) = -rk(17)				&		! NO+ + e- -> O   + N
			   -3.*Cij(indi_NO,indi_H)*Umnew(iz)*alt_geo_1(iz)	! divergence du tube

    chim(indi_NO,indi_O2) = rk(18)				&		! O2+ + N  -> NO+ + O
                           +rk(19)				&		! O2+ + N2 -> NO+ + NO
                           +rk(20)						! O2+ + NO -> NO+ + O2

  ! ion O2+

    y0  (indi_O2)	  = N6old(i)
    y1  (indi_O2)	  = (N6new(iz)-N6old(iz))/deltat_2+Po2(iz)

    chim(indi_O2,indi_O ) = rk(4)						! O+  + O2 -> O2+ + O

    chim(indi_O2,indi_H ) = 0.d0

    chim(indi_O2,indi_N ) = rk(11)						! N+  + O2 -> O2+ + N

    chim(indi_O2,indi_N2) = rk(15)						! N2+ + O2 -> O2+ + N2

    chim(indi_O2,indi_NO) = 0.d0

    chim(indi_O2,indi_O2) =-rk(18)				&		! O2+ + N  -> NO+ + O
                           -rk(19)				&		! O2+ + N2 -> NO+ + NO
                           -rk(20)				&		! O2+ + NO -> NO+ + O2
                           -rk(21)				&		! O2+ + e- -> O   + O
			   -3.*Cij(indi_O2,indi_H)*U3new(iz)*alt_geo_1(iz)	! divergence du tube

    call solve_ODE(nb_ion,y0,y1,chim,mat,deltat_2)

    N1new(iz)=max(y0(indi_O ),r_min)
    N2new(iz)=max(y0(indi_H ),r_min)
    N3new(iz)=max(y0(indi_N ),r_min)
    N4new(iz)=max(y0(indi_N2),r_min)
    N5new(iz)=max(y0(indi_NO),r_min)
    N6new(iz)=max(y0(indi_O2),r_min)
    Nenew(iz)=N1new(iz)+N2new(iz)+N3new(iz)+N4new(iz)+N5new(iz)+N6new(iz)


    !************************************************************************
    !*                          FREQUENCES DE COLLISION                     *
    !************************************************************************


    coefnu(indi_H ,indi_H )  = (1.-.083*log10(Tr(indi_H,indi_H)))**2*sqrt(Tr(indi_H,indi_H))
    coefnu(indi_O ,indi_O )  = (1.-.064*log10(Tr(indi_O,indi_O)))**2*sqrt(Tr(indi_O,indi_O))
    coefnu(indi_O ,indi_O1)  = (1.-.064*log10(Tr(indi_O,indi_O1)))**2*sqrt(Tr(indi_O,indi_O1))
    coefnu(indi_O2,indi_O2)  =(1.-.073*log10(Tr(indi_O2,indi_O2)))**2*sqrt(Tr(indi_O2,indi_O2))
    coefnu(indi_H ,indi_O )  = (1.-.047*log10(temperature(indi_H)))**2*sqrt(temperature(indi_H))
    coefnu(indi_e ,indi_H )  = exp(-1.35e-4*temperature(indi_e))*Teh
    coefnu(indi_e ,indi_O )  = Teh+5.7e-4*Teh_3
    coefnu(indi_e ,indi_O1)  = Teh+5.7e-4*Teh_3
    coefnu(indi_e ,indi_O2)  = Teh+3.6e-2*Teh_3
    coefnu(indi_e ,indi_N2)  = exp(-1.21e-4*temperature(indi_e))*Teh

    lnCou=lnCou_o-.5*log(Nenew(iz))+1.5*log(Tenew(iz))

    do j_ion=1,nb_ion
      do indi_ion=1,nb_ion
        nuij(indi_ion,j_ion)= nuijo(indi_ion,j_ion)*lnCou*dens_i(j_ion)*Tij_r_3(indi_ion,j_ion)
      enddo
    enddo

    do indi_ion=1,nb_ion
      do j=1,nb_spn
        nuin(indi_ion,j_neutre)= nuino(indi_ion,j_neutre)*coefnu(indi_ion,j_neutre)*dens_n(j_neutre)
      end do
    end do
	  	

    Dcei=0.d0

    do i=1,nb_ion

      y1(i)=-Cij(iref,i)*G(iz)-nuij(i,indi_e)*Jpara(iz)/dens_i(indi_e)/Co(i)
      y1(i+nb_ion) = nuij(i,indi_e)*Dco(i,indi_e)*dens_i(i)*temperature(i)/dens_i(indi_e)*Jpara(iz)/Co(i)

      diff(i,i)=0.d0
      diff(i,i+nb_ion)=0.d0

      do j=1,nb_spn
        y1(i)=y1(i)+nuin(i,j)*Un(iz)/Co(i)

        diff(i,i)=diff(i,i)-nuin(i,j)
      end do

      do j=1,nb_ion
        if (j/=i) then
          diff(i,j       ) = (nuij(i,j)+nuij(i,indi_e)*dens_i(j)/dens_i(indi_e))*Cij(j,i)
          diff(i,i       ) = diff(i,i)-nuij(i,j)
          diff(i,i+nb_ion) = diff(i,i+nb_ion)+nuij(i,j)/Tij(i,j)*Ac(i,j)*dens_i_1(i)
          diff(i,j+nb_ion) = -Ac(j,i)*dens_i_1(j)*Cij(j,i)
        end if
      end do

      diff(i,i            ) = diff(i,i)-nuij(i,indi_e)*(1.-dens_i(i)/dens_i(indi_e))
      diff(i,i+nb_ion     ) = diff(i,i+nb_ion)+nuij(i,indi_e)/Tij(i,indi_e)*Ac(i,indi_e)*dens_i_1(i)
      diff(i,indi_e+nb_ion) = -Ac(indi_e,i)*dens_i_1(indi_e)*Cij(indi_e,i)

      diff(i+nb_ion,i)=0.d0
      diff(i+nb_ion,i+nb_ion)=-0.8*nuij(i,i)

      do j=1,nb_ion
        if (i/=j) then
          diff(i+nb_ion,i       ) = diff(i+nb_ion,i)+dens_i(i)*temperature(i)	&
          					*(nuij(i,j)*Dco(i,j)+nuij(i,indi_e)*Dco(i,indi_e)*(1.-dens_i(i)/dens_i(e)))
          diff(i+nb_ion,i+nb_ion)=diff(i+nb_ion,i+nb_ion)-nuij(i,j)*(Dc1(i,j)+Dco(i,j)*temperature(i)/Tij(i,j))

          diff(i+nb_ion,j       ) = -dens_i(i)*temperature(i)*Cij(j,i)		&
          				       *(nuij(i,j)*Dco(i,j)+nuij(i,indi_e)*Dco(i,indi_e)*dens_i(j)/dens_i(e))
          diff(i+nb_ion,j+nb_ion) = nuij(i,j)*dens_i(i)*dens_i_1(j)*Cij(j,i)*(Dc4(i,j)+Dco(j,i)*temperature(i)/Tij(i,j))
        end if
      end do

      do j=1,nb_spn
        diff(i+nb_ion,i+nb_ion) = diff(i+nb_ion,i+nb_ion)-nuin(i,j)*Dm1(i,j)
      end do

      diff(i+nb_ion,i+nb_ion) = diff(i+nb_ion,i+nb_ion)	&
      				-nuij(i,indi_e)*(Dc1(i,indi_e)+Dco(i,indi_e)*temperature(i)/Tij(i,indi_e))
      				
      diff(i+nb_ion,indi_e+nb_ion) = nuij(i,indi_e)*dens_i(i)*dens_i_1(indi_e)*Cij(indi_e,i)	&
				      *(Dc4(i,indi_e)+Dco(indi_e,i)*temperature(i)/Tij(i,indi_e))

      Dcei=Dcei+Dco(indi_e,i)*nuij(indi_e,i)

    end do

    Dcei=Dcei*temperature(indi_e)

    y1(indi_e+nb_ion) = -Dcei*Jpara(iz)/Co(indi_e)

    diff(indi_e+nb_ion,indi_e+nb_ion)=-0.8*nuij(indi_e,indi_e)

    do i=1,nb_ion
      diff(indi_e+nb_ion,i+nb_ion     ) = nuij(indi_e,i)*dens_i(indi_e)*dens_i_1(i)*Cij(i,indi_e)	&
                                           *(Dc4(indi_e,i)+Dco(i,indi_e)*temperature(indi_e)/Tij(indi_e,i))
      diff(indi_e+nb_ion,i            ) = (Dcei*dens_i(i)						&
                                          -nuij(indi_e,i)*Dco(indi_e,i)*dens_i(indi_e)*temperature(indi_e))*Cij(indi_e,i)

      diff(indi_e+nb_ion,indi_e+nb_ion)=diff(indi_e+nb_ion,indi_e+nb_ion)			&
					  -nuij(indi_e,i)*(Dc1(indi_e,i)+Dco(indi_e,i)*temperature(indi_e)/Tij(indi_e,i))
    end do

    do j=1,nb_spn
      diff(indi_e+nb_ion,indi_e+nb_ion)=diff(indi_e+nb_ion,indi_e+nb_ion)-nuin(indi_e,j)*Dm1(indi_e,j)
    end do

    y0  (indi_O)	  = U1old(iz)
    y1  (indi_O)	  = y1(indi_O)+(U1new(iz)-U1old(iz))/deltat_2
    y0  (indi_O+nb_ion)	  = q1old(iz)
    y1  (indi_O+nb_ion)	  = y1(indi_O+nb_ion)+(q1new(iz)-q1old(iz))/deltat_2

    y0  (indi_H)	  = U2old(iz)
    y1  (indi_H)	  = y1(indi_H)+(U2new(iz)-U2old(iz))/deltat_2
    y0  (indi_H+nb_ion)	  = q2old(iz)
    y1  (indi_H+nb_ion)	  = y1(indi_H+nb_ion)+(q2new(iz)-q2old(iz))/deltat_2

    y0  (indi_N)	  = U3old(iz)
    y1  (indi_N)	  = y1(indi_N)+(U3new(iz)-U3old(iz))/deltat_2
    y0  (indi_N+nb_ion)	  = q3old(iz)
    y1  (indi_N+nb_ion)	  = y1(indi_N+nb_ion)+(q3new(iz)-q3old(iz))/deltat_2

    y0  (indi_N2)	  = U4old(iz)
    y1  (indi_N2)	  = y1(indi_N2)+(U4new(iz)-U4old(iz))/deltat_2
    y0  (indi_N2+nb_ion)  = q4old(iz)
    y1  (indi_N2+nb_ion)  = y1(indi_N2+nb_ion)+(q4new(iz)-q4old(iz))/deltat_2

    y0  (indi_NO)	  = U5old(iz)
    y1  (indi_NO)	  = y1(indi_NO)+(U5new(iz)-U5old(iz))/deltat_2
    y0  (indi_NO+nb_ion)  = q5old(iz)
    y1  (indi_NO+nb_ion)  = y1(indi_NO+nb_ion)+(q5new(iz)-q5old(iz))/deltat_2

    y0  (indi_O2)	  = U6old(iz)
    y1  (indi_O2)	  = y1(indi_O2)+(U6new(iz)-U6old(iz))/deltat_2
    y0  (indi_O2+nb_ion)  = q6old(iz)
    y1  (indi_O2+nb_ion)  = y1(indi_O2+nb_ion)+(q6new(iz)-q6old(iz))/deltat_2

    y0  (indi_e+nb_ion)	  = qeold(iz)
    y1  (indi_e+nb_ion)	  = y1(indi_e+nb_ion)+(qenew(iz)-qeold(iz))/deltat_2

    call solve_ODE(nb_ion+indi_e,y0,y1,colli,mat,deltat_2)

    U1new(iz) = y0(indi_O)
    q1new(iz) = y0(indi_O+nb_ion)
    U2new(iz) = y0(indi_H)
    q2new(iz) = y0(indi_H+nb_ion)
    U3new(iz) = y0(indi_N)
    q3new(iz) = y0(indi_N+nb_ion)
    U4new(iz) = y0(indi_N2)
    q4new(iz) = y0(indi_N2+nb_ion)
    U5new(iz) = y0(indi_NO)
    q5new(iz) = y0(indi_NO+nb_ion)
    U6new(iz) = y0(indi_O2)
    q6new(iz) = y0(indi_O2+nb_ion)
    qenew(iz) = y0(indi_e+nb_ion)


    vmax=max(abs(U1new(iz)*Co(1)),abs(U2new(iz)*Co(2)))
    vmax=max(vmax,abs(U3new(iz)*Co(3)))
    vmax=max(vmax,abs(U4new(iz)*Co(4)))
    vmax=max(vmax,abs(U5new(iz)*Co(5)))
    vmax=max(vmax,abs(U6new(iz)*Co(6)))

  enddo

  if (isnant(N1new,nb_alt)) then
    print*,'probleme lors du calcul de N1new dans la boucle 1 apres traitement de la chimie'
    goto 246
  endif
  if (isnant(N2new,nb_alt)) then
    print*,'probleme lors du calcul de N2new dans la boucle 1 apres traitement de la chimie'
    goto 246
  endif
  if (isnant(N3new,nb_alt)) then
    print*,'probleme lors du calcul de N3new dans la boucle 1 apres traitement de la chimie'
    goto 246
  endif
  if (isnant(N4new,nb_alt)) then
    print*,'probleme lors du calcul de N4new dans la boucle 1 apres traitement de la chimie'
    goto 246
  endif
  if (isnant(N5new,nb_alt)) then
    print*,'probleme lors du calcul de N5new dans la boucle 1 apres traitement de la chimie'
    goto 246
  endif
  if (isnant(N6new,nb_alt)) then
    print*,'probleme lors du calcul de N6new dans la boucle 1 apres traitement de la chimie'
    goto 246
  endif

  if (isnant(U1new,nb_alt)) then
    print*,'probleme lors du calcul de U1new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(U2new,nb_alt)) then
    print*,'probleme lors du calcul de U2new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(U3new,nb_alt)) then
    print*,'probleme lors du calcul de U3new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(U4new,nb_alt)) then
    print*,'probleme lors du calcul de U4new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(U5new,nb_alt)) then
    print*,'probleme lors du calcul de U5new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(U6new,nb_alt)) then
    print*,'probleme lors du calcul de U6new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif

  if (isnant(q1new,nb_alt)) then
    print*,'probleme lors du calcul de q1new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(q2new,nb_alt)) then
    print*,'probleme lors du calcul de q2new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(q3new,nb_alt)) then
    print*,'probleme lors du calcul de q3new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(q4new,nb_alt)) then
    print*,'probleme lors du calcul de q4new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(q5new,nb_alt)) then
    print*,'probleme lors du calcul de q5new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(q6new,nb_alt)) then
    print*,'probleme lors du calcul de q6new dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(qenew,nb_alt)) then
    print*,'probleme lors du calcul de qenew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif

  if (isnant(T1pnew,nb_alt)) then
    print*,'probleme lors du calcul de T1pnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T1tnew,nb_alt)) then
    print*,'probleme lors du calcul de T1tnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T2pnew,nb_alt)) then
    print*,'probleme lors du calcul de T2pnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T2tnew,nb_alt)) then
    print*,'probleme lors du calcul de T2tnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T3pnew,nb_alt)) then
    print*,'probleme lors du calcul de T3pnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T3tnew,nb_alt)) then
    print*,'probleme lors du calcul de T3tnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T4pnew,nb_alt)) then
    print*,'probleme lors du calcul de T4pnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T4tnew,nb_alt)) then
    print*,'probleme lors du calcul de T4tnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T5pnew,nb_alt)) then
    print*,'probleme lors du calcul de T5pnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T5tnew,nb_alt)) then
    print*,'probleme lors du calcul de T5tnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T6pnew,nb_alt)) then
    print*,'probleme lors du calcul de T6pnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(T6tnew,nb_alt)) then
    print*,'probleme lors du calcul de T6tnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(Tepnew,nb_alt)) then
    print*,'probleme lors du calcul de Tepnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif
  if (isnant(Tetnew,nb_alt)) then
    print*,'probleme lors du calcul de Tetnew dans la boucle 1 apres traitement des collisions'
    goto 246
  endif

enddo


return
end subroutine modele
