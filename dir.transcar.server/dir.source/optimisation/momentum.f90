Real		:: lTr,str,lpr
do i=1,nbpoint

!----------------------------------!
! equation d'impulsion de l'ion H+ !
!----------------------------------!

  Tijr(indi_O,indi_H)=mirij(indi_O,indi_H)*T2new(i)+mirij(indi_H,indi_O)*T1new(i)
  Tijr(indi_O,indi_N)=mirij(indi_O,indi_N)*T3new(i)+mirij(indi_N,indi_O)*T1new(i)
  Tijr(indi_O,indi_N2)=mirij(indi_O,indi_N2)*Tmnew(i)+mirij(indi_N2,indi_O)*T1new(i)
  Tijr(indi_O,indi_NO)=mirij(indi_O,indi_NO)*Tmnew(i)+mirij(indi_NO,indi_O)*T1new(i)
  Tijr(indi_O,indi_O2)=mirij(indi_O,indi_O2)*Tmnew(i)+mirij(indi_O2,indi_O)*T1new(i)
  Tijr(indi_O,indi_H)=mirij(indi_O,indi_H)*T2new(i)+mirij(indi_H,indi_O)*T1new(i)

  Tijr(indi_H,indi_O)=Tijr(indi_O,indi_H)
  Tijr(indi_H,indi_N)=mirij(indi_H,indi_N)*T3new(i)+mirij(indi_N,indi_H)*T2new(i)
  Tijr(indi_H,indi_N2)=mirij(indi_H,indi_N2)*Tmnew(i)+mirij(indi_N2,indi_H)*T2new(i)
  Tijr(indi_H,indi_NO)=mirij(indi_H,indi_NO)*Tmnew(i)+mirij(indi_NO,indi_H)*T2new(i)
  Tijr(indi_H,indi_O2)=mirij(indi_H,indi_O2)*Tmnew(i)+mirij(indi_O2,indi_H)*T2new(i)

  Tijr(indi_N,indi_O)=Tijr(indi_O,indi_N)
  Tijr(indi_N,indi_H)=Tijr(indi_H,indi_N)
  Tijr(indi_N,indi_N2)=mirij(indi_N,indi_N2)*Tmnew(i)+mirij(indi_N2,indi_N)*T3new(i)
  Tijr(indi_N,indi_NO)=mirij(indi_N,indi_NO)*Tmnew(i)+mirij(indi_NO,indi_N)*T3new(i)
  Tijr(indi_N,indi_O2)=mirij(indi_N,indi_O2)*Tmnew(i)+mirij(indi_O2,indi_N)*T3new(i)

  Tijr(indi_N2,indi_O)=Tijr(indi_O,indi_N2)
  Tijr(indi_N2,indi_H)=Tijr(indi_H,indi_N2)
  Tijr(indi_N2,indi_N)=Tijr(indi_N,indi_N2)
  Tijr(indi_N2,indi_N2)=Tmnew(i)
  Tijr(indi_N2,indi_NO)=Tmnew(i)
  Tijr(indi_N2,indi_O2)=Tmnew(i)

  Tr(i,j)=(Ti_reel(iz,i)+Tn_reel(iz,j))/2.
  Tin(i,j)=mirin(i,j)*Tn(iz,j)+mnrin(i,j)*Temperature(iz,i)
  Tr1=T_o*T2new(i)
  Tr=(Tn(i)+Tr1)/2.

  T1_15=T1new(i)**(-1.5)
  T2_15=T2new(i)**(-1.5)
  T3_15=T3new(i)**(-1.5)
  Tm_15=Tmnew(i)**(-1.5)
  Te_15=Tenew(i)**(-1.5)
	
  thermacc(i)=thermodiff*qenew(i)/Tenew(i)**2.5

        !************************************************************************
	!*                   FREQUENCES DE COLLISION				*
	!************************************************************************

  T2r=.5*(T2new(i)*T_o+Tn(i))
  ltr=log10(T2r)
  str=sqrt(T2r)
  nuin(indi_H,indn_N2)=nuino(indi_H,indn_N2)*Nn2(i)
  nuin(indi_H,indn_O2)=nuino(indi_H,indn_O2)*No2(i)
  nuin(indi_H,indn_O)=nuino(indi_H,indn_O)*No(i)*(1.-.047*ltr)**2*str
  nuin(indi_H,indn_H)=nuino(indi_H,indn_H)*Nh(i)*(1.-.083*ltr)**2*str

  T1r=.5*(T1new(i)*T_o+Tn(i))
  nuin(indi_O,indn_N2)=nuino(indi_O,indn_N2)*Nn2(i)
  nuin(indi_O,indn_O2)=nuino(indi_O,indn_O2)*No2(i)
  nuin(indi_O,indn_O)=nuino(indi_O,indn_O)*No(i)*(1.-.064*alog10(T1r))**2*sqrt(T1r)
  nuin(indi_O,indn_H)=nuino(indi_O,indn_H)*Nh(i)


  nuin(indi_N,indn_N2)=nuino(indi_N,indn_N2)*Nn2(i)
  nuin(indi_N,indn_O2)=nuino(indi_N,indn_O2)*No2(i)
  nuin(indi_N,indn_O)=nuino(indi_N,indn_O)*No(i)
  nuin(indi_N,indn_H)=nuino(indi_N,indn_H)*Nh(i)

  Tmr=.5*(Tmnew(i)*T_o+Tn(i))
  ltr=log10(Tmr)
  str=sqrt(Tmr)
  nuin(indi_N2,indn_N2)=nuino(indi_N2,indn_N2)*Nn2(i)*(1.-.069*ltr)**2*str
  nuin(indi_N2,indn_O2)=nuino(indi_N2,indn_O2)*No2(i)
  nuin(indi_N2,indn_O)=nuino(indi_N2,indn_O)*No(i)
  nuin(indi_N2,indn_H)=nuino(indi_N2,indn_H)*Nh(i)

  nuin(indi_NO,indn_N2)=nuino(indi_NO,indn_N2)*Nn2(i)
  nuin(indi_NO,indn_O2)=nuino(indi_NO,indn_O2)*No2(i)
  nuin(indi_NO,indn_O)=nuino(indi_NO,indn_O)*No(i)
  nuin(indi_NO,indn_H)=nuino(indi_NO,indn_H)*Nh(i)

  nuin(indi_O2,indn_N2)=nuino(indi_O2,indn_N2)*Nn2(i)
  nuin(indi_O2,indn_O2)=nuino(indi_O2,indn_O2)*No2(i)*(1.-.073*ltr)**2*str
  nuin(indi_O2,indn_O)=nuino(indi_O2,indn_O)*No(i)
  nuin(indi_O2,indn_H)=nuino(indi_O2,indn_H)*Nh(i)



	Tij=Tij_reel/T_o
   	Teh=sqrt(Ti_reel(iz,e))
        Vi=sqrt(Vz(iz,:)**2+Vy(iz,:)**2)
        Tij_reel_1=1./Tij_reel
	Tij_reel_3=Tij_reel_1*sqrt(Tij_reel_1)




  nuij(indi_H,indi_H)=nuijo(indi_H,indi_O)*N1new(i)*Tij_reel_3(i,j)


            lpr=log(Nenew(i)*Tepnew(i))
	    C2am1(i)=-Tepnew(i)
	    D2am1(i)=lpr
            C2bm1(i)=-T1pnew(i)*N1new(i)/xn2(i)
	    D2bm1(i)=log(xn1(i)*T1pnew(i))
	
	    C2am2(i)=-Tepnew(i)
	    D2am2(i)=lpr
            C2bm2(i)=-T2pnew(i)*N2new(i)/xn2(i)
	    D2bm2(i)=log(xn2(i)*T2pnew(i))
	
	    C2am3(i)=-Tepnew(i)
	    D2am3(i)=lpr
            C2bm3(i)=-T3pnew(i)*N3new(i)/xn3(i)
	    D2bm3(i)=log(xn3(i)*T3pnew(i))
	
	    C2amm(i)=-Tepnew(i)
	    D2amm(i)=lpr
            C2bmm(i)=-Tmpnew(i)*Nmnew(i)/xnm(i)
	    D2bmm(i)=log(xnm(i)*Tmpnew(i))
	
enddo

call velocity(Vel1m,Ipos1,Iposnp,deltat_2)
	
D2al=(D2a(1)+D2a(2))/2.
D2bl=(D2b(1)+D2b(2))/2.
	
D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
D2br=.5*log(xn1(np)*T1pnew(np)*xn1(nx)*T1pnew(nx))

call sources(Ipos1,Iposn,deltat_2,2,C2am1,D2am1,D2al,D2ar)
call sources(Ipos1,Iposn,deltat_2,2,C2bm1,D2bm1,D2bl,D2br)
call sources(Ipos1,Iposn,deltat_2,8,zero,D3m1,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7m1,0.,0.)
lbc=1.
call lcpfct(U1old,U1new,Ipos1,Iposn,lbc,0.,0.,U1new(np),.false.,0)

if (isnant(U1new,nx)) then
  print*,'probleme lors du calcul de U1new dans la boucle ',iboucle
  goto 246
endif

call velocity(Vel2m,Ipos1,Iposnp,deltat_2)
	
D2bl=(D2b(1)+D2b(2))/2.
	
D2br=.5*log(xn2(np)*T2pnew(np)*xn2(nx)*T2pnew(nx))

call sources(Ipos1,Iposn,deltat_2,2,C2am2,D2am2,D2al,D2ar)
call sources(Ipos1,Iposn,deltat_2,2,C2bm2,D2bm2,D2bl,D2br)
call sources(Ipos1,Iposn,deltat_2,8,zero,D3m1,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7m1,0.,0.)
lbc=1.
call lcpfct(U2old,U2new,Ipos1,Iposn,lbc,0.,0.,U2new(np),.false.,0)

if (isnant(U2new,nx)) then
  print*,'probleme lors du calcul de U2new dans la boucle ',iboucle
  goto 246
endif

call velocity(Vel3m,Ipos1,Iposnp,deltat_2)
	
D2bl=(D2b(1)+D2b(2))/2.
	
D2br=.5*log(xn3(np)*T3pnew(np)*xn3(nx)*T3pnew(nx))

call sources(Ipos1,Iposn,deltat_2,2,C2am3,D2am3,D2al,D2ar)
call sources(Ipos1,Iposn,deltat_2,2,C2bm3,D2bm3,D2bl,D2br)
call sources(Ipos1,Iposn,deltat_2,8,zero,D3m3,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7m3,0.,0.)
lbc=1.
call lcpfct(U3old,U3new,Ipos1,Iposn,lbc,0.,0.,U3new(np),.false.,0)

if (isnant(U3new,nx)) then
  print*,'probleme lors du calcul de U3new dans la boucle ',iboucle
  goto 246
endif
	
call velocity(Velmm,Ipos1,Iposnp,deltat_2)
	
D2bl=(D2b(1)+D2b(2))/2.
	
D2br=.5*log(xnm(np)*Tmpnew(np)*xnm(nx)*Tmpnew(nx))

call sources(Ipos1,Iposn,deltat_2,2,C2amm,D2amm,D2al,D2ar)
call sources(Ipos1,Iposn,deltat_2,2,C2bmm,D2bmm,D2bl,D2br)
call sources(Ipos1,Iposn,deltat_2,8,zero,D3mm,0.,0.)
call sources(Ipos1,Iposn,deltat_2,7,zero,D7mm,0.,0.)
lbc=1.
call lcpfct(Umold,Umnew,Ipos1,Iposn,lbc,0.,0.,Umnew(np),.false.,0)

if (isnant(Umnew,nx)) then
  print*,'probleme lors du calcul de Umnew dans la boucle ',iboucle
  goto 246
endif

do i=1,npoint

  Uenew(i)=(Co(1)*N1new(i)*U1new(i)+Co(2)*N2new(i)*U2new(i)	&
          +Co(3)*N3new(i)*U3new(i)+Co(4)*Nmnew(i)*Umnew(i)	&
          -JJ(i)/N_o-0.*Jes(i)/N_o)/Nenew(i)/Co(0)
enddo

	
	
	
	    D3(i)=-G(i)+thermacc(i)/mi/Ci0
     &            -3.*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*(T2pnew(i)-T2tnew(i))*alt_geo_1(i)
     &		  +nuij(i)*Cji*U1new(i)
     &		  +(nuik(i)*Czi+nuil(i)*Czi+nuim(i)*Czi)*Umnew(i)
     &		  +nuin(i)*Cni*U3new(i)

	    D3(i)=D3(i)+.6*q2new(i)/xn2(i)
     &		       *(28.*nuiN2(i)/(Tn(i)/T_0+28.*T2new(i))
     &                  +32.*nuiO2(i)/(Tn(i)/T_0+32.*T2new(i))
     &                  +16.*nuiO(i)/(Tn(i)/T_0+16.*T2new(i))
     &                  +    nuiH(i)/(Tn(i)/T_0+T2new(i))
     &                  +16.*nuij(i)/(T1new(i)+16.*T2new(i))
     &                  +28.*nuik(i)/(Tmnew(i)+28.*T2new(i))
     &                  +32.*nuil(i)/(Tmnew(i)+32.*T2new(i))
     &                  +30.*nuim(i)/(Tmnew(i)+30.*T2new(i))
     &                  +14.*nuin(i)/(T3new(i)+14.*T2new(i))
     &                  +5.46e-4*nuie(i)/(Tenew(i)+5.46e-4*T2new(i)))

	    D3(i)=D3(i)-.6*q1new(i)/xn1(i)*Cji
     &			  *nuij(i)/(T1new(i)+16.*T2new(i))

	    D3(i)=D3(i)-.6*qenew(i)/xne(i)*Cei
     &			  *nuie(i)/(Tenew(i)+5.46e-4*T2new(i))

	    D3(i)=D3(i)-.6*q3new(i)/xn3(i)*Cni
     &			  *nuin(i)/(T3new(i)+14.*T2new(i))
	
	    D3(i)=D3(i)+(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiH(i))
     &                    *(Un(i)/Ci0)
     &                 -.6*nuiN2(i)/(Tn(i)/T_0+28.*T2new(i))
     &                    *(q_Nn2(i)*Niqi0)
     &                 -.6*nuiO2(i)/(Tn(i)/T_0+32.*T2new(i))
     &                    *(q_No2(i)*Niqi0)
     &                 -.6*nuiO(i)/(Tn(i)/T_0+16.*T2new(i))
     &                    *(q_No(i)*Niqi0)
     &                 -.6*nuiH(i)/(Tn(i)/T_0+T2new(i))
     &                    *(q_Nh(i)*Niqi0)
     &                 +nuie(i)*Cei*Uenew(i)

	    D7(i)=-(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiH(i)
     &              +nuij(i)+nuik(i)+nuil(i)+nuim(i)
     &              +nuin(i)+nuie(i))

	  enddo



C]]]

C[[[    O+ momentum equation resolution

	  do i=1,nx



	    D3(i)=-Cij*G(i)+thermacc(i)/mj/Cj0
     &            -3.*Cji*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*Cji*(T1pnew(i)-T1tnew(i))*alt_geo_1(i)
     &		  +nuji(i)*Cij*U2new(i)

	    D3(i)=D3(i)+(nujk(i)+nujl(i)+nujm(i))*Czj*Umnew(i)
     &		       +nujn(i)*Cnj*U3new(i)

	    D3(i)=D3(i)-.6*q2new(i)/xn2(i)*Cij*16.
     &			*nuji(i)/(16.*T2new(i)+T1new(i))

	    D3(i)=D3(i)+.6*q1new(i)/xn1(i)
     &			*(28.*nujN2(i)/(16.*Tn(i)/T_0+28.*T1new(i))
     &                   +32.*nujO2(i)/(16.*Tn(i)/T_0+32.*T1new(i))
     &                   +16.*nujO(i)/(16.*Tn(i)/T_0+16.*T1new(i))
     &                   +    nuji(i)/(16.*T2new(i)+T1new(i))
     &                   +28.*nujk(i)/(16.*Tmnew(i)+28.*T1new(i))
     &                   +32.*nujl(i)/(16.*Tmnew(i)+32.*T1new(i))
     &                   +30.*nujm(i)/(16.*Tmnew(i)+30.*T1new(i))
     &                   +14.*nujn(i)/(16.*T3new(i)+14.*T1new(i))
     &                   +    nuje(i)/(16.*Tenew(i)+5.46e-4*T1new(i))
     &                        *5.46e-4)

	    D3(i)=D3(i)-.6*qenew(i)/xne(i)*Cej*16.
     &			*nuje(i)/(16.*Tenew(i)+5.46e-4*T1new(i))

	    D3(i)=D3(i)-.6*q3new(i)/xn3(i)*Cnj*16.
     &			*nujn(i)/(16.*T3new(i)+14.*T1new(i))

	    D3(i)=D3(i)+(nujN2(i)+nujO2(i)+nujO(i))*(Un(i)/Cj0)
     &                 -.6*nujN2(i)/(16.*Tn(i)/T_0+28.*T1new(i))
     &                  *(16.*q_Nn2(i)*Njqj0)
     &                 -.6*nujO2(i)/(16.*Tn(i)/T_0+32.*T1new(i))
     &                  *(16.*q_No2(i)*Njqj0)
     &                 -.6*nujO(i)/(16.*Tn(i)/T_0+16.*T1new(i))
     &                  *(16.*q_No(i)*Njqj0)
     &                 +nuje(i)*Cje*Uenew(i)

	    D7(i)=-(nujN2(i)+nujO2(i)+nujO(i)
     &              +nuji(i)+nujn(i)
     &              +nuje(i)+nujk(i)+nujl(i)+nujm(i))

	  enddo




C]]]

C[[[    heavy ions momentum equation resolution

	  do i=1,nx



          D3(i)=-Ciz*G(i)+thermacc(i)/mz/Cz0
     &            -3.*Czi*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*Czi*(Tmpnew(i)-Tmtnew(i))*alt_geo_1(i)
     &		 +(nukj(i)+nulj(i)+numj(i))*Cjz/3.*U1new(i)
     &		 +(nuki(i)+nuli(i)+numi(i))*Ciz/3.*U2new(i)
     &		 +(nukn(i)+nuln(i)+numn(i))*Cnz/3.*U3new(i)

	  D3(i)=D3(i)-.2*q2new(i)/xn2(i)*Ciz
     &		      *(28.*nuki(i)/(28.*T2new(i)+Tmnew(i))
     &                 +32.*nuli(i)/(32.*T2new(i)+Tmnew(i))
     &                 +30.*numi(i)/(30.*T2new(i)+Tmnew(i)))

	  D3(i)=D3(i)-.2*q1new(i)/xn1(i)*Cjz
     &		      *(28.*nukj(i)/(28.*T1new(i)+16.*Tmnew(i))
     &                 +32.*nulj(i)/(32.*T1new(i)+16.*Tmnew(i))
     &                 +30.*numj(i)/(30.*T1new(i)+16.*Tmnew(i)))

	  D3(i)=D3(i)-.2*qenew(i)/xne(i)*Cez
     &		      *(28.*nuke(i)/(28.*Tenew(i)+5.46e-4*Tmnew(i))
     &                 +32.*nule(i)/(32.*Tenew(i)+5.46e-4*Tmnew(i))
     &                 +30.*nume(i)/(30.*Tenew(i)+5.46e-4*Tmnew(i)))

	  D3(i)=D3(i)-.2*q3new(i)/xn3(i)*Cnz
     &		      *(28.*nukn(i)/(28.*T3new(i)+14.*Tmnew(i))
     &                 +32.*nuln(i)/(32.*T3new(i)+14.*Tmnew(i))
     &                 +30.*numn(i)/(30.*T3new(i)+14.*Tmnew(i)))

	  D3(i)=D3(i)+(nukN2(i)+nukO2(i)+nukO(i)
     &               +nulN2(i)+nulO2(i)+nulO(i)
     &               +numN2(i)+numO2(i)+numO(i))*Un(i)/Cz0/3.
     &           -.2*
     &               (nukN2(i)/(28.*Tn(i)/T_0+28.*Tmnew(i))
     &           *(28.*q_Nn2(i)*Nkqk0)
     &               +nukO2(i)/(28.*Tn(i)/T_0+32.*Tmnew(i))
     &           *(28.*q_No2(i)*Nkqk0)
     &               +nukO(i)/(28.*Tn(i)/T_0+16.*Tmnew(i))
     &           *(28.*q_No(i)*Nkqk0))
     &           -.2*
     &               (nulN2(i)/(32.*Tn(i)/T_0+28.*Tmnew(i))
     &           *(32.*q_Nn2(i)*Nlql0)
     &               +nulO2(i)/(32.*Tn(i)/T_0+32.*Tmnew(i))
     &           *(32.*q_No2(i)*Nlql0)
     &               +nulO(i)/(32.*Tn(i)/T_0+16.*Tmnew(i))
     &           *(32.*q_No(i)*Nlql0))
     &           -.2*
     &               (numN2(i)/(30.*Tn(i)/T_0+28.*Tmnew(i))
     &           *(30.*q_Nn2(i)*Nmqm0)
     &                +numO2(i)/(30.*Tn(i)/T_0+32.*Tmnew(i))
     &           *(30.*q_No2(i)*Nmqm0)
     &                +numO(i)/(30.*Tn(i)/T_0+16.*Tmnew(i))
     &           *(30.*q_No(i)*Nmqm0))
     &           +Uenew(i)*Cez*(nuke(i)+nule(i)+nume(i))/3.

          D7(i)=-(nukN2(i)+nukO2(i)+nukO(i)
     &            +nuki(i)+nuke(i)+nukj(i)+nukn(i)
     &            +nulN2(i)+nulO2(i)+nulO(i)
     &            +nuli(i)+nule(i)+nulj(i)+nuln(i)
     &            +numN2(i)+numO2(i)+numO(i)
     &            +numi(i)+nume(i)+numj(i)+numn(i))/3.
	

	  enddo

	
	  do i=1,nx




	    D3(i)=-Cin*G(i)+thermacc(i)/mn/Cn0
     &            -3.*Cni*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*Cni*(T3pnew(i)-T3tnew(i))*alt_geo_1(i)
     &		  +nuni(i)*Cin*U2new(i)+nunj(i)*Cjn*U1new(i)
     &		  +(nunk(i)+nunl(i)+nunm(i))*Czn*Umnew(i)

	  D3(i)=D3(i)-.6*q2new(i)/xn2(i)*Cin*14.
     &			*nuni(i)/(14.*T2new(i)+T3new(i))

	  D3(i)=D3(i)-.6*q1new(i)/xn1(i)*Cjn*14.
     &			*nunj(i)/(14.*T1new(i)+16.*T3new(i))

	  D3(i)=D3(i)-.6*qenew(i)/xne(i)*Cen*14.
     &			*nune(i)/(14.*Tenew(i)+5.46e-4*T3new(i))

	  D3(i)=D3(i)+.6*q3new(i)/xn3(i)
     &			*(28.*nunN2(i)/(14.*Tn(i)/T_0+28.*T3new(i))
     &                  +32.*nunO2(i)/(14.*Tn(i)/T_0+32.*T3new(i))
     &                  +16.*nujO(i)/(14.*Tn(i)/T_0+16.*T3new(i))
     &                  +    nuni(i)/(14.*T2new(i)+T3new(i))
     &                  +28.*nunk(i)/(14.*Tmnew(i)+28.*T3new(i))
     &                  +32.*nunl(i)/(14.*Tmnew(i)+32.*T3new(i))
     &                  +30.*nunm(i)/(14.*Tmnew(i)+30.*T3new(i))
     &                  +    nune(i)/(14.*Tenew(i)+5.46e-4*T3new(i))
     &                       *5.46e-4)

	  D3(i)=D3(i)+(nunN2(i)+nunO2(i)+nunO(i))*Un(i)/Cn0
     &           -.6*nunN2(i)/(14.*Tn(i)/T_0+28.*T3new(i))
     &           *(14.*q_Nn2(i)*Nnqn0)
     &           -.6*nunO2(i)/(14.*Tn(i)/T_0+32.*T3new(i))
     &           *(14.*q_No2(i)*Nnqn0)
     &           -.6*nunO(i)/(14.*Tn(i)/T_0+16.*T3new(i))
     &           *(14.*q_No(i)*Nnqn0)
     &           +nune(i)*Cen*Uenew(i)

	    D7(i)=-(nunN2(i)+nunO2(i)+nunO(i)
     &             +nuni(i)+nunj(i)+nune(i)
     &             +nunk(i)+nunl(i)+nunm(i))



	  enddo


C[[[    Velocities corrections

