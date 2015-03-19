	real chim(nb_ion,nb_ion),y0(nb_ion),y1(nb_ion)

	r4_9=4./9.
	r8_9=8./9.
	r11_9=11./9.
	r17_9=17./9.

!-------------------------------------!
!     equations du flux de chaleur    !
!-------------------------------------!

!	ion O+
	
    call velocity(Vel1,Ipos1,Iposnp,deltat_2)

    do iz=1,nb_alt

      C2a(iz)=-2.2*Cji*q1new(iz)
      D2a(iz)=U1new(iz)
      C2b(iz)=-Cji*(r11_18*T1pnew(iz)+r8_9*T1tnew(iz))*N1new(iz)
      D2b(iz)=T1pnew(iz)
      C2c(iz)=-Cji*(r17_9*T1pnew(iz)-r8_9*T1tnew(iz))*N1new(iz)
      D2c(iz)=T1tnew(iz)
      C2d(iz)=r4_9*Cji*(T1pnew(iz)-T1tnew(iz))**2
      D2d(iz)=N1new(iz)
    enddo

    D2al=(D2a(1)+D2a(2))/2.
    D2ar=ylimd(Radn(np),D2a,extra,nx)
    D2bl=(D2b(1)+D2b(2))/2.
    D2br=ylimd(Radn(np),D2b,extra,nx)
    D2cl=(D2c(1)+D2c(2))/2.
    D2cr=ylimd(Radn(np),D2c,extra,nx)
    D2dl=(D2d(1)+D2d(2))/2.
    D2dr=ylimd(Radn(np),D2d,extra,nx)
	
    D2ar=.5*(U1new(np)+U1new(nx))
    D2br=.5*(T1pnew(np)+T1pnew(nx))
    D2cr=.5*(T1tnew(np)+T1tnew(nx))
    D2dr=sqrt(N1new(np)*N1new(nx))
	
    call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
    call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
    call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
    call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
    lbc=1.
    rbc=bclimd(Radn(nx+1),q1old,extra,nx)
    rbc=1.
    call lcpfct(q1old,q1new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

    q1new(nx)=max(0.,q1new(nx))

          if (isnant(q1new,nx)) then
            print*,'probleme lors du calcul de q1new dans la boucle 1'
            goto 246
          endif
		

!	ion H+
	
	  call velocity(Veliq,Ipos1,Iposnp,deltat_2)

	 do iz=1,nb_alt

	    C2a(iz)=-2.2*Cji*q2new(iz)
	    D2a(iz)=U1new(iz)
	    C2b(iz)=-(r11_18*T2pnew(iz)+r8_9*T2tnew(iz))*N2new(iz)
	    D2b(iz)=T2pnew(iz)
	    C2c(iz)=-(r17_9*T2pnew(iz)-r8_9*T2tnew(iz))*N2new(iz)
	    D2c(iz)=T2tnew(iz)
	    C2d(iz)=r4_9*(T2pnew(iz)-T2tnew(iz))**2
	    D2d(iz)=N2new(iz)
	enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
          D2cl=(D2c(1)+D2c(2))/2.
          D2cr=ylimd(Radn(np),D2c,extra,nx)
          D2dl=(D2d(1)+D2d(2))/2.
          D2dr=ylimd(Radn(np),D2d,extra,nx)
	
	  D2ar=.5*(U2new(np)+U2new(nx))
	  D2br=.5*(T2pnew(np)+T2pnew(nx))
	  D2cr=.5*(T2tnew(np)+T2tnew(nx))
	  D2dr=sqrt(N2new(np)*N2new(nx))
	
	  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
	  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
	  call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
	  call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
          lbc=1.
          rbc=bclimd(Radn(nx+1),q2old,extra,nx)
	  rbc=1.
	  call lcpfct(q2old,q2new,Ipos1,Iposn,
     &		      lbc,0.,rbc,0.,.false.,0)

          q2new(nx)=max(0.,q2new(nx))

          if (isnant(q2new,nx)) then
            print*,'probleme lors du calcul de q2new dans la boucle 1'
            goto 246
          endif

          	
!	ion N+
	
	  call velocity(Velnq,Ipos1,Iposnp,deltat_2)

	 do iz=1,nb_alt

	    C2a(iz)=-2.2*Cji*q3new(iz)
	    D2a(iz)=U3new(iz)
	    C2b(iz)=-Cni*(r11_18*T3pnew(iz)+r8_9*T3tnew(iz))*N3new(iz)
	    D2b(iz)=T3pnew(iz)
	    C2c(iz)=-Cni*(r17_9*T3pnew(iz)-r8_9*T3tnew(iz))*N3new(iz)
	    D2c(iz)=T3tnew(iz)
	    C2d(iz)=r4_9*Cni*(T3pnew(iz)-T3tnew(iz))**2
	    D2d(iz)=N3new(iz)
	enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
          D2cl=(D2c(1)+D2c(2))/2.
          D2cr=ylimd(Radn(np),D2c,extra,nx)
          D2dl=(D2d(1)+D2d(2))/2.
          D2dr=ylimd(Radn(np),D2d,extra,nx)
	
	  D2ar=.5*(U3new(np)+U3new(nx))
	  D2br=.5*(T3pnew(np)+T3pnew(nx))
	  D2cr=.5*(T3tnew(np)+T3tnew(nx))
	  D2dr=sqrt(N3new(np)*N3new(nx))
	
	  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
	  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
	  call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
	  call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
          lbc=1.
          rbc=bclimd(Radn(nx+1),q3old,extra,nx)
	  rbc=1.
	  call lcpfct(q3old,q3new,Ipos1,Iposn,
     &		      lbc,0.,rbc,0.,.false.,0)

          q3new(nx)=max(0.,q3new(nx))

          if (isnant(q3new,nx)) then
            print*,'probleme lors du calcul de q3new dans la boucle 1'
            goto 246
          endif
		


!-------------------------------------!
!        equations d'impulsion        !
!-------------------------------------!


!       calcul des collisions
	
	  do iz=1,nb_alt

          Tr=(T_0*T1new(iz)+Tn(iz))/2.
	  T1_15(iz)=T1new(iz)**(-1.5)/nu_0
	  T2_15(iz)=T2new(iz)**(-1.5)/nu_0
	  T3_15(iz)=T3new(iz)**(-1.5)/nu_0
	  Tm_15(iz)=Tmnew(iz)**(-1.5)/nu_0
	  Te_15(iz)=Tenew(iz)**(-1.5)/nu_0
          Tr1=T_0*T2new(iz)
          Tr=(Tn(iz)+Tr1)/2.
	
          Tr=(T_0*Tmnew(iz)+Tn(iz))/2.

          nujN2(iz)=6.82e-10*Nn2(iz)                              *t0
          nujO2(iz)=6.64e-10*No2(iz)                              *t0
          nujO (iz)=3.67e-11*No(iz)*sqrt(Tr)
     &             *(1.-.064*alog10(Tr))**2                     *t0

          nujj(iz)=.22*N1new(iz)
     &           /(T1new(iz))**1.5				/nu_0
          nuji(iz)=0.077*N2new(iz)
     &           /((T1new(iz)+16.*T2new(iz))/17.)**1.5		/nu_0
          nujn(iz)=0.22*N3new(iz)
     &           /((14.*T1new(iz)+16.*T3new(iz))/30.)**1.5 	/nu_0
          nujk(iz)=0.25*N4new(iz)
     &           /((28.*T1new(iz)+16.*Tmnew(iz))/44.)**1.5	/nu_0
          nujm(iz)=0.26*N5new(iz)
     &           /((30.*T1new(iz)+16.*Tmnew(iz))/46.)**1.5	/nu_0
          nujl(iz)=0.26*N6new(iz)
     &	         /((32.*T1new(iz)+16.*Tmnew(iz))/48.)**1.5	/nu_0
          nuje(iz)=1.87e-3*Nenew(iz)*Te_15(iz)


	  thermacc(iz)=thermodiff*qenew(iz)/Tenew(iz)**2.5

          nuiN2(iz)=3.36e-9*Nn2(iz)                               *t0

          nuiO2(iz)=3.20e-9*No2(iz)                               *t0
          nuiO (iz)=6.61e-11*No (iz)
     &             *(1.-.047*alog10(Tr1))**2*sqrt(Tr1)          *t0
          nuiH (iz)=2.65e-10*Nh (iz)
     &             *(1.-.083*alog10(Tr))**2*sqrt(Tr)            *t0
          nuij (iz) = 1.23*N1new(iz)
     &              /((16.*T2new(iz)+T1new(iz))/17.)**1.5		/nu_0
          nuii (iz) = .9*N2new(iz)
     &              /(T2new(iz))**1.5 				/nu_0
          nuin (iz) = 1.23*N3new(iz)
     &              /((14.*T2new(iz)+T3new(iz))/15.)**1.5 	/nu_0
          nuik (iz) = 1.25*N4new(iz)
     &              /((28.*T2new(iz)+Tmnew(iz))/27.)**1.5		/nu_0
          nuim (iz) = 1.25*N5new(iz)
     &              /((30.*T2new(iz)+Tmnew(iz))/31.)**1.5  	/nu_0
          nuil (iz) = 1.25*N6new(iz)
     &              /((32.*T2new(iz)+Tmnew(iz))/33.)**1.5		/nu_0
          nuie (iz)=0.03    *Nenew(iz)*Te_15(iz)



          nunH (iz)=0.145e-9*Nh(iz)                               *t0
          nunN2(iz)=0.747e-9*Nn2(iz)                              *t0
          nunO (iz)=0.442e-9*No (iz)                              *t0
          nunO2(iz)=0.725e-9*No2 (iz)                             *t0
          nunj(iz)=0.25*N1new(iz)
     &           / ((16.*T3new(iz)+14.*T1new(iz))/30.)**1.5	/nu_0
          nuni(iz)=0.088*N2new(iz)
     &            / (( T3new(iz)+14.*T2new(iz))/15.)**1.5		/nu_0
          nunn(iz)=0.24*N3new(iz)/ (T3new(iz))**1.5		/nu_0
          nunk(iz)=0.28*N4new(iz)
     &            / ((28.*T3new(iz)+14.*Tmnew(iz))/42.)**1.5	/nu_0
          nunm(iz)=0.28*N5new(iz)
     &            / ((30.*T3new(iz)+14.*Tmnew(iz))/44.)**1.5	/nu_0
          nunl(iz)=0.28*N6new(iz)
     &            / ((32.*T3new(iz)+14.*Tmnew(iz))/46.)**1.5	/nu_0
          nune(iz)=2.136e-3*Nenew(iz)*Te_15(iz)




          nukH(iz)=0.074e-9*Nh(iz)                                *t0
          nukO2(iz)=0.449e-9*No2(iz)                              *t0
          nukO (iz)=0.258e-9*No (iz)                              *t0
          nukN2 (iz)=5.14e-11*Nn2 (iz)
     &             *(1.-.069*alog10(Tr))**2*sqrt(Tr)            *t0

          nukj (iz)=0.15    *N1new(iz)*T1_15(iz)
          nuki (iz)=0.045   *N2new(iz)*T2_15(iz)
          nukn (iz)=0.14    *N3new(iz)*T3_15(iz)
          nukk (iz)=0.17    *N4new(iz)*Tm_15(iz)
          nukm (iz)=0.17    *N5new(iz)*Tm_15(iz)
          nukl (iz)=0.18    *N6new(iz)*Tm_15(iz)
          nuke (iz)=1.06e-3 *Nenew(iz)*Tm_15(iz)

          nulH(iz)=0.065e-9*Nh(iz)                                *t0
          nulN2(iz)=0.413e-9*Nn2(iz)                              *t0
          nulO (iz)=0.231e-9*No (iz)                              *t0
          nulO2 (iz)=2.59e-11*No2 (iz)
     &             *(1.-.073*alog10(Tr))**2*sqrt(Tr)            *t0
          nulj (iz)=0.13    *N1new(iz)*T1_15(iz)
          nuli (iz)=0.039   *N2new(iz)*T2_15(iz)
          nuln (iz)=0.12    *N3new(iz)*T3_15(iz)
          nulk (iz)=0.15    *N4new(iz)*Tm_15(iz)
          nulm (iz)=0.16    *N5new(iz)*Tm_15(iz)
          null (iz)=0.16    *N6new(iz)*Tm_15(iz)
          nule (iz)=9.276e-4*Nenew(iz)*Te_15(iz)

          numH(iz)=0.069e-9*Nh(iz)                                *t0
          numN2(iz)=0.434e-9*Nn2(iz)                              *t0
          numO (iz)=0.244e-9*No (iz)                              *t0
          numO2(iz)=0.427e-9*No2 (iz)                             *t0
          numj (iz)=0.14    *N1new(iz)*T1_15(iz)
          numi (iz)=0.042   *N2new(iz)*T2_15(iz)
          numn (iz)=0.13    *N3new(iz)*T3_15(iz)
          numk (iz)=0.16    *N4new(iz)*Tm_15(iz)
          numm (iz)=0.16    *N5new(iz)*Tm_15(iz)
          numl (iz)=0.17    *N6new(iz)*Tm_15(iz)
          nume (iz)=9.894e-4*Nenew(iz)*Tm_15(iz)

	enddo

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
    D2ar=ylimd(Radn(np),D2a,extra,nx)
    D2bl=(D2b(1)+D2b(2))/2.
    D2br=ylimd(Radn(np),D2b,extra,nx)
	
    call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
    call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
    lbc=1.
    rbc=1.
    call lcpfct(U1old,U1new,Ipos1,Iposn,lbc,0.,rbc,0.,.false.,0)

!	ion H+

	  call velocity(Velim,Ipos1,Iposnp,deltat_2)
	
	  do iz=1,nb_alt

	    C2a(iz)=-Tepnew(iz)
!	    D2a(iz)=log(Nenew(iz)*Tepnew(iz))
            C2b(iz)=-T2pnew(iz)*N2new(iz)/xn2(iz)
	    D2b(iz)=log(xn2(iz)*T2pnew(iz))

	  enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
	
!	  D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
	  D2br=.5*log(xn2(np)*T2pnew(np)*xn2(nx)*T2pnew(nx))
	  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
	  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
          lbc=1.
          rbc=1.
	  call lcpfct(U2old,U2new,Ipos1,Iposn,
     &		      lbc,0.,0.,U2new(np),.false.,0)

!	ion N+

	  call velocity(Velnm,Ipos1,Iposnp,deltat_2)
	
	  do iz=1,nb_alt

	    C2a(iz)=-Cni*Tepnew(iz)
!	    D2a(iz)=log(Nenew(iz)*Tepnew(iz))
            C2b(iz)=-Cni*T3pnew(iz)*N3new(iz)/xn3(iz)
            D2b(iz)=log(xn3(iz)*T3pnew(iz))

	  enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
	
!	  D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
	  D2br=.5*log(xn3(np)*T3pnew(np)*xn3(nx)*T3pnew(nx))
	
	  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
	  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U3old,extra,nx)
	  rbc=1.
	  call lcpfct(U3old,U3new,Ipos1,Iposn,
     &		      lbc,0.,0.,U3new(np),.false.,0)


!	ions lourds

	  call velocity(Velmm,Ipos1,Iposnp,deltat_2)
	
	  do iz=1,nb_alt

            C2a(iz)=-Czi*Tepnew(iz)
!  	    D2a(iz)=log(Nenew(iz)*Tepnew(iz))
            C2b(iz)=-Czi*Tmpnew(iz)*Nmnew(iz)/xnm(iz)
            D2b(iz)= log(xnm(iz)*Tmpnew(iz))
	  enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
	
!	  D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
	  D2br=.5*log(xnm(np)*Tmpnew(np)*xnm(nx)*Tmpnew(nx))
	
	  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
	  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U1old,extra,nx)
	  rbc=1.
	  call lcpfct(Umold,Umnew,Ipos1,Iposn,
     &		      lbc,0.,0.,Umnew(np),.false.,0)



















	  D3(iz)=Cji*alt_geo_1(iz)*N1new(iz)*(T1pnew(iz)-T1tnew(iz))
     &		                   *(T1pnew(iz)-4.*T1tnew(iz))
     &	         +N1new(iz)*T1new(iz)*U2new(iz)*Cij*
     &                   nuji(iz)*(2.412-2.5)

          D3(iz)=D3(iz)-N1new(iz)*T1new(iz)*U1new(iz)*          (
     &                         nujN2(iz)*1.545+nujO2(iz)*1.5
     &                        +nujO(iz)*1.75  +nuji(iz)*2.412
     &                        +nujn(iz)*1.8   +nujk(iz)*1.545
     &                        +nujl(iz)*1.5   +nujm(iz)*1.522
     &                        +nuje(iz)*2.5                 )
     &               + 2.5*N1new(iz)*T1new(iz)*U1new(iz)*(nujN2(iz)
     &                              +nujO2(iz)+nujO(iz)
     &                              +nuji(iz)+nujn(iz)
     &                              +nujk(iz)+nujl(iz)+nujm(iz)+nuje(iz))

          D3(iz)=D3(iz)+N1new(iz)*T1new(iz)*Umnew(iz)*Czj*                (
     &                         nujk(iz)*1.545+nujl(iz)*1.5+nujm(iz)*1.522
     &                        - 2.5*(nujk(iz)+nujl(iz)+nujm(iz))        )

          D3(iz)=D3(iz)-N1new(iz)*T1new(iz)*U3new(iz)*Cnj*.7*nujn(iz)

          D3(iz)=D3(iz)-N1new(iz)*q2new(iz)*Cij*xn2_1(iz)*nuji(iz)*
     &                   (1.262- 24.*T1new(iz)/(16.*T2new(iz)+T1new(iz)))

           D3(iz)=D3(iz)-N1new(iz)*qenew(iz)*Cej*xne_1(iz)*nuje(iz)*
     &             (1.5-24.*T1new(iz)/(16.*Tenew(iz)+5.46e-4*T1new(iz)))

           D3(iz)=D3(iz)-N1new(iz)*q3new(iz)*Cnj*xn3_1(iz)*nujn(iz)*
     &             (.128-24.*T1new(iz)/(16.*T3new(iz)+14.*T1new(iz)))

           D3(iz)=D3(iz)+N1new(iz)*(
     &                           nujN2(iz)*(.139*q_Nn2(iz)*Njqj0
     &                                    +1.545*T1new(iz)*Un(iz)/Cj0)
     &                          +nujO2(iz)*(.2*q_No2(iz)*Njqj0
     &                                    +1.5*T1new(iz)*Un(iz)/Cj0)
     &                          +nujO(iz)* (-.075*q_No(iz)*Njqj0
     &                                    +1.75*T1new(iz)*Un(iz)/Cj0)
     &                          )
     &           -2.5*T1new(iz)*N1new(iz)*(
     &                         (nujN2(iz)+nujO2(iz)+nujO(iz))*(Un(iz)/Cj0)
     &               -9.6*Njqj0*
     &                  (q_Nn2(iz)*nujN2(iz)/(16.*Tn(iz)/T_0+28.*T1new(iz))
     &                  +q_No2(iz)*nujO2(iz)/(16.*Tn(iz)/T_0+32.*T1new(iz))
     &                  +q_No(iz)*nujO (iz)/(16.*Tn(iz)/T_0+16.*T1new(iz))
     &                  )               )

          D7(iz)=-(nujN2(iz)*.34+nujO2(iz)*.27+nujO(iz)*.725
     &           +nuji(iz)*2.66+nujk(iz)*.34+nujl(iz)*.27
     &           +nujm(iz)*.3+nujn(iz)*.835+nujj(iz)*.8+nuje(iz)*3.)
     &          -1.5*T1new(iz)*(
     &                        nujN2(iz)*28./(16.*Tn(iz)/T_0+28.*T1new(iz))
     &                       +nujO2(iz)*32./(16.*Tn(iz)/T_0+32.*T1new(iz))
     &                       +nujO (iz)*16./(16.*Tn(iz)/T_0+16.*T1new(iz))
     &                       +nuji(iz)/(16.*T2new(iz)+T1new(iz))
     &                       +nujk(iz)*28./(16.*Tmnew(iz)+28.*T1new(iz))
     &                       +nujl(iz)*32./(16.*Tmnew(iz)+32.*T1new(iz))
     &                       +nujm(iz)*30./(16.*Tmnew(iz)+30.*T1new(iz))
     &                       +nujn(iz)*14./(16.*T3new(iz)+14.*T1new(iz))
     &                       +nuje(iz)*5.46e-4
     &                               /(16.*Tenew(iz)+5.46e-4*T1new(iz)))
     &          -4.2*Cji*U1new(iz)*alt_geo_1(iz)

	 enddo



C]]]

C[[[    H+ heat flow equation resolution

	  do iz=1,nb_alt
	
	  D3(iz)=alt_geo_1(iz)*N2new(iz)*(T2pnew(iz)-T2tnew(iz))
     &		                   *(T2pnew(iz)-4.*T2tnew(iz))
     &	         -N2new(iz)*T2new(iz)*U2new(iz)
     &                     *(nuij(iz)*18.5/17.
     &                      +nuik(iz)*30.5/29.
     &                      +nuil(iz)*34.5/33.
     &                      +nuim(iz)*32.5/31.
     &                      +nuin(iz)*16.5/15.
     &                      +nuie(iz)*2.5
     &                      +nuiN2(iz)*30.5/29.
     &                      +nuiO2(iz)*34.5/33.
     &                      +nuiO(iz)*18.5/17.
     &                      +nuiH(iz)*3.5/2.)
     &            + 2.5*N2new(iz)*T2new(iz)*U2new(iz)*(
     &                       nuiN2(iz)+nuiO2(iz)+nuiO(iz)+nuiH(iz)
     &                      +nuij(iz)
     &                      +nuik(iz)+nuil(iz)+nuim(iz)
     &                      +nuin(iz)
     &                      +nuie(iz))

            D3(iz)=D3(iz)+N2new(iz)*T2new(iz)*U1new(iz)*Cji*
     &                           (nuij(iz)*18.5/17.
     &                           -2.5*nuij(iz))

            D3(iz)=D3(iz)+N2new(iz)*T2new(iz)*Umnew(iz)*Czi*
     &                           (nuik(iz)*30.5/29.
     &                           +nuil(iz)*34.5/33.+nuim(iz)*32.5/31.
     &                           -2.5*(nuik(iz)+nuil(iz)+nuim(iz))     )

            D3(iz)=D3(iz)+N2new(iz)*T2new(iz)*U3new(iz)*Cni*
     &                           (nuin(iz)*16.5/15.
     &                           -2.5*nuin(iz)         )

            D3(iz)=D3(iz)+N2new(iz)*q1new(iz)*Cji*xn1_1(iz)*
     &            nuij(iz)*(.0612+1.5*T2new(iz)/(T1new(iz)+16.*T2new(iz)))

            D3(iz)=D3(iz)+N2new(iz)*qenew(iz)*Cei*xne_1(iz)*
     &            nuie(iz)*1.5*(T2new(iz)/(Tenew(iz)+5.46e-4*T2new(iz))-1)

            D3(iz)=D3(iz)+N2new(iz)*q3new(iz)*Cni*xn3_1(iz)*
     &            nuin(iz)*(.068+1.5*T2new(iz)/(T3new(iz)+14.*T2new(iz)))

            D3(iz)=D3(iz)+N2new(iz)*(
     &                            nuiN2(iz)*(1.07*q_Nn2(iz)*Niqi0
     &                                     +1.052*T2new(iz)*Un(iz)/Ci0)
     &                           +nuiO2(iz)*(1.084*q_No2(iz)*Niqi0
     &                                     +1.045*T2new(iz)*Un(iz)/Ci0)
     &                           +nuiO(iz) *(.98*q_No(iz)*Niqi0
     &                                     +1.088*T2new(iz)*Un(iz)/Ci0)
     &                           +nuiH(iz) *(-.075*q_Nh(iz)*Niqi0
     &                                     +1.75*T2new(iz)*Un(iz)/Ci0 )
     &                           )
     &          -2.5*T2new(iz)*N2new(iz)*        (
     &                    (nuiN2(iz)+nuiO2(iz)+nuiO(iz)+nuiH(iz))
     &                   *(Un(iz)/Ci0)
     &                   -.6*(nuiN2(iz)/(Tn(iz)/T_0+28.*T2new(iz))
     &                       *(q_Nn2(iz)*Niqi0)
     &                       +nuiO2(iz)/(Tn(iz)/T_0+32.*T2new(iz))
     &                       *(q_No2(iz)*Niqi0)
     &                       +nuiO (iz)/(Tn(iz)/T_0+16.*T2new(iz))
     &                       *(q_No(iz)*Niqi0)
     &                       +nuiH (iz)/(Tn(iz)/T_0+    T2new(iz))
     &                       *(q_Nh(iz)*Niqi0))  )

            D7(iz)=
     &            (
     &              nuiN2(iz)*.1795+nuiO2(iz)*.1824
     &              +nuiO(iz)*.1612-nuiH(iz)*.725
     &              -nuii(iz)*.8+nuij(iz)*.1612+nuik(iz)*.1795
     &              +nuil(iz)*.1824+nuim(iz)*.1811+nuin(iz)*.1547
     &              -nuie(iz)*3.
     &            )
     &            -1.5*T2new(iz)*(
     &                           nuiN2(iz)*28./(Tn(iz)/T_0+28.*T2new(iz))
     &                          +nuiO2(iz)*32./(Tn(iz)/T_0+32.*T2new(iz))
     &                          +nuiO (iz)*16./(Tn(iz)/T_0+16.*T2new(iz))
     &                          +nuiH (iz)/(Tn(iz)/T_0+T2new(iz))
     &                          +nuij(iz)*16./(T1new(iz)+16.*T2new(iz))
     &                          +nuik(iz)*28./(Tmnew(iz)+28.*T2new(iz))
     &                          +nuil(iz)*32./(Tmnew(iz)+32.*T2new(iz))
     &                          +nuim(iz)*30./(Tmnew(iz)+30.*T2new(iz))
     &                          +nuin(iz)*14./(T3new(iz)+14.*T2new(iz))
     &                          +nuie(iz)*5.46e-4/
     &                               (Tenew(iz)+5.46e-4*T2new(iz))     )
     &            -4.2*U2new(iz)*alt_geo_1(iz)

	  enddo



C]]]



	  call velocity(Velnq,Ipos1,Iposnp,deltat_2)

C[[[    N+ heat flow equation resolution

	 do iz=1,nb_alt



	  D3(iz)=Cni*alt_geo_1(iz)*N3new(iz)*(T3pnew(iz)-T3tnew(iz))
     &		                   *(T3pnew(iz)-4.*T3tnew(iz))
     &	        -.1*N3new(iz)*U2new(iz)*Cin*T3new(iz)*nuni(iz)
     &          -.8*N3new(iz)*U1new(iz)*Cjn*T3new(iz)*nunj(iz)
     &          -   N3new(iz)*Umnew(iz)*Czn*T3new(iz)
     &                        *(nunk(iz)+1.04*nunl(iz)+1.02*nunm(iz))

          D3(iz)=D3(iz)+N3new(iz)*T3new(iz)*U3new(iz)*               (
     &                         nunN2(iz)+nunO2(iz)*1.04+nunO(iz)*.8
     &                        +nuni(iz)*.1+nunj(iz)*.8+nunk(iz)
     &                        +nunl(iz)*1.04+nunm(iz)*1.02        )

          D3(iz)=D3(iz)+N3new(iz)*q2new(iz)*Cin*xn2_1(iz)*nuni(iz)*
     &               (21.*T3new(iz)/(14.*T2new(iz)+T3new(iz))-1.232)

          D3(iz)=D3(iz)+N3new(iz)*q1new(iz)*Cjn*xn1_1(iz)*nunj(iz)*
     &               (21.*T3new(iz)/(14.*T1new(iz)+16.*T3new(iz))-.028)

          D3(iz)=D3(iz)+N3new(iz)*qenew(iz)*Cen*xne_1(iz)*nune(iz)*
     &               (21.*T3new(iz)/(14.*Tenew(iz)+5.46e-4*T3new(iz))-1.5)

          D3(iz)=D3(iz)+N3new(iz)*  (
     &            nunN2(iz)* (.1*q_Nn2(iz)*Nnqn0
     &                      +1.5*T3new(iz)*Un(iz)/Cn0)
     &            +nunO2(iz)*(.115*q_No2(iz)*Nnqn0
     &                      +1.456*T3new(iz)*Un(iz)/Cn0)
     &            +nunO(iz)* (-.028*q_No(iz)*Nnqn0
     &                      +1.7*T3new(iz)*Un(iz)/Cn0)
     &                           )
     &          -2.5*T3new(iz)*N3new(iz)*(
     &                        (nunN2(iz)+nunO2(iz)+nunO(iz))*(Un(iz)/Cn0)
     &           -8.4*Nnqn0*(
     &                  nunN2(iz)*q_Nn2(iz)/(14.*Tn(iz)/T_0+28.*T3new(iz))
     &                 +nunO2(iz)*q_No2(iz)/(14.*Tn(iz)/T_0+32.*T3new(iz))
     &                 +nunO (iz)*q_No(iz)/(14.*Tn(iz)/T_0 +16.*T3new(iz)))
     &                                 )


          D7(iz)=-( nunN2(iz)*.27+nunO2(iz)*.2+nunO(iz)*.62+nuni(iz)*2.62
     &            +nunk(iz)*.27+nunl(iz)*.2 +nunm(iz)*.23
     &            +nunj(iz)*.62+nunn(iz)*.8 +nune(iz)*3.)
     &            -1.5*T3new(iz)*(
     &                 nunN2(iz)*28./(14.*Tn(iz)/T_0+28.*T3new(iz))
     &                +nunO2(iz)*32./(14.*Tn(iz)/T_0+32.*T3new(iz))
     &                +nunO (iz)*16./(14.*Tn(iz)/T_0+16.*T3new(iz))
     &                +nuni(iz)/(14.*T2new(iz)+T3new(iz))
     &                +nunk(iz)*28./(14.*Tmnew(iz)+28.*T3new(iz))
     &                +nunl(iz)*32./(14.*Tmnew(iz)+32.*T3new(iz))
     &                +nunm(iz)*30./(14.*Tmnew(iz)+30.*T3new(iz))
     &                +nunj(iz)*16./(14.*T1new(iz)+16.*T3new(iz))
     &                +nune(iz)*5.46e-4/(14.*Tenew(iz)+5.46e-4*T3new(iz)))
     &          -4.2*Cni*U3new(iz)*alt_geo_1(iz)

	 enddo



C]]]

	
	
C[[[    H+ momentum equation resolution

	  call velocity(Velim,Ipos1,Iposnp,deltat_2)
	
	  do iz=1,nb_alt


	    D3(iz)=-G(iz)+thermacc(iz)/mi/Ci0
     &            -3.*(Tepnew(iz)-Tetnew(iz))*alt_geo_1(iz)
     &            -3.*(T2pnew(iz)-T2tnew(iz))*alt_geo_1(iz)
     &		  +nuij(iz)*Cji*U1new(iz)
     &		  +(nuik(iz)*Czi+nuil(iz)*Czi+nuim(iz)*Czi)*Umnew(iz)
     &		  +nuin(iz)*Cni*U3new(iz)

	    D3(iz)=D3(iz)+.6*q2new(iz)*xn2_1(iz)
     &		       *(28.*nuiN2(iz)/(Tn(iz)/T_0+28.*T2new(iz))
     &                  +32.*nuiO2(iz)/(Tn(iz)/T_0+32.*T2new(iz))
     &                  +16.*nuiO(iz)/(Tn(iz)/T_0+16.*T2new(iz))
     &                  +    nuiH(iz)/(Tn(iz)/T_0+T2new(iz))
     &                  +16.*nuij(iz)/(T1new(iz)+16.*T2new(iz))
     &                  +28.*nuik(iz)/(Tmnew(iz)+28.*T2new(iz))
     &                  +32.*nuil(iz)/(Tmnew(iz)+32.*T2new(iz))
     &                  +30.*nuim(iz)/(Tmnew(iz)+30.*T2new(iz))
     &                  +14.*nuin(iz)/(T3new(iz)+14.*T2new(iz))
     &                  +5.46e-4*nuie(iz)/(Tenew(iz)+5.46e-4*T2new(iz)))

	    D3(iz)=D3(iz)-.6*q1new(iz)*xn1_1(iz)*Cji
     &			  *nuij(iz)/(T1new(iz)+16.*T2new(iz))

	    D3(iz)=D3(iz)-.6*qenew(iz)*xne_1(iz)*Cei
     &			  *nuie(iz)/(Tenew(iz)+5.46e-4*T2new(iz))

	    D3(iz)=D3(iz)-.6*q3new(iz)*xn3_1(iz)*Cni
     &			  *nuin(iz)/(T3new(iz)+14.*T2new(iz))
	
	    D3(iz)=D3(iz)+(nuiN2(iz)+nuiO2(iz)+nuiO(iz)+nuiH(iz))
     &                    *(Un(iz)/Ci0)
     &                 -.6*nuiN2(iz)/(Tn(iz)/T_0+28.*T2new(iz))
     &                    *(q_Nn2(iz)*Niqi0)
     &                 -.6*nuiO2(iz)/(Tn(iz)/T_0+32.*T2new(iz))
     &                    *(q_No2(iz)*Niqi0)
     &                 -.6*nuiO(iz)/(Tn(iz)/T_0+16.*T2new(iz))
     &                    *(q_No(iz)*Niqi0)
     &                 -.6*nuiH(iz)/(Tn(iz)/T_0+T2new(iz))
     &                    *(q_Nh(iz)*Niqi0)
     &                 +nuie(iz)*Cei*Uenew(iz)

	    D7(iz)=-(nuiN2(iz)+nuiO2(iz)+nuiO(iz)+nuiH(iz)
     &              +nuij(iz)+nuik(iz)+nuil(iz)+nuim(iz)
     &              +nuin(iz)+nuie(iz))

	  enddo


C]]]

C[[[    O+ momentum equation resolution

	
	  do iz=1,nb_alt

	    D3(iz)=-Cij*G(iz)+thermacc(iz)/mj/Cj0
     &            -3.*Cji*(Tepnew(iz)-Tetnew(iz))*alt_geo_1(iz)
     &            -3.*Cji*(T1pnew(iz)-T1tnew(iz))*alt_geo_1(iz)
     &		  +nuji(iz)*Cij*U2new(iz)

	    D3(iz)=D3(iz)+(nujk(iz)+nujl(iz)+nujm(iz))*Czj*Umnew(iz)
     &		       +nujn(iz)*Cnj*U3new(iz)

	    D3(iz)=D3(iz)-.6*q2new(iz)*xn2_1(iz)*Cij*16.
     &			*nuji(iz)/(16.*T2new(iz)+T1new(iz))

	    D3(iz)=D3(iz)+.6*q1new(iz)*xn1_1(iz)
     &			*(28.*nujN2(iz)/(16.*Tn(iz)/T_0+28.*T1new(iz))
     &                   +32.*nujO2(iz)/(16.*Tn(iz)/T_0+32.*T1new(iz))
     &                   +16.*nujO(iz)/(16.*Tn(iz)/T_0+16.*T1new(iz))
     &                   +    nuji(iz)/(16.*T2new(iz)+T1new(iz))
     &                   +28.*nujk(iz)/(16.*Tmnew(iz)+28.*T1new(iz))
     &                   +32.*nujl(iz)/(16.*Tmnew(iz)+32.*T1new(iz))
     &                   +30.*nujm(iz)/(16.*Tmnew(iz)+30.*T1new(iz))
     &                   +14.*nujn(iz)/(16.*T3new(iz)+14.*T1new(iz))
     &                   +    nuje(iz)/(16.*Tenew(iz)+5.46e-4*T1new(iz))
     &                        *5.46e-4)

	    D3(iz)=D3(iz)-.6*qenew(iz)*xne_1(iz)*Cej*16.
     &			*nuje(iz)/(16.*Tenew(iz)+5.46e-4*T1new(iz))

	    D3(iz)=D3(iz)-.6*q3new(iz)*xn3_1(iz)*Cnj*16.
     &			*nujn(iz)/(16.*T3new(iz)+14.*T1new(iz))

	    D3(iz)=D3(iz)+(nujN2(iz)+nujO2(iz)+nujO(iz))*(Un(iz)/Cj0)
     &                 -.6*nujN2(iz)/(16.*Tn(iz)/T_0+28.*T1new(iz))
     &                  *(16.*q_Nn2(iz)*Njqj0)
     &                 -.6*nujO2(iz)/(16.*Tn(iz)/T_0+32.*T1new(iz))
     &                  *(16.*q_No2(iz)*Njqj0)
     &                 -.6*nujO(iz)/(16.*Tn(iz)/T_0+16.*T1new(iz))
     &                  *(16.*q_No(iz)*Njqj0)
     &                 +nuje(iz)*Cje*Uenew(iz)

	    D7(iz)=-(nujN2(iz)+nujO2(iz)+nujO(iz)
     &              +nuji(iz)+nujn(iz)
     &              +nuje(iz)+nujk(iz)+nujl(iz)+nujm(iz))

	  enddo




C]]]

C[[[    heavy ions momentum equation resolution

	
	  do iz=1,nb_alt



          D3(iz)=-Ciz*G(iz)+thermacc(iz)/mz/Cz0
     &            -3.*Czi*(Tepnew(iz)-Tetnew(iz))*alt_geo_1(iz)
     &            -3.*Czi*(Tmpnew(iz)-Tmtnew(iz))*alt_geo_1(iz)
     &		 +(nukj(iz)+nulj(iz)+numj(iz))*Cjz/3.*U1new(iz)
     &		 +(nuki(iz)+nuli(iz)+numi(iz))*Ciz/3.*U2new(iz)
     &		 +(nukn(iz)+nuln(iz)+numn(iz))*Cnz/3.*U3new(iz)

	  D3(iz)=D3(iz)-.2*q2new(iz)*xn2_1(iz)*Ciz
     &		      *(28.*nuki(iz)/(28.*T2new(iz)+Tmnew(iz))
     &                 +32.*nuli(iz)/(32.*T2new(iz)+Tmnew(iz))
     &                 +30.*numi(iz)/(30.*T2new(iz)+Tmnew(iz)))

	  D3(iz)=D3(iz)-.2*q1new(iz)*xn1_1(iz)*Cjz
     &		      *(28.*nukj(iz)/(28.*T1new(iz)+16.*Tmnew(iz))
     &                 +32.*nulj(iz)/(32.*T1new(iz)+16.*Tmnew(iz))
     &                 +30.*numj(iz)/(30.*T1new(iz)+16.*Tmnew(iz)))

	  D3(iz)=D3(iz)-.2*qenew(iz)*xne_1(iz)*Cez
     &		      *(28.*nuke(iz)/(28.*Tenew(iz)+5.46e-4*Tmnew(iz))
     &                 +32.*nule(iz)/(32.*Tenew(iz)+5.46e-4*Tmnew(iz))
     &                 +30.*nume(iz)/(30.*Tenew(iz)+5.46e-4*Tmnew(iz)))

	  D3(iz)=D3(iz)-.2*q3new(iz)*xn3_1(iz)*Cnz
     &		      *(28.*nukn(iz)/(28.*T3new(iz)+14.*Tmnew(iz))
     &                 +32.*nuln(iz)/(32.*T3new(iz)+14.*Tmnew(iz))
     &                 +30.*numn(iz)/(30.*T3new(iz)+14.*Tmnew(iz)))

	  D3(iz)=D3(iz)+(nukN2(iz)+nukO2(iz)+nukO(iz)
     &               +nulN2(iz)+nulO2(iz)+nulO(iz)
     &               +numN2(iz)+numO2(iz)+numO(iz))*Un(iz)/Cz0/3.
     &           -.2*
     &               (nukN2(iz)/(28.*Tn(iz)/T_0+28.*Tmnew(iz))
     &           *(28.*q_Nn2(iz)*Nkqk0)
     &               +nukO2(iz)/(28.*Tn(iz)/T_0+32.*Tmnew(iz))
     &           *(28.*q_No2(iz)*Nkqk0)
     &               +nukO(iz)/(28.*Tn(iz)/T_0+16.*Tmnew(iz))
     &           *(28.*q_No(iz)*Nkqk0))
     &           -.2*
     &               (nulN2(iz)/(32.*Tn(iz)/T_0+28.*Tmnew(iz))
     &           *(32.*q_Nn2(iz)*Nlql0)
     &               +nulO2(iz)/(32.*Tn(iz)/T_0+32.*Tmnew(iz))
     &           *(32.*q_No2(iz)*Nlql0)
     &               +nulO(iz)/(32.*Tn(iz)/T_0+16.*Tmnew(iz))
     &           *(32.*q_No(iz)*Nlql0))
     &           -.2*
     &               (numN2(iz)/(30.*Tn(iz)/T_0+28.*Tmnew(iz))
     &           *(30.*q_Nn2(iz)*Nmqm0)
     &                +numO2(iz)/(30.*Tn(iz)/T_0+32.*Tmnew(iz))
     &           *(30.*q_No2(iz)*Nmqm0)
     &                +numO(iz)/(30.*Tn(iz)/T_0+16.*Tmnew(iz))
     &           *(30.*q_No(iz)*Nmqm0))
     &           +Uenew(iz)*Cez*(nuke(iz)+nule(iz)+nume(iz))/3.

          D7(iz)=-(nukN2(iz)+nukO2(iz)+nukO(iz)
     &            +nuki(iz)+nuke(iz)+nukj(iz)+nukn(iz)
     &            +nulN2(iz)+nulO2(iz)+nulO(iz)
     &            +nuli(iz)+nule(iz)+nulj(iz)+nuln(iz)
     &            +numN2(iz)+numO2(iz)+numO(iz)
     &            +numi(iz)+nume(iz)+numj(iz)+numn(iz))/3.
	

	  enddo



C[[[[   Boundaries conditions

C]]]

C[[[    N+ momentum equation resolution
	
	
	  do iz=1,nb_alt


	    D3(iz)=-Cin*G(iz)+thermacc(iz)/mn/Cn0
     &            -3.*Cni*(Tepnew(iz)-Tetnew(iz))*alt_geo_1(iz)
     &            -3.*Cni*(T3pnew(iz)-T3tnew(iz))*alt_geo_1(iz)
     &		  +nuni(iz)*Cin*U2new(iz)+nunj(iz)*Cjn*U1new(iz)
     &		  +(nunk(iz)+nunl(iz)+nunm(iz))*Czn*Umnew(iz)
	
	  D3(iz)=D3(iz)-.6*q2new(iz)*xn2_1(iz)*Cin*14.
     &			*nuni(iz)/(14.*T2new(iz)+T3new(iz))


	  D3(iz)=D3(iz)-.6*q1new(iz)*xn1_1(iz)*Cjn*14.
     &			*nunj(iz)/(14.*T1new(iz)+16.*T3new(iz))


	  D3(iz)=D3(iz)-.6*qenew(iz)*xne_1(iz)*Cen*14.
     &			*nune(iz)/(14.*Tenew(iz)+5.46e-4*T3new(iz))


	  D3(iz)=D3(iz)+.6*q3new(iz)*xn3_1(iz)
     &			*(28.*nunN2(iz)/(14.*Tn(iz)/T_0+28.*T3new(iz))
     &                  +32.*nunO2(iz)/(14.*Tn(iz)/T_0+32.*T3new(iz))
     &                  +16.*nujO(iz)/(14.*Tn(iz)/T_0+16.*T3new(iz))
     &                  +    nuni(iz)/(14.*T2new(iz)+T3new(iz))
     &                  +28.*nunk(iz)/(14.*Tmnew(iz)+28.*T3new(iz))
     &                  +32.*nunl(iz)/(14.*Tmnew(iz)+32.*T3new(iz))
     &                  +30.*nunm(iz)/(14.*Tmnew(iz)+30.*T3new(iz))
     &                  +    nune(iz)/(14.*Tenew(iz)+5.46e-4*T3new(iz))
     &                       *5.46e-4)

	  D3(iz)=D3(iz)+(nunN2(iz)+nunO2(iz)+nunO(iz))*Un(iz)/Cn0
     &           -.6*nunN2(iz)/(14.*Tn(iz)/T_0+28.*T3new(iz))
     &           *(14.*q_Nn2(iz)*Nnqn0)
     &           -.6*nunO2(iz)/(14.*Tn(iz)/T_0+32.*T3new(iz))
     &           *(14.*q_No2(iz)*Nnqn0)
     &           -.6*nunO(iz)/(14.*Tn(iz)/T_0+16.*T3new(iz))
     &           *(14.*q_No(iz)*Nnqn0)
     &           +nune(iz)*Cen*Uenew(iz)

	    D7(iz)=-(nunN2(iz)+nunO2(iz)+nunO(iz)
     &             +nuni(iz)+nunj(iz)+nune(iz)
     &             +nunk(iz)+nunl(iz)+nunm(iz))

	  enddo


C]]]
