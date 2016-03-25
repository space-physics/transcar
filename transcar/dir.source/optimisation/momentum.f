

!----------------------------------!
! equation d'impulsion de l'ion H+ !
!----------------------------------!

	  call velocity(Velim,Ipos1,Iposnp,deltat_2)

	  do i=1,nx

          Tr1=T_0*T2new(i)
          Tr=(Tn(i)+Tr1)/2.

	  T1_15(i)=T1new(i)**(-1.5)/nu_0
	  T2_15(i)=T2new(i)**(-1.5)/nu_0
	  T3_15(i)=T3new(i)**(-1.5)/nu_0
	  Tm_15(i)=Tmnew(i)**(-1.5)/nu_0
	  Te_15(i)=Tenew(i)**(-1.5)/nu_0

	  thermacc(i)=thermodiff*qenew(i)/Tenew(i)**2.5

          nuiN2(i)=3.36e-9*Nn2(i)                               *t0

          nuiO2(i)=3.20e-9*No2(i)                               *t0
          nuiO (i)=6.61e-11*No (i)
     &             *(1.-.047*alog10(Tr1))**2*sqrt(Tr1)          *t0
          nuiH (i)=2.65e-10*Nh (i)
     &             *(1.-.083*alog10(Tr))**2*sqrt(Tr)            *t0
          nuij (i) = 1.23*N1new(i)
     &              /((16.*T2new(i)+T1new(i))/17.)**1.5		/nu_0
          nuii (i) = .9*N2new(i)
     &              /(T2new(i))**1.5 				/nu_0
          nuin (i) = 1.23*N3new(i)
     &              /((14.*T2new(i)+T3new(i))/15.)**1.5 	/nu_0
          nuik (i) = 1.25*N4new(i)
     &              /((28.*T2new(i)+Tmnew(i))/27.)**1.5		/nu_0
          nuim (i) = 1.25*N5new(i)
     &              /((30.*T2new(i)+Tmnew(i))/31.)**1.5  	/nu_0
          nuil (i) = 1.25*N6new(i)
     &              /((32.*T2new(i)+Tmnew(i))/33.)**1.5		/nu_0
          nuie (i)=0.03    *Nenew(i)*Te_15(i)



	    C2a(i)=-Tepnew(i)
	    D2a(i)=log(Nenew(i)*Tepnew(i))
            C2b(i)=-T2pnew(i)*N2new(i)/xn2(i)
	    D2b(i)=log(xn2(i)*T2pnew(i))
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
          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)

	  D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
	  D2br=.5*log(xn2(np)*T2pnew(np)*xn2(nx)*T2pnew(nx))
	  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
	  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c	  call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
	  call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=1.
	  call lcpfct(U2old,U2new,Ipos1,Iposn,
     &		      lbc,0.,0.,U2new(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif

            U2new(i)=dexpnu*(U2new(i)-U2old(i)+D3(i)*deltat_2)
     &		     +expnu*U2old(i)
          enddo

          if (isnant(U2new,nx)) then
            print*,'probleme lors du calcul de U2new dans la boucle 1'
            goto 246
          endif


C]]]

C[[[    O+ momentum equation resolution

	  call velocity(Veljm,Ipos1,Iposnp,deltat_2)

	  do i=1,nx

          Tr=(T_0*T1new(i)+Tn(i))/2.

          nujN2(i)=6.82e-10*Nn2(i)                              *t0
          nujO2(i)=6.64e-10*No2(i)                              *t0
          nujO (i)=3.67e-11*No(i)*sqrt(Tr)
     &             *(1.-.064*alog10(Tr))**2                     *t0

          nujj(i)=.22*N1new(i)
     &           /(T1new(i))**1.5				/nu_0
          nuji(i)=0.077*N2new(i)
     &           /((T1new(i)+16.*T2new(i))/17.)**1.5		/nu_0
          nujn(i)=0.22*N3new(i)
     &           /((14.*T1new(i)+16.*T3new(i))/30.)**1.5 	/nu_0
          nujk(i)=0.25*N4new(i)
     &           /((28.*T1new(i)+16.*Tmnew(i))/44.)**1.5	/nu_0
          nujm(i)=0.26*N5new(i)
     &           /((30.*T1new(i)+16.*Tmnew(i))/46.)**1.5	/nu_0
          nujl(i)=0.26*N6new(i)
     &	         /((32.*T1new(i)+16.*Tmnew(i))/48.)**1.5	/nu_0
          nuje(i)=1.87e-3*Nenew(i)*Te_15(i)



	    C2a(i)=-Cji*Tepnew(i)
	    D2a(i)=log(Nenew(i)*Tepnew(i))
            C2b(i)=-Cji*T1pnew(i)*N1new(i)/xn1(i)
            D2b(i)=log(xn1(i)*T1pnew(i))
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

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)

	  D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
	  D2br=.5*log(xn1(np)*T1pnew(np)*xn1(nx)*T1pnew(nx))

	  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
	  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c	  call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
	  call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=1.
	  call lcpfct(U1old,U1new,Ipos1,Iposn,
     &		      lbc,0.,0.,U2new(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            U1new(i)=dexpnu*(U1new(i)-U1old(i)+D3(i)*deltat_2)
     &		     +expnu*U1old(i)

          enddo
c	  U1new(nx)=max(U1new(nx),-1000./Cj0)
c	  U1new(nx)=max(U1new(nx),0.)

          if (isnant(U1new,nx)) then
            print*,'probleme lors du calcul de U1new dans la boucle 1'
            goto 246
          endif




C]]]

C[[[    heavy ions momentum equation resolution

	  call velocity(Velmm,Ipos1,Iposnp,deltat_2)

	  do i=1,nx

          Tr=(T_0*Tmnew(i)+Tn(i))/2.

          nukH(i)=0.074e-9*Nh(i)                                *t0
          nukO2(i)=0.449e-9*No2(i)                              *t0
          nukO (i)=0.258e-9*No (i)                              *t0
          nukN2 (i)=5.14e-11*Nn2 (i)
     &             *(1.-.069*alog10(Tr))**2*sqrt(Tr)            *t0

          nukj (i)=0.15    *N1new(i)*T1_15(i)
          nuki (i)=0.045   *N2new(i)*T2_15(i)
          nukn (i)=0.14    *N3new(i)*T3_15(i)
          nukk (i)=0.17    *N4new(i)*Tm_15(i)
          nukm (i)=0.17    *N5new(i)*Tm_15(i)
          nukl (i)=0.18    *N6new(i)*Tm_15(i)
          nuke (i)=1.06e-3 *Nenew(i)*Tm_15(i)

          nulH(i)=0.065e-9*Nh(i)                                *t0
          nulN2(i)=0.413e-9*Nn2(i)                              *t0
          nulO (i)=0.231e-9*No (i)                              *t0
          nulO2 (i)=2.59e-11*No2 (i)
     &             *(1.-.073*alog10(Tr))**2*sqrt(Tr)            *t0
          nulj (i)=0.13    *N1new(i)*T1_15(i)
          nuli (i)=0.039   *N2new(i)*T2_15(i)
          nuln (i)=0.12    *N3new(i)*T3_15(i)
          nulk (i)=0.15    *N4new(i)*Tm_15(i)
          nulm (i)=0.16    *N5new(i)*Tm_15(i)
          null (i)=0.16    *N6new(i)*Tm_15(i)
          nule (i)=9.276e-4*Nenew(i)*Te_15(i)

          numH(i)=0.069e-9*Nh(i)                                *t0
          numN2(i)=0.434e-9*Nn2(i)                              *t0
          numO (i)=0.244e-9*No (i)                              *t0
          numO2(i)=0.427e-9*No2 (i)                             *t0
          numj (i)=0.14    *N1new(i)*T1_15(i)
          numi (i)=0.042   *N2new(i)*T2_15(i)
          numn (i)=0.13    *N3new(i)*T3_15(i)
          numk (i)=0.16    *N4new(i)*Tm_15(i)
          numm (i)=0.16    *N5new(i)*Tm_15(i)
          numl (i)=0.17    *N6new(i)*Tm_15(i)
          nume (i)=9.894e-4*Nenew(i)*Tm_15(i)


          C2a(i)=-Czi*Tepnew(i)
	  D2a(i)=log(Nenew(i)*Tepnew(i))

          C2b(i)=-Czi*Tmpnew(i)*Nmnew(i)/xnm(i)
          D2b(i)= log(xnm(i)*Tmpnew(i))

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

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)

	  D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
	  D2br=.5*log(xnm(np)*Tmpnew(np)*xnm(nx)*Tmpnew(nx))

	  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
	  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c	  call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
	  call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U1old,extra,nx)
	  rbc=1.
	  call lcpfct(Umold,Umnew,Ipos1,Iposn,
     &		      lbc,0.,0.,Umnew(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Umnew(i)=dexpnu*(Umnew(i)-Umold(i)+D3(i)*deltat_2)
     &		     +expnu*Umold(i)
          enddo

          if (isnant(Umnew,nx)) then
            print*,'probleme lors du calcul de Umnew dans la boucle 1'
            goto 246
          endif


C[[[[   Boundaries conditions

C]]]

C[[[    N+ momentum equation resolution

	  call velocity(Velnm,Ipos1,Iposnp,deltat_2)

	  do i=1,nx

          nunH (i)=0.145e-9*Nh(i)                               *t0
          nunN2(i)=0.747e-9*Nn2(i)                              *t0
          nunO (i)=0.442e-9*No (i)                              *t0
          nunO2(i)=0.725e-9*No2 (i)                             *t0
          nunj(i)=0.25*N1new(i)
     &           / ((16.*T3new(i)+14.*T1new(i))/30.)**1.5	/nu_0
          nuni(i)=0.088*N2new(i)
     &            / (( T3new(i)+14.*T2new(i))/15.)**1.5		/nu_0
          nunn(i)=0.24*N3new(i)/ (T3new(i))**1.5		/nu_0
          nunk(i)=0.28*N4new(i)
     &            / ((28.*T3new(i)+14.*Tmnew(i))/42.)**1.5	/nu_0
          nunm(i)=0.28*N5new(i)
     &            / ((30.*T3new(i)+14.*Tmnew(i))/44.)**1.5	/nu_0
          nunl(i)=0.28*N6new(i)
     &            / ((32.*T3new(i)+14.*Tmnew(i))/46.)**1.5	/nu_0
          nune(i)=2.136e-3*Nenew(i)*Te_15(i)


	    C2a(i)=-Cni*Tepnew(i)
	    D2a(i)=log(Nenew(i)*Tepnew(i))
	    C2b(i)=-Cni*T3pnew(i)*N3new(i)/xn3(i)
            D2b(i)=log(xn3(i)*T3pnew(i))


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

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)

	  D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
	  D2br=.5*log(xn3(np)*T3pnew(np)*xn3(nx)*T3pnew(nx))

	  call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
	  call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c	  call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
	  call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U3old,extra,nx)
	  rbc=1.
	  call lcpfct(U3old,U3new,Ipos1,Iposn,
     &		      lbc,0.,0.,U3new(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            U3new(i)=dexpnu*(U3new(i)-U3old(i)+D3(i)*deltat_2)
     &		     +expnu*U3old(i)
          enddo
c	  U3new(nx)=max(U3new(nx),0./Cn0)
c	flag=.false.
c	do i=1,nx
c	  U3new(i)=U1new(i)
c	enddo

          if (isnant(U3new,nx)) then
            print*,'probleme lors du calcul de U3new dans la boucle 1'
            goto 246
          endif


C]]]


C[[[    Velocities corrections

	  do i=1,nx

	 Uenew(i)=(Ci0*N2new(i)*U2new(i)
     &           +Cj0*N1new(i)*U1new(i)
     &           +Cn0*N3new(i)*U3new(i)
     &           +Cz0*(N4new(i)+N6new(i)+N5new(i))*Umnew(i)
     &            -JJ(i)/N_0-0.*Jes(i)/N_0)
     &               /Nenew(i)/Ce0

          enddo

