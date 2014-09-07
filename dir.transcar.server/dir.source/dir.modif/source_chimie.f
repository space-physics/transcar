	real chim(nb_ion,nb_ion),y0(nb_ion),y1(nb_ion)



!-------------------------------------!
!       equations de continuite       !
!-------------------------------------!

!	ion O+

	  call velocity(Veljc,Ipos1,Iposnp,deltat_2)
          lbc=sqrt(N1old(2)/N1old(3))
	  lbc=1.
	  rbc=(N1old(nx)/N1old(nx-1))
	  call lcpfct(N1old,N1new,Ipos1,Iposn,
     &		      lbc,0.,0.,N1new(np),.false.,1)
	
!	ion H+

	  call velocity(Velic,Ipos1,Iposnp,deltat_2)
          lbc=sqrt(N2old(2)/N2old(3))
	  lbc=1.
	  rbc=(N2old(nx)/N2old(nx-1))
	  call lcpfct(N2old,N2new,Ipos1,Iposn,
     &		      lbc,0.,0.,N2new(np),.false.,1)
	
!	ion N+

	  call velocity(Velnc,Ipos1,Iposnp,deltat_2)
          lbc=sqrt(N3old(2)/N3old(3))
	  lbc=1.
	  rbc=(N3old(nx)/N3old(nx-1))
	  call lcpfct(N3old,N3new,Ipos1,Iposn,
     &		      lbc,0.,0.,N3new(np),.false.,1)
	
!	ion N2+

	  call velocity(Velmc,Ipos1,Iposnp,deltat_2)
          lbc=sqrt(N4old(2)/N4old(3))
	  lbc=1.
	  rbc=(N4old(nx)/N4old(nx-1))
	  call lcpfct(N4old,N4new,Ipos1,Iposn,
     &		      lbc,0.,0.,N4new(np),.false.,1)
	
!	ion NO+

          lbc=sqrt(N5old(2)/N5old(3))
	  lbc=1.
	  rbc=(N5old(nx)/N5old(nx-1))
	  call lcpfct(N5old,N5new,Ipos1,Iposn,
     &		      lbc,0.,0.,N5new(np),.false.,1)
	
!	ion O2+

          lbc=sqrt(N6old(2)/N6old(3))
	  lbc=1.
	  rbc=(N6old(nx)/N6old(nx-1))
	  call lcpfct(N6old,N6new,Ipos1,Iposn,
     &		      lbc,0.,0.,N6new(np),.false.,1)
	
	  do i=1,nx

	
	    Tjr=T_0*T1new(i)
	    Tmr=T_0*Tmnew(i)
	    Tr=(Tn(i)+Tmr)/2.
	    Trh=sqrt(Tr)
	    lTr=log10(Tr)

            nukO2(i)=0.449e-9*No2(i)                            *t0
            nukO (i)=0.258e-9*No (i)                            *t0
            nukN2 (i)=5.14e-11*Nn2 (i)
     &             *(1.-.069*alog10(Tr))**2*sqrt(Tr)            *t0
            nulN2(i)=0.413e-9*Nn2(i)                            *t0
            nulO (i)=0.231e-9*No (i)                            *t0
            nulO2 (i)=2.59e-11*No2 (i)
     &             *(1.-.073*alog10(Tr))**2*sqrt(Tr)            *t0
            numN2(i)=0.434e-9*Nn2(i)                            *t0
            numO (i)=0.244e-9*No (i)                            *t0
            numO2(i)=0.427e-9*No2 (i)                           *t0

	    nu_omega=(nukN2(i)+nukO2(i)+nukO(i)
     &               +nulN2(i)+nulO2(i)+nulO(i)
     &               +numN2(i)+numO2(i)+numO(i))/omega(i)*omz
	    coef_conv=1./(1.+nu_omega**2)
	    Tperp(i)=cofterp*coef_conv*Vm_2(i)

	
	    TjO2=(Tjr*32.+Tn(i)*16.)/48.+(16.*32./48.*Tperp(i))
	    kjO2=TjO2/300.
	    TjN2=(Tjr*28.+Tn(i)*16.)/44.+(16.*28./44.*Tperp(i))
	    kjN2=TjN2/300.


	    TN2O=(Tmnew(i)*T_0*16.+Tn(i)*28.)/44.+(16.*28./44.*Tperp(i))
	    kN2O=TN2O/300.

	    TONO=(Tjr*30.+Tn(i)*16.)/46.+(16.*30./46.*Tperp(i))
	    kONO=TONO/300.

	    TlO2=(32.*Tmr+16.*Tn(i))/48.+(16.*32./48.*Tperp(i))
	    klO2=TlO2/300.

	    TN2O2=(Tmr*32.+Tn(i)*28.)/60.+(32.*28./60.*Tperp(i))

c****************************************************************
c O+ + N2 -> NO+ + N  in an O buffer: *
c**************************************
	if (TjN2 .ge. 100. .and. TjN2 .le. 6200.) then
	  ak1 = 1.248e-12 - 1.751e-13*kN2O
     &           - 5.101e-14*kN2O**2 + 1.345e-14*kN2O**3
     &           - 2.922e-16*kN2O**4
        elseif ( TjN2 .gt. 6200. .and. TjN2 .le. 22000.) then
	  ak1 = -9.626e-11 + 6.994e-12*kN2O
     &           - 2.315e-14*kN2O**2
        elseif ( TjN2 .gt. 22000.) then
	  ak1 = -9.626e-11 + 6.994e-12*kN2O
     &           - 2.315e-14*kN2O**2
        endif
c*****************************************
c O+ + N2 -> NO+ + N  in an N2 buffer:   *
c*****************************************
	if (TjN2 .ge. 100. .and. TjN2 .le. 5500.) then
	  ak1 = 1.417e-12 - 3.223e-13*kN2O - 2.362e-14*kN2O**2
     &		 + 1.247e-14*kN2O**3 - 3.030e-16*kN2O**4	
	elseif(TjN2 .gt. 5500. .and. TjN2 .le. 29000.) then
	  ak1 =-4.74e-11 + 3.74e-12*kN2O + 2.80e-14*kN2O**2
     &		 - 1.88e-16*kN2O**3
        elseif ( TjN2 .gt. 29000.) then
	  ak1 =-4.74e-11 + 3.74e-12*kN2O + 2.80e-14*kN2O**2
     &		 - 1.88e-16*kN2O**3
        endif


c****************************************
c O+ + O2 -> O2+ + O in an O buffer:    *
c****************************************
	if (TjO2 .ge. 100. .and. TjO2 .le. 6400.) then
	   ak2 = 2.836e-11 - 7.521e-12*kjO2 + 1.039e-12*kjO2**2
     &            - 4.981e-14*kjO2**3 + 9.087e-16*kjO2**4
	elseif (TjO2 .gt. 6400. .and. TjO2 .le. 22000.) then
	   ak2 =-3.42e-11 + 4.08e-12*kjO2 - 1.70e-14*kjO2**2
        elseif (TjO2 .gt. 22000.) then
           ak2 =-3.42e-11 + 4.08e-12*kjO2 - 1.70e-14*kjO2**2
        endif
c****************************************
c O+ + O2 -> O2+ + O in an N2 buffer:    *
c****************************************
	if (TjO2 .ge. 100. .and. TjO2 .le. 8400. ) then
	   ak2 = 2.763e-11 - 6.733e-12*kjO2 + 8.383e-13*kjO2**2
     &            - 3.317e-14*kjO2**3 + 4.805e-16*kjO2**4
	elseif (TjO2 .gt. 8400. .and. TjO2 .le. 31000.) then
	   ak2 =-2.57e-11 + 3.48e-12*kjO2 - 1.01e-14*kjO2**2
        elseif (TjO2 .gt. 31000.) then
           ak2 =-2.57e-11 + 3.48e-12*kjO2 - 1.01e-14*kjO2**2
	endif


c*************************
c  N2+ + O -> O+ + N2    *
c*************************
	    if (TN2O.le.1500.) then
	      ak3=1.e-11*kN2O**(-.23)
	      ak5=1.4e-10*kN2O**(-.44)
	    else
	      ak3=3.6e-12*kN2O**.41
	      ak5=5.2e-11*kN2O**.2
	    endif


c****************************************
c O+ + NO -> NO+ + O in an  O buffer:   *
c****************************************
	if (TONO .gt. 100. .and. TONO .le. 6300. ) then
	  ak4 = 5.974e-13 - 9.422e-14*kONO + 6.583e-14*kONO**2
     &           - 2.156e-15*kONO**3 + 3.957e-17*kONO**4
        elseif( TONO .gt. 6300. .and. TONO .le. 22000.) then
	  ak4 =-1.557e-11 + 1.397e-12*kONO + 2.461e-15*kONO**2
        elseif( TONO .gt. 22000.) then
          ak4 =-1.557e-11 + 1.397e-12*kONO + 2.461e-15*kONO**2
	endif
c****************************************
c O+ + NO -> NO+ + O in an  N2 buffer:   *
c****************************************
	if (TONO .gt. 100. .and. TONO .le. 8200. ) then
	  ak4 = 5.622e-13 - 6.094e-14*kONO + 5.74e-14*kONO**2
     &           - 1.399e-15*kONO**3 + 1.84e-17*kONO**4
        elseif( TONO .gt. 8200. .and. TONO .le. 30000.) then
	  ak4 =-2.22e-11 + 1.64e-12*kONO - 6.7e-17*kONO**2
        elseif( TONO .gt. 30000.) then
          ak4 =-2.22e-11 + 1.64e-12*kONO - 6.7e-17*kONO**2
	endif
      	
c****************************************
c  O+ + O2 -> O2+ + O in an N2 buffer:  *
c****************************************
	    if( TlO2 .ge. 100. .and. TlO2 .le. 8400. ) then
	      akl2 = 2.763e-11 - 6.733e-12*klO2 + 8.383e-13*klO2**2
     &               - 3.317e-14*klO2**3 + 4.805e-16*klO2**4
	    elseif( TlO2 .gt. 8400. .and. TlO2 .le. 31000. ) then
	      akl2 =-2.57e-11 + 3.48e-12*klO2 - 1.01e-14*klO2**2
            elseif( TlO2 .gt. 31000.) then
              akl2 =-2.57e-11 + 3.48e-12*klO2 - 1.01e-14*klO2**2
	    endif

      	
    	
      	    y0  (1  )=N1old(i)
      	    y1  (1  )=(N1new(i)-N1old(i))/deltat_2+Po(i)/N_0*t0
      	    chim(1,1)=-(ak1*Nn2(i)+ak2*No2(i)+ak4*N_0*Nnonew(i)
     &                +2.5e-11*sqrt(Tn(i)+T1new(i)*T_0/16.
     &                +1.2e-8*(U1new(i)*Cj0)**2)*Nh(i))		*t0
     &                -3.*Cji*U1new(i)*alt_geo_1(i)
      	    chim(1,2)=2.2e-11*sqrt(T_0*T2new(i)+Tn(i)/16.
     &    	     +1.2e-8*(U2new(i)*Ci0)**2)*No(i)*t0
      	    chim(1,3)=(3.7e-11*No2(i)+5.e-13*No(i))*t0
      	    chim(1,4)=ak3*No(i)*t0
      	    chim(1,5)=0.d0
      	    chim(1,6)=0.d0
      	
      	    y0  (2  )=N2old(i)
      	    y1  (2  )=(N2new(i)-N2old(i))/deltat_2+Ph(i)/N_0*t0
      	    chim(2,1)=2.5d-11*sqrt(Tn(i)+T1new(i)*T_0/16.+1.2e-8*
     &      	     (U1new(i)*Cj0)**2)*Nh(i)*t0
      	    chim(2,2)=-2.2d-11*sqrt(T_0*T2new(i)+Tn(i)/16.
     &		      +1.2d-8*(U2new(i)*Ci0)**2)*No(i)		*t0
     &		      -3.*U2new(i)*alt_geo_1(i)
      	    chim(2,3)=3.6d-12*Nh(i)*t0
      	    chim(2,4)=0.d0
      	    chim(2,5)=0.d0
      	    chim(2,6)=0.d0
      	
      	    y0  (3  )=N3old(i)
      	    y1  (3  )=(N3new(i)-N3old(i))/deltat_2+0.21*Pn2(i)*t0/N_0
      	    chim(3,1)=0.d0
      	    chim(3,2)=0.d0
      	    chim(3,3)=-(2.6e-10*No2(i)+3.1e-10*No2(i)+3.7e-11*No2(i)
     &                +2.e-11*N_0*Nnonew(i)+3.6e-12*Nh(i)
     &                +5.e-13*No(i))*t0
     &		      -3*Cni*U3new(i)*alt_geo_1(i)
      	    chim(3,4)=0.d0
      	    chim(3,5)=0.d0
      	    chim(3,6)=0.d0
      	
      	    y0  (4  )=N4old(i)
      	    y1  (4  )=(N4new(i)-N4old(i))/deltat_2+Pn2(i)*0.79*t0/N_0      	
	    chim(4,1)=0.d0
      	    chim(4,2)=0.d0
      	    chim(4,3)=0.d0
      	    chim(4,4)=-(5.e-11*(300./TN2O2)*No2(i)+ak3*No(i)+ak5*No(i)
     &                +3.3e-10*N_0*Nnonew(i)
     &                +1.8e-7*Ter_1(i)**.39*Nenew(i)*N_0)*t0
     &		      -3.*Czi*Umnew(i)*alt_geo_1(i)
      	    chim(4,5)=0.d0
      	    chim(4,6)=0.d0
      	
      	    y0  (5  )=N5old(i)
      	    y1  (5  )=(N5new(i)-N5old(i))/deltat_2+6.e-7*Nnonew(i)*t0      	
      	    chim(5,1)=(ak4*N_0*Nnonew(i)+ak1*Nn2(i))*t0
      	    chim(5,2)=0.d0
      	    chim(5,3)=(No2(i)*2.6e-10+2.e-11*N_0*Nnonew(i))*t0
      	    chim(5,4)=ak5*No(i)*t0
      	    chim(5,5)=-4.2e-7*Ter_1(i)**.85*Nenew(i)*N_0*t0
     &		      -3.*Czi*Umnew(i)*alt_geo_1(i)
      	    chim(5,6)=(5.e-16*Nn2(i)+1.2e-10*Nn(i)
     &                 +4.5e-10*N_0*Nnonew(i))*t0
      	
      	    y0  (6  )=N6old(i)
      	    y1  (6  )=(N6new(i)-N6old(i))/deltat_2+Po2(i)/N_0*t0
      	    chim(6,1)=akl2*No2(i)*t0
      	    chim(6,2)=0.d0
      	    chim(6,3)=3.1e-10*No2(i)*t0
      	    chim(6,4)=5.e-11*(300./TN2O2)*No2(i)*t0
      	    chim(6,5)=0.d0
      	    chim(6,6)=-(5.e-16*Nn2(i)+1.2e-10*Nn(i)
     &		      +4.5e-10*N_0*Nnonew(i)
     &                +1.6e-7*Ter_1(i)**.55*Nenew(i)*N_0)*t0
     &		      -3.*Czi*Umnew(i)*alt_geo_1(i)

	    call solve_ODE(y0,y1,chim,mat,deltat_2)
	
	    N1new(i)=max(y0(1),r_min)
	    N2new(i)=max(y0(2),r_min)
	    N3new(i)=max(y0(3),r_min)
	    N4new(i)=max(y0(4),r_min)
	    N5new(i)=max(y0(5),r_min)
	    N6new(i)=max(y0(6),r_min)

	  enddo
	
          N1new(np)=min(1.,N1new(nx-1)/N1new(nx-2))*N1new(nx)
          N2new(np)=min(1.,N2new(nx-1)/N2new(nx-2))*N2new(nx)
          N3new(np)=min(1.,N3new(nx-1)/N3new(nx-2))*N3new(nx)
          N4new(np)=min(1.,N4new(nx-1)/N4new(nx-2))*N4new(nx)
          N5new(np)=min(1.,N5new(nx-1)/N5new(nx-2))*N5new(nx)
          N6new(np)=min(1.,N6new(nx-1)/N6new(nx-2))*N6new(nx)

          if (isnant(N1new,nx)) then
            print*,'probleme lors du calcul de N1new dans la boucle 1'
            goto 246
          endif
          if (isnant(N2new,nx)) then
            print*,'probleme lors du calcul de N2new dans la boucle 1'
            goto 246
          endif
          if (isnant(N3new,nx)) then
            print*,'probleme lors du calcul de N3new dans la boucle 1'
            goto 246
          endif
          if (isnant(N4new,nx)) then
            print*,'probleme lors du calcul de N4new dans la boucle 1'
            goto 246
          endif
          if (isnant(N5new,nx)) then
            print*,'probleme lors du calcul de N5new dans la boucle 1'
            goto 246
          endif
          if (isnant(N6new,nx)) then
            print*,'probleme lors du calcul de N6new dans la boucle 1'
            goto 246
          endif
