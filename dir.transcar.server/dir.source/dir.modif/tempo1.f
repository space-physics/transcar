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

