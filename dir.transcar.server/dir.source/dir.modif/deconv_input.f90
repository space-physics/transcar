do iz=1,nb_alt
  ipos=(iz+1)*ncol
  alt(iz)      =buffer(ipos+ipos_i%z)

  N1new(iz)    =buffer(ipos+ipos_i%n1)/N_o
  N2new(iz)    =buffer(ipos+ipos_i%n2)/N_o
  N3new(iz)    =buffer(ipos+ipos_i%n3)/N_o
  N4new(iz)    =buffer(ipos+ipos_i%n4)/N_o
  N5new(iz)    =buffer(ipos+ipos_i%n5)/N_o
  N6new(iz)    =buffer(ipos+ipos_i%n6)/N_o
  U1new(iz)    =buffer(ipos+ipos_i%u1)/Co(1)
  U2new(iz)    =buffer(ipos+ipos_i%u2)/Co(2)
  U3new(iz)    =buffer(ipos+ipos_i%u3)/Co(3)
  U4new(iz)    =buffer(ipos+ipos_i%u4)/Co(4)
  U5new(iz)    =buffer(ipos+ipos_i%u5)/Co(5)
  U6new(iz)    =buffer(ipos+ipos_i%u6)/Co(6)
  Uenew(iz)    =buffer(ipos+ipos_i%ue)/Co(indi_e)
  T1pnew(iz)   =buffer(ipos+ipos_i%t1p)/T_o
  T1tnew(iz)   =buffer(ipos+ipos_i%t1t)/T_o
  T2pnew(iz)   =buffer(ipos+ipos_i%t2p)/T_o
  T2tnew(iz)   =buffer(ipos+ipos_i%t2t)/T_o
  T3pnew(iz)   =buffer(ipos+ipos_i%t3p)/T_o
  T3tnew(iz)   =buffer(ipos+ipos_i%t3t)/T_o
  T4pnew(iz)   =buffer(ipos+ipos_i%t4p)/T_o
  T4tnew(iz)   =buffer(ipos+ipos_i%t4t)/T_o
  T5pnew(iz)   =buffer(ipos+ipos_i%t5p)/T_o
  T5tnew(iz)   =buffer(ipos+ipos_i%t5t)/T_o
  T6pnew(iz)   =buffer(ipos+ipos_i%t6p)/T_o
  T6tnew(iz)   =buffer(ipos+ipos_i%t6t)/T_o
  Tepnew(iz)   =buffer(ipos+ipos_i%tep)/T_o
  Tetnew(iz)   =buffer(ipos+ipos_i%tet)/T_o
  q1new(iz)    =buffer(ipos+ipos_i%q1)/q_o(1)
  q2new(iz)    =buffer(ipos+ipos_i%q2)/q_o(2)
  q3new(iz)    =buffer(ipos+ipos_i%q3)/q_o(3)
  q4new(iz)    =buffer(ipos+ipos_i%q4)/q_o(4)
  q5new(iz)    =buffer(ipos+ipos_i%q5)/q_o(5)
  q6new(iz)    =buffer(ipos+ipos_i%q6)/q_o(6)
  qenew(iz)    =buffer(ipos+ipos_i%qe)/q_o(indi_e)
  Nnonew(iz)   =buffer(ipos+ipos_i%nno)/N_o
  Unonew(iz)   =buffer(ipos+ipos_i%uno)/Co(indi_NO)
  Po(iz)       =buffer(ipos+ipos_i%po)
  Ph(iz)       =buffer(ipos+ipos_i%ph)
  Pn(iz)       =buffer(ipos+ipos_i%pn)
  Pn2(iz)      =buffer(ipos+ipos_i%pn2)
  Po2(iz)      =buffer(ipos+ipos_i%po2)
  Heat(iz)     =buffer(ipos+ipos_i%heat)

  Nenew(iz)=N1new(iz)+N2new(iz)+N3new(iz)+N4new(iz)+N5new(iz)+N6new(iz)

  T1new(iz)=(T1pnew(iz)+2.*T1tnew(iz))/3.
  T2new(iz)=(T2pnew(iz)+2.*T2tnew(iz))/3.
  T3new(iz)=(T3pnew(iz)+2.*T3tnew(iz))/3.
  T4new(iz)=(T4pnew(iz)+2.*T4tnew(iz))/3.
  T5new(iz)=(T5pnew(iz)+2.*T4tnew(iz))/3.
  T6new(iz)=(T6pnew(iz)+2.*T4tnew(iz))/3.
  Tenew(iz)=(Tepnew(iz)+2.*Tetnew(iz))/3.

enddo
