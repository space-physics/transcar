
  Lenrot=(5.35e-39*Nn2(iz)+1.1e-38*No2(iz))*dTen/Teh

  ff=1.06e4+7.51e3*tanh(1.1e-3*temperature(indi_e)-1.98)
  gg=1245.+2.261*temperature(indi_e)-2.056e-4*temperature(indi_e)**2
  LeN2vib=4.78e-39*Nn2(iz)*exp(ff*(5.e-4-Te_1)*(exp(gg*dTen*Te_1*Tn_1)-1.)

  hh=3300.-839.*sin(1.91e-4*temperature(indi_e)-.5157))
  LeO2vib=8.3e-38*No2(iz)*exp(hh*(1.429e-3-Te_1))*(exp(2770.*dTen*Te_1*Tn_1)-1.)

  dd=-2.34e4*Te_1+1.55e-4*temperature(indi_e)-6.49e-9*temperature(indi_e)**2
  LeOexc=4.1e-34*No(iz)*exp(dd)*(exp(22713.*dTen*Te_1*Tn_1)-1.)
  !LeOexc=2.5e-37*No(iz)*exp(dd*(3.333e-4-Te_1))*(exp(22713.*dTen*Te_1*Tn_1)-1.)


  LeOfin=0.
  Zz=5.
  do j=1,3
    select case j
      case (1)
        Ex=exp(-Eq(1)*Te_1)
        Dx=exp(-Eq(1)*Tn_1)
        Zz=Zz+3.*Dx
      case (2)
        Ex=exp(-Eq(2)*Te_1)
        Dx=exp(-Eq(2)*Tn_1)
        Zz=Zz+Dx
      case default
        Ex=exp(-Eq(3)*Te_1-Eq(1)*Tn_1)
        Dx=exp(-Eq(2)*Tn_1)
    end select
    Ff=epsq(j)*(Dx-Ex)-5.91e-9*dTen*((1.+Bq(j))*Dx+Ex*(Eq(j)*Te_1+1.+Bq(j)))
    LeOfin=LeOfin+Aq(j)*Cq(j)*Ff*(temperature(indi_e))**(Bq(j)-.5)
  enddo
  LeOfin=1.38e-30*No(iz)/Zz*LeOfin

  Len=Len_o*(LeOfin+Lenrot+LeN2vib+LeO2vib+LeOexc+LenCO2vib)/dTen
  Len_o=deux_tiers*to*no/kb
  Qen=Heat(iz)+Len*Tn(iz)

  do i=1,indi_e
    ther(i,i)=0.d0
    ther(i,i+indi_e)=0.d0
    do j=1,indi_e
      if (j/=i) then
        ther(i,j) = nuij(i,j)*Dp1(i,j)
        ther(i,j+indi_e)=nuij(i,j)*Dp3(i,j)
        ther(i,i)=ther(i,i)-nuij(i,j)*Dp2(i,j)
        ther(i,i+indi_e)=nuij(i,j)*Dp3(j,i)

        ther(i+indi_e,j) = nuij(i,j)*Dt3(i,j)
        ther(i+indi_e,j+indi_e)=nuij(i,j)*Dt1(i,j)
        ther(i+indi_e,i)=nuij(i,j)*Dt3(j,i)
        ther(i+indi_e,i+indi_e)=ther(i+indi_e,i+indi_e)-nuij(i,j)*Dt2(i,j)
      else
        ther(i,i)=ther(i,i)-Dp4(i)*nui(i,i)
        ther(i,i+indi_e)=ther(i,i+indi_e)+Dp4*nuij(i,i)

        ther(i+indi_e,i)=ther(i+indi_e,i)+Dt4(i)*nui(i,i)
        ther(i+indi_e,i+indi_e)=ther(i+indi_e,i+indi_e)-Dt4*nuij(i,i)
      endif
    enddo
  enddo
  ther(indi_e,indi_e)=ther(indi_e,indi_e)-Len
  ther(indi_e+indi_e,indi_e+indi_e)=ther(indi_e+indi_e,indi_e+indi_e)-2.*Len


  y0     (1        ) = T1pold(iz)
  y1     (1        ) = (T1pnew(iz)-T1pold(iz))/deltat_2
  y0     (1+indi_e ) = T1told(iz)
  y1     (1+indi_e ) = (T1tnew(iz)-T1told(iz))/deltat_2

  y0     (2        ) = T2pold(iz)
  y1     (2        ) = (T2pnew(iz)-T2pold(iz))/deltat_2
  y0     (2+indi_e ) = T2told(iz)
  y1     (2+indi_e ) = (T2tnew(iz)-T2told(iz))/deltat_2

  y0     (3        ) = T3pold(iz)
  y1     (3        ) = (T3pnew(iz)-T3pold(iz))/deltat_2
  y0     (3+indi_e ) = T3told(iz)
  y1     (3+indi_e ) = (T3tnew(iz)-T3told(iz))/deltat_2

  y0     (4        ) = T4pold(iz)
  y1     (4        ) = (T4pnew(iz)-T4pold(iz))/deltat_2
  y0     (4+indi_e ) = T4told(iz)
  y1     (4+indi_e ) = (T4tnew(iz)-T4told(iz))/deltat_2

  y0     (5        ) = T5pold(iz)
  y1     (5        ) = (T5pnew(iz)-T5pold(iz))/deltat_2
  y0     (5+indi_e ) = T5told(iz)
  y1     (5+indi_e ) = (T5tnew(iz)-T5told(iz))/deltat_2

  y0     (6        ) = T6pold(iz)
  y1     (6        ) = (T6pnew(iz)-T6pold(iz))/deltat_2
  y0     (6+indi_e ) = T6told(iz)
  y1     (6+indi_e ) = (T6tnew(iz)-T6told(iz))/deltat_2

  y0(indi_e        ) = Tepold(iz)
  y1(indi_e        ) = (Tepnew(iz)-Tepold(iz))/deltat_2+Qen
  y0(indi_e+indi_e ) = Tetold(iz)
  y1(indi_e+indi_e ) = (Tetnew(iz)-Tetold(iz))/deltat_2+2.*Qen

  call solve_ODE(2*indi_e,y0,y1,ther,mat,deltat_2)

  T1pnew(iz) = y0(1)
  T1tnew(iz) = y0(1+indi_e)
  T2pnew(iz) = y0(2)
  T2tnew(iz) = y0(2+indi_e)
  T3pnew(iz) = y0(3)
  T3tnew(iz) = y0(3+indi_e)
  T4pnew(iz) = y0(4)
  T4tnew(iz) = y0(4+indi_e)
  T5pnew(iz) = y0(5)
  T5tnew(iz) = y0(5+indi_e)
  T6pnew(iz) = y0(6)
  T6tnew(iz) = y0(6+indi_e)
  Tepnew(iz) = y0(indi_e)
  Tetnew(iz) = y0(indi_+indi_e)

  T1new(iz)=r1_3*T1pnew(iz)+r2_3*T1tnew(iz)
  T2new(iz)=r1_3*T2pnew(iz)+r2_3*T2tnew(iz)
  T3new(iz)=r1_3*T3pnew(iz)+r2_3*T3tnew(iz)
  T4new(iz)=r1_3*T4pnew(iz)+r2_3*T4tnew(iz)
  T5new(iz)=r1_3*T5pnew(iz)+r2_3*T4tnew(iz)
  T6new(iz)=r1_3*T6pnew(iz)+r2_3*T4tnew(iz)
  Tenew(iz)=r1_3*Tepnew(iz)+r2_3*Tetnew(iz)
