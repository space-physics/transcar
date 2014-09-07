    Dcei=0.d0

    do i=1,nb_ion

      y1(i)=-Cij(iref,i)*G(iz)-nuij(i,indi_e)*Jpara(iz)/dens_i(indi_e))/Co(i)
      y1(i+nb_ion) = nuij(i,indi_e)*Dco(i,indi_e)*dens_i(i)*temperature(i)/dens_i(indi_e)*Jpara(iz)/Co(i)

      colli(i,i)=0.d0
      colli(i,i+nb_ion)=0.d0

      do j=1,nb_neutre
        y1(i)=y1(i)+nuin(i,j)*Un(iz)/Co(i)

        colli(i,i)=colli(i,i)-nuin(i,j)
      end do

      do j=1,nb_ion
        if (j/=i) then
          colli(i,j       ) = (nuij(i,j)+nuij(i,indi_e)*dens_i(j)/dens_i(indi_e))*Cij(j,i)
          colli(i,i       ) = colli(i,i)-nuij(i,j)
          colli(i,i+nb_ion) = colli(i,i+nb_ion)+nuij(i,j)/Tij(i,j)*Ac(i,j)*dens_i_1(i)
          colli(i,j+nb_ion) = -Ac(j,i)*dens_i_1(j)*Cij(j,i)
        end if
      end do

      colli(i,i            ) = colli(i,i)-nuij(i,indi_e)*(1.-dens_i(i)/dens_i(indi_e))
      colli(i,i+nb_ion     ) = colli(i,i+nb_ion)+nuij(i,indi_e)/Tij(i,indi_e)*Ac(i,indi_e)*dens_i_1(i)
      colli(i,indi_e+nb_ion) = -Ac(indi_e,i)*dens_i_1(indi_e)*Cij(indi_e,i)

      colli(i+nb_ion,i)=0.d0
      colli(i+nb_ion,i+nb_ion)=-0.8*nuij(i,i)

      do j=1,nb_ion
        if (i/=j) then
          colli(i+nb_ion,i       ) = colli(i+nb_ion,i)+dens_i(i)*temperature(i)	&
          					*(nuij(i,j)*Dco(i,j)+nuij(i,indi_e)*Dco(i,indi_e)*(1.-dens_i(i)/dens_i(e)))
          colli(i+nb_ion,i+nb_ion)=colli(i+nb_ion,i+nb_ion)-nuij(i,j)*(Dc1(i,j)+Dco(i,j)*temperature(i)/Tij(i,j))

          colli(i+nb_ion,j       ) = -dens_i(i)*temperature(i)*Cij(j,i)		&
          				       *(nuij(i,j)*Dco(i,j)+nuij(i,indi_e)*Dco(i,indi_e)*dens_i(j)/dens_i(e))
          colli(i+nb_ion,j+nb_ion) = nuij(i,j)*dens_i(i)*dens_i_1(j)*Cij(j,i)*(Dc4(i,j)+Dco(j,i)*temperature(i)/Tij(i,j))
        end if
      end do

      do j=1,nb_neutre
        colli(i+nb_ion,i+nb_ion) = colli(i+nb_ion,i+nb_ion)-nuin(i,j)*Dm1(i,j)
      end do

      colli(i+nb_ion,i+nb_ion) = colli(i+nb_ion,i+nb_ion)	&
      				-nuij(i,indi_e)*(Dc1(i,indi_e)+Dco(i,indi_e)*temperature(i)/Tij(i,indi_e))
      				
      colli(i+nb_ion,indi_e+nb_ion) = nuij(i,indi_e)*dens_i(i)*dens_i_1(indi_e)*Cij(indi_e,i)	&
				      *(Dc4(i,indi_e)+Dco(indi_e,i)*temperature(i)/Tij(i,indi_e))

      Dcei=Dcei+Dco(indi_e,i)*nuij(indi_e,i)

    end do

    Dcei=Dcei*temperature(indi_e)

    y1(indi_e+nb_ion) = -Dcei*Jpara(iz)/Co(indi_e)

    colli(indi_e+nb_ion,indi_e+nb_ion)=-0.8*nuij(indi_e,indi_e)

    do i=1,nb_ion
      colli(indi_e+nb_ion,i+nb_ion     ) = nuij(indi_e,i)*dens_i(indi_e)*dens_i_1(i)*Cij(i,indi_e)	&
                                           *(Dc4(indi_e,i)+Dco(i,indi_e)*temperature(indi_e)/Tij(indi_e,i))
      colli(indi_e+nb_ion,i            ) = (Dcei*dens_i(i)-nuij(indi_e,i)*Dco(indi_e,i)*dens_i(indi_e)*temperature(indi_e))*Cij(indi_e,i)

      colli(indi_e+nb_ion,indi_e+nb_ion)=colli(indi_e+nb_ion,indi_e+nb_ion)			&
					  -nuij(indi_e,i)*(Dc1(indi_e,i)+Dco(indi_e,i)*temperature(indi_e)/Tij(indi_e,i))
    end do

    do j=1,nb_neutre
      colli(indi_e+nb_ion,indi_e+nb_ion)=colli(indi_e+nb_ion,indi_e+nb_ion)-nuin(indi_e,j)*Dm1(indi_e,j)
    end do

