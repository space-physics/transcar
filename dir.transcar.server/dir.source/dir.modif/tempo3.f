




        Un(i)=-w(1)*sin(dipangle*deg2rad)	!	Un est la composante parallele du vent neutre
        Vn(i)=w(2)  				!	Vn est la composante est du vent neutre
        Wn(i)=w(1)*cos(dipangle*deg2rad)  	!       Wn est la composante nord du vent neutre
        Un (iz) = Un (iz-1)
        Vn (iz) = Vn (iz-1)
        Wn (iz) = Wn (iz-1)

  	  B=B0*(Re/(Re+z(i)))**3
  	  omega(i)=1.6e-19*B*1.e-4/1.667e-27*t0

    !C       Vm est la composante est de la derive magnetique

  	  Vm(i)=-Enord/Bmag*1000.-Vn(i)

    !C       Wm est la composante nord de la derive magnetique

  	  Wm(i)=Eest/Bmag*1000.-Wn(i)


  	Vm_2(i)=Vm(i)**2+Wm(i)**2

    !c	J0 est en ÊA/m-2
  	  JJ(i)=J0/1.6e-9*((800.+Re)/(z(i)+Re))**3

              nu_omega=1.e-9*Nn2(i)*t0/omega(i)*30.5
            coef_cour=1./(1.+nu_omega**2)

            JJ(i)=JJ(i)*coef_cour


  	enddo
