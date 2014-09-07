c
c--------------------------------------------------------------------
c
	SUBROUTINE LOW_PROTON(Eflux,Emoy_keV,nalt,ZPHT_,
     .		XNN2_,XNO2_,XNO_,QTI_,QIA_)
c
c  Computation of the electron and ion production rates
c  induced by a proton beam of a given mean energy 
c  and a given energy flux.
c
c  * INPUTS:
c
c  Eflux: energy flux of the incident proton flux, in mW m-2.
c  Emoy_keV  : mean energy of the incident proton flux, in keV.
c  	  Restriction : 2 < Emoy_keV < 40 keV
c	  If out of the range, set to closest value (2 keV or 40 keV)
c  	  Emoy_keV = 2*E0, with E0, the characteristic energy of 
c	  the incident flux(having assumed the distribution in 
c 	  energy to be Maxwellian)
c
c  nalt : number of altitude level
c  ZPHT_: altitude grid (in km) [decreasing altitude grid !]
c  XNN2_, XNO2_, XNO_: number density of N2, O2, and O, respectively (in m-3)
c
c  * OUTPUTS:
c
c  QIA_(4,nalt) : ion production rate (m-3 s-1)
c               j=1,2,3,4 -> N2+,O2+,O+,N+ 
c  QTI_ : electron production rate (m-3 s-1)
c
c  Notes: * if the altitude grid has only few levels, set nk to an integer
c  larger than 1 (-> addition of nk points between each two successive
c  altitude levels)
c         * otherwise, set nk to 1.
c	  * Number maximum of the altitude grid on which the computation is
c  performed: nbrz.
c
c  You can find the detailled explanation of the parameterization in
c  Galand et al., Ionization by energetic protons in Thermosphere-
c  Ionosphere Electrodynamics General Circulation Model,
c  Journal of Geophysical Research, 104, 27973-27989, 1999.
c--------------------------------------------------------------------

c Parameters
	integer nk,nbrz
	real*4 tab_EL(6),tab_KL(6),tab_E0(6)
	real*4 EminV,Mean_mass
	real*4 P_W1(7),P_W2(7),P_W3(7),P_W4(7),P_W5(7),P_W6(7)
	real*4 tab_shape(9)
	real*4 P_Kd1(3),P_Kd2(3),P_Kd3(3),P_Kd4(3)
	real*4 E_Kd(4)
	real*4 Mass(3)

	parameter(nk=1,nbrz=201)
        DATA tab_EL/0.87,   0.84,  0.883,   0.9,    0.92,   0.92/
	DATA tab_KL/4.5e-18,8.1e-18,7.0e-18,7.0e-18,7.0e-18,8.e-18/
	DATA tab_E0/13.5,   6.5,   3.5,    2.5,    1.5,    0/
	DATA EminV/100/
	DATA Mean_mass/4.23e-23/
c 	N2/H+
	DATA P_W1/37.90,  2.95e-1,2.03e-1, -1.57e-2,3.44e-4, 0,    0/
c 	N2/H
	DATA P_W2/156.68,-63.86,17.33,-2.32,1.64e-1,-5.83e-3,8.21e-5/
c 	O2/H+
	DATA P_W3/8.01e-1,2.09,  -2.91e-2,  1.58e-4,  0,     0,    0/
c 	O2/H
	DATA P_W4/14.26,  2.73,  -8.55e-2,  2.93e-3,  0,     0,    0/
c 	O/H+ ou O/H
	DATA P_W5/28.81, 32.69,  -1.83,     4.02e-2,  0,     0,    0/
c 	electrons
	DATA P_W6/18.94, 1.32,   -8.25e-2, 3.80e-3,-7.06e-5, 0    ,0/
	DATA tab_shape/1,0,0, 0,0.1,1, 0,0.7,1/
	DATA P_Kd1/3.85e-2,   3.03e-2,   -4.74e-4/
	DATA P_Kd2/4.56e-2,   1.27e-2,   -1.15e-4/
	DATA P_Kd3/8.90e-2,   3.95e-2,   -4.21e-4/
	DATA P_Kd4/3.15e-1,   2.19e-2,   -3.92e-4/
	DATA E_Kd/30, 50, 40, 30/
c Mass in g
c	DATA Mass/2*2.34*1.e-23,2*2.67*1.e-23,2.67*1.e-23/
	DATA Mass/4.68e-23,5.34e-23,2.67e-23/


c Input data
	real Eflux,Emoy_keV

c Complete set of input/output data (MG) 
	integer nalt,nz
	real*4 ZPHT_(nalt),XNN2_(nalt),XNO2_(nalt),XNO_(nalt)
	real*4 QTI_(nalt),QIA_(4,nalt)
	
c Local variables
	integer i,j,K,ik,iz
	real*4 ZPHT(nbrz),XNN2(nbrz),XNO2(nbrz),XNO(nbrz)
	real*4 QTI(nbrz),QIA(4,nbrz)
	real*4 RHO(nbrz)
	real*4 Emoy_eV,EL,KL,RANGE,F_E,R_norm,LAMBDA,ETA,W(6),shape(3)
	real*4 Kd(4),E,QIA_TOT

c INITIALISATION
c --------------
c 	Attention ! Emoy_keV est l'energie moyenne, egale a deux fois l'
c 	energie caracteristique de la maxwellienne.
c 	Passage en eV
	Emoy_eV=Emoy_keV*1.e3
	if (Emoy_eV.lt.2000.) then 
		Emoy_eV=2000.
	endif
	if (Emoy_eV.gt.40000.) then 
		Emoy_eV=40000.
	endif

c Increase number of altitude levels
	nz=nk*(nalt-1)+nalt
c	write(*,*)nz
	if (nz.gt.nbrz) then
		write(*,*)'Increase nbrz'
c		write(*,*)nz
		stop
	endif
	if (nk.ge.1) then
	  do iz=2,nalt
		ZPHT((iz-2)*(nk+1)+1)=ZPHT_(iz-1)
		XNN2((iz-2)*(nk+1)+1)=XNN2_(iz-1)
		XNO2((iz-2)*(nk+1)+1)=XNO2_(iz-1)
		XNO((iz-2)*(nk+1)+1)=XNO_(iz-1)
		do ik=1,nk
	ZPHT((iz-2)*(nk+1)+ik+1)=exp((log(ZPHT_(iz))-log(ZPHT_(iz-1)))
     .		/(1.+nk)*ik+log(ZPHT_(iz-1)))
	XNN2((iz-2)*(nk+1)+ik+1)=exp((log(XNN2_(iz))-log(XNN2_(iz-1)))
     .		/(ZPHT_(iz)-ZPHT_(iz-1))
     .		*(ZPHT((iz-2)*(nk+1)+ik+1)-ZPHT_(iz-1))
     .		+log(XNN2_(iz-1)))
	XNO2((iz-2)*(nk+1)+ik+1)=exp((log(XNO2_(iz))-log(XNO2_(iz-1)))
     .		/(ZPHT_(iz)-ZPHT_(iz-1))
     .		*(ZPHT((iz-2)*(nk+1)+ik+1)-ZPHT_(iz-1))
     .		+log(XNO2_(iz-1)))
	XNO((iz-2)*(nk+1)+ik+1)=exp((log(XNO_(iz))-log(XNO_(iz-1)))
     .		/(ZPHT_(iz)-ZPHT_(iz-1))
     .		*(ZPHT((iz-2)*(nk+1)+ik+1)-ZPHT_(iz-1))
     .		+log(XNO_(iz-1)))
		enddo
	  enddo
	  ZPHT(nz)=ZPHT_(nalt)
c 	  Passages des m-3 aux cm-3
	  XNN2(nz)=XNN2_(nalt)*1.e-06
	  XNO2(nz)=XNO2_(nalt)*1.e-06
	  XNO(nz)=XNO_(nalt)*1.e-06

	else
	  do iz=1,nz
		ZPHT(iz)=ZPHT_(iz)
		XNN2(iz)=XNN2_(iz)*1.e-06
		XNO2(iz)=XNO2_(iz)*1.e-06
		XNO(iz)=XNO_(iz)*1.e-06
	  enddo
	endif

c Mass density (g cm-3)
	do iz=1,nz
	     RHO(iz)=Mass(1)*XNN2(iz)+Mass(2)*XNO2(iz)+Mass(3)*XNO(iz)
	enddo


c CALCULATION
c -----------
c RANGE (g.cm-2). Attention ! Un seul range global pour une masse 
c 	moyenne... 78% N2, 11% O2, 11% O
	do i=6,1,-1
	  if (Emoy_keV.gt.2.*tab_E0(i)) then
		EL=tab_EL(i)
		KL=tab_KL(i)
	  endif
	enddo
c 	Equation 10
	RANGE = Mean_mass/KL/(1-EL)*(Emoy_eV**(1-EL)-EminV**(1-EL)) 
c 	Apres equation 11
	F_E   = Emoy_eV**(1-EL)/KL/(1-EL)/RANGE*Mean_mass

c  	Energy loss per electron or ionization type.
c 	Valeurs du tableau 2
	do i=1,6
	      W(i)=0
	enddo
	do i=1,7
	      W(1)=W(1)+P_W1(i)*(Emoy_keV/2)**(i-1)  ! N2/H+
	      W(2)=W(2)+P_W2(i)*(Emoy_keV/2)**(i-1)  ! N2/H
	      W(3)=W(3)+P_W3(i)*(Emoy_keV/2)**(i-1)  ! O2/H+
	      W(4)=W(4)+P_W4(i)*(Emoy_keV/2)**(i-1)  ! O2/H
	      W(5)=W(5)+P_W5(i)*(Emoy_keV/2)**(i-1)  ! O/H+ ou O/H
	      W(6)=W(6)+P_W6(i)*(Emoy_keV/2)**(i-1)  ! electrons
	enddo

c  	Dissociation ratio
c 	Tableau 3
	do isp=1,4
		Kd(isp)=0.
                if (Emoy_keV.gt.E_Kd(isp)) then
                        E=E_Kd(isp)
		else
			E=Emoy_keV
		endif
		if (isp.eq.1) then	!N2+H+ --> N + N+ +H+ +e-
			do i=1,3
			Kd(isp)=Kd(isp)+P_Kd1(i)*E**(i-1)
			enddo

		elseif (isp.eq.2) then	!N2+H --> N + N+ +H +e-
			do i=1,3
			Kd(isp)=Kd(isp)+P_Kd2(i)*E**(i-1)
			enddo

		elseif (isp.eq.3) then	!O2+H+ --> O + O+ +H+ +e-
			do i=1,3
			Kd(isp)=Kd(isp)+P_Kd3(i)*E**(i-1)
			enddo

		elseif (isp.eq.4) then	!O2+H+ --> O + O+ +H +e-
			do i=1,3
			Kd(isp)=Kd(isp)+P_Kd4(i)*E**(i-1)
			enddo
		endif
	enddo

c-----------------------------------------------------------------------
c Decreasing altitude grid
	DO iz=1,nz
c R_norm : Normalized atmospheric scattering depth (g.cm-2)
c Decreasing altitudes
	  if (iz.eq.1) then
	     R_norm=0
	  else   
c Decreasing altitudes
 	     R_norm=R_norm
     .	     +(RHO(iz-1)+RHO(iz))/2*(ZPHT(iz-1)-ZPHT(iz))*1.e5/RANGE
	  endif

c LAMBDA : Normalized energy deposition function
 	  if (abs(1-EL)/(1-EL)*(F_E-R_norm).gt.0) then
    		LAMBDA=KL*RANGE/Mean_mass/(Emoy_eV-EminV)
     .	    *(((1-EL)*KL*RANGE/Mean_mass*(F_E-R_norm))**(EL/(1-EL)))
	  else
	    	LAMBDA=0
          endif

c ETA : Energy deposition rate (eV.cm-3.s-1)
 	  ETA=Eflux*1.e-7/(1.602e-19)/RANGE*RHO(iz)*LAMBDA

c QTI : Electron production rate (cm-3.s-1)
	  QTI(iz)=ETA/W(6)

c QIA(j,iz) : Ion production rate (cm-3.s-1)
c j=1,2,3,4 -> N2+,O2+,O+,N+ 
	  shape(1)=tab_shape(1)*XNN2(iz)/(tab_shape(1)*XNN2(iz)
     .		+tab_shape(2)*XNO2(iz)+tab_shape(3)*XNO(iz))
	  shape(2)=tab_shape(5)*XNO2(iz)/(tab_shape(4)*XNN2(iz)
     .		+tab_shape(5)*XNO2(iz)+tab_shape(6)*XNO(iz))
	  shape(3)=tab_shape(9)*XNO(iz)/(tab_shape(7)*XNN2(iz)
     .		+tab_shape(8)*XNO2(iz)+tab_shape(9)*XNO(iz))

	  QIA(1,iz)=shape(1)*(ETA/W(1)/(1+Kd(1))+ETA/W(2)/(1+Kd(2)))
	  QIA(4,iz)=shape(1)*(ETA/W(1)*Kd(1)/(1+Kd(1))
     .	   +ETA/W(2)*Kd(2)/(1+Kd(2)))
	  QIA(2,iz)=shape(2)*(ETA/W(3)/(1+Kd(3))+ETA/W(4)/(1+Kd(4)))
	  QIA(3,iz)=shape(2)*(ETA/W(3)*Kd(3)/(1+Kd(3))
     .	   +ETA/W(4)*Kd(4)/(1+Kd(4)))+shape(3)*ETA/W(5)

c Normalization
	  QIA_TOT=QIA(1,iz)+QIA(2,iz)+QIA(3,iz)+QIA(4,iz) 
	  if ((QIA_TOT.gt.1.e-20).and.(QTI(iz).gt.1.e-20)) then
		do i=1,4
			QIA(i,iz)=QIA(i,iz)*QTI(iz)/QIA_TOT
		enddo
	  else
		do i=1,4
			QIA(i,iz)=0.
		enddo
		QTI(iz)=0.
	  endif

	ENDDO

	do iz=1,nalt
		QTI_(iz)=QTI((iz-1)*(nk+1)+1)*1.e6
		do i=1,4
		QIA_(i,iz)=QIA(i,(iz-1)*(nk+1)+1)*1.e6
		enddo
	enddo

	return
	end
c------------------------------------------------
