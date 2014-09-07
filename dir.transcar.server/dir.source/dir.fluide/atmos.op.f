	subroutine atmos(iyd,ces,stl,z,glat,glong,jpreci,f107,
     &			ap,Ne,Te,Tj,nx,kiappel,file_cond)
 



        include 'TRANSPORT.INC'

        logical isnan,isnant,flgnan
        common/nan/flgnan
        external isnan,isnant

C	SIGNIFICATION DES VARIABLES DE CE SOUS-PROGRAMME

C	liste d'appel
c       -------------
C	IYD	jour de l'annee sous la forme YYJJJ, 
c 		avec YY=annee et JJJ=jour dans l'annee
C	SEC	temps universel 				(s)
C	STL	temps local	 				(heure)
C	Z	tableau des altitudes 				(km)
C	GLAT	latitude geographique				(degres)
C	GLONG	longitude geographique				(degres)
C	F107A	flux solaire moyen sur les 81 derniers jours
C	F107	flux solaire instantane de la veille
C	AP	tableau des ap (convention MSIS86)
C	Ne	tableau des concentrations electroniques
c 		(unites normalisees)
C	Te	tableau des temperatures electroniques			
c 		(unites normalisees)
C	Tj	tableau des temperatures de l'ion O+			
c 		(unites normalisees)
C	NX	nombre d'altitudes resolues
c     jpreci = 0 if sun only, 1 if prec. only, 2 if both of them

c
C	atmosphere neutre
c       ----------------- 
C	NH	concentration de l'hydrogene atomique		(cm-3)
C	NO	concentration de l'oxygene atomique		(cm-3)
C	NO2	concentration de l'oxygene moleculaire		(cm-3)
C	NN2	concentration de l'azote moleculaire		(cm-3)
C	NN	concentration de l'azote atomique		(cm-3)
C	TN	temperature neutre				(K)
C	UN	vitesse neutre alignee a B			(cm.s-1)
C	VN	vitesse neutre vers l'est	 		(cm.s-1)
C	WN	vitesse neutre vers le nord		 	(cm.s-1)
C	Q_NH	flux de chaleur reduit de H		(erg.cm.s-1)
C	Q_NO	flux de chaleur reduit de O		(erg.cm.s-1)
C	Q_NO2	flux de chaleur reduit de O2		(erg.cm.s-1)
C	Q_NN2	flux de chaleur reduit de N2		(erg.cm.s-1)
C	Q_NN	flux de chaleur reduit de N		(erg.cm.s-1)

C	DTINF	ecart sur Texo donne par MSIS86			(K)
C	DTZ	gradient de Tn donne par MSIS86			(K.km-1)
C	GRADTN	gradient de Tn donne par MSIS86			(K.cm-1)

C	T	tableau pour appel de MSIS86
C	D	tableau pour appel de MSIS86
C	W	tableau pour appel de  HWM87
C	COFH	coefficient pour la concentration de H
C	COFO 	coefficient pour la concentration de O
C	COFN	coefficient pour les autres concentrations neutres
C	T69	intermediare de calcul pour les flux de chaleur
C	FACT	coefficient pour l'extrapolation des concentrations
C	ZCIRA	altitude limite pour le calcul par extrapolation
C	ICIRA	indice de l'altitude limite zcira

C	perturbations externes
c       ----------------------- 
C	PERTURB	routine decrivant la perturbation
C	PERT	coefficient de la perturbation		(0<pert<1)

C	B0	amplitude du champ magnetique au sol	(Gauss)
C	COSI0	cosinus du dip angle
C	SINI0	sinus du dip angle
C	J0	amplitude du courant aligne		(�.m-2)
C	Eest0	amplitude de la composante
c		est du champ electrique de convection	(mV.m-1)
C	Enord	amplitude de la composante
c		nord du champ electrique de convection	(mV.m-1)

C	B	tableau des intensites du champ magnetique	(Gauss)
C	VM	tableau de la composante est
c		de la derive magnetique				(cm.s-1)
C	WM	tableau de la composante nord
c		de la derive magnetique				(cm.s-1)
C	JJ	tableau du courant aligne applique		(�.m-2)

C	EE	energie moyenne des electrons precipitants	(eV)
C	FE	flux des electrons precipitants		(m-2.s-1)

C	conditions limites pour le flux de chaleur electronique
c       -------------------------------------------------------
C	VARTEMP	coefficient decrivant l'evolution temporelle
c		du flux de chaleur a l'altitude superieure
C	TEMPLOC	intermediaure de calcul pour vartemp

C	QETOPMIN	valeur minimale du flux de chaleur
c			electronique a l'altitude superieure
c			(erg.cm-2.s-1)
C	QETOPMAX	valeur maximale du flux de chaleur
c			electronique a l'altitude superieure		
c			(erg.cm-2.s-1)
C	QETOP	valeur du flux de chaleur
c		electronique a l'altitude superieure			
c			(erg.cm-2.s-1)

C	variables liees a l'evolution en temps
c       -------------------------------------- 
C	CHI	angle zenithal					(rad)
C	NHEU	heure universelle
C	MIN	minute universelle
C	NAN	annee
C	ajour	jour dans l'annee

C	parametres

C	INDLIM		indice de l'altitude limite utilisee pour les 
c 			calcul de Jean
C	INDLIM_1	indice de l'altitude avant l'altitude limite
C	ZLIM		valeur de l'altitude limite
C	ZLIM_1		veleur de l'altitude precedent l'altitude limite

C	R0	coefficient de normalisation des altitudes	(cm)
C	G0	coefficient de normalisation de la gravite    (cm-2.s-1)
C	T0	coefficient de normalisation du temps		(s)

C	N_0	coefficient de normalisation des concentrations	(cm-3)
C	T_0	coefficient de normalisation des temperatures	(K)
C	P_0	coefficient de normalisation des pressions    (erg.cm-3)
C	CI0	coefficient de normalisation de la vitesse de H+ 	
c		(cm.s-1)
C	CJ0	coefficient de normalisation de la vitesse de O+ 	
c		(cm.s-1)
C	CK0	coefficient de normalisation de la vitesse de N2+ 	
c		(cm.s-1)
C	CL0	coefficient de normalisation de la vitesse de O2+ 	
c		(cm.s-1)
C	CM0	coefficient de normalisation de la vitesse de NO+	
c		(cm.s-1)
C	CE0	coefficient de normalisation de la vitesse de e- 	
c		(cm.s-1)
C	QI0	coefficient de normalisation du flux de chaleur de H+   
c		(erg.cm-2.-s1)
C	QJ0	coefficient de normalisation du flux de chaleur de O+   
c		(erg.cm-2.-s1)
C	QE0	coefficient de normalisation du flux de chaleur de e-	
c		(erg.cm-2.-s1)

C	RE	rayon de la terre				(km)
C	M0	masse de H+ (amu)				(g)
C	ME	masse de le electron				(g)
C	KB	constante de Boltzmann			    (erg.K-1)

C	SRCFIC	repertoire ou se trouvent les fichiers temporaires de 
c		transport cinetique

C	PRECIPITATION	routine determinant les electrons precipitants 
c		        (+ les ions precipitants)
C	TRANSPORT	routine calculant le transport des electrons 
c 			suprathermaux

C       Added to simulate an assumed neutral Hot O profile
C       Matt Zettergren, 13/09/2004
C
C	CHKOHOT	Check the Hot O input file for info on simulating hot O
C       NOHOT   The neutral Hot O concentration (1/cm3)
C       Q_NOHOT The heat flow for neutral Hot O
C       TNOHOT  The neutral Hot O temperature (K)
C       PCTOHOT Percentage of Hot O compared with Cold O (NOHOT/NO * 100)
C	ZNOHOTREF Reference altitude for density
C	NOHOTREF Reference density
C       TNOHOTINF Exospheric Temperature of the Hot O
C	ZTNOHOTREF Reference Hot O altitude for temperature
C	TNOHTOREF Reference Hot O temperature
C	TNOHOTDECAY Decay rate for bates temperature profile
C       FLAGOHOT False to omit Hot O, True to simulate Hot O
C   
C       --MZ


	parameter (npt=500)

	integer file_cond
	real z(npt),Ne(npt),Te(npt),Tj(npt)
	real t(2),d(8),f107(3),ap(7),w(2)
	real sec,stl,ces
        real secref,stlref,zref,Tref
	real nu_omega

	real Nh(npt),No(npt),No2(npt),Nn2(npt),Nn(npt),Tn(npt)
	real Un(npt),Vn(npt),Wn(npt)
	common /neutral/ Nh,No,No2,Nn2,Nn,Tn,Un,Vn,Wn

	real q_Nh(npt),q_No(npt),q_No2(npt),q_Nn2(npt),q_Nn(npt)
	common/fluneut/q_Nh,q_No,q_No2,q_Nn2,q_Nn

	real dTinf,vmerid_inf,vzonal_inf,vmerid_200,vzonal_200
        common/ventexo/vmerid_inf,vzonal_inf,vmerid_200,vzonal_200

	common/exo/dTinf
	common/grad/dtz

        real zlim,zlim_1
        common/limite/indlim,indlim_1,zlim,zlim_1

	real Vm(npt),Wm(npt),Vm_2(npt),omega(npt)
	common/champ/Vm,Wm,omega,Vm_2

	real JJ(npt)
	common/J_aligne/JJ

	real vartemp,Qetop
	common/fluxtop/Qetop

	real N_0,T_0,Ci0,Cj0,Ck0,Cl0,Cm0,Ce0,Pi0,Pj0,Pe0,Qi0,Qj0,Qe0
        common/param/	N_0,T_0,P_0,Ci0,Cj0,Ck0,Cl0,Cm0,Ce0,
     &			Qi0,Qj0,Qe0
     &			
        real R0,g0,t0
        common/adim/R0,t0,G0
  
        real dt
        common/temps/dt

	real Re,zsup,zcira,m0,kb,me
	data Re,zsup,zcira/6378.,1000.,500./
	data m0,me,kb/1.667e-24,9.11e-28,1.38e-16/

	real J0,B0,sinI0,cosI0,chi
        data B0,sinI0,cosI0/.542,.9848,.1736/

	real deg2rad
	data deg2rad/0.01745329251994/

        real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0,chi0
        real Bmag,dipangle,Enord,Eest,vperpnord,vperpest,vhorizon,vpara
	real kp,ddp,Jtop
        integer ikp

        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi0,
     &                  Fe0,Ee0,Fi0,Ei0,
     &                  Bmag,dipangle,Enord,Eest,
     &                  vperpnord,vperpest,vhorizon,vpara,ddp,Jtop
        real cofn2,cofo2
        common/buff2/ cofn2,cofo2

	real coef_flux
        common/prec/coef_flux

	real Pr(npt,3),rhot(npt)
        real*8 tu,dlon,dlat
c
        real Ph(npt),Po(npt),Po2(npt),Pn2(npt),Pn(npt)
        common/prodion/Ph,Po,Po2,Pn2,Pn
        
        real Nes(npt),Jes(npt),Tes(npt),qes(npt)
        common /supra/Nes,Jes,Tes,qes

c       PIERRE LOUIS : SI TU METS CES MOMENTS DANS UN COMMON, NOUBLIE
C       PAS DE ME RENVOYER LE CODE MODIFIE.

        real Heat(npt)
        common/heate/Heat

        real colen
        common/coeflen/colen
        common/cool/coloss

        logical flg_fct
        common/fl/flag_fct

        real burn
        common/burn/burnside

	real temps,En,Ee,Eno,Eeo
cc 	 
        logical flag,flagatmos,flgini,flgread
	common /atm/flagatmos

        data flag/.true./
        data flgini/.true./
 	integer kiappel

	real perturb,pert

	real dTinf0(14),qetop0(14)
	real cofo0(14),cofh0(14)

	real vparaB
        common /avril/vparaB


	real rap
	common /rap_prec/rap


        logical precflag
        integer prect
        common/precstate/precflag,prect

        real timeser(256), edist(256)
        real fluxdist(256,256)
        integer ntimeser,nfluxdist,precint,precext
        common /precdist/ nfluxdist,ntimeser,timeser,edist,fluxdist,precint,precext

	real e1,f1
        common /E930216/e1,f1,coefchamp

	real*8 a0,c0,s0,x0,y0,xpos,ypos,pos0

	
	logical flg_err
	common/err/flg_err	

            
C       Added to simulate an assumed neutral Hot O profile
C       Matt Zettergren, 13/09/2004
        logical chkOHot,flagOHot
	common /flagsOHot/ chkOHot,flagOHot
        real NOHot(npt),TnOHot(npt),q_NOHot(npt),POHot(npt),fluxhot(npt)
	common /OHot/ NOHot,TnOHot,q_NOHot,POHot,fluxhot
	real secOHot
        real zNOHotref,pctOHot,NOHotRef
	real zTnOHotRef,TnOHotRef,TnOHotDecay,TnOHotInf
	common /OHotParams/ zNOHotref,pctOHot,NOHotRef,zTnOHotRef,
     &		TnOHotRef,TnOHotDecay,TnOHotInf
C       --MZ
	
C       Added in to account for RBR electron precipitation
        real Eprec,Fprec
        common /RBR/ Eprec,Fprec
        real tstartprec,tstopprec
        common /prectimes/ tstartprec,tstopprec
C       --Z


	sec=ces
        kp=ap2kp(ap(1))
        ikp=kp*3
	nan=iyd/1000
	njour=mod(iyd,1000)
	chi=acos(coskhi(glat,glong,sec/3600.,nan,njour,2))
	chi0=chi*180./3.14159265

	rap=1.
!	if (njour.eq.48) then
!	  open(83,file=data_path(1:lpath_data)
!     &                       //'dir.fluide/var_rap.dat',
!     &		form='formatted',status='unknown',iostat=ierr)
    
!	  if (ierr.eq.0) then
!	    tf=0.
!	    do while(ces.ge.tf) 
!	      read(83,*)td,tf,rd,rf
!	    enddo
!	    close(83)
!	    rap=(rf-rd)/(tf-td)*(ces-td)+rd
!	    rap=(rap-1.)/10.+1.
!	rap=1.
!	  endif
!	endif


        temploc=(91.-chi*180./3.1415926)
        vartemp=(1+tanh(temploc*3.))/2.
	J0=0.
!	Qetop  = -5.e-4-4.e-3*exp(-(chi0-91.)**2/30.)
!	Qetop=-2.e-3
	
	  vmerid_inf=0.
	  vzonal_inf=0.
	  vmerid_200=0.
	  vzonal_200=0.
	  dTinf=0.
!	  cofo=1.
!	  cofn=1.
!	  cofn2=cofn
!	  cofo2=cofn
!	  cofh=1.
        goto 987

987	continue	


	dlon=dble(15.*tmag)
	dlat=dble(latmag)
	tu=dble(sec)
	call courant(iyd,tu,kp,dlon,dlat,Jtop)

	
!       Added for simulation of electron precpitation
        if(sec.ge.tstartprec.and.sec.le.tstopprec) then
          precflag=.true.
          jpreci=2

          prect=1
          do while(sec .gt. timeser(prect+1) .and. prect .lt. ntimeser)
            prect=prect+1
          enddo
          if(prect .eq. ntimeser) then
            precflag=.false.
            jpreci=0
          endif
        else
          precflag=.false.
          jpreci=0
        endif
!       --MZ


!       Added to simulate an assumed neutral Hot O profile
!       Matt Zettergren, 21/09/2004
!       Unfortunately I can't just get rid of this because of the collision
!       terms I've added to the transport equations.  I'll just set everything
!       to zero instead for now.

!	!Get Hot O input information
!	if(chkOHot) then
!	  open(917,file='/home/mattz/dir.transcar/dir.input/DATOHOTnull')
!	  read(917,*) flagOHot

!	  if(flagOHot) then
!	    read(917,*) zNOHotRef
!	    read(917,*) pctOHot
!	    read(917,*) zTnOHotRef
!	    read(917,*) TnOHotRef
!	    read(917,*) TnOHotDecay
!	    read(917,*) TnOHotInf
!	  endif
!
!	  close(917)
!	endif

        flagOHot=1
        zNOHotRef=400.
        pctOHot=0.
        zTnOHotRef=90.
        TnOHotRef=500.
        TnOHotDecay=0.005
        TnOHotInf=4000.

	if(flagOHot) then
	  !Find the reference hot O density
	  secOHot=sec
	  call gtd6(iyd,secOHot,zNOHotRef,glat,glong,stl,
     &		f107(3),f107(2),ap,48,d,t)
	  NOHotRef=pctOHot*0.01*d(2)

	  !Calculate the Hot O density over the desired altitude range
          do i=1,nx
	    !Parameterized Bates profile
!	    TnOHot(i)=TnOHotInf-(TnOHotInf-TnOHotRef)*
!     &		exp(-TnOHotDecay*(z(i)-zTnOHotRef))

            TnOHot(i)=4000.

	    !Diffusive equilibrium profile
	    NOHot(i)=NOHotRef*exp(-16.*m0*g0*1.e5/kb/TnOHot(i)*(z(i)-zNOHotRef)
     &			/(z(i)/Re+1.)**2)
            
	    !Heat Flow is negligible
	    q_NOHot(i)=0.
	    fluxhot(i)=0.!sqrt(TnOHot(i)/T_0)
          enddo
        endif
!       --MZ


	do i=1,nx

     
	if (flagatmos) then
c	if (flgini) then

	  if (z(i).le.zsup) indlim=i

C       Added to make sure that neutral heat flow is zero
C       Matt Zettergren, 21/09/2004
	  q_NH(i)=0
	  q_NO(i)=0
	  q_NN(i)=0
	  q_NO2(i)=0
	  q_NN2(i)=0
C       --MZ

	  if (z(i).le.zcira) then
	    icira=i
c 	    ... msis90
 	    call gtd6(iyd,sec,z(i),glat,glong,stl,
     &		f107(3),f107(2),ap,48,d,t)
	    Nh(i)=d(7)*cofh
	    No(i)=d(2)*cofo
	    No2(i)=d(4)*cofo2
	    Nn2(i)=d(3)*cofn2
            Nn(i)=d(8)*cofn
	    Tn(i)=t(2)
	  else
	    fact=exp(-m0*g0*Re**2/kb/Tn(icira)*(z(i)-z(icira))
     &			/(z(i)+Re)/(z(icira)+Re)*1.e5)
	    Nh(i)=Nh(icira)*fact
	    No(i)=No(icira)*fact**16
	    No2(i)=No2(icira)*fact**32
	    Nn2(i)=Nn2(icira)*fact**28
            Nn(i)=Nn(icira)*fact**14
	    Tn(i)=Tn(icira)
	  endif



!------Neutral wind calculation
!	  call GWS5(iyd,sec,z(i),glat,glong,stl,
!     &		f107(3),f107(2),ap,w)

!	Un est la composante parallele du vent neutre
!        Un(i)=-w(1)*1e2*sin(dipangle*deg2rad)


C	Vn est la composante est du vent neutre

c	  Vn(i)=w(2)*1e2

C       Wn est la composante nord du vent neutre

c          Wn(i)=w(1)*1e2*cos(dipangle*deg2rad)

	endif
!-------MZ

	  B=B0*(Re/(Re+z(i)))**3
	  omega(i)=1.6e-19*B*1.e-4/1.667e-27*t0

C       Vm est la composante est de la derive magnetique

	  Vm(i)=-Enord/Bmag*1000.-Vn(i)

C       Wm est la composante nord de la derive magnetique

	  Wm(i)=Eest/Bmag*1000.-Wn(i)


	Vm_2(i)=Vm(i)**2+Wm(i)**2

c	J0 est en �/m-2
	  JJ(i)=J0/1.6e-9*((800.+Re)/(z(i)+Re))**3

          nu_omega=1.e-9*Nn2(i)*t0/omega(i)*30.5
          coef_cour=1./(1.+nu_omega**2)

          JJ(i)=JJ(i)*coef_cour

	enddo

	if (flagatmos) then
          indlim_1=indlim-1
          zlim=z(indlim)
          zlim_1=z(indlim_1)
          
c
c
c 	Transport cinetique. Kinetic transport
          if (kiappel.eq.1 .or. kiappel.eq.2)then

	    if (flgini) then
	    flgnan=.true.
               if (isnan(glat)) then
                 print*,'probleme avec glat'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnan(glong)) then
                 print*,'probleme avec glong'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnan(sec)) then
                 print*,'probleme avec sec'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnan(stl)) then
                 print*,'probleme avec stl'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnant(ne,nx)) then
                 print*,'probleme avec ne'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnant(te,nx)) then
                 print*,'probleme avec te'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnant(tj,nx)) then
                 print*,'probleme avec tj'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnant(nh,nx)) then
                 print*,'probleme avec nh'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnant(no,nx)) then
                 print*,'probleme avec no'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnant(no2,nx)) then
                 print*,'probleme avec no2'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnant(nn2,nx)) then
                 print*,'probleme avec nn2'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnant(nn,nx)) then
                 print*,'probleme avec nn'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               elseif (isnant(tn,nx)) then
                 print*,'probleme avec tn'
                 print*,npt,iyd,sec,glat,glong,stl,f107,ap,chi
                 stop 'erreur NaN'
               endif
               flgnan=.false.
               flg_err=.false.
c               if (sec.ge.280..and.sec.lt.3000.) then
c                 flg_err=.true.
c 		endif

               call transelec(npt,iyd,sec,glat,glong,stl,f107,
     .          ap,chi,Ne,Te,Tj,nx,Nh,No,No2,Nn2,Nn,Tn,indlim,jpreci,
     .          N_0,T_0,kiappel,zlim,zlim_1,z,Heat,Ph,Po,Po2,Pn2,Pn,
     .          Nes,Jes,Tes,qes)

	    endif

          endif

	endif

	do i=indlim+1,nx
	  Po(i) =Po(indlim)/No(indlim)*No(i)
	  Pn2(i)=Pn2(indlim)/Nn2(indlim)*Nn2(i)
	  Po2(i)=Po2(indlim)/No2(indlim)*No2(i)
	  Pn(i)=Pn(indlim)/Nn(indlim)*Nn(i)
	  Ph(i)=Ph(indlim)/Nh(indlim)*Nh(i)
	enddo

	  do i=1,nx								!MZ
	    POHot(i)=0.!Po(i)*NOHot(i)/No(i)					!MZ
	  enddo									!MZ

c	if (flgini) flgini=.false.

	return
	end
