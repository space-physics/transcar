	program trace_fte

	integer transcar_dat,kiappel
	parameter(transcar_dat=10)
	Integer itube,iyd,iydfin,iyddeb,iydtmp,ipos
	Real longeo,latgeo,tempsdeb,postint,tempsfin
	real*8 dlontmp,dlattmp
	Real*8 dlongeo,dlatgeo
	Real temps,tps,dtconv,dto,dtref, postinto
	real latgeo_ini,longeo_ini
	Real*8 pot
	Logical flgdt,flgpot



        real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0
	real Bmag,dipangle,Enord,Eest,vperpnord,vperpest,vhorizon,vpara
	real orient,chi
	real B,dip,or,ddp,Jtop
        integer ikp

        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi,
     &			Fe0,Ee0,Fi0,Ei0,
     &			Bmag,dipangle,Enord,Eest,
     &			vperpnord,vperpest,vhorizon,vpara,ddp,Jtop




	open(transcar_dat, file='dir.input/DATCAR.fte')
        rewind(transcar_dat)
        read(transcar_dat,*)kiappel
        read(transcar_dat,'(a)')

	open(31,file='trace_fte.dat',form='formatted',status='unknown')

	read(transcar_dat,*)dto     	! pas d'integration numerique
	read(transcar_dat,*)sortie      ! intervalle de temps entre deux sorties 
	read(transcar_dat,*)iyd_ini     ! date de la periode simulee
	read(transcar_dat,*)tempsini    ! UT de debut (en heure)
	read(transcar_dat,*)tempslim    ! UT limite (en heure)
	read(transcar_dat,*)f107
	read(transcar_dat,*)f107a
	read(transcar_dat,*)
	

        read(transcar_dat,*)jpreci

	do ipos=1,4
	
	 read(transcar_dat,*)latgeo_ini(ipos),longeo_ini(ipos)

	end do
	
	tempsconv_1=0.
	read(transcar_dat,*)tempsconv_1	! duree de la convection en amont (<= 0 si pas de convection)

	tempsconv=0.
	read(transcar_dat,*)tempsconv	! duree de la convection en aval (<= 0 si pas de convection)

	read(transcar_dat,*)step        ! intervalle de temps entre deux tubes
	read(transcar_dat,*)postinto    ! intervalle de temps entre deux appel a transelec

	close(transcar_dat)
	
--------------------------------	
Debut boucle sur les 4 positions
--------------------------------

	do ipos=1,4
	

	tempsini=tempsini*3600.
	
	if (tempsconv_1.gt.0.) then
	  tempsconv_1=tempsconv_1*3600.
	endif
	
	if (tempsconv.gt.0.) then
	  tempsconv=tempsconv*3600.
	endif
	
	  itube=-1
	  ntubmax=int(tempslim/step)+1


	do while (itube.le.ntubmax)	! debut de la boucle sur les lignes de champs

	  itube=itube+1

cccccccccccccccccccccccccccccccccccccccccccccccccccc
C[      Densities and Velocities profiles initialisation
cccccccccccccccccccccccccccccccccccccccccccccccccccc


c determination du point de depart de la convection

	iw=0
	  dlongeo=longeo_ini(ipos)
	  dlatgeo=latgeo_ini(ipos)
	  longeo=longeo_ini(ipos)
	  latgeo=latgeo_ini(ipos)
	  iyd=iyd_ini
	  iydfin=iyd_ini
	  iydtube=iydfin

c le temps est fice a tempsini

c	  tempstube=tempsini+itube*step 		! temps de resolution du tube
	  tempstube=tempsini

	  do while (tempstube.ge.86400.)		!
	    tempstube=tempstube-86400.			! on ramene tempstube entre 0 et 24 heures
	    iydtube=iydtube+1				!
	  enddo						!

c	  tempsfin=tempstube+tempsconv		 	! temps desire en fin de convection (version mesure radar)
c	  iydfin=iydtube
	  do while (tempsfin.ge.86400.)			!
	    tempsfin=tempsfin-86400.			! on ramene tempsfin entre 0 et 24 heures
	    iydfin=iydfin+1				!
	  enddo						!

	  tempsdeb=tempstube-tempsconv_1		! temps desire en debut de convection (version mesure radar)
	  iyddeb=iydtube
	  do while (tempsdeb.lt.0.)			!
	    tempsdeb=tempsdeb+86400.			! on ramene tempsdeb entre 0 et 24 heures
	    iyddeb=iyddeb-1				!
	  enddo						!


	  temps=tempsfin
	  iyd=iydfin
	  flgpot=.true.
	  dtref=0.
	      call convec(iyd,temps,dlongeo,dlatgeo,0.,pot,flgpot)
	      
	      write(31,100) temps,tmag,latmag,dlongeo,dlatgeo
	  do while (temps.gt.tempsdeb.or.iyd.ne.iyddeb)
	    tps=temps
	    iydtmp=iyd
	    flgdt=.true.
	    do while (flgdt)
	      dlontmp=dlongeo
	      dlattmp=dlatgeo
	      call pas_de_temps(iydtmp,tps,dtconv,postint,dto,postinto)
c	dtconv=dto
	      call convec_1(iydtmp,tps,dlontmp,dlattmp,dtconv,pot,flgpot)
	      tps=tps-dtconv
	      if (tps.lt.0.) then			!
	        tps=tps+86400.			! on ramene temps entre 0 et 24 heures
	        iydtmp=iydtmp-1			!
	      endif				!
	      flgdt=.not.(dtref.eq.dtconv)
	      dtref=dtconv
	    enddo
c	    flgpot=.false.
	    temps=tps
	    iyd=iydtmp
	    dlatgeo=dlattmp
	    dlongeo=dlontmp
	    iw=iw+1
  	    if (mod(iw,300).eq.0) then
	      write(31,100) temps,tmag,latmag,dlongeo,dlatgeo
	    endif
100	format(f6.0,8(1x,g12.5))
	  enddo
	      write(31,100) temps,tmag,latmag,dlongeo,dlatgeo
	enddo
	
	
	enddo

------------------------------
Fin boucle sur les 4 positions
------------------------------
	
	close (31)

	end
