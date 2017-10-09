	program tube
	
	include 'TRANSPORT.INC'
	
	parameter(nbuf=8)

        parameter (nb_position_max=100)         ! Modif DA 02/02 2001
        logical multi_position                  ! Modif DA 02/02 2001
        real longeo_position(nb_position_max)   ! Modif DA 02/02 2001
        real latgeo_position(nb_position_max)   ! Modif DA 02/02 2001
        character lecture_lat_lon*80            ! Modif DA 02/02 2001


	integer kiappel
	Integer itube,iyd,iydfin,iyddeb,iydtmp
	Real longeo,latgeo,tempsdeb,postint,tempsfin,buffer(nbuf)
	real*8 dlontmp,dlattmp
	Real*8 dlongeo,dlatgeo,dlonmag,dlatmag,dlonref
	Real*8 temps,tps,dtconv,dto,dtref, postinto
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


        open(transcar_dat, 
     &       file='dir.input/DATCAR.tube',
     &       status='old')
        rewind(transcar_dat)

	open(31,file='trace.dat',form='unformatted',
     &                access='direct',status='new',recl=4*nbuf)
     	
     	irec=0

	read(transcar_dat,*)dto     	! pas d'integration numerique
	read(transcar_dat,*)sortie      ! intervalle de temps entre deux sorties
	read(transcar_dat,*)iyd_ini     ! date de la periode simulee
	read(transcar_dat,*)tempsini    ! UT de debut (en secondes)
	read(transcar_dat,*)tempslim    ! UT limite (en secondes)
	read(transcar_dat,*)iar    	! (iar=1, convection_1+convection; ar = 0,convection_1)
	

	read(transcar_dat,*)latgeo_ini,longeo_ini

	tempsconv_1=0.
	read(transcar_dat,*)tempsconv_1	! duree de la convection en amont (<= 0 si pas de convection)

	tempsconv=0.
	read(transcar_dat,*)tempsconv	! duree de la convection en aval (<= 0 si pas de convection)

	read(transcar_dat,*)step        ! intervalle de temps entre deux tubes

        if (latgeo_ini.gt.90.) then                                                ! Modif DA 0202 2001 - debut
           nb_position=0                                                ! Modif DA 0202 2001
           multi_position=.true.                                        ! Modif DA 0202 2001
           read(transcar_dat,*)latgeo_ini,longeo_ini                    ! Modif DA 0202 2001
           do while (latgeo_ini.le.90.)                                 ! Modif DA 0202 2001
              nb_position=nb_position+1                                 ! Modif DA 0202 2001
              latgeo_position(nb_position)=latgeo_ini                   ! Modif DA 0202 2001
              longeo_position(nb_position)=longeo_ini                   ! Modif DA 0202 2001
              lecture_lat_lon=' '                                       ! Modif DA 0202 2001
              read(transcar_dat,*)latgeo_ini,longeo_ini                 ! Modif DA 0202 2001
           enddo                                                        ! Modif DA 0202 2001
       else                                                             ! Modif DA 0202 2001
           nb_position=1                                                ! Modif DA 0202 2001
           multi_position=.false.                                       ! Modif DA 0202 2001
       endif                                                            ! Modif DA 0202 2001 - fin
	close(transcar_dat)

	
	  itube=0
	  ntubmax=int(tempslim/step)

	do while (itube.le.ntubmax)	! debut de la boucle sur les lignes de champs

          do i_position=1,nb_position             ! Modif DA 0202 2001

cccccccccccccccccccccccccccccccccccccccccccccccccccc
C[      Densities and Velocities profiles initialisation
cccccccccccccccccccccccccccccccccccccccccccccccccccc

c determination du point de depart de la convection

         if (multi_position) then                         ! Modif DA 0202 2001 - debut
            dlongeo=longeo_position(i_position)         ! Modif DA 0202 2001
            dlatgeo=latgeo_position(i_position)         ! Modif DA 0202 2001
            longeo=longeo_position(i_position)          ! Modif DA 0202 2001
            latgeo=latgeo_position(i_position)          ! Modif DA 0202 2001
         else                                           ! Modif DA 0202 2001
             dlongeo=longeo_ini
             dlatgeo=latgeo_ini
             longeo=longeo_ini
             latgeo=latgeo_ini
         endif                                          ! Modif DA 0202 2001 - fin
	  iyd=iyd_ini
	  iydfin=iyd_ini
	  iydtube=iydfin

	  tempstube=tempsini+itube*step 		! temps de resolution du tube
	  do while (tempstube.ge.86400.)		!
	    tempstube=tempstube-86400.			! on ramene tempstube entre 0 et 24 heures
	    iydtube=iydtube+1				!
	  enddo						!

	  tempsfin=tempstube+tempsconv		 	! temps desire en fin de convection (version mesure radar)
	  iydfin=iydtube
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
	  call geo2mag(dlatgeo,dlongeo,dlatmag,dlonmag,dlonref)
	      call convec(iyd,temps,kp,dlongeo,dlatgeo,
     &			    dlonmag,dlatmag,dlonref,
     &			  dto,pot,flgpot)
c	      write(31,100) temps,tmag,latmag,dlongeo,dlatgeo,pot
c	      write(*,100) temps,tmag,latmag,dlongeo,dlatgeo,pot
c	write(*,*) itube,i_position,tempsdeb,tempsfin,tempstube	

c	  write(31,100) itube,temps,tmag,latmag,dlongeo,dlatgeo,
c     &		pot,i_position
	  buffer(1)=itube
	  buffer(2)=temps
	  buffer(3)=tmag
	  buffer(4)=latmag
	  buffer(5)=dlongeo
	  buffer(6)=dlatgeo
	  buffer(7)=pot
	  buffer(8)=i_position
	  irec=irec+1
	  write(31,rec=irec) buffer
	  tps0=temps
	  do while (temps.gt.tempsdeb.or.iyd.ne.iyddeb)
	    tps=temps
	    iydtmp=iyd
	      dlontmp=dlongeo
	      dlattmp=dlatgeo
	      call convec_1(iydtmp,tps,kp,dlontmp,dlattmp,
     &			    dlonmag,dlatmag,dlonref,
     &			    dto,pot,flgpot)
	      tps=tps-dto
	      if (tps.lt.0.) then			!
	        tps=tps+86400.			! on ramene temps entre 0 et 24 heures
	        iydtmp=iydtmp-1			!
	        tps0=tps0+86400.
	      endif				!
	    temps=tps
	    iyd=iydtmp
	    dlatgeo=dlattmp
	    dlongeo=dlontmp
  	    if (tps0-tps.ge.sortie) then
c	      write(31,100) itube,temps,tmag,latmag,dlongeo,dlatgeo,
c     &		pot,i_position
	      buffer(1)=itube
	      buffer(2)=temps
	      buffer(3)=tmag
	      buffer(4)=latmag
	      buffer(5)=dlongeo
	      buffer(6)=dlatgeo
	      buffer(7)=pot
	      buffer(8)=i_position
	      irec=irec+1
	      write(31,rec=irec) buffer
c	      write(*,100) itube,temps,tmag,latmag,dlongeo,dlatgeo,pot,
c     &		i_position
	      tps0=tps
	    endif
100	format(i4,1x,f6.0,5(1x,g12.5),1x,i4)
	  enddo
c	      write(31,100) itube,temps,tmag,latmag,dlongeo,dlatgeo,
c     &		pot,i_position
	      irec=irec+1
	      write(31,rec=irec) buffer
c	      write(*,100) itube,temps,tmag,latmag,dlongeo,dlatgeo,pot,
c     &		i_position

	if (iar.gt.0) then
	tps0=tps-sortie
	  do while (temps.lt.tempsfin.or.iyd.ne.iydfin)
	    tps=temps
	    iydtmp=iyd
	      dlontmp=dlongeo
	      dlattmp=dlatgeo
	      call convec(iydtmp,tps,kp,dlontmp,dlattmp,
     &			    dlonmag,dlatmag,dlonref,
     &			    dto,pot,flgpot)
	      tps=tps+dto
	      if (tps.gt.86400.) then			!
	        tps=tps-86400.			! on ramene temps entre 0 et 24 heures
	        iydtmp=iydtmp+1			!
	        tps0=tps0-86400.
	      endif				!
	    temps=tps
	    iyd=iydtmp
	    dlatgeo=dlattmp
	    dlongeo=dlontmp
  	    if (tps-tps0.ge.sortie) then
c	      write(31,100) itube,temps,tmag,latmag,dlongeo,dlatgeo,
c     &		pot,i_position
	      buffer(1)=itube
	      buffer(2)=temps
	      buffer(3)=tmag
	      buffer(4)=latmag
	      buffer(5)=dlongeo
	      buffer(6)=dlatgeo
	      buffer(7)=pot
	      buffer(8)=i_position
	      irec=irec+1
	      write(31,rec=irec) buffer
c	      write(*,100) itube,temps,tmag,latmag,dlongeo,dlatgeo,
c     &		pot,i_position
	      tps0=tps
	    endif

	  enddo
c	      write(31,100) itube,temps,tmag,latmag,dlongeo,dlatgeo,
c     &		pot,i_position
	      irec=irec+1
	      write(31,rec=irec) buffer
c	      write(*,100) itube,temps,tmag,latmag,dlongeo,dlatgeo,
c     &		pot,i_position
          endif


        enddo           ! Modif DA 0202 2001: fin de boucle sur i_position
	  itube=itube+1
	enddo
	close (31)

	end
