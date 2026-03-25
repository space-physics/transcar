	programme conjonction


	Implicit none
	integer nbouclemax
	parameter(nbouclemax=10)

	real temps,delai,tempsdeb,tempsfin,dtconv,deltat
	integer iyd,iydfin,iyddeb,nboucle
	Real longeo,latgeo,dlon,dlat,lon_ref,lat_ref
	real longeo_tro,latgeo_tro,longeo_esr,latgeo_esr
	real lonmag_tro,latmag_tro,lonmag_esr,latmag_esr,lonref
	Real pot
	Logical flgpot,conjonction,sens
	real*8 dlatgeo,dlongeo,dlatmag,dlonmag,dlonref

        real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0
        real Bmag,dipangle,Enord,Eest,vperpnord,vperpest,vhorizon,vpara
        real orient,chi
        real B,dip,or
        integer ikp

        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi,
     &                  Fe0,Ee0,Fi0,Ei0,
     &                  Bmag,dipangle,Enord,Eest,
     &                  vperpnord,vperpest,vhorizon,vpara


	data longeo_tro,latgeo_tro/19.20,69.66/
	data longeo_esr,latgeo_esr/16.03,78.15/
	data dlon,dlat/5.,1./
	data flgpot/.true./
	data dtconv/1./
	data deltat/1800./
	data delai/5./

	dlatgeo=latgeo_esr
	dlongeo=longeo_esr
	call geo2mag(dlatgeo,dlongeo,dlatmag,dlonmag,dlonref)
	latmag_esr=dlatmag
	lonmag_esr=dlonmag
	dlatgeo=latgeo_tro
	dlongeo=longeo_tro
	call geo2mag(dlatgeo,dlongeo,dlatmag,dlonmag,dlonref)
	latmag_tro=dlatmag
	lonmag_tro=dlonmag
	lonref=dlonref

1	iyd=98001
	print*,'jour pour la conjonction (dans le format julien AAJJJ) : '
	read*,iyd
	ikp=1
	print*,'Type d''activite (kp = 1, 3 ou 5) : '
	read*,ikp
	sens=.true.
	print*,'sens de la conjonction (T = EISCAT -> ESR)'
	read*,sens
	if (sens) then
	  iyddeb=iyd
	  tempsdeb=0.
	else
	  iyddeb=iyd
	  tempsdeb=15.*3600.
	endif

	delai=delai*3600.
	do while (tempsdeb.lt.86400)
	  print*,'temps de depart = ',tempsdeb
	  tempsfin=tempsdeb+delai
	  iydfin=iyddeb
	  if (tempsfin.ge.86400.) then	!
	    tempsfin=tempsfin-86400.		! on ramene tempsfin entre 0 et 24 heures
	    iydfin=iydfin+1			!
	  endif				!
	  temps=tempsdeb
	  iyd=iyddeb
	  flgpot=.true.
	  if (sens) then
	    latgeo=latgeo_tro
	    longeo=longeo_tro
	    lat_ref=latmag_esr
	    lon_ref=lonmag_esr
	  else
	    latgeo=latgeo_esr
	    longeo=longeo_esr
	    lat_ref=latmag_tro
	    lon_ref=lonmag_tro
	  endif
	  do while (temps.le.tempsfin.or.iyd.ne.iydfin)
	    call convec(temps,latgeo,longeo,dtconv,pot,flgpot)
	    temps=temps+dtconv
	    flgpot=.false.
	    if (temps.ge.86400.) then			!
	      temps=temps-86400.			! on ramene temps entre 0 et 24 heures
	      iyd=iyd+1				!
	    endif					!

	if (mod(int(temps),300).eq.0)
     &		 print*,longeo,latgeo,tmag,latmag,pot

	    if (abs(latmag-lat_ref).le.dlat) then
	        print*,'heure de conjonction a EISCAT = ',
     &		         tempsdeb/3600.,' TU'
	        print*,'heure de conjonction a ESR = ',
     &		         temps/3600.,' TU'
	        print*,'temps MLT a la conjonction ESR / EISCAT = ',
     &		        tmag,mod(temps/3600.+(lon_ref+lonref)/15.+24.,24.)
     	        pause
	        exit
c	      else
c	        exit
c	      endif
	    endif
	  enddo
	  tempsdeb=tempsdeb+deltat
	enddo
c
	end
