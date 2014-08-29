	subroutine conv_ar(itube,dlongeo,dlatgeo,
     &			   dlonmag,dlatmag,dlonref,
     &			   iyddeb,tempsdeb,
     &                     iydfin,tempsfin,
     &			   iyd   ,temps   ,dto,postinto)

	Implicit none

	Integer itube,iyd,iydfin,iyddeb,iydtmp
	Real longeo,latgeo,kp,ap(7)
	Real*8 tempsdeb,postint,tempsfin
	real*8 dlontmp,dlattmp
	Real*8 dlongeo,dlatgeo,dlonmag,dlatmag,dlonref
	Real*8 temps,tps,dtconv,dto,dtref, postinto
	Real*8 pot
	Logical flgdt,flgpot

	temps=tempsfin
	iyd=iydfin
	flgpot=.true.
	dtref=0.
	do while (temps.gt.tempsdeb.or.iyd.ne.iyddeb)
	  tps=temps
	  iydtmp=iyd
	  flgdt=.true.
	  do while (flgdt)
	    dlontmp=dlongeo
	    dlattmp=dlatgeo
	    call pas_de_temps(iydtmp,tps,dtconv,postint,dto,postinto)
c	dtconv=dto
c	    call aptime(iydtmp,tps,kp,ap)
	    kp=1
c	    print*,dlontmp,dlattmp,
c     &			  dlonmag,dlatmag
	    call convec_1(iydtmp,tps,kp,dlontmp,dlattmp,
     &			  dlonmag,dlatmag,dlonref,dtconv,pot,flgpot)
	    tps=tps-dtconv
	    if (tps.lt.0.d0) then			!
	      tps=tps+86400.d0			! on ramene temps entre 0 et 24 heures
	      iydtmp=iydtmp-1			!
	    endif				!
	    flgdt=.not.(dtref.eq.dtconv)
	    dtref=dtconv
	  enddo
c	  flgpot=.false.
	  temps=tps
	  iyd=iydtmp
	  dlatgeo=dlattmp
	  dlongeo=dlontmp
c	print*,temps,dlatgeo,dlongeo,pot
	enddo
	return
	end
