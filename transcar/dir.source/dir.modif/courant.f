	subroutine courant(iyd,tu,kp,tmag,lat,Jsup)
c	implicit none

        integer npt
        parameter(npt=1000)
        integer iyd,ndeg,mdeg,ierr
	real*8 tu
	real kp,Jsup
	real*8 tmag,lon,lat
	complex*16 phicourant(npt)
        real*8 latequi,Lmin,Lmax
	
	save ierr

	call coef_cour(iyd,tu,ndeg,mdeg,phicourant,
     &		       Lmin,Lmax,latequi,ierr)
        Jsup=0.d0
	if (ierr.eq.0) then
          lon=15.d0*tmag
	  Jsup=val_fit(lon,lat,ndeg,mdeg,phicourant,
     &		       Lmin,Lmax,latequi)

c	dans le fichier initial le courant est en 10^-7 A/m2
	  Jsup=Jsup/10.

	endif
	
	return
	end
	
	
