	subroutine potentiel(iyd,tu,kp,tmag,lat,Eest,Enord,pot,ddp)
	
	implicit none

        integer npt
        parameter(npt=1000)
        integer iyd,ndeg,mdeg,ierr
	real*8 tu
	real kp,ddp
	real*8 tmag,lon,lat
	complex*16 phipot(npt)
        real*8 latequi,Lmin,Lmax
        real*8 Enord,Eest,pot

	call coef_pot(iyd,tu,kp,ndeg,mdeg,phipot,
     &		       Lmin,Lmax,latequi,ddp,ierr)

        lon=15.d0*tmag

	call val_pot(lon,lat,ndeg,mdeg,phipot,
     &			   Lmin,Lmax,latequi,
     &			   Enord,Eest,pot)

	return
	end
