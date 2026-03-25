	subroutine pas_de_temps(iyd,temps,dt,postint,dto,postinto)
c
c       Permet de modifier le pas de temps d'integration.
c
	integer iyd
	real*8 dt,temps,duree,postint,dto,postinto
	real Eperp

        real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0,chi0
        real Bmag,dipangle,Enord,Eest,vperpnord,vperpest,vhorizon,vpara
	real ddp,Jtop
        integer ikp
        real log2
        data log2_1/1.4427/
        real dt_max
        common/CFL/dt_max


        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi0,
     &                  Fe0,Ee0,Fi0,Ei0,
     &                  Bmag,dipangle,Enord,Eest,
     &                  vperpnord,vperpest,vhorizon,vpara,ddp,Jtop


	dt=dto
	postint=postinto
        goto 999
	Eperp=sqrt(Enord**2+Eest**2)
        if (Eperp.gt.75.) then
          dt=dto/10.
        elseif (Eperp.gt.50.) then
          dt=dto/5.
        elseif (Eperp.gt.30.) then
          dt=dto/2.
        endif
999	continue

	if (dt_max.gt.0..and.dt.gt.0.d0) then
	  n=int(max(0.d0,log(dt*dt_max)*log2_1))
          if (n.gt.0) dt=dt*2.d0**(-n)
        endif

	return
	end
