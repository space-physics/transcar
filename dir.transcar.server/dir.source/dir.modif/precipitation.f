	subroutine precipitation(iyd,tu,kp,tmag,lat,Energie,fluxE)
c	implicit none

        integer npt
        parameter(npt=1000)
        integer iyd,ndeg,mdeg,ierr
	real*8 tu
	real kp,Energie,fluxE
	complex*16 phienerg(npt),phifluxE(npt)
        real*8 tmag,lon,lat,latequi,Lmin,Lmax
	
	save ierr

	call coef_prec(iyd,tu,kp,ndeg,mdeg,phienerg,phifluxE,
     &		       Lmin,Lmax,latequi,ierr)
	
	if (ierr.gt.0) then

          call hardelec(latmag,tmag,ikp,Flux_elec_int,fluxE)
c
c         Les sorties de hardelec sont en log decimal :
c         Flux de particules integre [cm-2.s-1.sr-1]
          Flux_elec_int=10.**Flux_elec_int
c         Flux d'energie integre [keV.cm-2.s-1.sr-1]
          fluxE=10.**fluxE

c
c         Calcule l'energie moyenne en eV
c         On assume que la forme du flux est une maxwellienne, bien que
c         hardy et al,(JGR Mai90, p4229-4248, 1985) ont publie des
c 	  formes de flux dont ils precisent autant que faire se peut
c 	  qu'elles ne sont justement pas maxwelliennes. Helas, ils ne
c 	  donnent pas une loi generale de representation de leurs
c 	  precipitations. Que faire ?
c 	  Pour une maxwellienne F = k E exp(-E/Eo), on a
c 	  Flux de particules = integrale de  k F(E) dE
c 	  Flux d'energie     = integrale de  k E F(E) dE
c 	  Le rapport des 2 donne
c 	      Eo = (Flux d'energie/2*Flux de particules)
          Energie=(1000.*fluxE)/(2.*Flux_elec_int)
          fluxE=fluxE/7.



        else
	
          lon=15.d0*tmag
	  energie=val_fit(lon,lat,ndeg,mdeg,phienerg,
     &		       Lmin,Lmax,latequi)
	  fluxE=val_fit(lon,lat,ndeg,mdeg,phifluxE,
     &		       Lmin,Lmax,latequi)


c	dans le fichier initial l'energie est en keV et on la veut en eV
	  energie=exp(energie)*1000./2.
c	dans le fichier initial le flux en energie est en W/m^2 et on le veut en  keV/(cm^2.s.sr)	
	  fluxE=exp(fluxE)*9.9472e10

c	facteur multiplicatif pour le 930216
c	  fluxE=fluxE*6.2832

	endif
	
	return
	end
	
	
