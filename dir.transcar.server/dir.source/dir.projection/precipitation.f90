subroutine precipitation(iyd,tu,kp,lon,lat,Energie,fluxE)

implicit none

integer,parameter :: npt=10000
integer ::iyd,ndeg,mdeg,ierr,ikp
real*8 :: tu
real :: kp,Energie,fluxE,flux_elec_int,rtmag,rlat
real*8 :: phienerg(npt),phifluxE(npt),psi
real*8 :: tmag,lon,lat,latequi,Lmin,Lmax
	
save ierr

interface

  subroutine val_fit(lon,lat,ndeg,mdeg,coef_psi,latmin,latmax,latequi,psi,psi_est,psi_nord)
    Integer	::	ndeg,mdeg
    Real*8	::	latmin,latmax,latequi
    Real*8	::	lon,lat,coef_psi(:)
    Real*8	::	psi
    Real*8, optional :: psi_est,psi_nord
  end subroutine val_fit

end interface



Lmin=0.d0
Lmax=0.d0	
!call coef_prec(iyd,tu,kp,ndeg,mdeg,phienerg,phifluxE,	&
!	       Lmin,Lmax,latequi,ierr)

!call IMM_PREC(ndeg,mdeg,phienerg,phifluxE,Lmin,Lmax,latequi)
if (Lmin.lt.Lmax) then
  ierr=0
else
  ierr=1
endif
	
if (ierr.gt.0) then
  rtmag=lon/15.d0
  rlat=lat
  ikp=kp
  call hardelec(rlat,rtmag,ikp,Flux_elec_int,fluxE)
!
!         Les sorties de hardelec sont en log decimal :
!         Flux de particules integre [cm-2.s-1.sr-1]
          Flux_elec_int=10.**Flux_elec_int
!         Flux d'energie integre [keV.cm-2.s-1.sr-1]
          fluxE=10.**fluxE

!
!         Calcule l'energie moyenne en eV
!         On assume que la forme du flux est une maxwellienne, bien que
!         hardy et al,(JGR Mai90, p4229-4248, 1985) ont publie des
! 	  formes de flux dont ils precisent autant que faire se peut
! 	  qu'elles ne sont justement pas maxwelliennes. Helas, ils ne
! 	  donnent pas une loi generale de representation de leurs
! 	  precipitations. Que faire ?
! 	  Pour une maxwellienne F = k E exp(-E/Eo), on a
! 	  Flux de particules = integrale de  k F(E) dE
! 	  Flux d'energie     = integrale de  k E F(E) dE
! 	  Le rapport des 2 donne
! 	      Eo = (Flux d'energie/2*Flux de particules)
          Energie=(1000.*fluxE)/(2.*Flux_elec_int)
          fluxE=fluxE/7.


else
	
  call val_fit(lon,lat,ndeg,mdeg,phienerg,Lmin,Lmax,latequi,psi)
  energie=psi
  call val_fit(lon,lat,ndeg,mdeg,phifluxE,Lmin,Lmax,latequi,psi)
  fluxE=psi


!	dans le fichier initial l'energie est en keV et on la veut en eV
!	  energie=exp(energie)*1000./2.
!	dans le fichier initial le flux en energie est en W/m^2 et on le veut en  keV/(cm^2.s.sr)	
!	  fluxE=exp(fluxE)*9.9472e10

!	dans le fichier initial l'energie est en eV et on la veut en eV
	  energie=exp(energie)
!	dans le fichier initial le flux en energie est en mW/m^2 et on le veut en  keV/(cm^2.s.sr)	
	  fluxE=exp(fluxE)*9.9472e7
	if (energie<50.) then
	  energie=50.
	  fluxE=0.
	endif

!	facteur multiplicatif pour le 930216
!	  fluxE=fluxE*6.2832

!	print*,Energie,fluxE
endif
	
return
end subroutine precipitation
