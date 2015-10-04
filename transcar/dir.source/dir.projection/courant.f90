subroutine courant(iyd,tu,kp,lon,lat,Jsup)


implicit none

integer,parameter :: npt=10000
integer :: iyd,ndeg,mdeg,ierr
real*8 :: tu,lon,lat
real :: kp,Jsup
real*8 :: phicourant(npt),psi
real*8 :: latequi,Lmin,Lmax
	
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


	call coef_cour(iyd,tu,kp,ndeg,mdeg,phicourant,	&
      		       Lmin,Lmax,latequi,ierr)
Lmin=0.d0
Lmax=0.d0	
!call IMM_COUR(ndeg,mdeg,phicourant,Lmin,Lmax,latequi)

if (Lmin.lt.Lmax) then
  ierr=0
else
  ierr=1
endif
 
Jsup=0.
if (ierr.eq.0) then
  call val_fit(lon,lat,ndeg,mdeg,phicourant,Lmin,Lmax,latequi,psi)
  Jsup=psi
!	dans le fichier initial le courant est en 10^-7 A/m2
  Jsup=Jsup/10.
endif
	
return
end subroutine courant
	
	
