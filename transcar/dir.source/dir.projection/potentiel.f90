subroutine potentiel(iyd,tu,kp,lon,lat,Eest,Enord,pot,ddp)
	
implicit none

!f2py intent(out) :: pot, ddp

integer,intent(in) :: iyd
real,intent(in) :: kp
double precision,intent(in) :: tu,lon,lat
double precision,intent(out) :: Eest,Enord,pot
real,intent(out) :: ddp

integer,parameter :: npt=10000
integer :: ndeg,mdeg,ierr,i,j
double precision :: phipot(npt), latequi,Lmin,Lmax, fit(3)

save ierr

!TODO segfault if this interface is not defined 

interface

  subroutine val_fit(lon,lat,ndeg,mdeg,coef_psi,latmin,latmax,latequi,psi,psi_est,psi_nord)
    Integer	::	ndeg,mdeg
    Real*8	::	latmin,latmax,latequi
    Real*8	::	lon,lat,coef_psi(:)
    Real*8	::	psi
    Real*8, optional :: psi_est,psi_nord
 end subroutine val_fit

end interface
!Lmin=0.d0
!Lmax=0.d0

print*,'call coef_pot'
call coef_pot(iyd,tu,kp,ndeg,mdeg,phipot,	&
	      Lmin,Lmax,latequi,ddp,ierr)
!call IMM_POT(ndeg,mdeg,phipot,Lmin,Lmax,latequi)
if (Lmin.lt.Lmax) then
  print*,'call val_fit'
  call val_fit(lon,lat,ndeg,mdeg,phipot,Lmin,Lmax,latequi,pot,Eest,Enord)
  ddp=pot
else
  Enord=0.d0
  Eest=0.d0
  pot=0.d0
endif

end subroutine potentiel
