subroutine potentiel(iyd,tu,kp,lon,lat,Eest,Enord,pot,ddp)
use comm, only: dp,tic
implicit none


integer,intent(in) :: iyd
real,intent(in) :: kp
real(dp),intent(in) :: tu,lon,lat
real(dp),intent(out) :: Eest,Enord,pot
real,intent(out) :: ddp

integer,parameter :: npt=10000
integer :: ndeg,mdeg,i,j
real(dp) :: phipot(npt), latequi,Lmin,Lmax, fit(3)

!TODO segfault if this interface is not defined 

interface

  subroutine val_fit(lon,lat,ndeg,mdeg,coef_psi,latmin,latmax,latequi,psi,psi_est,psi_nord)
    import dp
    Real(dp),intent(in)	::	lon,lat,coef_psi(:)
    Integer,intent(in)	::	ndeg,mdeg
    Real(dp)	::	latmin,latmax,latequi

    Real(dp)	::	psi
    Real(dp), optional :: psi_est,psi_nord
 end subroutine val_fit

end interface
!Lmin=0.d0
!Lmax=0.d0

call cpu_time(tic)
print *,tic,'potentiel.f: call coef_pot'
call coef_pot(iyd,tu,kp,&
	                     ndeg,mdeg,phipot,Lmin,Lmax,latequi,ddp)
!call IMM_POT(ndeg,mdeg,phipot,Lmin,Lmax,latequi)
if (Lmin < Lmax) then
  call cpu_time(tic)
  print *,tic,'potentiel.f: call val_fit:  lon,lat  ',lon,lat
  call val_fit(lon,lat,ndeg,mdeg,phipot,Lmin,Lmax,latequi,pot,Eest,Enord)
  ddp=pot
else
  Enord=0.d0
  Eest=0.d0
  pot=0.d0
endif

end subroutine potentiel
