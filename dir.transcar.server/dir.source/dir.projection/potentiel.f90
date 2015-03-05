subroutine potentiel(iyd,tu,kp,lon,lat,Eest,Enord,pot,ddp)
	
implicit none

!f2py intent(in) :: iyd,tu,kp,lon,lat
!f2py intent(out) :: Eest,Enord,pot, ddp

integer,parameter :: npt=10000

integer:: iyd,ndeg,mdeg,ierr,i,j
real*8 :: tu
real   :: kp,ddp
real*8 :: lon,lat
real*8 :: phipot(npt)
real*8 :: latequi,Lmin,Lmax
real*8 :: fit(3),Enord,Eest,pot

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

!Lmin=0.d0
!Lmax=0.d0
call coef_pot(iyd,tu,kp,ndeg,mdeg,phipot,	&
	      Lmin,Lmax,latequi,ddp,ierr)
!call IMM_POT(ndeg,mdeg,phipot,Lmin,Lmax,latequi)
if (Lmin.lt.Lmax) then
  call val_fit(lon,lat,ndeg,mdeg,phipot,Lmin,Lmax,latequi,pot,Eest,Enord)
  ddp=pot
else
  Enord=0.d0
  Eest=0.d0
  pot=0.d0
endif

return
end subroutine potentiel
