subroutine val_fit(lon,lat,ndeg,mdeg,coef_psi,latmin,latmax,latequi,psi,psi_est,psi_nord)
use comm,only : dp, deg2rad,rad2deg,pi
implicit none


real(dp),intent(in) :: lon,lat,coef_psi(:),latmin,latmax,latequi
integer,intent(in) :: ndeg,mdeg
Real(dp),intent(out), optional :: psi,psi_est,psi_nord
!
Integer,parameter::	npt=100000

Real(dp),parameter ::     coeffc=6._dp,re=6.378_dp,lat_top=89.9_dp
Integer	::	i,j,k,l,rang

Logical :: mode
Real(dp)	::  LM,DL
Real(dp)	::	coef_est,coef_nord
Real(dp)	::	psi0
Real(dp)  ::      alpha,xmin,xequi,coslinf,dxmin,d12,d11,d210,d21,d22,det_1,xmax,coslsup
Real(dp)  ::      anm,bnm,cnm(0:npt),cosphi,sinphi
Real(dp)  ::	A(npt),A_est(npt),A_nord(npt),Pn(0:npt),dPn(0:npt)
Real(dp)  ::	xinf,xfit,xsup,phi,expdx
Real(dp)  ::	dxequi,dxequi2,expjdx,expcdx,expcdx_1,Fnm_inf,Fnm_sup,dFnm_inf,dFnm_sup,sum

mode=present(psi_est)

rang=(ndeg+1)*(2*mdeg+1)
LM=(latmax+latmin)/2.0_dp
DL=2.0_dp/(latmax-latmin)
alpha=-DL*rad2deg


phi=deg2rad*lon

coef_est=-1.0_dp/re/cos(deg2rad*min(lat,lat_top))



if (lat<latmin.and.lat>=latequi) then
  xinf=log(tan(pi/4.0_dp-lat*deg2rad/2.0_dp))

  coef_nord=coef_est

  xmin =log(tan(pi/4.0_dp-latmin*deg2rad/2.0_dp))
  xequi=log(tan(pi/4.0_dp-latequi*deg2rad/2.0_dp))
  coslinf=cos(latmin*deg2rad)
  dxequi=(xinf-xequi)
  dxequi2=dxequi*dxequi
  dxmin =xmin-xequi
  d12=dxmin*dxmin
  d11=xmin*d12
  d210=2.0_dp*dxmin*xmin+d12

  k=0
  do i=0,ndeg
    k=k+1
    call plegendr(i,1.0_dp,Pn(i),dPn(i))
    cnm(i)=-alpha*dPn(i)*coslinf
    A(k)=Pn(i)+cnm(i)*(xinf-xmin)
    if (mode) then
      A_est(k)=0.0_dp
      A_nord(k)=-cnm(i)
    endif
  enddo
  do j=1,mdeg
    d21=d210+j*d11
    d22=2.0_dp*dxmin+j*d12
    det_1=1.0_dp/(d11*d22-d12*d21)
    cosphi=cos(real(j,dp)*phi)
    sinphi=sin(real(j,dp)*phi)
    do i=0,ndeg
      k=k+1
      anm=(d22-cnm(i)*d12)*det_1
      bnm=(-d21+cnm(i)*d11)*det_1
      bnm=anm*xinf+bnm
      Fnm_inf=bnm*dxequi2
      A(k)=Fnm_inf*cosphi
      if (mode) then
        dFnm_inf=2.0_dp*dxequi*bnm+anm*dxequi2+j*Fnm_inf
        A_est(k)=-Fnm_inf*j*sinphi
        A_nord(k)=-dFnm_inf*cosphi
      endif

      k=k+1
      A(k)=Fnm_inf*sinphi
      if (mode) then
        A_est(k)=Fnm_inf*j*cosphi
        A_nord(k)=-dFnm_inf*sinphi
      endif
    enddo
  enddo

elseif(lat>=latmin.and.lat<=latmax) then
  xfit=(LM-lat)*DL
  coef_nord=alpha/re

  k=0
  do i=0,ndeg
    k=k+1
    call plegendr(i,xfit,Pn(i),dPn(i))
    A(k)=Pn(i)
    if (mode) then
      A_est(k)=0.0_dp
      A_nord(k)=-dPn(i)
    endif
  enddo
  do j=1,mdeg
    cosphi=dcos(real(j,dp)*phi)
    sinphi=dsin(real(j,dp)*phi)
    do i=0,ndeg
      k=k+1
      A(k)=Pn(i)*cosphi
      if (mode) then
        A_est(k)=-Pn(i)*j*sinphi
        A_nord(k)=-dPn(i)*cosphi
      endif

      k=k+1
      A(k)=Pn(i)*sinphi
      if (mode) then
        A_est(k)=Pn(i)*j*cosphi
        A_nord(k)=-dPn(i)*sinphi
      endif
    enddo
  enddo
elseif(lat>latmax) then

  xmax =log(tan(pi/4.0_dp-latmax*deg2rad/2.0_dp))

!  if (lat<lat_top) then
    xsup=log(tan(deg2rad*(45.0_dp-min(lat,lat_top)/2.0_dp)))
    expcdx_1=xsup-xmax
    expdx=exp(expcdx_1)
    expjdx=1.0_dp
    if (coeffc/=0.0_dp) then
      expcdx=exp(coeffc*expcdx_1)
      expcdx_1=(expcdx-1.0_dp)/coeffc
    else
      expcdx=1.0_dp
    endif
!  else
!    xsup=xmax
!    expcdx=1.0_dp
!    expcdx_1=0.0_dp
!    expdx=0.0_dp
!    expjdx=1.0_dp
!  endif
  coef_nord=coef_est

  coslsup=cos(latmax*deg2rad)

  k=0
  do i=0,ndeg
    k=k+1
    call plegendr(i,-1.0_dp,Pn(i),dPn(i))
    cnm(i)=-alpha*dPn(i)*coslsup
    A(k)=Pn(i)+cnm(i)*expcdx_1
    if (mode) then
      dFnm_sup=cnm(i)*expcdx
      A_est(k)=0.0_dp
      A_nord(k)=-dFnm_sup
    endif
  enddo
  do j=1,mdeg
    cosphi=cos(real(j,dp)*phi)
    sinphi=sin(real(j,dp)*phi)
    expjdx=expdx*expjdx
    do i=0,ndeg
      k=k+1
      anm=(cnm(i)-j*Pn(i))
      Fnm_sup=(Pn(i)+anm*expcdx_1)*expjdx
      A(k)=Fnm_sup*cosphi
      if (mode) then
        dFnm_sup=(anm*(expcdx+j*expcdx_1)+j*Pn(i))*expjdx
        A_est(k)=-Fnm_sup*j*sinphi
        A_nord(k)=-dFnm_sup*cosphi
      endif

      k=k+1
      A(k)=Fnm_sup*sinphi
      if (mode) then
        A_est(k)=Fnm_sup*j*cosphi
        A_nord(k)=-dFnm_sup*sinphi
      endif
    enddo
  enddo

endif


psi=dot_product(A(1:rang),coef_psi(1:rang))
if (mode) then
  psi_est =coef_est*dot_product(A_est(1:rang),coef_psi(1:rang))
  psi_nord=coef_nord*dot_product(A_nord(1:rang),coef_psi(1:rang))
endif


end subroutine val_fit
