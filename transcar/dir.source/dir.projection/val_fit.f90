subroutine val_fit(lon,lat,ndeg,mdeg,coef_psi,latmin,latmax,latequi,psi,psi_est,psi_nord)

implicit none

real*8,intent(in) :: lon,lat,coef_psi(:),latmin,latmax,latequi
integer,intent(in) :: ndeg,mdeg
Real*8,intent(out), optional :: psi,psi_est,psi_nord
!
Integer,parameter::	npt=100000
Real*8,parameter ::	pi=3.14159265358979d0,deg2rad=1.745329251994330d-2,	&
			rad2deg=57.295779513082320d0
Real*8,parameter ::     coeffc=6.d0,re=6.378d0,lat_top=89.9d0
Integer	::	i,j,k,l,rang

Logical :: mode
Real*8	::  LM,DL
Real*8	::	coef_est,coef_nord
Real*8	::	psi0
Real*8  ::      alpha,xmin,xequi,coslinf,dxmin,d12,d11,d210,d21,d22,det_1,xmax,coslsup
Real*8  ::      anm,bnm,cnm(0:npt),cosphi,sinphi
Real*8  ::	A(npt),A_est(npt),A_nord(npt),Pn(0:npt),dPn(0:npt)
Real*8  ::	xinf,xfit,xsup,phi,expdx
Real*8  ::	dxequi,dxequi2,expjdx,expcdx,expcdx_1,Fnm_inf,Fnm_sup,dFnm_inf,dFnm_sup,sum

mode=present(psi_est)

rang=(ndeg+1)*(2*mdeg+1)
LM=(latmax+latmin)/2.d0
DL=2.d0/(latmax-latmin)
alpha=-DL*rad2deg


phi=deg2rad*lon

coef_est=-1.d0/re/cos(deg2rad*min(lat,lat_top))



if (lat<latmin.and.lat>=latequi) then
  xinf=log(tan(pi/4.d0-lat*deg2rad/2.d0))

  coef_nord=coef_est

  xmin =log(tan(pi/4.d0-latmin*deg2rad/2.d0))
  xequi=log(tan(pi/4.d0-latequi*deg2rad/2.d0))
  coslinf=cos(latmin*deg2rad)
  dxequi=(xinf-xequi)
  dxequi2=dxequi*dxequi
  dxmin =xmin-xequi
  d12=dxmin*dxmin
  d11=xmin*d12
  d210=2.d0*dxmin*xmin+d12

  k=0
  do i=0,ndeg
    k=k+1
    call plegendr(i,1.d0,Pn(i),dPn(i))
    cnm(i)=-alpha*dPn(i)*coslinf
    A(k)=Pn(i)+cnm(i)*(xinf-xmin)
    if (mode) then
      A_est(k)=0.d0
      A_nord(k)=-cnm(i)
    endif
  enddo
  do j=1,mdeg
    d21=d210+j*d11
    d22=2.d0*dxmin+j*d12
    det_1=1.d0/(d11*d22-d12*d21)
    cosphi=dcos(j*phi)
    sinphi=dsin(j*phi)
    do i=0,ndeg
      k=k+1
      anm=(d22-cnm(i)*d12)*det_1
      bnm=(-d21+cnm(i)*d11)*det_1
      bnm=anm*xinf+bnm
      Fnm_inf=bnm*dxequi2
      A(k)=Fnm_inf*cosphi
      if (mode) then
        dFnm_inf=2.d0*dxequi*bnm+anm*dxequi2+j*Fnm_inf
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
      A_est(k)=0.d0
      A_nord(k)=-dPn(i)
    endif
  enddo
  do j=1,mdeg
    cosphi=dcos(j*phi)
    sinphi=dsin(j*phi)
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

  xmax =log(tan(pi/4.d0-latmax*deg2rad/2.d0))

!  if (lat<lat_top) then
    xsup=log(tan(deg2rad*(45.d0-min(lat,lat_top)/2.d0)))
    expcdx_1=xsup-xmax
    expdx=exp(expcdx_1)
    expjdx=1.d0
    if (coeffc/=0.d0) then
      expcdx=exp(coeffc*expcdx_1)
      expcdx_1=(expcdx-1.d0)/coeffc
    else
      expcdx=1.d0
    endif
!  else
!    xsup=xmax
!    expcdx=1.d0
!    expcdx_1=0.d0
!    expdx=0.d0
!    expjdx=1.d0
!  endif
  coef_nord=coef_est

  coslsup=cos(latmax*deg2rad)

  k=0
  do i=0,ndeg
    k=k+1
    call plegendr(i,-1.d0,Pn(i),dPn(i))
    cnm(i)=-alpha*dPn(i)*coslsup
    A(k)=Pn(i)+cnm(i)*expcdx_1
    if (mode) then
      dFnm_sup=cnm(i)*expcdx
      A_est(k)=0.d0
      A_nord(k)=-dFnm_sup
    endif
  enddo
  do j=1,mdeg
    cosphi=cos(j*phi)
    sinphi=sin(j*phi)
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
