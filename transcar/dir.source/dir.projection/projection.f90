subroutine geometrie(nb_point,nb_data,ind_data,lon,lat,latmin,latmax,latequi,mode,ndeg,mdeg,rang,A,B)

!    mode=0 on prend toutes les données entre latequi et 90°
!    mode=1 on prend toutes les données entre latequi et latmax
!     mode=2 on prend toutes les données entre latmin et 90°
!    mode=3 on prend toutes les données entre latmin et latmax


implicit none

Integer,parameter::    npt=100000
Real*8,parameter ::    pi=3.14159265358979d0,deg2rad=1.745329251994330d-2,    &
            rad2deg=57.295779513082320d0
Real*8,parameter ::     coeffc=6.d0
Integer    ::    i,j,k,l,longueur,sta_inf,end_inf,sta_fit,end_fit,sta_sup,end_sup

Integer    ::    nb_point,nb_data,ind_data(npt),mode,ndeg,mdeg,rang,nb_inf,nb_fit,nb_sup
Real*8    ::    latmin,latmax,latequi,LM,DL
Real*8    ::    lon(nb_point),lat(nb_point)
Real*8  ::      alpha,xmin,xequi,coslinf,dxmin,d12,d11,d210,d21,d22,det_1,xmax,coslsup
Real*8  ::      anm,bnm,cnm
Real*8  ::  A(nb_point*rang),B(rang*nb_point),x1(1)=1.d0,x_1(1)=-1.d0
Real*8,allocatable  :: xinf(:),xfit(:),xsup(:),phi(:),Pleg(:,:),dPleg(:,:),C(:,:),S(:,:)
Real*8,allocatable  :: dxequi2(:),expcdx_1(:),Pinf(:,:),dPinf(:,:),Psup(:,:),dPsup(:,:),Fnm_inf(:),Fnm_sup(:)
Real*8,allocatable  :: At(:,:),mat(:,:),mat_1(:,:),mat1(:,:),A1(:,:),B1(:,:)
Integer,allocatable :: ind_pos(:)
Logical,allocatable :: flg_inf(:),flg_fit(:),flg_sup(:)

interface

  subroutine polyleg(x,P,dP)
    Real*8  :: x(:),P(:,:),dP(:,:)
  end subroutine polyleg

  subroutine solve_sys(y,mat,x)
    Real*8  ::      y(:),x(:)
    Real*8  ::      mat(:,:)
  end subroutine solve_sys

  subroutine inverse(mat,mat_1)
    Real*8  ::      mat(:,:),mat_1(:,:)
  end subroutine inverse

  subroutine coefficient(ndeg,mdeg,nb_data,ind_data,B,pot,coef)
!  subroutine coefficient(ndeg,mdeg,nb_data,ind_data,B,pot,coef_r,coef_i)
    Integer ::      ndeg,mdeg,nb_data,ind_data(:)
    Real*8  ::      B(:,:)
    Real*8  ::      pot(:)
    Real*8  ::        coef(:)
!    Real*8  ::        coef_r(:),coef_i(:)
  end subroutine coefficient

  subroutine projection(nb_data,ind_data,rang,A,B,pot,pot_fit)
    Integer :: nb_data,ind_data(:),rang
    Real*8  ::      A(:,:),B(:,:)
    Real*8  ::      pot(:),pot_fit(:)
  end subroutine projection

end interface


LM=(latmax+latmin)/2.d0
DL=2.d0/(latmax-latmin)
alpha=-DL*rad2deg


allocate(ind_pos(nb_point))
allocate(flg_inf(nb_point))
allocate(flg_sup(nb_point))
allocate(flg_fit(nb_point))

ind_pos=(/(i,i=1,nb_point)/)
flg_inf=.false.
flg_fit=.false.
flg_sup=.false.

select case (mode)
  case (1)
    where (lat<latmin.and.lat>=latequi)
      flg_inf=.true.
    endwhere
  case (2)
    where(lat>latmax.and.lat<90.d0)
      flg_sup=.true.
    endwhere
  case (3)

  case default
    where (lat<latmin.and.lat>=latequi)
      flg_inf=.true.
    endwhere
    where(lat>latmax.and.lat<90.d0)
      flg_sup=.true.
    endwhere
endselect

where(lat>=latmin.and.lat<=latmax)
  flg_fit=.true.
endwhere

nb_inf=count(flg_inf)
nb_fit=count(flg_fit)
nb_sup=count(flg_sup)


nb_data=nb_inf+nb_sup+nb_fit

if (nb_inf>0) then
  allocate(xinf(nb_inf))
  allocate(dxequi2(nb_inf))
  allocate(Fnm_inf(nb_inf))
  allocate(Pinf(1,0:ndeg))
  allocate(dPinf(1,0:ndeg))

  sta_inf=1
  end_inf=nb_inf
  ind_data(sta_inf:end_inf)=pack(ind_pos,flg_inf)

  xinf=log(tan(pi/4.d0-lat(ind_data(sta_inf:end_inf))*deg2rad/2.d0))
  xmin =log(tan(pi/4.d0-latmin*deg2rad/2.d0))
  xequi=log(tan(pi/4.d0-latequi*deg2rad/2.d0))
  coslinf=cos(latmin*deg2rad)
  dxequi2=(xinf-xequi)**2
  dxmin =xmin-xequi
  d12=dxmin*dxmin
  d11=xmin*d12
  d210=2.d0*dxmin*xmin+d12
  call polyleg(x1,Pinf,dPinf)

endif

if (nb_fit>0) then
  allocate(xfit(nb_fit))
  allocate(Pleg(nb_fit,0:ndeg))
  allocate(dPleg(nb_fit,0:ndeg))

  sta_fit=nb_inf+1
  end_fit=sta_fit-1+nb_fit
  ind_data(sta_fit:end_fit)=pack(ind_pos,flg_fit)

  xfit=(LM-lat(ind_data(sta_fit:end_fit)))*DL
  call polyleg(xfit,Pleg,dPleg)
endif

if (nb_sup>0) then
  allocate(xsup(nb_sup))
  allocate(expcdx_1(nb_sup))
  allocate(Fnm_sup(nb_sup))
  allocate(Psup(1,0:ndeg))
  allocate(dPsup(1,0:ndeg))

  sta_sup=nb_inf+nb_fit+1
  end_sup=sta_sup-1+nb_sup
  ind_data(sta_sup:end_sup)=pack(ind_pos,flg_sup)

  xsup=log(tan(pi/4.d0-lat(ind_data(sta_sup:end_sup))*deg2rad/2.d0))
  xmax =log(tan(pi/4.d0-latmax*deg2rad/2.d0))
  coslsup=cos(latmax*deg2rad)
  expcdx_1=xsup-xmax
  if (coeffc/=0.d0) expcdx_1=(exp(coeffc*expcdx_1)-1.d0)/coeffc
  call polyleg(x_1,Psup,dPsup)
endif

allocate(phi(nb_data))
allocate(C(nb_data,mdeg))
allocate(S(nb_data,mdeg))

phi=deg2rad*lon(ind_data(1:nb_data))
do j=1,mdeg
  C(:,j)=cos(j*phi)
  S(:,j)=sin(j*phi)
enddo


allocate(A1(nb_data,rang))
allocate(At(rang,nb_data))
allocate(B1(rang,nb_data))
allocate(mat(rang,rang))
allocate(mat_1(rang,rang))
allocate(mat1(rang,rang))


if (nb_fit>0) A1(sta_fit:end_fit,1:ndeg+1)=Pleg(:,0:ndeg)

k=0
do i=0,ndeg
  k=k+1
  if (nb_inf>0) A1(sta_inf:end_inf,k)=Pinf(1,i)-alpha*dPinf(1,i)*coslinf*(xinf-xmin)
  if (nb_sup>0) A1(sta_sup:end_sup,k)=Psup(1,i)-alpha*dPsup(1,i)*coslsup*expcdx_1
enddo

do j=1,mdeg
  if (nb_inf>0) then
    d21=d210+j*d11
    d22=2.d0*dxmin+j*d12
    det_1=1.d0/(d11*d22-d12*d21)
  endif
  do i=0,ndeg
    k=k+1
    if (nb_inf>0) then
      cnm=-alpha*dPinf(1,i)*coslinf
      anm=(d22-cnm*d12)*det_1
      bnm=(-d21+cnm*d11)*det_1
      Fnm_inf=(anm*xinf+bnm)*dxequi2
      A1(sta_inf:end_inf,k)=Fnm_inf*C(sta_inf:end_inf,j)
    endif

    if (nb_fit>0) A1(sta_fit:end_fit,k)=Pleg(:,i)*C(sta_fit:end_fit,j)

    if (nb_sup>0) then
      Fnm_sup=(Psup(1,i)-(alpha*dPsup(1,i)*coslsup+j*Psup(1,i))*expcdx_1)*exp(j*(xsup-xmax))
      A1(sta_sup:end_sup,k)=Fnm_sup*C(sta_sup:end_sup,j)
    endif

    k=k+1
    if (nb_inf>0) A1(sta_inf:end_inf,k)=Fnm_inf*S(sta_inf:end_inf,j)
    if (nb_fit>0) A1(sta_fit:end_fit,k)=Pleg(:,i)*S(sta_fit:end_fit,j)
    if (nb_sup>0) A1(sta_sup:end_sup,k)=Fnm_sup*S(sta_sup:end_sup,j)
  enddo
enddo


At=transpose(A1)
mat=matmul(At,A1)
mat1=mat
call inverse(mat,mat_1)
mat=matmul(mat1,mat_1)
B1=matmul(mat_1,At)




longueur=rang*nb_data
A(1:longueur)=transfer(A1,A1(1,1),longueur)
B(1:longueur)=transfer(B1,B1(1,1),longueur)

deallocate(ind_pos)
deallocate(At)
deallocate(mat)
deallocate(mat_1)
deallocate(mat1)
deallocate(A1)
deallocate(B1)

deallocate(phi)
deallocate(C)
deallocate(S)

if (nb_inf>0) then
  deallocate(flg_inf)
  deallocate(xinf)
  deallocate(Pinf)
  deallocate(dPinf)
  deallocate(Fnm_inf)
  deallocate(dxequi2)
endif

if (nb_fit>0) then
  deallocate(flg_fit)
  deallocate(xfit)
  deallocate(Pleg)
  deallocate(dPleg)
endif

if (nb_sup>0) then
  deallocate(flg_sup)
  deallocate(xsup)
  deallocate(expcdx_1)
  deallocate(Psup)
  deallocate(dPsup)
  deallocate(Fnm_sup)
endif


end subroutine geometrie

Subroutine polyleg(x,pn,dpn)

implicit none

integer,parameter :: npt=100000
integer    :: n,i,nb_point,ndeg
real*8  :: x(:),pn(:,0:),dpn(:,0:),xm
real*8    :: x2(npt),x3(npt),x4(npt),x5(npt)
real*8    :: pn_1(npt),pn_2(npt),dpn_1(npt),dpn_2(npt)


ndeg=size(pn,2)-1
nb_point=size(x)
do n=0,ndeg
  select case (n)
    case(0)
      pn (1:nb_point,n)= 1.d0
      dpn(1:nb_point,n)= 0.d0
    case(1)
      pn (1:nb_point,n) = x(:)
      dpn(1:nb_point,n) = 1.d0
    case(2)
      x2(1:nb_point)=x(:)*x(:)
      pn (1:nb_point,n)=1.5d0*x2(1:nb_point)-.5d0
      dpn(1:nb_point,n)=3.d0*x(:)
    case(3)
      x2(1:nb_point)=x(:)*x(:)
      x3(1:nb_point)=x2(1:nb_point)*x(:)
      pn(1:nb_point,n)=2.5d0*x3(1:nb_point)-1.5d0*x(:)
      dpn(1:nb_point,n)=7.5d0*x2(1:nb_point)-1.5d0
    case(4)
      x2(1:nb_point)=x(:)*x(:)
      x3(1:nb_point)=x2(1:nb_point)*x(:)
      x4(1:nb_point)=x3(1:nb_point)*x(:)
      pn (1:nb_point,n)=4.375d0*x4(1:nb_point)-3.75d0*x2(1:nb_point)+.375d0
      dpn(1:nb_point,n)=17.5d0*x3(1:nb_point)-7.5d0*x(:)
    case(5)
      x2(1:nb_point)=x(:)*x(:)
      x3(1:nb_point)=x2(1:nb_point)*x(:)
      x4(1:nb_point)=x3(1:nb_point)*x(:)
      x5(1:nb_point)=x4(1:nb_point)*x(:)
      pn (1:nb_point,n)=7.875d0*x5(1:nb_point) - 8.75d0*x3(1:nb_point)+1.875d0*x(:)
      dpn(1:nb_point,n)=39.375d0*x4(1:nb_point)-26.25d0*x2(1:nb_point)+1.875d0
    case default
      x2(1:nb_point)=x(:)*x(:)
      x3(1:nb_point)=x2(1:nb_point)*x(:)
      x4(1:nb_point)=x3(1:nb_point)*x(:)
      x5(1:nb_point)=x4(1:nb_point)*x(:)
      pn_2 (1:nb_point)=4.375d0*x4(1:nb_point)-3.75d0*x2(1:nb_point)+.375d0
      dpn_2(1:nb_point)=17.5d0*x3(1:nb_point)-7.5d0*x(:)
      pn_1 (1:nb_point)=7.875d0*x5(1:nb_point) - 8.75d0*x3(1:nb_point)+1.875d0*x(:)
      dpn_1(1:nb_point)=39.375d0*x4(1:nb_point)-26.25d0*x2(1:nb_point)+1.875d0
      xm=5.d0
      do i=1,n-5
        xm=xm+1.d0
        pn (1:nb_point,n)=(2.d0-1.d0/xm)*x(:)*pn_1(1:nb_point)-(1.d0-1.d0/xm)*pn_2(1:nb_point)
        dpn(1:nb_point,n)=(2.d0-1.d0/xm)*(pn_1(1:nb_point)+x(:)*dpn_1(1:nb_point))-(1.d0-1.d0/xm)*dpn_2(1:nb_point)
        pn_2 (1:nb_point)=pn_1 (1:nb_point)
        dpn_2(1:nb_point)=dpn_1(1:nb_point)
        pn_1 (1:nb_point)=pn   (1:nb_point,n)
        dpn_1(1:nb_point)=dpn  (1:nb_point,n)
      enddo
  end select
enddo

end subroutine polyleg

subroutine solve_sys(y,mat,x)

implicit none

Integer ::      i,j,k,l,rang

Real*8  ::      y(:),x(:)
Real*8  ::      mat(:,:)
Real*8  ::      xnorm

rang=size(y)

do i=1,rang
  xnorm=1.d0/mat(i,i)
  mat(i,i)=1.d0
  y(i)=y(i)*xnorm
  mat(i,i+1:rang)=mat(i,i+1:rang)*xnorm
  do k=i+1,rang
    xnorm=mat(k,i)
    mat(k,i)=0.d0
    mat(k,i+1:rang)=mat(k,i+1:rang)-mat(i,i+1:rang)*xnorm
    y(k)=y(k)-xnorm*y(i)
  enddo
enddo

x(rang)=y(rang)/mat(rang,rang)

do i=rang-1,1,-1
  x(i)=y(i)-dot_product(mat(i,i+1:rang),x(i+1:rang))
  x(i)=x(i)/mat(i,i)
enddo

end subroutine solve_sys

subroutine inverse(mat,mat_1)

implicit none

Integer,parameter :: npt=100000
Integer ::      i,j,k,l,rang,ind(npt),signe

Real*8  ::      mat(:,:),mat_1(:,:)
Real*8  ::      xnorm,temp(npt)

rang=size(mat,1)
mat_1=0.d0
do i=1,rang
  mat_1(i,i)=1.d0
enddo

do i=1,rang
  ind(1:rang-i+1)=maxloc(abs(mat(i:rang,i)))
  k=i-1+ind(1)
  signe=(-1)**(i+k)
  temp(1:rang)=mat(i,:)
  mat(i,:)=signe*mat(k,:)
  mat(k,:)=signe*temp(1:rang)
  temp(1:rang)=mat_1(i,:)
  mat_1(i,:)=signe*mat_1(k,:)
  mat_1(k,:)=signe*temp(1:rang)
  xnorm=1.d0/mat(i,i)
  mat(i,i)=1.d0
  mat(i,i+1:rang)=mat(i,i+1:rang)*xnorm
  mat_1(i,:)=mat_1(i,:)*xnorm
  do k=1,rang
    if (k/=i) then
      xnorm=mat(k,i)
      mat(k,i)=0.d0
      mat(k,i+1:rang)=mat(k,i+1:rang)-mat(i,i+1:rang)*xnorm
      mat_1(k,:)=mat_1(k,:)-mat_1(i,:)*xnorm
    endif
  enddo
enddo

end subroutine inverse

subroutine coefficient(ndeg,mdeg,rang,nb_point,nb_data,ind_data,B,pot,coef)

implicit none

Integer,parameter :: npt=100000
Integer :: ndeg,mdeg,rang,nb_point,nb_data,ind_data(nb_data)
Integer :: i,j,k,l

Real*8  ::      B(rang,nb_data)
Real*8  ::      pot(npt)
Real*8  ::      coef(rang)

coef=matmul(B,pot(ind_data))

end subroutine coefficient


subroutine projection(nb_data,ind_data,rang,A,B,pot,pot_fit)

implicit none

Integer,parameter :: npt=100000
Integer :: nb_data,ind_data(nb_data),rang,i


Real*8  ::      A(nb_data,rang),B(rang,nb_data)
Real*8  ::      pot(npt),pot_fit(npt)
Real*8  ::      coef(npt)


coef(1:rang)=matmul(B,pot(ind_data(1:nb_data)))
pot_fit(ind_data(1:nb_data))=matmul(A,coef(1:rang))



end subroutine projection

