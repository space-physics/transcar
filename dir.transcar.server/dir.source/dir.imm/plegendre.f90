Subroutine polyleg(x,pn,dpn)

implicit none

integer,parameter :: npt=100000
integer	:: n,i,nb_point,ndeg
real*8  :: x(:),pn(:,0:),dpn(:,0:),xm
real*8	:: x2(npt),x3(npt),x4(npt),x5(npt)
real*8	:: pn_1(npt),pn_2(npt),dpn_1(npt),dpn_2(npt)


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
        dpn_2(1:nb_point)=dpn_1(1:nb_point,n)
        pn_1 (1:nb_point)=pn   (1:nb_point)
        dpn_1(1:nb_point)=dpn  (1:nb_point,n)
      enddo
  end select
enddo

end subroutine polyleg
