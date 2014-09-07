Subroutine stabenerg(npoint,dens,told,tnew,qold,qnew,D7e,D7q,Cji,dt)

implicit none

Integer,parameter	:: npt=500
Integer			:: npoint,i
Real			:: Cji,dt
Real,parameter		:: omega=.9
Real, dimension(npoint)	:: dens,told,tnew,qold,qnew,d7e,d7q
Real, dimension(npt)	:: a1,a2,d1,d2
Real, dimension(npt)	:: x1,x2,y1,y2,cf

Real, dimension(npt)	::LO,LN(NPT),AH,RLN,LH,RLH,ROH,RNH,ADUGTH
Common  /FCT_GRID/ LO, LN, AH, RLN, LH, RLH, ROH, RNH, ADUGTH

cf(1:npoint)=Cji*dt*omega*RLN(1:npoint)


y1(1:npoint:2)=(tnew(1:npoint:2)-told(1:npoint:2))*dens(1:npoint:2)*omega
!y1(1:npoint:2)=(tnew(1:npoint:2)-told(1:npoint:2))*omega
y1(2:npoint:2)=(qnew(2:npoint:2)-qold(2:npoint:2))*omega

d1(1:npoint:2)=1.-dt*D7e(1:npoint:2)
d1(2:npoint:2)=1.-dt*D7q(2:npoint:2)
a1(1:npoint:2)=cf(1:npoint:2)/3.
a1(2:npoint:2)=1.25*told(2:npoint:2)*cf(2:npoint:2)
!a1(1:npoint:2)=cf(1:npoint:2)/3./dens(1:npoint:2)
!a1(2:npoint:2)=1.25*told(2:npoint:2)*dens(2:npoint:2)*cf(2:npoint:2)

y2(1:npoint:2)=(qnew(1:npoint:2)-qold(1:npoint:2))*omega
y2(2:npoint:2)=(tnew(2:npoint:2)-told(2:npoint:2))*dens(2:npoint:2)*omega
!y2(2:npoint:2)=(tnew(2:npoint:2)-told(2:npoint:2))*omega

d2(1:npoint:2)=1.-dt*D7q(1:npoint:2)
d2(2:npoint:2)=1.-dt*D7e(2:npoint:2)
a2(1:npoint:2)=1.25*told(1:npoint:2)*cf(1:npoint:2)
a2(2:npoint:2)=cf(2:npoint:2)/3.
!a2(1:npoint:2)=1.25*told(1:npoint:2)*dens(1:npoint:2)*cf(1:npoint:2)
!a2(2:npoint:2)=cf(2:npoint:2)/3./dens(2:npoint:2)

y1(2:npoint-1)=d1(2:npoint-1)*y1(2:npoint-1)
y2(2:npoint-1)=d2(2:npoint-1)*y2(2:npoint-1)

d1(1)=1.
a1(1)=0.

d2(1)=1.
a2(1)=0.

do i=2,npoint-1
  d1(i)=d1(i)+a1(i)*a1(i-1)
  y1(i)  =(y1(i)  +a1(i)*y1(i-1))/d1(i)
  a1(i)=a1(i)/d1(i)

  d2(i)=d2(i)+a2(i)*a2(i-1)
  y2(i)  =(y2(i)  +a2(i)*y2(i-1))/d2(i)
  a2(i)=a2(i)/d2(i)
enddo

x1(1:npoint)=y1(1:npoint)
x2(1:npoint)=y2(1:npoint)

do i=npoint-1,2,-1
  x1(i)=y1(i)-a1(i)*x1(i+1)
  x2(i)=y2(i)-a2(i)*x2(i+1)
enddo

y2=x2
x2(2:npoint:2)=x1(2:npoint:2)
x1(2:npoint:2)=y2(2:npoint:2)


tnew(2:npoint-1)=told(2:npoint-1)+x1(2:npoint-1)/dens(2:npoint-1)/omega
!tnew(2:npoint-1)=told(2:npoint-1)+x1(2:npoint-1)/omega
qnew(2:npoint-1)=qold(2:npoint-1)+x2(2:npoint-1)/omega

return

End subroutine stabenerg
