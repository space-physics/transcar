subroutine solve_ODE(nli,y1,y2,a,b,dt)

implicit none

integer		:: nli
real*8		:: y1(nli),y2(nli),a(nli,nli),b(nli,nli)
real		:: dt

integer,parameter :: npt=20
complex*16	:: mat(npt,npt),mat_p1(npt,npt),mat_p2(npt,npt),x1(npt),x2(npt),y11(npt),y21(npt)
real*8		:: fv1(npt),wr(npt),wi(npt),rnorm
complex*16	:: cnorm,lambda,expnu,dexpnu,temp
integer		:: ierr,i,j,k,is1,is2,iv1(npt),ind_permut

interface

  subroutine  balanc(nm,n,a,low,igh,scale)
      integer	:: n,nm,igh,low
      real*8	:: a(nm,n),scale(n)
  end subroutine  balanc

  subroutine balbak(nm,n,low,igh,scale,m,z)
      integer	:: n,m,nm,igh,low
      real*8	:: scale(n),z(nm,m)
  end subroutine balbak

  subroutine  elmhes(nm,n,low,igh,a,int)
      integer	:: n,nm,low,igh
      real*8	:: a(nm,n)
      integer	:: int(igh)
  end subroutine  elmhes

  subroutine  eltran(nm,n,low,igh,a,int,z)
      integer	:: n,nm,igh,low
      real*8	:: a(nm,igh),z(nm,n)
      integer	:: int(igh)
  end subroutine  eltran

  subroutine  hqr2(nm,n,low,igh,h,wr,wi,z,ierr)
      integer	:: n,nm,igh,low,ierr
      real*8	:: h(nm,n),wr(n),wi(n),z(nm,n)
  end subroutine  hqr2


end interface

!write(*,'(6(g12.5,1x))')y1
!write(*,'(6(g12.5,1x))')y2

!do i=1,nli
!  write(*,'(i2,1x,6(g12.5,1x))')i,(a(i,j),j=1,nli)
!enddo



!nli=size(y1,1)

call  balanc(nli,nli,a,is1,is2,fv1)
call  elmhes(nli,nli,is1,is2,a,iv1)
call  eltran(nli,nli,is1,is2,a,iv1,b)
call  hqr2(nli,nli,is1,is2,a,wr,wi,b,ierr)
if (ierr .eq. 0) call  balbak(nli,nli,is1,is2,fv1,nli,b)


do i=1,nli
  y11(i)=dcmplx(y1(i),0.d0)
  y21(i)=dcmplx(y2(i),0.d0)
  j=1
  do while (j<=nli)
    if (wi(j)==0.d0) then
      mat(i,j)=dcmplx(b(i,j),0.d0)
      j=j+1
    else
      mat(i,j)=dcmplx(b(i,j),b(i,j+1))
      mat(i,j+1)=dcmplx(b(i,j),-b(i,j+1))
      j=j+2
    endif
  enddo
enddo
do j=1,nli
  rnorm=dot_product(mat(1:nli,j),mat(1:nli,j))
  mat(1:nli,j)=mat(1:nli,j)/dsqrt(rnorm)
enddo
!do i=1,nli
!  write(*,'(i2,1x,6(g12.5,1x))')i,(real(mat(i,j)),j=1,nli)
!enddo

do j=1,nli
  lambda=dcmplx(wr(j),wi(j))
!  print*,'--->',j,lambda,lambda*dt,dt
  expnu=exp(lambda*dt)
  if (abs(lambda*dt)<1.d-16) then
    dexpnu=dt
  else
    dexpnu=(expnu-1.d0)/lambda
  endif
  do i=1,nli
    mat_p1(i,j)=mat(i,j)*expnu
    mat_p2(i,j)=mat(i,j)*dexpnu
  enddo
!  print*,'--->',j,expnu,dexpnu,lambda
enddo
!do i=1,nli
!  write(*,'(i2,1x,6(g12.5,1x))')i,(real(mat_p1(i,j)),j=1,nli)
!  write(*,'(i2,1x,6(g12.5,1x))')i,(real(mat_p2(i,j)),j=1,nli)
!enddo

do i=1,nli
  rnorm=abs(mat(i,i))
  ind_permut=i
  do k=i+1,nli
    if (abs(mat(k,i))>rnorm) then
      ind_permut=k
      rnorm=abs(mat(k,i))
    endif
  enddo
  if (ind_permut>i) then
    do j=i,nli
      temp=mat(i,j)
      mat(i,j)=mat(ind_permut,j)
      mat(ind_permut,j)=temp
    enddo
    temp=y11(i)
    y11(i)=y11(ind_permut)
    y11(ind_permut)=temp
    temp=y21(i)
    y21(i)=y21(ind_permut)
    y21(ind_permut)=temp
  endif
  cnorm=1.d0/mat(i,i)
  mat(i,i)=1.d0
  y11(i)=y11(i)*cnorm
  y21(i)=y21(i)*cnorm
  do j=i+1,nli
    mat(i,j)=mat(i,j)*cnorm
  enddo
  do k=i+1,nli
    cnorm=mat(k,i)
    if (cnorm.ne.0.d0) then
      mat(k,i)=0.d0
      y11(k)=y11(k)-cnorm*y11(i)
      y21(k)=y21(k)-cnorm*y21(i)
      do j=i+1,nli
        mat(k,j)=mat(k,j)-cnorm*mat(i,j)
      enddo
    endif
  enddo
enddo


x1(nli)=y11(nli)
x2(nli)=y21(nli)
do i=nli-1,1,-1
  x1(i)=y11(i)
  x2(i)=y21(i)
  do j=i+1,nli
    x1(i)=x1(i)-mat(i,j)*x1(j)
    x2(i)=x2(i)-mat(i,j)*x2(j)
  enddo
enddo

do i=1,nli
  cnorm=0.d0
  do j=1,nli
    cnorm=cnorm+mat_p1(i,j)*x1(j)+mat_p2(i,j)*x2(j)
  enddo
  y1(i)=cnorm
enddo
!write(*,'(6(g12.5,1x))')y1

end subroutine solve_ODE
