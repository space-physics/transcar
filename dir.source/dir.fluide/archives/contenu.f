	subroutine contenu(nx,iyd,sec,glat,glong,stl,f107a,f107,
     &				ap,z,rhot)

	parameter (npt=500,Ngauss=6)

	real z(npt),zlim,zlim_1,zi,rho,rhot(npt)
	real sec,glat,glong,stl,f107a,f107,ap(7),d(8),t(2)
	real amu,mn(3),hauteur
	integer index(3)

	real Nh(npt),No(npt),No2(npt),Nn2(npt),Tn(npt)
	real Un(npt),Vn(npt),Wn(npt)

	common /neutral/ Nh,No,No2,Nn2,Tn,Un,Vn,Wn

	common/limite/indlim,indlim_1,zlim,zlim_1

	real wi(Ngauss),ci(Ngauss)


	data wi/ .171324492379170, .360761573048139, .467913934572691,
     &           .467913934572691, .360761573048139, .171324492379170/
	data ci/-.932469514203152,-.661209386466265,-.238619186083197,
     &           .238619186083197, .661209386466265, .932469514203152/

	data amu/1.667e-27/
	data mn/16.,32.,28./
	data index/2,4,3/

        real Re,kb,g0
        data Re,g0/6.378e6,9.80665/
        data kb/1.38e-23/

	hauteur(y1,y2,x1,x2)=(x1-x2)/(log(y2)-log(y1))

c	contenu entre zlim et l'infini

	  HO =hauteur(No (nx-1),No (nx),z(nx-1),z(nx))
	  HO2=hauteur(No2(nx-1),No2(nx),z(nx-1),z(nx))
	  HN2=hauteur(Nn2(nx-1),Nn2(nx),z(nx-1),z(nx))

	  rhoH=HO *mn(1)*No (nx)
     &	      +HO2*mn(2)*No2(nx)
     &	      +HN2*mn(3)*Nn2(nx)
	  rhoH=rhoH*amu*1.e8

c	contenu entre z et z(nx)

	do ialt=1,nx

	  rhot(ialt)=rhoH

	  alpha=(z(nx)-z(ialt))/2.
	  beta =(z(nx)+z(ialt))/2.
	  do i=1,Ngauss
	    zi=alpha*ci(i)+beta
	    if (zi.le.z(indlim)) then
	      call gtd6(iyd,sec,zi,glat,glong,stl,f107a,f107,ap,48,d,t)
	    else
	      fact=exp(-amu*g0*Re**2/kb/Tn(indlim)*(zi-z(indlim))
     &             /(zi+Re)/(z(indlim)*1000.+Re))
              d(2)=No(indlim)*fact**mn(1)
              d(4)=No2(indlim)*fact**mn(2)
              d(3)=Nn2(indlim)*fact**mn(3)
            endif

	    rho=0.
	    do j=1,3
	      rho=rho+mn(j)*d(index(j))
	    enddo
	    rho=rho*amu*1.e8
	    rhot(ialt) =rhot(ialt) +alpha*wi(i)*rho
	  enddo

	enddo

	return
	end
