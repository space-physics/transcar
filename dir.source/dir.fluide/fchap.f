      subroutine fchap(iyd,tu,z,lat,lon,f107,ap,chi,chapesp)
        
        implicit none

        integer,parameter :: nlaguer=15 , ngauss=8 , nmax=31,nesp=7
        real,intent(in) :: iyd,tu,z,lat,lon,f107(3),ap(7),chi
        real,intent(out) :: chapesp(nesp)

        
        real zi(nmax),lati(nmax),loni(nmax),wi(nmax),decli
        real zj(nlaguer),wj(nlaguer)
        integer ni,nj,i,j
        real dn(8),tn(2),sec,stl,chap1(nesp),chap2(nesp),indesp(nesp)
            

        data indesp/3,4,2,7,1,8,5/
        common /decl/decli
        
        stl=tu + lon/15. !wow big mistake in prior code due to no implicit none
        call grille(z,lat,lon,tu,chi,decli,zi,lati,loni,wi,ni,zj,wj,nj)
    
        if (ni.gt.0) then
          do i=1,nesp
            chap1(i)=0.
          enddo
           do i = 1,ni
            call gtd6(iyd,tu,zi(i),lati(i),loni(i),stl,
     &            f107(3),f107(2),ap,48,dn,tn)
            do j=1,nesp
               chap1(j)=chap1(j)+wi(i)*dn(indesp(j))
            enddo
          enddo
          do i=1,nesp
            chap2(i)=0.
          enddo
           do i = 1,nj
            call gtd6(iyd,tu,zj(i),lat,lon,stl,
     &              f107(3),f107(2),ap,48,dn,tn)
            do j=1,nesp
               chap2(j)=chap2(j)+wj(i)*dn(indesp(j))
            enddo
          enddo
          do i=1,nesp
            chapesp(i)=chap1(i)/chap2(i)
          enddo
        else
          do i=1,nesp
            chapesp(i)=1.e30
          enddo
        endif

        end subroutine fchap



      subroutine grille(z,lat,lon,tu,chi,decli,
     &                   zi,lati,loni,wi,ni,zj,wj,nj)
      use comm, only: dp,pi,deg2rad,rad2deg
      implicit none

     
      integer,parameter :: nlaguer=15 , ngauss=16 , nmax=31
     
      real,intent(in) :: z,lat,lon,tu,chi,decli
      real,intent(out) :: zi(nmax),lati(nmax),loni(nmax), wi(nmax),
     &                 zj(nlaguer),wj(nlaguer)
      integer,intent(out) :: ni,nj
     
      
      real latref,lonref
      real(dp) :: xl(nlaguer),wl(nlaguer)
      real absc(nmax),cfi,cki,clat,cosdec,dl,h,hver,rlim,rmax,rmin,rz,
     &  sfi,sindec,ski,slat,zmax,zmin
      integer i
      
      
    
      real,parameter :: pis2=pi/2, re=6356.766
    
c    abscisses et poids de Gauss ( a 16 points )
      real, parameter:: xg(*)=
     &            [-0.989400934991649932596,-0.944575023073232576078,
     &            -0.865631202387831743880,-0.755404408355003033895,
     &          -0.617876244402643748447,-0.458016777657227386342,
     &        -0.281603550779258913230,-0.095012509837637440185,
     &             0.095012509837637440185, 0.281603550779258913230,
     &           0.458016777657227386342, 0.617876244402643748447,
     &             0.755404408355003033895, 0.865631202387831743880,
     &         0.944575023073232576078, 0.989400934991649932596]
     
      real, parameter :: wg(*)=
     &        [0.027152459411754094852, 0.062253523938647892863,
     &         0.095158511682492784810, 0.124628971255533872052,
     &          0.149595988816576732081, 0.169156519395002538189,
     &             0.182603415044923588867, 0.189450610455068496285,
     &             0.189450610455068496285, 0.182603415044923588867,
     &          0.169156519395002538189, 0.149595988816576732081,
     &         0.124628971255533872052, 0.095158511682492784810,
     &         0.062253523938647892863, 0.027152459411754094852]

c    abscisses et poids de Laguerre ( a 15 points )



      data xl/  .093307812017,  .492691740302, 1.215595412071,
     &         2.269949526204, 3.667622721751, 5.425336627414,
     &         7.565916226613,10.120228568019,13.130282482176,
     &        16.654407708330,20.776478899449,25.623894226729,
     &        31.407519169754,38.530683306486,48.026085572686/

      data wl/  .239578170311,  .560100842793,  .887008262919,
     &         1.22366440215 , 1.57444872163 , 1.94475197653 ,
     &         2.34150205661 , 2.77404192683 , 3.25564334640 ,
     &         3.80631171423 , 4.45847775384 , 5.27001778443 ,
     &         6.35956346973 , 8.03178763212 ,11.5277721009  /
        latref=lat*deg2rad
        lonref=(lon+15.*tu)*deg2rad
        zmax=800.
        if (z.gt.zmax) zmax=z+100.
            cki = cos(chi)
        ski = sin(chi)
        clat=cos(latref)
        slat=sin(latref)
        cfi =cos(lonref)
        sfi =sin(lonref)
        sindec=sin(decli)
        cosdec=cos(decli)
        rmax=re+zmax
            rlim=re+z
             h=(sqrt(rmax**2-rlim**2*ski**2)-rlim*abs(cki))/xl(nlaguer)
        hver=(zmax-z)/xl(nlaguer)

        ni=0

        if (chi.gt.pis2) then
               cki=-cki

C  grille de gauss

       rmin=rlim*ski
       zmin=rmin-re
       if (zmin .lt. 90.) then
        ni = 0
            nj = 0
C   ds le prog d'appel faire  chap = 1.e30 p ex.
        return
       endif
       dl=rlim*cki

c calcul des abscisses abs des points de Gauss

       do i = 1,ngauss
          absc(i)=xg(i)*cki
          wi(i)=wg(i)*cki*rlim
       enddo

       ni=ngauss
       do i = 1,ni
              rz= sqrt(ski**2 +absc(i)**2)
              zi(i)=rz*rlim -re
c             lati(i)=asin((slat+absc(i)*sindec)/rz)
c             loni(i) =atan2(clat*sfi,clat*cfi-absc(i)*cosdec)
              lati(i)=asin(min(1.,(slat+(absc(i)+cki)*sindec)/rz))
              loni(i) =atan2(clat*sfi,clat*cfi-(absc(i)+cki)*cosdec)
              enddo
       latref=asin(min(1.,slat+2.*cki*sindec))
       lonref=atan2(clat*sfi,clat*cfi-2.*cosdec*cki)
       clat=cos(latref)
       slat=sin(latref)
       cfi =cos(lonref)
           sfi =sin(lonref)
      endif

c calcul des abscisses abs des points de Laguerre


        do i = 1,nlaguer
       absc(ni+i)= xl(i)*h/rlim
       wi(ni+i) = wl(i)*h
      enddo

C  Calcul des coordonnes geographiques aux points de grille
C         sur la direction solaire

      do i =ni+1,nlaguer+ni
           rz= sqrt(1+absc(i)**2 + 2.*absc(i)*cki)
           zi(i)=rz*rlim -re
           lati(i)=asin(min(1.,(slat+absc(i)*sindec)/rz))
           loni(i) =atan2(clat*sfi,clat*cfi-absc(i)*cosdec)
          enddo

        ni=ni+nlaguer

        do i = 1,ni
           lati(i)=lati(i)*rad2deg
           loni(i)= loni(i)*rad2deg - 15.*tu/3600.
        enddo
C  Calcul de l'altitude aux points de grille
C         sur la direction verticale



        do i = 1,nlaguer
               zj(i)=rlim + xl(i)*hver -re
           wj(i) = wl(i)*hver
              enddo
        nj = nlaguer

        end subroutine grille
