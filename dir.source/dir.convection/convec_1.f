        subroutine convec_1(iyd,tu,kp,dlongeo,dlatgeo,
     &			    dlonmag,dlatmag,dlonref,dt,psi0,flgpot)
        implicit none

        real*8 tu,dt
        real kp
        real*8 dlonmlt,lat,lon,dlat,dlon,dlon1,dlat1
        real*8 EE(2),v0
        real*8 ca,sa,re,ctet,dx,dy,pideg
        real*8 dlatgeo,dlongeo,dlonref,dtmag
        real*8 dlatgeo0,dlongeo0,dlonmlt0
        real*8 err,dcoef
        real*8 dlonmag0,dlatmag0
        real*8 dpsi0,dpsi10,dpsi20,dtheta,dtheta_1
        real*8 psi,psi0,psi1,psi2
        real*8 cor_cnv_1,cor_cnv,err1
        real*8 vh,vp,vpnord,vpest,dlonmag,dlatmag
        real*8 dlatmag1,dlonmlt1,dlonmag1,dlongeo1,dlatgeo1
        real*8 dlonmlt2,dlatmag2,dlonmlt3,dlatmag3
        real*8 dlonmlt4,dlatmag4,dlon0,dlat0,delta
        real*8 dt1,dt2,dt20,distance,transit_1,errdt,tol
        complex*16 cpsi
        real*8 x1,x2,x3,y1,y2,y3,norm,ax,bx,cx,a1,a2,a3
        real zref,year
        integer i,j,ii,nboucle,nbou1,nbou2,iyd

        real latgeo,longeo,pot
        real latgeo0,longeo0,latgeo1,longeo1
        real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0
        real Bmag,dipangle,Enord,Eest,vperpnord,vperpest,vhorizon,vpara
        real orient,chi
        real B,dip,or,ddp,Jtop
                integer ikp
        logical flgpot

        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi,
     &                  Fe0,Ee0,Fi0,Ei0,
     &                  Bmag,dipangle,Enord,Eest,
     &                  vperpnord,vperpest,vhorizon,vpara,ddp,Jtop


        real*8 rad2deg,deg2rad,pi
        data deg2rad/1.745329251994330d-002/
        data rad2deg/57.295779513082320d0/
        data pi/3.141592653589793d0/

        data tol/1.d-5/

        real*8 lat_top
        data lat_top/89.9d0/

        re=6.378d6
        year=1995.
        zref=300.

c   dlatgeo=latgeo
c   dlongeo=longeo

        latgeo=dlatgeo
        longeo=dlongeo

        latgeo0=latgeo
        longeo0=longeo

c	print*,'entree convec_1',tu,dlatgeo,dlongeo,dlatmag,dlonmag,dlonmlt

c        call geo2mag(dlatgeo,dlongeo,dlatmag,dlonmag,dlonref)

c          dtmag=(dlonmag+dlonref)/15.d0+ (tu-dt)/3600.d0
        dtmag=(dlonmag+dlonref)/15.d0+ tu/3600.d0
        dtmag=mod(dtmag+24.d0,24.d0)
        tmag=dtmag
        dlonmlt0=dtmag*15.d0

        dlatmag0=dlatmag
        dlonmag0=dlonmag

c	print*,'calcul',tu,dlatgeo,dlongeo,dlatmag,dlonmag,dlonmlt0

        call potentiel(iyd,tu-dt,kp,
     &                   dlonmlt0,dlatmag0,EE(1),EE(2),psi,ddp)
        Eest=EE(1)
        Enord=EE(2)
c	print*,'YY',tu,dlonmlt0,dlatmag,psi0,psi
        if (flgpot) psi0=psi


        pot=psi0

        call magfild(latgeo,longeo,zref,year,Bmag,dipangle,orient)

        vpest =-Enord/Bmag*10.d0
        vh= Eest/Bmag*10.d0/dcos(dipangle*deg2rad)

        v0   = dt/re
        dlon = vpest*v0/dcos(min(dlatmag,lat_top)*deg2rad)
        dlat = vh*v0

	do nboucle=1,3
c        err=1.
c        nboucle=0
c        do while (err.gt.tol.and.nboucle.lt.4)
c          nboucle=nboucle+1

          dlatmag=dlatmag0
          dlonmlt=dlonmlt0

          dtheta_1=cor_cnv_1(iyd,tu-dt,kp,
     &                       dlonmlt,dlatmag,dlat,dlon,psi0)

c          print*,'convec_1',tu,dtheta,dtheta_1

          dlonmag=dlonmlt-(tu-dt)/240.d0-dlonref
c            dlonmag=dlonmlt-tu/240.d0-dlonref
          dlonmag=mod(dlonmag+360.d0,360.d0)
          call mag2geo(dlatmag,dlonmag,dlatgeo,dlongeo)

          longeo=dlongeo
          latgeo=dlatgeo

          call magfild(latgeo,longeo,zref,year,Bmag,dipangle,orient)

          call potentiel(iyd,tu-dt,kp,
     &                     dlonmlt,dlatmag,EE(1),EE(2),psi1,ddp)
          Eest=EE(1)
          Enord=EE(2)

          vpest =-Enord/Bmag*10.d0
          vh        = Eest/Bmag*10.d0/dcos(dipangle*deg2rad)

          dlon = vpest*v0/dcos(min(dlatmag,lat_top)*deg2rad)
          dlat = vh*v0
c   	write(*,100)tu,dlonmlt,dlatmag,dlon*rad2deg,dlat*rad2deg,
c     &dtheta_1,err
c   	write(*,100)tu,dlonmlt,dlatmag,psi1,psi0,psi1-psi0

100	format(9(1x,g15.8))
c	  dlongeo0=dlongeo
c	  dlatgeo0=dlatgeo
c	  dlonmag0=dlonmag
c	  dlatmag0=dlatmag
c          call convec(iyd,tu,kp,dlongeo0,dlatgeo0,
c     &			    dlonmag0,dlatmag0,dlonref,dt,psi1,flgpot)
c          err=abs(psi1-psi0)

        enddo

        lonmag=dlonmag
        latmag=dlatmag

        dlatgeo1=dlatgeo
        dlongeo1=dlongeo
c	call convec(iyd,tu,dlongeo1,dlatgeo1,dt,psi1,flgpot)

c   print*,tu,flgpot,psi0,psi,pot

c        write(*,100) tu,tu,psi0,psi1,
c     &dlatgeo,dlongeo,dlatgeo1,dlongeo1
c        write(*,100)'YY',tu,tu,psi0,
c     &dlatmag,dlonmlt,dlatgeo,dlongeo
c	print*,'sortie',tu,dlatmag,dlonmag,psi0,EE

c	print*,'sortie convec_1',tu,dlatgeo,dlongeo,dlatmag,dlonmag,dlonmlt
        return
        end


        subroutine integ_1(lat,lon,dlat,dlon)
        implicit none
        real*8 lat,lon,x,x2,lon1,lat1
        real*8 clon,slon,clat,slat,det,y,z,dlat2
        real*8 dlat,dx,dlon,dy,ca,sa,cb,sb
        real*8 rad2deg,deg2rad,pi
        data rad2deg/57.2957795130823229d0/
        data deg2rad/0.174532925199432955d-01/
        data pi/3.14159265358979312d0/



	lon=lon-dlon*rad2deg
	lat=lat-dlat*rad2deg
	if (lat.ge.90.d0) then
	  lat=180.d0-lat
	  lon=lon+180.d0
	endif
	return


c                dx=dlon
c                x=dcos(lat*deg2rad)
c                x2=x*x
c                if (x2.lt.dx*dx) then
c                  dy=0.d0
c                else
c                  dy=dsqrt(x2-dx**2)
c                endif
c   if (dlat.lt.dy) dy=-dy
c                lat1=dacos(dlat-dy)*rad2deg
c                lon=lon*deg2rad
c                ca=dcos(lon)
c                sa=dsin(lon)
c                cb=dx*sa-dy*ca
c                sb=-dx*ca-dy*sa
c                lon1=datan2(sb,cb)*rad2deg
c                lon=mod(lon1+360.d0,360.d0)
c   lat=lat1

        clat=dcos(lat*deg2rad)
        slat=dsin(lat*deg2rad)
        clon=dcos(lon*deg2rad)
        slon=dsin(lon*deg2rad)
        x=clat*clon
        y=clat*slon
        z=slat

        dlat2=dlat*dlat
        det=dsqrt(1.d0+dlat2-z**2)

        lat=dacos((dlat*z+det)/(1+dlat2))*rad2deg

        lon=datan2(y*det-x*dlon,x*det+y*dlon)*rad2deg
                lon=mod(lon+360.d0,360.d0)

                return
                end

        real*8 function cor_cnv_1(iyd,tu,kp,lonmlt,latmag,
     &                                dlat,dlon,psi0)

        implicit none

        real*8 tu
        real kp
        real*8 dlon,dlat,latmag,lonmlt,psi0
        real*8 ca,sa,dlon1,dlat1,xa,lat,lon
        real*8 rad2deg,deg2rad,pi
        INTEGER ITMA        ,iyd
        REAL*8 tol,x1,x2,func,EPS,potar_1,deupi
        PARAMETER (ITMA        =100,EPS=3.d-8)
        INTEGER iter
        REAL*8 a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
        real*8 b0,fb0
        real*8 angle_ref,angle_max,ds,a1,a2,delta
        real*8 dlo,dla,coef_ds
        logical flgini
        data tol/1.d-8/
        data deg2rad/1.745329251994330d-002/
        data rad2deg/57.295779513082320d0/
        data pi/3.141592653589793d0/
        data angle_ref,angle_max,coef_ds/0.d0,0.1d0,1.d-2/


        deupi=2.d0*pi
        a=angle_ref
        lat=latmag
        lon=lonmlt
        dla=dlat*coef_ds
        dlo=dlon*coef_ds
        ds=dsqrt(dlo**2+dla**2)
        fa=potar_1(iyd,tu,kp,lonmlt,latmag,dlat,dlon,a,psi0)
        b=pi/2.d0
        fb=potar_1(iyd,tu,kp,lonmlt,latmag,dla,dlo,b,psi0)
        c=-b
        fc=potar_1(iyd,tu,kp,lonmlt,latmag,dla,dlo,c,psi0)
        a2=(fb+fc)/2.d0/ds**2
        a1=(fb-fc)/2.d0/ds
        if (a2.ne.0.d0) then
          delta=a1*a1+4.d0*fa*a2
          if (delta.gt.0.d0) then
            if (fb.ge.fc) then
              b=(a1-dsqrt(delta))/2.d0/a2
            else
              b=(a1+dsqrt(delta))/2.d0/a2
            endif
            b=2.d0*b/ds*coef_ds
          else
            b=angle_max
          endif
        else
          b=2.d0*fa/a1/ds*coef_ds
        endif
        fb=potar_1(iyd,tu,kp,lonmlt,latmag,dlat,dlon,b,psi0)
        flgini=.true.
        do while(fa*fb.gt.0.d0)
          if (flgini) then
            b=1.1*b
            flgini=.false.
          else
            b=-b
            flgini=.true.
          endif
          fb=potar_1(iyd,tu,kp,lonmlt,latmag,dlat,dlon,b,psi0)
        enddo
        c=b
        fc=fb
        do 12 iter=1,ITMA
          if((fb.gt.0.d0.and.fc.gt.0.d0)
     &           .or.(fb.lt.0.d0.and.fc.lt.0.d0))then
                c=a
                fc=fa
                d=b-a
                e=d
          endif
          if(abs(fc).lt.abs(fb)) then
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
          endif
          tol1=2.d0*EPS*abs(b)+0.5d0*tol
          xm=.5d0*(c-b)
          if(abs(xm).le.tol1 .or. fb.eq.0.) goto 999
          if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
                s=fb/fa
                if(a.eq.c) then
                  p=2.d0*xm*s
                  q=1.d0-s
                else
                  q=fa/fc
                  r=fb/fc
                  p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
                  q=(q-1.d0)*(r-1.d0)*(s-1.d0)
                endif
                if(p.gt.0.d0) q=-q
                p=abs(p)
                if(2.d0*p .lt. min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
                  e=d
                  d=p/q
                else
                  d=xm
                  e=d
                endif
          else
                d=xm
                e=d
          endif
          a=b
          fa=fb
          if(abs(d) .gt. tol1) then
                b=b+d
          else
                b=b+sign(tol1,xm)
          endif
          b=mod(b,deupi)
          fb=potar_1(iyd,tu,kp,lonmlt,latmag,dlat,dlon,b,psi0)
12	continue
999	continue
        cor_cnv_1=b
        ca=dcos(cor_cnv_1)
        sa=dsin(cor_cnv_1)
        dlon1=dlon*ca-dlat*sa
        dlat=dlon*sa+dlat*ca
        dlon=dlon1
        call integ_1(latmag,lonmlt,dlat,dlon)

        return
        end

        real*8 function potar_1(iyd,tu,kp,lonmlt,latmag,
     &                          dlat,dlon,dtheta,psi0)

        implicit none

        real*8 tu
        real kp,ddp
        integer iyd
        real*8 ca,sa,dx,dy,dlon,dlat,dtheta,x,y
        real*8 lat,lon,latmag,lonmlt,psi,psi0

        ca=dcos(dtheta)
        sa=dsin(dtheta)
        dx=dlon*ca-dlat*sa
        dy=dlon*sa+dlat*ca
        lat=latmag
        lon=lonmlt
        call integ_1(lat,lon,dy,dx)
        call potentiel(iyd,tu,kp,lon,lat,x,y,psi,ddp)
        potar_1=psi-psi0

c   print*,dtheta,psi0,psi,dx,dy,lonmlt,latmag,lon,lat,dy/dx,x/y
        return
        end
