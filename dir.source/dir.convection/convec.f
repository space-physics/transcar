       subroutine convec(iyd,tu,kp,dlongeo,dlatgeo,
     &                   dlonmag,dlatmag,dlonref,dt,psi0,flgpot)

       include 'comm.f'

       implicit none
       
       include 'comm_dp.f'

       real(dp),intent(inout) :: dlonmag

       real(dp) dlonmlt,lat,lon,dlat,dlon,dlon1,dlat1
       real(dp) tu,dt
       real kp
       real(dp) EE(2),Ex(2)
       real(dp) ca,sa,re,ctet,dx,dy,clat,transit
       real(dp) dlatgeo,dlongeo,dlonref,dtmag,distance
       real(dp) dlonmlt1,lattransi,lateps,dlonmlt0,dcoef
       real(dp) dpsi0,dpsi10,dpsi20,dtheta,dtheta1,dtheta2
       real(dp) psi,psi0,psi1,psi2,dlon0,dlat0
       real(dp) cor_cnv
       real(dp) dt1,dt2,pideg
       real(dp) vh,vp,vpnord,vpest,dlatmag
       real(dp) loc(2,2),dpot(2)
       complex(zp) cpsi
       real zref,year
       integer i,j,iyd
       logical flgpot,flgini
       data flgini/.true./

       real latgeo,longeo,pot
       real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0
       real Bmag,dipangle,Enord,Eest,vperpnord,vperpest,vhorizon,vpara
       real B,dip,or,ddp,Jtop
       real orient,chi
       integer ikp

        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi,
     &       		Fe0,Ee0,Fi0,Ei0,
     &			Bmag,dipangle,Enord,Eest,
     &			vperpnord,vperpest,vhorizon,vpara,ddp,Jtop

       real(dp) dlonmag0,dlatmag0,dlo0,dla0


        real(dp) lat_top
        data lat_top/89.9_dp/

       save dtmag,dlonmag0,dlatmag0
       
       if (flgini) then
       open(56,file='trace_conv',form='formatted',status='replace')
       flgini=.false.
       endif

       re=6.378d6
       year=1995.
       zref=300.

c       dlatgeo=latgeo
c       dlongeo=longeo
       latgeo=dlatgeo
       longeo=dlongeo

c       print*,'entree convec',longeo,latgeo
c        call geo2mag(dlatgeo,dlongeo,dlatmag,dlonmag,dlonref)
       dlo0=dlonmag-dlonmag0
       dla0=dlatmag-dlatmag0
        dtmag=(dlonmag+dlonref)/15._dp + tu/3600._dp
       dtmag=mod(dtmag+24._dp,24._dp)
       tmag=dtmag
       dlonmlt=dtmag*15._dp

       write(stdout,*),'convec.f: call potentiel, dlonmag,dlonmlt',
     &   dlonmag,dlonmlt
       call potentiel(iyd,tu,kp,dlonmlt,dlatmag ,EE(1),EE(2),psi,ddp)
       
       if (flgpot) psi0=psi
       
       print*,'convec.f: call magfild'  
       call magfild(latgeo,longeo,zref,year,Bmag,dipangle,orient)
       pot=psi0
       Enord=EE(2)
       Eest=EE(1)
       vpest =-EE(2)/Bmag*10._dp
       vperpest=vpest
       vpnord= EE(1)/Bmag*10._dp
       vperpnord=vpnord
       vh =vpnord/cos(dipangle*deg2rad)
       vhorizon=vh
       vp    =vpnord*tan(dipangle*deg2rad)*orient
       vpara=vp
       
        dlon = vpest*dt/re/cos(min(dlatmag,lat_top)*deg2rad)
        dlat =  vh*dt/re
       if (dlon.ne.0.d0 .and. dlat.ne.0.d0) then
         write(stdout,*),'convec.f: call cor_cnv   dlonmlt,dlatmag ',
     &                    dlonmlt,dlatmag
     
         dtheta=cor_cnv(iyd,tu,kp,dlonmlt,dlatmag,
     &			   dlat,dlon,psi0)
         dlonmag=dlonmlt-(tu+dt)/240.d0-dlonref
         dlonmag=mod(dlonmag+360.d0,360.d0)
      write(stdout,*),'convec.f: call mag2geo,',
     &                ' dlatmag,dlonmag,dlatgeo,dlongeo:',
     &                   dlatmag,dlonmag,dlatgeo,dlongeo
         call mag2geo(dlatmag,dlonmag,dlatgeo,dlongeo)
         print*,'mag2geo outputs:',dlatmag,dlonmag,dlatgeo,dlongeo
       endif


       lonmag=dlonmag
       latmag=dlatmag
       longeo=dlongeo
       latgeo=dlatgeo
       dlonmag0=dlonmag
       dlatmag0=dlatmag
c       write(*,100)tu,dlonmlt,dlatmag,Enord,Eest,
c     &dtheta*rad2deg

c        write(*,100)'XX',tu,psi0,
c     &dlatmag,dlonmlt,dlatgeo,dlongeo


c       print*,tu,dlonmlt,dlatmag,dlon*rad2deg,dlat*rad2deg,
c     &dtheta,psi0
!       print*,'attempting trace_conv write'
       write(56,*) tu,dlonmlt,dlatmag,psi0,dlo0,dla0
c100       format(a2,10(1x,g15.8))
100       format(9(1x,g15.8))

       end Subroutine convec


      subroutine integ(lat,lon,dlat,dlon)
        include 'comm.f'        
        implicit none
        include 'comm_dp.f'
        real(dp),intent(inout) :: lat,lon
        real(dp),intent(in) :: dlat,dlon

        real(dp) :: x,y,z,
     & clon,slon,clat,slat,
     & dx,dy,ca,sa,cb,sb

        lon=lon+dlon*rad2deg
        lat=lat+dlat*rad2deg
        if (lat.ge.90.d0) then
            lat=180._dp-lat
            lon=lon+180._dp
        endif
        return


c        dx=dlon
c        dy=(dlat-dcos(lat*deg2rad))
c        lat=dacos(dsqrt(dx**2+dy**2))*rad2deg
c        ca=dcos(lon*deg2rad)
c        sa=dsin(lon*deg2rad)
c        cb=-dx*sa-dy*ca
c        sb=dx*ca-dy*sa
c        lon=datan2(sb,cb)*rad2deg
        clat=cos(lat*deg2rad)
        slat=sin(lat*deg2rad)
        clon=cos(lon*deg2rad)
        slon=sin(lon*deg2rad)
        x=(clat-dlat*slat)*clon-dlon*slon
        y=(clat-dlat*slat)*slon+dlon*clon
        z=slat+dlat*clat

        lat=atan2(z,sqrt(x**2+y**2))*rad2deg
        lon=atan2(y,x)*rad2deg
        lon=mod(lon+360._dp,360._dp)

        end subroutine integ



      double precision function cor_cnv(iyd,tu,kp,lonmlt,latmag,
     &              dlat,dlon,psi0)
       include 'comm.f'
       implicit none
       include 'comm_dp.f'

       integer, intent(in) :: iyd
       real(dp), intent(in) :: tu,psi0
       real, intent(in) :: kp
       
       real(dp), intent(inout) :: dlat,dlon,lonmlt,latmag

       real(dp) ca,sa,dlon1,lat,lon,dlat1,xa
       integer itmax
       real(dp) x1,x2,func,eps,potar,deupi
       parameter (itmax=100,eps=3.d-8)
       integer iter
       real(dp) a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
        real(dp) b0,fb0
        real(dp) ds,a1,a2,delta
        real(dp) dlo,dla
        logical flgini
        real(dp),parameter :: tol=1.d-8,angle_ref=0._dp,
     &    angle_max=0.1_dp,coef_ds=1.d-2



        deupi=2.d0*pi
        a=angle_ref
        lat=latmag
        lon=lonmlt
        dla=dlat*coef_ds
        dlo=dlon*coef_ds
        ds=sqrt(dlo**2+dla**2)
        fa=potar(iyd,tu,kp,lonmlt,latmag,dlat,dlon,a,psi0)
        b=pi/2._dp
        fb=potar(iyd,tu,kp,lonmlt,latmag,dla,dlo,b,psi0)
        c=-b
        fc=potar(iyd,tu,kp,lonmlt,latmag,dla,dlo,c,psi0)
        a2=(fb+fc)/2.d0/ds**2
        a1=(fb-fc)/2.d0/ds
        delta=a1*a1+4.d0*fa*a2
        if (a2.ne.0.d0) then
          delta=a1*a1+4.d0*fa*a2
          if (delta.gt.0.d0) then
            if (fb.ge.fc) then
              b=(a1-sqrt(delta))/2.d0/a2
            else
              b=(a1+sqrt(delta))/2.d0/a2
            endif
            b=2.d0*b/ds*coef_ds
          else
            b=angle_max
          endif
        else
          b=2.d0*fa/a1/ds*coef_ds
        endif
    	fb=potar(iyd,tu,kp,lonmlt,latmag,dlat,dlon,b,psi0)
        flgini=.true.
        do while(fa*fb.gt.0.d0)
          if (flgini) then
            b=1.1*b
            flgini=.false.
          else
            b=-b
            flgini=.true.
          endif
          fb=potar(iyd,tu,kp,lonmlt,latmag,dlat,dlon,b,psi0)
        enddo
       
       ca=cos(b)
       sa=sin(b)
       dlon1=dlon*ca-dlat*sa
       dlat1=dlon*sa+dlat*ca
       lat=latmag
       lon=lonmlt
      
      if (debug) then
       write(stdout,*),'convec.f: cor_cnv 1st call:',
     & ' call integ  lat,lon,dlat,dlon'
     &, lat,lon,dlat,dlon
      endif
       
       call integ(lat,lon,dlat1,dlon1)
       c=b
       fc=fb
       do 12 iter=1,itmax
         if((fb.gt.0.d0.and.fc.gt.0.d0)
     &	     .or.(fb.lt.0.d0.and.fc.lt.0.d0))then
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
         tol1=2.d0*eps*abs(b)+0.5d0*tol
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
         fb=potar(iyd,tu,kp,lonmlt,latmag,dlat,dlon,b,psi0)
12       continue
999       continue
       cor_cnv=mod(b,deupi)
       ca=cos(cor_cnv)
       sa=sin(cor_cnv)
       dlon1=dlon*ca-dlat*sa
       dlat=dlon*sa+dlat*ca
       dlon=dlon1
       lat=latmag
       lon=lonmlt
       
      if (debug) then
      call cpu_time(tic)
      write(stdout,*),tic,' convec.f: cor_cnv 2nd call:',
     & ' call integ  lat,lon,dlat,dlon'
     &, lat,lon,dlat,dlon
      endif
      
       call integ(lat,lon,dlat,dlon)
       latmag=lat
       lonmlt=lon

       end function cor_cnv


       double precision function potar(iyd,tu,kp,lonmlt,latmag,
     &      dlat,dlon,dtheta,psi0)
        include 'comm.f'
       implicit none
       
       integer,intent(in) :: iyd
       real(dp),intent(in) :: tu,lonmlt,latmag,dlat,dlon,dtheta,psi0
       real, intent(in) :: kp
       real ddp
       real(dp) ca,sa,dx,dy,x,y,lat,lon,psi

       ca=cos(dtheta)
       sa=sin(dtheta)
       dx=dlon*ca-dlat*sa
       dy=dlon*sa+dlat*ca
       lat=latmag
       lon=lonmlt
       
       call integ(lat,lon,dy,dx)
       call potentiel(iyd,tu,kp,lon,lat,x,y,psi,ddp)
       potar=psi-psi0

       end function potar
