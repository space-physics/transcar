	subroutine prec_time(iyd,tu,Energie,fluxE)
c	implicit none

	integer iyddeb,iydfin,iyd
	real temps,tempsdeb,tempsfin
	logical flgini
	real tu,Energie,fluxE,Flux_elec_int,lonmlt
	real*8 Lmin,Lmax,DL,LM,deg2rad,coeffc,alfa
	parameter (deg2rad=.01745329251994,coeffc=6.)
	real*8 fx,d11,d12,d21,d22,det
	complex*16 phienerg(0:5,0:5),phiflux(0:5,0:5)
	complex*16 expphi(0:5),xcomp,ejphi
	real*8 x,pn,dp,xt,xtd,xtf,x1,x2
	integer ikp,m,n

        real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0
        real Bmag,dipangle,Enord1,Eest1
        real vperpnord,vperpest,vhorizon,vpara
        real orient,chi
        real B,dip,or,ddp,Jtop

        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi,
     &                  Fe0,Ee0,Fi0,Ei0,
     &                  Bmag,dipangle,Enord1,Eest1,
     &                  vperpnord,vperpest,vhorizon,vpara,ddp,Jtop


	data iyddeb  ,tempsdeb,iydfin  ,tempsfin,flgini
     &	    /3000365.,24.     ,1980001.,0.      ,.true./



        include 'TRANSPORT.INC'

	save ierr

	if (flgini) then
	  ierr=1	
	  flgini=.false.
	  open(77,file=data_path(1:lpath_data)
     &                       //'dir.cine/varenerg.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)
	  open(78,file=data_path(1:lpath_data)
     &                       //'dir.cine/varfluxe.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)
	endif

999	continue
	close(77)
	close(78)
	

	if (ierr.gt.0) then

          call hardelec(latmag,tmag,ikp,Flux_elec_int,fluxE)
c
c         Les sorties de hardelec sont en log decimal :
c         Flux de particules integre [cm-2.s-1.sr-1]
          Flux_elec_int=10.**Flux_elec_int
c         Flux d'energie integre [keV.cm-2.s-1.sr-1]
          fluxE=10.**fluxE

c
c         Calcule l'energie moyenne en eV 
c         On assume que la forme du flux est une maxwellienne, bien que
c         hardy et al,(JGR Mai90, p4229-4248, 1985) ont publie des 
c 	  formes de flux dont ils precisent autant que faire se peut 
c 	  qu'elles ne sont justement pas maxwelliennes. Helas, ils ne 
c 	  donnent pas une loi generale de representation de leurs 
c 	  precipitations. Que faire ?
c 	  Pour une maxwellienne F = k E exp(-E/Eo), on a 
c 	  Flux de particules = integrale de  k F(E) dE
c 	  Flux d'energie     = integrale de  k E F(E) dE
c 	  Le rapport des 2 donne 
c 	      Eo = (Flux d'energie/2*Flux de particules)
          Energie=(1000.*fluxE)/(2.*Flux_elec_int)
          fluxE=fluxE/7.
          


        else

	  open(77,file=data_path(1:lpath_data)
     &                       //'dir.cine/varenerg.dat',
     &		form='formatted',status='old',iostat=ierr)
	  open(78,file=data_path(1:lpath_data)
     &                       //'dir.cine/varfluxe.dat',
     &		form='formatted',status='old',iostat=ierr)

	  xt=iyd+tu/360000.d0
	if (iyd.lt.1900000) xt=xt+1900000.
	  xtd=iyddeb+tempsdeb/100.d0
	  xtf=iydfin+tempsfin/100.d0

	  dowhile (xtd.gt.xt.or.xtf.le.xt)

	    read(77,*) iyddeb,tempsdeb,iydfin,tempsfin,Lmin,Lmax
	    read(78,*) iyddeb,tempsdeb,iydfin,tempsfin,Lmin,Lmax

	    xtd=iyddeb+tempsdeb/100.d0
	    xtf=iydfin+tempsfin/100.d0

	    read(77,*) phienerg
	    read(78,*) phiflux
	  enddo
	  close(77)
	  close(78)

	  lonmlt=15.*tmag
	  DL=Lmax-Lmin
	  LM=Lmax+Lmin
	  alfa=-2./DL/deg2rad
	  latequi= 60.

	  ejphi=cmplx(dcos(deg2rad*lonmlt),dsin(deg2rad*lonmlt))
	  expphi(0)=(1.d0,0.d0)
	  do m=1,5
	    expphi(m)=expphi(m-1)*ejphi
	  enddo


	  energie=0.
	  fluxE=0.

          if (latmag.ge.Lmin.and.latmag.le.Lmax) then
	    x=(LM-2.*latmag)/DL
	    do n=0,5
	      call plegendr(n,x,pn,dp)
	      do m=0,5
	        energie=energie-dreal(phienerg(m,n)*pn*expphi(m))
	        fluxE=fluxE-dreal(phiflux(m,n)*pn*expphi(m))
	      enddo
	    enddo

	  else if (latmag.gt.Lmax) then

	    x=dlog(dtan(deg2rad*(45.-latmag/2.)))
	    x0=dlog(dtan(deg2rad*(45.-Lmax/2.)))
	    dxx=x-x0
	    x1=(LM-2.*Lmax)/DL
	    coslsup=dcos(deg2rad*Lmax)
	    do n=0,5
	      call plegendr(n,x1,pn,dp)
	      expcdx=exp(coeffc*dxx)
	      bnm=pn
	      do m=0,5
	        expmdx=exp(m*dxx)
	        anm=-(alfa*dp*coslsup+m*pn)/coeffc
	        Fnm=(anm*(expcdx-1.)+bnm)*expmdx*expphi(m)
                energie=energie-dreal(phienerg(m,n)*Fnm)
                fluxE=fluxE-dreal(phiflux(m,n)*Fnm)
	      enddo
	    enddo

	  else if (latmag.lt.Lmin.and.latmag.ge.latequi) then

	    xequi=log(dtan(deg2rad*(45.-latequi/2.)))
	    x=log(dtan(deg2rad*(45.-latmag/2.)))
	    x0=log(dtan(deg2rad*(45.-Lmin/2.)))
	    dx0=x0-xequi
	    dxx=x-xequi
	    dx2=dxx*dxx
	    x2=(LM-2.*Lmin)/DL
	    coslinf=dcos(deg2rad*Lmin)
	    do n=0,5
	      call plegendr(n,x2,pn,dp)
	      cnm=-alfa*dp*coslinf
	      d12=dx0*dx0
	      d11=x0*d12
	      do m=0,5
	        expmdx=exp(m*dxx)
	        d21=2.*dx0*x0+d12+m*d11
	        d22=2.*dx0+m*d12
	        det=d11*d22-d12*d21
	        anm=(pn*d22-cnm*d12)/det
	        bnm=(-pn*d21+cnm*d11)/det
	        fx=anm*x+bnm
	        Fnm=dx2*fx
	        if (m.eq.0) then
	         Fnm=cnm*(x-x0)+pn
	        endif
	        Fnm=Fnm*expphi(m)
                energie=energie-dreal(phienerg(m,n)*Fnm)
                fluxE=fluxE-dreal(phiflux(m,n)*Fnm)
	      enddo
	    enddo

	  endif

c	dans le fichier initial l'energie est en keV et on la veut en eV
	  energie=exp(energie)*1000./2.
c	dans le fichier initial le flux en energie est en W/m^2 et on le veut en  keV/(cm^2.s.sr)	
	  fluxE=exp(fluxE)*9.9472e10

c	facteur multiplicatif pour le 930216
c	  fluxE=fluxE*6.2832


	endif
	return
	end
