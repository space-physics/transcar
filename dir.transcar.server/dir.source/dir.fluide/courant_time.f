	subroutine courant_time(iyd,tu,Jsup)
c	implicit none

	integer iyddeb,iydfin,iyd
	real temps,tempsdeb,tempsfin
	logical flgini
	real tu,Jsup
	real*8 Lmin,Lmax,DL,LM,deg2rad,coeffc,alfa
	parameter (deg2rad=.01745329251994,coeffc=6.)
	real*8 fx,d11,d12,d21,d22,det
	complex*16 phicourant(0:5,0:5)
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



	real rap
	common /rap_prec/rap

        include 'TRANSPORT.INC'

	save ierr

	if (flgini) then
	  ierr=1	
	  flgini=.false.
	  open(67,file=chemin(1:lpath)
     &                       //'dir.source/dir.fluide/varcourant.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)
	endif

999	continue
	  close(67)
	

	if (ierr.gt.0) then

          Jsup=0.
          
        else

	  open(67,file=chemin(1:lpath)
     &                       //'dir.source/dir.fluide/varcourant.dat',
     &		form='formatted',status='old',iostat=ierr)

	  xt=iyd+tu/360000.d0
	if (iyd.lt.1900000) xt=xt+1900000.
	  xtd=iyddeb+tempsdeb/100.d0
	  xtf=iydfin+tempsfin/100.d0

	  dowhile (xtd.gt.xt.or.xtf.le.xt)

	    read(67,*) iyddeb,tempsdeb,iydfin,tempsfin,Lmin,Lmax

	    xtd=iyddeb+tempsdeb/100.d0
	    xtf=iydfin+tempsfin/100.d0

	    read(67,*) phicourant
	  enddo
	  close(67)

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


	  Jsup=0.

          if (latmag.ge.Lmin.and.latmag.le.Lmax) then
	    x=(LM-2.*latmag)/DL
	    do n=0,5
	      call plegendr(n,x,pn,dp)
	      do m=0,5
	        Jsup=Jsup-dreal(phicourant(m,n)*pn*expphi(m))
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
                Jsup=Jsup-dreal(phicourant(m,n)*Fnm)
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
                Jsup=Jsup-dreal(phicourant(m,n)*Fnm)
	      enddo
	    enddo

	  endif

c	dans le fichier initial le courant est en 10^-7 A/m2
	  Jsup=Jsup/10.
	  Jtop=Jsup

	endif
	return
	end
