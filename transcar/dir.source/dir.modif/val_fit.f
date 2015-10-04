	real function val_fit(lonmlt,latmag,ndeg,mdeg,phi,
     &			      Lmin,Lmax,latequi)

        implicit none

        integer npt
        parameter(npt=1000)

        integer ndeg,mdeg
	real*8 lonmlt,latmag,Lmin,Lmax,latequi
	complex*16 phi(0:mdeg,0:ndeg)
	
	real*8 DL,LM,deg2rad,coeffc,alfa
	parameter (deg2rad=.01745329251994d0,coeffc=6.d0)
	complex*16 expphi(0:npt),ejphi
	real*8 x,pn,dp,val
	real*8 x0,x1,x2,dx0,dxx,dx2,expcdx,expmdx
	real*8 anm,bnm,cnm,Fnm
	real*8 fx,d11,d12,d21,d22,det
	real*8 xequi,coslsup,coslinf
	integer m,n
	
	  	
	DL=Lmax-Lmin
	LM=Lmax+Lmin
	alfa=-2.d0/DL/deg2rad

	ejphi=cmplx(dcos(deg2rad*lonmlt),dsin(deg2rad*lonmlt))
	expphi(0)=(1.d0,0.d0)
	do m=1,mdeg
	  expphi(m)=expphi(m-1)*ejphi
	enddo


	val=0.d0

        if (latmag.ge.Lmin.and.latmag.le.Lmax) then
	  x=(LM-2.d0*latmag)/DL
	  do n=0,ndeg
	    call plegendr(n,x,pn,dp)
	    do m=0,mdeg
	      val=val-dreal(phi(m,n)*pn*expphi(m))
	    enddo
	  enddo
        else if (latmag.gt.Lmax) then
	  x=dlog(dtan(deg2rad*(45.d0-latmag/2.d0)))
	  x0=dlog(dtan(deg2rad*(45.d0-Lmax/2.d0)))
	  dxx=x-x0
	  x1=(LM-2.d0*Lmax)/DL
	  coslsup=dcos(deg2rad*Lmax)
	  do n=0,ndeg
	    call plegendr(n,x1,pn,dp)
	    expcdx=exp(coeffc*dxx)
	    bnm=pn
	    do m=0,mdeg
	      expmdx=exp(m*dxx)
	      anm=-(alfa*dp*coslsup+m*pn)/coeffc
	      Fnm=(anm*(expcdx-1.d0)+bnm)*expmdx*expphi(m)
              val=val-dreal(phi(m,n)*Fnm)
	    enddo
	  enddo
        else if (latmag.lt.Lmin.and.latmag.ge.latequi) then
	  xequi=log(dtan(deg2rad*(45.d0-latequi/2.d0)))
	  x=log(dtan(deg2rad*(45.d0-latmag/2.d0)))
	  x0=log(dtan(deg2rad*(45.d0-Lmin/2.d0)))
	  dx0=x0-xequi
	  dxx=x-xequi
	  dx2=dxx*dxx
	  x2=(LM-2.d0*Lmin)/DL
	  coslinf=dcos(deg2rad*Lmin)
	  do n=0,ndeg
	    call plegendr(n,x2,pn,dp)
	    cnm=-alfa*dp*coslinf
	    d12=dx0*dx0
	    d11=x0*d12
	    do m=0,mdeg
	      expmdx=exp(m*dxx)
	      d21=2.d0*dx0*x0+d12+m*d11
	      d22=2.d0*dx0+m*d12
	      det=d11*d22-d12*d21
	      anm=(pn*d22-cnm*d12)/det
	      bnm=(-pn*d21+cnm*d11)/det
	      fx=anm*x+bnm
	      Fnm=dx2*fx
	      if (m.eq.0) then
	        Fnm=cnm*(x-x0)+pn
	      endif
	      Fnm=Fnm*expphi(m)
              val=val-dreal(phi(m,n)*Fnm)
	    enddo
	  enddo
	endif
	val_fit=val
	end

	subroutine val_pot(lonmlt,latmag,ndeg,mdeg,phi,
     &			   Lmin,Lmax,latequi,
     &			   Enord,Eest,pot)

        implicit none

        integer npt
        parameter(npt=1000)

        integer ndeg,mdeg
	real*8 lonmlt,latmag,Lmin,Lmax,latequi
	real*8 Enord,Eest,pot
	complex*16 phi(0:mdeg,0:ndeg)
	
	real*8 DL,LM,deg2rad,coeffc,alfa,re
	parameter (deg2rad=.01745329251994d0,coeffc=6.d0,re=6.378d0)
	complex*16 expphi(0:npt),ejphi,psi,xcomp
	real*8 x,pn,dp
	real*8 x0,x1,x2,dx0,dxx,dx2,expcdx,expmdx
	real*8 anm,bnm,cnm,Fnm,DFnm
	real*8 fx,d11,d12,d21,d22,det
	real*8 xequi,coslsup,coslinf
	integer m,n

	DL=Lmax-Lmin
	LM=Lmax+Lmin
	alfa=-2.d0/DL/deg2rad


	Eest=0.d0
	Enord=0.d0
	psi=0.d0

	ejphi=cmplx(dcos(deg2rad*lonmlt),dsin(deg2rad*lonmlt))
	expphi(0)=(1.d0,0.d0)
	do m=1,mdeg
	  expphi(m)=expphi(m-1)*ejphi
	enddo

        if (latmag.ge.Lmin.and.latmag.le.Lmax) then
	  x=(LM-2.d0*latmag)/DL
	  do n=0,ndeg
	    call plegendr(n,x,pn,dp)
	    do m=0,mdeg
	      psi=psi-phi(m,n)*pn*expphi(m)
	      Eest =Eest +dimag(phi(m,n)*pn*m*expphi(m))
	      Enord=Enord+dreal(phi(m,n)*dp*expphi(m))
	    enddo
	  enddo
	  Eest =-Eest/dcos(deg2rad*latmag)/re
	  Enord=-2.d0/DL*Enord/deg2rad/re

	else if (latmag.gt.Lmax) then

	  x=dlog(dtan(deg2rad*(45.d0-latmag/2.d0)))
	  x0=dlog(dtan(deg2rad*(45.d0-Lmax/2.d0)))
	  dxx=x-x0
	  x1=(LM-2.d0*Lmax)/DL
	  coslsup=dcos(deg2rad*Lmax)
	  do n=0,ndeg
	    call plegendr(n,x1,pn,dp)
	    expcdx=exp(coeffc*dxx)
	    bnm=pn
	    do m=0,mdeg
	      expmdx=exp(m*dxx)
	      anm=-(alfa*dp*coslsup+m*pn)/coeffc
	      Fnm=(anm*(expcdx-1.)+bnm)*expmdx
	      DFnm=(anm*(coeffc+m)*expcdx+m*(bnm-anm))*expmdx
	      xcomp=phi(m,n)*expphi(m)
              psi=psi-xcomp*Fnm
	      Eest =Eest+dimag(xcomp)*m*Fnm
	      Enord=Enord+dreal(xcomp)*DFnm
	    enddo
	  enddo
	  Eest=-Eest/re/dcos(deg2rad*latmag)
	  Enord=-Enord/re/dcos(deg2rad*latmag)

	else if (latmag.lt.Lmin.and.latmag.ge.latequi) then

	  xequi=log(dtan(deg2rad*(45.d0-latequi/2.d0)))
	  x=log(dtan(deg2rad*(45.d0-latmag/2.d0)))
	  x0=log(dtan(deg2rad*(45.d0-Lmin/2.d0)))
	  dx0=x0-xequi
	  dxx=x-xequi
	  dx2=dxx*dxx
	  x2=(LM-2.d0*Lmin)/DL
	  coslinf=dcos(deg2rad*Lmin)
	  do n=0,ndeg
	    call plegendr(n,x2,pn,dp)
	    cnm=-alfa*dp*coslinf
	    d12=dx0*dx0
	    d11=x0*d12
	    do m=0,mdeg
	      expmdx=exp(m*dxx)
	      d21=2.d0*dx0*x0+d12+m*d11
	      d22=2.d0*dx0+m*d12
	      det=d11*d22-d12*d21
	      anm=(pn*d22-cnm*d12)/det
	      bnm=(-pn*d21+cnm*d11)/det
	      fx=anm*x+bnm
	      Fnm=dx2*fx
	      DFnm=2.d0*dxx*fx+anm*dx2+m*Fnm
	      if (m.eq.0) then
	       Fnm=cnm*(x-x0)+pn
	       DFnm=cnm
	      endif
	      xcomp=phi(m,n)*expphi(m)
              psi=psi-xcomp*Fnm
	      Eest =Eest+dimag(xcomp)*m*Fnm
	      Enord=Enord+dreal(xcomp)*DFnm
	    enddo
	  enddo
	  Eest=-Eest/re/dcos(deg2rad*latmag)
	  Enord=-Enord/re/dcos(deg2rad*latmag)

	endif

	pot=dreal(psi)
	
	return
	end
	
	