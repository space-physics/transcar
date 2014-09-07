	subroutine coef_pot(iyd,tu,kp,phipot,Lmin,Lmax,latequi,ddp,ierr)
	
	implicit none

        include 'TRANSPORT.INC'

        integer npt,ndeg0,mdeg0
        parameter(ndeg0=5,mdeg0=5,npt=(ndeg0+1)*(mdeg0+1))
        integer ndeg,mdeg,ierr,i
	integer iyddeb,iydfin,iyd
	real*8 tu
	real kp,ddp
	real*8 latequi,Lmin,Lmax
	complex*16 phipot(1),phi(npt,3)
	
	real*8 tempsdeb,tempsfin,xt,xtd,xtf
	real latlim,Linf,Lsup
	logical flgini

	data iyddeb    ,tempsdeb  ,iydfin    ,tempsfin  ,flgini
     &	    /3000365.d0,24.d0     ,1980001.d0,0.d0      ,.true./


	data phi /
c       	Kp = 0 a 2-
     &	( 0.293,-0.073),(-0.810,1.432),(-0.144,1.318),( 0.051, 0.513),
     & 	( 0.029, 0.100),(-0.035,-0.014),
     &  (-0.226, 0.098),(1.103,-1.740),( 0.620,-0.723),(0.241,-0.579),
     &	(0.290,-0.046),(-0.025,0.017),
     &	(-0.097,-0.030),(-0.462,0.807),(-0.422,0.399),(-0.102,0.042),
     &	(-0.164,-0.006),(-0.211,-0.063),
     & 	(0.026,0.004),(-0.147,-0.068),(-0.109,-0.089),(-0.215,0.078),
     &	(-0.071,0.025),(0.034,0.118),
     &  (-0.014,0.005),(0.009,-0.075),(0.127,-0.095),(0.086,-0.020),
     &	(0.035,-0.068),(0.052,-0.006),
     &	(0.017,-0.005),(0.048,0.042),(0.003,0.050),(0.042,0.012),
     &	(0.014,0.026),(0.003,-0.036),
c       	Kp = 2 a 4-
     &	(0.366,-0.304),(-0.475,7.321),(-1.470,2.696),(-0.170,0.358),
     &	(0.072,0.239),(0.070,0.034),
     &	(-0.145,0.319),(4.626,-8.764),(2.234,-1.445),(0.010,0.262),
     & 	(0.387,0.330),(0.117,0.175),
     &	(-0.316,0.022),(-1.718,1.601),(0.248,-0.838),(-0.038,-0.467),
     &	(-0.430,-0.427),(-0.364,0.097),
     &	(-0.005,-0.048),(-0.009,0.531),(-0.618,0.809),(0.055,0.118),
     &	(0.210,-0.163),(0.279,-0.266),
     &	(0.086,0.008),(0.175,-0.272),(0.078,-0.220),(-0.122,0.069),
     & 	(0.063,0.108),(0.171,-0.010),
     &	(0.014,0.003),(-0.025,0.014),(0.021,0.022),(0.014,-0.083),
     &	(-0.047,-0.010),(-0.145,0.114),
c       	Kp = 4 a 6-
     &	(1.242,-0.308),(-1.137,17.243),(-3.430,2.003),(0.453,-0.254),
     &	(0.259,-0.101),(-0.081,-0.492),
     &	(-0.382,0.311),(10.731,-11.999),(4.326,0.486),(-1.069,1.401),
     &	(0.003,0.336),(0.296,-0.287),
     &	(-1.096,0.028),(-1.825,-0.792),(2.459,-1.157),(1.278,-0.034),
     &	(-0.286,0.315),(-0.426,0.611),
     &	(0.199,-0.062),(-0.769,1.431),(-1.488,0.317),(0.103,0.003),
     &	(0.346,0.240),(-0.011,-0.135),
     &	(0.159,0.016),(0.092,-0.130),(-0.100,-0.009),(-0.188,0.085),
     &	(0.004,-0.055),(-0.004,-0.265),
     &	(-0.122,0.015),(0.076,-0.281),(0.155,-0.096),(0.057,-0.134),
     &	(-0.096,0.020),(-0.059,0.152)/

       data Linf,Lsup,latlim /61.5d0,72.5d0,60.d0/
	
	if (flgini) then
	  ierr=1	
	  flgini=.false.
	  open(87,file=chemin(1:lpath)
     &                       //'dir.source/dir.convection/varpot.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)
	endif

999	continue
	close(87)

	if (ierr.gt.0) then
	  ndeg=ndeg0
	  mdeg=mdeg0
	  Lmin=Linf
	  Lmax=Lsup
	  latequi=latlim
	  if(kp.le.1.) then
	    do i=1,(ndeg+1)*(mdeg+1)
	      phipot(i)=phi(i,1)
	    enddo
	  else if(kp.le.3.) then
	    do i=1,(ndeg+1)*(mdeg+1)
	      phipot(i)=phi(i,2)
	    enddo
	  else
	    do i=1,(ndeg+1)*(mdeg+1)
	      phipot(i)=phi(i,3)
	    enddo
	  endif
	else
	  open(87,file=chemin(1:lpath)
     &                       //'dir.source/dir.convection/varpot.dat',
     &		form='formatted',status='old',iostat=ierr)
	
	  xt=iyd+tu/360000.d0
          if (iyd.lt.1900000) xt=xt+1900000.

	  xtd=iyddeb+tempsdeb/100.d0
	  xtf=iydfin+tempsfin/100.d0
	  dowhile (xtd.gt.xt.or.xtf.le.xt)

	    read(67,*) iyddeb,tempsdeb,iydfin,tempsfin,
     &                 ndeg,mdeg,Lmin,Lmax,latequi,ddp
	
	    xtd=iyddeb+tempsdeb/100.d0
	    xtf=iydfin+tempsfin/100.d0

	    read(87,*) (phipot(i),i=1,(mdeg+1)*(ndeg+1))
	  enddo
	  close(87)
	endif
	
	return
	end
