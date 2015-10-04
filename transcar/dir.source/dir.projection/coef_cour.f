	subroutine coef_cour(iyd,tu,kp,ndeg,mdeg,phicour,
     &			     Lmin,Lmax,latequi,ierr)

	implicit none

        include 'TRANSPORT.INC'

        integer ndeg,mdeg,ierr,i
	integer iyddeb,iydfin,iyd
	real*8 tu
	real kp,ddp
	real*8 latequi,Lmin,Lmax
	complex*16 phicour(1)

	real*8 temps,tempsdeb,tempsfin,xt,xtd,xtf
	logical flgini

	data iyddeb    ,tempsdeb  ,iydfin    ,tempsfin  ,flgini
     &	    /3000365.d0,24.d0     ,1980001.d0,0.d0      ,.true./


	if (flgini) then
	  ierr=1
	  flgini=.false.
	  open(67,file='dir.data/dir.linux/dir.projection/varcourant.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)
	endif

999	continue
	  close(67)

        ndeg=0
        mdeg=0
	if (ierr.eq.0) then

	  open(67,file='dir.data/dir.linux/dir.projection/varcourant.dat',
     &		form='formatted',status='old',iostat=ierr)

	  xt=iyd+tu/360000.d0
	  if (iyd.lt.1900000) xt=xt+1900000.d0
	  xtd=iyddeb+tempsdeb/100.d0
	  xtf=iydfin+tempsfin/100.d0

	  dowhile (xtd.gt.xt.or.xtf.le.xt)

	    read(67,*) iyddeb,tempsdeb,iydfin,tempsfin,
     &                 ndeg,mdeg,Lmin,Lmax,latequi,ddp

	    xtd=iyddeb+tempsdeb/100.d0
	    xtf=iydfin+tempsfin/100.d0

	    read(67,*) (phicour(i),i=1,(2*mdeg+1)*(ndeg+1))
	  enddo
	  close(67)
	endif
	return
	end
