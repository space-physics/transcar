      subroutine coef_prec(iyd,tu,kp,ndeg,mdeg,phienerg,phifluxE,
     &			     Lmin,Lmax,latequi,ierr)

      implicit none

        integer ndeg,mdeg,ierr,i
        integer iyddeb,iydfin,iyd
        real*8 tu
        real kp,ddp
        real*8 latequi,Lmin,Lmax
        complex*16 phifluxE(1),phienerg(1)

	real*8 temps,tempsdeb,tempsfin,xt,xtd,xtf
	logical flgini

	data iyddeb    ,tempsdeb  ,iydfin    ,tempsfin  ,flgini
     &	    /3000365.d0,24.d0     ,1980001.d0,0.d0      ,.true./

	if (flgini) then
	  ierr=1
	  flgini=.false.
	  open(77,file='dir.data/dir.linux/dir.projection/varenerg.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)
	  open(78,file='dir.data/dir.linux/dir.projection/varfluxe.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)
	endif

999	continue
	close(77)
	close(78)


        ndeg=0
        mdeg=0
	if (ierr.eq.0) then

	  open(77,file='dir.data/dir.linux/dir.projection/varenerg.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)
	  open(78,file='dir.data/dir.linux/dir.projection/varfluxe.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)

	  xt=iyd+tu/360000.d0
	  if (iyd.lt.1900000) xt=xt+1900000.d0
	  xtd=iyddeb+tempsdeb/100.d0
	  xtf=iydfin+tempsfin/100.d0

	  dowhile (xtd.gt.xt.or.xtf.le.xt)

	    read(77,*) iyddeb,tempsdeb,iydfin,tempsfin,
     &                 ndeg,mdeg,Lmin,Lmax,latequi,ddp
	    read(78,*) iyddeb,tempsdeb,iydfin,tempsfin,
     &                 ndeg,mdeg,Lmin,Lmax,latequi,ddp


	    xtd=iyddeb+tempsdeb/100.d0
	    xtf=iydfin+tempsfin/100.d0

	    read(77,*) (phienerg(i),i=1,(2*mdeg+1)*(ndeg+1))
	    read(78,*) (phifluxE(i),i=1,(2*mdeg+1)*(ndeg+1))
	  enddo
	  close(77)
	  close(78)
	endif

	end subroutine coef_prec
