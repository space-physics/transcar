      subroutine coef_cour(iyd,tu,kp,ndeg,mdeg,phicour,
     &			     Lmin,Lmax,latequi,ierr)
       use comm, only : dp
        implicit none
 
        integer ndeg,mdeg,ierr,i,u
	    integer iyddeb,iydfin,iyd
	    real(dp) :: tu,latequi,Lmin,Lmax
	    real kp,ddp
	    complex(dp) :: phicour(1)

	    real(dp) :: temps,tempsdeb,tempsfin,xt,xtd,xtf
        logical flgini

      data iyddeb    ,tempsdeb  ,iydfin    ,tempsfin  ,flgini
     &	    /3000365.d0,24.d0     ,1980001.d0,0.d0      ,.true./


      if (flgini) then
          ierr=1
          flgini=.false.
        open(newunit=u,
     &    file='dir.data/dir.linux/dir.projection/varcourant.dat',
     &		form='formatted',status='old',iostat=ierr,err=999)
      endif

999   continue
      close(u)

        ndeg=0
        mdeg=0
      if (ierr==0) then

       open(newunit=u,
     &     file='dir.data/dir.linux/dir.projection/varcourant.dat',
     &		form='formatted',status='old',iostat=ierr)

	  xt=iyd+tu/360000.d0
	  if (iyd < 1900000) xt=xt+1900000.d0
	  xtd=iyddeb+tempsdeb/100.d0
	  xtf=iydfin+tempsfin/100.d0

	  dowhile (xtd.gt.xt.or.xtf.le.xt)

	    read(u,*) iyddeb,tempsdeb,iydfin,tempsfin,
     &                 ndeg,mdeg,Lmin,Lmax,latequi,ddp

	    xtd=iyddeb+tempsdeb/100.d0
	    xtf=iydfin+tempsfin/100.d0

	    read(u,*) (phicour(i),i=1,(2*mdeg+1)*(ndeg+1))
	  enddo
	  close(u)
      endif
      end subroutine coef_cour
