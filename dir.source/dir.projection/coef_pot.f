      subroutine coef_pot(iyd,tu,kp,
     &              ndeg,mdeg,phipot,Lmin,Lmax,latequi,ddp)

      include 'comm.f'

      implicit none
      
      include 'comm_sp.f'

      include 'TRANSPORT.INC'
      

      real(dp),intent(in) :: tu
      real,intent(in) :: kp
      integer,intent(in) :: iyd
      integer,intent(out) :: ndeg,mdeg
      real(dp),intent(out) :: phipot(*),Lmin,Lmax,latequi
      real,intent(out) :: ddp

      integer,save :: ierr  ! file read error
      
      real(dp),parameter :: Linf=61.5_dp,Lsup=72.5_dp,latlim=55._dp
      integer, parameter :: ndeg0=5, mdeg0=5, npt=(ndeg0+1)*(2*mdeg0+1)

      integer len_coef,len_rec,len_buf,irec
      
      integer iyddeb,iydfin,i
      real buffer(1000)
      real(dp) phi(3*npt)

      real(dp) tempsdeb,tempsfin,xt,xtd,xtf
      logical flgini

      data iyddeb    ,tempsdeb  ,iydfin    ,tempsfin  ,irec,flgini
     &        /3000365.d0,86400.d0  ,1980001.d0,0.d0      ,0   ,.true./


      data phi/
!        0<=Kp<=2-
     &   -2.930d-01,  2.260d-01,  9.700d-02, -2.600d-02,
     &    1.400d-02, -1.700d-02,  8.100d-01,  1.432d+00,
     &   -1.103d+00, -1.740d+00,  4.620d-01,  8.070d-01,
     &    1.470d-01, -6.800d-02, -9.000d-03, -7.500d-02,
     &   -4.800d-02,  4.200d-02,  1.440d-01,  1.318d+00,
     &   -6.200d-01, -7.230d-01,  4.220d-01,  3.990d-01,
     &    1.090d-01, -8.900d-02, -1.270d-01, -9.500d-02,
     &   -3.000d-03,  5.000d-02, -5.100d-02,  5.130d-01,
     &   -2.410d-01, -5.790d-01,  1.020d-01,  4.200d-02,
     &    2.150d-01,  7.800d-02, -8.600d-02, -2.000d-02,
     &   -4.200d-02,  1.200d-02, -2.900d-02,  1.000d-01,
     &   -2.900d-01, -4.600d-02,  1.640d-01, -6.000d-03,
     &    7.100d-02,  2.500d-02, -3.500d-02, -6.800d-02,
     &   -1.400d-02,  2.600d-02,  3.500d-02, -1.400d-02,
     &    2.500d-02,  1.700d-02,  2.110d-01, -6.300d-02,
     &   -3.400d-02,  1.180d-01, -5.200d-02, -6.000d-03,
     &   -3.000d-03, -3.600d-02,
!        2<=Kp<=4-
     &   -3.660d-01,  1.450d-01,  3.160d-01,  5.000d-03,
     &   -8.600d-02, -1.400d-02,  4.750d-01,  7.321d+00,
     &   -4.626d+00, -8.764d+00,  1.718d+00,  1.601d+00,
     &    9.000d-03,  5.310d-01, -1.750d-01, -2.720d-01,
     &    2.500d-02,  1.400d-02,  1.470d+00,  2.696d+00,
     &   -2.234d+00, -1.445d+00, -2.480d-01, -8.380d-01,
     &    6.180d-01,  8.090d-01, -7.800d-02, -2.200d-01,
     &   -2.100d-02,  2.200d-02,  1.700d-01,  3.580d-01,
     &   -1.000d-02,  2.620d-01,  3.800d-02, -4.670d-01,
     &   -5.500d-02,  1.180d-01,  1.220d-01,  6.900d-02,
     &   -1.400d-02, -8.300d-02, -7.200d-02,  2.390d-01,
     &   -3.870d-01,  3.300d-01,  4.300d-01, -4.270d-01,
     &   -2.100d-01, -1.630d-01, -6.300d-02,  1.080d-01,
     &    4.700d-02, -1.000d-02, -7.000d-02,  3.400d-02,
     &   -1.170d-01,  1.750d-01,  3.640d-01,  9.700d-02,
     &   -2.790d-01, -2.660d-01, -1.710d-01, -1.000d-02,
     &    1.450d-01,  1.140d-01,
!        4<=Kp<=6-
     &   -1.242d+00,  3.820d-01,  1.096d+00, -1.990d-01,
     &   -1.590d-01,  1.220d-01,  1.137d+00,  1.724d+01,
     &   -1.073d+01, -1.200d+01,  1.825d+00, -7.920d-01,
     &    7.690d-01,  1.431d+00, -9.200d-02, -1.300d-01,
     &   -7.600d-02, -2.810d-01,  3.430d+00,  2.003d+00,
     &   -4.326d+00,  4.860d-01, -2.459d+00, -1.157d+00,
     &    1.488d+00,  3.170d-01,  1.000d-01, -9.000d-03,
     &   -1.550d-01, -9.600d-02, -4.530d-01, -2.540d-01,
     &    1.069d+00,  1.401d+00, -1.278d+00, -3.400d-02,
     &   -1.030d-01,  3.000d-03,  1.880d-01,  8.500d-02,
     &   -5.700d-02, -1.340d-01, -2.590d-01, -1.010d-01,
     &   -3.000d-03,  3.360d-01,  2.860d-01,  3.150d-01,
     &   -3.460d-01,  2.400d-01, -4.000d-03, -5.500d-02,
     &    9.600d-02,  2.000d-02,  8.100d-02, -4.920d-01,
     &   -2.960d-01, -2.870d-01,  4.260d-01,  6.110d-01,
     &    1.100d-02, -1.350d-01,  4.000d-03, -2.650d-01,
     &    5.900d-02,  1.520d-01/


      if (flgini) then
          ierr=1
          flgini=.false.
          
          open(87,file='dir.data/dir.linux/dir.projection/varpot.dat',
     &       form='unformatted',access='direct',status='old',recl=40,
     &         iostat=ierr,err=999)

          read(87,rec=1)(buffer(i),i=1,10)
          ndeg=buffer(5)
          mdeg=buffer(6)
          len_coef=(2*mdeg+1)*(ndeg+1)
          write(stdout,*),'len_coef=',len_coef
          len_buf=len_coef+10
          len_rec=4*len_buf
      endif

999    continue
      close(87)
      if (ierr.gt.0) then 
      
        call cpu_time(tic)
        print *,tic,
     &  ' Error reading varpot.dat, fallback to defaults  kp=', kp
    !if there was an error reading varpot.dat, fallback to default param
          ndeg=ndeg0
          mdeg=mdeg0
          Lmin=Linf
          Lmax=Lsup
          latequi=latlim

          if(kp.le.1.) then
            do i=1,(ndeg+1)*(mdeg+1)
              phipot(i)=phi(i)
            enddo
          elseif(kp.le.3.) then
            do i=1,(ndeg+1)*(mdeg+1)
              phipot(i)=phi(i+npt)
            enddo
          else
            do i=1,(ndeg+1)*(mdeg+1)
              phipot(i)=phi(i+2*npt)
            enddo
          endif
          
      else
! continue to read varpot.dat
           open(87,file='dir.data/dir.linux/dir.projection/varpot.dat',
     &        form='unformatted',access='direct',status='old',
     &         recl=len_rec,iostat=ierr)

           xt=iyd+tu*1.d-6
              if (iyd.lt.1900000) xt=xt+1900000.d0

           xtd=iyddeb+tempsdeb*1.d-6
           xtf=iydfin+tempsfin*1.d-6
           do while (xtd.gt.xt.or.xtf.le.xt)
                if (xtf.gt.xt) irec=max(irec-1,1)
                if (xtf.le.xt) irec=irec+1

                read(87,rec=irec) (buffer(i),i=1,len_buf)
                iyddeb=buffer(1)
                tempsdeb=buffer(2)
                iydfin=buffer(3)
                tempsfin=buffer(4)
                Lmin=buffer(7)
                Lmax=buffer(8)
                latequi=buffer(9)
                ddp=buffer(10)
!          read(87,*) iyddeb,tempsdeb,iydfin,tempsfin,
!    &                 ndeg,mdeg,Lmin,Lmax,latequi,ddp
                xtd=iyddeb+tempsdeb*1.d-6
                xtf=iydfin+tempsfin*1.d-6

!        read(87,*) (phipot(i),i=1,(2*mdeg+1)*(ndeg+1))
                do i=1,len_coef
                  phipot(i)=buffer(10+i)
                end do
           end do
           close(87)

      End If

      end subroutine coef_pot
