      real function dE_dt(E,Ne,Te)
        implicit none

          real,intent(in) :: E,Ne,Te

          real*8 v,ve,x,Lncoul
          real*8 :: fct_F,fct_G

          real*8,parameter :: Lncoul0=16.47435303047670d0,
     &      const=8.308124389799321d-19, const1=5.185521124430914d-6,
     &      cve=5.505711034864502d3,     cv=5.930968920239551d5

! debye=sqrt(eps0*kb*te/ne/qe**2)
! Lncoul=log(8*pi*eps0*kb*te/qe**2/gamma*debye)=log(8*pi*(eps0*kb)**(3/2)/qe**3/gamma)+3/2*log(Te)-1/2*log(Ne)
      Lncoul=Lncoul0+1.5*log(Te)-.5*Log(Ne)

!      ve=sqrt(2*kb/me*Te)
      ve=cve*sqrt(Te)
!      v=sqrt(2*qe/me*E)
      v=cv*sqrt(E)

      x=v/ve

!   omega2= qe**2/eps0/me*Ne
!  dE_dx=1e6/2/pi**1.5/eps0*omega2*qe**2/ve*(Lncoul*(sqrt(pi)/2*erf_d(x)./x-2*exp(-x.^2))+G)/v
!  dE_dx=1e6/2/pi**1.5/eps0**2*qe**4/me*Ne/ve*(Lncoul*(sqrt(pi)/2*erf_d(x)./x-2*exp(-x.^2))+G)/v
      dE_dt=max(0.d0,const*Ne*(Lncoul*fct_F(x)+fct_G(x))/ve/v)

!      dE_dt=1/2/pi**1.5/eps0**2*qe**3/me*Ne/ve*(Lncoul*(sqrt(pi)/2*erf_d(x)./x-2*exp(-x.^2))+G)
!      dE_dt=max(0.d0,const1*Ne*(Lncoul*fct_F(x)+fct_G(x))/ve)

      end function dE_dt


       real*8 function fct_G(x)
        implicit none
        real*8,intent(in) :: x

        integer i
        integer,parameter :: n1=16,n2=9,n3=10
        real*8,parameter :: x1=4.d0,x2=10.d0
        real*8,parameter:: sqrt_pi3_2=2.65868077635827d0,
     & sqrt2_3=1.12246204830937d0
        real*8 x_1
       real*8,parameter ::
     & p_G1(n1)=(/ -0.00002048624814d0,  0.00069046243775d0,
     &              -0.01043544684716d0,  0.09303613548637d0,
     &              -0.54086367860165d0,  2.13603106149447d0,
     &              -5.77733935730334d0, 10.50893913422101d0,
     &              -12.34413794791237d0,  9.24815805250947d0,
     &              -6.08981347536415d0,  5.07702299840575d0,
     &              -0.54863842550586d0, -3.12958804967662d0,
     &              -0.00398280232498d0,  1.07727401230231d0 /)
      real*8,parameter ::
     & p_G2(n2)=(/-0.00000033809574d0, 0.00002062408644d0,
     &           -0.00054924836885d0, 0.00835350728231d0,
     &           -0.07952212228605d0, 0.48638957620716d0,
     &           -1.86840995049315d0, 4.05569648338900d0,
     &            -2.77845422485329d0 /)
      real*8,parameter ::
     & p_G3(n3)=(/-0.00000056862955d0, 0.00003089756689d0,
     &           -0.00072947821499d0, 0.00982615690014d0,
     &           -0.08342833064046d0, 0.46621855866187d0,
     &           -1.74464065804196d0, 4.41458930978275d0,
     &            -8.05865486553700d0,-0.70055029305823d0 /)

        if (x.le.x1) then
          fct_G=p_G1(1)
          do i=2,n1
            fct_G=fct_G*x+p_G1(i)
          enddo
        elseif (x.le.x2) then
          fct_G=p_G2(1)
          do i=2,n2
            fct_G=fct_G*x+p_G2(i)
          enddo
        else
          x_1=x/10.
          fct_G=p_G3(1)
          do i=2,n3
            fct_G=fct_G*x_1+p_G3(i)
          enddo
          fct_G=log(sqrt2_3*x)/x*sqrt_pi3_2-exp(fct_G)
        endif

      end function fct_G


      real*8 function fct_F(x)

          real*8,intent(in) :: x

          real*8 :: erf_d
          real*8,parameter :: sqrtpi_2=0.88622692545276d0


          if (x.ne.0.d0) then
            fct_F=sqrtpi_2*erf_d(x)/x-2.d0*exp(-x**2)
          else
            fct_F=-1.d0
          endif

      end function fct_F


      real function fzero(Ne,Te)

    !  fonction calculant la racine de dE_dt=0 pour 5<=Lncoul<=25

          real, intent(in) :: Ne,Te

          integer i
          integer,parameter :: n0=16
          real*8 x
          real*8,parameter :: Lncoul0=16.47435303047670d0,
     &    kb_qe=8.617385692256675d-5,
     &    p_zero(n0)=(/-0.00000001747884d0,0.00000070827925d0,
     &                  -0.00001324818763d0,0.00015159023895d0,
     &                  -0.00118445573920d0,0.00667249752758d0,
     &                  -0.02783941459867d0,0.08672516052003d0,
     &                  -0.19951250008517d0,0.32512866229200d0,
     &                  -0.33163827045288d0,0.10773434490681d0,
     &                   0.20475850700617d0,-0.23539720452491d0,
     &                  -0.10085290532879d0,1.27475858388761d0/)


          x=log(Lncoul0+1.5d0*log(Te)-.5d0*log(Ne))

            fzero=p_zero(1)
            do i=2,n0
              fzero=fzero*x+p_zero(i)
            enddo
          fzero=kb_qe*fzero**2*Te

      end function fzero


      pure real*8 function erf_d(x)
          implicit none

          real*8,intent(in) :: x

          erf_d=x
      end function erf_d

