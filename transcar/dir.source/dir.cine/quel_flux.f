c
c------------------------- fluxkappa ----------------------------
c
        subroutine fluxkappa (nango2,nen,centE,isotro,
     .                  gmu,fluxdown,fluxup,E0,FE)
c
c
c
        implicit none
        include 'TRANSPORT.INC'
c
        real centE(nbren),gmu(nbrang),flux_val
        real fluxup(nbren,nbrango2),fluxdown(nbren,nbrango2)
       integer nen,isotro,nango2,ind
       real x(4),kappa,Ecutoff,xE,E0,FE
c
       real phi0,fludesc
       integer ien,iang
c
       call zeroit(fluxup,nbren*nbrango2)
       call zeroit(fluxdown,nbren*nbrango2)
c
!        open(44,file='fluxtop.dat',status='replace')
        do ien=1,nen
         fludesc=flux_val(centE(ien),E0,FE)
          if(isotro.eq.2)then
c           On est totalement forward :
            fluxdown(ien,1)=fludesc
            fluxup(ien,1)=1.0e-33
            do iang=2,nango2
              fluxdown(ien,iang)=1.0e-33
              fluxup(ien,iang)=1.0e-33
            enddo
          else
            do iang=1,nango2
              if(isotro.eq.1)then
c               The integral is over mu.dmu =1/2
                fluxdown(ien,iang)=fludesc*2.
                fluxup(ien,iang)=fludesc*2.
              else
c               The integral is over mu**2.dmu =1/3.
                fluxdown(ien,iang)=fludesc*gmu(iang)*3.
                fluxup(ien,iang)=abs(fludesc*gmu(iang+nango2))*3.
              endif
            enddo
          endif
!          write(44,*) centE(ien),fludesc
        enddo
!        close(44)
c
        return
        end

c
c------------------------- inmaxwl ------------------------------
c
        subroutine inmaxwl (Fe,Eave,nango2,nen,centE,isotro,
     .                  gmu,fluxdown,fluxup)
c
c       Calcul du flux maxwellien
c       isotro = parametre d'isotropie du flux precipite :
c       1 = flux isotrope
c       0 = distribution gaussienne
c       2 = dirac
c
        implicit none
        include 'TRANSPORT.INC'
c
        real centE(nbren),gmu(nbrang)
        real fluxup(nbren,nbrango2),fluxdown(nbren,nbrango2)
       integer nen,isotro,nango2
       real Eave,Fe
c
       real phi0,fludesc
       integer ien,iang
c
        phi0 = Fe*1.e+03/(2.*Eave**3)
       call zeroit(fluxup,nbren*nbrango2)
       call zeroit(fluxdown,nbren*nbrango2)
c
        do ien=1,nen
          fludesc=phi0*centE(ien)*exp(-centE(ien)/Eave)
          fludesc=max(fludesc/3.,1.e-33)
          if(isotro.eq.2)then
c           On est totalement forward :
            fluxdown(ien,1)=fludesc
            fluxup(ien,1)=1.0e-33
            do iang=2,nango2
              fluxdown(ien,iang)=1.0e-33
              fluxup(ien,iang)=1.0e-33
            enddo
          else
            do iang=1,nango2
              if(isotro.eq.1)then
c               The integral is over mu.dmu =1/2
                fluxdown(ien,iang)=fludesc*2.
                fluxup(ien,iang)=fludesc*2.
              else
c               The integral is over mu**2.dmu =1/3.
                fluxdown(ien,iang)=fludesc*gmu(iang)*3.
                fluxup(ien,iang)=abs(fludesc*gmu(iang+nango2))*3.
              endif
            enddo
          endif
        enddo
c
        return
        end
c
c------------------------- indirac ------------------------------
c
       subroutine indirac (Fe,Eave,nango2,nen,centE,isotro,
     .                  gmu,fluxdown,fluxup)
c
c       Calcul du flux monoenergetique
c       isotro = parametre d'isotropie du flux precipite :
c       1 = flux isotrope
c       0 = distribution gaussienne
c       2 = dirac
c
        implicit none
       include 'TRANSPORT.INC'
c
        real centE(nbren),gmu(nbrang)
        real fluxup(nbren,nbrango2),fluxdown(nbren,nbrango2)
       integer nen,isotro,nango2
       real Eave,Fe
c
        real phi0,fludesc
        integer ien,iang,ndirac
c
       phi0 = Fe*1.00e+03
       call zeroit(fluxup,nbren*nbrango2)
       call zeroit(fluxdown,nbren*nbrango2)
c
c       recherche de l'energie la plus proche de Eave (inferieure).
       if(centE(nen).le.Eave)then
         ndirac = nen
         go to 20
        endif
       do ien=2,nen
         if(centE(ien).gt.Eave)then
           ndirac=ien-1
           go to 20
         endif
        enddo
20       continue
c
       fludesc = phi0
       if(isotro.eq.2)then
c         On est totalement forward :
          fluxdown(ndirac,1)=phi0
         fluxup(ndirac,1)=0.
       else
        do iang=1,nango2
          if(isotro.eq.1)then
c             On est isotrope
              fluxdown(ndirac,iang)=phi0
             fluxup(ndirac,iang)=phi0
          elseif (isotro.eq.0)then
c            On se sert des poids gaussiens
              fluxdown(ndirac,iang)=phi0*gmu(iang)
             fluxup(ndirac,iang)=abs(phi0*gmu(iang+nango2))
          endif
          enddo
       endif
c
       return
       end
c
c------------------------- ingauss ------------------------------
c
       subroutine ingauss (Fe,Eave,nango2,nen,centE,isotro,
     .                  gmu,fluxdown,fluxup)
c
c       Calcul du flux gaussien
c
        implicit none
        include 'TRANSPORT.INC'
c
c       Calcul du flux gaussien
c       isotro = parametre d'isotropie du flux precipite :
c       1 = flux isotrope
c       0 = distribution gaussienne
c       2 = dirac
c
        real centE(nbren),gmu(nbrang)
        real fluxup(nbren,nbrango2),fluxdown(nbren,nbrango2)
       integer nen,isotro,nango2
       real Eave,Fe
c
        real phi0,fludesc,ehalf
        integer ien,iang
c
       ehalf = 0.1*Eave
       phi0 = Fe*1.00e+03/(ehalf*Eave)
       call zeroit(fluxup,nbren*nbrango2)
       call zeroit(fluxdown,nbren*nbrango2)
c
       do ien=1,nen
        fludesc=phi0*exp(-((centE(ien)-Eave)/ehalf)**2)
         if(isotro.eq.2)then
c           On est totalement forward :
            fluxdown(ien,1)=fludesc
           fluxup(ien,1)=0.
          do iang=2,nango2
              fluxdown(ien,iang)=0.
             fluxup(ien,iang)=0.
           enddo
         else
          do iang=1,nango2
            if(isotro.eq.1)then
c             Distribution isotrope
                fluxdown(ien,iang)=fludesc
               fluxup(ien,iang)=fludesc
            elseif (isotro.eq.0)then
c               On se sert des poids gaussiens.
                fluxdown(ien,iang)=fludesc*gmu(iang)
               fluxup(ien,iang)=abs(fludesc*gmu(iang+nango2))
            endif
            enddo
         endif
       enddo
c
       return
       end
c
c ------------------------ normflux ---------------------------------
c
       subroutine normflux(Fe,nango2,nen,centE,ddeng,
     .                  gmu,gwt,fluxdown,fluxup)
c
c       Normalise le flux d'entree a une valeur en energie Fe donnee.
c
        implicit none
        include 'TRANSPORT.INC'
c
        real centE(nbren),ddeng(nbren),gmu(nbrang),gwt(nbrang)
        real fluxup(nbren,nbrango2),fluxdown(nbren,nbrango2)
        integer nen,isotro,nango2,iang,ien
        real Fe,qsum,qtot,xnorm
c
c       Compute input energy in eV/cm2/sec/sr
       qsum = Fe*1.00e+03      ! total energy in eV/cm2/sec/sr
      qtot = 0.            ! total flux energy in eV/cm2/sec/sr
       do iang=1,nango2
         do ien=1,nen
           qtot=qtot+fluxdown(ien,iang)*gwt(iang)*
     .            gmu(iang)*centE(ien)*ddeng(ien)
         enddo
        enddo

       xnorm = qsum/qtot
c
      do iang=1,nango2
        do ien=1,nen
          fluxdown(ien,iang) = fluxdown(ien,iang)*xnorm
          fluxup(ien,iang)   = fluxup(ien,iang)*xnorm
         enddo
       enddo
c       Initialisation des flux trop petits a 1.e-05
      do iang=1,nango2
        do ien=1,nen
          fluxdown(ien,iang) = max(fluxdown(ien,iang),1.e-05)
          fluxup(ien,iang)   = max(fluxup(ien,iang),1.e-05)
         enddo
       enddo
c
c       Verification ...
      qtot = 0.
       do iang=1,nango2
         do ien=1,nen
           qtot=qtot+fluxdown(ien,iang)*gwt(iang)*
     .            gmu(iang)*centE(ien)*ddeng(ien)
         enddo
        enddo
       qtot=qtot            ! total energy input in eV/cm2/sec/sr
      write(6,1000)Fe,qtot/1.00e+03,xnorm
1000  format ('Energie integree :',/,9x,1pe10.2,'keV/cm2/s/sr desire',
     &    /,9x,1pe10.2,' keV/cm2/s/sr calcule  ',
     &    '(facteur de normalisation :',0p1f10.6,')')
c
c
       return
       end
