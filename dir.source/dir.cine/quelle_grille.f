      subroutine quelle_grille(emax,nen,centE,botE,ddeng,
     &                          nang,nango2,gmu,gwt,angzb)

      use, intrinsic :: iso_fortran_env, only : error_unit

      implicit none

      include 'TRANSPORT.INC'
c
c    Entree du programme : La valeur caracteristique de l'energie
c.   des precipitations
c
c    Sorties : nombres d'energie, d'angles
c    grille d'energie (centre de grille, bas, largeur)
c    grille d'angle (angle, cosinus, poids)

      integer :: nen,nang,nango2
      integer :: i1,i2,i3
      real :: centE(nbren),botE(nbren),ddeng(nbren)
      real :: angzb(2*nbrango2),gmu(2*nbrango2),gwt(2*nbrango2)

c     Parametres internes :
      real :: spfac,emax,emin    !calcul eventuel d'une grille propre
      integer :: ien,iang,iost
      character(*), parameter ::
     &        datdegfn = 'dir.data/dir.linux/dir.cine/DATDEG'


c     Pour transferer via DATDEG les caracteristiques de grille a
c     degrad.f
      integer :: type_de_grille,idess(2),iprint1,iprint2,iprint3,
     &  iprint4

      logical :: logint,lt1,lt2,lt3,lt4,lopal

      character(80) :: crsinput,crs,rdt,crsfn
      namelist /DATDEG/ type_de_grille,crsinput,crs,rdt,idess,
     &    iprint1,iprint2,iprint3,iprint4,logint,lt1,lt2,lt3,lt4,
     &    lopal
c
      real,parameter :: pideg=57.29578
c
c     On stoke les parametres de calcul de degrad.
      open(fic_datdeg, file=datdegfn, status='old')

      rewind(fic_datdeg)
      read(fic_datdeg,nml=DATDEG)

      if(type_de_grille == 0) then
        print '(a)', 'DATDEG specifies grid type 0'
        close(fic_datdeg)
        return
      end if

      close(fic_datdeg)


!-MZ
        emax=70000.
!-MZ

      if (type_de_grille == 1) then
c       On utilise les grilles detaillees
!     "use detailed grid"
      if (emax.le.355.8)then
        crs = 'dir.cine/dir.seff/crsa1'
        rdt = 'dir.cine/dir.seff/rdta1'
      else if (emax.gt.355.8 .and. emax.le.720.6)then
        crs = 'dir.cine/dir.seff/crsa2'
        rdt = 'dir.cine/dir.seff/rdta2'
      else if (emax.gt.720.6 .and. emax.le.2193.)then
        crs = 'dir.cine/dir.seff/crsa3'
        rdt = 'dir.cine/dir.seff/rdta3'
      else if (emax.gt.2193. .and. emax.le.4394.)then
        crs = 'dir.cine/dir.seff/crsa4'
        rdt = 'dir.cine/dir.seff/rdta4'
      else if (emax.gt.4394. .and. emax.le.7057.)then
        crs = 'dir.cine/dir.seff/crsa5'
        rdt = 'dir.cine/dir.seff/rdta5'
      else if (emax.gt.7057. .and. emax.le.21700.)then
        crs = 'dir.cine/dir.seff/crsa6'
        rdt = 'dir.cine/dir.seff/rdta6'
      else if (emax.gt.21700. .and. emax.le.43410.)then
        crs = 'dir.cine/dir.seff/crsa7'
        rdt = 'dir.cine/dir.seff/rdta7'
      else if (emax.gt.43410. .and. emax.le.70440.)then
        crs = 'dir.cine/dir.seff/crsa8'
        rdt = 'dir.cine/dir.seff/rdta8'
      else
        nen = 400
        crs = 'dir.cine/dir.seff/crsa'
        rdt = 'dir.cine/dir.seff/rdta'
        emin = 1.e-01
!     On ne va pas au dessus de 140 keV...
!     "do not go above 140keV (?) "
        emax = min(emax,140000.0)
              print*,'call gridpolo'
            call gridpolo (nen,emin,emax,centE,ddeng,spfac)
c           Calcul des energies de bas de grille
            botE(1) = max(centE(1) - ddeng(1)/2.,1.e-03)
            do ien = 1,nen-1
              botE(ien+1) = botE(ien)+ddeng(ien+1)
            enddo
      print*,'using detail grid type 1, crs,rdt: ',crs,rdt
      endif

      else if (type_de_grille.eq.2)then
      if (emax.le.355.8)then
        crs = 'dir.cine/dir.seff/crsb1'
        rdt = 'dir.cine/dir.seff/rdtb1'
      else if (emax.gt.355.8 .and. emax.le.720.6)then
        crs = 'dir.cine/dir.seff/crsb2'
        rdt = 'dir.cine/dir.seff/rdtb2'
      else if (emax.gt.720.6 .and. emax.le.2193.)then
        crs = 'dir.cine/dir.seff/crsb3'
        rdt = 'dir.cine/dir.seff/rdtb3'
      else if (emax.gt.2193. .and. emax.le.4394.)then
        crs = 'dir.cine/dir.seff/crsb4'
        rdt = 'dir.cine/dir.seff/rdtb4'
      else if (emax.gt.4394. .and. emax.le.7057.)then
        crs = 'dir.cine/dir.seff/crsb5'
        rdt = 'dir.cine/dir.seff/rdtb5'
      else if (emax.gt.7057. .and. emax.le.21700.)then
        crs = 'dir.cine/dir.seff/crsb6'
        rdt = 'dir.cine/dir.seff/rdtb6'
      else if (emax.gt.21700. .and. emax.le.43410.)then
        crs = 'dir.cine/dir.seff/crsb7'
        rdt = 'dir.cine/dir.seff/rdtb7'
      else if (emax.gt.43410. .and. emax.le.70440.)then
        crs = 'dir.cine/dir.seff/crsb8'
        rdt = 'dir.cine/dir.seff/rdtb8'
      else
        nen = 150
        crs = 'dir.cine/dir.seff/crsb'
        rdt = 'dir.cine/dir.seff/rdtb'
        emin = 1.e-01
c      On ne va pas au dessus de 140 keV...
        emax = min(emax,140000.0)
              call gridpolo (nen,emin,emax,centE,ddeng,spfac)
c           Calcul des energies de bas de grille         :
              botE(1) = max(centE(1) - ddeng(1)/2.,1.e-03)
              do ien = 1,nen-1
                botE(ien+1) = botE(ien)+ddeng(ien+1)
              enddo
       endif
      print '(a,1x,a,1x,a)','using detail grid type 2, crs,rdt:',crs,rdt
      endif
c
c       Lecture des grilles d'energie.
c goto 314
      crsfn = 'dir.data/dir.linux/'//crs
      print '(a,1x,a)', 'attempting to open',crsfn
      open(icrsin,file=crsfn, status='OLD',form='UNFORMATTED',
     &     iostat=iost, err=992)
      rewind icrsin
      read(icrsin) nen,i1,i2,i3
      read(icrsin) (centE(ien),ien=1,nen)
      read(icrsin) (botE(ien),ien=1,nen)
      read(icrsin) (ddeng(ien),ien=1,nen)
      close(icrsin)

c314      nen = 40
c      crs = 'dir.cine/dir.seff/crsb1'
c      rdt = 'dir.cine/dir.seff/rdtb1'
c      emin = 1.e-01
cc      On ne va pas au dessus de 140 keV...
c      emax = 355.8
c            call gridpolo (nen,emin,emax,centE,ddeng,spfac)
cc           Calcul des energies de bas de grille
c            botE(1) = max(centE(1) - ddeng(1)/2.,1.e-03)
c            do ien = 1,nen-1
c              botE(ien+1) = botE(ien)+ddeng(ien+1)
c            enddo
c
c  On previent maintenant degrad de ou il faut lire les sections
c  efficaces.
      print '(a,1x,a)', 'attempting to open',datdegfn
            open(fic_datdeg,file=datdegfn, status='replace',err=993,
     &           delim='quote')
      print '(a,1x,a)', 'beginning to rewrite',datdegfn
            write(fic_datdeg,nml=DATDEG)

      close(fic_datdeg)
c
c
        nang=8
        nango2=nang/2

        gmu(1)= .9305681586266
        gmu(2)= .6699905395508
        gmu(3)= .3300094604492
        gmu(4)= .0694318413734
        gmu(5)=-.0694318413734
        gmu(6)=-.3300094604492
        gmu(7)=-.6699905395508
        gmu(8)=-.9305681586266

        gwt(1)=.1739274263382
        gwt(2)=.3260725736618
        gwt(3)=.3260725736618
        gwt(4)=.1739274263382
        gwt(5)=.1739274263382
        gwt(6)=.3260725736618
        gwt(7)=.3260725736618
        gwt(8)=.1739274263382

        do iang=1,nang
          angzb(iang)=pideg*acos(gmu(iang))
        enddo
c
      return

992   write(error_unit,'(a,1x,a,1x,i0)') 'Cross-section file',crsfn,
     &      'is in error. Status=',iost
      error stop

993   write(error_unit,'(a,1x,a)') 'trouble writing crs file',datdegfn
      error stop

      end subroutine quelle_grille
