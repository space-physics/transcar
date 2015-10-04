        subroutine quelle_grille(emax,nen,centE,botE,ddeng,
     &                          nang,nango2,gmu,gwt,angzb)

 	    implicit none

        include 'TRANSPORT.INC'
c
c 	Entree du programme : La valeur caracteristique de l'energie
c 	des precipitations
c
c 	Sorties : nombres d'energie, d'angles
c 	grille d'energie (centre de grille, bas, largeur)
c 	grille d'angle (angle, cosinus, poids)

     	integer nen,nang,nango2
     	integer i1,i2,i3
     	real centE(nbren),botE(nbren),ddeng(nbren)
     	real angzb(2*nbrango2),gmu(2*nbrango2),gwt(2*nbrango2)
c
c 	Parametres internes :
     	real spfac,emax,emin	   !calcul eventuel d'une grille propre
     	integer ien,iang,iost
      character(len=80),parameter ::
     &        datdegfn = 'dir.data/dir.linux/dir.cine/DATDEG'



c 	Pour transferer via DATDEG les caracteristiques de grille a
c 	degrad.f
 	    integer type_de_grille,lenc,idess(2),iprint1,iprint2,iprint3,
     &		iprint4
 	    logical logint,lt1,lt2,lt3,lt4,lopal
        character*80 crsin,crs,rdt,crsfn
c
 	    real,parameter :: pideg=57.29578
c
c 	On stoke les parametres de calcul de degrad.
        open(fic_datdeg, file=datdegfn, status='old')

        rewind(fic_datdeg)
     	read(fic_datdeg,*) type_de_grille
     	if(type_de_grille .eq.0) then
        print*,'DATDEG specifies grid type 0'
     	  close(fic_datdeg)
     	  return
     	endif
        read(fic_datdeg,1000) crsin
        crsin = crsin(1:lenc(crsin))
        read(fic_datdeg,1000) crs
        crs = crs(1:lenc(crs))
        read(fic_datdeg,1000) rdt
        rdt = rdt(1:lenc(rdt))
        call xline(4,fic_datdeg) !skip 4 lines
        read (fic_datdeg,*) idess(1),idess(2)
        call xline(2,fic_datdeg) !skip 2 lines
        read (fic_datdeg,*) iprint1,iprint2,iprint3,iprint4
        read(fic_datdeg,1000) logint
        read(fic_datdeg,1000) lt1
        read(fic_datdeg,1000) lt2
        read(fic_datdeg,1000) lt3
        read(fic_datdeg,1000)lt4
        read(fic_datdeg,1000)lopal
 	    close(fic_datdeg)


!-MZ
        emax=70000.
!-MZ

 	   if (type_de_grille.eq.1)then
c 	  On utilise les grilles detaillees
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
!	    On ne va pas au dessus de 140 keV...
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
c 	    On ne va pas au dessus de 140 keV...
   	    emax = min(emax,140000.0)
              call gridpolo (nen,emin,emax,centE,ddeng,spfac)
c           Calcul des energies de bas de grille         :
              botE(1) = max(centE(1) - ddeng(1)/2.,1.e-03)
              do ien = 1,nen-1
                botE(ien+1) = botE(ien)+ddeng(ien+1)
              enddo
   	   endif
      print*,'using detail grid type 2, crs,rdt: ',crs,rdt
 	    endif
c
c       Lecture des grilles d'energie.
c	goto 314
        crsfn = 'dir.data/dir.linux/'//crs
        print*,'attempting to open ',crsfn
        open(icrsin,file=crsfn,
     &       status='OLD',form='UNFORMATTED',
     &          iostat=iost,err=992)
        rewind icrsin
        read(icrsin) nen,i1,i2,i3
        read(icrsin) (centE(ien),ien=1,nen)
        read(icrsin) (botE(ien),ien=1,nen)
        read(icrsin) (ddeng(ien),ien=1,nen)
        close(icrsin)

c314 	    nen = 40
c 	    crs = 'dir.cine/dir.seff/crsb1'
c 	    rdt = 'dir.cine/dir.seff/rdtb1'
c 	    emin = 1.e-01
cc 	    On ne va pas au dessus de 140 keV...
c 	    emax = 355.8
c            call gridpolo (nen,emin,emax,centE,ddeng,spfac)
cc           Calcul des energies de bas de grille
c            botE(1) = max(centE(1) - ddeng(1)/2.,1.e-03)
c            do ien = 1,nen-1
c              botE(ien+1) = botE(ien)+ddeng(ien+1)
c            enddo
c
c 	On previent maintenant degrad de ou il faut lire les sections
c 	efficaces.
        print*,'attempting to open ',datdegfn
        open(fic_datdeg,file=datdegfn,
     &       status='unknown',err=993) !yes unknown in case not all values rewritten
!        rewind(fic_datdeg)
       print*,'beginning to rewrite',datdegfn
     	write(fic_datdeg,1010)type_de_grille
        write(fic_datdeg,1020) crsin
        write(fic_datdeg,1030) crs
        write(fic_datdeg,1040) rdt
        write(fic_datdeg,1050) idess(1),idess(2)
        write(fic_datdeg,1060) iprint1,iprint2,iprint3,iprint4

 	    if(logint)then
          write(fic_datdeg,1070)
 	    else
          write(fic_datdeg,1075)
 	    endif
1070 	format('.true',11x,
     &    'logint: interpolation type for cross sections')
1075 	format('.false.',10x,
     &    'logint: interpolation type for cross sections')
c
     	if(lt1)then
          write(fic_datdeg,1080)
     	else
          write(fic_datdeg,1085)
     	endif
1080 	format('.true.',11x,'lt1   : test without ionization')
1085 	format('.false.',10x,'lt1   : test without ionization')

     	if(lt2)then
              write(fic_datdeg,1090)
     	else
              write(fic_datdeg,1095)
     	endif
1090 	format('.true.',11x,'lt2   : test only with ionization')
1095 	format('.false.',10x,'lt2   : test only with ionization')

     	if(lt3)then
              write(fic_datdeg,1100)
     	else
              write(fic_datdeg,1105)
     	endif
1100 	format('.true.',11x,'lt3   : test without degradation')
1105 	format('.false.',10x,'lt3   : test without degradation')

     	if(lt4)then
              write(fic_datdeg,1110)
     	else
              write(fic_datdeg,1115)
     	endif
1110 	format('.true.',11x,'lt4   : test without excitation')
1115 	format('.false.',10x,'lt4   : test without excitation')

     	if(lopal)then
              write(fic_datdeg,1120)
     	else
              write(fic_datdeg,1125)
     	endif
1120 	format('.true.',11x,
     &    'lopal : sec. distrib. de Opal (Rees sinon)')
1125 	format('.false.',10x,
     &    'lopal : sec. distrib. de Opal (Rees sinon)')

 	  close(fic_datdeg)
c
1000    format(a)
1010 	format(i2,10x,' for iniflu : 1 use detailed grid, 2 coarse grid')
c1020 	format(a,3x,' Fichier d''entree des sections efficaces')
c1030 	format(a,3x,' Fichier de sortie des seff(E)')
c1040 	format(a,3x,' Fic. de sortie des seff diff. (Eprim-->E)')
1020 	format(a)
1030 	format(a)
1040 	format(a)
1050 	format('idess =-1 no plot.',/,
     &    'idess = 0 metacode created, without interactive session.',/,
     &    'idess = 1 metacode not created, but screened plots.',/,
     &    'idess = 2 metacode created, with interactive session.',/,
     &   2i5,5x,'(seff, seff detailled)')
1060 	format('iprint = 1 if print (whatever number otherwise)',/,
     &    8x,'sigel  siginel  rmatrix diff-sig-test ',/,4i10)
1130 	format(a7)
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

992     print*,' Cross-section file is in error. Status=',iost
     	stop

993     print*,' trouble writing crs file'
        stop

        end subroutine quelle_grille
