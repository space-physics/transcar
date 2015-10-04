c
c-------------------------- ecr ---------------------------------
c
        subroutine ecr(nspec,nalt,zbot,ztop,UT,hrloc,day,nan,
     &      jpreci,tempexo,f107,f107a,Apind,fctemp,fcdens,glat,
     &      glong,modatmos,albedo,altkm,nen,botE,centE,ddeng,knm,eave,
     &      alpha,nang,angzb,gmu,gwt,fluxdown,fluxup,denelc,dne,
     &      temelc,dte,temion,dti,derivte,comp,z50,densneut,tneutre,
     &      xmasdens,colden,mgdpa,smgdpa,chideg,icont,iprt,spfac,
     &      moddene,modtemp,modcomp)
c
 	    include 'TRANSPORT.INC'
c
 	    real altkm(nbralt),denelc(nbralt),dne(nbralt),temelc(nbralt),
     &		dte(nbralt),temion(nbralt),dti(nbralt),derivte(nbralt),
     &		tneutre(nbralt),gwt(nbrang),comp(nbralt)
c
 	    real angzb(nbrang),gmu(nbrang)
c
     	real botE(nbren),centE(nbren),ddeng(nbren)
c
 	    real fcdens(8),densneut(8,nbralt),xmasdens(nbralt)
      	real colden(8,nbralt),mgdpa(nbralt),smgdpa(nbralt)
c
 	    real fluxup(nbren,nbrango2),fluxdown(nbren,nbrango2)

     	integer moddene,modtemp,md2
     	real chideg
c
     	integer icont(12),iprt(9)
     	real travalt(nbralt)
c
     	md2 = nang/2
     	pi = 4.*atan(1.)
     	pio180 = pi/180.
c
c-------------------- neutraldata ----------------------------
c
 	    if (iprt(1).eq.1)then
c 	write(6,*)'Writting NEUTRAL'
c	write(6,*)'Writting NEUTRAL                             [A'
 3000 format (
     &'      ----- neutraldata, reecrit apres iniflu.f  ---',/,
     &'      Nombre d''especes a considerer (nspec), Id #',/,i1,10x,i10)
3016 	format(
     &'       altitude grid parameters', /,
     &'      nalt: # altitude grid points pour photo ',/,
     &'      zbot: minimum altitude. ztop: maximum altitude.   ',/,
     &'      zbot,ztop definissent la grille de travail global')
 3001 format (i10,40x,'nalt')
 3002 format (3f10.2,20x,'zbot/ztop')
 3003 format ( '1  neutral atmosphere parameters',/,
     &'2  ----------------------------- ',/,
     &'3     hrloc: Local hour angle, UT:Universal time',/,
     &'4     day: julian day (1 to 365); Annee.',/,
     &'5     jpreci=0sol,1el,2sol/el,3prot,4prot/el,5sol/prot,6tout',/,
     &'6     tempexo: exospheric temperature (calculee si =1.)',/, 
     &'7     f107/bar: solar 10.7-cm line output (this day/average).',/,
     &'8     chideg : solar zenith angle (degrees)',/,
     &'9     Ap: Magnetic index')
 3004 format ('      fcdens,fctemp:Density and Tn scaling factors',/,
     &'1     glat, glong:  Lat. (-90 to +90). Long. (degree)',/,
     &'2     modatmos : 1 = Atmosphere imposee ',/,
     &'3                2 = MSIS,iterations sur f107',/,
     &'4                3 = MSIS,iterations sur f107 et Ap',/,
     &'5     Albedo at the lowest altitude.')
3005  format(1f8.4,1f8.1,1i8,1i8,1f10.4,3x,'hr loc/day/annee/jpreci/UT')
 3006 format (f8.2,3(f8.1),f10.2,3x,'tempexo/f107/f107a/Ap/chideg')
 3007 format (9f6.2,5x,'fctemp/fcdens')
 3008 format (2f10.2,i10,f8.2,7x,'lat,long,modatmos,albedo')
 3011 format('1  no en. (nen) and energies for transport:',/,
     & '2  in eV (increasing order):')
 3012 format(i10)
 3013 format(5f10.2)
 3014 format('alt.(z) for computation of production(s)(1ry and 2ry)')
 3015 format(5f10.2)
3050    format ('Neutral temperature')
3055    format (5f14.6)
3057    format(1i3,20x,'Number of neutral species in atmosphere')
3060    format (a5, ' neutral density [cm-3]')
3065    format (5(1pe14.6))
3066    format (a5,' Column density (Warning! non divided by ',
     &    'sin(magn.dip angle) ')
c
      open (ineutr, file=data_path(1:lpath_data)
     &                         //'dir.cine/NEUTRAL_VERIF')
      rewind ineutr
      write (ineutr,3000)nspec,knm
      write (ineutr,3016)
      write (ineutr,3001) nalt
      write (ineutr,3002) zbot, ztop
      write (ineutr,3003)
      write (ineutr,3004)
      write (ineutr,3005) hrloc,day,nan,jpreci,UT
      write (ineutr,3006) tempexo, f107, f107a,Apind,chideg
      write (ineutr,3007) fctemp,(fcdens(i),i=1,8)
      write (ineutr,3008) glat,glong,modatmos,albedo
      write (ineutr,3014)
      write (ineutr,3012)nalt
      write (ineutr,3015)(altkm(i),i=1,nalt)
      write(ineutr,3050)
      write(ineutr,3055)(tneutre(i),i=1,nalt)
      write(ineutr,3057)nspec
      do jsp=1,nspec
 	  write(ineutr,3060)specie(jsp)
        write(ineutr,3065)(densneut(jsp,ialt),ialt=1,nalt)
 	  write(ineutr,3066)specie(jsp)
        write(ineutr,3065)(colden(jsp,ialt),ialt=1,nalt)
      enddo
c
      close(ineutr)
 	  endif
c
c----------------- elecprofil ---------------------------
c
 	  if (iprt(2).eq.1)then
c	write(6,*)'Writting ELEC                                [A'
c 	write(6,*)'Writting ELEC'
 4000 format('  ----- elecprofil apres iniflu.f --------',/)
 4001 format(i10,10x,'Nombre d''energies (nen)')
 4035 format(i10,10x,f10.7,9x,'nen, facteur de croissance')
 4002 format(5(1pe14.7))
 4003 format('1  ------  flux descendant (cm-2.s-1.eV-1.sr-1)',/,
     &'2   angle, cos(angle),gaussian weight and fluxes')
 4004 format(i10,10x,'Nombre d''angles (nang)')
 4006 format('1  pitch angles ')
 4007 format(5(1pe10.2))
 4008 format ('1  ------ densite electronique ',
     &      /,'2     ambient electron parameters ',
     &      /,'3     ne: number of alt.',
     &      /,'4     alt: altitudes. (km)       ',
     &      /,'5     denelc: number densities.(/cm3)',
     &      /,'6 ') 
 4009 format ('1  ------ densite electronique ',
     &      /,'2     ambient electron parameters ',
     &      /,'3     ne: number of alt.',
     &      /,'4     alt: altitudes. (km)       ',
     &      /,'5     denelc: number densities.(/cm3) (from IRI)',
     &      /,'6 ') 
 4010 format(5f14.6)
 4011 format(5(1pe14.6))
 4012 format(i10)
 4013 format('1  ------  flux montant (cm-2.s-1.eV-1.sr-1)',/,
     &  '2   angle, cos(angle),gaussian weight and fluxes')
 4014 format('Compositions [O+]/Ne (%)Mesures EISCAT. z50 =',f7.2,'km')
 4015 format(1f20.4,2f20.13)
 4016 format('Angles en ordre croissant.{0->90=fd},{90->180=fu}')
 4017 format ('Id. No (knm), mean E(maxw) (0 if no),     alpha'
     &       ,/,47('-'),/,i10,10x,f15.2,f10.2)
 4018 format ('Energie (eV)(en ordre decroissant) ')
 4030 format ('E(eV) (ordre decr.), loi de puissance')
 4019 format (4f18.13)
 4020 format ('Corresponding cosines')
 4021 format ('Corresponding gaussian weight factors')
 4022 format (4f18.4)
 4023 format (5f13.5)
 4024 format ('% O+/Ne Modele EISCAT-Ch.L. (dependant du temps). z50 =',
     &  f7.2,'km')
 4025 format ('% O+/Ne Modele Millestone')
c
      open (ielec,file=data_path(1:lpath_data)
     &                       //'dir.cine/ELEC_VERIF')
      rewind ielec
      write (ielec,4000)
      write (ielec,4017)knm,eave,alpha
      if(icont(3).eq.5)then
 	  write(ielec,4030)
        write (ielec,4035)nen,spfac
      else
        write (ielec,4018)
        write (ielec,4001)nen
      endif
      write(ielec,*)'Bottom energy'
      write (ielec,4002) (botE(nen+1-ien),ien=1,nen)
      write(ielec,*)'Center energy'
      write (ielec,4002) (centE(nen+1-ien),ien=1,nen)
      write(ielec,*)'Energy grid width'
      write (ielec,4002) (ddeng(nen+1-ien),ien=1,nen)
      write (ielec,4016)
      write (ielec,4004)nang
      write (ielec,4022) (angzb(iang),iang=1,nang)
      write (ielec,4020) 
      write (ielec,4019) (gmu(iang),iang=1,nang)
      write (ielec,4021) 
      write (ielec,4019) (gwt(iang),iang=1,nang)
      write (ielec,4003)
      do iang=1,md2
	  write (ielec,4015) angzb(iang),gmu(iang),gwt(iang)
        write (ielec,4002) (fluxdown(nen+1-ien,iang), ien = 1,nen)
      enddo
      write (ielec,4013)
      do iang=1,md2
        write (ielec,4015)angzb(md2+iang),gmu(md2+iang),gwt(md2+iang)
        write (ielec,4002) (fluxup(nen+1-ien,iang), ien = 1,nen)
      enddo
      if(moddene.eq.2)then
 	  write (ielec,4009)
      else
 	  write (ielec,4008)
      endif
      write (ielec,4012) nalt
c     In order to avoid division by 0, it is necessary to reduce
c 	the alt. extrema.
      travalt(1)=altkm(1)+0.01
     	do 10 ialt=2,nalt-1
     	  travalt(ialt)=altkm(ialt)
10    continue
      travalt(nalt)=altkm(nalt)-0.01
      write (ielec,4010) (travalt(i), i = 1,nalt)
      write (ielec,4011) (denelc(i), i = 1,nalt)
c     write (ielec,*)'Erreurs sur les densites electroniques'
c     write (ielec,4011) (dne(i), i = 1,nalt)
      if (modtemp.eq.2)then
     	write (ielec,*)'Temperatures electroniques from IRI'
          else
     	write (ielec,*)'Temperatures electroniques'
      endif
      write (ielec,4010) (temelc(i), i = 1,nalt)
c     write (ielec,*)'Erreurs sur les temperatures electroniques'
c     write (ielec,4011) (dte(i), i = 1,nalt)
      if (modtemp.eq.2)then
     	write (ielec,*)'Temperatures ioniques from IRI'
          else
     	write (ielec,*)'Temperatures ioniques'
      endif
      write (ielec,4010) (temion(i), i = 1,nalt)
c     write (ielec,*)'Erreurs sur les temperatures ioniques'
c     write (ielec,4011) (dti(i), i = 1,nalt)
      write (ielec,*)'Derivee de la temperature electronique'
      write (ielec,4011) (derivte(i), i = 1,nalt)
      if(modcomp.eq.0)write (ielec,4014)z50
      if(modcomp.eq.1)write (ielec,4024)z50
      if(modcomp.eq.2)write (ielec,4025)
      write (ielec,4002) (comp(i),i=1,nalt)
      write(ielec,*)'Magnetic deep angle (mgdpa), [degres]'
      write (ielec,4023) (mgdpa(i)/pio180,i=1,nalt)
      write(ielec,*)'sin{Magnetic deep angle (mgdpa)}'
      write (ielec,4023) (smgdpa(i),i=1,nalt)
 
      close(ielec)
 	  endif
c
      	return
      	end
c
c----------------------------------------------------------------------
c
        subroutine ecrit_DATDEG
c
 	    implicit none

        include 'TRANSPORT.INC'
c
c 	Le but de ce sous programme est de palier le probleme suivant :
c 	lorsqu'on appelle transcar, il laisse dans DATDEG la trace du
c 	fichier d'entree des sections efficaces utilise. Si (comme ca
c 	arrive souvent), on appelle transsolo ensuite, avec des 
c 	conditions qui n'ont rien a voir, on ecrase les sections
c 	efficaces pour en recalculer d'autres. Pour eviter ca, ce
c 	sous programme, a n'appeler qu'apres transsolo (ce serait meme
c 	mieux a la fin de transcar) reecrit dans DATDEG des noms de
c 	fichier de sections efficaces differents.

        integer i
        character*70 chaine(18),crs,rdt
1000 	format(a)
c
        open(fic_datdeg,file=data_path(1:lpath_data)
     &                             //'dir.cine/DATDEG')
	    rewind(fic_datdeg)
     	do i = 1,18
     	  read(fic_datdeg,1000)chaine(i)
     	enddo
     	close(fic_datdeg)
c
     	chaine(3) = 'dir.cine/dir.seff/crs'
     	chaine(4) = 'dir.cine/dir.seff/rdt'
c
        open(fic_datdeg,file=data_path(1:lpath_data)
     &                             //'dir.cine/DATDEG')
	    rewind(fic_datdeg)
     	do i = 1,18
     	  write(fic_datdeg,1000)chaine(i)
     	enddo
     	close(fic_datdeg)
c
 	    return
    	end
