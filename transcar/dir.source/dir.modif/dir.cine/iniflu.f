c
c 
  	subroutine iniflu(npt,iyd,UTsec,z,glat,glong,f107,
     .                ap,chi,Ne,Te,Tj,indlim,jpreci,
     .   	      Nh,No,No2,Nn2,Nn,Tn,N_0,T_0,Po,Po2,Pn2,Ph,Pn,Heat,
     .                  nspec,knm,nen,nang,nango2,nalt,
     .                  ddeng,botE,centE,gmu,gwt,angzb,altkm,altcm,
     .                  dipang,smgdpa,
     .			fluxup,fluxdown,densneut,tneutre,tempexo,
     .                  albedo,chideg,denelc,temelc,temion,UT,hrloc,
     .                  ztop,zbot,colden,day,year)
c
c 	Initialisation des parametres necessaires au programme cinetique
c 	lorsque celui-ci est appele par le programme fluide.
c
c 	ENTREES DE INIFLU ISSUES DU PROGRAMME FLUIDE 
c 	--------------------------------------------
c 	Passees par liste d'appel :
c  	---------------------------
C       IYD     jour de l'annee sous la forme YYJJJ, 
c               avec YY=annee et JJJ=jour dans l'annee
C       SEC     temps universel                                 (s)
C       Z       tableau des altitudes                           (km)
C       GLAT    latitude geographique                           (degres)
C       GLONG   longitude geographique                          (degres)
C       F107    flux solaire moyen
C       F107    flux solaire instantane de la veille
C       AP      tableau des ap (convention MSIS86)
c 	CHI 	Angle solaire zenithal (radian)
c
C       Ne      tableau des concentrations electroniques
c               (unites normalisees)
C       Te      tableau des temperatures electroniques        
c               (unites normalisees)
C       Tj      tableau des temperatures de l'ion O+        
c               (unites normalisees)
C       NX      nombre d'altitudes resolues
c 	indlim 	= Nombre d'altitudes de calcul du transport fluide
c 	npt 	= dimension des tableaux utilises par fluide.
c       Nh(npt),No(npt),No2(npt),Nn2(npt),Nn(npt),Tn(npt) 
c 		= densites et temperature neutres de fluide.f
C       N_0     coefficient de normalisation des concentrations (cm-3)
C       T_0     coefficient de normalisation des temperatures   (K)
c 	Px 	Productions, simplement initialisees a 0 ici
c 	Heat 	Chauffage, simplement initialisees a 0 ici
c	jpreci = 0 if sun only, 1 if prec. only, 2 if both
c
c 	PARAMETRES INTERNES A CE SPGMME
c 	-------------------------------
c 	En particulier pour l'ecriture eventuelle de ELEC et NEUTRAL
c 	index(7)= index des elements neutres de fluide.f
c 	ndelta 	= nombre d'altitudes a rajouter pour descendre jusqu'a
c 		90 km
c 	ztrav 	= altitudes a rajouter pour descendre jusqu'a 90 km
c 	dens(8,100),fcdens(8),Ttrav,xmasdens,fctemp
c 	modatmos, alpha (ici, param. bidon),dne(id.),dte,dti,spfac
c 	moddene,modtemp,modcomp,z50,derivte,comp
c 	
c 	sw : tableau de switches pour msis95
c
c 	SORTIES DE INIFLU INTERNES A CE SPGMME, POUR LE TRANS. CINETIQUE
c 	----------------------------------------------------------------
c 	A terme, tout ceci devra etre mieux calcule...
c	nspec, knm, nen Emax, nang, nango2, stepnrj, anormnrj, ddeng,
c	botE, centE, pi, pideg, albedo
c
c 	SORTIES DE INIFLUX CALCULEES A PARTIR DES ENTREES DE FLUIDE
c 	------------------------------------------------------------
c	altkm, altcm(ialt) = altitudes
c 	zbot, ztop
c 	fluxup, fluxdown
c 	smgdpa
c 	tempexo
c 	densneut,colden, tneutre
c	UT, hrloc, year, nan, day, Apind, chideg, jpreci,
c	Flux_elec_int = Flux d'electrons integre (cm-2.s-1.sr-1)
c 	Flux_ener_int = Flux d'energie integre (kev.cm-2.s-1.sr-1)
c	 	  par liste d'appel
c
 	implicit none
c
	logical flg_err
	common/err/flg_err	
	
	integer jpreci0
        logical flgini,flgchi
        data flgini/.true./
        data flgchi/.false./

        include 'TRANSPORT.INC'
c
 	integer npt,ind,indref,i
 	integer nspec,knm,nen,nang,nango2,nalt,indlim,index(7)
        data index/3,4,2,7,1,8,5/
	integer in
	data in/1/

c
        real Nh(npt),No(npt),No2(npt),Nn2(npt),Nn(npt),Tn(npt)
 	real Ph(npt),Po(npt),Po2(npt),Pn2(npt),Pn(npt),Heat(npt)
c
        real N_0,T_0
 	real Eave,Flux_elec_int,E,Flux_ener_int,coef_flux,fluxE
c
 	real z(npt),Ne(npt),Te(npt),Tj(npt)
c
 	real stepnrj,anormnrj
 	real ddeng(nbren),botE(nbren),centE(nbren),Emax
 	real angzb(2*nbrango2),gmu(2*nbrango2),gwt(2*nbrango2)
 	real altkm(nbralt),altcm(nbralt),dipang(nbralt),smgdpa(nbralt)
        real fluxdown(nbren,nbrango2),fluxup(nbren,nbrango2)
 	real tneutre(nbralt),densneut(8,nbralt),colden(8,nbralt)
 	real denelc(nbralt),temelc(nbralt),temion(nbralt)
 	real dne(nbralt),dte(nbralt),dti(nbralt)
 	real tempexo,albedo
 	real UT,UTsec,hrloc,year,day,chi,chideg,f107(3),
     .		ap(7),glong,glat
 	real ztop,zbot
 	integer iyd,nan,jpreci
	integer ialt,ien,iang,isp,ii
 	real pi,d(8),t(2),sw(26)
 	integer ndelta
 	real ztrav(50),Ttrav(50)
  	real dens(8,50),fcdens(8),xmasdens(50),fctemp,comp(nbralt)
 	integer isotro,modatmos,moddene,modtemp,modcomp
c 	Pour IRI :
 	integer nziri,jmag,iddd
 	real ztopiri,zbotiri,zstepiri,ziri(5),Neiri(50),xind
        logical jf(12)
        real outiri1(11,50),outiri2(30),traviri(5)
        real date,xbid,ybid,zbid,bf,decbid
        real mgdpa(nbralt),field(nbralt)
c
        real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0,chi0
        real Bmag,dipangle,Enord,Eest,vperpnord,vperpest,vhorizon,vpara
	real ddp,Jtop
        integer ikp

        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi0,
     &                  Fe0,Ee0,Fi0,Ei0,
     &                  Bmag,dipangle,Enord,Eest,
     &                  vperpnord,vperpest,vhorizon,vpara,ddp,Jtop

     	real*8 dtu,dlon,dlat
     	real kp
	common/prec/coef_flux

 	real alpha,spfac,z50,derivte(nbralt)
 	integer icont(12),iprt(9)
c
 	data fcdens /1.,1.,1.,1.,1.,1.,1.,1./
 	data fctemp /1./,modatmos /2/,moddene /1/,modtemp/1/,modcomp/1/
 	data icont /0,0,0,0,0,0,0,0,0,0,0,0/
 	data iprt /1,1,0,0,0,0,0,0,0/
 	integer fanion
 	data fanion /0/
c       Les options sont les options standard de IRI
        data jf /.false.,.false.,.false.,.true.,.true.,.true.,.true.,
     .           .true.,.true.,.true.,.true.,.true./
c
        real PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),PMA(100,10)
        COMMON /PARM6/     PT,PD,PS,PDL,PTL,PMA


	integer ndeg_pn,ndeg_pf

c
c
c ====================================================================
c 	definition de la grille d'altitude
c ====================================================================
c
 	ndelta = 0
	do ialt=1,indlim
	  altkm(ialt)=z(indlim+1-ialt)
	  altcm(ialt)=1.e5*altkm(ialt)
	enddo
 	nalt = indlim
c 	Les altitudes doivent descendre a 90 km au moins.
c 	On genere eventuellement les altitudes supplementaires pour ca.
  	if(z(1).gt.90.)then
c 	  On va calculer une grille avec un pas moyen de 2 km
 	  ndelta = ifix((z(1) - 90.)/2.)
 	  ndelta = min(ndelta,50)
 	  if (ndelta.gt.1) then
 	    call gridexp (ndelta,90.,z(1),ztrav)
 	    do ialt = 2,ndelta
 	      nalt = nalt + 1
 	      altkm(nalt) = ztrav(ialt)
	      altcm(nalt)=1.e5*altkm(nalt)
 	    enddo
 	  endif
 	endif
	ztop=altkm(1)
	zbot=altkm(nalt)
c
c ====================================================================
c 	Chargement de l'atmosphere neutre
c ====================================================================
c
c
c 	Pour le moment, le code fluide ne se sert pas de he.
 	nspec = 4
c
	day=mod(float(iyd),1000.)
 	ikp = 2*max(0.,log((ap(1)+2.5)/4.))+.5
	UT=UTsec/3600.
	hrloc=UT+glong/15.
	year=float(iyd/1000)
	if (iyd.lt.1900000) then
	  year=year+1900.
	endif
	nan=ifix(year)
	tempexo=Tn(indlim)
c
	do ialt=1,indlim
	  densneut(1,ialt)=Nn2(indlim+1-ialt)/1.e6
	  densneut(2,ialt)=No2(indlim+1-ialt)/1.e6
	  densneut(3,ialt)=No (indlim+1-ialt)/1.e6
	  densneut(4,ialt)=Nh (indlim+1-ialt)/1.e6
	  tneutre(ialt)=Tn(indlim+1-ialt)
 	enddo
c 	Si le transport fluide ne descend pas assez bas, il fo modeliser
c 	les densites neutres en bas
  	if(z(1).gt.90. .and. ndelta.gt.1)then
 	  do ialt = 2,ndelta
c 	    ... msis90
  	    if (flg_err) then
  	      print*,iyd,UTsec,ztrav(ialt),glat,glong,
     &	hrloc,f107(3),f107(2),ap
            endif
  	
  	    call msis90(iyd,UTsec,ztrav(ialt),glat,glong,
     .			f107,ap,48,d,t,w)
c
 	    densneut(1,indlim+ialt-1)=d(3)/1.e6		! N2
 	    densneut(2,indlim+ialt-1)=d(4)/1.e6		! O2
 	    densneut(3,indlim+ialt-1)=d(2)/1.e6		! O
 	    densneut(4,indlim+ialt-1)=d(7)/1.e6		! H
 	    tneutre(indlim+ialt-1)=t(2)
  	  enddo
  	endif
c
c 	Calcul de la densite colonne (cm-2).
	call column(nbralt,nalt,altcm,nspec,densneut,colden)
c
c
c ====================================================================
c 	Chargement des conditions ionospheriques
c ====================================================================
c
c 	Attention ! Le programme fluide sort des te/ti a basse alt.
c 	qui peuvent etre < Tn ...
	do ialt=1,indlim
	  denelc(ialt)=N_0*Ne(indlim+1-ialt)
	  temelc(ialt)=max(T_0*Te(indlim+1-ialt),tneutre(ialt))
	  temion(ialt)=max(T_0*Tj(indlim+1-ialt),tneutre(ialt))
	enddo
c 	Si le transport fluide ne descend pas assez bas, il fo modeliser
c 	Ne, Te, Ti
c
        if(z(1).gt.90. .and. ndelta.gt.1)then
          nziri = 5
          ztopiri = altkm(indlim+1)
          zbotiri = altkm(nalt)
          zstepiri = (ztopiri-zbotiri)/4.
          call gridcst (nziri,zbotiri,ztopiri,ziri,traviri)
c         Calcule Ne ...
          jf(1) = .true.
c 	  IRI n'estime pas les Te Ti en dessous de 120 km.
          jf(2) = .false.  
c         Les latitudes et longitudes sont geographiques.
          jmag = 0
c         La date est en jours : negative pour IRI
          iddd = -ifix(day)
c         L'indice est f107 : negatif pour IRI
          xind = -f107(3)
          call iris12(jf,jmag,glat,glong,xind,iddd,hrloc,zbotiri,
     .              ztopiri,zstepiri,outiri1,outiri2)
c 	  Il reste a interpoler depuis la grille reguliere de iri
c 	  vers la grille exp de calcul :
 	  do ialt = 1,nziri
 	    traviri(ialt) = outiri1(1,ialt)/1.e+06
 	  enddo
 	  call intlin(nziri,ziri,traviri,ndelta,ztrav,Neiri)
 	  do ialt = 2,ndelta
	    denelc(indlim+ialt-1)=Neiri(ndelta+2-ialt)
 	  enddo
c 	
c 	  iri ne modelise pas les temperatures en dessous de 120 km-->
c 	  on prend Te=Ti=Tn
 	  do ialt = indlim+1,nalt
	    temelc(ialt)=tneutre(ialt)
	    temion(ialt)=tneutre(ialt)
 	  enddo
        endif
c
c 	Initialisation a 0 des productions
 	call zeroit(Pn2,npt)
 	call zeroit(Po2,npt)
 	call zeroit(Po,npt)
 	call zeroit(Ph,npt)
 	call zeroit(Pn,npt)
 	call zeroit(Heat,npt)
c
c 
c ====================================================================
c 	definition des grilles de calcul
c ====================================================================
c
c 	Numero d'identification du run
	knm=ifix(UTsec)
c 	S'il y a des precipitations, elles sont isotropes
 	isotro = 1
 	
c 	Si il n'y a pas de precipitations, on garde la 1ere grille.
 	if (jpreci.eq.0)then
 	  Emax = 300.
 	  call quelle_grille(Emax,nen,centE,botE,ddeng,
     .		nang,nango2,gmu,gwt,angzb)
 	else
c
c 	  Il faut connaitre l'energie caracteristique des precipitations
c 	  pour definir sur quelle grille il faut travailler. On  
c 	  calcule le flux donne par le modele statistique de Hardy. 
c 	  Ce programme donne  
c 		le flux total Flux_elec_int EN ELECTRONS/CM2.S.SR 
c 		le flux total d'energie Flux_ener_int EN KEV/CM2.S.SR 
c
c	  call prec_time(iyd,UTsec,Eave,Flux_ener_int)
	  dlon=dble(15.*tmag)
	  dlat=dble(latmag)
	  dtu=dble(UTsec)
	  kp=ikp/3.
	  call precipitation(iyd,dtu,kp,dlon,dlat,Eave,Fluxe)
	  Flux_ener_int=Fluxe

c 	  Pour une maxwellienne, 99% du signal est contenu avant 7Eo
          emax =10.*Eave
CME	  correction CME
	  Emax=min(40000.,10.*Eave)
c	Emax=40000.
c
c 	  write(6,*)'Flux_ener_int,Eave',Flux_ener_int,Eave
 	  call quelle_grille(Emax,nen,centE,botE,ddeng,
     .		nang,nango2,gmu,gwt,angzb)
c
c
c         ============================================================
c         Initialisation du flux precipite
c         ============================================================
c
c 	  On a suppose donc (a l'encontre de ce que precise Hardy) que
c 	  le flux est maxwellien

c	if (latmag.ge.0.) then
c	if (latmag.le.0.) then
c	  isotro=1
  	   call fluxkappa(nango2,nen,centE,isotro,
     .                  gmu,fluxdown,fluxup,Eave,Flux_ener_int)
          call flux_integre(nango2,nen,centE,ddeng,
     .                    gmu,gwt,fluxdown,fluxup,
     .			Eave,Flux_ener_int)
	  Fe0=Flux_ener_int*6.2832
	  Ee0=Eave
c
c	else
c	  isotro=1
c  	  call inmaxwl(Flux_ener_int,Eave,nango2,nen,centE,
c     .			isotro,gmu,fluxdown,fluxup)	
c          call normflux(Flux_ener_int,nango2,nen,centE,ddeng,gmu,gwt,
c     .			fluxdown,fluxup)
c	  Fe0=Flux_ener_int*6.2832
c	  Ee0=Eave
	endif

	Fe0=Fe0*1.6e-12
c 	endif
c
c
c ====================================================================
c       Initialisation d'un flux precipite de perturbation
c ====================================================================
c
c 	Connaissant le flux de particules integre [cm-2.s-1.sr-1]
c       et le flux d'energie integre [keV.cm-2.s-1.sr-1], on a le 
c 	choix entre diverses possibilites. 
c 	Il faut appeler un programme qui en fonction de ces 
c 	parametre estime un flux. En routine, il existe un flux 
c 	maxwellien, gaussien ou un flux monoenergitique (dirac)
c 	On les appelle par
c
c 	call inmaxwl(Flux_ener_int,Eave,nango2,nen,centE,
c    .			isotro,gmu,fluxdown,fluxup)	
c
c       call ingauss (Flux_ener_int,Eave,nango2,nen,centE,
c    .			isotro,gmu,fluxdown,fluxup)
c
c       call indirac (Flux_ener_int,Eave,nango2,nen,centE,
c    .			isotro,gmu,fluxdown,fluxup)
c
c       isotro = parametre d'isotropie du flux precipite :
c       0 = distribution gaussienne 
c       1 = flux isotrope
c       2 = dirac
c
c 	La discretisation impose dans tous les cas de normaliser a 
c 	Flux_ener_int.
c       call normflux(Flux_ener_int,nango2,nen,centE,ddeng,gmu,gwt,
c    .		fluxdown,fluxup)
c
c 	Pour mettre des perturbations, on peut changer brutalement la 
c	forme et l'intensite du flux. Il suffit d'initialiser 
c 	Flux_ener_int,Eave,isotro et de retirer les c de commentaire
c 	a l'appel souhaite. Il faut eventuellement redefinir la grille
c 	d'energie.
c 	
c 	Remarque : 1 Erg = 6.25e+11 eV...
c
c ---- 	ICI COMMENCENT LES CHANGEMENTS EVENTUELS POUR UNE PERTURBATION
c
c 	emax = 400.
c	call quelle_grille(Emax,nen,centE,botE,ddeng,
c    .		nang,nango2,gmu,gwt,angzb)
c	isotro = 1
c	Flux_ener_int = 6.25e+11 * 3.
c 	Eave = 200.
c
c 	call inmaxwl(Flux_ener_int,Eave,nango2,nen,centE,
c    .			isotro,gmu,fluxdown,fluxup)	
c       call normflux(Flux_ener_int,nango2,nen,centE,ddeng,gmu,gwt,
c    .			fluxdown,fluxup)
c
c	isotro = 1
c       call ingauss (Flux_ener_int,Eave,nango2,nen,centE,
c    .			isotro,gmu,fluxdown,fluxup)
c       call normflux(Flux_ener_int,nango2,nen,centE,ddeng,gmu,gwt,
c    .			fluxdown,fluxup)
c
c	isotro = 1
c       call indirac (Flux_ener_int,Eave,nango2,nen,centE,
c    .			isotro,gmu,fluxdown,fluxup)
c       call normflux(Flux_ener_int,nango2,nen,centE,ddeng,gmu,gwt,
c    .			fluxdown,fluxup)
c
c ----	FIN DE CES CHANGEMENTS.
c
c
c ====================================================================
c 	Chargement des conditions du run
c ====================================================================
c
c
c       Calcul du champ : mgdpa = magnetic dip angle (radians)
        date=float(nan)+day/100.
        if(nan.lt.100)date=1900.+float(nan)+day/100.
        do ialt = 1,nalt
          call bfield(date,altkm(ialt),glat,glong,xbid,ybid,zbid,bf,
     .              mgdpa(ialt),decbid)
          mgdpa(ialt)=abs(mgdpa(ialt))
          smgdpa(ialt)=sin(mgdpa(ialt))
c 	  La valeur du champ pourra dans le futur etre a garder...
c         Conversion de B, des gammas en Tesla
c         field(ialt)=bf*1.0e-09
 	enddo
	albedo=1.
c
        chideg = chi * 180.0 / (4.*atan(1.))

        if (chideg.ge.108.) then
          if (.not.flgchi) jpreci0=jpreci
          flgchi=.true.
          if (jpreci.eq.2) jpreci=1
          if (jpreci.eq.6) jpreci=4
        else
          if (flgchi) then
            jpreci=jpreci0
            flgchi=.false.
          else
            jpreci0=jpreci
          endif
        endif

c        if (Fluxe.lt.1.e-30) then
c          if (.not.flgchi) jpreci0=jpreci
c          flgchi=.true.
c          if (jpreci.eq.2) jpreci=1
c          if (jpreci.eq.6) jpreci=4
c        else
c          if (flgchi) then
c            jpreci=jpreci0
c            flgchi=.false.
c          else
c            jpreci0=jpreci
c          endif
c        endif


c
c ====================================================================
c 	Eventuellement, ecriture pour verifications.
c ====================================================================
c
c 	Pour mettre off ces ecritures, mettre iprt(1) et iprt (2) a 
c 	0 dans la ligne data en debut de ce sous programme
 	if (iprt(1).eq.1 .or. iprt(2).eq.1) then
 	  fanion = fanion + 1
 	  do ialt = 1,nalt
 	    dne (ialt) = denelc(ialt)/10.
 	    dte (ialt) = temelc(ialt)/10.
 	    dti (ialt) = temion(ialt)/10.
 	  enddo
          call ecr(nspec,nalt,zbot,ztop,UT,hrloc,day,nan,
     .      jpreci,tempexo,f107(2),f107(3),ap(1),fctemp,fcdens,glat,
     .      glong,modatmos,albedo,altkm,nen,botE,centE,ddeng,knm,eave,
     .      alpha,nang,angzb,gmu,gwt,fluxdown,fluxup,denelc,dne,
     .      temelc,dte,temion,dti,derivte,comp,z50,densneut,tneutre,
     .      xmasdens,colden,mgdpa,smgdpa,chideg,icont,iprt,spfac,
     .      moddene,modtemp,modcomp)
 	endif
c
 	return
 	end
c
c------------------------ column -------------------------------
c 
	subroutine column(nbralt,nalt,altcm,nspec,densneut,colden)
*
*	Calculates the column density HTAU
*
 	implicit none
c
c 	Parametres d'entree/sortie
 	integer nbralt,nalt,nspec
	real densneut(8,nbralt),altcm(nbralt),colden(8,nbralt)
c 	Parametres internes
 	real epsilon,d
 	integer i1,i2,i3,i4,isp,ialt
c
	epsilon=1.e-10
	if(altcm(nalt).lt.altcm(1)) then
		i1=1
		i2=2
		i3=nalt
		i4=1
	else
		i1=nalt
		i2=nalt-1
		i3=1
		i4=-1
	end if
 	do isp=1,nspec
	  d=0.
	  if(densneut(isp,i2).gt.0..and.densneut(isp,i1).gt.0.)
     .		d=abs(log(densneut(isp,i2)/densneut(isp,i1)))
	  if(abs(d).ge.epsilon) 
     .		colden(isp,i1)=densneut(isp,i1)/d*(altcm(i1)-altcm(i2))
	  do ialt=i2,i3,i4
	  d=0.
	    if(densneut(isp,ialt-i4).gt.0..and.densneut(isp,ialt).gt.0.)
     .		d=abs(log(densneut(isp,ialt-i4)/densneut(isp,ialt)))
	    if(abs(d).ge.epsilon) then
	      colden(isp,ialt)=colden(isp,ialt-i4)+
     .	      		abs((altcm(ialt)-altcm(ialt-i4))/d*
     .			(densneut(isp,ialt)-densneut(isp,ialt-i4)))
	    else
	      colden(isp,ialt)=colden(isp,ialt-i4)+
     .			densneut(isp,ialt)*(altcm(ialt-i4)-altcm(ialt))
	    end if
 	  enddo
 	enddo
c
	return
	end
