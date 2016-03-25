      subroutine transelec(npt,iyd,UTsec,glat,glong,stl,f107,
     .          ap,chi,Ne,Te,Tj,nx,Nh,No,No2,Nn2,Nn,Tn,indlim,jpreci0,
     .          N_0,T_0,kiappel,zlim,zlim_1,z,Heat,Ph,Po,Po2,Pn2,Pn,
     .          Ne_sup,courant_sup,Te_sup,Chaleur_sup)
c
c 	Ce programme est un driver.
c       This program is a driver.
c
c 	- Lit toutes les donnees necessaires a une execution.
c         Reads all the required data to an execution.
c
c		--> si kiappel = 1, ce sont les donnees contenue dans ELEC et NEUTRAL pour un run unique du programme.
c		    if kiappel = 1, the data are contained in ELEC and NEUTRAL for a single run of the program.
c
c 		--> si kiappel = 2, ce sont les donnees necessaires au transport fluide.
c		    if kiappel = 2, these are the necessary data transmission fluid
c
c 	- appelle le programme de degradation d'energie.
c	  program called Energy degradation.
c
c 	- appelle si necessaire celui de calcul de la photoionisation.
c         call if necessary the calculation photoionization
c 		--> si jpreci = 0, 2, 5 ou 6. if jpreci = 0, 2, 5 or 6.
c
c 	- appelle le code de transport
c	  called the transport code
c
c 	- eventuellement, interpolle sur les grilles du transfluide
c	  possibly interpolated on the fence transfluide
c
c       The necessary data for a code to another by the intermediary are non-formatted files.
c       This method is costly in time MAYBE, but very practical.
c	Each program runs with its own data files that says what to print, draw, or find the data of cross sections, solar flux ...
c	Again, this way of proceeding can be costly in time, but more comfortable to use than go through commons which regulate the switches.
c	DO NOT MODIFY THIS WAY TO MAKE if you want an easy followed this 'black box' of transport (in particular releases).
c
c 	Les donnees necessaires d'un code a l'autre se font par l'inter
c 	mediaire de fichiers non formattes. Cette methode est peut
c 	etre couteuse en temps, mais tres pratique.
c 	Chaque programme tourne avec son propre fichiers de donnees
c	qui dit quoi imprimer, dessiner, ou trouver les donnees de
c	sections efficaces, de flux solaire... Encore une fois, cette
c 	facon de proceder peut etre couteuse en temps, mais beaucoup
c 	plus commode a l'utilisation que passer par des commons qui
c 	reglent les switches.
c 	NE PAS MODIFIER CETTE FACON DE FAIRE si on veut un suivit
c 	facile de cette 'boite noire' du transport (en particulier
c 	des releases).
c
c ====================================================================
c Computational parameters :
c ====================================================================
c nbren    =    IPN: max. length of the energy grid
c nen      =    IPN: actual length of the energy grid
c nbralt   =    IPM:(MMAX) length of altitude grid (number of layers +1)
c nalt     =    actual # of altitudes
c nbrang   =    maximum number of streams
c nbrango2 =    one half of the maximum number of streams
c nang     =    actual number of stream (must be less or equal 16;
c               should preferably be a multiple of two, if more than 16
c               streams are required,recompile all subroutines with
c               larger arrays)
c nango2   =    one half of the actual number of streams
c nbrsp    =    max # of species (1=N2 , 2=O2 , 3=O , 4=H , 5=He)
c nspec    =    actual # of species dans l'ordre(imperatif) N2,O2,O,H,He
c nbrionst =    number of ionized states
c nbrexc   =    maximum number of excitation states
c
c       Les angles descendants vont de 0 a
c       90 degres et les angles montant de 90 a 180 degres.
c       gwt : gaussian weigths (         1 ...    nango2 : up,
c                               nbrango2+1 ...  2*nango2 : down)
c       gmu : gaussian angles
c
c       The up and downward intensity, flux etc is defined in terms of
c       positive and negative mu. My definition is, however, independent
c       of the direction of the magnetic field and thus a bit arbitrary.
c       On the northern hemisphere the field is pointing down, and
c       thus a pitch angle of 0 degrees (mu=1) would be down and a pitch
c       angle of 180 degrees (mu=-1) would be pointing up.  However, I
c       define the pitch angle 90-180 (mu=[0,-1]) to be down, and the
c       pitch angle of 0-90 (mu=[0,+1]) up.  This way the transport code
c       doesn't care about the hemisphere anymore.
c       The program is built after the Northern hemisphere assumption,
c       i.e. {0,90} degrees --> downward, and {90,180} degrees -->
c       upward.
c
c ====================================================================
c
        implicit none
c
        integer npt
c
        include 'TRANSPORT.INC'
c
        integer nspec,nalt,jpreci0,jpreci,modatmos,neutspe
        real zbot,ztop,hrloc,ut,year,tempexo,f107(3),Apind,day,
     . 	 	glat,glong,albedo,altcm(nbralt),altkm(nbralt),
     .		tneutre(nbralt),densneut(8,nbralt),
     .     	colden(8,nbralt)
      integer knm,nang,nango2,nen
      real botE(nbren),centE(nbren),ddeng(nbren)
      real angzb(2*nbrango2),gmu(2*nbrango2),gwt(2*nbrango2)
      real fluxdown(nbren,nbrango2),fluxup(nbren,nbrango2)
      real denelc(nbralt),temelc(nbralt),dipang(nbralt),
     .	       smgdpa(nbralt),temion(nbralt),chideg
        real prodiontot(nbrsp*2,nbralt),chaufelec(nbralt)
        real Ne_supra(nbralt),courant_supra(nbralt),Te_supra(nbralt),
     .          Chaleur_supra(nbralt)
c
           integer i,j,iyd,nx
        real z(npt),Ne(npt),Te(npt),ap(7),stl
        real Tj(npt),chi,pideg,dTinf,stepnrj,anormnrj
        real UTsec,dz
        real HQe_1,HPo_1,HPo2_1,HPn2_1,HPh_1,flux
        real Eave
         real zlim,zlim_1

        real Nh(npt),No(npt),No2(npt),Nn2(npt),Nn(npt),Tn(npt)
        real Ne_sup(npt),courant_sup(npt),Te_sup(npt),Chaleur_sup(npt)

        !Variables for photodissociation of O2
        integer npt2
        parameter (npt2=500)
        real phdisso2(500),pfluxsr(8,201)
        real Po1sdisso2(npt2)
        common/photodissociation/phdisso2,Po1sdisso2
        common/kininds/nalt,nspec,nen
        !-----MZ

         integer kiappel
        integer ien,ialt,iang,isp

        real coef
c
        integer indlim
c
        real Ph(npt),Po(npt),Po2(npt),Pn2(npt),Pn(npt)
        real Heat(npt)

        real N_0,T_0
        common/exo/dTinf

c
       data pideg/57.29578/
c
1000   format (a)
c
c ====================================================================
c 	Lecture des parametres d'entree pour le transport
c ====================================================================
c ---- 	Identification des grilles
c 	knm 	= identification du run(sert aux impressions ou dessins)
c 	zbot 	= altitude min (km)
c 	ztop 	= altitude max (km)
c 	altkm 	= grille d'altitudes, sens decroissant (km)
c 	alt 	= grille d'altitudes, sens decroissant (cm)
c 	botE 	= Energie basse de la grille d'energie (eV).
c 		  L'energie la plus basse doit etre = a 0 eV.
c 		  Sert au rangement des photoelectrons dans la grille.
c 		  Les energies doivent etre en ordre croissant.
c 	centE 	= Energie centrale de la grille d'energie (eV).
c 		  L'energie la plus basse doit etre <= a 1 eV.
c 	ddeng	= largeur de chaque grille d'energie.
c 	angzb	= grille d'angles en degres (en nombre nang)
c 		  Doit etre une grille gaussienne pour les integrations
c 		  Ordre croissant de 0 a 180 degres.
c 		  {0->90=down},{90->180=up}
c 	gmu	= cosinus (grille d'angles) (en nombre nang)
c 	gwt	= Poids gaussiens de la grille d'angles (en nombre nang)
c ----	Parametres physiques
c 	hrloc	= heure decimale locale
c 	ut	= heure universelle
c 	day	= jour
c 	year 	= annee
c 	tempexo = temperature exospherique
c 	f107, f107a = ...
c 	Apind	= indice Ap
c 	chideg  = angle solaire zenithal en degres (inutile si on
c 		  n'appelle pas la photoionisation)
c 	glat  	= latitude (Tromso = 69.5)
c 	glong  	= longitude (Tromso = 19.2)
c 	tneutre = temperature neutre sur la grille altkm
c 	densneut= densite des neutres en nombre neutspe (cm-3)
c 		  ATTENTION A LA DIMENSION ! densneut est dimensionne
c 		  en (8,nbralt). 8 correspond au nombre max de neutres
c 		  dans msis, plus NO que je peux calculer par ailleurs.
c 		  Le transport se calcul sur un nombre max. de 5
c 		  especes (= nbrsp), mais pour d'autres applications,
c 		  il est important de conserver ce dimensionnement.
c 	colden  = densite colonne des neutres en nombre neutspe (cm-2)
c 		  ATTENTION A LA DIMENSION ! colden est dimensionne
c 		  en (8,nbralt) comme densneut.
c 	denelc	= densite electronique (sert pour le calcul de la
c 		  friction, cm-3)
c 	temelc	= temperature electronique (sert pour le calcul de la
c 		  friction, K)
c 	temion	= temperature ionique (sert pour le calcul des
c 	 	  conductivites (dans les frequences de collisions).
c 	dipang	= Magnetic dip angle (xmgdpa), [degres]
c 	smgdpa 	= sin{Magnetic deep angle (xmgdpa)}
c 	fluxdown= Flux descendant (cm-2.s-1.eV-1.sr-1), dans le sens
c 	 	  des energies croissantes.
c 	fluxup  = Flux montant (cm-2.s-1.eV-1.sr-1). Ne sert pas au
c 		  calcul lui meme, mais pour un test de comparaison
c 		  entre cette mesure eventuelle et le calcul
c 	Pour fluxdown et fluxup, on imagine un satellite qui tourne
c 	et qui mesure en fluxdown(1) le flux descendant //B, puis
c 	en fluxdown(nang/2) le flux descendant proche de perp. B,
c 	puis en fluxup(1) le flux montant proche de perp. B et enfin
c 	en fluxup(nang/2) le flux montant //B
c 	L'indice des angles (ainsi que leurs cosinus et poids) va
c 	pendant ce temps de 1 a nang), et les angles de pres de 0 a
c 	pres de 180 degres.
c 	Le programme de transport remet ces tableau dans l'ordre qui
c 	lui est digerable (different). Cela prend un peu de temps,
c 	mais aide dans le cas (comme pour Viking) ou on a a lire
c 	directement des donnees satellites.
c ----	Conditions du run
c 	jpreci = 0 if sun only
c 	       = 1 if electron precipitation only
c 	       = 2 if sun + electron precipitations
c	       = 3 if proton precipitation only
c	       = 4 if proton and electron precipitations
c	       = 5 if surn + proton precipitations
c	       = 6 if proton, electron precipitations + sun
c 	albedo  = charge la condition limite basse (la condition haute
c 		  etant le flux precipite). Dit quelle proportion de
c 		  particules remontent de la couche au dessous de l'
c 		  altitude la plus basse. Quand cette altitude zbot
c 		  est en region E, albedo peut etre n'importe quoi entre
c 		  0 et 1 (tout est local). Mais si zbot est en region
c 		  F, il faut faire differents essais pour en voir l'
c 		  influence.
c 		  Albedo conseille : 1
c 	 	  Si le flux descendant par l'altitude la plus basse est
c 		  trop eleve, il apparait un warning.
c
c 	Sorties du transport :
c 	----------------------
c
c 	prodiontot(sp,z) =  prodionprim(sp,z) + prodionsec(sp,z)
c                       prodiontot(ionspe=1,iz)--->N2+
c                       prodiontot(ionspe=2,iz)--->O2+
c                       prodiontot(ionspe=3,iz)--->O+
c                       prodiontot(ionspe=4,iz)--->N+
c                       prodiontot(ionspe=5,iz)--->H+
c                       prodiontot(ionspe=6,iz)--->He+
c 	chaufelec = chauffage electronique par perte continue (Heat dans
c 	       trans.f
c 	Ne_supra(z) = 1er moment : Ne suprathermique [cm-3]
c	courant_supra(z) = 2eme moment : courant suprathermique
c 			[cm-2.s-1]
c	Te_supra(z) = 3eme moment : Te suprathermique [K]
c	Chaleur_supra(z) = 4eme moment : flux d'energy suprathermique
c		           (Heat flux) [eV.cm-2.s-1]
c
c ====================================================================
c 	Initialisation des entrees (selon qui appelle)
c ====================================================================
c
c
c	if(kiappel.eq.1)write(6,*)'transelec, appele par transsolo' 'transelec called by transsolo'
c	if(kiappel.eq.2)write(6,*)'transelec, appele par transcar'  'transelec called by transcar'
c	write(6,*)
c
       if (kiappel.eq.1) then
c
c 	  On est appele par un transsolo.f. Les entrees sont lues dans
c 	  ELEC et NEUTRAL
c	  It is called by transsolo.f. The entries are read and ELEC NEUTRAL.

          write(6,*)
          write(6,*)'lect.f : reading the inflow of transport.f'
          write(6,*)'-------'
          print*,'call lect'
          call lect (nspec,knm,nen,nalt,zbot,ztop,hrloc,UT,day,year,
     .     jpreci,tempexo,f107(2),f107(3),Apind,chi,chideg,glat,glong,
     .          albedo,altkm,altcm,tneutre,densneut,colden,botE,centE,
     .          ddeng,nang,nango2,angzb,gmu,gwt,fluxdown,fluxup,
     .          denelc,temelc,temion,dipang,smgdpa)
         ap(1)=Apind
c
       elseif (kiappel.eq.2) then
c 	  On est appele par le transport fluide
         jpreci = jpreci0
            print*,'call iniflu'
            call iniflu(npt,iyd,UTsec,z,glat,glong,stl,f107,
     .                ap,chi,Ne,Te,Tj,indlim,jpreci,
     .                Nh,No,No2,Nn2,Nn,Tn,N_0,T_0,Po,Po2,Pn2,Ph,Pn,Heat,
     .
     .                  nspec,knm,nen,nang,nango2,nalt,
     .                  ddeng,botE,centE,gmu,gwt,angzb,altkm,altcm,
     .                  dipang,smgdpa,
     .                  fluxup,fluxdown,densneut,tneutre,tempexo,
     .                  albedo,chideg,denelc,temelc,temion,UT,hrloc,
     .                  ztop,zbot,colden,day,year)
c
        else
          write(6,*)'kiappel is false'
          stop
c
      endif            ! endif kiappel

c
c ====================================================================
c 	Calcul de la degradation en energie
c	Calculation of degradation in energy
c ====================================================================
      print*,'call degrad'
        call degrad(knm,nen,centE,botE,ddeng,nspec,kiappel)
c
c ====================================================================
c 	Calcul de la photoproduction primaire
c	Calculation of primary photoproduction
c ====================================================================
c
       if(jpreci.ne.1 .and. jpreci.ne.3 .and. jpreci.ne.4)
     .  call felin(knm,nspec,hrloc,day,year,UT,
     .  tempexo,f107,ap,glat,glong,nen,botE,centE,
     .  ddeng,nalt,altkm,tneutre,densneut,colden,chi,chideg,
     .          kiappel,phdisso2,pfluxsr,Po1sdisso2)
c
c ====================================================================
c 	Transport : calcul du flux stationnaire d'electrons
c	Transport: Calculation of the steady flow of electrons
c ====================================================================
c
c
        call trans(knm,nspec,nalt,zbot,ztop,hrloc,day,year,jpreci,
     .        tempexo,f107,ap,chideg,glat,glong,albedo,
     .        altkm,altcm,tneutre,densneut,colden,nang,nango2,
     .        angzb,gmu,gwt,nen,centE,botE,ddeng,fluxdown,fluxup,denelc,
     .        temelc,temion,smgdpa,prodiontot,chaufelec,kiappel,
     .        Ne_supra,courant_supra,Te_supra,Chaleur_supra,ut)
c
c ====================================================================
c 	Integration des productions et flux chaleur
c	Integration of production and heat flow
c ====================================================================
c
c 	Uniquement pour le transport fluide
c	Solely for fluid transportation
      if (kiappel.eq.2) then
            print*,'call cineout'
           call cineout(nalt,chaufelec,denelc,prodiontot,
     .          Ne_supra,courant_supra,Te_supra,Chaleur_supra,
     .      	Ne,npt,indlim,nx,zlim,zlim_1,z,Heat,Ph,Po,Po2,Pn2,Pn,
     .          Ne_sup,courant_sup,Te_sup,Chaleur_sup)

      end if
      end subroutine transelec
