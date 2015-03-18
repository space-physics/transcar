        subroutine trans(knm,nspec,nalt,zbot,ztop,hrloc,day,year,jpreci,&
     &        tempexo,f107,ap,chideg,glat,glong,albedo,                 &
     &        altkm,alt,tneutre,densneut,colden,nang,nango2,            &
     &        angzb,gmu,gwt,nen,centE,botE,ddeng,fluxdown,fluxup,denelc,&
     &        temelc,temion,smgdpa,prodiontot,chaufelec,kiappel,        &
     &        Ne_supra,courant_supra,Te_supra,Chaleur_supra,ut)

!
!-----------------------------------------------------------------------
!	This code solves the transport equation with the multistream
!	approach using the modified DISORT subroutine
!	(C) Copyright by D.Lummerzheim, 1985
!	                 J.Lilensten  , 1993
!-----------------------------------------------------------------------
!	Version with anisotropic boundary condition
!-----------------------------------------------------------------------
! 	jl, Avril 1995
!            Version de trans originale ou les sorties sont 
!	ecrites dans un format lisible, et on ne recompile pas quand on 
!	veut changer un nbre d'energies, angles ou altitudes. 
!	Dans le processus d'indentation et de lisibilite du pgmme,
! 	quelques bugs ont ete corriges, sur les coefficients d'emission
!	et des dimensions mal faites. A ces bugs pres, memes resultats
!	que l'original.
! 	Programme de dessins inclus.
! 	La photoproduction primaire est incluse.
! 	Ce programme peut etre appele indifferemment par le couplage
! 	fluide/cinetique ou independemment.
! 	Correction de la normalisation des flux (nov 2000)
!
!-----------------------------------------------------------------------
	implicit logical (l)
	include 'TRANSPORT.INC'
!
! In the following comments
! 	 z = altitude, E = energy, A = angle, sp = species
! 		 exc = excitation, ionst = ion state
!
! ====================================================================
! Computational parameters :
! ====================================================================
! nbren    =	IPN: max. length of the energy grid
! nen      =	IPN: actual length of the energy grid
! nbralt   =	IPM:(MMAX) length of altitude grid (number of layers +1)
! nalt	   =	actual # of altitudes
! nbrsp    =	IPJ: max # of species (1=N2 , 2=O2 , 3=O , 4=H , 5=He)
! nspec    =	IPJ: actual # of species 
! nbrexc   =	IPJS: maximum number of excitation states
! nbrionst =	IPJSP: number of ionized states
! nbrango2 =	IPK :  one half of the maximum number of streams
! nango2   =    kstr2 : one half of the actual number of streams
! nang     = 	actual number of stream (must be less or equal 16; 
! 	        should preferably be a power of two, if more than 16 
! 		streams are required,recompile all subroutines with 
! 		larger arrays)
! nbrlayer =	nbralt-1: number of layers
! ncalc  =  2*nbralt-1: number of layer boundaries (computational layer)
! jsg (sp) : number of states
! jsp (sp) : number of excited ion states
!
! e(E) 	   = energy grid [eV]
! engdd(E) = energy width grid [eV]
!
! alt(z)   = altitude grid [cm]
! altkm(z) = altitude grid [km]
!
! tau(0:z-1) = collision depth [1]
!
! pitchang(2*nbrango2) = angle grid [sr]
! weitang(2*nbrango2) = gaussian weigths     (1 ...    nbrango2 : up, 
!			             nbrango2+1 ...  2*nbrango2 : down)
! cosang(2*nbrango2) = cos(pitchang)
!
! ====================================================================
! Atmosphere and ionosphere parameters :
! ====================================================================
!
! Cross sections
! --------------
! ethres(sp,exc,ionst) = threshold energy [eV]
! bratio(ionst,sp) = branching ratio
! ntherm(z) = energy index for thermal electrons
! cel(sp,E) = elastic cross sections [cm2]
! cin(sp,E) = inelastic cross sections [cm2]
! cinex(sp,exc,E) = excitation cross sections [cm2]
! omdeg(E,sp) = differential cross section for primary (exc. and ion.)
! 		[cm2.eV-1]
! omsec(E,sp) = differential cross section for secondary [cm2.eV-1]
! elosse(E,z) = loss function (continuous slowing down approximation) 
! 		due to electron-electron and Coulomb interaction
! 		[eV.cm2]
!
! Neutral atmosphere
! ------------------
! densneut(sp,z) = Neutral density [cm-3] for N2,O2,O,H,He ...
! dentot(z) = Total neutral density [cm-3]
! colden(sp,z) = Column density [cm-2]
! denmass(z) = Mass of the neutral gaz
! press(z) = Atmospheric pressure [g.cm-2]
! yyddd = date (integer)
! tneutre(z) = Neutral temperature [K]
! smgdpa(z) = sin(magnetic dip angle)
! chideg,year,day,hrloc,tempexo,glat,Apind,f107,f107bar
!
! Ionosphere
! ----------
! denelc(z) = Electron density [cm-3]
! temelc(z) = Electron temperature [K]
! temion(z) = Ion temperature [K]
!
! ====================================================================
! Transport parameters :
! ====================================================================
!
!-------------------------------------------------------------------
! Q**	--> 	INTENSITIES 			[cm-2 s-1 eV-1 sr-1]  
!-------------------------------------------------------------------
!
! qxdown(E,-A/2 -->1)= Electron precipitations (ex qext)
! 	        -nango2 = vertical
! qxup(E,1 -->A/2)   = Measured e- up intensity (created)
! qdwn (-A/2 -->1)   = Electron precipitations (ex fiso)
! qntsty(E,z,A)      = Computed intensity (ex fint), 3 niveau/couche
! intensite 	     = Computed intensity (cree), 1 niveau/couche
! qgaupin(z,A)       = Computed intensity (ex gaupin)
! qint (z,-A-->A)    = Source from higher energies in the E loop 
! qprimpHot(E,z,A)   = Source function from primary photo electrons 
! qprimpRot(E,z,A)   = Source function from primary proto electrons 
! qprim(E,z,A)       = Total source function 

!
!-------------------------------------------------------------------
! F**	-->	FLUXES
!-------------------------------------------------------------------
!
! Projection en angle d'ordre 1 : flux hemispherique [eV-1.cm-2.s-1]
! ------------------------------------------------------------------
! fluxprim (z,E) = Primary electron flux from felin.f
! fhemu (E,z)    = Up  projected flux (2pi.sum(mu.qntsty.dmu)) (ex fup)
! fhemd (E,z)    = dwn projected flux (2pi.sum(mu.qntsty.dmu)) (ex fdwn)
! fhemtot(E,z)   = fhemd + fhemu
!
! Projection on angles and integration in energy : 
! ------------------------------------------------
! 			Hemispherical particle fluxes [cm-2.s-1]
!       		----------------------------------------
! fpartup(z)  = up particle flux (sum(dE.fhemu)) 
! fpartdwn(z) = down particle flux (sum(dE.fhemd))
! fpartsum(z) = Total particle flux (sum(dE.fhemtot))
!
! Projection on angles and energy : 
! ---------------------------------
!			Hemispherical energy fluxes [eV.cm-2.s-1]
!			-----------------------------------------
! feup(z) = up energy flux (sum(E.dE.fhemu))
! fedwn(z) = down energy flux (sum(E.dE.fhemd))
! fesum(z)= Total energy flux
!
!-------------------------------------------------------------------
! ENERGY
!-------------------------------------------------------------------
!
! engdep (z,2) = Energy deposition per layer
! edep         = Energy depositon per cm at layer boundaries
! Qeflux (z)   =  heating by fast electrons [eV cm-3 sec-1]
! Qetherm (z)  =  thermalized heating electron flux [eV cm-3 sec-1]
! chaufelec (z)       = Qeflux + Qetherm [eV cm-3 sec-1]
! enrate(sp+1,exc+1,z+1) = E deposited in ionisation and excitation.
!	 		 [eV.cm-3.s-1]
! 	  enrate(ist<jsg(isp)) --> exc. for specie isp, at each alt.
!	 		 [eV.cm-3.s-1]
! 	  enrate(jsg(isp)) --> ion. for specie isp, at each alt.
!	 		 [eV.cm-3.s-1]
! 	  enrate(jsg(isp)+1) --> total for specie isp, at each alt.
!	 		 [eV.cm-3.s-1]
!	  enrate(nalt+1)     --> height integrated for each specie.
!	 		 [eV.cm-2.s-1]
! enrion(z)    = ionization energy for all species.
!
! qphelerg  = Energy input due to photoelectrons [erg/cm2/sec]
! qpheleV   = Energy input due to photoelectrons [eV/cm2/sec]
! qprecerg  = Energy input due to precipitation [erg/cm2/sec]
! qpreceV   = Energy input due to precipitation [eV/cm2/sec]
! qtoterg   = Total energy input [erg/cm2/sec]
! qtoteV    = Total energy input [eV/cm2/sec]
! qsump     = Reflected energy [eV/cm2/sec]
! entot     = Total absorbed energy [eV/cm2/sec]
! shsum(sp)    = Total energy deposition [eV cm-2 sec-1]
!                         through inelastic coll. [eV/cm2/sec]
! shsumtot  = total E. absorbed through inelastic coll. [eV/cm2/sec]
! elhsum    = E. absorbed through heating  [eV/cm2/sec]
! elratio(z,2) = % of energy due to heating
! relerr    = Rate of energy conservation
! elos2     = Average energy loss per ion pair [eV] from qtot
! elos      = Average energy loss per ion pair [eV] from entot
! elos1     = Average energy loss per ion pair [eV] from {qtot-qsump}
!
!-------------------------------------------------------------------
! PRODUCTION
!-------------------------------------------------------------------
!
! photelec(z) = primary pHotoelectron production [cm-3.s-1]
!
! protelec(z) = primary pRotoelectron production [cm-3.s-1]
! prodionprot (isp,ialt) primary ion production rate due to protons  :
! 	 cm-3.s-1
!   1,2,3,6 -> N2+,O2+,O+,N+ (H+ et He+ ne sont pas calcules)

!
! Par rapport a felin, je prend l'ordre N2+,O2+,0+,H+,He+,N+ ,
! qui est celui dans lequel sont ecrites les especes dans 
! TRANSPORT.INC. Cela necessite de refaire passer N+ de la 4eme 
! position a la 6eme. Cela se fait dans phel.f
! prodionphot(sp+1,z) = Ion production due to photoionization [cm-3.s-1]
! 			sp + 1 tient compte de N+ dissocie par photons.
!               	[cm-3.s-1]
!               	prodionphot(ionspe=1,iz)--->N2+ 
!               	prodionphot(ionspe=2,iz)--->O2+ 
!               	prodionphot(ionspe=3,iz)--->O+ 
!               	prodionphot(ionspe=4,iz)--->H+ 
!               	prodionphot(ionspe=5,iz)--->He+ 
!               	prodionphot(ionspe=6,iz)--->N+ 
!
! primelec(z,E) = primary photoelectron production/eV  [cm-3.s-1.eV-1]
!
! primprotelec(z,E)=primary protoelectron production/eV  [cm-3.s-1.eV-1]
!
! prodionsec(sp,z) = Ion production due to transport [cm-3.s-1]
!                       prodionsec(ionspe=1,iz)--->N2+ 
!                       prodionsec(ionspe=2,iz)--->O2+ 
!                       prodionsec(ionspe=3,iz)--->O+ 
!                       prodionsec(ionspe=4,iz)--->H+ 
!                       prodionsec(ionspe=5,iz)--->He+ 
!                       prodionsec(ionspe=6,iz)--->N+ 
! prodiontot(sp,z) =  prodionprim(sp,z) + prodionsec(sp,z)
!                       prodiontot(ionspe=1,iz)--->N2+ 
!                       prodiontot(ionspe=2,iz)--->O2+ 
!                       prodiontot(ionspe=3,iz)--->O+ 
!                       prodiontot(ionspe=4,iz)--->H+ 
!                       prodiontot(ionspe=5,iz)--->He+ 
!                       prodiontot(ionspe=6,iz)--->N+ 
! prodelsec(z) = Electron production due to transport [cm-3.s-1]
! prate (sp,exc+1,z+1) = taux de production d'excitation ou ionisation
! 		         [cm-3.s-1]
! 	prate (ialt<=nalt)  = [cm-3.s-1]
! 	prate (ialt=nalt+1) = [cm-2.s-1] : integrale en colonne
!
! cemis (z) = integrated emission rate [cm-2.s-1]
! 
!-------------------------------------------------------------------
! SUPRATHERMAL MOMENTUM 
!-------------------------------------------------------------------
!
! First momentum : suprathermal electron density [cm-3] 
!       Ne_supra        = 1er  moment   [cm-3]
! Second momentum : suprathermal courant [cm-2.s-1] 
!       courant_supra   = 2eme moment   [cm-2.s-1]
! Third momentum : suprathermal temperature [K]
!       Te_supra        = 3eme moment   [K]
! Fourth momentum : Energy flux (Heat flux) [eV.cm-2.s-1]
!       Chaleur_supra   = 4eme moment   [eV.cm-2.s-1]
! 
!-------------------------------------------------------------------
! MISCELANEOUS
!-------------------------------------------------------------------
!
! conductivity calculation
! ------------------------
! denselcalc(z) = Electron density in the E region out of trans [cm-3]
! cped(z)       = Pedersen conductivity [MHO.m-1]
! cpedsum       = Integrated Pedersen conductivity [MHO]
! chal(z)       = Hall conductivity [MHO.m-1]
! chalsum       = Integrated Hall conductivity [MHO]
! ratHoP(z)     = Ratio Hall/Pedersen conductivity
! ratHoPsum     = Integrated Ratio Hall/Pedersen conductivity
! cpedCS        = Integrated Pedersen conductivity from Senior [MHO]
! chalCS        = Integrated Hall conductivity from Senior [MHO]
! gyreave       = Average electron gyrofrequency [s-1]
! gyriave       = Average ion gyrofrequency [s-1]
! collionSN(z)  = ion/n collision frequency [s-1], Schunk&Nagy
! collionRG(z)  = ion/n collision frequency [s-1], Rishbeth&Garriot
! colle(z)      = electron/neutral collision frequency [s-1]
! collOp(z)     = O+/Neutral collision frequency [s-1]
! collNOp(z)    = NO+/Neutral collision frequency [s-1]
! collO2p(z)    = O2+/Neutral collision frequency [s-1]
! collen2(z)    = electron/N2 collision frequency [s-1]
! colleo2(z)    = electron/O2 collision frequency [s-1]
! colleo1(z)    = electron/O1 collision frequency [s-1]
 
! range calculation
! ------------------
! range(z)      = distance in g/cm2 for range computation.
! zprodmax(sp)  = Altitude of electron prod. max [km]
! prodmax(sp)   = Electron production at this height [cm-3]
 
!
! optical emission rates
! ----------------------
! flbh (photon intensity and filter counts)
!
! libre parcours moyen 
! --------------------
! lpm(z,E)
!
!-------------------------------------------------------------------
! SWITCHES
!------------------------------------------------------------------- 
! lmaxw/lgauss/lpower : shape of input spectrum
! qinput	: total energy input in erg/cm2/sec
! epeak,epara : peak energy and parameter (power)
! ehalf   : halfwidth of Gaussian spectrum = EPEAK*EHALF
! eshift  : param. to define the amplitude of the added powerlaw
!		  I(eshift)=I(epeak)
! lpower1 : ADDITIONAL powerlaw to be added on maxwell or gauss
! lpower2 : ADDITIONAL powerlaw at high en. (above epower, with 
!	 	  apower exp.)
! switches for turning on/off different options (these should
!	generally left at the default setting):
! lruther : sceening selection (Wedde-Strand vs. Rutherford)
! ldelta  : Wiscombe's delta-M method applied
! iphase,gphase	: phase function selection
! lmono,kmono	: monodirectional input spectrum
! lporter : low energy phase function from Porter et al., 1987
! albedo   : particle albedo at bottom of atmosphere
! switches for testing of the correct preformance of the code
! nstop	: test-stop after NSTOP times through the energy loop
! external/internal source on file (if LEXT = true, set LMAXW
!          to false, otherwise the Maxwellian will be superposed
!          to the external source)(specify as Hollerith constant '...'H)
!
!-------------------------------------------------------------------
! INTERNAL WORKING ARRAYS 
!------------------------------------------------------------------- 
! ctot(2*z-1) = working array for energy 
! utau(2*z-1) = differential in tau
! gls(0-->A) = phase function moments
! eup = up energy flux (sum(E.dE.fhemu))
! edwn = down energy flux (sum(E.dE.fhemd))
! partup  = up particle flux (sum(dE.fhemu)) 
! partdwn = down particle flux (sum(dE.fhemd))
! fldn(z) = work array for sub. MSTREAM
! flup(z) = work array for sub. MSTREAM
! temp(2),dens(8) = arrays for MSIS
! twork(2*nbralt-1) = array for layers
! zwork(nbralt) = array for altitudes
!              --> si kiappel = 1, ce sont les donnees 
!                   contenue dans ELEC et NEUTRAL pour un run unique
!                   du programme.
!	if kiappel = 1, the data are contained in ELEC and NEUTRAL for a single run of the program.
!
!               --> si kiappel = 2, ce sont les donnees necessaires
!                   au transport fluide
!	if kiappel = 2, these are the necessary data transmission fluid
	
	
	
	logical flg_err
	common/err/flg_err	
	
	
	
	
	
	real e(nbren),engdd(nbren),Ebot(nbren)
 
	real ethres(nbrsp,nbrexc,nbrionst),bratio(nbrionst,nbrsp)
	integer ntherm(nbralt)	! energy index for thermal electrons
 	integer kiappel,icolin
 
	real  cel(nbrsp,nbren),cin(nbrsp,nbren),
     .		cinex(nbrsp,nbrexc,nbren),ctot(2*nbralt-1),
     .		elosse(nbren,nbralt)
        
	real omdeg(nbren,nbrsp),omsec(nbren,nbrsp)
 
 	real fluxdown(nbren,nbrango2),fluxup(nbren,nbrango2)
	real qdwn(-nbrango2:-1),qxdown(nbren,-nbrango2:-1),
     .	   qxup(nbren,nbrango2),qint(2*nbralt-1,-nbrango2:nbrango2)
 	
	real densneut(8,nbralt),temelc(nbralt),denelc(nbralt),
     .		dentot(nbralt),alt(nbralt),colden(8,nbralt),
     .		tau(0:nbralt-1),altkm(nbralt),temion(nbralt),
     .		denmass(nbralt),press(nbralt),utau(2*nbralt-1)

! 	densig = somme des seff + L(E)
 	real densig(nbralt,nbren)
 	real Ne(nbralt),lpm(nbralt,nbren)
!
	integer jsg(nbrsp*2+1),jsp(nbrsp*2+1)
*	JSG : number of states
*	JSP : number of excited ion states
 
	real gls(0:2*nbrango2),ssalb(nbralt-1)
*	IPLYR is the number of layers

	real f107(3),ap(7)
 	real photelec(nbralt),prodionprim(nbrsp*2,nbralt)
 	real prodionphot(nbrsp*2,nbralt)
 	real protelec(nbralt),prodionprot(6,nbralt)
 	real prodiontot(nbrsp*2,nbralt),prodeltot(nbralt)
       	real primelec(nbralt,nbren),fluxprim(nbralt,nbren)
       	real primprotelec(nbralt,nbren)
       	real fluxprimprot(nbralt,nbren)
	real qprim(nbren,nbralt,-nbrango2:nbrango2)
	real qprimpHot(nbren,nbralt,-nbrango2:nbrango2)
	real qprimpRot(nbren,nbralt,-nbrango2:nbrango2)
!
	real qntsty(nbren,2*nbralt-1,-nbrango2:nbrango2),
     .	  qgaupin(2*nbralt-1,-nbrango2:nbrango2),
     .    fhemu(nbren,nbralt),fhemd(nbren,nbralt),
     .	  fhemtot(nbren,nbralt)
	real intensite(nbren,nbralt,-nbrango2:nbrango2)
 	real feup(nbralt),fedwn(nbralt),fesum(nbralt),eup,edwn
 	real fpartup(nbralt),fpartdwn(nbralt),fpartsum(nbralt),
     .	     partup,partdwn
!
	real fldn(nbralt),flup(nbralt)    ! work array for sub. MSTREAM
 
	real enrate(nbrsp*2,nbrexc+1,nbralt+1),chaufelec(nbralt),
     .		shsum(nbrsp), prate(nbrsp*2,nbrexc+1,nbralt+1),
     .		prodionsec(nbrsp*2,nbralt),prodelsec(nbralt),
     .		enrion(nbralt),engdep(nbralt,2),elratio(nbralt,2),
     .		cemis(nbralt),Qeflux(nbralt),Qetherm(nbralt),
     .		edep(nbralt)
!
! 	Next arrays are for conductivity calculation.
      	real denselcalc(nbralt)
 	real cped(nbralt),chal(nbralt),ratHoP(nbralt)
 	real ratHoPsum,cpedsum,chalsum,cpedCS,chalCS
 	real gyreave,gyriave
 	real collionSN(nbralt),colle(nbralt),collionRG(nbralt)
 	real collOp(nbralt),collNOp(nbralt),collO2p(nbralt)
 	real collen2(nbralt),colleo2(nbralt),colleo1(nbralt)
 
	real weitang(2*nbrango2),cosang(2*nbrango2),pitchang(2*nbrango2)
 
*	for LBH and NI (photon intensity)
        real wavelbh(105),flbh(105),lbhint
        real waveni(9),fni(9),tni(9),niint,niextint
! 
	character word*4
*	real t0,t1
	character *80 crsin,option,rdtin,exfile
	integer file,ien1,isp1
	equivalence (file,exfile)
 	common /mess/ messag
	character headline*80,messag*130,sun*3
	character*9 title(nbrexc,nbrsp)
	character*30 istdate
	integer yyddd
	real temp(2),dens(8),tneutre(nbralt),smgdpa(nbralt)
 
	real twork(2*nbralt-1),zwork(nbralt)
 
	integer mcount(5),ncountE,ncountA
 	real prodmax(nbrsp),zprodmax(nbrsp),range(nbralt)
        real Ne_supra(nbralt),courant_supra(nbralt),Te_supra(nbralt),
     .          Chaleur_supra(nbralt)


!-------Commonize these variables for calculations of optical emissions
!-------in main TRANSCAR code
        common/pexc/prate,fhemd,fhemu,e
        common/excinds/jsg
!-------MZ
 
*	Default values
!
	data  lchem/.false./
*	LPRINT	: switch to turn printed output on/off
*	LPLOT	: switch to turn output for plot routines on/off
*	LPLOTM  : create cumulative plotfile (for PLTMUL)
*	NCOUNT  : print every NCOUNT energy step
*	MCOUNT  : print every MCOUNT altitude step
*	    (1) : intensity/density	(2) : energy deposition
*	    (3) : excitation rates	(4) : electron heating
*	LCHEM   : generate output for ionchemistry code
! 	Switches for prints :
*	IPRT allows to turn individual print statements on or off
*	(1): input spectrum		(2): energy deposition and error
*	(3): energy grid		(4): electron density
*	(5): electron heating (plot)	(6): lossfunction
*	(7): integrated sources		(8): scale heights
*	(9): table (for plotting)	(10): model atmosphere
*	(11): intensities		(12): integrated int. 
*					      by E-sections
*	(13): integrated flux		(14): ambient electron heating
*	(15): energy dep in neutals	(16): excitation rates
*	(17): MSTREAM output (short)	(18): MSTREAM output (long)
*	(19): table of emissions	(20): optical emissions
! 	idess(1)  --> atmosphere     			
! 	idess(2)  --> production 				
! 	idess(3)  --> flux prim,hemispherical tot,totd,totu = f( E)   
! 	idess(4)  --> flux prim,hemispherical tot,totd,totu = f(mu)     
! 	idess(5)  --> flux prim,hemispherical tot,totd,totu = f(z)    
! 	idess(6)  --> Hemispherical total flux = f(z) (zoom)
! 	idess(7)  --> E-region conductivities
! 	idess(8)  --> flux totd,totu a l'alt. la + hte, avec entrees. 
! 	idess(9)  --> chaufelec
! 	idess(10) --> header
! 	idess(11) --> Ne
! 	idess(12) --> Range computation
! 	idess(13) --> Suprathermal momentums
 	integer iprt(40),idess(20),izplt(4),ieplt(4)
 	real eplt(4)
 	real centE(nbren),botE(nbren),ddeng(nbren)
 	real angzb(2*nbrango2),gmu(2*nbrango2),gwt(2*nbrango2)
! 
	data lmaxw/.false./lgauss/.false./lpower/.false./
     .	      lpower1/.false./epeak/300./epara/2./
     .	      eshift/1./ehalf/.1/lpower2/.false./epower/1000./apower/3./
*	LMAXW/LGAUSS/LPOWER : shape of input spectrum
*	QINPUT	: total energy input in erg/cm2/sec
*	EPEAK,EPARA : peak energy and parameter (power)
*	ehalf   : halfwidth of Gaussian spectrum = EPEAK*EHALF
*	eshift  : param. to define the amplitude of the added powerlaw
*		  I(eshift)=I(epeak)
*	lpower1 : ADDITIONAL powerlaw to be added on maxwell or gauss
*	lpower2 : ADDITIONAL powerlaw at high en. (above epower, with 
*	 	  apower exp.)
!
 	logical onlyfl,usrang,exsorc,usrtau
 	data lruther/.true./iphase/0/gphase/.5/
!
*	switches for turning on/off different options (these should
*		generally left at the default setting):
*	LRUTHER : sceening selection (Wedde-Strand vs. Rutherford)
*	LDELTA  : Wiscombe's delta-M method applied
*	IPHASE,GPHASE	: phase function selection
*	LMONO,KMONO	: monodirectional input spectrum
*	lporter : low energy phase function from Porter et al., 1987
*	albedo   : particle albedo at bottom of atmosphere
! 
	data lt1/.true./lt3/.true./lt4/.true./consto/0.d0/
     .		constb/0.d0/lt5/.true./etest/100./nstop/0/
!	switches for testing of the correct preformance of the code
!	NSTOP	: test-stop after NSTOP times through the energy loop
! 
1000   format('No.=',i10,';',i3,' E;',i2,' A;',i3,' Alt; Sun=',
     .        a3,';yr=',i2,';day=',i3,';LT=',f5.2,';lat=',f6.2,
     .        ';Ap=',f5.1,';f10.7=',f6.2,';fav=',f6.2,
     .        ';texo=',f7.2,2x)
1001   format('No.=',i10,';',i3,' E;',i2,' A;',i3,' Alt; Sun=',
     .        a3,';yr=',i2,';day=',i3,';LT=',f5.2,';lat=',f6.2,
     .        ';Ap=',f5.1,';f10.7=',f6.2,';sza=',f6.2,
     .        ';texo=',f7.2,2x)
1010 	format('libre parcours moyen a ',1pe10.2,' eV et a ',1pe10.2,
     .         'eV',/,'altitude lpm 1   lpm2')
1020 	format(f10.2,2(1pe10.2))
1305 	format('   Energy    Qxdown   Computed    Qxup    Computed')
1300	format(5(1pe10.2))
!
2002	format(/1x,'The energy grid covers the range from :',f6.3,
     .		' eV up to :',f8.1,' eV')
2850    format (/,'atmosphere neutre, modele MSIS 90',/,
     .       ' no.   height    n(n2)     '
     .      ,' n(02)      n(01)',/,8x,' km   ',
     .       '   /cm3       /cm3       /cm3')
2855    format (/,' no.   height     n(He)      n(Ar)       ',
     .		'n(H)       n(N)', /,'         km        /cm3       ',
     .		'/cm3       /cm3       /cm3')
2860    format (i3,f10.2,4(1pe11.3))
2865    format (i3,2f10.2,3(1pe11.3))
2870    format (' no.   height  T neutral Mass-dsty   Pressure',/,
     .          '         km       [K]      [mbar]      [g]  ')
2875    format (' no.   height     Te         Ne ',/,
     .          '         km       [K]      [cm-3]')
2880    format (/,'COLUMN DENSITY',/, ' no.   height    n(n2)     '
     .      ,' n(02)      n(01)',/,8x,' km   ',
     .       '   /cm2       /cm2       /cm2')
2885    format (/,' no.   height     n(He)      n(Ar)       ',
     .		'n(H)       n(N)', /,'         km        /cm2       ',
     .		'/cm2       /cm2       /cm2')
!
4012	format(' *****************    Wiscombe''s delta-M method',
     .		' is NOT applied    ************')
4013	format(' Wiscombe''s delta-M method is applied')
4014	format(' *****************    Phase function is not',
     .		' Rutherford !!!!',
     .		t70,'IPHASE =',i3,',  GPHASE =',f8.4)
4016	format(' The screening parameter for the Rutherford phase',
     .		' function is calculated',/,
     .		' from the Bethe/Moliere formula for E > 500 eV')
!
5000	format('Date:',i6.5,' time:',f6.2,' LT'/'Longitude :',f6.1,
     .	  ' Latitude :',f6.1,/,'10.7 cm flux :',f6.1,
     .	  ' average :',f6.1,' Ap',f6.1,/,'MSIS exospheric',
     .	  ' temperature is ',f8.0,' K')
5007	format(
     .	'The total E. input from precipitation source is : ', 1pe10.4,
     .	' eV/cm2/s',/,45x,'or : ',1pe10.4,' erg/cm2/s',/,18x,
     .	'and from photoelectron source : ',1pe10.4,' eV/cm2/s',/,45x,
     .	'or : ',1pe10.4,' erg/cm2/s')
5008    format(
     .	'Input eV/cm2/s (erg/cm2/s): prec:',1pe9.2,' (',1pe8.2,')',
     .	' nue-: ',1pe9.2,' (',1pe8.2,')')
!
6002	format(' ',-5pf10.2,4(1pe12.4))
6006	format
     .	 ('Boundary conditions (intensity) (eV-1 cm-2 s-1 ster-1)  :',/,
     .		4(1pe12.5))
6007	format('Phase function moments :',/,f12.4,(4f12.4))
6008	format(/'    Altitude ',t15,' Flux (+) ',
     .	  t26,' Flux (-) ',t38,'Source (+)',t50,'Source (-)',/,
     .	  '    km',t17,'eV-1 cm-2 sec-1',t40,'eV-1 cm-2 sec-1'/)
6009	format
     .	 ('Calculated Intensity at boundary (eV-1 cm-2 s-1 ster-1) :',/,
     .		(4(1pe12.5)))
6010	format('Gaussian angle  : ',(4f12.5))
6011	format('Gaussian weight : ',(4f12.5))
6012 	format (/,'  Altitude',t21,'tau')
6013 	format (-5pf10.2,2(1pe20.8))
! 
7000	format(-5pf5.1,2x,7(1pe9.2))
7001 	format(/,' Alt.    ------ Energy Flux ------  E-depos. ',
     .   '---- Particle flux ----',/,' km            eV cm-2 sec-1',
     .   '       /eVcm3sec       cm-2 sec-1',/,'            (-)',
     .   '     (+)     total              (-)     (+)     total')
7002	format(/'Energy albedo at highest altitude          :',f11.5,/,
     .		'Particle albedo at highest altitude        :',f11.5/
     .          'Total energy deposition (qdep = Edwn-Eup)  :',1pe11.3,
     .		' eV cm-2 sec-1',//,
     .    'Characteristic energy of incident spectrum:',0pf12.3,
     .    ' keV ',/,12x,'or (ignoring flux below 50 eV):',f12.3,' keV'/
     .     7x,'of upward spectrum:       ',10x,f12.3,' keV ',/,12x,
     .    'or (ignoring flux below 50 eV):',f12.3,' keV')
7003	format(//' Energy deposited in heating of the ambient electrons'
     .  	,' [eV cm-3 sec-1]',/,'Alt. [km]',6x
     .		,'heating thermalization',3x,'total','  crossE [eV]')
!    .  	 '  % of e-dep'/)
7004	format(-5pf7.2,3x,3(1pe12.3),0pf7.2,2x,1pe8.1,2x,f7.3)
7005	format(//t10,' Energy deposition in ',a13,' (eV cm-3 sec-1)'//
     .		'Altitude   total ',5a10/)
7006	format(-5pf7.2,6(1pe10.2))
7008	format(/' Total Heating :',1pe11.3,' eV cm-2 sec-1'/)
7007	format(/' Total energy deposition in ',a5,' is ',1pe9.3,
     .	  ' eV cm-2 sec-1 (shsum)'/)
7009	format(//t10,' Excitation rates for ',a13,' (cm-3 sec-1)'//
     .		'Altitude  total   ',5a10/)
7010  	format(
     .  /'Total energy input of :',1pe9.3,' eV.cm-2.s-1 (qtoteV)',
     .  /,'             of which ',1pe9.3,' is reflected (qsump)',
     .	/,'                  and ',1pe9.3,' is absorbed (entot) ',
     .  /,25x,'through inelastic collisions (shsum =',1pe9.3,')',
     .	/,25x,'                  + heating (elhsum =',1pe9.3,')')
7011	format(/'total: ',10x,7(1pe10.2))
7012	format(//' The average E. loss per ion pair is :',
     .	         f7.2,' eV (from qtoteV)',/,34x,
     .	  'or :',f7.2,' eV (from (entot))',/,34x,
     .	  'or :',f7.2,' eV (from (qtoteV-qsump))')
7013	format(/t10,'Electron-electron interaction loss function :',
     .		' L/dE * [e]'//'   Energy   ',11(-5pf6.2,' km ')/)
7014	format(' ',12(1pe9.2,1x))
7015	format(1x,-5pf7.2,12x,2(1pe13.3))
7016	format(/' Height integrated :',2(1pe13.3))
7017	format(/' The average e. loss per ion pair of ',a13,
     .		  ' is',f7.2,' eV'/)
7018    format('Energy ',a4,' of',0pf7.2,' %')
7019	format(/'Warning: excessive downward energy flux thru lower border'/
     .    'Suggest: set albedo to 1. , increase Ne, reduce input flux'/)
7020	format(/'Integrated intensity between',f10.1,' and',f10.1,
     .    ' eV',/,' compared to integral from',f10.1,' eV to infinity'/)
7030	format(/,'Altitude ',7a10/)
7033	format(/'total: ',10(1pe10.2))
7050 	format(/,'Electron density',/,' altitude ',
     .	      ' e- prod.  Computed [Ne]   Input [Ne]')
7051 	format(1f10.2,3(1pe12.3))
7055 	format('Alt. of electron prod. max. from',a6,' : ',f10.2,' km',
     .	       /,'Value at this height  ',16x,' : ',1pe10.2,' cm-3')
7057 	format('Theoretical  range for [N2]  ',1pe10.2,' g/cm2')
!
7100 	format(/,10x,'Hemispherical fluxes [cm-2 s-1 eV-1]',/,10x,
     .	       37('-'),/,'fhemtot(E,z) = (2pi.sum(mu.qntsty.dmu))')
7110 	format(f10.2,' km       ',5(1pe10.2))
7120 	format(/,' Energies (eV) -->  ',5(1pe10.2),/,70('-'))
7130 	format(/,10x,'Upward hemispherical fluxes [cm-2 s-1 eV-1]',/,10x
     .	       ,37('-'),/,'fhemu(E,z) = (2pi.sum(mu.qntsty.dmu))')
7140 	format(/,10x,'Downward hemispherical fluxes [cm-2 s-1 eV-1]',/,
     .	       10x,37('-'),/,'fhemd(E,z) = (2pi.sum(mu.qntsty.dmu))')
!
8000	format(//15x,'Total Energy Deposition per layer',
     .		' (eV cm-2 sec-1)'/
     .		13x,'Altitude',3x,'E-flux diff.',
     .		9x,'E-dep.',5x,'rel. diff.'/)
8001	format(5x,-5pf15.2,2(2(1pe15.3),0pf15.3,2x))
!
!----------------------------------------------------------------------
!
        write(6,*)'    trans.f : Kinetic transport. jpreci=',
     .		jpreci,'ut=',UT
!


!       write(6,*)'Reed...                                      [A'
        call zeroit(qdwn,nbrango2)
	call zeroit(e,nbren)
	call zeroit(Qeflux,nbralt)
	call zeroit(Qetherm,nbralt)
	call zeroit(chaufelec,nbralt)
	call zeroit(engdd,nbren)
        call reed (iprt,idess,mcount,ncountE,ncountA,
     .	   linear,ldeltam,lporter,nspec,e,Ebot,engdd,nen,
     .	   nang,nango2,pitchang,cosang,weitang,angzb,gmu,gwt,
     .	   f107(3),f107(1),smgdpa,day,year,glong,alt,altkm,nalt,
     .	   tneutre,densneut,denelc,colden,title,jpreci,
     .	   zbot,ztop,hrloc,tempexo,knm,glat,ap(1),albedo,
     .	   qxdown,qxup,fluxup,fluxdown,temelc,temion,ezero,izplt,ieplt,
     .	   eplt,nke,jsg,jsp,ethres,bratio,cel,cin,cinex,
     .	   lamber,onlyfl,exsorc,usrang,usrtau,ddeng,centE,botE,icolin)


!
! 	Initialise phase function
!
	if(gphase.lt.0.) then
	  lgphase=.true.
	else
	  lgphase=.false.
	end if
!
	if(ldeltam) then
	  write(fic_transout,4013)
	else
	  write(fic_transout,4012)
	end if
	if(iphase.ne.0) write(fic_transout,4014) iphase,gphase
	if(lruther) write(fic_transout,4016)
!
! 	Initialise some values.
	nlayer=nalt-1
	mmax=2*nlayer+1
	pi=atan(1.)*4.
	boltz=8.6173e-05
	ergev=1.6022e-12	! convert from erg to eV
 	if (year.gt.99) then
 	  yyddd = int(year-1900.)*1000+int(day)
 	else
	  yyddd=int(year)*1000+int(day)
 	endif
!
	eng0=e(1)-engdd(1)/2.
	engmax=e(nen)+engdd(nen)/2.
	write(fic_transout,2002) eng0,engmax
! 
	do m=1,nalt
 	  denmass(m)=0.
 	  dentot(m)=0.
 	  do isp=1,nspec
	    denmass(m)=denmass(m)+densneut(isp,m)*atomas(isp)
   	    dentot(m)=dentot(m)+densneut(isp,m)
 	  enddo
	  denmass(m)=denmass(m)*1.6726e-24	! mass density in gramm
 	enddo
!
	do m=1,nalt
 	  press(m) = 0.
 	  do isp = 1,nspec
 	    press(m) = press(m)+colden(isp,m)*atomas(isp)
 	  enddo
	  press(m)=press(m)*1.6726e-24*0.981	! pressure in mbar
   	enddo
 
	write(fic_transout,5000) yyddd,hrloc,glong,glat,
     .			f107(1),f107(3),ap(1),tempexo
 	if(iprt(10).eq.1)then
          write(fic_transout,2850)
	  if(nspec.le.3)then
	    iwrite=nspec
	    do i=1,nalt,mcount(1)
              write(fic_transout,2860)i,altkm(i),
     .			(densneut(isp,i),isp=1,iwrite)
 	    enddo
	  else
	    iwrite=3
	    do i=nalt,1,-1
              write(fic_transout,2860)i,alt(i),
     .			(densneut(isp,i),isp=1,iwrite)
 	    enddo
            write(fic_transout,2855)
	    do i=nalt,1,-1
              write(fic_transout,2860)i,alt(i),
     .			(densneut(isp,i),isp=4,nspec)
 	    enddo
	  endif
 	  write(fic_transout,2870)
	  do i=1,nalt,mcount(1)
            write(fic_transout,2865)i,altkm(i),tneutre(i),
     .			denmass(i),press(i)
 	  enddo
 	  write(fic_transout,2875)
	  do i=1,nalt,mcount(1)
            write(fic_transout,2860)i,altkm(i),temelc(i),denelc(i)
 	  enddo
!
! 	  Densite colonne
          write(fic_transout,2880)
	  if(nspec.le.3)then
	    iwrite=nspec
	    do i=1,nalt,mcount(1)
              write(fic_transout,2860)i,altkm(i),
     .			(colden(isp,i),isp=1,iwrite)
 	    enddo
	  else
	    iwrite=3
	    do i=nalt,1,-1
              write(fic_transout,2860)i,alt(i),
     .			(colden(isp,i),isp=1,iwrite)
 	    enddo
            write(fic_transout,2885)
	    do i=nalt,1,-1
              write(fic_transout,2860)i,alt(i),
     .			(colden(isp,i),isp=4,nspec)
 	    enddo
	  endif
        endif
!

 	do ien = 1,nen
	  call eloss(e(ien),ien,temelc,denelc,elosse,
     .		nbren,nen,nbralt,nalt)
	  do ialt=1,nalt
!           densig(ialt,ien) = 0.
            densig(ialt,ien) = elosse(ien,ialt)/ddeng(ien)
 	    do isp = 1,nspec
              densig(ialt,ien) = densig(ialt,ien) + 
     .			(cin(isp,ien)+cel(isp,ien))*densneut(isp,ialt)
 	    enddo
 	  enddo
 	enddo
!       read in primary pHoto electron production rates. 
! 		photelec(z) 		: cm-3.s-1
! 	        prodionphot(isp+1,z)	: cm-3.s-1
! 		primelec(z,E) 		: cm-3.s-1.eV-1
! 		fluxprim(z,E)		: cm-2.s-1.eV-1
! 		qprimpHot(E,z,A)	: cm-2.s-1.eV-1.sr-1
        call phel(e,alt,altkm,iprt,fluxprim,
     .	      primelec,qprimpHot,photelec,nspec,nen,nang,nango2,nalt,
     .	      jpreci,mcount,prodionphot,densig)

!
!       read in primary pRoto electron production rates. 
! 		protelec(z) 		: cm-3.s-1
! 		primprotelec(z,E) 	: cm-3.s-1.eV-1
! 		fluxprimprot(z,E) 	: cm-2.s-1.eV-1
! 		qprimpRot(E,z,A)	: cm-2.s-1.eV-1.sr-1
!

        call prot(e,Ebot,ddeng,alt,altkm,iprt,
     .        nspec,nen,nang,nango2,nalt,jpreci,mcount,
     .        protelec,primprotelec,fluxprimprot,qprimpRot,
     .	      prodionprot,gmu,gwt,densig,densneut)


        do iang = -nango2,nango2
          do ialt = 1,nalt
            do ien = 1,nen
              qprim(ien,ialt,iang) = qprimpHot(ien,ialt,iang)+
     .                                 qprimpRot(ien,ialt,iang)
            enddo
          enddo
        enddo
c
c 	Compute input energy due to precipitations
 	if (jpreci.ne.0)then
  	  qpreceV = 0.
	  do iang=-nango2,-1,1
	    do ien=1,nen
	      qpreceV = qpreceV+qxdown(ien,iang)*weitang(iabs(iang))*
     .		      cosang(iabs(iang))*e(ien)*engdd(ien)
 	    enddo
 	  enddo
 	endif
	qpreceV = qpreceV*2.*pi
c
c 	Compute input energy due to photoionization
 	if (jpreci.ne.1 .and. jpreci.ne.3 .and. jpreci.ne.4)then
c 	  Integration over the energies
 	  do ialt = 1,nalt
 	    zwork(ialt) = 0.
	    do ien=1,nen
 	      zwork(ialt) = zwork(ialt)+
     .			    e(ien)*primelec(ialt,ien)*engdd(ien)
 	    enddo
 	  enddo
 	  call hint(nalt,alt,zwork,qpheleV)
 	endif
 	qtoteV   = qpheleV + qpreceV 
c
	qphelerg = qpheleV*ergev
	qprecerg = qpreceV*ergev
 	qtoterg  = qphelerg + qprecerg
c 
	write(fic_transout,5007) qpreceV,qprecerg,qpheleV,qphelerg
c
	write(6,5008) qpreceV,qprecerg,qpheleV,qphelerg
c	write(6,5009) yyddd,ut,qphelerg
5009    format(i10,f10.2,1pe12.5)

c
c--------------------------------------------------------------------
c
c	Solve the transport equation, starting with the highest energy
c       --------------------------------------------------------------
 
	if(iprt(11).eq.1.and.iprt(10).eq.1) write(fic_transout,*)
	read(irdtin) (bidon,j=1,nspec)	! read first line
	read(irdtin) (bidon,j=1,nspec)	! and forget it
 
*	Find the normalization factor for first energy step
 
	do n=1,nen
	  call eloss(e(n),n,temelc,denelc,elosse,nbren,nen,nbralt,nalt)
 	enddo
	do ialt=1,nalt
	  mm=2*ialt-1
	  ctot(mm)= elosse(nen,ialt)/engdd(nen)
	  do isp=1,nspec
	    ctot(mm)=ctot(mm)+
     .		(cin(isp,nen)+cel(isp,nen))*densneut(isp,ialt)
 	  enddo
 	enddo
	do ilyr=1,nlayer
	  mm=2*ilyr
	  ctot(mm)=(elosse(nen,ilyr)+elosse(nen,ilyr+1))*0.5/engdd(nen)
	  do isp=1,nspec
	    ctot(mm)=ctot(mm) + (cin(isp,nen)+cel(isp,nen))*
     .			(densneut(isp,ilyr)+densneut(isp,ilyr+1))*0.5
 	  enddo
 	enddo
	do m=1,mmax
	  if(ctot(m).le.1e-20) then
	    print*,'Message issued from /dir.cine/trans.f'
	    print*,nen,' Something is wrong with CTOT :',
     .	            (ctot(mm),mm=1,mmax)
	    stop 'error'		! call abort
	  end if
 	enddo
c
c	print*,'    Begin of energy loop'
 
	nloop=0
	call zeroit(qint,(2*nbralt-1)*(nbrango2*2+1))
	if(.not.lt4) nstop=1
	if(.not.lt3) then
	  ncountE=1
	  nstop=2
	end if
	if(.not.lt5) ntest=nlevtrans(e,engdd,etest,nen)
	nprt=mod(nen,ncountE)
	englimit=engmax
c
************************ energy loop ***********************
c
	do n=nen,2,-1
 	  lpr11=.false.
c	  Pour inhiber les impressions de disort
 	  lpr17=.false.
 	  lpr18=.false.
c
c 	  Si on veut les sorties de dissort mettre :
c	  lpr17=.true. 		! pour sorties courtes
c	  lpr18=.true. 		! pour sorties detaillees
c
 	  if(mod(n,ncountE).eq.nprt.and.iprt(11).eq.1)lpr11=.true.
	  if(n-1.eq.nen-nstop) then
	    print*,' Test stop after',nstop,' runs through the loop'
	    stop 'error'		! call abort
	  end if
	  eng = e(n)
 
*	  generate tau
 
	  call zeroit(twork,nlayer)
	  do m=1,nlayer
   	    twork(m+1)=twork(m)+(elosse(n,m)+elosse(n,m+1))/engdd(n)
     .		*(alt(m)-alt(m+1))/2.
 	  enddo
c
	  do ialt=1,nalt
	    tau(ialt-1)=twork(ialt)
	    do isp=1,nspec
   	      tau(ialt-1)=tau(ialt-1)+
     .		          (cin(isp,n)+cel(isp,n))*colden(isp,ialt)
 	    enddo
 	  enddo
	  if(.not.lt3) then
	    do m=1,nlayer
	      tau(m)=tau(m-1)+.5
 	    enddo
	  end if
	  do m=nlayer,0,-1
	    tau(m)=tau(m)-tau(0)
	    tau(m)=tau(m)/smgdpa(m+1)
 	  enddo
	  do m=1,nlayer
	    utau(2*m-1)=tau(m-1)
	    utau(2*m)=(tau(m-1)+tau(m))/2.
 	  enddo
	  utau(mmax)=tau(nlayer)
c 	  Pour eviter les erreurs d'arrondi, on retire a la derniere
c 	  difference une valeur petite (sinon, il arrive que la somme
c 	  des utau soit superieure a tau(nlayer), ce qui provoque un
c 	  arret dans disort.
	  utau(mmax)=utau(mmax) - tau(nlayer-1)/1000.
*	  generate phase function
 
	  if(eng.gt.12.) then
	    eps=(12./min(eng,englimit))**.75
	    eps=eps/(1.-eps)/2.
	    if(eng.gt.447..and.lruther) then
	      erel=eng/511000.
	      eps=6.22e-5/(2.+erel)/erel
	    end if
	  end if
	  if(lporter.and.eng.lt.447.) then     !alternative phase fction
	    call porter(gls,eng,nang)
	  else
	    call ruther(gls,eng,nang)
	  end if
	  if(lgphase) gphase=gls(1)
 
	  call zeroit(ssalb,nlayer)
	  do m=1,nlayer
	    do j=1,nspec
	      ssalb(m)=ssalb(m)+cel(j,n)*(densneut(j,m)+
     .			densneut(j,m+1))*0.5/ctot(2*m)
c     	if (n.eq.34) then
c     	print*,n,m,j,altkm(m),cel(j,n),densneut(j,m),densneut(j,m+1),
c     &ctot(2*m)
c     	endif
 	    enddo
   	    ssalb(m)=min(ssalb(m),0.99999)
 	  enddo
 
*	  solve equations

	if (flg_err) then
	  do m=1,nlayer
	    print*,n,altkm(m),tau(m),ssalb(m)
	  enddo
	endif
	  do k=-nango2,-1,1
	    qdwn(k)=qxdown(n,k)
 	  enddo
	  call mstream(eng,n,nang,nlayer,qdwn,tau,utau,ssalb,albedo,
     .		gls,eps,ldeltam,lpr17,lpr18,linear,iphase,gphase,
     .		qint,qgaupin,fldn,flup,nalt,lamber,onlyfl,exsorc,usrang,
     .		usrtau,qprim)
     	
	  flux=0.
	  do m=1,nalt
	    fhemd(n,m)=fldn(m)
	    fhemu(n,m)=flup(m)
 	    fhemtot(n,m) = fhemd(n,m)+fhemu(n,m)
   	  enddo
	  do k=-nango2,nango2,1
	    do m=1,mmax
   	      qntsty(n,m,k)=qgaupin(m,k)
 	    enddo
 	  enddo
 
*	  integrate intensity over dmu, dtau, and EdE
	  do m=1,nlayer
	    flux=flux + (qgaupin(2*m-1,0)+qgaupin(2*m+1,0))*
     .		(tau(m)-tau(m-1))/2.
 	  enddo
 
*	  Optional printing
 
	  if(lpr11) then
	    write(fic_transout,6007) (gls(k),k=0,nang)
	    write(fic_transout,6010) (cosang(k),k=nango2,1,-1)
	    write(fic_transout,6011) (weitang(k),k=nango2,1,-1)
	    write(fic_transout,6006) (qdwn(k),k=-nango2,-1,1)
	    write(fic_transout,6009) (qntsty(n,1,k),k=-nango2,-1,1)
c
	    write(fic_transout,6012)
	    do m=1,nalt,mcount(1)
   	      write(fic_transout,6013) alt(m),tau(m-1)
 	    enddo
	    write(fic_transout,6008)
	    do m=1,nalt,mcount(1)
	      qqdown=0.
	      qqup=0.
	      do k=-nango2,-1,1
   	        qqdown=qqdown+qint(2*m-1,k)*weitang(iabs(k))*
     .			cosang(iabs(k))
 	      enddo
	      do k=1,nango2,1
   	        qqup=qqup+qint(2*m-1,k)*weitang(iabs(k))*cosang(iabs(k))
 	      enddo
   	      write(fic_transout,6002) alt(m),fhemu(n,m),fhemd(n,m),
     .			qqup,qqdown
 	    enddo
	  end if
 
*	  Find the normalization factor for next energy step
 
	  if(n.gt.1) then
	    do m=1,nalt
	      ctot(2*m-1)= elosse(n-1,m)/engdd(n-1)
	      do j=1,nspec
	        ctot(2*m-1)=ctot(2*m-1)+(cin(j,n-1)+cel(j,n-1))*
     .			    densneut(j,m)
 	      enddo
 	    enddo
	    do m=1,nlayer
	      ctot(2*m)=(elosse(n-1,m)+elosse(n-1,m+1))*0.5/engdd(n-1)
	      do j=1,nspec
	        ctot(2*m)=ctot(2*m) + (cin(j,n-1)+cel(j,n-1))*
     .		(densneut(j,m)+densneut(j,m+1))*0.5
 	      enddo
 	    enddo
	    do m=1,mmax
	      if(ctot(m).le.1e-20) then
	        print*,n-1,' Something is wrong with CTOT :',
     .	                    (ctot(mm),mm=1,nalt)
	        stop 'error'		! call abort
	      end if
 	    enddo
 
c	    print*,n,'=n, calling qmstr'
	    if(lt1) call qmstr(n,mmax,nspec,engdd,densneut,ctot,elosse,
     .	        qint,qntsty,omdeg,omsec,weitang,nbren,nen,nbralt,
     .		nalt,nbrsp,nbrango2,nango2,irdtin)
	  end if
 
	  if(.not.lt3) then
	    do k=-nango2,nango2,1
	      do m=1,mmax
	        qint(m,k)=(1.-consto)/cosang(1)/2./pi
	        if(k.eq.0) qint(m,k)=0.
 	      enddo
 	    enddo
	  end if
	  do k=-nango2,nango2,1
   	    qint(mmax,k)=0.
 	  enddo
 
	  nloop=nloop+1
c         if (kiappel.eq.1)
c    .	  write(6,*)'Through with loop ',nloop+1,'/',nen,'         [A'
c         write(6,*)'Through with loop ',nloop+1,'/',nen,'         [A'
c
 	enddo
c   	print*,'    Done with the energy loop                    '
c
********************** end of energy loop *****************************
c 
c 	Pour plus de simplicite, recentre l'intensite.
 	do ien = 1,nen
 	  do ialt = 1,nalt
 	    ilyr = 2.*ialt-1
 	    do iang =-nango2,nango2
 	      intensite(ien,ialt,iang) = qntsty(ien,ilyr,iang)
 	    enddo
 	  enddo
 	enddo
c
        if(iprt(17).eq.1) then
c         Intensites calculees a toutes les altitudes, energies, et tous angles 
c 	  d'attaque
          do ialt = 1,nalt
            write(fic_transout,7150)altkm(ialt)
            do ien = 2,nen,5*ncountE
              write(fic_transout,7120) (e(ien+i*ncountE),i=0,4)
c             Downward ...
              do iang = 1,nango2
                write(fic_transout,7160)angzb(iang),
     .             (max(0.,intensite(ien+i*ncountE,ialt,nango2+1-iang)),
     .             i=0,4)
              enddo
              do iang = 1,nango2,1
                write(fic_transout,7170)angzb(iang+nango2),
     .             (max(0.,intensite(ien+i*ncountE,ialt,-iang)),
     .             i=0,4)
              enddo
            enddo
          enddo
        endif
c
 	if(iprt(8).eq.1) then
c 	  Intensites calculees aux altitudes programmees dans DATTRANS
 	  do ialt = 1,4
 	    iz = izplt(ialt)
 	    write(fic_transout,7150)altkm(iz)
 	    do ien = 2,nen,5*ncountE
 	      write(fic_transout,7120) (e(ien+i*ncountE),i=0,4)
c 	      Downward ...
 	      do iang = 1,nango2
                write(fic_transout,7160)angzb(iang),
     .             (max(0.,intensite(ien+i*ncountE,iz,nango2+1-iang)),
     .		   i=0,4)
 	      enddo
 	      do iang = 1,nango2,1
                write(fic_transout,7170)angzb(iang+nango2),
     .             (max(0.,intensite(ien+i*ncountE,iz,-iang)),
     .		   i=0,4)
 	      enddo
 	    enddo
  	  enddo
 	endif


7150 	format(/,10x,'Electron stationary flux [cm-2 s-1 eV-1 sr-1]',/,
     .		10x,47('-'),/,'qntsty(E,z,mu) at z =',f10.2,' km',
     .		/,'(0 deg --> fully down, 180 deg --> fully up)')
7160 	format('Down angle',f8.2,2x,5(1pe10.2))
7170 	format('Up   angle',f8.2,2x,5(1pe10.2))

	if(iprt(7).eq.1) then
c	  Flux hemispherique total
 	  write(fic_transout,7100)
 	  do ien = 2,nen,5*ncountE
 	    write(fic_transout,7120) (e(ien+i*ncountE),i=0,4)
 	    do ialt = 1,nalt,mcount(1)
 	      write(fic_transout,7110)altkm(ialt),
     .		   (max(0.,fhemtot(ien+i*ncountE,ialt)),i=0,4)
 	    enddo
  	  enddo
c	  Flux hemispherique montant
 	  write(fic_transout,7130)
 	  do ien = 2,nen,5*ncountE
 	    write(fic_transout,7120) (e(ien+i*ncountE),i=0,4)
 	    do ialt = 1,nalt,mcount(1)
 	      write(fic_transout,7110)altkm(ialt),
     .		   (max(0.,fhemu(ien+i*ncountE,ialt)),i=0,4)
 	    enddo
  	  enddo
c	  Flux hemispherique descendant
 	  write(fic_transout,7140)
 	  do ien = 2,nen,5*ncountE
 	    write(fic_transout,7120) (e(ien+i*ncountE),i=0,4)
 	    do ialt = 1,nalt,mcount(1)
 	      write(fic_transout,7110)altkm(ialt),
     .		   (max(0.,fhemd(ien+i*ncountE,ialt)),i=0,4)
 	    enddo
  	  enddo
 	endif
	if(iprt(13).eq.1) write(fic_transout,7001)
c
	call zeroit(engdep,nbralt*2)
	call zeroit(edep,nbralt)
	d1fe=0.
	do ialt=1,nalt
	  etherm=temelc(ialt)*boltz	! start at thermal energy
	  ntherm(ialt)=
     .		ncross(e,engdd,qntsty,ialt,etherm,denelc(ialt),nen)
c	  cross over of thermal/streaming electrons
	  estart=e(ntherm(ialt))	
	  call intgrl(estart,engmax,fhemd,e,engdd,partdwn,edwn,nen,ialt)
	  call intgrl(estart,engmax,fhemu,e,engdd,partup,eup,nen,ialt)
 	  feup(ialt) = eup
 	  fedwn(ialt) = edwn
 	  fpartup(ialt) = partup
 	  fpartdwn(ialt) = partdwn
	  fesum(ialt)=feup(ialt) + fedwn(ialt)
	  fpartsum(ialt)=fpartdwn(ialt) + fpartup(ialt)
	  d0fe=d1fe
	  d1fe=edwn - eup
	  ddz=1.
	  if(ialt.ge.2) then
	    engdep(ialt,1)=d0fe - d1fe
	    ddz=alt(ialt-1)-alt(ialt)
	  end if
 	  if(iprt(13).eq.1)then
	    if(mod(ialt,mcount(1)).eq.1.or.mcount(1).eq.1)
     .	       write(fic_transout,7000) alt(ialt),edwn,eup,fesum(ialt),
     .			engdep(ialt,1)/ddz,partdwn,partup,fpartsum(ialt)
 	  endif
	  if(ialt.eq.1) then
	    qdep=fedwn(ialt)-feup(ialt)
	    qsump=feup(ialt)
	    ealb=feup(ialt)/fedwn(ialt)
	    falb=fpartup(ialt)/fpartdwn(ialt)
	    fpartnet=fpartdwn(ialt)-fpartup(ialt)
	    echar=fedwn(ialt)/fpartdwn(ialt)/2000.
	    echarp=feup(ialt)/fpartup(ialt)/2000.
	  end if
 	enddo
	if(edwn/qtotev.ge.1.e-2) then
 	  write(fic_transout,7019)
 	  write(6,7019)
 	endif
	e1=50.
	call intgrl(e1,engmax,fhemu,e,engdd,sumf,sume,nen,1)
	echarp2=sume/sumf/2000.
	call intgrl(e1,engmax,fhemd,e,engdd,sumf,sume,nen,1)
c	Bob Robinson's definition of E_char
	if(iprt(12).eq.1) then
	  s1=0.
	  s2=0.
	  dz=0.
	  e1=10.
	  e2=20.
	  write(fic_transout,7020) e1,e2,e2
	  do m=1,nalt
	    if(m.ge.2) dz=(alt(m-1)-alt(m))/2./smgdpa(m)
	    su1=sum1
	    call intgrl(e1,e2,fhemd,e,engdd,sum1,sum,nen,m)
	    call intgrl(e1,e2,fhemu,e,engdd,sum2,sum,nen,m)
	    sum1=sum1+sum2
	    s1=s1+(su1+sum1)*dz
	    su2=sum2
	    call intgrl(e2,engmax,fhemd,e,engdd,sum2,sum,nen,m)
	    call intgrl(e2,engmax,fhemu,e,engdd,sum3,sum,nen,m)
	    sum2=sum2+sum3
	    s2=s2+(su2+sum2)*dz
	    if(mod(m,mcount(1)).eq.1)write(fic_transout,7015) alt(m),
     .			sum1,sum2
 	  enddo
	  write(fic_transout,7016) s1,s2
	end if
c 
c	Calculate production and energy deposition profiles
c
	call endep(nspec,nalt,e,cinex,engdd,ethres,bratio,qntsty,
     .		alt,smgdpa,densneut,jsg,jsp,shsum,enrate,prate,ntherm,
     .		nen)


c
c 	Calcul des moments suprathermiques
c
        call moments(nalt,alt,altkm,nen,e,engdd,
     .          nang,nango2,weitang,cosang,intensite,iprt,ntherm,
     .          Ne_supra,courant_supra,Te_supra,Chaleur_supra)

c
	call zeroit(prodelsec,nalt)
	call zeroit(enrion,nalt)
	do isp=1,nspec
*	    get deposition profile for species J
	    if (kiappel.eq.1)
     .      print*,'     ',specie(isp),' deposition profile done','[A'
	    do m=1,nalt
 	      prodionsec(isp,m) = 0.
 	    enddo
	    do m=1,nalt
	      enrion(m)=enrion(m)+enrate(isp,jsg(isp),m)
 	      prodionsec(isp,m)=prodionsec(isp,m)+prate(isp,jsg(isp),m)
   	      prodelsec(m)=prodelsec(m)+prate(isp,jsg(isp),m)
 	      if (isp.eq.1)then
c 		On tient compte de N+ dissocie
c 		Il ne faut pas incrementer enrion car l'energie
c 	        des ions N+ a deja ete comptee dans les etats excites.
c	 	Parmi les sections efficaces prises en compte, seule
c 	 	B'1Su est predissociative. Il s'agit de l'etat 12
 	        prodionsec(6,m)=prodionsec(6,m)+prate(1,12,m)
   	        prodelsec(m)=prodelsec(m)+prate(1,12,m)
 	     endif
 	    enddo
	    edep(1)=enrate(isp,jsg(isp)+1,1)
	    do m=2,nalt
	      edep(m)=edep(m)+enrate(isp,jsg(isp)+1,m)
	      engdep(m,2)=engdep(m,2)+
     .		(enrate(isp,jsg(isp)+1,m-1)+enrate(isp,jsg(isp)+1,m))*
     .		(alt(m-1)-alt(m))/2./smgdpa(m)
 	    enddo
	    if(iprt(15).eq.1.and.shsum(isp).ne.0.)  then
	      write(fic_transout,7005)specie(isp),(title(js,isp),js=1,5)
	      do m=1,nalt,mcount(2)
	        write(fic_transout,7006)alt(m),enrate(isp,jsg(isp)+1,m),
     .	   		(enrate(isp,js,m),js=1,5)
   	      enddo
	      write(fic_transout,7011) (enrate(isp,js,nalt+1),js=1,5)
 	      if (jsg(isp).gt.5)then
	        write(fic_transout,7030) (title(js,isp),js=6,jsg(isp))
	        do m=1,nalt,mcount(2)
	          write(fic_transout,7006) alt(m),
     .				(enrate(isp,js,m),js=6,jsg(isp))
   	        enddo
	        write(fic_transout,7033) (enrate(isp,js,nalt+1),
     .			js=6,jsg(isp))
 	      endif
	    end if
	    write(fic_transout,7007) specie(isp),shsum(isp)
 	enddo
	do m=1,nalt
   	if(prodelsec(m).le.0.) prodelsec(m)=1.e-36
 	enddo
 
	if(iprt(16).eq.1) then		! optional printing
	  do isp=1,nspec
	      write(fic_transout,7009)specie(isp),(title(js,isp),js=1,5)
	      do m=1,nalt,mcount(3)
   	        write(fic_transout,7006) alt(m),prate(isp,jsg(j)+1,m),
     .	  	              (prate(isp,js,m),js=1,5)
 	      enddo
	      write(fic_transout,7011) (prate(isp,js,nalt+1),js=1,5)
 	      if (jsg(isp).gt.5)then
	        write(fic_transout,7030) (title(js,isp),js=6,jsg(isp))
	        do m=1,nalt,mcount(3)
   	          write(fic_transout,7006)alt(m),
     .			(prate(isp,js,m),js=6,jsg(isp))
 	        enddo
	        write(fic_transout,7033) 
     .			(prate(isp,js,nalt+1),js=6,jsg(isp))
 	      endif
	      elos=shsum(isp)/prate(isp,jsg(isp),nalt+1)
	      write(fic_transout,7017) specie(isp),elos
 	  enddo
	end if
c 
	if(iprt(6).eq.1) then
	  write(fic_transout,7013) (alt(m),m=1,nalt,nalt/5)
	  do n=1,nlevtrans(e,engdd,160.,nen),1
   	    write(fic_transout,7014) e(n),
     .		(elosse(n,m)/engdd(n),m=1,nalt,nalt/5)
 	  enddo
	end if
	elhsum=0.
	do m=1,nalt
	  nn=ntherm(m)
c	  if(nn.eq.0)nn=1
 	  if(nn.le.1)nn=2
	  etherm=temelc(m)*boltz
	  Qetherm(m)=elosse(nn,m)*(e(nn)-1.5*etherm)*
     .	      qntsty(nn,2*m-1,0)
	  do n=ntherm(m),nen,1
	    Qeflux(m)=Qeflux(m)+elosse(n,m)*qntsty(n,2*m-1,0)*engdd(n)
   	    chaufelec(m)=Qeflux(m)+Qetherm(m)
 	  enddo
	  if(m.ge.2) then
	    ddh=(chaufelec(m-1)+chaufelec(m))/2.*(alt(m-1)-alt(m))/
     .			smgdpa(m)
	    engdep(m,2)=engdep(m,2)+ddh
	    edep(m)=edep(m)+chaufelec(m)
	    if(engdep(m,1).gt.0.) elratio(m,1)=ddh/engdep(m,1)
	    if(engdep(m,2).gt.0.) elratio(m,2)=ddh/engdep(m,2)
	  end if
 	enddo
	do m=2,nalt
	  ddz=(alt(m-1)-alt(m))/smgdpa(m)
   	  elhsum=elhsum+(chaufelec(m-1)+chaufelec(m))*ddz/2.
 	enddo
	
	  if(iprt(14).eq.1) then
	    write(fic_transout,7003)
	    do m=1,nalt,mcount(4)
	      write(fic_transout,7004) alt(m),Qeflux(m),Qetherm(m),
     .			chaufelec(m),e(ntherm(m))
c    . 		,elratio(m,1)*100.,elratio(m,2)*100.
   	    enddo
	  end if
	  write(fic_transout,7008) elhsum
 
	if(iprt(2).eq.1) then
	  write(fic_transout,8000)
	  do m=1,nalt,mcount(2)
	    if(engdep(m,1).gt.0.)
     .		edeperr=(engdep(m,1)-engdep(m,2))/engdep(m,1)
	    write(fic_transout,8001) 
     .			alt(m),engdep(m,1),engdep(m,2),edeperr
*    .	       prodelsec(m),press(m),engdep(m,2)/1e6/fpartnet/denmass(m)
 	  enddo
	end if
c
C 	if (iprt(20).eq.1)then
c 	  Emission lines intensities.
c
c 	  Calcul en fonction des longueurs d'onde (source ?)
C          call lbh(nalt,prate)
C          call ni(nalt,prate,colden,alt,smgdpa)
c
c 	  Calculs integres (jl, 1994/1995)
c 	  calcule les emissions de l'azote suivantes :
c         em3914,em4278,emMeinel,em2p,em1p,emVK,emLBH,emBH,em5001,
c         emNI,emNII,emRy
C          call emisazote (nalt,alt,nen,e,engdd,bratio,qntsty,cinex,
C     .          jsg,enrate,prate,prodionprim,ntherm)
c
c 	  1304 et 1356 A dues a O et O2
C	  if(nspec.ge.2)call oi(prate,colden,alt,smgdpa,nalt)
c
c	  D'autres emissions a verifier un jour...
C          call emission(cemis,prate,alt,smgdpa,densneut,nalt,jsg,
C     .                  emis4278,emis7320,prodoii2p,lossoii2p,
C     .                  densoii2p)
C
C	endif
c
*	calculate the energy balance
 
	entot=elhsum
	shsumtot=0.
	prtot=0.
 	if (jpreci.ne.1 .and. jpreci.ne.3 .and. jpreci.ne.4)
     .		call hint(nalt,alt,photelec,prtot)
	do j=1,nspec
	  prtot=prtot + prate(j,jsg(j),nalt+1)     ! total ionization
   	  entot=entot+shsum(j)		! total energy deposition
   	  shsumtot=shsumtot+shsum(j)	! Energy deposition/coll.
 	enddo
c 	average energyloss per ion pair
c 	Il faut tenir compte de la stationnarite : l'energie des photons
c 	au fur et a mesure qu'elle cree des electrons primaires puis
c 	secondaires est renouvelee : il faut donc la recompter dans 
c 	l'energie par paire.
	elos=(qpheleV+entot)/prtot     
	elos1=(qpheleV+qtoteV-qsump)/prtot
	elos2=(qpheleV+qtoteV)/prtot
 	
	if(qtoteV.ne.0.) relerr=((qtoteV-qsump-entot)/qtoteV)*100.
	word='loss'
	if(relerr.lt.0.) word='gain'
c
c	Computes the total ion and electron production
 	do ialt = 1,nalt
c 	  electrons
 	  prodeltot(ialt) = photelec(ialt) + prodelsec(ialt)+
     .			    protelec(ialt)
c 	  N2+
 	  prodionprim(1,ialt) = prodionphot(1,ialt)+ prodionprot(1,ialt)
 	  prodiontot(1,ialt) = prodionprim(1,ialt) + prodionsec(1,ialt)
c 	  O2+
 	  prodionprim(2,ialt) = prodionphot(2,ialt)+ prodionprot(2,ialt)
 	  prodiontot(2,ialt) = prodionprim(2,ialt) + prodionsec(2,ialt)
c 	  O+
 	  prodionprim(3,ialt) = prodionphot(3,ialt)+ prodionprot(3,ialt)
 	  prodiontot(3,ialt) = prodionprim(3,ialt) + prodionsec(3,ialt)
c 	  H+
 	  prodionprim(4,ialt) = prodionphot(4,ialt)+ prodionprot(4,ialt)
 	  prodiontot(4,ialt) = prodionprim(4,ialt) + prodionsec(4,ialt)
c 	  He+
 	  prodionprim(5,ialt) = prodionphot(5,ialt)+ prodionprot(5,ialt)
 	  prodiontot(5,ialt) = prodionprim(5,ialt) + prodionsec(5,ialt)
c 	  N+
 	  prodionprim(6,ialt) = prodionphot(6,ialt)+ prodionprot(6,ialt)
 	  prodiontot(6,ialt) = prodionprim(6,ialt) + prodionsec(6,ialt)
 	enddo
2010 	format(a6,' ion production : ',/,
     .		'Altitude   Primary    Secondary  Ratio S/P, Total',/,
     .          ' [km]     [cm-3s-1]   [cm-3s-1]          [cm-3s-1]')
2015 	format(' Electron production : ',/,
     .		'Altitude   Primary    Secondary  Ratio S/P, Total',/,
     .          ' [km]     [cm-3s-1]   [cm-3s-1]          [cm-3s-1]')
2020 	format(1f10.2,4(1pe10.2))
2030 	format(1f10.2,2(1pe10.2),10x,1pe10.2)
 	if(iprt(19).eq.1) then
 	  do isp = 1,nspec
 	    write(fic_transout,2010)specie(isp)
 	    do ialt = 1,nalt
 	      if(prodionprim(isp,ialt).ne.0.)then
 	        rapport = prodionsec(isp,ialt)/prodionprim(isp,ialt)
 	        write(fic_transout,2020)
     .		  altkm(ialt),prodionprim(isp,ialt),
     .		  prodionsec(isp,ialt),rapport,prodiontot(isp,ialt)
 	       else
 	        write(fic_transout,2030)
     .		  altkm(ialt),prodionprim(isp,ialt),
     .		  prodionsec(isp,ialt),prodiontot(isp,ialt)
 	       endif
 	    enddo
 	  enddo
c	  On ecrit le N+ qui vient de la dissociation de N2
          isp = 6
            write(fic_transout,2010)specie(isp)
            do ialt = 1,nalt
              if(prodionprim(isp,ialt).ne.0.)then
                rapport = prodionsec(isp,ialt)/prodionprim(isp,ialt)
                write(fic_transout,2020)
     .            altkm(ialt),prodionprim(isp,ialt),
     .            prodionsec(isp,ialt),rapport,prodiontot(isp,ialt)
               else
                write(fic_transout,2030)
     .            altkm(ialt),prodionprim(isp,ialt),
     .            prodionsec(isp,ialt),prodiontot(isp,ialt)
               endif
            enddo
 	  continue
c
 	  write(fic_transout,2015)
 	  if(jpreci.eq.0)then
 	    do ialt = 1,nalt
 	      rapport=0.
  	      if(photelec(ialt).ne.0.)
     .			rapport = prodelsec(ialt)/photelec(ialt)
 	      write(fic_transout,2020)altkm(ialt),photelec(ialt),
     .		prodelsec(ialt),rapport,prodeltot(ialt)
 	    enddo
 	  elseif (jpreci.eq.3) then
 	    do ialt = 1,nalt
 	      rapport=0.
  	      if(protelec(ialt).ne.0.)
     .			rapport = prodelsec(ialt)/protelec(ialt)
 	      write(fic_transout,2020)altkm(ialt),protelec(ialt),
     .		prodelsec(ialt),rapport,prodeltot(ialt)
 	    enddo
 	  endif
 	endif
c
c 	Computes the electrons and ion densities in the E region
 	if(iprt(4).eq.1 .or. idess(10).ne.-1)then
          call densout(hrloc,nalt,altkm,prodeltot,denelc,denselcalc)
c	  if (iprt(4).eq.1)then
c
c47 	format(1i4,4f10.2,2(1pe12.3))
 	    write(fic_transout,7050)
 	    do ialt = 1,nalt
 	      write(fic_transout,7051)altkm(ialt),prodeltot(ialt),
     .		  denselcalc(ialt),denelc(ialt)
 	    enddo
 	endif
!
! 	Computes the E-region conductivities
 	if(iprt(23).eq.1 .or. iprt(24).eq.1 .or. idess(11).eq.1) then
 	  do ialt = 1,nalt
! 	    Il se peut qu'on ait produit entre 2 appels de atmos...
 	    Ne(ialt) = max(denelc(ialt),denselcalc(ialt))
 	  enddo
       	  call conductivite(nalt,altkm,Ne,densneut,
     .		tneutre,temelc,temion,chideg,glat,glong,year,day,
     .		hrloc,cped,chal,ratHoP,cpedsum,chalsum,cpedCS,chalCS,
     .		ratHoPsum,gyreave,gyriave,collOp,collNOp,collO2p,
     .		collionSN,collionRG,collen2,colleo2,colleo1,colle,iprt,
     .		f107(1),icolin)
 	endif
!
 	if(iprt(9).eq.1)then
 	  write(fic_transout,*)
 	  write(fic_transout,1305)
  	  do ien=1,nen,ncountE
 	    write(fic_transout,1300) e(ien),qxdown(ien,-nango2),
     .	     qntsty(ien,1,-nango2),qxup(ien,nango2),qntsty(ien,1,nango2)
 	  enddo
 	endif
!
! 	Computes the theoretical range.
 	avogadro = 6.023e+23
 	if (ezero.ne.0.)then
 	  E0 = ezero/1000.
      	  thrange = 4.3e-7+5.36e-6*E0**1.67-3.8e-9*E0**(-0.7)
 	  write(6,7057)thrange
 	  write(fic_transout,7057)thrange
 	endif
! 	calculates the altitude in g/cm2
 	do iz = 1,nalt
 	  range(iz)=(alt(1)-alt(iz))*densneut(1,iz)*atomas(1)/avogadro
 	enddo
!
 	if(iprt(21).eq.1 .or. iprt(22).eq.1)then
! 	  Search for the max. production altitude.
 	  do isp = 1,nspec
 	    do iz=1,nalt
 	      twork(iz) = prodionsec(isp,iz)
 	    enddo
 	    call mnmxi(twork,nalt,iprodmax,imin)
 	    if(iprodmax.gt.0)then
 	      prodmax(isp)=twork(iprodmax)
 	      zprodmax(isp)=altkm(iprodmax)
 	      write(fic_transout,7055) specie(isp),zprodmax(isp),
     .			prodmax(isp)
 	    endif
 	  enddo
 	endif
 
*	final printing
! 	calcul du libre parcours moyen
 	if(iprt(26).eq.1)then
 	write(fic_transout,1010)e(5),e(nen)
        do ialt=1,nalt
          do ien = 1,nen
            lpm(ialt,ien) = 0.
            do isp = 1,nspec
              lpm(ialt,ien) = lpm(ialt,ien) + 
     .                  (cin(isp,ien)+cel(isp,ien))*densneut(isp,ialt)
            enddo
            lpm(ialt,ien)= 1./lpm(ialt,ien) 
          enddo
 	  write(fic_transout,1020)altkm(ialt),lpm(ialt,5),lpm(ialt,nen)
        enddo
	endif

 
! 	Il faut revoir le calcul de l'energie qui ne veut plus rien 
! 	dire apres toutes les modifs faites
!	if(kiappel.eq.1) write(6,7010) qtoteV,qsump,entot,shsumtot,
!    .		      elhsum
! 	write(6,7018)word,abs(relerr)
  	write(6,*)
  	write(6,*)'Kinetic Transport execute success'

	write(fic_transout,*)'Il faut revoir le calcul de l''energie'
    	write(fic_transout,*)'qui ne veut plus rien dire apres '
 	write(fic_transout,*)'toutes les modifs'
	write(fic_transout,7010)qtoteV,qsump,entot,shsumtot,elhsum
 	write(fic_transout,7018)word,abs(relerr)
	write(fic_transout,7012) elos2,elos,elos1
!
 	  iyd = ifix(year)
 	  if(iyd.gt.99)iyd=iyd-1900
	  if(jpreci.eq.1 .or. jpreci.eq.3 .or. jpreci.eq.4)then
            sun='off'
            write(messag,1000)knm,nen,nang,nalt,sun,iyd,ifix(day),
     .                      hrloc,glat,ap(1),f107(1),f107(3),tempexo
	  else
            sun=' on'
            write(messag,1001)knm,nen,nang,nalt,sun,iyd,ifix(day),
     .                      hrloc,glat,ap(1),f107(1),chideg,tempexo
	  endif
!
 	close (ifeltrans)
 	close (fic_transout)
 	close (irdtin)
!
! 	Remise a zero des conditions initiales de sections efficaces
! 	differentielles.
        if (kiappel.eq.1) call ecrit_DATDEG

 	return
!
993	print*,' Diff.cross-sect. file is in error. Status=',iost
	stop 'error'		! call abort
992	print*,' Cross-section file is in error. Status=',iost
 	print*,' Message de trans.f'
	stop 'error'		! call abort
997	print*,' Option file is in error. Status=',iost
	stop 'error'		! call abort
996	print*, ' Error with assigning input file for external source.'
	stop 'error'		! call abort
!
	end
! 
! ----------------------------------------------------------------------
!
**********************************************************************
*	Subroutines for the electron transport code TRANS
*	in addition the these routines the DISORT and MSIS packages
*	are needed. DISORT makes reference to some LINPACK routines
*	which should reside on the system. If they are not available,
*	the needed routines from LINPACK are extracted from the LINPACK
*	library and reside on a separate file of the TRANS package.
*	Version for Cray (single precision)
**********************************************************************
! 
!----------------------------------------------------------------------
! 
	subroutine eloss(eng,n,temelc,denelc,elosse,nbren,nen,nbralt,
     .			 nalt)
*
*	This subroutine calculates the loss function (continuous slowing
*	down approximation) due to electron-electron and Coulomb
* 	interaction, (See Stamnes,Rees: Heating of thermal ionospheric
*	electrons....,JRL, 10, 309-312,1983) and the backscatter ratio	
*	for e-e interaction						
* 	The variable ELOSSE contains the loss function : ([e]=e-density)
*		ELOSSE = [e] * L(E)					
*
	real temelc(nbralt),denelc(nbralt),elosse(nbren,nbralt)
        real dE_dt
 
*	compute electron-electron interaction loss function		
 
	boltz=8.6173e-05    ! en. (eV) associated with 1 degree Kelvin *
	do m=1,nalt
	  elosse(n,m)=0.
	  eth=boltz*temelc(m)
	  if(eng.le.1500.e3.and.eng.gt.eth) then
	    fac=((eng-eth)/(eng-.53*eth))**2.36
	    elosse(n,m)=fac*3.37e-12*denelc(m)**.97/eng**.94
	    if(fac.le.0.) elosse(n,m)=0.
	  end if
!          elosse(n,m)=dE_dt(eng,denelc(m)*1.e6,temelc(m))/1.6e-11
 	enddo
! 
	return
	end
!
!--------------------------------------------------------------------
!
	subroutine endep(nspec,nalt,e,cinex,engdd,ethres,bratio,qntsty,
     .	    alt,smgdpa,densneut,jsg,jsp,shsum,enrate,prate,ntherm,
     .	    nen)
 
!
 	include 'TRANSPORT.INC'
!
	real cinex(nbrsp,nbrexc,nbren),ethres(nbrsp,nbrexc,nbrionst),
     .	     engdd(nbren),bratio(nbrionst,nbrsp),
     .       qntsty(nbren,2*nbralt-1,-nbrango2:nbrango2),alt(nbralt),
     .       densneut(8,nbralt),shsum(nbrsp),smgdpa(nbralt),
     .	     enrate(nbrsp*2,nbrexc+1,nbralt+1),e(nbren),
     .       prate(nbrsp*2,nbrexc+1,nbralt+1)
	integer ntherm(nbralt),jsg(nbrsp*2+1),jsp(nbrsp*2+1)
 	integer nspec,isp
!
! 	Calculate excitation and ionisation rates (PRATE) and energy 
! 	deposision in excited states (ENRATE)
! 	the "state" JSG(isp) contains the ionization rate
! 	the "state" JSG(isp)+1 contains the sum of all states (including
! 	ionization).
! 	the "altitude level" nalt+1 contains the integral over altitude
 
      call zeroit(enrate,nbrsp*2*(nbrexc+1)*(nbralt+1))
      call zeroit(prate,nbrsp*2*(nbrexc+1)*(nbralt+1))
      do isp = 1,nspec
	shsum(isp)=0.
!
! 	prate = production rate
! 	     1= N2
! 	     2= O2
! 	     3= O
! 	     4= H
! 	     5= He
! 	     6= N
! 	 
	do ialt=1,nalt
	  do istate=1,jsg(isp)
	    do ien=ntherm(ialt),nen
! 	      ethres(isp,ist,1)is the smallest excitation threshold.
	      cros=cinex(isp,istate,ien)*ethres(isp,istate,1)
	      if(istate.eq.jsg(isp)) then
	        cros=0.
	        do jp=1,jsp(isp)
   	          cros=cros+cinex(isp,jsg(isp),ien)*bratio(jp,isp)
     .			*ethres(isp,jsg(isp),jp)
 	        enddo
	      end if
	      prate(isp,istate,ialt)=prate(isp,istate,ialt) +
     .	         qntsty(ien,2*ialt-1,0)*cinex(isp,istate,ien)*engdd(ien)
   	      enrate(isp,istate,ialt)=enrate(isp,istate,ialt) +
     .		 qntsty(ien,2*ialt-1,0)*cros*engdd(ien)
 	    enddo		! boucle sur les energies
	    prate(isp,istate,ialt)=
     .			prate(isp,istate,ialt)*densneut(isp,ialt)
   	    enrate(isp,istate,ialt)=
     .			enrate(isp,istate,ialt)*densneut(isp,ialt)
 	  enddo			! boucle sur les etats d'excitation
 	  if (isp.eq.1)then
	    do ien=ntherm(ialt),nen
!
! 	      On tient compte des ionisations pre-dissociatives.
!             NI = N+. Issue de 50% B'1Su(#12) + 50% sum-sngl (#8) + 
! 	      5% B1pu (#11). Il manque !'4S+u. Est ce sum.sngl  
! 	      (= Rydberg ?)?
 	      prate (6,jsg(6),ialt) = prate(1,12,ialt)*0.5 + 
     .               prate(1,8,ialt) * 0.5 + prate(1,11,ialt)*0.05
 	      enrate (6,jsg(6),ialt) = enrate(1,12,ialt)*0.5 + 
     .               enrate(1,8,ialt) * 0.5 + enrate(1,11,ialt)*0.05

 	    enddo		! boucle sur les energies
 	  endif			! If l'azote
 	enddo			! boucle sur les altitudes
!
! 	Computes integral over altitudes. Stored in nalt+1
	do ialt=2,nalt
	  ddz=(alt(ialt-1)-alt(ialt))/smgdpa(ialt)
	  do istate=1,jsg(isp)
	    prate(isp,istate,nalt+1)=prate(isp,istate,nalt+1)+
     .	      (prate(isp,istate,ialt-1)+prate(isp,istate,ialt))*ddz/2.
	    enrate(isp,istate,nalt+1)=enrate(isp,istate,nalt+1)+
     .	      (enrate(isp,istate,ialt-1)+enrate(isp,istate,ialt))*ddz/2.
   	  enddo
 	enddo
	do istate=1,jsg(isp)
   	  shsum(isp)=shsum(isp)+enrate(isp,istate,nalt+1)
 	enddo
!
! 	Computes sum over all states. Stored in jsp(isp)+1
	do ialt=1,nalt
	  enrate(isp,jsg(isp)+1,ialt)=0.
	  prate(isp,jsg(isp)+1,ialt)=0.
	  do istate=1,jsg(isp)
	    prate(isp,jsg(isp)+1,ialt)=
     .		prate(isp,jsg(isp)+1,ialt)+prate(isp,istate,ialt)
   	    enrate(isp,jsg(isp)+1,ialt)=
     .		enrate(isp,jsg(isp)+1,ialt)+enrate(isp,istate,ialt)
 	  enddo
 	enddo
      enddo			! Boucle sur les especes
 
      return
      end
 
!----------------------------------------------------------------------
 
	subroutine intgrl(e1,e2,f,e,de,sum,sume,nen,m)
!
 	include 'TRANSPORT.INC'
	real f(nbren,nbralt),e(nbren),de(nbren)
!
	n2=nlevtrans(e,de,e2,nen)
	n1=nlevtrans(e,de,e1,n2)
	sum=0.
	sume=0.
	do n=n1,n2,1
	  sum=sum+f(n,m)*de(n)
   	  sume=sume+f(n,m)*e(n)*de(n)
 	enddo
!
	return
	end
!----------------------------------------------------------------------
!
 	subroutine moments(nalt,alt,altkm,nen,e,engdd,
     .		nang,nango2,weitang,cosang,intensite,iprt,ntherm,
     .		Ne_supra,courant_supra,Te_supra,Chaleur_supra)
!
 	implicit none
!
        integer npt
!
 	include 'TRANSPORT.INC'
!
! 	Entrees :
! 	---------
	real engdd(nbren),alt(nbralt),altkm(nbralt),e(nbren)
        real weitang(2*nbrango2),cosang(2*nbrango2)
 	real intensite(nbren,nbralt,-nbrango2:nbrango2)
	integer nalt,nen,iprt(40),nang,nango2
	integer ntherm(nbralt)
	logical flux_flg

! 	Sorties :
! 	---------
! 	Ne_supra 	= 1er  moment 	[cm-3]
! 	courant_supra 	= 2eme moment 	[cm-2.s-1]
! 			  POSITIF VERS LE BAS 
! 	Te_supra 	= 3eme moment 	[K] 
! 	Chaleur_supra 	= 4eme moment 	[eV.cm-2.s-1]
! 			  POSITIVE VERS LE BAS 
! 	
     	real Ne_supra(nbralt),courant_supra(nbralt),Te_supra(nbralt),
     .		Chaleur_supra(nbralt)

! 	Parametres internes :
! 	---------------------
 	integer ialt,ien,iang
 	real trav(nbren),me,boltz,pi
 	real vmoy(nbralt) 		!vitesse moyenne = courant/Ne
! 	Flux stationnaire signe : + vers le haut, - vers le bas
! 	Projections en angle d'ordre 0, 1 et 2 [eV-1.cm-2.s-1] :
! 	fproj0(E,z) = 2pi*sum(intensite*dmu)
! 	fproj1(E,z) = 2pi*sum(mu*intensite*dmu)
! 	fproj2(E,z) = 2pi*sum(mu**2*intensite*dmu)
 	real fproj0(nbren,nbralt),fproj1(nbren,nbralt),
     .		fproj2(nbren,nbralt)
!
        real Qetopcalc
        common/fluxtopcalc/Qetopcalc
        real Qealt(nbralt)
!
! 	Attention aux unites ! La masse electronique est exprimee
! 	en eV.cm-2.s2
! 	La constante de Boltzmann est en eV.K-1
 	data me /5.69E-16/
 	data boltz /8.6173E-05/
!
 	pi = 4.*atan(1.)
!
! 	calcul des flux signes projetes :
!       fproj0 = sum(intensite*dmu)
!       fproj1 = sum(mu*intensite*dmu)
!       fproj2 = sum(mu**2*intensite*dmu)
!

	call zeroit(Ne_supra,nalt)
	call zeroit(Te_supra,nalt)
	call zeroit(courant_supra,nalt)
	call zeroit(chaleur_supra,nalt)
	
	open(76,file=data_path(1:lpath_data)
     &                     //'dir.cine/flux.flag',
     &		form='formatted',status='unknown')
	read(76,*)flux_flg
	close(76)
	if (flux_flg) then
	  open(77,file=chemin(1:lpath)
     &                       //'dir.output/flux.output',
     &		form='formatted',status='unknown')

          do ialt = 1,nalt
	    write(77,*) nalt,nen,nango2,ntherm(ialt)
	    write(77,*) alt(ialt)
	    write(77,*) (cosang(iang),iang=1,nango2),
     &			(weitang(iang),iang=1,nango2)
	    do ien=1,nen
	      write(77,*)e(ien),
     &			 (intensite(ien,ialt,iang),iang=-nango2,-1),
     &		         (intensite(ien,ialt,iang),iang=1,nango2)
	    enddo
	  enddo
	  close(77)
	endif


        do ialt = 1,nalt
          do ien = ntherm(ialt),nen
            fproj0(ien,ialt) = 0.
            fproj1(ien,ialt) = 0.
!            fproj2(ien,ialt) = 0.
            do iang = -nango2,-1,1
! 	      Descendant : cosinus positif
              fproj0(ien,ialt)=fproj0(ien,ialt)+intensite(ien,ialt,iang)
     .			*weitang(iabs(iang))
              fproj1(ien,ialt)=fproj1(ien,ialt)-intensite(ien,ialt,iang)
     .                  *weitang(iabs(iang))*cosang(iabs(iang))
!              fproj2(ien,ialt)=fproj2(ien,ialt)+intensite(ien,ialt,iang)
!     .                  *weitang(iabs(iang))*cosang(iabs(iang))**2
            enddo
            do iang = 1,nango2,1
! 	      montant : cosinus negatif
              fproj0(ien,ialt)=fproj0(ien,ialt)
     .                  +max(0.,intensite(ien,ialt,iang))
     .			*weitang(iabs(iang))
              fproj1(ien,ialt)=fproj1(ien,ialt)
     .                  +max(0.,intensite(ien,ialt,iang))
     .                  *weitang(iabs(iang))*cosang(iabs(iang))
!              fproj2(ien,ialt)=fproj2(ien,ialt)
!     .                  +max(0.,intensite(ien,ialt,iang))
!     .                  *weitang(iabs(iang))*cosang(iabs(iang))**2
            enddo
            fproj0(ien,ialt)=fproj0(ien,ialt)*2.*pi
            fproj1(ien,ialt)=fproj1(ien,ialt)*2.*pi
            fproj2(ien,ialt)=fproj2(ien,ialt)*2.*pi
          enddo
        enddo
!
	do ialt = 1,nalt
  	  do ien = ntherm(ialt),nen
 	    Ne_supra(ialt) = Ne_supra(ialt) +
     .		fproj0(ien,ialt)/sqrt(e(ien)) *engdd(ien)
 	    Te_supra(ialt) = Te_supra(ialt) +
     .		fproj0(ien,ialt)*sqrt(e(ien)) *engdd(ien)
 	    courant_supra(ialt) = courant_supra(ialt) +
     .		fproj1(ien,ialt)*engdd(ien)
 	    chaleur_supra(ialt) = chaleur_supra(ialt) +
     .		fproj1(ien,ialt)*engdd(ien)*e(ien)
 	  enddo
!
 	  Ne_supra(ialt) = Ne_supra(ialt)*sqrt(me/2.)
 	  vmoy(ialt) = courant_supra(ialt)/Ne_supra(ialt)
!
 	enddo
!
 	do ialt = 1,nalt
! 	  Attention ! a basse altitude, on peut avoir des valeurs 
! 	  de Te_supra negatives (faibles).... 
! 	  Quelle conduite pour le code couple ?
!
 	  Te_supra(ialt) = Te_supra(ialt)/Ne_supra(ialt)
     .			*sqrt(2.*me)/3./boltz
!
 	enddo
! 	  
 	if(iprt(25).eq.1)then
 	  write (fic_transout,1000)
 	  do ialt = 1,nalt
 	    write(fic_transout,1010)altkm(ialt),Ne_supra(ialt),
     .		courant_supra(ialt),Te_supra(ialt),Chaleur_supra(ialt)
 	  enddo
 	endif
!
!       division par 2 pour isotropie
        Qetopcalc = -1.* Chaleur_supra(1)*1.6022e-12/2.
        Qetopcalc = Qetopcalc/2.
!
1000	format(10x,'Moments suprathermiques',/,10x,23('-'),/,
     .	'Altitude    Ne_supra courant_supra Te_supra Chaleur_supra',/,
     .  '  [km]       [cm-3]    [cm-2.s-1]    [K]    [eV.cm-2.s-1]')
1010 	format(f10.2,4(1pe10.2))
 	  
!
       	return
       	end
!
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
	subroutine mstream(eng,ien,nang,maxlyr,fiso,tau,utau,ssalb,
     .		albedo,gls,eps,ldeltam,lpr1,lpr2,linear,iphase,gphase,
     .		qint,gaupin,fldn,flup,nalt,lamber,onlyfl,exsorc,usrang,
     .		usrtau,qprim)
!
!	"driver" for the 1990 DISORT version.  This subroutine sets up
!	the input for disort and reorders the output to make it 
!	compatible
!	with the electron transport code.
!
 	include 'TRANSPORT.INC'
	parameter (maxphi=1)

*	input/output variables
 
	real fiso(-nbrango2:-1), tau(0:nbralt-1), ssalb(nbralt-1),
     .	     gls(0:2*nbrango2), qint(2*nbralt-1,-nbrango2:nbrango2),
     .	     gaupin(2*nbralt-1,-nbrango2:nbrango2),fldn(nbralt),
     .	     flup(nbralt)
	real qprim(nbren,nbralt,-nbrango2:nbrango2)
	logical lpr1,lpr2,ldeltam
 
*	internal (active) variables
 
	real dtauc(nbralt-1), pmom(0:2*nbrango2,nbralt-1), 
     .	      umu(2*nbrango2),uou(2*nbrango2,2*nbralt-1),
     .	      uavg(2*nbralt-1),utau(2*nbralt-1),wkfldn(2*nbralt-1), 
     .        wkflup(2*nbralt-1)
	logical prnt(7),exsorc
 
*	internal (dummy) variables
 
	real hl(0:2*nbrango2),phi(maxphi),src(3*(nbralt-1),2*nbrango2),
     .		srcu(3*(nbralt-1),2*nbrango2),rfldir(2*nbralt-1), 
     .	        dfdt(2*nbralt-1), uu(2*nbrango2,2*nbralt-1,maxphi)
 
	character*8 header
	logical lamber, onlyfl, usrang, usrtau, lpass, linear
	save lpass
	data lpass /.false./
	data pi4/12.56637061435917295385057/		! 4*pi
 
	call zeroit(hl,2*nbrango2+1)
	call zeroit(phi,maxphi)
	ntau=2*maxlyr+1
	numu=0
	call zeroit(umu,2*nbrango2)
	nphi=0
	fbeam=0.
	umu0=1.
	phi0=0.
	accur=0.
	if(lpass) then
	  prnt(1)=lpr1	! print "input" only the first time
	  lpass=.false.
	else
	  prnt(1)=.false.
	end if
	prnt(2)=lpr1
	prnt(3)=lpr2
	do i=4,7
   	 prnt(i)=.false.
 	enddo
	header=' '
	do ilyr=1,maxlyr
   	  dtauc(ilyr)=tau(ilyr)-tau(ilyr-1)		! set up tau
 	enddo
 
*	set up phase function

        if(iphase.eq.-1) then           ! Rutherford phase function
          pmom(0,1)=1.
          pmom(1,1)=1. + 2.*eps*(1.-(1.+eps)*log(1.+1./eps))
          do k=2,nang,1
            pmom(k,1)=((2.*k-1.)*(1.+2.*eps)*pmom(k-1,1) -
     .                  k*pmom(k-2,1))/(k-1.)
          enddo
        else if(iphase.eq.2) then               ! Henyey Greenstein
          pmom(0,1)=1.
          pmom(1,1)=gphase
          do k=2,nang,1
            pmom(k,1)=pmom(k-1,1)*gphase
          enddo
        else if(iphase.ge.1) then               ! isotropic
          pmom(0,1)=1.
          do k=1,nang,1
            pmom(k,1)=0.
          enddo
        else if(iphase.eq.0) then      ! iphase=0 => phase fct. supplied
          pmom(0,1)=1.
          do k=1,nang,1
            pmom(k,1)=gls(k)
          enddo
        else
          print*,'*** error in IPHASE: ',iphase
          stop '***error'
        end if
!
        do  k=0,nang,1
          do ilyr=2,maxlyr
            pmom(k,ilyr)=pmom(k,1)      ! same phase fct for every alt.
          enddo
        enddo
 
 
	nango2=nang/2
! 	The photoelectron source function is considered as regularly
! 	distributed on the considered layer : One half is attibuted
! 	to the top (from top to middle), one half to the middle (from
! 	middle to bottom)
	do m=1,maxlyr			! maxlyr = nalt-1
	  ms=3*(m-1)
	  do k=1,nango2,1
! 	    Top of layer
 	    quelle = qint(2*m-1,k)+qprim(ien,m,k)
	    src(ms+1,nango2+k)=max(quelle,0.)
 	    quelle = qint(2*m-1,-k)+qprim(ien,m,-k)
	    src(ms+1,nango2+1-k)=max(quelle,0.)
! 	    Center of layer
!	    quelle = qint(2*m,k)+qprim(ien,m,k)/2.
!    	    src(ms+2,nango2+k)=max(quelle,0.)
!	    quelle = qint(2*m,-k)+qprim(ien,m,-k)/2.
!    	    src(ms+2,nango2+1-k)=max(quelle,0.)
! 	    Bottom of layer
 	    quelle = qint(2*m+1,k)+qprim(ien,m+1,k)
	    src(ms+3,nango2+k)=max(quelle,0.)
 	    quelle = qint(2*m+1,-k)+qprim(ien,m+1,-k)
	    src(ms+3,nango2+1-k)=max(quelle,0.)
 	  enddo
 	enddo
	scru=0.
	bsrc=0.
	tsrc=0.
 
	call disort(ien,maxlyr,dtauc,ssalb,pmom,
     .		usrtau, ntau, utau, nang, usrang,
     .		numu, umu, nphi, phi, fbeam, umu0, src, srcu,
     .		phi0, fiso, lamber, albedo, hl, bsrc,
     .		tsrc, ldeltam, exsorc, onlyfl,
     .		accur, prnt, header, nbralt-1, 2*nbralt-1,
     .		2*nbrango2, 2*nbrango2, maxphi, rfldir, wkfldn,
     .		wkflup, dfdt, uavg, uu, uou, nbrango2, linear, 
     .		fic_transout)
!
	do mlyr=1,maxlyr
	  mt=2*mlyr-1				! top of the layer
	  mc=2*mlyr				! center of the layer
! 	  Dans DISORT, uavg est divise par 4pi (et multiplie par 2pi).
! 	  On retabli la multipication par 2pi.
	  gaupin(mt,0)=uavg(mt)*pi4
	  gaupin(mc,0)=uavg(mc)*pi4
	  fldn(mlyr)=wkfldn(mt)
	  flup(mlyr)=wkflup(mt)
	  do k=1,nango2,1
	    gaupin(mt,k)=uou(nango2+k,mt)	
	    gaupin(mt,-k)=uou(nango2+1-k,mt)
	    gaupin(mc,k)=uou(nango2+k,mc)	
	    gaupin(mc,-k)=uou(nango2+1-k,mc)
 	  enddo
 	enddo
!
	mb=2*maxlyr+1				! bottom of last layer
	gaupin(mb,0)=uavg(mb)*pi4
	fldn(maxlyr+1)=wkfldn(mb)
	flup(maxlyr+1)=wkflup(mb)
	do 402 k=1,nango2,1
	gaupin(mb,k)=uou(nango2+k,mb)	
	gaupin(mb,-k)=uou(nango2+1-k,mb)
402	continue
	return
	end
!
!--------------------------- phel ---------------------------------
! 
        subroutine phel(e,alt,altkm,iprt,fluxprim,
     .	      primelec,qprim,photelec,nspec,nen,nang,nango2,nalt,
     .	      jpreci,mcount,prodionprim,densig)
!
	include 'TRANSPORT.INC'
!
!       read in primary electron production rates.  determine angular
!       dependence of primaries.  calculate source function qprim.
! 	On interpole la production primaire prophel(z,ephel) sur la 
! 	production primaire primelec(alt,e). Puis primelec est 
! 	transforme en flux fluxprim(alt,e) qui est distribue
! 	sur la fonction source. Reecrit (jl,12/1989).
!
! 	Lus dans feltrans :
! 		prophel (z,E)		: cm-3.s-1.eV-1
!		proelec(z)		: cm-3.s-1
!		prodion(z,isp+1)	: cm-3.s-1 (4 = N+)
! 	Outputs interpolees :
! 		primelec(z,E) 		: cm-3.s-1.eV-1
! 		photelec(z) 		: cm-3.s-1
! 		prodionprim(isp+1,z)	: cm-3.s-1 (4 = N+)
! 	Outputs deduites des interpolations :
! 		fluxprim(z,E)		: cm-2.s-1.eV-1
! 		qprim(E,z,A)		: cm-2.s-1.eV-1.sr-1
!
! 	INPUTS
 	integer iprt(40),mcount(5)
        real e(nbren),alt(nbralt),altkm(nbralt),densneut(8,nbralt)
	real  cel(nbrsp,nbren),cin(nbrsp,nbren)
        integer nen,nang,nalt,nango2
! 	INTERNAL PARAMETERS
	real pi,densig(nbralt,nbren)
       	real ephel(nbren),prophel(nbralt,nbren),z(nbralt)
       	real finput(nbren),foutput(nbren),finter(nbren,nbralt)
 	real proelec(nbralt),prodion(nbralt,nbrsp*2)
! 	OUTPUTS
 	real photelec(nbralt),prodionprim(nbrsp*2,nbralt)
       	real primelec(nbralt,nbren),fluxprim(nbralt,nbren)
	real qprim(nbren,nbralt,-nbrango2:nbrango2)
!
900  	format (5f10.2)
910  	format (5(1pe10.2))
920  	format (20x,'Altitude:',f10.2,' km')
1000 	format (/,'Input photoelectron production [cm-3.s-1.eV-1]')
1010 	format (/,'# of Energies and altitudes for primary ',/,
     .        'photoelectron fluxes (nnen,nz):',2i10)
1040 	format (/,' Prod. (prophel) issus de felin (cm-3.s-1.eV-1):')
1050 	format (/,'# of Energies and altitudes for secondary ',/,
     .        'electon fluxes (nen,nalt):',2i10)
1060 	format (/,' Primary photoelectron flux : (cm-2.s-1.eV-1) :')
1070 	format (/,' no of e. and alt. for sec. ',
     .        'el. fluxes (nen,nalt):',2i10)
1080 	format ('Primary interpolated prod. :(cm-3.s-1.eV-1):')
!
 	if(jpreci.eq.1 .or. jpreci.eq.3 .or. jpreci.eq.4)then
 	  do ialt = 1,nalt
	    photelec(ialt) = 0.
 	    do ien = 1,nen
 	      primelec(ialt,ien) = 0.
 	      fluxprim(ialt,ien) = 0.
 	      do iang = -nango2,nango2
  		qprim(ien,ialt,iang) = 0.
 	      enddo
 	    enddo
 	  enddo
 	  return
 	endif
      	open (ifeltrans,file=
     .	data_path(1:lpath_data)
     &        //'dir.cine/FELTRANS',
     .	form='unformatted')
      	rewind(ifeltrans)
	pi=atan(1.)*4.
      	r2pi = 1./(pi*2.)
!
!  *  	read in primary photoelectron flux data
!     	nnen   : # energies.
!     	nz     : # altitudes.
!     	ephel   : energies.
!     	z      : altitudes.
!     	prophel    : primary production,cm-3.s-1.eV-1,issus de prime .
!
      	read  (ifeltrans) nnen, nz,nnspec
      	read  (ifeltrans) (ephel(j), j = 1,nnen)
      	read  (ifeltrans) (z(i), i = 1,nz)
      	read  (ifeltrans) ((prophel(i,j),j=1,nnen),i=1,nz)
      	read  (ifeltrans) (proelec(i), i = 1,nz)
      	read  (ifeltrans) ((prodion(ialt,isp),ialt=1,nz),isp=1,nnspec+1)
      	if(iprt(22).eq.1)then
      	  write (fic_transout,1000)
          write (fic_transout,*) 'Altitudes (z)'
          write (fic_transout,900) (z(i), i = 1,nz)
          write (fic_transout,*) 'Energies (ephel)'
          write (fic_transout,910) (ephel(j), j = 1,nnen)
          write (fic_transout,1040)
          do i=1,nz
            write(fic_transout,920)  z(i)
            write (fic_transout,910)(prophel(i,j), j = 1,nnen)
!           [prophel]=cm-3.s-1.eV-1
 	  enddo
      	endif
!
! 	Interpolations en altitude
! 	Production totale
 	call intlin(nz,z,proelec,nalt,altkm,photelec)
! 	Production d'ions par especes.
 	do isp = 1,nnspec+1
 	  do ialt = 1,nz
 	    finput(ialt) = prodion(ialt,isp)
 	  enddo
 	  call intlin(nz,z,finput,nalt,altkm,foutput)
! 	  On va repasser N+ en position 6, comme pour les especes
! 	  neutres
 	  if (isp.le.3)then
 	    do ialt = 1,nalt
 	      prodionprim(isp,ialt) = foutput(ialt)
 	    enddo
 	  else if (isp.eq.4)then
! 	    On avait N+, on le passe en position 6
 	    do ialt = 1,nalt
 	      prodionprim(6,ialt) = foutput(ialt)
 	    enddo
 	  else if (isp.eq.5)then
! 	    On avait H+, on le passe en position 4
 	    do ialt = 1,nalt
 	      prodionprim(4,ialt) = foutput(ialt)
 	    enddo
 	  else if (isp.eq.6)then
! 	    On avait He+, on le passe en position 5
 	    do ialt = 1,nalt
 	      prodionprim(5,ialt) = foutput(ialt)
 	    enddo
 	  endif
 	enddo
! 
!
! 	Interpolation en energie:
	do iz=1,nz
	  do ien=1,nnen
	    finput(ien)=prophel(iz,ien)
 	  enddo
 	  call intlin(nnen,ephel,finput,nen,e,foutput)
	  do ien=1,nen
	    finter(ien,iz)=foutput(ien)
 	  enddo
 	enddo
! 	Interpolation en altitude:
 	do ien=1,nen
	  do iz=1,nz
	    finput(iz)=finter(ien,iz)
 	  enddo
 	  call intlin(nz,z,finput,nalt,altkm,foutput)
	  do ialt=1,nalt
	    primelec(ialt,ien)=foutput(ialt)
 	  enddo
 	enddo
!
      if(iprt(22).eq.1)then
!       print primary production rate.
        write (fic_transout,1070) nen, nalt
        write (fic_transout,*) 'Altitudes (alt)'
        write (fic_transout,900) (altkm(ialt), ialt = 1,nalt)
        write (fic_transout,*) 'Energies (e)'
        write (fic_transout,910) (e(j), j = 1,nen)
        write (fic_transout,1080)
        do ialt=1,nalt,mcount(5)
          write(fic_transout,920)  altkm(ialt)
          write (fic_transout,910)(primelec(ialt,k), k = 1,nen)
 	enddo
      endif
!
! 	On passe en flux.
 	do ien = 1,nen
	  do ialt=1,nalt
  	    fluxprim(ialt,ien)= primelec(ialt,ien)/densig(ialt,ien) 
            do iang=-nango2,nango2
! 	      La division par 2 pi vient de l'integration polaire
              qprim(ien,ialt,iang)= fluxprim(ialt,ien)* r2pi
! 	      On redivise par 2 parce que le flux primaire est 
! 	      isotrope (la moitie vers le haut et la moitie vers le
! 	      bas)
              qprim(ien,ialt,iang)= qprim(ien,ialt,iang)/2. 
 	    enddo
 	  enddo
 	enddo
!
        if(iprt(21).eq.1 .or. iprt(22).eq.1)then
!         print normalized source function.
          write (fic_transout,1050) nen, nalt
          write (fic_transout,*) 'Energies (e)'
          write (fic_transout,910) (e(j), j = 1,nen)
          write (fic_transout,*) 'Altitudes (alt)'
          write (fic_transout,900) (altkm(i), i = 1,nalt)
          write (fic_transout,*) 
     .		'Electron primary production (photelec) [cm-3.s-1]'
          write (fic_transout,900) (photelec(ialt), ialt = 1,nalt)
          write (fic_transout,1060)
          do ialt=1,nalt,mcount(5)
            write(fic_transout,920)  altkm(ialt)
            write(fic_transout,910)(fluxprim(ialt,ien), ien = 1,nen)
 	  enddo
        endif
!
        return
        end
!
!--------------------------------------------------------------------
!
        subroutine qmstr(n,mmax,nspec,engdd,densneut,ctot,elosse,
     .           qint,fint,omdeg,omsec,weitang,nbren,nen,nbralt,nalt,
     .		 nbrsp,nbrango2,nango2,irdtin)
!
!		Accelerated Cray version
!       This subroutine calculates the sources due to degradation 
!       of electrons in energy. (See equ. 2)
!       The second half of the derivatives of the loss functions is 
!       also calculated here.
!
        implicit logical (l)
        real engdd(nbren),densneut(8,nbralt),ctot(2*nbralt-1),
     .	  elosse(nbren,nbralt),
     .    qint(2*nbralt-1,-nbrango2:nbrango2),
     .		 fint(nbren,2*nbralt-1,-nbrango2:nbrango2),
     .    omdeg(nbren,nbrsp),omsec(nbren,nbrsp),weitang(nbrango2)
	real rdp(2*nbralt-1),rds(2*nbralt-1)
	real fsec(2*nbralt-1),fac(2*nbralt-1)
 
        call zeroit(qint,(2*nbrango2+1)*(2*nbralt-1))
	maxlyr=(mmax-1)/2
 
	read(irdtin) ((omdeg(nn,j),nn=nen,n-1,-1),j=1,nspec)
	read(irdtin) ((omsec(nn,j),nn=nen,n-1,-1),j=1,nspec)
 
        do 100 in=nen,n,-1
          dde = engdd(in)
	  call zeroit(rdp,2*nbralt-1)
	  call zeroit(rds,2*nbralt-1)
          do 103 j=1,nspec
            do 109 m=1,maxlyr+1
              rdp(2*m-1) = rdp(2*m-1) + omdeg(in,j)*densneut(j,m)
              rds(2*m-1) = rds(2*m-1) + omsec(in,j)*densneut(j,m)
109	    continue
	  do 103 m=1,maxlyr
	    rdp(2*m) = rdp(2*m) + omdeg(in,j)*
     .			(densneut(j,m)+densneut(j,m+1))*0.5
	    rds(2*m) = rds(2*m) + 
     .			omsec(in,j)*(densneut(j,m)+densneut(j,m+1))*0.5
103       continue
	  do 102 m=1,mmax
            rdp(m) = rdp(m)/ctot(m)*dde    ! definition of Rjm an equ. 9
102         rds(m) = rds(m)/ctot(m)*dde
 
            call zeroit(fsec,2*nbralt-1)
*				! isotropic distribution of secondaries
            do 104 k=-nango2,nango2,1
              if(k.ne.0) then
		do 105 m=1,mmax
	  	  fsec(m)=fsec(m)+fint(in,m,k)*weitang(iabs(k))/2.
105		continue
	      end if
104         continue
 
 	    testmax = 1.e+33
 	    testmin = 1.e-33
	    do 106 m=1,mmax
	      toto = rds(m)*fsec(m)/2.
	      if(toto.lt.testmax .and. toto.gt.testmin)
     .   	    qint(m,0)=qint(m,0) + toto
106	    continue
 
            do 107 k=-nango2,nango2,1
            if(k.ne.0) then
		do 108 m=1,mmax
                toto = rdp(m)*fint(in,m,k)+ rds(m)*fsec(m)
	        if(toto.lt.testmax .and. toto.gt.testmin)
     .	  	  qint(m,k)=qint(m,k) + toto
108		continue
	    end if
107       continue
100     continue
 
	do 201 m=1,maxlyr+1
	  fac(2*m-1)=elosse(n,m)/ctot(2*m-1)/engdd(n-1)
201	continue
	do 203 m=1,maxlyr
	  fac(2*m)=(elosse(n,m)+elosse(n,m+1))*0.5/ctot(2*m)/engdd(n-1)
203	continue
	do 200 k=-nango2,nango2,1
	  if(k.ne.0) then
	    do 202 m=1,mmax
 	      toto = fac(m)*fint(n,m,k)
	      if(toto.lt.testmax .and. toto.gt.testmin)
     .	       		qint(m,k)=qint(m,k) + toto
202	    continue
	  end if
200	continue
 
        return
        end
 
!----------------------------------------------------------------------
 
	function ncross(e,engdd,qntsty,ialt,Etherm,denelc,nen)
!
!	calculates the crossover energy of the thermal and streaming 
!	electrons
!
 	implicit none
!
        integer npt
!
 	include 'TRANSPORT.INC'
!
 	real etherm,emass
 	integer ialt,nen,nmax,nn,nlevtrans,n,ncross
	real qntsty(nbren,2*nbralt-1,-nbrango2:nbrango2),e(nbren)
 	real engdd(nbren),denelc
 	real pi,a,b,c,flx
 
        real eth,fzero

        eth=fzero(denelc*1.e6,Etherm*1.16046e4)
        ncross=nlevtrans(e,engdd,5.*Etherm,nen)
        return

	pi=atan(1.)*4.
	if(denelc.le.1.e-30.or.Etherm.le.0.)then
!		ncross=1
		return
	end if
	emass=5.69e-16
	nmax=nlevtrans(e,engdd,10.,nen)
	nn=ncross
	a=log(emass/(2.*pi*Etherm))*1.5
	do 100 n=nmax,nn,-1
	  b=log(denelc*2.) - 2.*log(emass) - e(n)/Etherm
	  c=a + b + log(e(n))
	  flx=exp(c)
	  if(flx.ge.qntsty(n,2*ialt-1,0)) then
	    ncross=n
	    return
	  end if
100	continue
!	print*,'>> ncross error.',ialt,flx,qntsty(nmax,2*ialt-1,0),
!     .				nn,nmax
!	ncross=1
	return
	end
!
!----------------------------------------------------------------------
!
	subroutine porter(g,eng,nang)
	real g(0:nang),eng
	real e(0:16),gama(16),beta(16),alfa(16)
	data e/0.,2.3,2.9,3.9,5.,10.,15.,20.,25.,30.,50.,90.,
     .		100.,300.,500.,800.,1000./
	data gama/2.68e-1,3.32e-1,2.47e-1,3.81e-1,2.28e-1,1.27e-1,
     .		7.59e-2,6.23e-2,5.55e-2,2.95e-2,1.76e-2,1.61e-2,
     .		8.60e-3,6.70e-3,5.09e-3,4.21e-3/
	data beta/5.24e1,9.13e-1,5.39e-1,7.34e-1,4.66e-1,6.47e-1,
     .		6.98e-1,6.74e-1,7.11e-1,1.06,6.12,8.64,2.37e1,8.32,
     .		3.34,1.97/
	data alfa/2.7,4.62e-1,3.17e-1,4.77e-1,1.92e-1,1.32e-1,1.32e-1,
     .		1.24e-1,8.07e-2,1.19e-1,3.78e-1,4.73e-1,-1.71,-1.27,
     .		-1.12,-1.08/
!
!	alternate phase function to replace Rutherford for low energy
!	see Porter et al., 1987
!
	if(nang.gt.32) then
		print*, 'Arrays in PORTER too small'
		stop '***error'
	end if
	do 400 i=1,16
	if(e(i-1).lt.eng.and.eng.le.e(i)) then
		if(nang.gt.2) then
			call recdown(g,nang,alfa(i),beta(i),gama(i))
		else
			call recup(g,nang,alfa(i),beta(i),gama(i))
		end if
		return
	end if
400	continue
	do 410 k=1,nang		! outside the energy range: isotropic
410	g(k)=0.
	g(0)=1.
	return
	end
 
!----------------------------------------------------------------------
	subroutine recdown(g,nang,alfa,beta,gama)
	parameter (maxstr=32,madd=9,mdim=maxstr+madd)
	real g1(0:mdim),g2(0:mdim),g(0:nang),alfa,beta,gama
!
!	downward recursion for Porter's phase function
!
	if(beta.eq.0.) then
		if(gama.eq.0.) then		! isotropic
			g(0)=1.
			do 100 k=1,nang
100			g(k)=0.
			return
		end if				! Rutherford
		fnorm=1./(4.*gama*(1.+gama))
		fnorm1=1.
		fnorm2=0.
	else					! Porter
		fnorm=0.25/(gama*(1.+gama))+
     .			0.25*beta/(alfa*(1.+alfa))
		fnorm1=1./(4.*gama*(1.+gama))/fnorm
		fnorm2=beta/(4.*alfa*(1.+alfa))/fnorm
	end if
	infty=mdim
!	infty=nang+madd	      ! is this better? faster? still accuarate?
	del=1./(nang+madd)**4
	call zeroit(g2,mdim+1)
430	continue
	g1(infty)=0.
	g1(infty-1)=g1(infty)+del
	g2(infty)=0.
	if(beta.ne.0.) g2(infty-1)=g2(infty)+del
	do 410 k=infty-1,1,-1
	g1(k-1)=((2.*k+1.)*(1.+2.*gama)*g1(k)-k*g1(k+1))/(k+1.)
	if(beta.ne.0.) g2(k-1)=(-(2.*k+1.)*(1.+2.*alfa)*g2(k)-
     .			k*g2(k+1))/(k+1.)
	if(g1(k-1).gt.1e10.or.g2(k-1).gt.1e10) then
		infty=infty-1
		if(infty.le.1) then
	 	   print*,'error in recursion >>> too many iterations'
		   stop '***error'
		end if
		goto 430
	end if
410 	continue
	do 400 k=nang,0,-1
	g1(k)=g1(k)*fnorm1/g1(0)
	if(fnorm2.ne.0.) then
		g2(k)=g2(k)*fnorm2/g2(0)
	else
		g2(k)=0.
	end if
400	continue
	do 420 k=0,nang
	g(k)=g1(k)+g2(k)
420	continue
	return
	end
 
!----------------------------------------------------------------------
 
	subroutine recup(g,nang,alfa,beta,gama)
	parameter (maxstr=32,mdim=maxstr)
	real g1(0:mdim),g2(0:mdim),g(0:nang),alfa,beta,gama
!
! 	use recursion relation upwards
!
	if(beta.eq.0.) then
		if(gama.eq.0.) then			! isotropic
			g(0)=1.
			do 100 k=1,nang
100			g(k)=0.
			return
		end if
		fnorm=1./(4.*gama*(1.+gama))   		! rutherford
	else
		fnorm=(1./(gama*(1.+gama))+    		! porter
     .			beta/(alfa*(1.+alfa)))/4.
	end if
	g1(0)=1./(4.*gama*(1.+gama))/fnorm
	g1(1)=(1.+2.*gama)/(4.*gama*(1.+gama))-
     .		.5*log(1.+1./gama)
	g1(1)=g1(1)/fnorm
	call zeroit(g2,mdim+1)
	if(beta.ne.0.) then
		g2(0)=beta/(4.*alfa*(1.+alfa))/fnorm
		g2(1)=(1.+2.*alfa)/(4.*alfa*(1.+alfa))-
     .			.5*log(1.+1./alfa)
		g2(1)=-beta*g2(1)/fnorm
	end if
	do 410 k=2,nang,1
	g1(k)=((2.*k-1.)*(1.+2.*gama)*g1(k-1)-k*g1(k-2))/(k-1.)
	if(beta.ne.0.) g2(k)=(-(2.*k-1.)*(1.+2.*alfa)*g2(k-1)-
     .			k*g2(k-2))/(k-1.)
410	continue
	do 420 k=0,nang,1
	g(k)=g1(k)+g2(k)
420	continue
	return
	end
 
!----------------------------------------------------------------------
 
	subroutine ruther(g,eng,nang)
!
!	screened Rutherford phase function
!	corr is an adjustment to the screening parameter eps
!
	real g(0:nang),eps,eng
	if(eng.ge.1000.) then
		corr=0.6*(1000./eng)**0.09
		eps=6.22e-5/(2+eng/511000.)/(eng/511000.)*corr
		call recup(g,nang,0.,0.,eps)
		return
	else if(eng.ge.477.) then
		corr=.6
		eps=6.22e-5/(2+eng/511000.)/(eng/511000.)*corr
		call recdown(g,nang,0.,0.,eps)
		return
	else if(eng.gt.12.) then
		eps=0.5/((eng/12.)**0.75-1.)
		call recdown(g,nang,0.,0.,eps)
		return
	else
		do 100 k=1,nang
100		g(k)=0.
		g(0)=1.
		return
	end if
	end
! 
!----------------------------------------------------------------------
! 
	integer function nlevtrans(e,engdd,eng,nen)
 
!	Function to find the grid number to which a given energy belongs
!	input   ENG : energy
!	output  NLEV : nearest cell to energy ENG
 
 	real e(nen),engdd(nen)
 
	do n=nen,1,-1
	  de=engdd(n)/2.
	  if(eng.ge.e(n)-de) then
	    nlevtrans=n
	    return
	  end if
 	enddo
	if(eng.le.0.) then
	  nlevtrans=1
	  return
	end if
	nlevtrans=1
 	return
!
	end
!
!---------------------------- densout -------------------------------
! 
        subroutine densout(hrloc,nalt,altkm,prodeltot,denelc,denselcalc)
!
        include 'TRANSPORT.INC'
!
        real prodeltot(nbralt),altkm(nbralt),denelc(nbralt)
        real denselcalc(nbralt)
        real prod,hrloc
        integer ialt,nalt
!
!       calcul de densites
        do ialt=1,nalt
          if (altkm(ialt).lt.85.)then
            alphaeff=7.30e+04*exp(-altkm(ialt)/3.3)
            prod = prodeltot(ialt)
            denselcalc(ialt)=sqrt(prod/alphaeff)      
          elseif (altkm(ialt).ge.85. .and. altkm(ialt).le.180.)then
            alphaeff=2.50e-06*exp(-altkm(ialt)/51.2)
            prod = prodeltot(ialt)
            denselcalc(ialt)=sqrt(prod/alphaeff)      
          else
            denselcalc(ialt) = denelc(ialt)
          endif
        enddo
!
        return
        end
!
!-------------------------- reed ----------------------------
!
      subroutine reed (iprt,idess,mcount,ncountE,ncountA,
     .	   linear,ldeltam,lporter,nspec,e,Ebot,engdd,nen,
     .	   nang,nango2,pitchang,cosang,weitang,angzb,gmu,gwt,
     .	   f107a,f107,smgdpa,day,year,glong,alt,altkm,nalt,
     .	   tneutre,densneut,denelc,colden,title,jpreci,
     .	   zbot,ztop,hrloc,tempexo,knm,glat,Apind,albedo,
     .	   qxdown,qxup,fluxup,fluxdown,temelc,temion,ezero,izplt,ieplt,
     .	   eplt,nke,jsg,jsp,ethres,bratio,cel,cin,cinex,
     .	   lamber,onlyfl,exsorc,usrang,usrtau,ddeng,centE,botE,icolin)
!
	include 'TRANSPORT.INC'
!
 	logical lporter,linear,ldeltam,lamber,onlyfl,exsorc,usrang,
     .		usrtau
 	character*80 crsin,rdtin,crsinput
 	integer iprt(40),mcount(5),ncountE,ncountA,idess(20),
     .		izplt(4),ieplt(4),icolin
 	real eplt(4)
!
 	integer nspec,knmneutral,nalt,jpreci
 	integer isp,ialt,iexc,ilinear,iporter,ien,iost,i,j,iang,jp
	character*30 istdate
 	integer lenc
 	real bid
 	real zbot,ztop,hrloc,year,tempexo,f107,f107a,Apind,day,
     . 	 	glat,glong,albedo,alt(nbralt),
     .		altkm(nbralt),tneutre(nbralt),densneut(8,nbralt),
     .     	colden(8,nbralt)
!
 	character headline*80,rtext*20,phasfct*6
!
 	integer iopal,nen,knmsig
 	real centE(nbren),botE(nbren),ddeng(nbren)
 	real e(nbren),Ebot(nbren),engdd(nbren),esig(nbren),bsig(nbren),
     .  	ethres(nbrsp,nbrexc,nbrionst),bratio(nbrionst,nbrsp),
     .  	cel(nbrsp,nbren), cin(nbrsp,nbren),
     .  	cinex(nbrsp,nbrexc,nbren)
!
 	integer knm,nensig,ne,nang,nango2,jsg(nbrsp*2+1),jsp(nbrsp*2+1)
 	real ezero,trav(nbren)
 	real pitchang(2*nbrango2),cosang(2*nbrango2),weitang(2*nbrango2)
 	real angzb(2*nbrango2),gmu(2*nbrango2),gwt(2*nbrango2)
 	real qxdown(nbren,-nbrango2:-1),qxup(nbren,nbrango2)
 	real fluxdown(nbren,nbrango2),fluxup(nbren,nbrango2)
 	real denelc(nbralt),temelc(nbralt),
     .	       zel(nbralt),smgdpa(nbralt),temion(nbralt)
	character*9 title(nbrexc,nbrsp)
!
	open(fic_transout,file=
     .	data_path(1:lpath_data)
     &        //'dir.cine/TRANSOUT'
     .	,status='unknown')
	rewind fic_transout
! 
!-------------------------------------------------------------
!-------------------- lecture DATTRANS ------------------------
!
	open(fic_datdeg,file=
     .	data_path(1:lpath_data)
     &        //'dir.cine/DATDEG',
     .	status='OLD',iostat=iost,err=997)
 	rewind(fic_datdeg)
 	read (fic_datdeg,*)ibid
 	read(fic_datdeg,'(a)') crsinput
	read(fic_datdeg,'(a)') crsin
 	crsin = crsin(1:lenc(crsin))
	read(fic_datdeg,'(a)') rdtin
 	rdtin = rdtin(1:lenc(rdtin))
 	close(fic_datdeg)
!
	open(fic_dattrans,file= data_path(1:lpath_data)
     &                                //'dir.cine/DATTRANS',
     .	status='OLD', iostat=iost,err=997)
	rewind fic_dattrans
 	call xline(15,fic_dattrans)
 	read(fic_dattrans,*)(iprt(i),i=1,26)
!	MCOUNT  : print every MCOUNT altitude step
!	    (1) : intensity/density	(2) : energy deposition
!	    (3) : excitation rates	(4) : electron heating
 	call xline(3,fic_dattrans)
 	read(fic_dattrans,*)(mcount(i),i=1,5)
! 	ncount = Printing energy and angle steps
 	read(fic_dattrans,*) ncountE,ncountA
 	call xline(19,fic_dattrans)
 	read(fic_dattrans,*)(idess(i),i=1,13)
 	call xline(3,fic_dattrans)
 	read(fic_dattrans,*)(izplt(ialt),ialt=1,4)
 	read(fic_dattrans,*)(ieplt(ien),ien=1,4)
 	read(fic_dattrans,*)(eplt(ien),ien=1,4)
 	linear = .false.
	read(fic_dattrans,*)ilinear
 	if (ilinear.eq.1)linear=.true.
 	lporter = .true.
	read(fic_dattrans,*)iporter
 	if (iporter.eq.0)lporter=.false.
 	ldeltam = .true.
	read(fic_dattrans,*)ideltam
 	if (ideltam.eq.0)ldeltam=.false.
 	lamber = .true.
	read(fic_dattrans,*)ilamber
 	if (ilamber.eq.0)lamber=.false.
 	onlyfl = .true.
 	read(fic_dattrans,*)ionlyfl
 	if (ionlyfl.eq.0)onlyfl=.false.
 	exsorc = .true.
 	read(fic_dattrans,*)iexsorc
 	if (iexsorc.eq.0)exsorc=.false.
 	usrang = .true.
 	read(fic_dattrans,*)iusrang
 	if (iusrang.eq.0)usrang=.false.
 	usrtau = .true.
 	read(fic_dattrans,*)iusrtau
 	if (iusrtau.eq.0)usrtau=.false.
 	read(fic_dattrans,*)icolin
 	close (fic_dattrans)
!
!
!----------------------------------------------------------------
!----------------- lecture de sections efficaces ----------------
!
!  *    read in ionization and excitation cross-sections vs energy.
!
!	print*,' open file: ',data_path(1:lpath_data)
!     &                              //crsin
	open(icrsin,file=data_path(1:lpath_data)
     &                         //crsin,status='OLD',form='UNFORMATTED',
     .		iostat=iost,err=992)
	rewind icrsin
!	print*,' open file: ',data_path(1:lpath_data)
!     &                              //rdtin
	open(irdtin,file=data_path(1:lpath_data)
     &                         //rdtin,status='OLD',form='UNFORMATTED',
     .		iostat=iost,err=993)
	rewind irdtin
	read(irdtin) headline
	iopal=index(headline,'Opal')
	if(iopal.ne.0) then
	  rtext='  '//headline(iopal:iopal+4)
	else
	  rtext='       .'
	end if
	read(irdtin) nensig,nspecsig
 	if (nspecsig.lt.nspec)then
 	  write(6,*)'number of specie in rdtin   =',nspecsig
 	  write(6,*)'number of specie in NEUTRAL =',nspec
 	  write(6,*)'impossible to continue'
 	  stop
 	endif
	read(icrsin) ien,isp,iexc,jp	! test for correct inputfile
 	if(ien.ne.nensig)then
	   print*,rdtin,' does not match ',crsin
	   print*,'#(energie) in ',rdtin,' = ',nensig 
	   print*,'#(energie) in ',crsin,' = ',ien
	   print*,'#(specie) in  ',rdtin,' = ',nspec
	   print*,'#(specie) in  ',crsin,' = ',isp
 	  stop
 	endif

	if(nensig.gt.nbren.or.nspec.gt.nbrsp.or.iexc.gt.nbrexc.or.
     .		jp.gt.nbrionst) then
 	  print*,'actual and maximum values'
	  print*,'nensig,nbren :',nensig,nbren
	  print*,'nspec,nbrsp :',nspec,nbrsp
	  print*,'iexc,nbrexc :',iexc,nbrexc
	  print*,'jp,nbrionst :',jp,nbrionst
	  stop 'error'		! call abort
 	endif
!
	read(icrsin) (esig(ien),ien=1,nensig)
	read(icrsin) (bsig(ien),ien=1,nensig)
	read(icrsin) (engdd(ien),ien=1,nensig)
!	Read cross-sections
	read(icrsin) (jsg(isp),isp=1,nspecsig),(jsp(isp),isp=1,nspecsig)
	read(icrsin) ((title(iexc,isp),isp=1,nspecsig),iexc=1,nbrexc)
  
!Write out the titles for each species and excitation
!        open(45,file='titles.dat',status='replace')
!        do isp=1,nspecsig
!                write(45,*) isp,jsg(isp),jsp(isp)
!                do iexc=1,nbrexc
!                        write(45,*) title(iexc,isp)
!                enddo
!                write(45,*) ''
!        enddo
!        close(45)
        !
  
	read(icrsin)(((ethres(isp,iexc,jp),isp=1,nspecsig),
     .		iexc=1,nbrexc),jp=1,nbrionst)
	read(icrsin) ((bratio(jp,isp),jp=1,nbrionst),isp=1,nspecsig)
	read(icrsin) ((cel(isp,ien),isp=1,nspecsig),ien=1,nensig)
	read(icrsin) ((cin(isp,ien),isp=1,nspecsig),ien=1,nensig)
	read(icrsin) (((cinex(isp,iexc,ien),isp=1,nspecsig),
     .		iexc=1,nbrexc),ien=1,nensig)
 	close (icrsin)

!
! 	On tient compte de l'azote NI (N+)
 	jsg(6) = 1		! On ne produit que N+
 	jsp(6) = 1		! et dans un seul etat
!
! 	The differential cross sections are read in the transport
! 	equation resolution process.
!
!       Met les energies en ordre croissant
        if(centE(1).lt.centE(nen))then
          do ien = 1,nen
            e(ien) = centE(ien)
 	    Ebot(ien) = botE(ien)
            engdd(ien) = ddeng(ien)
          enddo
        else
          do ien = 1,nen
            e(nen+1-ien) = centE(ien)
            Ebot(nen+1-ien) = botE(ien)
            engdd(nen+1-ien) = ddeng(ien)
          enddo
        endif

!       Print out the cross section for a few states of interest
!       O(1D), O(1S),O(3p3P)
        open(955,file='cross.dat',status='replace')
        do ien=1,nen
          write(955,*) e(ien), cinex(3,2,ien), cinex(3,3,ien), 
     &          cinex(3,5,ien),cinex(3,4,ien),cinex(3,jsg(3),ien),
     &          cinex(1,3,ien),cinex(1,4,ien),cinex(1,5,ien),
     &          cinex(1,7,ien),
     &          cinex(1,jsg(1),ien),cinex(2,jsg(2),ien)
        enddo
        close(955)
!       -MZ

!
 	if (e(1) .ne. esig(1)) then
 	  write(6,*)'erreur:sig et trans ne tournent pas sur la '
 	  write(6,*)'meme grille d''entree des energies'
   	  write(6,*)'grille sig:'
 	  write(6,*)nensig
      	  write(6,1000)(esig(ien),ien=1,nensig)
   	  write(6,*)'grille trans:'
 	  write(6,*)nen
      	  write(6,1000)(e(ien),ien=1,nen)
 	  stop
 	endif 
1000    format(5(1pe10.2))
!
	if(nango2.gt.nbrango2) then
	  print*,' Too many streams: nbrango2 = ', nbrango2,
     .		 '                   nango2   = ', nango2
	  stop 'error'		! call abort
	end if
	if(nang.eq.2) lporter=.false.
	if(lporter) then
	  phasfct='Porter'
	else
	  phasfct='Wedde '
	end if
!
!  	Oriente les angles comme il faut pour disort.
! 	qxdown (-nango2--->-1)=flux down 
! 	qxup   (1--->nango2)=flux up 
!
 	do iang = 1,nango2
 	  pitchang(iang) = angzb(nango2+1-iang)
 	  pitchang(nango2+iang) = angzb(iang+nango2)
 	  cosang(iang) = gmu(nango2+1-iang)
 	  cosang(nango2+iang) = gmu(iang+nango2)
 	  weitang(iang) = gwt(nango2+1-iang)
 	  weitang(nango2+iang) = gwt(iang+nango2)
 	  do ien = 1,nen
            qxdown(ien,iang-nango2-1)= fluxdown(ien,iang)
            qxup(ien,nango2+1-iang) = fluxup(ien,iang)
 	  enddo
 	enddo
!
! 	Writting
!
	write(fic_transout,4002) istdate
	write(fic_transout,4003)nbren,nen,nbralt,nalt,2*nbrango2,nang,
     .	  nbrsp,nspec,nbrexc,nbrionst
	write(fic_transout,4004) crsin,rdtin
	if(linear) then
	  write(fic_transout,4009)
	else 
	  write(fic_transout,4008)
	end if
	write(fic_transout,4018) zbot,albedo,phasfct
	write(fic_transout,4019) headline
4002	format(' Electron transport code',/,a30)
4003	format(/' The input parameter statement gives:'/
     .	' max. (nbren)  and actual (nen)  length of E-grid     :',2i5,/,
     .	' max. (nbralt) and actual (nalt) length of alt-grid   :',2i5,/,
     .	' max. (nbrang) and actual (nang) length of angle-grid :',2i5,/,
     .	' max. (nbrsp)  and actual (nspec)number of species    :',2i5,/,
     .	' number of excitation states (nbrexc) :',i5,/,
     .  ' number of ion states (nbrionst) :',i5)
4004	format(/' Cross-section input files              : ',a20,/,
     .          ' Differential cross-section input files : ',a20)
4008	format(/' The energy degradation source',
     .	  ' term is interpolated using an exp-lin',/,' function in tau')
4009	format(' The energy degradation source',
     .		'term is interpolated linearly with tau.')
4018	format(' The albedo at ',f4.0,' km is set to ',f4.2/
     .		' The low energy phase function is taken from ',a6)
4019	format(1x,a70)
!
      	write (fic_transout,5030) zbot, ztop
      	write (fic_transout,5040) hrloc, glat, day
      	write (fic_transout,5050) tempexo, f107a, f107,albedo
      	write (fic_transout,5070) year,nalt
      	write (fic_transout,5080) nen, nang, knm
!
      	write (fic_transout,5090) (altkm(ialt), ialt = 1,nalt)
      	write (fic_transout,5020) (denelc(ialt), ialt = 1,nalt)
      	write (fic_transout,5100) (e(n), n = 1,nen)
      	write (fic_transout,5105) (engdd(n), n = 1,nen)
      	write (fic_transout,5110) (pitchang(iang), iang = 1,nang)
	write (fic_transout,5112) (cosang(iang),iang=1,nang)
	write (fic_transout,5114) (weitang(iang),iang=1,nang)
	write (fic_transout,*)
!
	if(iprt(3).eq.1) then
 	  do isp=1,nspec
	    write(fic_transout,2003) specie(isp)
 	    do ien=1,nen
	      write(fic_transout,2004) e(ien),engdd(ien),
     .	             cin(isp,ien)+cel(isp,ien),cin(isp,ien),cel(isp,ien)
 	    enddo
 	  enddo
	endif
 	write(fic_transout,*)
2003	format(/,4x,'Energy',2x,'bin width',9x,' cross sect',a5,/
     .		4x,' (eV) ',2x,'  (eV)   ',3x,'total',3x,'inelastic',
     .		4x,'elastic')
2004	format(1x,0pf9.2,f9.3,3(1pe11.3))
1080    format('Elastic cross-sections (cel(specie,energy)):',/,
     .		3x,'Energy',6x,5(4x,a4,4x))
1090    format(1x)
1100    format(1pe10.2,2x,5(1pe12.3))
1110	format(/,t20,'Inelastic cross-sections of ',a5,
     .	       /,t20,'---------------------------------',
     .	       /,'    Energy ','   total   ',5a10)
1120	format( /,'    Energy ',10a11)
1130	format(1pe10.2,10e10.3)
1140	format(/' Threshold in eV :')
1150	format(3(1a11,1f6.2,' | '))
        if (iprt(5) .ne. 0) then
	  write(fic_transout,1080) (specie(isp),isp=1,nspec)
	  write(fic_transout,1090)
	  do ien=1,nen
            write(fic_transout,1100) e(ien),(cel(isp,ien),isp=1,nspec)
 	  enddo
	  do isp=1,nspec
	    write(fic_transout,1110) specie(isp),(title(i,isp),i=1,5)
	    do ien=1,nen
	      write(fic_transout,1130) e(ien),cin(isp,ien),
     .			(cinex(isp,iexc,ien),iexc=1,5)
 	    enddo
 	    if(jsg(isp).gt.5)then
	      write(fic_transout,1120) (title(i,isp),i= 6,jsg(isp))
	      do ien=1,nen
	        write(fic_transout,1130) e(ien),
     .		  (cinex(isp,iexc,ien),iexc=6,jsg(isp))
 	      enddo
 	    endif
   	    write(fic_transout,1140)
   	    write(fic_transout,1150)
     .		(title(i,isp),ethres(isp,i,1),i=1,jsg(isp))
 	  enddo
        end if
 	if (iprt(1).eq.1)then
 	  write(fic_transout,*)
          do 130 iang=1,nango2
 	    write(fic_transout,*)'Downward flux. Angle =',pitchang(iang)
            write (fic_transout,5190) (qxdown(ien,-iang),ien=1,nen)
130 	  continue
          do 140 iang=1,nango2
 	    write(fic_transout,*)
     .		'Upward flux.   Angle =',pitchang(nango2+iang)
            write (fic_transout,5190) (qxup(ien,iang),ien=1,nen)
140 	  continue
	endif
5000 	format (7(1x,1pe10.2))
5020 	format (' Electron dens. (denelc)[cm-3]:',/,5(1x,1pe10.2))
5030 	format (' Min. and max. alt. (zbot,ztop)       :',2f10.2)
5040 	format (' Hour angle (hrloc)                       :',f10.2,/,
     .          ' Latitude (glat)                          :',f10.2,/,
     .          ' Julian day (day)                         :',f10.2)
5050 	format (' Desired exospheric temp. (tempexo)       :',f10.2,/,
     .          ' Solar 10.7-cm line output- ave.(f107a)  :',f10.2,/,
     .          ' Solar 10.7-cm line output- today (f107)  :',f10.2,/,
     .	        ' Albedo at lowest altitude (albedo)       :',f10.2)
5070 	format (' Year                                     :',1f10.2,/,
     .		' Number of altitude grid points (nalt)    :',i10)
5080 	format (' Number of energy grid points (nen)       :',i10,/,
     .          ' Number of angle grid points (nang)       :',i10,/,
     .          ' Experiment id. number (knm)              :',i10)
5090 	format (' Working altitude grid points (alt):',/,5(1x,1f10.2))
5100 	format (/,' Energie grid points (e)  :',/,5(1x,1f10.2))
5105 	format (/,' dE grid points (e)  :',/,5(1x,1f10.2))
5110 	format (/,' Angle grid points (pitchang):',/,4(1f16.4,2x))
5112 	format (/,' Corresponding cosines (cosang):',/,4(1f16.13,2x))
5114 	format (/,' Corresponding weight (weitang):',/,4(1f16.13,2x))
5190 	format (' Incident flux vs energy :',/,5(1x,1pe10.2))
!
      return
!
993	print*,' Diff.cross-sect. file is in error. Status=',iost
	stop 'error'		! call abort
992	print*,' Cross-section file is in error. Status=',iost
	stop 'error'		! call abort
997	print*,' Option file is in error. Status=',iost
	stop 'error'		! call abort
      end
