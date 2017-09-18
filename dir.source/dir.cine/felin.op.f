        subroutine felin(knm,nspec,hrloc,day,year,UT,
     &          tempexo,f107,ap,glat,glong,botE,centE,
     &          ddeng,nalt,altkm,tneutre,densneut,colden,chi,chideg,
     &		kiappel,phdisso2,pfluxsr,Po1sdisso2)

c	Computation of the primary photoelectrons, i.e. the electrons
c 	created by the solar photon flux.(jl,1993).
c
c	The photon input flux may be :
c 	-----------------------------
c 		. computed by interpolating on f10.7 between the 2
c		  fluxes for sunspot max and min. (Torr and Torr)
c		. given by the Tobiska's code (ref. in the subroutine)
c 		  for a given date and solar activity
c 	In any case, this solar flux is computed on the enlarged
c 	wavelength grid : 39 Tobiska points instead of the initial
c 	Hintereger, Torr and Torr' 37. The 39 energies are the same
c 	as the 37, plus two higher energies.
c 	The choice of the solar flux is made with the switch iflux in
c 	the file DATFEL
c 	iflux = 0 --> photon flux computed from f107.
c             = 1 --> photon flux solar min with requiered f107.
c             = 2 --> photon flux solar max with requiered f107.
c             = 3 --> photon flux from Tobiska for given date
c 		      (68172 --> 88366)
c             = 4 --> photon flux from EUVAC
c
c 	Neutrals are computed after N2,O2,O,H,He
c
c	Branching ratio
c 	---------------
c 	Rees, Physics and chemistry of the upper atmosphere, cambridge,
c 	1989.
c
c 	Cross sections
c 	--------------
c 	from Torr and Torr (1985), or Fennelly and Torr (1992).
c 	The choice of the set is made with the switch iseff in the file
c 	DATFEL.
c 	1 = Torr et Torr, 1985, 2 = Fennelly and Torr, 1992
c
c 	There are 2 computation modes :
c 	-------------------------------
c	  1) On the 39 energy points
c	  2) On a 1 eV width grid, made out of dividing the continuum
c 	     intervals.
c 	The choice of the set is made with the switch imod in the file
c 	DATFEL.
c 	Important remark : the bad discretization of the wavelength
c 	grid from Hintereger and Torr leads to a bad discretization of
c 	the energies at which the photoelectrons are created.
c 	For example, a wavelength 2.327 nm corresponds to an energy
c 	533 eV. The ionization threshold is about 15 eV. All created
c 	photoelectrons will have the energy 533-15 = 518 eV, and will
c 	fall down in the same energy box (the one that includes 518 eV).
c 	The next wavelength 3.750 nm, corresponds to 331 eV. All next
c 	created photoelectron will have energy 331-15 = 316 eV, and
c 	fall down in this energy box. Between 518 and 316 eV, no
c 	electron can be created ! One could compute the energy distribu-
c 	tion on the 39 Tobiska /Hintereger ' steps, but the cross
c 	sections vary greatly between one step to the other, so that the
c 	best way is to use the computation mode that divides the
c 	continuum in short intervals, when the energy spectrum is
c 	requiered.
c
c 	The Chapman function
c 	--------------------
c 	Given by Stan Solomon, who got it from Nagy. It comes from:
c 	"Numerical evaluation of Chapman's grazing incidence integral
c 	ch(X,x)", Smithe and Smith, JGR, vol 77, 1972, p3592,3597
c
c               --> si kiappel = 1, ce sont les donnees
c                   contenue dans ELEC et NEUTRAL pour un run unique
c                   du programme.
c               --> si kiappel = 2, ce sont les donnees necessaires
c                   au transport fluide
        include 'comm.f'
        implicit none
        include 'TRANSPORT.INC'

        real, intent(in) :: knm, tempexo, glong, altkm(nbralt),
     &    tneutre(nbralt),densneut(8,nbralt),colden(8,nbralt)
        integer, intent(in) :: nspec
	
      	common /bloc/ threshold,nbseff,eVseff,seffion,sefftot,pfluxmin,
     .  	      pfluxmax,wave,eV,wavemin,wavemax,eVmin,eVmax,
     .               nwave,ns,nns,f107min,f107max,iseff,wnmseff,
     .               lambdasr,sigsro2,Isr,lineflux,sigabso2,qyield,
     .               Isr2
        integer nalt,ns,nns,nwave,iseff
      real f107max,f107min,wavemin(39),wavemax(39),eVmin(39),eVmax(39)
      	real f107(3),day,ap(7),hrloc,UT,glat,chi,chideg
c 	wwt (specie, excit, energy) = branching ratio
       	real threshold(7),wwt(7,6,39)
        integer number(7),nbseff
c     	Ecent = energy boxes (eV).
c     	Ebot  = lower boundary of energy boxes (eV).
c     	engdd = energy boxes width(eV).
c 	Remark: If the energy grid on which the production is computed
c 	has its max below the ionization max (250eV), one can create
c 	electrons at higher energies than the max of the grid.
c
       	real eV(39),centE(nbren),ddeng(nbren),botE(nbren)
       	real Ecent(nbren),engdd(nbren),Ebot(nbren)
       	real proelec(nbralt),
     .		  prodion(nbralt,nbrsp*2),proneut(nbralt,nbrsp),
     . 	  	  prophel(nbralt,nbren),produc(nbralt,nbrsp,nbren)
        real prodcont(nbralt),prodraie(nbralt)
      	real seffion(2000,7),sefftot(2000,5),eVseff(2000),wnmseff(2000)
      	real sigi(39,7),sigt(39,5)         ! Seff sur les en(Torr)
      	real pflux(39),wave(39)
        real pfluxmin(39),pfluxmax(39)
        real altregE(nbralt),denselregE(nbralt)
        real xchap(nbralt,nbrsp)
        integer iprt(12),idess(9),kiappel
        real trav(nbren),sum,year
        real sflux(nbralt,nbren),phdisso2(500),sigsro2(8),lambdasr(9)
        real lineflux(7),sigabso2(7),qyield(7)
        real Isr(8),pfluxsr(8,nbralt),chapsp(8)
        real Isr2(8)
        real Po1sdisso2(500)
        integer iz,iwave,ichapman,iflux,ilambda,imod,nrege
      integer ifeldat,ifelprt,ifeltrans
        real maxchi,nan
        

        print *,'felin.f: nen=',nen

c ---	Open up files to be used in this program.
      open (newunit=ifeldat,file='dir.data/dir.linux/dir.cine/DATFEL',
     &          status='old')

      open (newunit=ifelprt,file='dir.data/dir.linux/dir.cine/FELOUT',
     &          status='unknown')
        rewind(ifelprt)

        open (newunit=ifeltrans,
     .		file='dir.data/dir.linux/dir.cine/FELTRANS',	access='stream',
     &          status='old')
        rewind(ifeltrans)

      	write(ifelprt,*) 'felin1.f'

c --- 	Resets all output arrays.

       	call zeroit(proelec,nbralt)
      	call zeroit(prophel,nbralt*nbren)
        call zeroit(produc,nbralt*nbrsp*nbren)
     	call zeroit(prodion,nbralt*(nbrsp*2))
        call zeroit(proneut,nbralt*nbrsp)
 	    call zeroit(prodcont,nbralt)
 	    call zeroit(prodraie,nbralt)
c
c --- 	Set the input parameters.
c
        call setpar(kiappel,imod,nspec,idess,hrloc,UT,day,nan,tempexo,
     .	   f107,ap,glat,glong,nalt,year,
     .	   altkm,centE,botE,ddeng,Ecent,Ebot,engdd,iprt,pflux,knm,
     .	   tneutre,densneut,colden,iflux,wwt,number,xchap,chi,chideg,
     .	   sigi,sigt,ichapman)
c
c --- 	Calculates primary photoelectron production.
        call depo(kiappel,nalt,ns,nns,nwave,chi,densneut,colden,altkm,
     &		  threshold,eV,Ecent,Ebot,engdd,produc,
     &		  proelec,prodion,proneut,prophel,prodcont,prodraie,
     &    	  sigi,sigt,pflux,tneutre,iprt,number,wwt,xchap,
     &		  eVmin,eVmax,imod,nspec,sflux)

!       Compute production of O(1D) from photodissociation of O2
        maxchi=1.7453

          do iz=1,nalt
                phdisso2(iz)=0.

                if(chi .lt. maxchi) then
                call fchap(day,UT,altkm(iz),
     &           glat,glong,f107,ap,chi,chapsp)

                do ilambda=1,8
                        pfluxsr(ilambda,iz)=Isr2(ilambda)*
     &                   exp(-sigsro2(ilambda)*chapsp(2)*colden(2,iz))

                        phdisso2(iz)=phdisso2(iz)+sigsro2(ilambda)
     &                   *pfluxsr(ilambda,iz)*densneut(2,iz)
                enddo
                endif
         enddo

!       Compute O(1S) production
         do iz=1,nalt
                Po1sdisso2(iz)=0.

                if(chi .lt. maxchi) then
!               Compute production of O(1S) from continuum dissociation of O2
                Po1sdisso2(iz)=densneut(2,iz)*(
     &                  sflux(iz,34)*12.81e-18*0.029+
     &                  sflux(iz,36)*21.108e-18*0.0575+
     &                  sflux(iz,39)*1.346e-18*0.0865)

!               Compute production of O(1S) from spectral line dissociation of O2
                call fchap(day,UT,altkm(iz),
     &           glat,glong,f107,ap,chi,chapsp)

                do ilambda=1,7
                        pfluxsr(ilambda,iz)=lineflux(ilambda)*
     &                   exp(-sigabso2(ilambda)*chapsp(2)*colden(2,iz))

                        Po1sdisso2(iz)=Po1sdisso2(iz)+sigabso2(ilambda)*
     &                    pfluxsr(ilambda,iz)*densneut(2,iz)*
     &                    qyield(ilambda)
                enddo

!               Compute production of O(1S) from Ly alpha dissociation of O2
                pfluxsr(1,iz)=251.e9*
     &           exp(-1.e-20*chapsp(2)*colden(2,iz))

                Po1sdisso2(iz)=Po1sdisso2(iz)+1.e-20*pfluxsr(1,iz)*
     &           densneut(2,iz)*0.012
                endif
!                print*,1.e-20*pfluxsr(1,iz)*densneut(2,iz)*0.012
!                print*,altkm(iz),pfluxsr(1,iz),densneut(2,iz)
!                print*,1.e-20*0.012,pfluxsr(1,iz)*densneut(2,iz)
         enddo
!       -MZ

c
c --- 	Calculates electron density.
        call densout1(nalt,altkm,proelec,denselregE,nregE,altregE,
     .		    densneut)
c
c ---	Writes the different productions.
        write(stdout,*),'felin.f: call prodprt nen=',nen
        call prodprt(ns,nalt,altkm,Ecent,engdd,produc,proelec,
     .  	      prodion,proneut,prophel,iprt,ichapman)

      end subroutine felin
c
c---------------------- function chap ---------------------------------
c
      real function chapsmith (CHI, ZCM, T, AM)
        implicit none

        real, intent(in) :: chi, zcm, t, am 
c 	The Chapman function was given by Stan Solomon, who got it from
c  	Nagy. Comes from:
c 	"Numerical evaluation of Chapman's grazing incidence integral
c 	ch(X,x)", Smithe and Smith, JG, vol 77, 1972, p3592,3597
c
       common /const/ pi,re,recm,bolt,gzero,amu
       real a,b,c,d
       real pi,re,recm,bolt,gzero,amu
       real gr,hf,hg,hn,sqhf
       real,external :: sperfc

        GR=GZERO*(RECM/(RECM+ZCM))**2
        HN=1.38E-16*T/(AM*1.662E-24*GR)
        HG=(RECM+ZCM)/HN
        HF=0.5*HG*(COS(CHI)**2)
        SQHF=SQRT(HF)
        if(chi.le.pi/2.)then
      	  chapsmith=SQRT(0.5*PI*HG)*SPERFC(SQHF)
        else
          a = sqrt(0.5*pi*hg)
          c = sqrt(sin(chi))*exp(hg*(1.-sin(chi)))
          d = 0.5*sperfc(sqhf)
          b = c - d
          chapsmith = 2.*a * b
        endif

      END function chapsmith
c
c------------------- function chapgreen -------------------------------
c
        real function chapgreen (CHI, ZCM, Tneutre, atomas)
c
        implicit none
        real, intent(in) :: chi,zcm,Tneutre,atomas
c
c       "Molecular absorption in planetary atmosphere"
c       Green, Lindenmeyer and Griggs, JGR, vol 69, 1964, p493-504
c
        real pi,re,recm,bolt,gzero,amu
        common /const/ pi,re,recm,bolt,gzero,amu
c
        real gr,hn,X,alpha,c
c
c       gr = gravite a l'atitude zcm [cm.s-2]
        GR=GZERO*(RECM/(RECM+ZCM))**2
        c = pi/2.
c
c       hn = kT/mg = hauteur d'echelle [cm]
        HN=bolt*Tneutre/(atomas*amu*GR)
c       X = hauteur au centre de la Terre/hauteur d'echelle
        X=(RECM+ZCM)/HN

        if(chi <= pi/2.) then
          alpha = 1./c**4 - 0.115/c**2 - 0.5/(c**2*log(sqrt(c*X)))
          chapgreen = exp(0.5*chi**2/1.-0.115*chi**2-alpha*chi**4)
        endif

        end function chapgreen
c
c----------------------------------------------------------------------
c
        subroutine branchratio(nwave,eV,threshold,num,wwt)

        implicit none
        include 'TRANSPORT.INC'

        integer, intent(in) :: nwave
        real, intent(in) :: eV(nwave),threshold(7)
        integer,intent(out) :: num(7)
        real,intent(out)::wwt(7,6,39)
 	
c
c 	Calcule les branching ratio sur les energies des ondes
c 	d'entree.
c
c 	INPUTS
c
        real eVbrO(10),brO(5,10),eVbrN2(13),brN2(5,13),
     .  	eVbrO2(20),brO2(3,20),eVbrO2dis(14),brO2dis(6,14),
     .		eVbrN2dis(2),brN2dis(1,2),eVbrH(2),brH(1,2),
     .		eVbrHe(2),brHe(1,2)
        integer neVbrO,nstO,neVbrN2,nstN2,
     .  	neVbrO2,nstO2,neVbrO2dis,nstO2dis,neVbrN2dis,nstN2dis,
     .  	neVbrH,nstH,neVbrHe,nstHe
        common /branching/
     .		eVbrN2,    brN2,    neVbrN2,    nstN2,
     .		eVbrO2,    brO2,    neVbrO2,    nstO2,
     .		eVbrO,     brO,     neVbrO,     nstO,
     .		eVbrO2dis, brO2dis, neVbrO2dis, nstO2dis,
     .		eVbrN2dis, brN2dis, neVbrN2dis, nstN2dis,
     .		eVbrH,     brH,     neVbrH,     nstH,
     .		eVbrHe,    brHe,    neVbrHe,    nstHe

        integer iprt(12)

c 	OUTPUTS
c 	wwt (iontype,state,energy) is the normalized branching ratio
c 	array.

c 	INTERNAL
 	integer ist,iwave,iontype
 	real tabin(500),tabout(500),xnorm
c
c ---- 	N2 -->N2+
c
 	iontype = 1
 	do ist = 1,nstN2
 	  do iwave = 1,neVbrN2
 	    tabin(iwave) = brN2(ist,iwave)
 	  enddo
 	  call intlin(neVbrN2,eVbrN2,tabin,nwave,eV,tabout)
 	  do iwave = 1,nwave
 	    wwt(iontype,ist,iwave) = max(0.,tabout(iwave))
 	    if(ev(iwave).lt.threshold(iontype))wwt(iontype,ist,iwave)=0.
 	  enddo
 	enddo
c 	Normalization
 	do iwave = 1,nwave
 	  if(ev(iwave).ge.threshold(iontype))then
 	    xnorm = 0.
 	    do ist = 1,nstN2
 	      xnorm = xnorm + wwt(iontype,ist,iwave)
 	    enddo
 	    do ist = 1,nstN2
 	      wwt(iontype,ist,iwave) = wwt(iontype,ist,iwave)/xnorm
 	    enddo
 	  endif
 	enddo
 	num(iontype) = nstN2
c
c ---- 	O2 -->O2+
c
 	iontype = 2
 	do ist = 1,nstO2
 	  do iwave = 1,neVbrO2
 	    tabin(iwave) = brO2(ist,iwave)
 	  enddo
 	  call intlin(neVbrO2,eVbrO2,tabin,nwave,eV,tabout)
 	  do iwave = 1,nwave
 	    wwt(iontype,ist,iwave) = max(0.,tabout(iwave))
 	    if(ev(iwave).lt.threshold(iontype))wwt(iontype,ist,iwave)=0.
 	  enddo
 	enddo
c 	Normalization
 	do iwave = 1,nwave
 	  if(ev(iwave).ge.threshold(iontype))then
 	    xnorm = 0.
 	    do ist = 1,nstO2
 	      xnorm = xnorm + wwt(iontype,ist,iwave)
 	    enddo
 	    do ist = 1,nstO2
 	      wwt(iontype,ist,iwave) = wwt(iontype,ist,iwave)/xnorm
 	    enddo
 	  endif
 	enddo
 	num(iontype) = nstO2
c
c ---- 	O -->O+
c
 	iontype = 3
 	do ist = 1,nstO
 	  do iwave = 1,neVbrO
 	    tabin(iwave) = brO(ist,iwave)
 	  enddo
 	  call intlin(neVbrO,eVbrO,tabin,nwave,eV,tabout)
 	  do iwave = 1,nwave
 	    wwt(iontype,ist,iwave) = max(0.,tabout(iwave))
 	    if(ev(iwave).lt.threshold(iontype))wwt(iontype,ist,iwave)=0.
 	  enddo
 	enddo
c 	Normalization
 	do iwave = 1,nwave
 	  if(ev(iwave).ge.threshold(iontype))then
 	    xnorm = 0.
 	    do ist = 1,nstO
 	      xnorm = xnorm + wwt(iontype,ist,iwave)
 	    enddo
 	    do ist = 1,nstO
 	      wwt(iontype,ist,iwave) = wwt(iontype,ist,iwave)/xnorm
 	    enddo
 	  endif
 	enddo
 	num(iontype) = nstO
c
c ---- 	N2 -->N+ + N
c
 	iontype = 4
 	do ist = 1,nstN2dis
 	  do iwave = 1,neVbrN2dis
 	    tabin(iwave) = brN2dis(ist,iwave)
 	  enddo
 	  call intlin(neVbrN2dis,eVbrN2dis,tabin,nwave,eV,tabout)
 	  do iwave = 1,nwave
 	    wwt(iontype,ist,iwave) = max(0.,tabout(iwave))
 	    if(ev(iwave).lt.threshold(iontype))wwt(iontype,ist,iwave)=0.
 	  enddo
 	enddo
c 	Normalization
 	do iwave = 1,nwave
 	  if(ev(iwave).ge.threshold(iontype))then
 	    xnorm = 0.
 	    do ist = 1,nstN2dis
 	      xnorm = xnorm + wwt(iontype,ist,iwave)
 	    enddo
 	    do ist = 1,nstN2dis
 	      wwt(iontype,ist,iwave) = wwt(iontype,ist,iwave)/xnorm
 	    enddo
 	  endif
 	enddo
 	num(iontype) = nstN2dis
c
c ---- 	O2 -->O+ + O
c
 	iontype = 5
 	do ist = 1,nstO2dis
 	  do iwave = 1,neVbrO2dis
 	    tabin(iwave) = brO2dis(ist,iwave)
 	  enddo
 	  call intlin(neVbrO2dis,eVbrO2dis,tabin,nwave,eV,tabout)
 	  do iwave = 1,nwave
 	    wwt(iontype,ist,iwave) = max(0.,tabout(iwave))
 	    if(ev(iwave).lt.threshold(iontype))wwt(iontype,ist,iwave)=0.
 	  enddo
 	enddo
c 	Normalization
 	do iwave = 1,nwave
 	  if(ev(iwave).ge.threshold(iontype))then
 	    xnorm = 0.
 	    do ist = 1,nstO2dis
 	      xnorm = xnorm + wwt(iontype,ist,iwave)
 	    enddo
 	    do ist = 1,nstO2dis
 	      wwt(iontype,ist,iwave) = wwt(iontype,ist,iwave)/xnorm
 	    enddo
 	  endif
 	enddo
 	num(iontype) = nstO2dis
c
c ---- 	H -->H+
c
 	iontype = 6
 	do ist = 1,nstH
 	  do iwave = 1,neVbrH
 	    tabin(iwave) = brH(ist,iwave)
 	  enddo
 	  call intlin(neVbrH,eVbrH,tabin,nwave,eV,tabout)
 	  do iwave = 1,nwave
 	    wwt(iontype,ist,iwave) = max(0.,tabout(iwave))
 	    if(ev(iwave).lt.threshold(iontype))wwt(iontype,ist,iwave)=0.
 	  enddo
 	enddo
c 	Normalization
 	do iwave = 1,nwave
 	  if(ev(iwave).ge.threshold(iontype))then
 	    xnorm = 0.
 	    do ist = 1,nstH
 	      xnorm = xnorm + wwt(iontype,ist,iwave)
 	    enddo
 	    do ist = 1,nstH
 	      wwt(iontype,ist,iwave) = wwt(iontype,ist,iwave)/xnorm
 	    enddo
 	  endif
 	enddo
 	num(iontype) = nstH
c
c ---- 	He -->He+
c
 	iontype = 7
 	do ist = 1,nstHe
 	  do iwave = 1,neVbrHe
 	    tabin(iwave) = brHe(ist,iwave)
 	  enddo
 	  call intlin(neVbrHe,eVbrHe,tabin,nwave,eV,tabout)
 	  do iwave = 1,nwave
 	    wwt(iontype,ist,iwave) = max(0.,tabout(iwave))
 	    if(ev(iwave).lt.threshold(iontype))wwt(iontype,ist,iwave)=0.
 	  enddo
 	enddo
c 	Normalization
 	do iwave = 1,nwave
 	  if(ev(iwave).ge.threshold(iontype))then
 	    xnorm = 0.
 	    do ist = 1,nstHe
 	      xnorm = xnorm + wwt(iontype,ist,iwave)
 	    enddo
 	    do ist = 1,nstHe
 	      wwt(iontype,ist,iwave) = wwt(iontype,ist,iwave)/xnorm
 	    enddo
 	  endif
 	enddo
 	num(iontype) = nstHe
c
 	end subroutine branchratio
c
c---------------------------- densout1 -------------------------------
c
      subroutine densout1(nalt,altkm,proelec,denselregE,nregE,altregE,
     &		 densneut)

      include 'TRANSPORT.INC'

      real, intent(in) :: altkm(nbralt),densneut(8,nbralt)
      integer,intent(out) :: nregE

      real proelec(nbralt),denselregE(nbralt)

      real altregE(nbralt)
      integer ialt,nalt
c
c     	calcul de densites
 	nregE = 0
        do ialt=1,nalt
          if (altkm(ialt).lt.85.)then
 	    nregE = nregE+1
 	    altregE(nregE) = altkm(ialt)
	    alphaeff=7.30e+04*exp(-altkm(ialt)/3.3)
 	    denselregE(nregE)=sqrt(proelec(ialt)/alphaeff)
	  elseif (altkm(ialt).ge.85. .and. altkm(ialt).le.200.)then
 	    nregE = nregE+1
 	    altregE(nregE) = altkm(ialt)
	    alphaeff=2.50e-06*exp(-altkm(ialt)/51.2)
 	    denselregE(nregE)=sqrt(proelec(ialt)/alphaeff)
	  endif
 	enddo
 	write(ifelprt,7050)
 	do ialt = 1,nregE
 	  write(ifelprt,7051)altregE(ialt),proelec(ialt),
     .		denselregE(ialt)
 	enddo
7050 	format(/,'Electron density in the E region',/,' altitude ',
     &    ' e- prod.  Computed [Ne]')
7051 	format(1f10.2,2(1pe12.3))
c
 	return
 	end
c
c------------------------- depo ----------------------------------
c
      subroutine depo(kiappel,nalt,ns,nns,nwave,chi,densneut,colden,
     .        altkm,threshold,eV,Ecent,Ebot,engdd,produc,
     .	      proelec,prodion,proneut,prophel,prodcont,prodraie,
     .        sigi,sigt,pflux,tneutre,iprt,num,wwt,xchap,
     .	      eVmin,eVmax,imod,nspec,sflux)
        include 'comm.f'
        implicit none

        include 'TRANSPORT.INC'

        integer, intent(in) :: nspec,imod
	      real, intent(in) :: altkm(nbralt),tneutre(nbralt),
     &    densneut(8,nbralt),colden(8,nbralt)

c
c 	Computes the different productions :
c     	produc  = electron prod. at alt no iz,neutral specie j,
c 	 	  box iener  [cm-3.s-1.ev-1]
c     		produc(iz,neutspe=1,iener)--->N2
c     		produc(iz,neutspe=2,iener)--->O2
c     		produc(iz,neutspe=3,iener)--->O
c     		produc(iz,neutspe=4,iener)--->H
c     		produc(iz,neutspe=5,iener)--->He
c     	proneut = Ion prod. at alt no iz,due to neutral specie j.
c     	        [cm-3.s-1]
c     		proneut(iz,neutspe=1)--->N2
c     		proneut(iz,neutspe=2)--->O2
c     		proneut(iz,neutspe=3)--->O
c     		proneut(iz,neutspe=4)--->H
c     		proneut(iz,neutspe=5)--->He
c 	Prodstion is only computed because it is a "natural" output of
c 	the code. Should it be used sometimes, it should be included
c 	in the subroutine calling list
c     	prodstion = Ion prod. at alt no iz,ion specie j, state ist.
c     	        [cm-3.s-1]
c     		prodstion(iz,iontype=1,ist)--->N2+
c     		prodstion(iz,iontype=2,ist)--->O2+
c     		prodstion(iz,iontype=3,ist)--->O+ (from O --> O+)
c     		prodstion(iz,iontype=4,ist)--->N+
c     		prodstion(iz,iontype=5,ist)--->O+ (from O2 --> O + O+)
c     		prodstion(iz,iontype=6,ist)--->H+
c     		prodstion(iz,iontype=7,ist)--->He+
c     	prodion = Ion prod. at alt no iz,ion specie j.
c     	        [cm-3.s-1]
c     		prodion(iz,ionspe=1)--->N2+
c     		prodion(iz,ionspe=2)--->O2+
c     		prodion(iz,ionspe=3)--->O+
c     		prodion(iz,ionspe=4)--->N+
c     		prodion(iz,ionspe=5)--->H+
c     		prodion(iz,ionspe=6)--->He+
c     	prophel = electron prod. at alt no iz,box iener, all species.
c     	        [cm-3.s-1.ev-1]
c     	proelec = electron prod. at alt no iz, all species mixed.
c     	        [cm-3.s-1]
c     	prodraie= electron prod. at alt no iz due to discrete
c     	        lines [cm-3.s-1]
c     	prodcont = electron prod. at alt no iz due to the
c     	        solar continuum [cm-3.s-1]

       integer,intent(in) :: kiappel,nalt,ns,nns,nwave
       real,intent(out) :: sflux(nbralt,nbren)        !I(z,lambda)

       	real chi
       	real threshold(7),proelec(nbralt),prodion(nbralt,nbrsp*2),
     . 		proneut(nbralt,nbrsp),prophel(nbralt,nbren)
        real prodcont(nbralt),prodraie(nbralt)
       	real produc(nbralt,nbrsp,nbren),pflux(nwave),
     .            sigi(nwave,7),eV(nwave),eVmin(nwave),eVmax(nwave)
        real Ecent(nbren),engdd(nbren),Ebot(nbren)
c 	wwt = branching ratio sur les energies de Torr
       	real wwt(7,6,nwave)
c 	wt = branching ratio sur les energies des petites boites
       	real wt(6),dev,flux,exa
        integer num(7),iontype,nexcit
       	integer iprt(12),iwave,iz,isp,ist,ninter,ne,ionspe,i
      real xchap(nbralt,nbrsp),sigt(nwave,5),dele,deltae,amu,gzero,bolt,
     &  recm,re,pi,sefint,energy
     


 	common /const/ pi,re,recm,bolt,gzero,amu
c
c 	DATFEL = Input file to be used by felin.f (formatted).
c 	FELTRANS = Result file to be used by trans.f (unformatted).
c 	FELOUT   = Result file to be read by user (formatted).
 	real sigtot(5),sigion(7)
c
!	write(6,*)'Computes the primary photoelectron production.'
!      write(stdout,*),'felin.f:depo   nen,nwave,nalt,ns,nns',
!     & nen,nwave,nalt,ns,nns
     
      !write(stdout,*),shape(sflux),sizeof(sflux)

        do iwave=1,nwave
!     .		write(6,*)'Wavelength ',iwave,'/',nwave
 	  deV = eVmax(iwave)-eVmin(iwave)
 	  if(deV .le. 1. .or. imod==1)then
c 	    On est sur une ligne discrete .ou. on ne coupe pas les
c 	    intervalles en petite boites.
 	    flux = pflux(iwave)
      	do iz = 1,nalt
 	      exa = 0.
 	      do isp = 1,ns
                exa = exa + xchap(iz,isp)*sigt(iwave,isp)*colden(isp,iz)
 	      enddo
              if (exa.gt.35.0) cycle
 	      do iontype = 1,nns
!               Computes the productions :
 	        sigion(iontype) = sigi(iwave,iontype)
 	        nexcit = num(iontype)
 	        dele = eV(iwave)-threshold(iontype)
     	 	do ist = 1,num(iontype)
     	  	  wt(ist) = wwt(iontype,ist,iwave)
     		enddo
!       write(stdout,*),'felin.f:depo iwave,iontype,iz,nen',
!     & iwave,iontype,iz,nen
       	        call eval (dev,dele,flux,iontype,sigion(iontype),
     &		           nexcit,wt,Ebot,engdd,densneut,exa,
     &	                   iz,proelec,prodion,produc,proneut,prophel,
     &			   prodcont,prodraie)
 	      enddo		! boucle sur les especes
              sflux(iz,iwave) = pflux(iwave)*exp(-exa)
!              write(stdout,*),exa,sflux(iz,iwave)
       end do !iz
 	  else
c	    On est sur un continuum. On divise la boite en plusieurs
c 	    petits intervalles de largeur a peu pres egale a 1 eV
c 	    Leur nombre est ninter :
 	    ninter = ifix(eVmax(iwave)-eVmin(iwave)) - 1
c 	    Leur largeur est delta :
  	    deltaE = (eVmax(iwave)-eVmin(iwave))/float(ninter)
c 	    Le flux dans chaque petite boite est :
 	    flux = pflux(iwave)/float(ninter)
c
 	    do i = 1,ninter
              write(6,*)'Wavelength ',iwave,'/',nwave,
     .              ' interval ',i,'/',ninter,'  [A'
c
c 	      On place les energies au milieu de chaque petite boite
 	      energy = eVmin(iwave) + deltaE/2. + float(i-1)*deltaE
c
 	      do ionspe = 1,ns
c 	        On calcule la section efficace totale d'absorption a
c	        cette energie, pour le coeff d'absorption exa.
 	        sigtot(ionspe) = sefint(energy,ionspe,1)
 	      enddo
 	      do iontype = 1,nns
c 	        On calcule la section efficace d'ionisation (diss. ou
c 		non) a cette energie.
 	        sigion(iontype) = sefint(energy,iontype,2)
c 	        Le rapport de branchement a cette energie est (en
c 	        gros) celui de l'energie iwave
 	 	do ist = 1,num(iontype)
 	  	  wt(ist) = wwt(iontype,ist,iwave)
 		enddo
 	      enddo
 	      
 	   ! write(stdout,*),'felin.f:depo "else" eval  nen=',nen
      	      do iz = 1,nalt
 	            exa = 0.
 	            do isp = 1,ns
                  exa = exa + xchap(iz,isp)*sigtot(isp)*colden(isp,iz)
 	            enddo
                if (exa.gt.35.0) cycle
                
 	           do iontype = 1,nns
c                 Computes the productions :
 	            nexcit = num(iontype)
 	            dele = energy-threshold(iontype)
       	          call eval (dev,dele,flux,iontype,sigion(iontype),
     .		           nexcit,wt,Ebot,engdd,densneut,exa,
     .		           iz,proelec,prodion,produc,proneut,prophel,
     .			   prodcont,prodraie)
 	           enddo		! iontype species
                sflux(iz,iwave)=pflux(iwave)*exp(-exa)
             enddo  !iz
 	    enddo		! boucle sur les intervalles
 	  endif
c
      enddo ! iwave

      	end subroutine depo
c
c--------------------------- eval -----------------------------------
c
      subroutine eval (deV,dele,flux,iontype,sigion,nexcit,wt,
     .	      Ebot,engdd,densneut,exa,iz,
     .	      proelec,prodion,produc,proneut,prophel,
     .	      prodcont,prodraie)
c
c 	INPUTS
c 	------
c 	flux   	: Incident photon flux at energy
c 	iontype	: considered ion type (considering the type of reaction)
c 	sigion 	: ionization or dissociative ionization cross section
c 	nexcit 	: # of excitation states
c 	wt 	: weitghing factors (branching ratio)
c 	dele 	: Energy in eV of the incident wavelength minus the
c   		  ionization threshold (eV)
c 	 	  i.e. the energy that is available for the electron
c 	nen, Ecent, engdd, Ebot 	: Ouput energy grid for
c 					  created electrons
c 	densneut	: Neutral density
c 	exa 	: Sum(xchap(iz,isp)*sigt(iwave,isp)*colden(isp,iz))
c 	iz 	: Altitude number
c
c 	OUTPUTS
c 	-------
c 	Production:
c 	Prim. photoel. prod. = Integral(N(specie j).sigma(E).I(z,E).dE)
c                       I = flux at infinity.exp(-exa).wt
c                         = ---------dflx------------.wt
c                         = -------------depr-----------
c 	wt = Normalized weight factors. Determine which ion state is
c 	excited.
        include 'comm.f'
        implicit none
        include 'TRANSPORT.INC'
        
        integer,intent(in) :: nexcit,iontype
        real,intent(in) :: deV,dele,exa,wt(6),sigion,flux,
     &   densneut(8,nbralt)

     	real engdd(nbren),Ebot(nbren)


      	real proelec(nbralt),prodion(nbralt,nbrsp*2),
     .		produc(nbralt,nbrsp,nbren),proneut(nbralt,nbrsp),
     .		prophel(nbralt,nbren),prodstion(nbralt,7,6)
        real prodcont(nbralt),prodraie(nbralt)
c       nbrelec = nbre d'electrons crees dans l'ionisation
        real nbrelec
c
        integer iz,ist,neutspe,ionspe,iener
        real psurde,depr

        !write(stdout,*),'felin.f:eval   nen=',nen

        if(dele.gt.0.)then
          depr = flux*exp(-exa)
          if(iontype==1) then
c           case N2 --> N2+  (neutral 1 -->ion 1)
 	    ionspe=1
 	    neutspe=1
 	    nbrelec = 1.
          elseif(iontype==2) then
c           case O2 --> O2+  (neutral 2 -->ion 2)
 	    ionspe=2
 	    neutspe=2
 	    nbrelec = 1.
          elseif(iontype==3) then
c           case O --> O+  (neutral 3 -->ion 3)
 	    ionspe=3
 	    neutspe=3
 	    nbrelec = 1.
          elseif(iontype==4) then
c           Case N2 -->N+ + N  (neutral 1 -->ion 4)
            ionspe=4
 	    neutspe=1
 	    nbrelec = 1.
          elseif(iontype==5) then
c           Case O2 --> O+ + O  (neutral 2 -->ion 3)
            ionspe=3
 	    neutspe=2
 	    nbrelec = 1.
          elseif(iontype==6) then
c           Case H --> H+  (neutral 4 -->ion 5)
            ionspe=5
 	    neutspe=4
 	    nbrelec = 1.
          elseif(iontype==7) then
c           Case He --> He+  (neutral 5 -->ion 6)
            ionspe=6
 	    neutspe=5
 	    nbrelec = 1.
          endif
          depr=depr*sigion*densneut(neutspe,iz)
 	  do ist = 1,nexcit
            prodstion(iz,iontype,ist)=prodstion(iz,iontype,ist)+
     .					depr*wt(ist)
 	  enddo
          prodion(iz,ionspe)=prodion(iz,ionspe)+ depr
          proneut(iz,neutspe)=proneut(iz,neutspe)+ depr
          proelec(iz)=proelec(iz)+ depr*nbrelec
 	  if(deV.lt.1.e-05)then
 	    prodraie(iz) = prodraie(iz) + depr*nbrelec
 	  else
 	    prodcont(iz) = prodcont(iz) + depr*nbrelec
 	  endif
c
c         search in which box iener to put the created electron.
          call search (Ebot,engdd,dele,iener)
c 	  Makes it in [cm-3.s-1.ev-1].
 	  psurde = depr*nbrelec/engdd(iener)
          produc(iz,neutspe,iener)=produc(iz,neutspe,iener) + psurde
          prophel(iz,iener) = prophel(iz,iener) + psurde
        endif

      	end subroutine eval
c
c ---------------------------- prodprt -----------------------------
c
        subroutine prodprt(ns,nalt,altkm,Ecent,engdd,produc,proelec,
     .  	      prodion,proneut,prophel,iprt,ichapman)

        include 'comm.f'
        implicit none


        include 'TRANSPORT.INC'


        integer,intent(in) :: ns,nalt,iprt(12)
       	real,intent(in) :: altkm(nbralt),Ecent(nbren),engdd(nbren),
     & prophel(nbralt,nbren)
       	
       	real proneut(nbralt,nbrsp),proelec(nbralt),
     .		prodion(nbralt,nbrsp*2),
     .		produc(nbralt,nbrsp,nbren)
        real zwork(nbralt),qphelev,qphelerg
        
        integer :: ien1,ien2,neutspe,naltO6,ialt1,ialt2,ii,ialt,i,isp,
     &  ichapman,ien,iz,m,nnz,nalt06,j

c
1000  	format(/,'Electron and ion production(/cm3.s)',/)
1010  	format(5x,'alt',7x,'Ne',8x,'Nn2+',7x,'No2+',7x,'No+',6x,'Nn+')
1015  	format(/,5x,'alt',7x,'H+',8x,' He+')
1020  	format(1h ,6(1pe10.2))
1040  	format( /,'Production (/cm3.s),at 5 alt., in',i5,
     .    ' energy boxes,due to N2,O2,O,H,He(produc(iz,neutspe,ien))')
1050    format(1x,/,'alt=',f6.1)
1060  	format(2(1x,i2,3(1pe10.2)))
1070  	format (/,5x,'Produced Photoel. (/cm3.s.eV) at each altitude',
     .    ' (prophel(ialt,ien))',
     .    / 12x, 'Altitude  ',
     .    / , 5x, ' E ', 1x, 7f10.2)
1080  	format (1x, f8.2, 8(1pe10.2))
1100  	format(1h ,7(1pe10.2))
1110  	format('Electron and ionized neutral production (cm-3.s-1)',
     .    ' (proneut(iz,neutspe))',/,
     .    5x,'alt',7x,'Ne',9x,'Nn2',7x,'No2 ',6x,'No1 ',
     .    7x,'H ',7x,'He')

        write(stdout,*),'felin.f:prodprt  nen=',nen


        if(iprt(8)==1)then
      	  write(ifelprt,1000)
      	  write(ifelprt,1010)
      	  do iz=1,nalt
           write(ifelprt,1020)altkm(iz),proelec(iz),
     .		(prodion(iz,m),m=1,4)
          enddo
      	  write(ifelprt,1015)
      	  do iz=1,nalt
           write(ifelprt,1020)altkm(iz),(prodion(iz,m),m=5,6)
          enddo
        endif
c
c 	Different productions
c     	produc  (alt,neutral specie,energy box)
c     	        = prod. at alt no iz,neutral specie j,box iener.
c     	        [cm-3.s-1.ev-1]
c     	proneut(alt,neutral specie)
c     	        = prod. at alt no iz,due to neutral specie j.
c     	        [cm-3.s-1]
c     	prodion (alt,ionized specie) = nbre d'ion de chaque espece .
c     	        = prod. at alt no iz,ion specie j.
c     	        [cm-3.s-1]
c     	prophel    (alt,boite d'energie)
c     	        = prod. at alt no iz,box iener, all species.
c     	        [cm-3.s-1.ev-1]
c     	proelec (alt)
c     	        = prod. at alt no iz, all species mixed.
c     	        [cm-3.s-1]
c
        if (iprt(9)==1)then
      	  write(ifelprt,1040)nen
	  nnz=(nalt-1)/4
	  iz=1
      	  do i=1,5
            write(ifelprt,1050) altkm(iz)
            do ien1=1,nen,2
              ien2=ien1+1
              write(ifelprt,1060)ien1,(produc(iz,neutspe,ien1),
     &         neutspe=1,3 ),ien2,(produc(iz,neutspe,ien2),neutspe=1,3)
            enddo
            iz=iz+nnz
          enddo
        endif

        if(iprt(10)==1)then
          nalto6 = (nalt+1)/6
          ialt1 = 1
          ialt2 = ialt1 + 5
          do ii = 1,nalto6
            write(ifelprt, 1070) (altkm(ialt),ialt = ialt1,ialt2)
            do ien = 1,nen
              write(ifelprt,1080)Ecent(ien),
     &		(prophel(ialt,ien),ialt=ialt1,ialt2)
 	    enddo
            ialt1 = ialt1 + 6
            ialt2 = ialt2 + 6
 	  enddo
c
      	  ialt2 = ialt2 - 5
      	  write(ifelprt, *)
      	  write(ifelprt, 1070) (altkm(ialt), ialt = ialt2,nalt)
      	  do ien = 1,nen
       	    write(ifelprt,1080) Ecent(ien),
     &		(prophel(ialt,ien),ialt=ialt2,nalt)
 	  enddo
        endif
c
        if(iprt(11)==1)then
      	  write(ifelprt,1110)
      	  do iz=1,nalt
      	   write(ifelprt,1100)altkm(iz),proelec(iz),
     &			      (proneut(iz,m),m=1,5)
          enddo
        endif
c
c 	calculates and prints the energy through primary production
c 	Integration over the energies
      do ialt = 1,nalt
 	    zwork(ialt) = 0.
	    do ien=1,nen
 	      zwork(ialt) = zwork(ialt)+
     &		                     Ecent(ien)*prophel(ialt,ien)*engdd(ien)
 	    enddo
      enddo
c
 	    call hint(nalt,altkm,zwork,qpheleV)
c
c 	Dans ce pgme, les altitudes sont en km :
 	    qpheleV = qpheleV*1.e+05
	    qphelerg = qpheleV*1.6022e-12
	    write(ifelprt,5007) qpheleV,qphelerg
c 	write(6,5007) qpheleV,qphelerg
5007	format(/,
     .  'The total E. from photoelectron prod.  is : ', 1pe11.4,
     .  ' eV/cm2/s',/,45x,'or : ',1pe11.4,' erg/cm2/s')
c
c     	write in file ifeltrans for transport program.
      	write(ifeltrans) nen,nalt,ns
      	write(ifeltrans) (Ecent(i),i=1,nen)
      	write(ifeltrans) (altkm(i),i=1,nalt)
      	write(ifeltrans) ((prophel(i,j),j=1,nen),i=1,nalt)
 	    write(ifeltrans) (proelec(i),i=1,nalt)
      	write(ifeltrans) ((prodion(ialt,isp),ialt=1,nalt),isp=1,ns+1)
      	close(ifeltrans)
      	close(ifelprt)

	  end subroutine prodprt
c
c------------------------ setpar ------------------------------
c
      subroutine setpar(kiappel,imod,nspec,idess,hrloc,UT,day,nan,
     .	   tempexo,f107,ap,glat,glong,nalt,year,
     .	   altkm,centE,botE,ddeng,Ecent,Ebot,engdd,iprt,pflux,knm,
     .	   tneutre,densneut,colden,iflux,wwt,number,xchap,chi,chideg,
     .	   sigi,sigt,ichapman)

      include 'comm.f'

        include 'TRANSPORT.INC'

        real, intent(in) :: f107(3),ap(7), knm, glong,altkm(nbralt),
     &     tneutre(nbralt),densneut(8,nbralt),colden(8,nbralt)
        integer, intent(in) :: nspec
        integer, intent(out) :: iflux,imod
        real, intent(out) :: nan

      	common /bloc/ threshold,nbseff,eVseff,seffion,sefftot,pfluxmin,
     .  	      pfluxmax,wave,eV,wavemin,wavemax,eVmin,eVmax,
     .		      nwave,ns,nns, f107min,f107max,iseff,wnmseff, lambdasr,
     .                sigsro2,Isr,lineflux,sigabso2,qyield,Isr2
        common /const/ pi,re,recm,bolt,gzero,amu
 	real wavemin(39),wavemax(39),eVmin(39),eVmax(39)
 	integer index(39)
       	real Ecent(nbren),engdd(nbren),Ebot(nbren)
       	real centE(nbren),ddeng(nbren),botE(nbren)
        real xchap(nbralt,nbrsp)
 	integer iprt(12),idess(9)
      	real seffion(2000,7),sefftot(2000,5),eVseff(2000),wnmseff(2000)
      	real sigi(39,7),sigt(39,5)
 	real wave(39),Pflux(39),Eflux(39),eV(39)
 	real pfluxmin(39),pfluxmax(39)
 	real wwt(7,6,39),threshold(7)
 	integer number(7),iseff
 	real trav(39),eVmid(39),work(39),place,year
 	integer yyddd,nbseff
 	real enflux
        real chapesp(8)
c
1000    format (1x,/,'hour=',f6.2, 8x,'chi =',f6.2, 3x,'degrees',
     .    6x,'exo temp =',f8.2,/)
1010  	format ('     felin.f : Production of primary photoelectrons',
     .    '. Sza (deg) :',f10.2)
c    .		'prod., Sza (deg) :',f10.2,/,' --------')
1020   	format (' Spectral lines :',/,'      Wavelength (nm)     ',7x,
     .    'energy(eV)        incident flux (photons/s.cm2)',/,
     .    51x,'SC#21REFW  Computed  F79059N',/,48x,3f10.2,/,
     .    1x,23('-'),2x,23('-'),2x,28('-'))
1030 	format(3f8.3,1x,3f8.3,3(1pe10.2))
1050  	format (' Ionization potentials (threshold, [eV]), ',
     .    'N2,O2,O,N2d,O2d,H,He ')
1060  	format( /,' Comparison 1/cos(chi) and',
     .          ' Chapman''s function (o,n2,o2)')
1070    format (1x, 7(1pe15.7) )
1080  	format (1x, 10(1pe10.2))
1090  	format (1x, f7.2, 6(1x,f5.2))
1130    format(/,70('-'),/,
     .  ' Lambda |      N_2 (*1.e18cm2)     ',
     .  '|    O_2   (*1.e18cm2)      |', /,70('-'),/,
     .  'ANGSTROMS  ABS   N_2+   N+     ION ',
     .   '|  ABS    O_2+    O+    ION', /,70('-'))
1140   format(12f7.3)
1150   format(/,28('-'),/,' Lambda  |  O   (*1.e18cm2) |',
     .    /,28('-'),/,'ANGSTROMS| ABS    ION       |', /,28('-'))
1160   format(/,70('-'),/,
     . ' Lambda |      H   (*1.e18cm2)     |    He  (*1.e18cm2)     |',
     .  /,70('-'),/,
     .  'ANGSTROMS      ABS        ION   |  ABS         ION',/,70('-'))
c
2010    format (/,'atmosphere neutre',/,f10.2, 20x, 'tempexo')
2020    format (' lat=' ,f7.2,3x,'long=',f7.2,3x,'ap=',f7.2)
2023    format (' jour=',f10.1,3x,'annee=',i7)
2030    format('flux moyen ',f7.2,3x,'flux de la veille ',f7.2)
2050    format (/,' no.   height    temp        n(N2)       '
     .      ,' n(02)      n(01)',/,8x,' km      deg k    ',
     .       '   /cm3       /cm3       /cm3')
2055    format (/,' no.   height     n(H)       n(He)       ',
     .    'n(N)       n(A)', /,'         km        /cm3       ',
     .    '/cm3       /cm3       /cm3')
2060    format (i3,2f10.2,3(1pe11.3))
2065    format (i3,f10.2,4(1pe11.3))
c
c     	ns is Number of neutral Species in atmosphere
c     	nns = ns+2 allows for diss ioniz of o2 and n2
 	    ns  = nspec
 	    nns = ns + 2

 	    nan = ifix(year)
c
c       Met les energies en ordre decroissant
        print *,'felin:setpar: nen=',nen
        if(centE(1).gt.centE(nen))then
          do ien = 1,nen
            Ebot(ien)  = botE(ien)
            Ecent(ien) = centE(ien)
            engdd(ien) = ddeng(ien)
          enddo
        else
          do ien = 1,nen
            Ebot(nen+1-ien)  = botE(ien)
            Ecent(nen+1-ien) = centE(ien)
            engdd(nen+1-ien) = ddeng(ien)
          enddo
        endif

        print *,'hi'
        call xline(13,ifeldat)
        print *,'bye'
        read(ifeldat,*) (idess(i),i=1,8)
c	iprt(1) = 1 if print Energy boxes definition [eV, angstrom, cm].
c	iprt(2) = 1 if print Spectral lines.
c	iprt(3) = 1 if print of cross sections
c	iprt(4) = 1 if print Ionization  potentials (threshold, [eV]).
c	iprt(5) = 1 if print Normalized branching ratio.
c	iprt(6) = 1 if print Neutral atmosphere.
c	iprt(7) = 1 if print Comparison 1/cos(chi) and Chapman's fction.
c	iprt(8) = 1 if print Electron and ion production(/cm3).
c	iprt(9) = 1 if print Production (/cm3),at 5 alt. and each en.
c	iprt(10)= 1 if print Produced Photoel. (/cm3.sec) at each alt.
c	iprt(11)= 1 if print Electron and ion production(/cm3).
c	iprt(12)= available.
        call xline(13,ifeldat)
        read(ifeldat,*) (iprt(i),i=1,12)
        write(ifelprt,*)
        write(ifelprt,*)'fcwave = fudge factor to multiply the flux'
        write(ifelprt,*)'values below 257 angstrom (Published by Torr):'
        call xline(1,ifeldat)
        read(ifeldat,*) fcwave
        write(ifelprt,*)fcwave
        call xline(5,ifeldat)
        read(ifeldat,*) iflux
        call xline(1,ifeldat)
        read(ifeldat,*) iseff
        call xline(3,ifeldat)
        read(ifeldat,*) imod
        call xline(1,ifeldat)
        read(ifeldat,*) ichapman

      close(ifeldat)


c----   Lecture des sections efficaces
c 	sefftot = sections eff. totale, N2,O2,O
c 	sefftion = sections eff. d'ionisation et d'ionisation dissocia
c 	tive, N2,O2,O,Ndiss,Odiss.
 	  if(iseff==1)then
c 	  Fichier Torr et Torr de 39 energies (2 energies rajoutees
c 	  pour tenir compte du flux solaire de Tobiska, jl 1993)
 	  open(icrsphot,
     .	    file='dir.data/dir.linux/dir.cine/dir.seff/crsphot1.dat',
     &      status='old')
      	  call xline(6,icrsphot)
 	  read(icrsphot,*)nbseff
 	  do i = 1,nbseff
 	    read(icrsphot,*)wnmseff(i),sefftot(i,1),seffion(i,1),
     .		seffion(i,4),sigion1,sefftot(i,2),seffion(i,2),
     .		seffion(i,5),sigion2
 	    eVseff(i) = 1239.8/wnmseff(i)
 	  enddo
      	  call xline(6,icrsphot)
 	  do i = 1,nbseff
 	    read(icrsphot,*)wnm,sefftot(i,3),seffion(i,3)
 	  enddo
      	  call xline(10,icrsphot)
 	  do i = 1,nbseff
 	    read(icrsphot,*)wnm,sefftot(i,4),seffion(i,6),
     .		          sefftot(i,5),seffion(i,7)
 	  enddo
 	  close(icrsphot)
 	  elseif(iseff==2)then
c 	  Fichier Fennely et Torr de 1946 energies (dont 2 rajoutees
c 	  pour tenir compte du flux solaire de Tobiska, jl 1993)
 	  open(icrsphot,
     .	   file = 'dir.data/dir.linux/dir.cine/dir.seff/crsphot2.dat',
     &     status ='old')
      	  call xline(12,icrsphot)
 	  read(icrsphot,*)nbseff
 	  do i = 1,nbseff
 	    read(icrsphot,*)wA,sefftot(i,1),seffion(i,1),seffion(i,4),
     .  	  sigion1,sefftot(i,2),seffion(i,2),seffion(i,5),sigion2
 	    wnmseff(i) = wA/10.
 	    eVseff(i) = 1239.8/wnmseff(i)
c	    Bizzarement, certaines seff d'abs. sont inf a la seff
c 	    d'ionis. totale...
 	    sefftot(i,1) = max(sefftot(i,1),(seffion(i,1)+seffion(i,4)))
 	    sefftot(i,2) = max(sefftot(i,2),(seffion(i,2)+seffion(i,5)))
 	  enddo
      	  call xline(17,icrsphot)
 	  read(icrsphot,*)nbseff
 	  do i = 1,nbseff
 	    read(icrsphot,*)wnm,sefftot(i,3),x1,x2,x3,x4,x5,x6,x7,
     .		seffion(i,3)
 	    sefftot(i,3) = max(sefftot(i,3),seffion(i,3))
 	  enddo
      	  call xline(14,icrsphot)
 	  read(icrsphot,*)nbseff
 	  do i = 1,nbseff
 	    read(icrsphot,*)wnm,sefftot(i,4),seffion(i,6),
     . 	                  sefftot(i,5),seffion(i,7)
 	    sefftot(i,4) = max(sefftot(i,4),seffion(i,6))
 	    sefftot(i,5) = max(sefftot(i,5),seffion(i,7))
 	  enddo
 	  close(icrsphot)
 	  endif
c
c ----	Estimation du flux solaire sur les 39 energies initiales,
c 	en photons/cm2/s
 	  if(iflux == 0)then
	  do i=1,nwave
 	    pflux(i)= (f107(1)-f107max)*(pfluxmin(i)-pfluxmax(i))/
     . 	   		(f107min-f107max) + pfluxmax(i)
c	    pflux(i)= (130.5-f107max)*(pfluxmin(i)-pfluxmax(i))/
c    . 	   		(f107min-f107max) + pfluxmax(i)

 	  enddo
 	  elseif (iflux==1)then
 	  do i=1,nwave
 	    pflux(i)= pfluxmin(i)
 	  enddo
 	  elseif (iflux==2)then
 	  do i=1,nwave
 	    pflux(i)= pfluxmax(i)
 	  enddo
 	  elseif (iflux==3)then
 	  if (nan > 99)then
 	    yyddd = (nan-1900)*1000+ifix(day)
 	  else
 	    yyddd = nan*1000+ifix(day)
 	  endif
 	  write(6,*)'yyddd = ',yyddd
 	  if(yyddd.lt.68172. or. yyddd.gt.88366)then
 	    write(6,*)'No source flux in Tobiska file for this date'
 	    write(6,*)'yyddd must be 68172<= yyddd <= 88366'
 	    write(6,*)'Program stopped'
 	    stop
 	  endif
 	  numday = 1
 	  call euv91(ifelprt,yyddd,numday,wavemin,eVmin,wavemax,
     .		     eVmax,eV,wave,Pflux,Eflux,iprt(1))

        elseif (iflux==4)then
c   	  ici le flux EUVAC94 peut etre obtenu.
c   	  ref: Richards et al, JGR,99,8981,1994
c   	  modification : OW, 12 dec 1997
          open (47,file=
     .    'dir.data/dir.linux/dir.cine/dir.euvac/EUVAC.dat',
     &    status='old')
 	  do i=1,nwave
            p=(f107(1)+f107(3))/2.
            read(47,*,end=333)n,wave1,wave2,fbase,a
            pflux(i)=fbase*(1.+a*(p-80))*1e9
          enddo
 	  close(47)
333       continue

 	  endif
c 	Prend en compte la distance au soleil
 	do i=1,nwave
	  pflux(i) = pflux(i)/(rayonUA**3)
 	enddo
c
      	write(6,1010)chideg
  	write(ifelprt, 1000) hrloc, chideg, tempexo 
c
c ---- 	Computes the energy at which to put the solar flux
c   	eV = energy of uv line [eV]  . eV=(hc)/(q.line.1E-9)
c	wavemin, wavemax = longueur d'onde en nm.
 	place = 0.5
        write(ifelprt,*)'Place de l''energie dans la grille ',place
 	do ien = 1,39
 	  eVmin(ien) = 1239.8/wavemax(ien)
 	  eVmax(ien) = 1239.8/wavemin(ien)
  	  eVmid(ien)=eVmin(ien)+(eVmax(ien)-eVmin(ien))*place
 	enddo
 	call clasdesc(eVmid,39,index)
c 	En sortie de clasdesc, les energies sont reordonnees en ordre
c 	descendant.
 	do ien = 1,39
  	  eV(ien)=eVmid(ien)
          wave(ien)=1239.8/eV(ien)
 	enddo
c 	On reordonne tout ce qui depend des energies.
 	call reord(wavemin,39,index,trav)
 	call reord(wavemax,39,index,trav)
 	call reord(evmin,39,index,trav)
 	call reord(evmax,39,index,trav)
 	call reord(pfluxmin,39,index,trav)
 	call reord(pfluxmax,39,index,trav)
 	call reord(pflux,39,index,trav)
        enflux = 0.
        do ien = 1,39
          enflux = enflux + pflux(ien)*eV(ien)
        enddo
c
c 	Takes the Torr fudge factor into account
 	do i=1,nwave
	  if(wave(i).lt.25.7)pflux(i) = pflux(i)*fcwave
 	enddo
c
c ----	Computes the normalized branching ratio
 	    call  branchratio(39,eV,threshold,number,wwt)
c
c 	Interpole les sections efficaces sur la grille d'energie de
c 	travail du flux solaire
 	    do iwave = 1,nwave
 	      sigt(iwave,1) = sefint(eV(iwave),1,1)
     	  sigt(iwave,2) = sefint(eV(iwave),2,1)
     	  sigt(iwave,3) = sefint(eV(iwave),3,1)
 	      sigt(iwave,4) = sefint(eV(iwave),4,1)
 	      sigt(iwave,5) = sefint(eV(iwave),5,1)
 	      sigi(iwave,1) = sefint(eV(iwave),1,2)
 	      sigi(iwave,2) = sefint(eV(iwave),2,2)
 	      sigi(iwave,3) = sefint(eV(iwave),3,2)
 	      sigi(iwave,4) = sefint(eV(iwave),4,2)
 	      sigi(iwave,5) = sefint(eV(iwave),5,2)
 	      sigi(iwave,6) = sefint(eV(iwave),6,2)
 	      sigi(iwave,7) = sefint(eV(iwave),7,2)
 	    enddo
c
c
c ----	Ecritures
c
        if (iprt(6)==1) then
          write(ifelprt,2010)tempexo
          write(ifelprt,2020)glat,glong,ap(1)
  	      write(ifelprt,2023)day,nan
          write(ifelprt,2030)f107(3),f107(1)
          write(ifelprt,2050)
	  if(nspec <= 3)then
	    iwrite=nspec
	    do 130 i=1,nalt
              write(ifelprt,2060)i,altkm(i),tneutre(i),
     .                          (densneut(isp,i),isp=1,iwrite)
130 	    continue
	  else
	    iwrite=3
	    do 140 i=1,nalt
              write(ifelprt,2060)i,altkm(i),tneutre(i),
     .                          (densneut(isp,i),isp=1,iwrite)
140 	    continue
            write(ifelprt,2055)
	    do 150 i=1,nalt
              write(ifelprt,2065)i,altkm(i),
     .			(densneut(isp,i),isp=4,nspec)
150 	    continue
	  endif
        endif
c
2070    format('total energy contained in the solar flux',
     .          ' [eV/s/cm2]:',1pe10.2)
        if(kiappel==1)write(6,2070)enflux
        if(iprt(2)==1)then
  	  write(ifelprt,1020)f107min, f107(1), f107max
	  do 35 i=1,nwave
	    write(ifelprt,1030)wavemin(i),wave(i),wavemax(i),eVmin(i),
     .		eV(i),eVmax(i),pfluxmin(i),pflux(i),pfluxmax(i)
35	  continue
          write(ifelprt,2070)enflux
	endif
c
 	if (iprt(3)==1)then
 	  if(iseff==1) write(ifelprt,*)
     .			'Torr and Torr (1985) cross section set'
 	  if(iseff==2) write(ifelprt,*)
     .			'Fennelly and Torr (1992) cross section set'
 	  write(ifelprt,1130)
 	  do i = 1,39
 	    write(ifelprt,1140)wave(i),sigt(i,1)*1.e18,sigi(i,1)*1.e18,
     .	        sigi(i,4)*1.e18,(sigi(i,1)+sigi(i,4))*1.e18,
     . 		sigt(i,2)*1.e18,sigi(i,2)*1.e18,sigi(i,5)*1.e18,
     .		(sigi(i,2)+sigi(i,5))*1.e18
 	  enddo
 	  write(ifelprt,1150)
 	  do i = 1,39
 	    write(ifelprt,1140)wave(i),sigt(i,3)*1.e18,sigi(i,3)*1.e18
 	  enddo
 	  write(ifelprt,1160)
 	  do i = 1,39
 	    write(ifelprt,1140)wave(i),sigt(i,4)*1.e18,sigi(i,6)*1.e18,
     .			               sigt(i,5)*1.e18,sigi(i,7)*1.e18
 	  enddo
 	endif
c
 	if (iprt(4)==1)then
      	  write(ifelprt,1050)
      	  do 100 j = 1,nns
            write(ifelprt,1080)threshold(j)
100  	  continue
 	endif
c
 	if (iprt(5)==1)then
      	  write(ifelprt,*)
      	  write(ifelprt,*)'Normalized branching ratio'
      	  write(ifelprt,*)
c
	  write(ifelprt,*)' N2 --> N2+ '
 	  iontype = 1
	  write(ifelprt,*)' Energy  X2S+g A2Pu B2S+u F2Su  2S+g'
      	  do k = 1,nwave
            write(ifelprt,1090)eV(k),
     .		 (wwt(iontype,ist,k),ist =1,number(iontype))
 	  enddo
	  write(ifelprt,*)' O2 --> O2+ '
 	  iontype = 2
	  write(ifelprt,*)' Energy  X2Pg  a4Pu+A2Pu  b4S-g'
      	  do k = 1,nwave
            write(ifelprt,1090)eV(k),
     .		 (wwt(iontype,ist,k),ist =1,number(iontype))
 	  enddo
	  write(ifelprt,*)' O --> O+ '
 	  iontype = 3
	  write(ifelprt,*) ' Energy  2P3/   2P3/  2P3/  2P4/  2P4/'
	  write(ifelprt,*) '         4S0    2D0   2D0    4P    2P'
      	  do k = 1,nwave
            write(ifelprt,1090)eV(k),
     .		 (wwt(iontype,ist,k),ist =1,number(iontype))
 	  enddo
	  write(ifelprt,*)' N2 --> N + N+ '
 	  iontype = 4
	  write(ifelprt,*)' Energy  Fundamental state'
      	  do k = 1,nwave
            write(ifelprt,1090)eV(k),
     .		 (wwt(iontype,ist,k),ist =1,number(iontype))
 	  enddo
	  write(ifelprt,*)' O2 --> O + O+ '
 	  iontype = 5
	  write(ifelprt,*)
     .   ' Energy  B2S-g  2Pu  c4S-u  2S-u 2,4S-g 662A'
      	  do k = 1,nwave
            write(ifelprt,1090)eV(k),
     .		 (wwt(iontype,ist,k),ist =1,number(iontype))
 	  enddo
	  write(ifelprt,*)' H --> H+ '
 	  iontype = 6
	  write(ifelprt,*) ' Energy  One state'
      	  do k = 1,nwave
            write(ifelprt,1090)eV(k),
     .		 (wwt(iontype,ist,k),ist =1,number(iontype))
 	  enddo
	  write(ifelprt,*)' He --> He+ '
 	  iontype = 7
	  write(ifelprt,*) ' Energy  One state'
      	  do k = 1,nwave
            write(ifelprt,1090)eV(k),
     .		 (wwt(iontype,ist,k),ist =1,number(iontype))
 	  enddo
 	endif
c
c       Chapman:
c	--------
c     	if(iprt(7)==1)write(ifelprt,1060)
c
 	if(ichapman==0) then
      	  do iz = 1,nalt
            call fchap(iyd,UT,altkm(iz),glat,glong,f107,
     &                  ap,chi,chapesp)
            do isp=1,ns
              xchap(iz,isp)=chapesp(isp)
            enddo
 	  enddo
 	else if (ichapman==1)then
      	  do iz = 1,nalt
 	    altcm=altkm(iz)*1.e+05
 	    do  isp=1,ns
              xchap(iz,isp) =
     .		chapsmith(chi,altcm,tneutre(iz),atomas(isp))
 	    enddo
 	  enddo
        else if(ichapman==2)then
      	  do iz = 1,nalt
 	    altcm=altkm(iz)*1.e+05
 	    do  isp=1,ns
              xchap(iz,isp) =
     .		chapgreen(chi,altcm,tneutre(iz),atomas(isp))
 	    enddo
 	  enddo
 	endif
c       Comparison:
 	if(iprt(7)==1)then
	  if(cos(chi).ne.0.) xcos = 1/cos(chi)
 	  do iz=1,nalt
            write(ifelprt, 1070) altkm(iz),xcos,(xchap(iz,j), j = 1,ns)
 	  enddo
 	endif

	end subroutine setpar

c------------------------ search ------------------------------

        subroutine search (Ebot,engdd,dele,iener)
        implicit none
        
        include 'TRANSPORT.INC'
	
!       Search in which box ll to put the created electron.
        real,intent(in) :: Ebot(nbren),engdd(nbren),dele
        integer, intent(out) :: iener

        real sup,inf
        integer ien

      do ien=1,nen
         inf = Ebot(ien)
         sup = Ebot(ien)+engdd(ien)
         if(dele.ge.inf .and. dele.lt.sup) then
            iener=ien
            return
         endif
      enddo
      
      if (dele.ge.Ebot(1)+engdd(1)) iener = 1
      if (dele.lt.Ebot(nen)) iener = nen

      end subroutine search
c
c------------------------ Function sefint ---------------------------
c
      function sefint(ener,isp,type)
c
      implicit none
      include 'TRANSPORT.INC'
c
      	common /bloc/ threshold,nbseff,eVseff,seffion,sefftot,pfluxmin,
     &  	      pfluxmax,wave,eV,wavemin,wavemax,eVmin,eVmax,
     &		      nwave,ns,nns,f107min,f107max,iseff,wnmseff,lambdasr,
     &                sigsro2,Isr,lineflux,sigabso2,qyield,Isr2
      	real seffion(2000,7),sefftot(2000,5),eVseff(2000),threshold(7)
        real wnmseff(2000)
        real eV(39),wave(39),pfluxmin(39),pfluxmax(39)
        real wavemin(39),wavemax(39),eVmin(39),eVmax(39)
        real f107min,f107max,ener,sefint
        integer ns,nns,isp,nwave,nbseff,type,iseff
        real sigsro2(8),lambdasr(9),Isr(8)
        real lineflux(7),sigabso2(7),qyield(7),Isr2(8)

 	real tabin(2000),tabeV(2000),tabout(1),seffout(1)
 	integer i
c 	type = 1 --> total
c 	type = 2 --> ionisation
c
 	if(type==1)then
 	  do i = 1,nbseff
 	    tabin(i) = sefftot(nbseff+1-i,isp)
 	    tabeV(i) = eVseff(nbseff+1-i)
 	  enddo
 	elseif(type==2)then
 	  do i = 1,nbseff
 	    tabin(i) = seffion(nbseff+1-i,isp)
 	    tabeV(i) = eVseff(nbseff+1-i)
 	  enddo
 	endif
 	tabout(1) = ener
c
 	call intlin(nbseff,tabeV,tabin,1,tabout,seffout)
c
 	sefint = seffout(1)*1.e-18

      end function sefint

c
c------------------------ Function sperfc ---------------------------
c
      pure real function sperfc(dummy)
      implicit none
      real,intent(in) :: dummy
c
c 	Used to compute the Chapman function.
c
      	if (dummy <= 8.) then
          sperfc = (1.0606963+0.55643831*dummy) /
     &           (1.0619896+1.7245609*dummy+dummy*dummy)
        else
          sperfc=0.56498823/(0.06651874+dummy)
        endif

      end function sperfc
c
c----------------------------------------------------------------------
c
        subroutine reord(tab,ntab,indx,trav)
        implicit none
        integer,intent(in) :: indx(ntab),ntab
        real,intent(inout)   :: tab(ntab),trav(ntab)
        integer i

        do i = 1,ntab
            trav(i) = tab(indx(i))
        enddo
        do i = 1,ntab
           tab(i) = trav(i)
        enddo

        end subroutine reord
c
c-------------------- donnees -------------------------------
c
 	block data

      	common /bloc/ threshold,nbseff,eVseff,seffion,sefftot,pfluxmin,
     .  	      pfluxmax,wave,eV,wavemin,wavemax,eVmin,eVmax,
     .		      nwave,ns,nns,f107min,f107max,iseff,wnmseff,lambdasr,
     .                sigsro2,Isr,lineflux,sigabso2,qyield,Isr2
 	common /const/ pi,re,recm,bolt,gzero,amu
c
 	real wavemin(39),wavemax(39),eVmin(39),eVmax(39),f107max
      	real seffion(2000,7),sefftot(2000,5),eVseff(2000),wnmseff(2000)
       	dimension threshold(7),eV(39),wave(39)
 	dimension pfluxmin(39),pfluxmax(39)
        real sigsro2(8),lambdasr(9),Isr(8)
        real lineflux(7),sigabso2(7),qyield(7),Isr2(8)
c
 	real eVbrO(10),brO(5,10),eVbrN2(13),brN2(5,13),
     .  	eVbrO2(20),brO2(3,20),eVbrO2dis(14),brO2dis(6,14),
     .		eVbrN2dis(2),brN2dis(1,2),eVbrH(2),brH(1,2),
     .		eVbrHe(2),brHe(1,2)
 	integer neVbrO,nstO,neVbrN2,nstN2,
     .  	neVbrO2,nstO2,neVbrO2dis,nstO2dis,neVbrN2dis,nstN2dis,
     .  	neVbrH,nstH,neVbrHe,nstHe
 	common /branching/
     .		eVbrN2,    brN2,    neVbrN2,    nstN2,
     .		eVbrO2,    brO2,    neVbrO2,    nstO2,
     .		eVbrO,     brO,     neVbrO,     nstO,
     .		eVbrO2dis, brO2dis, neVbrO2dis, nstO2dis,
     .		eVbrN2dis, brN2dis, neVbrN2dis, nstN2dis,
     .		eVbrH,     brH,     neVbrH,     nstH,
     .		eVbrHe,    brHe,    neVbrHe,    nstHe
c

        data lambdasr /1350.,1400.,1450.,1500.,1550.,1600.
     .                  ,1650.,1700.,1750./
        data Isr /.0628e10,.0768e10,.1078e10,.1894e10,
     .            .2436e10,.3120e10,.6430e10,1.132e10/
        data sigsro2 /1.2e-17,1.5e-17,1.3e-17,1.0e-17,
     .                6.e-18,3.4e-18,1.5e-18,5.e-19/
        data lineflux /2.07e9,2.7e8,4.36e9,1.12e9,1.78e9,
     .                  8.92e8,3.48e9/
        data sigabso2 /4.0e-18,1.75e-18,1.58e-18,1.e-18,
     .                  7.8e-19,1.8e-18,1.3e-18/
        data qyield /0.066,0.074,0.094,0.097,0.100,0.101,
     .                  0.012/
        data Isr2 /6.1e9,9.5e9,16.2e9,25.2e9,35.6e9,56.e9,
     .                  121.5e9,225.e9/

 	integer ne
c 	WARNING! In a block data, an array T(i,j) is read by i begining
c 	and j ending.
c
c     	threshold(isp) = potentiel d'ionisation [eV], espece isp
c	(Rees, Physics and chemistry of the upper atmosphere,
c 	cambridge, 1989)
c 	           1   = N2 --> N2+
c 	           2   = 02 --> 02+
c 	           3   = 0  --> 0+
c 	           4   = N2 --> N+ + N
c 	           5   = O2 --> O+ + O
c 	           6   = H  --> H+
c 	           7   = He --> He+
 	data threshold /15.6, 12.08, 13.6, 24.3, 19.0,13.6,24.6/
c
      	data re/6370.00/,recm/6.37E8/
 	data bolt/1.38e-16/,gzero/980.665/,amu/1.662e-24/
c
      	data nwave /39/
c   	wave = Wavelength, in nm of the considered line.
         data wavemin /
     .     1.900,    3.000,    5.000,   10.000,   15.000,
     .    20.000,   25.630,   28.415,   25.000,   30.331,
     .    30.378,   30.000,   36.807,   35.000,   40.000,
     .    46.522,   45.000,   50.000,   55.437,   58.433,
     .    55.000,   60.976,   62.973,   60.000,   65.000,
     .    70.331,   70.000,   76.515,   77.041,   78.936,
     .    75.000,   80.000,   85.000,   90.000,   97.702,
     .    95.000,  102.572,  103.191,  100.000/
         data wavemax /
     .     3.000,    5.000,   10.000,   15.000,   20.000,
     .	  25.000,   25.630,   28.415,   30.000,   30.331,
     .    30.378,   35.000,   36.807,   40.000,   45.000,
     .    46.522,   50.000,   55.000,   55.437,   58.433,
     .    60.000,   60.976,   62.973,   65.000,   70.000,
     .    70.331,   75.000,   76.515,   77.041,   78.936,
     .    80.000,   85.000,   90.000,   95.000,   97.702,
     .   100.000,  102.572,  103.191,  105.000/
c   	pfluxmin = Solar input flux for sunspot minimum.
c 	2 lowest wavelength --> Tobiska, 1991
c 	Others --> Torr, 1985
	data f107min /68.00/
        data pfluxmin /
     .	      1.119e+07, 3.119e+07,
     .	      3.834e+08, 1.346e+08, 18.418e+08, 9.235e+08 , 2.713e+08,
     .	      1.000e+08, 8.405e+08, 2.350e+08 , 6.000e+09 , 8.661e+08,
     .	      7.394e+08, 2.121e+08, 3.926e+08 , 1.800e+08 , 3.063e+08,
     .	      5.085e+08, 7.992e+08, 1.580e+09 , 4.843e+08 ,4.500e+08,
     .	      1.500e+09, 1.746e+08, 2.223e+08 , 3.915e+08 , 1.667e+08,
     .	      1.997e+08, 2.425e+08, 7.931e+08 , 8.728e+08 , 19.311e+08,
     .	     44.325e+08, 42.170e+08,59.570e+08, 17.850e+08,
     .	     43.750e+08, 31.840e+08,36.401e+08/
c
c   	pfluxmax = Solar input flux for sunspot maximum.
	data f107max /243.00/
        data pfluxmax /
     .	     1.2330e+08,  2.4690e+08,
     .       1.1487E+09, 3.4330E+08, 4.8498E+09, 3.7013E+09, 5.9470E+08,
     .       3.1675E+09, 4.1358E+09, 2.4995E+09, 1.1280E+10, 5.6326E+09,
     .       1.3949E+09, 2.1965E+09, 9.9320E+08, 3.6210E+08, 1.6716E+09,
     .       1.5468E+09, 1.5904E+09, 4.8664E+09, 1.0213E+09, 1.4621E+09,
     .       3.0180E+09, 4.8200E+08, 4.5540E+08, 7.1650E+08, 4.2560E+08,
     .       4.3180E+08, 6.7090E+08, 1.5869E+09, 2.1809E+09, 5.0135E+09,
     .       1.3298E+10, 1.2035E+10, 1.3177E+10, 4.4204E+09,
     .       1.3125E+10, 9.0426E+09, 8.6669E+09/
c
c	Les branching ratio qui suivent viennent de Rees, Physics and
c 	chemistry of the upper atmosphere, cambridge, 1989.
	data neVbrO /10/
 	data eVbrO /247.96,39.48,39.36,28.50,28.43,18.64,18.62,
     .			16.94,16.91,13.60/
 	data nstO /5/
c 		  2P3-4S0  2P3-2D0   2P3-2D0   2P4-4P   2P4-2P
 	data brO / 0.25,     0.36,    0.23,     0.10,    0.06,
     . 	           0.25,     0.36,    0.23,     0.10,    0.06,
     .		   0.26,     0.41,    0.26,     0.07,    0.00,
     .		   0.26,     0.41,    0.26,     0.07,    0.00,
     .		   0.30,     0.45,    0.25,     0.00,    0.00,
     .		   0.30,     0.45,    0.25,     0.00,    0.00,
     .		   0.43,     0.57,    0.00,     0.00,    0.00,
     .		   0.43,     0.57,    0.00,     0.00,    0.00,
     .		   1.00,     0.00,    0.00,     0.00,    0.00,
     .		   1.00,     0.00,    0.00,     0.00,    0.00/
c
	data neVbrN2 /13/
 	data eVbrN2 /59.04,51.66,44.28,41.33,37.34,28.97,24.80,
     .		     20.66,18.79,18.78,17.22,16.60,15.6/
 	data nstN2 /5/
c 		   X2S+g  A2Pu   B2S+u   F2Su   2S+g
 	data brN2 /0.271, 0.275, 0.110, 0.064, 0.278,
     .		   0.271, 0.345, 0.110, 0.064, 0.210,
     .		   0.271, 0.470, 0.095, 0.040, 0.124,
     .		   0.271, 0.470, 0.110, 0.074, 0.075,
     .		   0.300, 0.520, 0.120, 0.060, 0.000,
     .		   0.460, 0.460, 0.080, 0.000, 0.000,
     .		   0.404, 0.506, 0.090, 0.000, 0.000,
     .		   0.308, 0.589, 0.103, 0.000, 0.000,
     .		   0.308, 0.589, 0.103, 0.000, 0.000,
     .		   0.308, 0.692, 0.000, 0.000, 0.000,
     .		   0.420, 0.580, 0.000, 0.000, 0.000,
     .		   1.000, 0.000, 0.000, 0.000, 0.000,
     .		   1.000, 0.000, 0.000, 0.000, 0.000/
c
	data neVbrO2 /20/
 	data eVbrO2 /40.78,38.38,27.31,26.89,24.60,23.09,22.30,21.64,
     .		     21.23,20.73,20.32,19.46,19.22,19.02,18.13,17.61,
     .		     17.22,16.82,16.02,12.1/
 	data nstO2 /3/
c	            X2Pg  a4Pu+A2Pu  b4S-g
   	data brO2 /0.365,   0.205,   0.125,
     .		   0.374,   0.210,   0.124,
     .		   0.432,   0.243,   0.120,
     .		   0.435,   0.245,   0.120,
     . 		   0.384,   0.270,   0.126,
     .		   0.345,   0.290,   0.130,
     .		   0.356,   0.230,   0.225,
     .		   0.365,   0.270,   0.216,
     .    	   0.306,   0.330,   0.210,
     .		   0.230,   0.295,   0.375,
     .		   0.235,   0.385,   0.305,
     .		   0.245,   0.350,   0.370,
     .		   0.340,   0.305,   0.330,
     .		   0.270,   0.385,   0.345,
     .		   0.482,   0.518,   0.000,
     .		   0.675,   0.325,   0.000,
     .		   0.565,   0.435,   0.000,
     .		   0.565,   0.435,   0.000,
     .		   1.000,   0.000,   0.000,
     .		   1.000,   0.000,   0.000/
c
c 	On cree N a un seul etat d'excitation.
	data neVbrN2dis /2/
 	data eVbrN2dis /500.,24.3/
 	data nstN2dis /1/
   	data brN2dis /1.,1./
c
	data neVbrO2dis /14/
 	data eVbrO2dis /40.78,38.38,27.31,26.89,24.60,23.09,22.30,21.64,
     .		     21.23,20.73,20.32,19.46,19.22,19.00/
 	data nstO2dis /6/
c 	           B2S-g    2Pu      c4S-u     2S-u     2,4S-g     662A
   	data brO2dis
     .		  /0.055,  0.060,    0.035,   0.030,    0.125,    0.000,
     .		   0.055,  0.060,    0.035,   0.030,    0.000,    0.112,
     .		   0.055,  0.060,    0.035,   0.000,    0.000,    0.055,
     .		   0.055,  0.060,    0.035,   0.000,    0.000,    0.050,
     .		   0.079,  0.026,    0.000,   0.000,    0.000,    0.115,
     .		   0.098,  0.000,    0.000,   0.000,    0.000,    0.137,
     .		   0.109,  0.000,    0.000,   0.000,    0.000,    0.080,
     .		   0.119,  0.000,    0.000,   0.000,    0.000,    0.030,
     .		   0.125,  0.000,    0.000,   0.000,    0.000,    0.030,
     .		   0.058,  0.000,    0.000,   0.000,    0.000,    0.045,
     .		   0.000,  0.000,    0.000,   0.000,    0.000,    0.075,
     .		   0.000,  0.000,    0.000,   0.000,    0.000,    0.036,
     .		   0.000,  0.000,    0.000,   0.000,    0.000,    0.025,
     .		   0.000,  0.000,    0.000,   0.000,    0.000,    0.000/
c
c 	on cree H et He a un seul etat d'
c 	excitation.
	data neVbrH /2/
 	data eVbrH /500.,13.6/
 	data nstH /1/
 	data brH / 1.,1./
c
	data neVbrHe /2/
 	data eVbrHe /500.,24.6/
 	data nstHe /1/
 	data brHe / 1.,1./
c

 	end block data
