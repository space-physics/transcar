	   subroutine degrad(knm,neninit,centE,botE,ddeng,nspec,kiappel)
c
c----------------------------------------------------------------------!
c 	version degrad3.f					       !
c 	(C) Copyright by D.Lummerzheim, (1984,...      		       !
c 	(C)          and J.Lilensten (1992...        		       !
c 	This is a pre-TRANS code which prepares the cross-section data !
c 	for the main TRANS code					       !
c 	Les sections efficaces sont lues espece par espece, et dans
c 	chaque cas, les grilles d'energie peuvent etre differents.
c	Cross sections are read species by species, and in each case, the energy grids can be different.
c----------------------------------------------------------------------!
c
c e(n)			= Energy grid, increasing order (e(1) = min)
c engdd 		= Total width of each cell
c cel(isp,ien)		= elastic cross section (not used here)
c cin(isp,ien)		= total inelastic cross section
c nexcst(isp) 		= excited state number of specie isp. THE LAST
c 			  STATE IS FOR IONIZATION.
c			  (ex jsg(isp))
c nionst(isp) 		= actual number of excited ion states
c 		 	  (ex jsp(isp))
c cinex(isp,js,ien) 	= inelastic cross section for each state.
c			  CINEX(isp,JSG,ien) --> ionization.
c omdeg(E1,E2,isp)	= redist. matrix.
c			  E1 > E2 ==> degraded
c			  E2 > E1 ==> secondary
********************* Single precision ! *******************************

	    implicit logical (l)
	    include 'TRANSPORT.INC'
c
        real e(nbren),engdd(nbren)
        real esig(nbren),bsig(nbren),delta(nbren)
        integer nen
        common /eng_/ e,engdd,nen
c
        real ethres(nbrsp,nbrexc,nbrionst),bratio(nbrionst,nbrsp)
        integer nexcst(nbrsp),nionst(nbrsp)
        common /eng2_/ ethres,bratio,nexcst,nionst
c
       real cel(nbrsp,nbren),cin(nbrsp,nbren),cinex(nbrsp,nbrexc,nbren)
        common /cros_/ cel,cin,cinex
c
 	    dimension centE(nbren),ddeng(nbren),botE(nbren)
c 	botE(nbren) ne sert a rien ici, si ce n'est a permettre au
c 	programme simplifie pour la raie de plasma de connaitre
c 	les bas de grille pour la photoionisation.
c 	Il faut donc le laisser...!!!...
c
	    common /rdist_/ omdeg(nbren,nbren,nbrsp)
	    common /distfct_/lopal

*	The redistribution function as calculated by the SWARTZ
*	code.
*	OMDEG(N1,N2,J) : cross-section for inelastic scattering
*			  from energy N1 to N2
*	OMDEG(N2,N1,J)  : cross-section for the production of a
*			  secondary with energy N2

	character*8 istdate,isttime
	character*9 title(nbrexc,nbrsp)
 	character*9 excion(nbrionst,nbrsp)
	character(len=80) crsin,crsout,rdtout,preout
	character(len=80) line
 	integer idess(2),lenc,kiappel,i1,i2
c
1000    format(a)
1010	format(' Pre-electron transport code version 3 ',/,a8,3x,a8/)
1020	format(' The maximum input parameter authorized are:',/,
     .    ' length of energy grid (nbren) :',i5,/,
     .    ' half length of angle grid (nbrango2) :',i5,/,
     .    ' number of species (nbrsp) :',i2,/,
     .    ' number of states (nbrexc) :',i2,/,
     .    ' number of ion states (nbrionst) :',i2,/)
1025	format(' The actual input parameter are:'/,
     .    ' length of energy grid (nen) :',i5,/,
     .    ' number of species (nspec) :',i2)
1030	format(' Cross-section input file :',a20/
     .    ' Output files are         :',a20,' and ',a20/)
1040	format(' Testrun without one or more out of the following :'/
     .    ' ionization-cross-sections :',l2,' excitaton-',
     .    'cross-sections :',l2,' , energy degradation :',l2/
     .    ' test with constant scatter-ratios :',l2/)
1050	format( ' Energy grid : Center cell energies [eV] (centE)',
     .    ' cell widths [eV] (engdd)'/)
1060	format(3(f10.2,f10.2,4x))
1070	format(/'The energy grid covers the range from :',f11.3,
     .    ' eV up to :',f8.1,' eV'/)
1080    format('Elastic cross-sections (cel(specie,energy)):',/,
     .    3x,'Energy',6x,5(4x,a4,4x))
1090    format(1x)
1100    format(1pe10.2,2x,5(1pe12.3))
1110	format(/,t20,'Inelastic cross-sections of ',a5,
     .    /,t20,'---------------------------------',
     .    /,'    Energy ','   total   ',5a10)
1120	format( /,'    Energy ',10a11)
1130	format(1pe10.2,10e10.3)
1140	format(/' Threshold in eV :')
1150	format(3(1a11,1f6.2,' | '))
1160	format(/,'                     Cross-section test for ',a5,//,
     .    19x, 'degadation                 seconday production'/,
     .    ' energy |  integral   input  rel. diff.|  ',
     .    'integral   input  rel. diff.',/,
     .    7x,' |                              |')
1170	format(1pe9.2,'|',3(1pe10.2),'|',3(1pe10.2))
1180 	format(/' Redistribution coefficients for ',a5,/)
1190 	format('   E(n)  -->  E(m)     E(n)-E(m)   ',
     .    'Degraded  Ioniz.sec.')
1200 	format(1pe9.2,'-->',1pe9.2,4x,1pe8.1,1x,2(1pe11.2))
1210    format('E # = ',i3,'  emin=',1pe10.2,'   emax=',1pe10.2,
     .    '  Number of species=',i3)
c
c
c       write(6,*)
        write(6,*)'    degrad.f : computation of the cross sections'
c       write(6,*)'    ---------'
c       write(6,*)

c 	Met les energies en ordre croissant
 	    nen = neninit
 	    if(centE(1).lt.centE(nen))then
 	      do ien = 1,nen
 	        e(ien) = centE(ien)
 	        engdd(ien) = ddeng(ien)
 	      enddo
 	    else
 	      do ien = 1,nen
 	       e(nen+1-ien) = centE(ien)
 	       engdd(nen+1-ien) = ddeng(ien)
 	     enddo
 	    endif

c 	Input control parameters
     	open(fic_datdeg,file='dir.data/dir.linux/dir.cine/DATDEG',
     &           status='old')
     	read(fic_datdeg,*)ibid
     	read(fic_datdeg,1000)crsin
     	crsin = crsin(1:lenc(crsin))
	    read(fic_datdeg,1000) crsout
     	crsout = crsout(1:lenc(crsout))
	    read(fic_datdeg,1000) rdtout
     	rdtout = rdtout(1:lenc(rdtout))
     	call xline(4,fic_datdeg)
     	read (fic_datdeg,*) idess(1),idess(2)
     	call xline(2,fic_datdeg)
     	read (fic_datdeg,*) iprint1,iprint2,iprint3,iprint4
     	read(fic_datdeg,1000)logint   ! logint: interpolation type for
c 				             cross sections
     	read(fic_datdeg,1000)lt1           ! lt1   : without ionization
     	read(fic_datdeg,1000)lt2           ! lt2   : only with ioniz.
     	read(fic_datdeg,1000)lt3           ! lt3   : without degradation
     	read(fic_datdeg,1000)lt4           ! lt4   : without excitation
     	read(fic_datdeg,1000)lopal    ! Opal sec. distrib. (Rees sinon)
	    close (fic_datdeg)
c
c	Check for old computation
c
	    open(unfic_crsout_degrad,
     &	     	file='dir.data/dir.linux/'
     &                     //crsout,status='OLD',form='UNFORMATTED',
     &		iostat=iost)
 	    if (iost.ne.0)then
   	    print*,' No old differential cross section file'
 	      print*,' Cross sections therefore computed'
c 	  Puisqu'il n'existe pas, cree le fichier de sortie
	      open(unfic_crsout_degrad,
     .	     	file='dir.data/dir.linux/'
     &                     //crsout,status='new',form='UNFORMATTED',
     .		iostat=iost)
 	      rewind(unfic_crsout_degrad)
  	    close(unfic_crsout_degrad)
 	      go to 10
 	    endif
c
c 	Teste si les sections efficaces ont deja ete calculees sur cette
c 	grille d'energie
c
c 	On teste si le fichier contient quelque chose
c 	end=20 renvoie a l'etiquette 20 en cas de non lecture
	    read(unfic_crsout_degrad,end = 20) ien,isp,i1,i2

c
     	if(ien.ne.nen .or. isp.lt.nspec)then
     	  if(kiappel.eq.1)then
       	    print*,'Old cross section file does not match new energy'
     	    print*,'grid requirement for nen and nspec'
     	    print*,'Cross sections therefore computed'
     	  endif
      	  close(unfic_crsout_degrad)
     	  go to 10
     	endif
c
	    read(unfic_crsout_degrad) (esig(ien),ien=1,nen)
	    read(unfic_crsout_degrad) (bsig(ien),ien=1,nen)
	    read(unfic_crsout_degrad) (delta(ien),ien=1,nen)
c
     	do ien = 1,nen
     	  if (e(ien).ne.esig(ien) .or. ddeng(ien).ne.delta(ien)) then
     	    if(kiappel.eq.1)then
       	      print*,'Old cross section file does not match new energy'
     .		,' grid requirement'
     	      print*,'Cross sections therefore computed'
     	    endif
      	    close(unfic_crsout_degrad)
     	    go to 10
     	  endif
     	enddo
c
c 	On teste si les sections efficaces sont vraiment calculees
c 	end=20 renvoie a l'etiquette 20 en cas de non lecture
        read(unfic_crsout_degrad,end=20)(nexcst(iisp),iisp=1,nspec),
     .		(nionst(iisp),iisp=1,nspec)

c
     	if(kiappel.eq.1)then
       	  print*,' Old cross section file matches new energy',
     . 	       ' grid requirement'
         	  print*,'No new cross sections computed'
         	  print*
     	endif
c
      	close(unfic_crsout_degrad)
     	return 		! cross sections already computed
c
20 	print*,'Old cross section file not complete'
 	print*,'Cross sections therefore computed'
c 	Puisque ce fichier contient quelque chose, on le rembobine
 	rewind(unfic_crsout_degrad)

10 	continue	! cross sections not yet computed
c
c----  	Open input files : seff file
	open(fic_crsin_degrad,file='dir.data/dir.linux/'
     &                                   //crsin,status='old')
 	rewind fic_crsin_degrad
c
c----  	Open formatted output files
 	open(fic_degout,file= 'dir.data/dir.linux/'
     &                              //'dir.cine/DEGOUT',
     .			status='unknown')
 	rewind fic_degout

c 	Open unformatted output cross section files
	open(unfic_crsout_degrad,
     .		file='dir.data/dir.linux/'
     &                     //crsout,status='old',form='UNFORMATTED')
 	rewind unfic_crsout_degrad
c
	write(fic_degout,1010) istdate,isttime
	write(fic_degout,1020) nbren,nbrango2,nbrsp,nbrexc,nbrionst
	write(fic_degout,1030) crsin,crsout,rdtout
	if(lt1.or.lt2.or.lt3.or.lt4)
     .		write(fic_degout,1040) lt1,lt2,lt3,lt4

c----   Read cross section file, generate energy grid and interpolate
	call sigma(nspec,logint,title,excion,cel,cin,cinex,
     .  	         ethres,bratio,nexcst,nionst,e,nen)
	write(6,1210)nen,e(1),e(nen),nspec

	write(unfic_crsout_degrad) nen,nspec,nbrexc,nbrionst
	write(unfic_crsout_degrad) (centE(n),n=1,nen)
	write(unfic_crsout_degrad) (botE(n),n=1,nen)
	write(unfic_crsout_degrad) (engdd(n),n=1,nen)
	write(fic_degout,1025) nen,nspec
	write(fic_degout,1050)
	write(fic_degout,1060) (centE(n),engdd(n),n=1,nen)
	e1=max(e(1)-engdd(1)/2.,0.)
	e2=e(nen)+engdd(nen)/2.
	write(fic_degout,1070) e1,e2

	if(iprint1.eq.1) then
	  write(fic_degout,1080) (specie(isp),isp=1,nspec)
	  write(fic_degout,1090)
	  do ien=1,nen
            write(fic_degout,1100) e(ien),(cel(isp,ien),isp=1,nspec)
 	  enddo
	end if
	print*,' Elastic cross-sections done'

	lfstr=.true.			! O-fine-structure
c
c----	Inelastic cross-sections
	call crosin(lfstr,lt1,lt2,title,nspec)
	print*,' Inelastic cross-sections done'
	if(iprint2.eq.1) then
	  do isp=1,nspec
	    write(fic_degout,1110) specie(isp),(title(i,isp),i=1,5)
	    do ien=1,nen
	      write(fic_degout,1130)
     .		e(ien),cin(isp,ien),(cinex(isp,js,ien),js=1,5)
 	    enddo
 	    if(nexcst(isp).gt.5)then
	      write(fic_degout,1120) (title(i,isp),i= 6,nexcst(isp))
	      do ien=1,nen
	        write(fic_degout,1130)e(ien),(cinex(isp,js,ien),
     .				js=6,nexcst(isp))
 	      enddo
 	    endif
   	    write(fic_degout,1140)
   	    write(fic_degout,1150)
     .		(title(i,isp),ethres(isp,i,1),i=1,nexcst(isp))
 	  enddo
	end if
	if(lt3) then
	  etest=15.
	  do isp=1,nbrsp
	    do iexc=1,nbrexc
	      ethres(isp,iexc,1)=e(nen)
	      ethres(isp,iexc,2)=e(nen)
	      ethres(isp,iexc,3)=e(nen)
	      ethres(isp,iexc,4)=e(nen)
	      ethres(isp,iexc,5)=e(nen)
	      do ien=1,nen
	        if(isp.eq.1.and.iexc.eq.9.and.e(ien).ge.etest) then
		  cinex(isp,iexc,ien)=1.
		  cin(isp,ien)=1.
		else
		  cinex(isp,iexc,ien)=0.
		  cin(isp,ien)=0.
		end if
 	      enddo
 	    enddo
 	  enddo
	  ethres(1,9,1)=etest
	end if
c
c----	degradation cross-sections and adjustment with the Swartz method
	line='Secondary distribution from Rees et al.,'//
     .		' PSS 17, 1997, 1969'
	if(lopal) line='Secondary distribution from Opal et al.,'//
     .		' J. Chem. Phys., 55, 4100, 1971'
	do isp=1,nspec
	  call redist(lt4,isp,kiappel)
 	enddo
	print*,' Degradation cross-sections done'
c	print*,line

c----	Test for the differential cross-sections

	do isp=1,nspec
	  if(iprint4.eq.1) write(fic_degout,1160) specie(isp)
	  do ien=1,nen
   	    tomdeg=0.
	    tomsec=0.
	    rsec=0.
	    rdeg=0.
	    do iens=1,ien
	      tomdeg=tomdeg+omdeg(ien,iens,isp)*engdd(iens)
   	      tomsec=tomsec+omdeg(iens,ien,isp)*engdd(iens)
 	    enddo
	    if(cin(isp,ien).ne.0.)
     .		rdeg=(tomdeg-cin(isp,ien))/cin(isp,ien)
	    if(cinex(isp,nexcst(isp),ien).ne.0.)
     .	        rsec= (tomsec-cinex(isp,nexcst(isp),ien))
     .		      / cinex(isp,nexcst(isp),ien)
	    if(iprint4.eq.1) write(fic_degout,1170) e(ien),tomdeg,
     .		cin(isp,ien),rdeg,tomsec,cinex(isp,nexcst(isp),ien),rsec
 	  enddo
  	enddo

	write(unfic_crsout_degrad)
     .		(nexcst(isp),isp=1,nspec),(nionst(isp),isp=1,nspec)
	write(unfic_crsout_degrad)
     .		((title(js,isp),isp=1,nspec),js=1,nbrexc)
	write(unfic_crsout_degrad)
     .	   (((ethres(isp,js,jp),isp=1,nspec),js=1,nbrexc),jp=1,nbrionst)
	write(unfic_crsout_degrad)
     .		((bratio(jp,isp),jp=1,nbrionst),isp=1,nspec)
	write(unfic_crsout_degrad) ((cel(j,n),j=1,nspec),n=1,nen)
	write(unfic_crsout_degrad) ((cin(j,n),j=1,nspec),n=1,nen)
	write(unfic_crsout_degrad)
     .		(((cinex(j,js,n),j=1,nspec),js=1,nbrexc),n=1,nen)
	write(unfic_crsout_degrad)
     .		(((cinex(j,js,n),j=1,nspec),js=1,nbrexc),n=1,nen)
	close (unfic_crsout_degrad)
c
c 	Open unformatted output differential cross section files
	open(unfic_rdtout_degrad,file='dir.data/dir.linux/'
     &                                      //rdtout,form='UNFORMATTED')
 	rewind unfic_rdtout_degrad
	write(unfic_rdtout_degrad) line
	write(unfic_rdtout_degrad) nen,nspec
c
	do nfix=nen,1,-1
	  nend=nfix
	  write(unfic_rdtout_degrad)
     .		((omdeg(n,nfix,isp),n=nen,nend,-1),isp=1,nspec)
	  write(unfic_rdtout_degrad)
     .		((omdeg(nfix,n,isp),n=nen,nend,-1),isp=1,nspec)
   	enddo
	close (unfic_rdtout_degrad)
c
 	if (iprint3.eq.1)then
 	  do isp = 1,nspec
 	    write(fic_degout,1180) specie(isp)
c 	    On va de n a nfix
 	    nfix = 1
c 	    do nfix=nen-1,1,-1
 	      write(fic_degout,1190)
 	      do n=nen,nfix,-1
 	        write(fic_degout,1200)e(n),e(nfix),e(n)-e(nfix),
     .  	    omdeg(n,nfix,isp),omdeg(nfix,n,isp)
 	      enddo
c	    enddo
 	  enddo
 	endif
c
	close (fic_degout)
c
 	return
	end

*-----------------------------------------------------------------------

        pure real function terplin(xs,f,mxloi,xval)
				implicit none
*	TERPLIN performs a linear interpolation
*	of the function F(XS). MXLOI is the length
*	of the arrays F and XS, XVAL is the
*	X-value to where the interpolation is required
*
			  integer, intent(in) :: mxloi
      	real,intent(in) :: xs(mxloi),f(mxloi),xval
				integer n

      	if(xs(1).gt.xs(mxloi)) goto 1
      	do n=2,mxloi
          if(xval.gt.xs(n-1).and.xval.le.xs(n)) go to 3
 	      end do
        if(xval.ge.xs(mxloi)) goto 3

1       continue
        do n=mxloi,2,-1
          if(xval.lt.xs(n-1).and.xval.ge.xs(n)) goto 3
 	      end do

3       continue
        terplin=f(n-1)+ (xval - xs(n-1))*(f(n) - f(n-1))/(xs(n)-xs(n-1))

        end function terplin
c
*-----------------------------------------------------------------------
c
	real function avf(n,np,w,emax)
*	calculates averages over energy

	include 'TRANSPORT.INC'
	common /eng_/ e(nbren),engdd(nbren),nen

*	secondary production  (ionization collision)
	dde=engdd(n)/2.
	e1=e(n)-dde
	e2=min(e(n)+dde,emax)
	if(e2.le.e1) then
	  avf=0.
	  return
	end if
	avf=fnorm(e1,e2,e(np),w)/(e2-e1)
	return
*
C       primary degradation (ionization collision)
	entry avg(n,np,w,emax)
C       (this is F backwards)
	dde=engdd(n)/2.
	e1=e(n)-dde
	e2=min(e(n)+dde,emax)
	if(e2.le.e1) then
	  avg=0.
	  return
	end if
	avg=fnorm(e(np)-e2-w,e(np)-e1-w,e(np),w)/(e2-e1)

	end function avf

*-----------------------------------------------------------------------

	subroutine crosin(lfstr,lt1,lt2,title,jspec)

C	Generates inelastic cross-sections by reading a	data file

***********************  Single precision  ****************************

	implicit logical (l)
	include 'TRANSPORT.INC'
c
        real e(nbren),engdd(nbren)
        integer nen
        common /eng_/ e,engdd,nen
c
        real ethres(nbrsp,nbrexc,nbrionst),bratio(nbrionst,nbrsp)
        integer nexcst(nbrsp),nionst(nbrsp)
        common /eng2_/ ethres,bratio,nexcst,nionst
c
        real cel(nbrsp,nbren),cin(nbrsp,nbren),cinex(nbrsp,nbrexc,nbren)
        common /cros_/ cel,cin,cinex
c
	character*9 title(nbrexc,nbrsp)
C       Thresholds for different states of the produced ion
c 	are ethresh(nbrsp,nbrexc,nbrionst),bratio(nbrionst,nbrsp)
c
	do isp=1,jspec
	  do n=1,nen
	    if(12.0.le.e(n).and.e(n).le.23.1) then
	      crosp=0.
	      do ionst=1,nionst(isp)
	        if(e(n).ge.ethres(isp,nexcst(isp),ionst))crosp=crosp+
     .		            bratio(ionst,isp)*cinex(isp,nexcst(isp),n)
 	      enddo
	      cinex(isp,nexcst(isp),n)=crosp
	    end if
	    if(lt1.or.lt2) then   ! test with one cross section only !!!
	      do js=1,nexcst(isp)
	        jstest=nexcst(isp)
	        if(lt1.and.js.eq.jstest) cinex(isp,js,n)=0.
	        if(lt2.and.js.ne.jstest) cinex(isp,js,n)=0.
 	      enddo
	    end if
 	  enddo
c
 	enddo

C	inelastic cross-section is the sum over all states
	if(.not.lfstr.or.lt2) goto 506
C	include O-fine-structure in inelastic cross-sections
C	move all other cross-sections one index up

	do js=nexcst(3),1,-1
 	  title(js+1,3)=title(js,3)
	  do jp=1,nbrionst
   	    ethres(3,js+1,jp)=ethres(3,js,jp)
 	  enddo
	  do n=1,nen
   	    cinex(3,js+1,n)=cinex(3,js,n)
 	  enddo
 	enddo
	nexcst(3)=nexcst(3)+1
	ethres(3,1,1)=1.78e-02
 	title(1,3) = 'fne struc'
	do n=1,nen
	if(e(n).gt.5.) then
	  cinex(3,1,n)=0.
	else
	  cinex(3,1,n)=finestr(e(n))
	end if
 	enddo
506	continue

	do isp=1,jspec
	  do n=1,nen
	    cin(isp,n)=0.
	    do js=1,nexcst(isp)
   	      cin(isp,n)=cin(isp,n)+cinex(isp,js,n)
 	    enddo
 	  enddo
   	enddo

	return
	end

*-----------------------------------------------------------------------

	subroutine redist(ltst,isp,kiappel)

*	This subroutine uses the formula given in Rees, Steward and
*	Walker, PSS, 17, 1997-2008, 1969 for calculating the cross-
* 	sections for the production of secondary electrons in
*	ionization processes
*	Note that these cross-sections have the units cm**2/eV
*	normalized so that the integral over all secondary energies
*	reproduces the original cross-section (ionization as well as
*	exitation). The method used to calculate the degradation
*	cross-sections differs slightly from the Swartz method.
* 	The Rees-formula is used directly instead.
*	If LTST=.TRUE. do only degradation due to ionization.
*	If LSEC=.TRUE., the secondary production cross section is not
*		Swartz adjusted.

******************** Single precision **********************************

	implicit logical (l)
	include 'TRANSPORT.INC'
c
        real e(nbren),engdd(nbren)
        integer nen
        common /eng_/ e,engdd,nen
c
        real ethres(nbrsp,nbrexc,nbrionst),bratio(nbrionst,nbrsp)
        integer nexcst(nbrsp),nionst(nbrsp)
        common /eng2_/ ethres,bratio,nexcst,nionst
c
        real cel(nbrsp,nbren),cin(nbrsp,nbren),cinex(nbrsp,nbrexc,nbren)
        common /cros_/ cel,cin,cinex
c
	common /rdist_/ omdeg(nbren,nbren,nbrsp)

!	First treat ionization: loss of primaries and production
!	of secondaries
c
4790 	format('   E(n)  -->  E(m)     E(n)-E(m)   ',
     .    'Degraded  Ioniz.sec.')
4700 	format(1pe9.2,'-->',1pe9.2,4x,1pe8.1,1x,2(1pe11.2))
c
	  lsec=.true.
        if (kiappel.eq.1)
     .  write(6,*)'REDIST started with specie',isp,'            [A'
c       write(6,*)'REDIST started with specie',isp
c 	On part de l'energie la plus elevee (nen) vers la plus basse.
c	We start from the highest energy (nen) to the lowest.
	  do np=nen,2,-1
	    dde=engdd(np)
	   crosp=0.
	   if(lsec) crosec=0.
c 	  On examine tous les ions excites
c	  All excited ions were examined
	   do jp=1,nionst(isp)
c 	    On regarde le seuil pour l'ionisation excitative.
c	    We look at the ionization threshold for the excitative
	    seuil=ethres(isp,nexcst(isp),jp)
c 	    Ce seuil doit etre plus petit que l'energie de depart !
c 	    et de meme, l'energie de depart moins le seuil doit etre
c 	    plus grande que l'energie la plus basse
c	    This threshold should be smaller than the energy of departure! and similarly, the energy of departure less the threshold must be higher that the lowest energy.
 	    if(seuil.lt.e(np) .and. e(np)-seuil .ge.e(1)) then
c             sink energy < than source
	      nk=nlev(e(np)-seuil,np-1,isp,1)
c
	      if(e(np)-seuil.lt.e(nk)) then  ! take nearest cell bound.
	        ek=e(nk)-engdd(nk)/2.
	      else
	        ek=e(nk)+engdd(nk)/2.
	      end if
c
c 	      The next line might be in  c mode, BUT
c 	      When the next line is active, the results look better :
c 	      the wiggle around the energy where the cell size equals
c 	      the energy loss goes almost away.
 	      ek=e(np)-max(engdd(np)/2.,seuil)
	      ek=max(0.,ek)
	      es=max(0.,e(np)-seuil)
	      wi=e(np)-ek         ! de facto min. energy loss of primary
c 	      adjust if necessary
	      cros=cinex(isp,nexcst(isp),np)*bratio(jp,isp)*seuil/wi
	      crosp=crosp+cros
	      if(lsec) crosec=cinex(isp,nexcst(isp),np)*bratio(jp,isp)
	      fac=fnorm(0.,ek,e(np),wi)
	      if(lsec) facs=fnorm(0.,es,e(np),seuil)
	      ffac=0.
	      ffacs=0.
*	      cros=1.				! test
*	      crosec=1.				! test
	      do n=1,nk,1
	        ffac=ffac+avg(n,np,wi,ek)*engdd(n)
	        ffacs=ffacs+avf(n,np,seuil,es)*engdd(n)
 	      enddo
	      fac=ffac
	      facs=ffacs
	      if(fac.eq.0.) then
 	         write(6,1000)np,e(np),nk,e(nk),jp,seuil
1000 	  	 format('Warning ! Energy step too big from e(',
     .        i3,') = ',1pe10.2,' to e(',i3,') = ',1pe10.2,/,
     .       'for inelastic phenomenum # ',i2,' (Threshold = )',
     .        0pf10.2)
c 	         print*,'NP,E(NP),WI,NK,EK',np,e(np),wi,nk,ek
 	      endif
	      if(fac.ne.0.) fac=cros/fac
	      if(facs.ne.0..and.lsec) facs=crosec/facs
	      do ns=1,nk		   ! NS : sink energy cell
	        if(lsec) then
c 	          Secondaries :
	          omdeg(ns,np,isp)=omdeg(ns,np,isp)+
     .				  avf(ns,np,seuil,es)*facs
	        else
c 	          Secondaries :
	          omdeg(ns,np,isp)=omdeg(ns,np,isp)+avf(ns,np,wi,ek)*fac
	        end if
c 	        Primaries :
	        omdeg(np,ns,isp)=omdeg(np,ns,isp)+avg(ns,np,wi,ek)*fac
 	      enddo
 	    endif
 	  enddo
c
*	  "adjusted" ionization cross-section
	  cinex(isp,nexcst(isp),np)=crosp
	  cin(isp,np)=crosp
 	 enddo

	 if(ltst) return

*	Degrade primaries due to exitation. For reference see Swartz
*	"Optimization of discrete energy degradation", JGR 1985, p.6587
c
	jsend=nexcst(isp)-1		   ! don't do ionization again
c
c 	On part cette fois de l'energie la plus basse
c 	The smallest energy is taken below usual thermal energy
c 	(0.01eV = 100 K)
 	einit = max(e(1) - engdd(1)/2.,0.01)
	do n=2,nen
  	  dde=engdd(n)			   ! source cellwidth
	  dde2=dde/2.
 	  do js=1,jsend
	    seuil=ethres(isp,js,1)	   ! energyloss
 	    ew=max(seuil,dde)
 	    emin=e(n)-dde2-ew		   ! spread primaries over sink
 	    emax=e(n)+dde2-ew		   ! cells so that sum of sink
c
  	    if(emin .gt. einit) then
c
	      nk=nlev(e(n)-seuil,n-1,isp,2)
c
	      ntop=nlev(emax,n-1,isp,3)
c
	      nbot=nlev(emin,ntop,isp,4)
c
	      nbot=max(min(nbot,ntop-1),1)
	      cros=cinex(isp,js,n)*seuil/ew
	      do i=nbot+1,ntop-1,1
   	        omdeg(n,i,isp)=omdeg(n,i,isp)+cros/dde
 	      enddo
         fac=(emax-e(ntop)+engdd(ntop)/2.)/engdd(ntop)
	      omdeg(n,ntop,isp)=omdeg(n,ntop,isp)+cros*fac/dde
	      fac=(-emin+e(nbot)+engdd(nbot)/2.)/engdd(nbot)
	      omdeg(n,nbot,isp)=omdeg(n,nbot,isp)+cros*fac/dde
	      cinex(isp,js,n)=cros        ! restore modified cross-sect.
	      cin(isp,n)=cin(isp,n)+cinex(isp,js,n)
 	    endif
 	  enddo
 	enddo
c
c ----	Fin de modif
	return
	end

*-----------------------------------------------------------------------

	function finestr(e)
*
*	This function returns the cross-section for the oxygene 3P2,1,0
*	fine structure excitation. The collision strengths are taken
*	from Dourneuf and Nesbet, J.Phys.B,Vol.9,No.9,L241,1976
*	(formula 2)
* 	The cross-section is then calculated using Hoegy,GRL,Vol.8,541,
*	1976 (formula 1). The cross-sections for the several levels is
*	then added to give a total cross-section
*
**************************  Single precision  **************************

	g2=5.
	g1=3.
	pia=8.7974e-17
	om21=0.
	om20=0.
	om10=0.
	t=e/8.617e-5
	if(t.gt.228.) om21=(9.76e-6+3.46e-11*(t-228.))*(t-228.)
	if(t.gt.326.) om20=(3.29e-6-2.90e-11*(t-326.))*(t-326.)
	if(t.gt.98. ) om10=(1.89e-6+8.00e-11*(t-98. ))*(t-98. )
	sqk=t*6.3335e-6
	sig21=pia/g2/sqk*om21
	sig20=pia/g2/sqk*om20
	sig10=pia/g1/sqk*om10
	finestr=sig21+sig20+sig10
	return
	end

*-----------------------------------------------------------------------

	integer function nlev(eng,ntop,isp,imod)

*	Function to find the grid number to which a given energy belongs
*	input 	ENG : energy 	,	ntop : guess for an upper limit
*	output 	NLEV : cell that contains energy ENG
c 	imod = 1 pour ionisation, 2 pour excitation (juste pour info).

*************************  Single precision  ***************************

	include 'TRANSPORT.INC'
	common /eng_/e(nbren),engdd(nbren),nen
c
	if(eng.le.0..or.ntop.le.1) then
	  nlev=1
	  return
	end if
c
	do n=ntop,1,-1
	  de=engdd(n)/2.
	  if(eng.ge.e(n)-de) then
	    nlev=n
	    return
	  end if
 	enddo
c
 	write(6,*)
 	if(imod.eq.1)then
	  write(6,1000)specie(isp),imod
 	  write(6,1020)ntop+1,e(ntop+1),e(ntop+1)-eng
 	else if(imod.eq.2)then
 	  write(6,1010)specie(isp),imod
 	  write(6,1020)ntop+1,e(ntop+1),e(ntop+1)-eng
 	else if(imod.eq.3)then
 	  write(6,1010)specie(isp),imod
 	  write(6,1030)ntop+1,e(ntop+1),e(ntop+1)-eng
 	else if(imod.eq.4)then
 	  write(6,1010)specie(isp),imod
 	  write(6,1040)eng,ntop
 	endif
	write(6,*)
1000 	format('Error stop in NLEV for specie ',a4,
     .    ' in ionization redistribution, imod =',i2)
1010 	format('Error stop in NLEV for specie ',a4,
     .    ' in excitation redistribution, imod =',i2)
1020 	format('e(',i3,') = ',1pe10.2,' eV. Seuil =',1pe10.2,' eV')
1030 	format('e(',i3,') = ',1pe10.2,' eV. ',
     .    'Max(seuil,ddeng(n)/2) - ddeng(n)/2 =',1pe10.2,' eV')
1040 	format('emin = ',1pe10.2,' eV. ',
     .    'Estimation impossible a partir de l''energie # ',i3)
        stop 'arret dans nlev'
c
	end

*-----------------------------------------------------------------------

	real function fnorm(a,b,yy,z)

*	Integration with Takahashi's method to find normalization
*	for secondary electron distribution function

	real ww(-17:17),xx(-17:17)
	data ww
     .	/1.36435219e-07, 1.64996788e-06, 1.40066322e-05, 8.77183047e-05,
     .	 4.23141464e-04, 1.63230451e-03, 5.20187337e-03, 1.40839620e-02,
     .	 3.31747346e-02, 6.93263561e-02, 1.30518615e-01, 2.23881721e-01,
     .	 3.52479190e-01, 5.11343718e-01, 6.84302330e-01, 8.44149709e-01,
     .	 9.58392859e-01, 1.00000000e+00, 9.58392859e-01, 8.44149709e-01,
     .	 6.84302330e-01, 5.11343718e-01, 3.52479190e-01, 2.23881721e-01,
     .	 1.30518615e-01, 6.93263561e-02, 3.31747346e-02, 1.40839620e-02,
     .	 5.20187337e-03, 1.63230451e-03, 4.23141464e-04, 8.77183047e-05,
     .	 1.40066322e-05, 1.64996788e-06, 1.36435219e-07/
	data xx
     . /-1.00000000e+00,-9.99999821e-01,-9.99998450e-01,-9.99988973e-01,
     .	-9.99938786e-01,-9.99728441e-01,-9.99006689e-01,-9.96921241e-01,
     .	-9.91719127e-01,-9.80284095e-01,-9.57760513e-01,-9.17497337e-01,
     .	-8.51593614e-01,-7.52291977e-01,-6.14237905e-01,-4.37139213e-01,
     .	-2.27766961e-01, 0.00000000e+00, 2.27766961e-01, 4.37139213e-01,
     .	 6.14237905e-01, 7.52291977e-01, 8.51593614e-01, 9.17497337e-01,
     .	 9.57760513e-01, 9.80284095e-01, 9.91719127e-01, 9.96921241e-01,
     .	 9.99006689e-01, 9.99728441e-01, 9.99938786e-01, 9.99988973e-01,
     .	 9.99998450e-01, 9.99999821e-01, 1.00000000e+00/
        data h/ 2.30999470e-01/zero/0./

	t(x)=((bb-aa)*x + aa+bb)/2.	! transformation [a,b] to [-1,1]
	data nw/0/
c
	if(a.ge.b) then
	  print*,'WARNING FROM FNORM',a,b,yy,z
	  nw=nw+1
	  if(nw.gt.10) stop 'arret dans fnorm'
	  fnorm=0.
	  return
	end if
c
	ntop=17
	delmax=100.
	fnorm=0.
	nn=(b-a)/delmax
	d=(b-a)/(nn+1)
	do i=1,nn+1,1
	  aa=a+(i-1)*d
	  bb=a+i*d
	  fn=0.
	  do n=-ntop,ntop,1
	    fn=fn + f(t(xx(n)),yy,z)*ww(n)
 	  enddo
	  fn=d/2.*h*fn
	  fnorm=fnorm+fn
 	enddo

	end function fnorm

*-----------------------------------------------------------------------

      real function f(x,y,z)
*	secondary electron distribution function
*	select from various sources (Rees et al, Opal et al, etc)
			logical lopal
			common /distfct_/lopal
			data zero/0./elimit/5600./
			if(x.le.0..or.x.ge.y-z) then
			  f=0.
			  return
			end if
			if(lopal) then					! Opal et al
			  f=1.+(x/18.)**2.1
			  f=1./f
			  if(x.gt.(y-z)/2.) f=0.
			  return
			else						! Rees et al
			  if((x+z)/2.49.le.elimit) then
			      f=exp(-(x+z)/31.5-339.*exp(-(x+z)/2.49))/(x+z)
     .	            *log((sqrt(y)+sqrt(max(zero,y-x-z)))
     .	              /(sqrt(y)-sqrt(max(zero,y-x-z))))
			  else
			    f=0.
			  end if
			end if

      end function f
c
c---------------------- sigma ---------------------------------------
c
      subroutine sigma(nspec,logint,title,excion,cel,cin,cinex,
     .  	         ethres,bratio,nexcst,nionst,e,nen)
c
c	This program reads the datafile with the cross-section data
c	and converts it to the format required by the electron transport
c	code.
c	The variable CRSIN contains the name of this input data file.
c
         implicit none
	     include 'TRANSPORT.INC'

c 	INPUTS
 	real e(nbren)
 	integer nen,nspec
 	logical logint
c
c 	OUTPUTS
	real ethres(nbrsp,nbrexc,nbrionst),bratio(nbrionst,nbrsp)
 	real e1min,e1max
        integer	nexcst(nbrsp),nionst(nbrsp)
	real cel(nbrsp,nbren),cin(nbrsp,nbren),cinex(nbrsp,nbrexc,nbren)
	character*9 title(nbrexc,nbrsp),excion(nbrionst,nbrsp)
c
c 	INTERNAL
	real e1(nbren),work(nbren),e1log(nbren),elog(nbren)
 	real xout(1),yout(1)
	real celmul(nbrsp),cinelmul(nbrsp,nbrexc),terplin,w1,w2,fac
 	integer ien,nnmax,isp,n400,n2200,ist,ionst
c
1000 	format(a)
1010 	format(10(a9))
c
  	do ien=1,nen
 	  elog(ien)=log(e(ien))
 	enddo
c
 	do isp = 1,nspec
c
c 	  Lecture des energies pour les seff elastiques
 	  call xline(3,fic_crsin_degrad)
	  read(fic_crsin_degrad,*) nnmax
	  read(fic_crsin_degrad,*) (e1(ien),ien=1,nnmax)
          call mnmx(e1,nnmax,e1min,e1max,0)
 	  do ien=1,nnmax
 	    e1log(ien)=log(e1(ien))
 	  enddo
c	  Read input elastic cross-sections and interpolate to the
c		output energy-grid
	  read(fic_crsin_degrad,*) celmul(isp)
	  read(fic_crsin_degrad,*) (work(ien),ien=1,nnmax)
c 	  Multiplication with fudge factors.
 	  do ien = 1,nnmax
 	    work(ien) = work(ien)*celmul(isp)
	    if(logint) work(ien)=log(work(ien)*celmul(isp))
 	  enddo
 	  w1 = work(1)
 	  w2 = work(2)
	  do ien=1,nen
            if(e(ien).lt.e1min) then
              call intlin(nnmax,e1log,work,1,log(e(ien)),yout)
              cel(isp,ien) = yout(1)
            elseif(e(ien).le.e1max) then
	      cel(isp,ien)=terplin(e1,work,nnmax,e(ien))
	      if(logint) cel(isp,ien)=exp(cel(isp,ien))
	    end if
 	  enddo
c
c ----    Changes from Dirk, 05/05/1992, de-bugged (jl, sept 92)
*	  elastic cross section above 400 eV are proportional to E^-0.65
 	  if(e(nen).ge.400.)then
	    do ien=1,nen
	      if(e(ien).ge.400.) then
	        n400=ien
	        goto 400
	      endif
	    enddo
 	  endif
400	  continue
 	  if(e(nen).ge.2200.)then
	    do ien=1,nen
	      if(e(ien).ge.2200.) then
	        n2200=ien
	        goto 401
	      endif
	    enddo
 	  else
 	    n2200 = nen
 	  endif
401	  continue
c
 	  if(e(nen).ge.400.)then
 	    fac=e(n400)**0.65
	    do ien=n400,n2200
	      cel(1,ien)=cel(1,n400)*fac/e(ien)**0.65
	      cel(2,ien)=cel(2,n400)*fac/e(ien)**0.65
	      cel(3,ien)=cel(2,ien)*0.5
	    enddo
	  endif
c
*	  elastic cross sections above 2200 eV are proportional to E^-1
 	  if(e(nen).ge.2200.)then
   	    fac=e(n2200)
	    do ien=n2200,nen
	      cel(1,ien)=cel(1,n2200)*fac/e(ien)
	      cel(2,ien)=cel(2,n2200)*fac/e(ien)
	      cel(3,ien)=cel(2,ien)*0.5
	    enddo
 	  endif
c ---- 	  End of changes 05/05/92
c
c	  Read input inelastic cross-sections, interpolate, extrapolate
c 	  at high energies
 	  call xline(1,fic_crsin_degrad)
	  read(fic_crsin_degrad,*) nexcst(isp)		!  states #
 	  do ist = 1,nexcst(isp)
 	    read (fic_crsin_degrad,1000)title(ist,isp)
	    read(fic_crsin_degrad,*) ethres(isp,ist,1)  ! exc. thresh.
	    read(fic_crsin_degrad,*) nnmax
 	    call xline(1,fic_crsin_degrad)
	    read(fic_crsin_degrad,*) (e1(ien),ien=1,nnmax)
 	    do ien=1,nnmax
 	      e1log(ien)=log(e1(ien))
 	    enddo
	    read(fic_crsin_degrad,*) cinelmul(isp,ist)
 	    call xline(1,fic_crsin_degrad)
	    read(fic_crsin_degrad,*) (work(ien),ien=1,nnmax)
c 	    Multiplication with fudge factor.
 	    do ien =1,nnmax
 	      work(ien) = work(ien)*cinelmul(isp,ist)
	      if(logint) work(ien)=log(work(ien)*cinelmul(isp,ist))
 	    enddo
c
c 	    Interpolation sur la grille d'energie e(ien)
	    do ien=1,nen
 	      if(e(ien).lt.ethres(isp,ist,1)) then
		cinex(isp,ist,ien)=0.
 	      elseif(e(ien).ge.ethres(isp,ist,1).and.
     .		     e(ien).lt.e1(1)) then
c 		extrapolation.
 	        xout(1)=elog(ien)
 		call intlin(nnmax,e1log,work,1,xout,yout)
 		cinex(isp,ist,ien)=yout(1)
		if(logint) cinex(isp,ist,ien)=exp(cinex(isp,ist,ien))
	      elseif(e(ien).ge.e1(1) .and. e(ien).le.e1(nnmax)) then
		cinex(isp,ist,ien)=terplin(e1,work,nnmax,e(ien))
		if(logint) cinex(isp,ist,ien)=exp(cinex(isp,ist,ien))
 	      elseif(e(ien).gt.e1(nnmax)) then
c 		extrapolation.
 	        xout(1)=elog(ien)
 		call intlin(nnmax,e1log,work,1,xout,yout)
 		cinex(isp,ist,ien)=yout(1)
		if(logint) cinex(isp,ist,ien)=exp(cinex(isp,ist,ien))
	      end if
   	    enddo		! fin boucle sur energies
   	  enddo			! fin boucle sur etats inelastiques.
c 	  Read # of ion states (nionst), thresholds and branching ratios
 	  call xline(1,fic_crsin_degrad)
	  read(fic_crsin_degrad,*) nionst(isp)		! ion state #
c 	  ion excited states
	  read(fic_crsin_degrad,1010)(excion(ist,isp),ist=1,nionst(isp))
c 	  ionisation thresholds
c 	  L'etat ionise est l'etat inelastique nexcst(isp)
	  read(fic_crsin_degrad,*)
     .		(ethres(isp,nexcst(isp),ionst),ionst=1,nionst(isp))
c 	  branching ratio
	  read(fic_crsin_degrad,*)
     .		(bratio(ionst,isp),ionst=1,nionst(isp))

 	enddo			! fin boucle sur especes.
c
	close (fic_crsin_degrad)

	end subroutine sigma
c
c------------------------------------------------------------------
c
