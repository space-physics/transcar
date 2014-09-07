c
c----------------------------------------------------------------------
c
 	subroutine solarflux (iout,yyddd,nen,e,Pflx,Eflx,impress)
c
c 	Cette subroutine calcule un flux solaire a une date donnee 
c 	sur une grille d'energie e.
c 	La source est :
c     	Tobiska, W.K., Revised solar extreme ultraviolet flux 
c       model, JATP, (accepted) 1991. The preprint was XXVIII 
c       COSPAR, The Hague, paper STP I.2-P1, 1990.
c 	Le modele Tobiska est appele, puis les energies sont re-
c 	ordonnees et le flux est interpole sur la grille d'entree.
c
c 	INPUTS
c 	iout  = numero logique du fichier de sortie (integer)
c 	impress = 1 si on imprime le flux , 0 sinon
c 	yyddd = date (integer) between 68172 - 88366: '
c 	nen   = nombre d'energies de la grille d'entree
c 	e     = grille d'entree des energies [eV]
c
c	OUTPUTS
c 	Pflx  = flux de photons sur la grille d'entree [cm-2.s-1]
c 	Eflx  = flux d'energie sur la grille d'entree [eV.cm-2.s-1]
c
 	implicit none
c
 	integer iout,yyddd,numday,nen,ien,impress
c
C 	Sorties de euv91
 	real w1nm(39),w2nm(39),wmid(39)
 	real eVinf(39),eVsup(39)
 	real eVmid(39),Eflux(39),Pflux(39)
c
c 	Sorties reordonnees
 	real w1nmor(39),w2nmor(39),wmidor(39)
 	real eVinfor(39),eVsupor(39)
 	real eVmidor(39),Efluxor(39),Pfluxor(39)
c
c 	Energies d'entree et longueur d'onde correspondantes.
 	real e(400),ondenm(400)
c 	Flux interpolle sur les energies.
 	real PflxE(400),EflxE(400)
c 	Flux interpolle sur les longueurs d'onde.
 	real PflxW(400),EflxW(400)
c 	Facteurs de normalisation sur le nombre total de photons.
 	real totin,totout,rapp
c 	Flux renvoye par le programme
 	real Pflx(400),Eflx(400)

c 	Tableau de travail.
 	integer index(39)
 	double precision eVmiddp(39)
 	real trav(400),tabx(400),taby(400),xmin,xmax,ymin,ymax
c
c 	Estime le flux solaire pour le jour requis.
 	numday = 1
 	call euv91(iout,yyddd,numday,w1nm,eVinf,w2nm,eVsup,eVmid,wmid,
     .		Pflux,Eflux,impress)
c 	Reordonne les tableaux
 	do ien=1,39
 	  eVmiddp(ien) = dble(eVmid(ien))
 	enddo
 	do ien = 1,39
 	  w1nmor (ien) = w1nm (index(ien))
 	  w2nmor (ien) = w2nm (index(ien))
 	  wmidor (ien) = wmid (index(ien))
 	  eVinfor(ien) = eVinf(index(ien))
 	  eVsupor(ien) = eVsup(index(ien))
 	  eVmidor(ien) = eVmid(index(ien))
 	  Efluxor(ien) = Eflux(index(ien))
 	  Pfluxor(ien) = Pflux(index(ien))
 	enddo
c
c 	Interpolation en lin/log du flux solaire sur les longueurs d'
c 	onde.
 	do ien=1,39
 	  tabx(ien) = wmidor(ien)
 	  taby(ien) = log10(Pfluxor(ien))
 	enddo
 	do ien = 1,nen
 	  ondenm(ien) = 1240./e(ien)
 	enddo
 	call intlin (39,tabx,taby,nen,ondenm,PflxW)
 	do ien = 1,nen
 	  PflxW(ien) = 10.**PflxW(ien)
 	enddo
c
c 	Interpolation en log/log du flux solaire sur les energies
 	do ien=1,39
 	  tabx(ien) = log10(eVmidor(ien))
 	  taby(ien) = log10(Pfluxor(ien))
 	enddo
 	do ien = 1,nen
 	  trav(ien) = log10(e(ien))
 	enddo
 	call intlin (39,tabx,taby,nen,trav,PflxE)
 	do ien = 1,nen
 	  PflxE(ien) = 10.**PflxE(ien)
 	enddo
c
c 	Quelle est la meilleure interpolation? En lin/log sur les 
c 	longueurs d'onde ou en log/log sur les energies?
c 	Mystere... 
 	do ien = 1,nen
 	  Pflx(ien) = PflxW(ien)
 	enddo
c
c 	Que faut-il conserver? L'energie dans chaque boite, l'energie
c 	totale, le nombre de photons dans chaque boite, le nombre 
c	de photons total?
c 	Conservation du nombre de photons total:
 	totin = 0.
 	do ien = 1,39
 	  totin = totin + Pfluxor(ien)
 	enddo
 	totout = 0.
 	do ien = 1,nen
 	  totout = totout + Pflx(ien)
 	enddo
 	if (totout.ne.0.)then
 	  rapp = totin/totout
 	  write(6,*)'rapport de normalisation : ',rapp
 	  write(iout,*)'rapport de normalisation : ',rapp
 	  do ien = 1,nen
 	    Pflx(ien) = Pflx(ien)*rapp
 	  enddo
 	endif
c 	
 	write(iout,*)
 	write(iout,*)'Apres interpolation et normalisation du nombre',
     . 	 	     ' total de particules'
 	write(iout,*)'Energies [eV]'
 	write(iout,1040)(e(ien),ien=1,nen)
 	write(iout,*)'Longueurs d''ondes [nm]'
 	write(iout,1040)(ondenm(ien),ien=1,nen)
 	write(iout,*)'Flux [Photons.cm-2.s-1.nm-1]'
 	write(iout,1040)(Pflx(ien),ien=1,nen)
1040 	format(7(1pe10.2))
c
 	end
c
c-----------------------------------------------------------------------
c
 	subroutine euv91(iout,yyddd,numday,w1nm,eVinf,w2nm,eVsup,
     .			eVmid,wmid,Pflx,Eflx,impress)
c
c 	Re-ecrit en subroutine 07/93 (jl)
c 	PARAMETRES D'ENTREE : 
c	yyddd :  Enter YYDDD date between 68172 - 88366: '
c	numday
c 	PARAMETRES DE SORTIE
c 		TABLEAUX REELS DE 39 ELEMENTS
c 	w1nm : borne inf de la boite de longueur d'onde [nm] pour le 
c 		flux
c 	eVinf : borne inf de la boite d'energie [eV] pour le flux
c 	w2nm : borne sup de la boite de longueur d'onde [nm] pour le 
c 		flux
c 	eVsup : borne inf de la boite d'energie [eV] pour le flux
c 	eVmid : milieu de la boite d'energie [eV].
c 	wmid  : milieu de la boite de longueur d'ondes [nm].
C	Eflx  :	energy flux vector (39) = Eflux pour eviter pb. common
C	Pflx  :	photon flux vector (39) = Pflux pour eviter pb. common
c
c	program euv91
C********1*********2*********3*********4*********5*********6*********7**
C     This program creates the solar EUV irradiance for 18-1050 A based
C     on a multiple linear regression model for solar EUV flux. 
C     Reference: Tobiska, W.K., Revised solar extreme ultraviolet flux 
C     model, JATP, (accepted) 1991. The preprint was XXVIII 
C     COSPAR, The Hague, paper STP I.2-P1, 1990. (See accompanying file 
C     EUV91.TXT)
C     Files required for model to run:
C	     euv91coef.txt	= ASCII file of model coefficients
C	     euv91index2.dat	= ASCII file of proxy values
C     The model values produced in this program include: energy and 
C	photon flux for 39 EUV wavelength intervals or lines
C     Files created by running this model: fyyddd.dat (ASCII)
C      Note: DEC equipment requires IDINT function to be JINT
C
C     W. Kent Tobiska       Baseline version: 10-16-90    Rev: 4-30-91
C********1*********2*********3*********4*********5*********6*********7**
C			VARIABLES AND DEFINITIONS
C********1*********2*********3*********4*********5*********6*********7**
C	i,j,num,sign	loop and index control
C	frcnum,frcsin   decade counter
C	expont		flag for 2 digit ASCII exponent characters
C	expsin,expnum	sign and number of exponent (ASCII chars)
C	yyddd		date in YYDDD format
C	numday		number of days for calculating flux
C	fname		filename for flux data file
C	c		floating point array (12,52) of model coefficients
C	number		temporary number constructed from ASCII coef file
C	indice		floating point array (10,3750) of proxy values
C	sumflx		temporary variable for summing photon flux
C	Eflux		energy flux vector (39)
C	Pflux		photon flux vector (39)
C	date		date
C	model		model photon flux vector (52)
C	W		missing data scaling function vector (52)
C	S		missing data step function vector (4)
C	line1		character string for ASCII line
C	line2		character string for ASCII line
C	w1		beginning wavelengths vector (39) in A
C	w2		ending wavelengths vector (39) in A
C	w1nm		beginning wavelengths vector (39) in nm (JL, 93)
C	w2nm		ending wavelengths vector (39) in nm (JL, 93)
C	wave		average wavelengths vector (39)
C	itemp		temporary vector (6)
C	temp		temporary vector (6)
C	row		loop control for indices row number on a date
C	colum		loop control for indices column number on a date
C********1*********2*********3*********4*********5*********6*********7**
C	implicit	none	
	include 'TRANSPORT.INC'
	integer		yyddd,i,j,numday,num
	real		c(12,52),indice(10,3750),S(4),w1(39),w2(39)
	real		date,Eflux(39),Pflux(39),wave(39)
	real		Eflx(39),Pflx(39)
	real		w1nm(39),w2nm(39),wmid(39),eVinf(39),
     .			eVsup(39),eVmid(39)
	character*10	fname
c	data		numday /1/
c	data		fname(1:1),fname(7:10) /'f','.dat'/
 	integer impress
	common	/coeff/ c,indice,date,wave,w1,w2,S,Eflux,Pflux

C********1*********2*********3*********4*********5*********6*********7**
C	Initialize variables
C********1*********2*********3*********4*********5*********6*********7**
	write (*,*) ' Reading the coefficient and indices tables'
		call rddata
c	write (*,'(A,$)') ' Enter YYDDD date between 68172 - 88366: '
c		read  (*,'(I5)') yyddd
c	write (*,'(A,I5,A,$)') ' Enter # of days including ',yyddd,': '
c		read  (*,'(I10)') numday
 	date = real(yyddd)

C********1*********2*********3*********4*********5*********6*********7**
C	Calculate the flux for each day by first determining the date 
C	and proxy values for that date then calculating the flux.
C********1*********2*********3*********4*********5*********6*********7**
	do 10 i=0,numday-1
	  date = real(yyddd+i)
	  call getind
	  call flux
	  do 20 j=1,5
c	    num = jint(date*(10.**(j-5)))-jint(date*(10.**(j-6)))*10
c 	    num = idint(date*(10.**(j-5)))-idint(date*(10.**(j-6)))*10
  	    num = int(date*(10.**(j-5)))-int(date*(10.**(j-6)))*10
c	    fname(j+1:j+1) = char(num+ichar('0'))
20	  continue

C********1*********2*********3*********4*********5*********6*********7**
C	Write the results to an output file
C********1*********2*********3*********4*********5*********6*********7**

1000 	  format('  EUV 91 Energy and Photon Flux for ',i10,/)
1010 	format (5x,'Wavelength',5x,'Middle',6x,'Energy',6x,'Middle',8x,
     .		'Photon',5x,'Energy',/,58x,'flux',8x,'flux',/,
     .		7x,'[nm]',10x,'[nm]',8x,'[eV]',8x,'[eV]',9x,'[cm-2.',
     . 		5x,'[eV.cm-2.',/,3x,13('-'),12x,14('-'),15x,'s-1.nm-1]',
     .		3x,'s-1.nm-1]',/)
1020 	format
     .	    (2(1F8.1,1x),1f8.1,1x,2(1F7.2,1x),1f7.2,6x,2(1p1E9.3,3x))
c
 	  if (impress.eq.1)then
	    write (iout,1000)yyddd+i 
	    write (iout,1010) 
 	  endif
c 	  w1 et w2 sont en angstrom. On multiplie par 10 pour avoir des
c 	  nm.
	  do 30 j=1,39
	    w1nm(j) = w1(j)/10.
	    w2nm(j) = w2(j)/10.
	    wmid(j) = wave(j)/10.
 	    eVsup(j)= 12400./w1(j)
 	    eVinf(j)= 12400./w2(j)
 	    Pflx(j) = Pflux(j)
c
c 	    wave(j) est defini par Tobiska comme le milieu de la 
c 	    boite en longueur d'onde. Le flux d'energie est alors
c 	    calcule selon pflux*12400/wave. Dans ce cas, on a :
c 	    eVmid(j)= 12400./wave(j)
c 	    On peut prendre aussi comme energie pour le calcul du
c 	    flux d'energie le milieu de la boite en energie (ou un
c 	    point quelconque dans l'intervalle) :
  	    eVmid(j)=eVinf(j)+(eVsup(j)-eVinf(j))*0.5
c 	    Dans ce cas, le flux d'energie est :
 	    Eflux(j) = eVmid(j) * Pflux(j)
 	    Eflx(j) = Eflux(j)
 	    if (impress.eq.1) write(iout,1020) w1nm(j),w2nm(j),wmid(j),
     .		     eVinf(j),eVsup(j),eVmid(j),Pflux(j),Eflux(j)
30	  continue
10	continue
	end
c
c-----------------------------------------------------------------------
c
	subroutine rddata
C********1*********2*********3*********4*********5*********6*********7**
C	This subroutine reads in the coefficients and indices into the 
C	common block COEFF. It may seem awkward; however, it is useful 
C	to keep the coefficient and indices ASCII files in their 
C	present format for easy transfer to other users.
C********1*********2*********3*********4*********5*********6*********7**
C	implicit	none	
	include 'TRANSPORT.INC'
	integer		i,j,num,sign,itemp(6)
	integer		expsin,expnum,frcsin,frcnum
	real		c(12,52),indice(10,3750),S(4),temp(6),number
	real		date,Eflux(39),Pflux(39)
	real		w1(39),w2(39),wave(39)
	character*1	expont
	character*80	line1,line2
	common	/coeff/ c,indice,date,wave,w1,w2,S,Eflux,Pflux
 	real w1dat(39),w2dat(39)
C********1*********2*********3*********4*********5*********6*********7**
C	Default wavelength lines and intervals
C********1*********2*********3*********4*********5*********6*********7**
	data	    w1dat /18.62,30.02,50.52,100.54,150.10,200.02,
     1			256.32,284.15,251.10,303.31,303.78,303.31,
     2			368.07,356.01,401.14,465.22,453.00,500.00,
     3			554.37,584.33,554.37,609.76,629.73,609.76,
     4			650.30,703.36,701.00,765.15,770.41,787.71,
     5			750.01,801.00,851.00,901.00,977.02,951.00,
     6			1025.72,1031.91,1001.00/
	data	    w2dat /29.52,49.22,99.99,148.40,198.58,249.18,
     1			256.32,284.15,299.50,303.31,303.78,349.85,
     2			368.07,399.82,436.70,465.22,499.37,550.00,
     3			554.37,584.33,599.60,609.76,629.73,644.10,
     4			700.00,703.36,750.00,765.15,770.41,790.15,
     5			800.00,850.00,900.00,950.00,977.02,1000.00,
     6			1025.72,1031.91,1050.00/
c
 	do i=1,39
 	  w1(i) = w1dat(i)
 	  w2(i) = w2dat(i)
 	enddo

C********1*********2*********3*********4*********5*********6*********7**
C	Get average wavelength and read in the coefficients 
C********1*********2*********3*********4*********5*********6*********7**
	do 10 i = 1,39
	  wave(i) = (w1(i)+w2(i))/2.
10	continue

	open (unit=fic_euvcoeff,
     .		file=data_path(1:lpath_data)
     &                     //'euv91coef.txt',
     .		status='old',form='formatted')
	do 20 i = 1,3 
	  read(fic_euvcoeff,'(A)') line1
20	continue
	do 30 i = 1,52
	  frcnum = 0
	  frcsin = 1
	  expont = 'n'
	  expnum = 0
	  expsin = 1
	  num = 0
	  number = 0.
	  sign = 1
  	  read (fic_euvcoeff,'(2(I4,1x),I1,1x,G12.5E4,1x,A51)') 
     1		itemp(1),itemp(2),itemp(3),c(4,i),line2
	  c(1,i) = dble(itemp(1))
	  c(2,i) = dble(itemp(2))
	  c(3,i) = dble(itemp(3))
	  do 40 j = 1,51
	    if (line2(j:j) .ne. ' ') then 
	      if (line2(j:j) .eq. '-') then 
		if (j .eq. 1) then 
		  sign = -1
		else if (line2(j-1:j-1) .eq. ' ') then
		  sign = -1
		else if (expont .eq. 'y') then
		  expsin = -1
		endif
	      else if (line2(j:j) .eq. '.') then 
		frcsin = -1
	        frcnum = 1
	      else if (expont .eq. 'y') then 
		if (line2(j:j) .eq. '+') then
		  expsin = 1
		else if (line2(j:j) .eq. '0') then
		  continue
		else if (line2(j:j) .ne. '-') then
		  expnum = ichar(line2(j:j))-ichar('0')
		endif
	      else if (line2(j:j) .eq. 'e') then
		expont = 'y'
	      else 
		if (frcsin .eq. -1) then 
		  number = (ichar(line2(j:j))-ichar('0')) *
     1			(10.**(frcsin*frcnum)) + number
		  frcnum = frcnum+1
		else if (frcsin .eq. 1) then 
		  number = (ichar(line2(j:j))-ichar('0')) + 
     1			(10.**(frcsin*frcnum)) * number
		  frcnum = 1
	        endif
	      endif
	    else
	      c(5+num,i) = sign*number*(10.**(expsin*expnum))
	      num = num+1
	      frcnum = 0
	      frcsin = 1
	      number = 0.
	      sign = 1
	      expont = 'n'
	      expnum = 0
	      expsin = 1
	      if (num .eq. 8) goto 30
	    endif
40	  continue
30	continue
	close (unit=fic_euvcoeff)
C********1*********2*********3*********4*********5*********6*********7**
C	Read in the indices
C********1*********2*********3*********4*********5*********6*********7**
	open(unit=fic_euvindex,file=
     .	data_path(1:lpath_data)
     &        //'euv91index2.dat'
     .		,status='old',form='formatted')
  	do 50 i = 1,2
	  read(fic_euvindex,'(A)') line1
50	continue
  	do 60 i = 1,3750
	  read (fic_euvindex,
     1 '(1x,F5.0,2(1x,G11.5E4),2(1x,F3.0),3x,F5.0,2(1x,G11.5E4),
     2   2(1x,F3.0))')temp(1),indice(2,i),indice(3,i),temp(2),temp(3),
     3	    temp(4),indice(7,i),indice(8,i),temp(5),temp(6)
	  indice( 1,i) = dble(temp(1))
	  indice( 4,i) = dble(temp(2))
	  indice( 5,i) = dble(temp(3))
	  indice( 6,i) = dble(temp(4))
	  indice( 9,i) = dble(temp(5))
	  indice(10,i) = dble(temp(6))
60	continue
	close (unit=fic_euvindex)

	return
	end
c
c-----------------------------------------------------------------------
c

	subroutine flux
C********1*********2*********3*********4*********5*********6*********7**
C	This subroutine uses the coefficients and indices in the 
C	common block COEFF to calculate the EUV energy flux for a date.
C	For mid 1968 through early 1977 when neither Lyman alpha nor 
C	He I 10830 were present, the model calculates the chromospheric
C	emission using F10.7 linear relationship with Lyman alpha based
C	upon Barth et al., GRL, 17, 571-574, 1990.
C********1*********2*********3*********4*********5*********6*********7**
C	implicit	none	
	integer		i,j
	real		c(12,52),indice(10,3750),sumflx,Eflux(39)
	real		Pflux(39),wave(39),w1(39),w2(39)
	real		date,model(52),W(52),S(4),near0
	data		W,model,sumflx,near0 /52*1.,52*0.,0.,1.e-5/
	common	/coeff/ c,indice,date,wave,w1,w2,S,Eflux,Pflux

C********1*********2*********3*********4*********5*********6*********7**
C	Replace the missing proxies on a date with surrogate value
C********1*********2*********3*********4*********5*********6*********7**
	if ((S(1) .lt. near0).and.(S(2) .gt. near0)) S(1) = S(2)
	if  (S(1) .lt. near0) S(1) = (8.7e8)*S(3) + 1.9e11
	if  (S(2) .lt. near0) S(2) = S(1)

C********1*********2*********3*********4*********5*********6*********7**
C	Create the empirically modeled flux from Equations 3, 4, 5 in 
C	JATP paper. Note that the "W" scaling function is now unused 
C	and all elements are set to 1. in this version.
C********1*********2*********3*********4*********5*********6*********7**
	do 10 j=1,52
	  sumflx = 0.
	  do 20 i=1,4
	    if (i .le. 2) W(j) = (1.+c(8+i,j)*exp(-S(i)*1.e-10))*W(j)
	    if (i .gt. 2) W(j) = (1.+c(8+i,j)*exp(-S(i)*.7))*W(j)
	    sumflx = c(i+4,j)*S(i) + sumflx
20	  continue
	  model(j) = (c(4,j) + sumflx) * W(j)
10	continue
C********1*********2*********3*********4*********5*********6*********7**
C	Combine the chromospheric and coronal flux into intervals
C********1*********2*********3*********4*********5*********6*********7**
	Pflux( 1)=model( 1)
	Pflux( 2)=model( 2)+model( 3)
	Pflux( 3)=model( 4)+model( 5)
	Pflux( 4)=model( 6)+model( 7)
	Pflux( 5)=model( 8)+model( 9)
	Pflux( 6)=model(10)+model(11)
	Pflux( 7)=model(12)
	Pflux( 8)=model(13)
	Pflux( 9)=model(14)+model(15)
	Pflux(10)=model(16)
	Pflux(11)=model(17)
	Pflux(12)=model(18)
	Pflux(13)=model(19)
	Pflux(14)=model(20)+model(21)
	Pflux(15)=model(22)+model(23)
	Pflux(16)=model(24)
	Pflux(17)=model(25)+model(26)
	Pflux(18)=model(27)+model(28)
	Pflux(19)=model(29)
	Pflux(20)=model(30)
	Pflux(21)=model(31)
	Pflux(22)=model(32)
	Pflux(23)=model(33)
	Pflux(24)=model(34)+model(35)
	Pflux(25)=model(36)+model(37)
	Pflux(26)=model(38)
	Pflux(27)=model(39)
	Pflux(28)=model(40)
	Pflux(29)=model(41)
	Pflux(30)=model(42)
	Pflux(31)=model(43)+model(44)
	Pflux(32)=model(45)
	Pflux(33)=model(46)
	Pflux(34)=model(47)
	Pflux(35)=model(48)
	Pflux(36)=model(49)
	Pflux(37)=model(50)
	Pflux(38)=model(51)
	Pflux(39)=model(52)

C********1*********2*********3*********4*********5*********6*********7**
C	Make an energy flux from the photon flux
C********1*********2*********3*********4*********5*********6*********7**
	do 30 i=1,39
c 	  L'expression initiale est exprimee en unites de 10-7 Joules
c	  Eflux(i) = Pflux(i)*(12400.*1.6022e-12)/wave(i)
c 	  J'exprime un flux d'energie en eV
 	  Eflux(i) = Pflux(i)*12400./wave(i)
30	continue

	return
	end
c
c-----------------------------------------------------------------------
c

	subroutine getind
C********1*********2*********3*********4*********5*********6*********7**
C	This subroutine gets the row and column indicators for a given
C	date.
C********1*********2*********3*********4*********5*********6*********7**
C	implicit	none	
	integer		i,row,colum
	real		c(12,52),indice(10,3750),Eflux(39)
	real		Pflux(39),wave(39),w1(39),w2(39),date,S(4)
	common	/coeff/ c,indice,date,wave,w1,w2,S,Eflux,Pflux

C********1*********2*********3*********4*********5*********6*********7**
C	Find the row and column numbers for indices on a given date then
C	place the proxy values for that date into "S" variable. Stop on 
C	invalid date.
C********1*********2*********3*********4*********5*********6*********7**
	do 10 i=1,3750
	  if (ifix(date) .eq. ifix(indice(1,i))) then
	    row = i
	    colum = 1
	    goto 20
	  else if (ifix(date) .eq. ifix(indice(6,i))) then
	    row = i
	    colum = 6
	    goto 20
	  endif
10	continue

20	continue
	S(1)  = indice(colum+1,row)
	S(2)  = indice(colum+2,row)
	S(3)  = indice(colum+3,row)
	S(4)  = indice(colum+4,row)

	if (ifix(date) .ne. ifix(indice(colum,row))) then
	  write(*,*) ' Invalid date - program terminated'
	  stop
	endif

	return
	end
