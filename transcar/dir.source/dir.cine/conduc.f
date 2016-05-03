c
c-----------------------------------------------------------------
c
        subroutine conductivite(nalt,altkm,denelc,densneut,
     &		tneutre,temelc,temion,chideg,glat,glong,year,day,
     &		hrloc,cped,chal,ratHoP,cpedsum,chalsum,cpedCS,chalCS,
     &		ratHoPsum,gyreave,gyriave,collOp,collNOp,collO2p,
     &		collionSN,collionRG,collen2,colleo2,colleo1,colle,iprt,
     &		f107,icolin,fic_transout)
c
        implicit none
        include 'TRANSPORT.INC'

        integer, intent(in) :: fic_transout

        integer nalt,ialt,iprt(40),icolin
        integer file_cond
        real altkm(nbralt),denelc(nbralt),O1prate(nbralt),
     &	     densneut(8,nbralt)
        real tneutre(nbralt),temelc(nbralt),temion(nbralt)
        real cped(nbralt),chal(nbralt),ratHoP(nbralt)
 	real omegae(nbralt), omegai(nbralt),gyreave,gyriave
 	real collionSN(nbralt),colle(nbralt),collionRG(nbralt)
 	real collOp(nbralt),collNOp(nbralt),collO2p(nbralt)
 	real collen2(nbralt),colleo2(nbralt),colleo1(nbralt)
 	real ratHoPsum,cpedsum,chalsum,cpedCS,chalCS
 	real glat,glong,day,date,chideg,f107
 	real denne,denn2,deno2,deno1
      	real aionmas,gyre,gyri
 	real fcOp,fcO2p,fcNOp,fci,fciRG,fciSN,fcen2,fceo2,fceo1,fce
 	real ti,te,tn,hrloc,year
 	real pi,cx,conp,conh,field,field120,uma,z50
        real xbid,ybid,zbid,xmbid,decbid
 	real dOp,dO2p,dNOp
 	real xmO1,xmO2,xmNO
c

c 	ENTREES :
c 	---------
c 	nalt,altkm = nombre d'altitudes et altitudes en km
c 	denelc(altkm) = densite electronique cm-3
c 	densneut(8,altkm) = densites neutres. Seules les 3 premieres
c 	 		(n2,o2,o) sont utilisees.
c	tneutre(altkm) = temperature neutre
c 	temelc(altkm) = temperature electronique (pour le calcul des 
c 	 		frequences de collision)
c 	temion(altkm) = temperature ionique (pour le calcul des 
c 	 		frequences de collision)
c 	chideg = angle zenithal solaire (pour le calcul des valeurs
c 	 		theoriques de catherine)
c	glat,glong,year,day,hrloc = pour le calcul de B et de la 
c 			composition ionique.
c 	iprt = 1 si on imprime les resultats dans le fichier file_cond
c 	file_cond = numero du fichier de sortie
c 	f107 = n'est la que pour info, pour ne pas melanger les fichiers
c 	icolin = 0 si freq. de coll. ions-neutres de Risbeth and Garriot
c 	       = 1 si freq. de coll. ions-neutres de Schunk et Nagy
c
c 	SORTIES :
c 	---------
c       cped(altkm),chal(altkm),ratHoP(altkm) = conductivites Pedersen,
c 	  	Hall et leur rapport a chaque altitude
c 	cpedsum,chalsum,ratHoPsum = leurs valeurs integrees.
c
c 	cpedCS,chalCS = valeurs theoriques du modele de Catherine
c
c 	On sort aussi les parametres qui servent au calcul des cond :
c
c     	gyreave,gyriave = gyrofrequences moyennes  eB/m
c 	collOp(altkm),collNOp(altkm),collO2p(altkm) = frequences de 
c 	 	collisions ion/neutre calculees par le pgmme de Vincent,
c 		selon les formules de Schunk et Nagy.
c	collionSN(altkm),collionRG(altkm) = frequence de collision 
c 		moyenne ion/neutre selon Schunk et Nagy ou selon 
c 		Rishbeth et Garriot (qui est ce que Catherine avait 
c 		utilise).
c 	collen2(altkm),colleo2(altkm),colleo1(altkm),colle(altkm) = 
c 		frequences de collisions electrons/neutres et 
c 		totale electron/neutre. Ouf, tout le monde utilise la
c 		meme formule...
c
c 	Pour comparer avec l'article de Catherine, on se donne un
c 	parametre "senior" qui vaut 1 si on se place dans les memes
c 	conditions qu'elle, et 0 sinon
c
 	integer senior
 	data senior/1/
 	file_cond = 48

        write(6,*)'Calcul de conductivite                     ' ,'[A'
c
        date= year+day/100.
        if(year.lt.100.)date=1900.+ year+day/100.
c
c 	On cherche la composition ionique a partir du modele de Chantal
 	call compos(hrloc,nalt,altkm,O1prate,z50)
c
c 	On cherche la valeur du champ a 120 km pour se placer dans les
c 	memes conditions que Catherine Senior.
        call bfield(date,120.,glat,glong,xbid,ybid,zbid,
     .               field120,xmbid,decbid)
 	field120=field120*1.0e-09
c
 	pi = 4.*atan(1.)
 	cx = cos(chideg*pi/180.)
c
	xmO1 = 15.9994*uma
 	xmO2 = 2.* xmO1
 	xmNO = (14.0067+15.9994)*uma
c
        write(fic_transout,990)chideg,cx
 	if (iprt(24).eq.1)then
 	  write(fic_transout,*)' Gyrofrequences :'
 	  write(fic_transout,*)'   altkm     %O+   %O2+   %NO+ ',
     .		       '   aionmas  field      gyre       gyri'
 	endif
c
 	gyreave = 0.
 	gyriave = 0.
        uma = 1.66056e-27
        do ialt=1,nalt
c 	  Densites electronique
          denne=denelc(ialt)
 	  dOp  = O1prate(ialt)
 	  dO2p = (1.-O1prate(ialt))/2.
 	  dNOp = (1.-O1prate(ialt))/2.
c
 	  if(senior.eq. 0)then
c 	    On utilise des parametres modernes.
c 	    Masse ionique : On considere que ce qui n'est pas O+ est un
c 	    melange a 50% O2+, % NO+
  	    aionmas = xmO1*dOp+xmO2*dO2p+xmNO*dNOp
c
c	    Neutral densities en cm-3
            denn2=densneut(1,ialt)
            deno2=densneut(2,ialt)
            deno1=densneut(3,ialt)
c
c 	    Temperatures
 	    tn = tneutre(ialt)
c
c 	    Intensite du champ magnetique.
            call bfield(date,altkm(ialt),glat,glong,xbid,ybid,zbid,
     .               field,xmbid,decbid)
c 	    Conversion de B, des gammas en Tesla
 	    field=field*1.0e-09
c
 	  else if(senior.eq.1) then
c 	    On se place dans les conditions de l'article de Catherine
 	    aionmas = 30*uma
 	    call vin2mil(altkm(ialt),tn,denn2,deno2,deno1)
 	    field = field120
c 	    On utilise les freq. de collision Rishbeth et Garriot
 	    icolin = 0
 	  endif
c
 	  ti = temion(ialt)
 	  te = temelc(ialt)
c
 	  if (denne.gt.0.)then
	    call conduct(denn2,deno2,deno1,denne,dOp,dO2p,
     .	        dNOp,te,ti,tn,conP,conH,field,aionmas,gyre,gyri,
     .	        fcOp,fcNOp,fcO2p,fci,fciRG,fciSN,
     .		fcen2,fceo2,fceo1,fce,icolin)
c
	    cped(ialt)    = conP
c 	    A haute altitude, les cond. Hall peuvent prendre des valeurs
c 	    negatives tres petites
 	    chal(ialt)    = max(conH,0.)
            ratHoP(ialt)  = chal(ialt)/cped(ialt)
c
 	    omegae(ialt)  = gyre
 	    omegai(ialt)  = gyri
c
 	    collOp(ialt)  = fcOp
 	    collNOp(ialt) = fcNOp
 	    collO2p(ialt) = fcO2p
 	    collionSN(ialt) = fci
 	    collionRG(ialt) = fciRG
c
 	    collen2(ialt) = fcen2
 	    colleo2(ialt) = fceo2
 	    colleo1(ialt) = fceo1
 	    colle  (ialt) = fce
 	  endif
 	  if (iprt(24).eq.1)write(fic_transout,2000)altkm(ialt),
     .			    dOp,dO2p,dNOp,aionmas/uma,field,gyre,gyri
 	  gyreave = gyreave + gyre
 	  gyriave = gyriave + gyri
 	enddo
 	gyreave = gyreave/nalt
 	gyriave = gyriave/nalt
 	if (iprt(24).eq.1)write(fic_transout,1060)gyreave,gyriave
c
c 	conductivites integrees:
        call hint (nalt,altkm,chal,chalsum)
c 	conversion des km vers les metres
        chalsum = chalsum*1.e+03
        call hint (nalt,altkm,cped,cpedsum)
c 	conversion des km vers les metres
        cpedsum = cpedsum*1.e+03
c 	
c 	Comparaison avec les valeurs theoriques de Catherine Senior.
 	cpedCS =  1.81 + 8.88 * cx
 	chalCS = 21.58 - 0.21 * chideg
c
        ratHoPsum = chalsum/cpedsum
c
 	if (iprt(24).eq.1) then
          write(fic_transout,1040)
          write(fic_transout,1045)
          do ialt=1,nalt
            write(fic_transout,1050)altkm(ialt),collen2(ialt),
     . 		colleo2(ialt),colleo1(ialt),colle(ialt)
          enddo
          write(fic_transout,1055)
          do ialt=1,nalt
            write(fic_transout,1050)altkm(ialt),collOp(ialt),
     . 	       collNOp(ialt),collO2p(ialt),collionSN(ialt),
     .		collionRG(ialt)
          enddo
c
          write(fic_transout,1000)
          do ialt=1,nalt
            write(fic_transout,1010)altkm(ialt),chal(ialt),cped(ialt),
     . 		ratHoP(ialt),denelc(ialt),temelc(ialt),temion(ialt)
          enddo
 	endif
c       write(6,*)'                                           ' ,'[A'
c
 	if (iprt(23).eq.1)then
          write(fic_transout,1020)
          write(fic_transout,1030) glat,f107,chideg,chalsum,cpedsum,
     .			 ratHoPsum,chalCS,cpedCS
 	endif
c
c 	L'ouverture en mode 'old' d'un fichier suivit d'une ecriture
c 	positionne l'ecriture en fin de fichier (ajout des nouvelles
c 	donnees). En lecture, la ligne lue est la premiere.
c	Ceci est une specificite aix, car la norme ANSII ne specifie
c 	rien.
 	open(file_cond,file = chemin(1:lpath)//'dir.output/conduc.res',
     .		status = 'old')
        write(file_cond,1030)glat,hrloc,chideg,chalsum,cpedsum,
     .		ratHoPsum,chalCS,cpedCS
c       do ialt=1,nalt
c         write(file_cond,1010)altkm(ialt),chal(ialt),cped(ialt),
c    . 		ratHoP(ialt),denelc(ialt),temelc(ialt),temion(ialt)
c       enddo
 	close(file_cond)
c
c       write(6,1030)glat,hrloc,chideg,chalsum,cpedsum,ratHoPsum,
c    .		chalCS,cpedCS
c
990   	format (60('-'),/,'chideg,coschi=',2(1p1e10.2))
1000  	format (/,'conductivites :',/,
     &    '   Altitude       Hall     Pedersen      H/P     Ne[cm-3]',
     &    '   Te [K]         Ti[K]')
1010  	format (f10.2,3x,2(1p1e11.2),1x,1p1e10.2,1x,1p1e10.2,
     &    2(0p1f10.2))
1020    format(/,'COND. INTEGREES [MHO]:',/,
     &    'Latitude  f107    Chideg',
     &    '     Hall   Pedersen    H/P   Hall.C.S. Ped.C.S.')
1030    format (1f6.2,7f9.3)
1040 	format(/,'Frequences de collision')
1045    format(/,
     .   '   Altitude   e/N2      e/O2       e/O     Total e')
1050 	format(1f10.2,5(1pe10.2))
1055    format(/,
     .   '   Altitude  O+/Ntres NO+/Ntres O2+/Ntres Total ions NuionCS')
1060  	format('Gyrofrequences moyennes      :',2(1pe10.2))
2000 	format(1f10.2,3f7.4,1f10.2,3(1p1e10.2))

      return
      end       
c
c---------------------------------------------------------------
c	
 	subroutine conduct(denn2,deno2,deno1,denne,dOp,dO2p,
     .		dNOp,te,ti,tn,conP,conH,field,aionmas,gyre,gyri,
     .		fcOp,fcNOp,fcO2p,fci,fciRG,fciSN,
     .		fcen2,fceo2,fceo1,fce,icolin)
c
 	implicit none
c
 	real e,elmas,aionmas,pi,esb,srt,gyre,gyri,field
 	real denn2,deno2,deno1,denne
 	real dOp,dO2p,dNOp
 	real fcOp,fcOpO,fcOpN2,fcOpO2,burnside
        real fcNOp,fcNOpO,fcNOpN2,fcNOpO2
 	real fcO2p,fcO2pO,fcO2pN2,fcO2pO2
 	real fci,fciRG,fciSN
 	real fcen2,fceo2,fceo1,fce
 	real te,ti,tn,ae,ai,conP,conH
 	integer icolin
c
	e=1.6021e-19
	elmas=9.1091e-31
c
 	pi = 4.*atan(1.)
c
	esb=e/field
	srt=sqrt(tn)
c
 	gyre=e*field/elmas
 	gyri=e*field/aionmas
c
 	burnside = 1.7
 	call nuoion(fcOp,fcOpO,fcOpN2,fcOpO2,ti,tn,denO1,denN2,denO2,
     .			burnside)
 	call noion(fcNOp,fcNOpO,fcNOpN2,fcNOpO2,denO1,denN2,denO2)
 	call o2ion(fcO2p,fcO2pO,fcO2pN2,fcO2pO2,ti,tn,denO1,denN2,
     .		 	denO2)
c 	On passe de frequences reduites aux frequences reelles.
 	fcOp  = dOp*fcOp
 	fcO2p = dO2p*fcO2p
 	fcNOp = dNOp*fcNOp
 	fciSN   = fcOp+fcO2p+fcNOp
c
 	fciRG=3.75e-10*(denn2+deno2+deno1)
c
 	if(icolin.eq.0)then
c 	  frequence de collision ions-neutres:
c 	  Formule Risbeth and Garriot, Introduction to ionospheric 
c 	  physics, intern geophys series, ed J van Mieghem, Academic 
c 	  Pess, NY, 1969, p 130
 	  fci=fciRG
 	else
c 	  On appelle les pgmmes de vincent. Les freq. de coll sont alors
c 	  calculees selon Schunk et Nagy, Rev Geophys, 16,355-399, 1978.
 	  fci   = fciSN
 	endif
c
c 	frequence de collision electrons-neutres from Banks et Kockarts,
c	1973, et Schunk et Walker, 1973 :
c 	Ce sont les memes que dans l'article de Catherine.
	fcen2=2.33e-11*denn2*(1.-1.21e-4*tn)*tn
	fceo2=1.22e-10*deno2*(1.+3.6e-2*srt)*srt
	fceo1=2.8e-10*deno1*srt
c 	On passe de frequences reduites aux frequences reelles. Les
c 	frequences e/neutre tiennent deja compte des densite neutres.
c 	(These Catherine, p23)
	fce=fcen2+fceo2+fceo1
c
 	ae=gyre**2+fce**2
 	ai=gyri**2+fci**2
 	conP=denne*1.e+06*esb*(fci*gyri/ai+fce*gyre/ae)
	conH=denne*1.e+06*esb*(gyre*gyre/ae-gyri*gyri/ai)
c
	return
	end
c
c----------------------------------------------------------------------
c
 	subroutine vin2mil(zkm,tn,denn2,deno2,deno1)
c
c 	Sous programme d'atmosphere neutre de Denis Alcayde. Pour le
c 	fonctionnement et tous les commentaires, voir
c 	An analytical static model of temperature and composition from
c 	20 to 2000 km altitude, annales geophys., vol 37, 515-528, 1981.
c
c
 	implicit none

 	real zkm,tn,denn2,deno2,deno1
c
 	real tdz,xn(5),xntot,ro,press,xm,h
 	common/outpt1/tdz,xn,xntot,ro,press,xm,h
c
 	real paresse,gtzbaz
 	common/outpt2/ paresse,gtzbaz
c 
 	real texo,tzbaz,tm,ts
 	common /modtemp/texo,tzbaz,tm,ts
c
 	real zbaz,zref,zm,zs
 	common /modalts/ zbaz,zref,zm,zs
c
 	real xno(5)
 	common /modconc/xno
c
 	real alfa(5),amas(5)
 	common /alfamas/alfa,amas
c
 	real zp,hp
 	common /recom/zp,hp
c
 	real amu,r,gzer,re,amsol
 	common /cstes/ amu,r,gzer,re,amsol
c
 	real fact,recomb,g,alt
 	real x,y,z
 	integer i
c
c
 	recomb(x,y) = (1.+exp(-(x-zp)/hp))/(1.+exp(-(y-zp)/hp))
c
 	alt = zkm
c
 	gtzbaz = 3.*(tzbaz-tm)/(zbaz-zm)
 	paresse = gtzbaz/(texo-tzbaz)
c
 	xntot = 0.
 	ro = 0.
 	xm = amsol
 	call fzrftoz(zref,xm,tdz,fact,g,1.)
 	do i = 1,5
 	  xn(i) = xno(i)
 	  xntot = xntot +xn(i)
 	  ro = ro +amas(i)*xn(i)
 	enddo
c
 	if(alt.eq.zref) then
 	  ro = ro*amu
 	  press = r*tdz/xm
 	  h = press/g
 	  press = press*ro
 	  denn2 = xn(1)/1.e+06		!conversion en cm-3
 	  deno2 = xn(2)/1.e+06		!conversion en cm-3
 	  deno1 = xn(3)/1.e+06		!conversion en cm-3
 	  tn = tdz
 	else if (alt.gt.zref) then
 	  z = amin1(alt,zbaz)
 	  xntot = 0.
 	  ro = 0.
  	  do i = 1,5
 	    call fzrftoz(z,amas(i),tdz,fact,g,alfa(i))
 	    xn(i) = xn(i)*fact
 	    if(i.eq.3)xn(i) = xn(i)/recomb(z,zref)
 	    xntot = xntot + xn(i)
 	    ro = ro+xn(i)*amas(i)
 	  enddo
 	  xm = ro/xntot
 	  ro = ro*amu
 	  press = r*tdz/xm
 	  h = press/g
 	  press = press*ro
 	  if(alt.le.zbaz) then
 	    denn2 = xn(1)/1.e+06		!conversion en cm-3
 	    deno2 = xn(2)/1.e+06		!conversion en cm-3
 	    deno1 = xn(3)/1.e+06		!conversion en cm-3
 	    tn = tdz
 	    return
 	  endif
c
 	  xntot = 0.
 	  ro = 0.
 	  z = alt
 	  do i = 1,5
 	    call fzbztoz(z,amas(i),tdz,fact,g,alfa(i))
 	    xn(i) = xn(i)*fact
 	    if(i.eq.3)xn(i) = xn(i)/recomb(z,zbaz)
 	    xntot = xntot + xn(i)
 	    ro = ro+xn(i)*amas(i)
 	  enddo
 	  xm = ro/xntot
 	  ro = ro*amu
 	  press = r*tdz/xm
 	  h = press/g
 	  press = press*ro
c
 	  denn2 = xn(1)/1.e+06		!conversion en cm-3
 	  deno2 = xn(2)/1.e+06		!conversion en cm-3
 	  deno1 = xn(3)/1.e+06		!conversion en cm-3
 	  tn = tdz
c
 	else if(alt.lt.zref) then
 	  xntot = 0.
 	  ro = 0.
 	  z = amax1(alt,zm)
 	  call fzrftoz(z,xm,tdz,fact,g,1.)
  	  do i = 1,5
 	    xn(i) = xn(i)*fact
 	    if(i.eq.3)xn(i) = xn(i)/recomb(z,zref)
 	    xntot = xntot + xn(i)
 	    ro = ro+xn(i)*amas(i)
 	  enddo
 	  ro = ro*amu
 	  press = r*tdz/xm
 	  h = press/g
 	  press = press*ro
c
  	  if(alt.lt.zm)then
 	    ro = 0.
 	    xntot = 0.
 	    z = alt
 	    xm = amsol
 	    call fzmstoz (z,xm,tdz,fact,g)
 	    do i =1,5
 	      xn(i) = xn(i)*fact
 	      if(i.eq.3)xn(i) = xn(i)/recomb(z,zm)
 	      xntot = xntot +xn(i)
 	      ro = ro+amas(i)*xn(i)
 	    enddo
 	    ro = ro*amu
 	    press = r*tdz/xm
 	    h = press/g
 	    press = press*ro
 	  endif
 	  denn2 = xn(1)/1.e+06		!conversion en cm-3
 	  deno2 = xn(2)/1.e+06		!conversion en cm-3
 	  deno1 = xn(3)/1.e+06		!conversion en cm-3
 	  tn = tdz
 	endif
c
 	return
 	end
c
c ----------------------------------------------------------------------
c
 	subroutine fzbztoz(alt,amas,tdz,fact,gg,alfa)
c
 	implicit none
 	real paresse,gtzbaz
 	common/outpt2/ paresse,gtzbaz
 	real amu,r,gzer,re,amsol
 	common /cstes/ amu,r,gzer,re,amsol
 	real texo,tzbaz,tm,ts
 	common /modtemp/texo,tzbaz,tm,ts
 	real zbaz,zref,zm,zs
 	common /modalts/ zbaz,zref,zm,zs
c
 	real x,y,z,g,gg,t,csi,amas,alt,tdz
 	real e,c,f,d,h,gam,alfa,sig,gzbaz,csy,rezbaz,fact
c
 	g(x) = gzer*(1.+x/re)**(-2)
 	t(x) = texo -(texo-tzbaz)*exp(-sig*x)
 	csi(x,y) = (x-y)*(re+y)/(re+x)
c
 	gg = g(alt)
 	gzbaz = g(zbaz)
 	rezbaz = re + zbaz
 	sig = paresse + 1./rezbaz
 	gam = amas*gzbaz/r/texo/sig*1.e+03
 	csy = csi(alt,zbaz)
 	tdz = t(csy)
 	e = -sig *csy
 	c = tzbaz / tdz
 	f = e*gam
 	d = exp(f)
 	h = gam + alfa
 	fact = d*c**h
c
 	return
 	end
c
c----------------------------------------------------------------------
c
 	subroutine fzrftoz(alt,amas,tdz,fact,gg,alfa)
c
 	implicit none
c
 	real pi,sq
 	common /circ/pi,sq
 	real amu,r,gzer,re,amsol
 	common /cstes/ amu,r,gzer,re,amsol
 	real texo,tzbaz,tm,ts
 	common /modtemp/texo,tzbaz,tm,ts
 	real zbaz,zref,zm,zs
 	common /modalts/ zbaz,zref,zm,zs
c
 	real g,x,cub,cubic,deltaz,a,tborn,binf,bsup
 	real alfa,alt,gg,tdz,amas,fact
c
 	g(x) = gzer*(1.+x/re)**(-2)
	cub(x) = x**3+1.
 	cubic(x) = atan((2.*x-1.)/sq)/sq - alog(1.-3.*x/(1.+x)**2)/6.
c
 	deltaz = zbaz - zm
 	a = (tzbaz/tm -1.)**(1./3.)
 	x = (zref - zm)/deltaz*a
 	tborn = tm*cub(x)
 	tdz = tborn
 	fact = 1.
 	if (alt.eq.zref)then
 	  gg = g(alt)
 	else
 	  binf = cubic(x)
 	  x = (alt-zm)/deltaz*a
 	  bsup = cubic(x) 	
 	  tdz = tm*cub(x)
 	  gg = g((alt+zref)/2.)
 	  fact = amas*gg*deltaz*(binf-bsup)/r/tm/a*1.e+03
 	  fact = exp(fact)*(tborn/tdz)**alfa
 	endif
c
 	return
 	end
c
c---------------------------------------------------------------------
c
 	subroutine  fzmstoz(alt,amas,tdz,fact,gg)
c
 	real alt,amas,tdz,fact,gg
c
 	real pi,sq
 	common /circ/pi,sq
 	real amu,r,gzer,re,amsol
 	common /cstes/ amu,r,gzer,re,amsol
 	real texo,tzbaz,tm,ts
 	common /modtemp/texo,tzbaz,tm,ts
 	real zbaz,zref,zm,zs
 	common /modalts/ zbaz,zref,zm,zs
c
 	real g,x,srta,arctg
c
 	g(x) = gzer*(1.+x/re)**(-2)
 	srta(x) = (tm+ts + (tm-ts)*cos(x))/2.
c
 	gg = g((alt+zm)/2.)
 	deltaz = (zm-zs)/pi
 	x = (zm-alt)/deltaz
 	fact = 2.*amas*gg*deltaz/r/sqrt(tm*ts)*1.e+03
 	if(abs(alt-zs) .ge. 1.e-04)then
 	  arctg = sqrt(ts/tm)*tan(x/2.)
 	  arctg = atan(arctg)
 	  if (alt.lt.zs)arctg = arctg+pi
  	else
 	  arctg = pi/2.
 	endif
 	fact = fact*arctg
 	tdz = srta(x)
 	fact = exp(fact)*tm/tdz
 	gg = g(alt)
c
 	return
 	end
c
c----------------------------------------------------------------------
c
 	block data dt22mil
c
 	implicit none
c
 	real pi,sq
 	common /circ/pi,sq
 	real amu,r,gzer,re,amsol
 	common /cstes/ amu,r,gzer,re,amsol
 	real texo,tzbaz,tm,ts
 	common /modtemp/texo,tzbaz,tm,ts
 	real zbaz,zref,zm,zs
 	common /modalts/ zbaz,zref,zm,zs
 	real alfa(5),amas(5)
 	common /alfamas/alfa,amas
 	real zp,hp
 	common /recom/zp,hp
 	real xno(5)
 	common /modconc/xno
c
 	data pi/3.1415927/,sq/1.7320508/
 	data amu,r,gzer,re /1.660531e-27,8.31434e+03,9.80665,6356.776/
 	data amas /28.,32.,16.,40.,4./,amsol /28.86/
 	data alfa /1.,1.,1.,1.17,0.62/
 	data zbaz,zref,zm,zs /120.,100., 85.,45./
 	data texo,tzbaz,tm,ts /1000.,330.,185.,270./
  	data xno /1.e+19,2.e+18,1.e+18,1.e+17,7.e+13/
 	data zp,hp /95.,2.5/
c
 	end
