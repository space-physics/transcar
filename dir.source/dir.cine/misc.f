c
c 			CONTENT OF MISC.F
c    
c ----  subroutine bfield(date,alt,la,lon,x,y,z,f,dip,dec)
c       subroutine igrf (date,itype,alt,colat,elong,x,y,z,f)
c       igs 1975 magnetic field model
c
c ---- 	subroutine clasdesc(tab,ntab,index)
c 	Classe un tableau en ordre descendant. 
c 	tab : en entree, tableau a classe
c 	      en sortie, tableau classe (et donc original ecrase)
c 	ntab : nbre de points a classer
c 	index : tableau d'entiers contenant en sortie l'ordre de 
c 	        classement.
c
c ----  subroutine intlin (nin,xin,yin,nout,xout,yout)
c 	Linear interpolation.
c 	No restriction on extrema. Extrapolation is performed if 
c       necessary.
c
c ----  subroutine intquad(nin,xin,yin,nout,xout,yout)
c  	Performs quadratic interpolations on array yin(i) vs xin(i)
c    	result is yout(i) at the x points given by xout(i)
c
c ----  subroutine gaussint (ngau,nx,x,fx,sum)
c    	Computes |integrale fx(i) dx(i)| using gaussian integration.
c 	done on ngauss points within each interval. To speed up the 
c 	program, the absissas and weight factors for ngauss=6, 12 or 24 
c 	are stored. Otherwise they are computed.
c
c---- 	subroutine hint(nin,z,fction,fint)
c    	integrate using trapezoid rule
c 	fint=integral(fction(z).dz) on nin points.
c
c ----  subroutine mnmx(tab,ntab,tmin,tmax,linlog)
c 	if linlog=0, finds min and max of tab, defined on ntab points.
c 	if linlog=1, finds the first min non = 0, and the max of tab
c
c ----  subroutine mnmxplt(tmin,tmax,linlog)
c       find best min and max for nice plot.
c 	linlog = 0 if linear axis.
c 	linlog = 1 if log axis.
c
c ----	subroutine mnmxi(tab,ntab,imax,imin)
c 	renvoie les indices du max et min du tableau tab.
c
c ----  subroutine qgauss(m,gmu,gwt)
c 	Computes the absissas gmu and weight factors gwt for a gaussian
c 	integration on m points
c
c ----  subroutine xline(nbline,ijfile)                   
c    	xline skips nbline in file ijfile.
c
c ---- 	subroutine nuoion (xnutot,xnuo,xnun2,xnuo2, ti,tn,deno,denn2,
c    . 			deno2,burnside)
C     	routine by vincent wickwar, sri, may 1975.
c     	routine to find the ion-neutral collision frequencies
c     	for o+, with o, n2, and o2.
c     	these are the *reduced* frequencies ADAPTED FROM schunk
c     	and walker (1973), schunk and nagy (1980)
c
c ---- 	subroutine noion (xnutot,xnuo, xnun2,xnuo2,ti,tn,deno,denn2,
c    . 			 deno2)
c     	routine by vincent wickwar, sri, may 1975.
c     	routine to find the ion-neutral collision frequencieS
c     	for no+ with o, n2, and o2.
c     	these are the *reduced* frequencies given by schunk and
c     	walker (1973), pss 21, page 1896.
c 
c ---- 	subroutine o2ion (xnutot,xnuo,xnun2,xnuo2,ti,tn,deno,denn2,
c    . 			 deno2)
c     	routine by vincent wickwar, sri, may 1975.
c     	routine to find the ion-neutral collision frequencies
c     	for o2 with o, n2 and o2.
c     	these are the *reduced* frequencies given by schunk and
c     	walker (1973), pss 21, page 1896.
c
c ---- 	subroutine compos(timeloc,nalt,alt,comp,z50)
c
c 	This subroutine computes the O+ rate [%] (i.e [O+]/Ne)
c 	using the "all season" coefficients fitted by Chantal
c 	Lathuillere (Personnal communication, 1991).
c 	These formulas have been calculated for a time which is
c 	local time + 1
c
c ----  subroutine gridcst (ntab,tabmin,tabmax,tab,dtab)
c       determine grid with equally spaced points
c
c ---- 	subroutine gridexp (ntab,tabmin,tabmax,tab)
c       determine grid with exponentially spaced points
c
c ----  subroutine gridpolo(ntab,tmin,tmax,tab,widthtab,spfac)
c       computes an array of ntab points between tmin and tmax, with
c       power law spacing
c       tab outputs are in increasing order
c       the growth factor is spfac
c
c----------------------------------------------------------------
c
      subroutine bfield(date,alt,la,lon,x,y,z,f,dip,dec)
c
c     cette sub. a pour objet, ici , de calculer B a chaque alt
c     igs 1975 magnetic field model				  *
c     input						          *
c     date	 =    epoch in years and decimals of the year 	  *
c 		      sous forme 19**				  *
c     itype	 = 1  geodetic coordinates			  *
c     itype	 = 2  geocentric coordinates			  *
c     alt	 =    height in km above sea level if itype = 1
c     alt	 =    distance from earths center in km if itype = 2
c     la	 =    latitude in degrees. positive north
c     lon	 =    longitude in degrees. positive east
c output
c     x	 =    north component of field in gammas
c     y	 =    east component of field in gammas
c     z	 =    vertical downward component of field in gammas
c     f	 =    total field in gammas
c     dip	 =    dip angle in radians. negative when z is positive
c     dec	 =    declination in radians. positive east
c observe that alt, la  and lon as well as x, y, and z are given
c in the coordinate system chosen by itype.
c
      real la,lon
c
      colat=90.-la
      elong=lon
      if(elong.lt.0.)elong=elong+360.
      itype=1
      call igrf(date,itype,alt,colat,elong,x,y,z,f)
      a=sqrt(x*x+y*y)
      dip=-atan(z/a)
      dec=atan2(y,x)
c
      return
      end
c
c---------------------------------------------------------------
c
      subroutine igrf (date,itype,alt,colat,elong,x,y,z,f)
c
c     cette sub. est appelee par bfield
c
      dimension agh(168),dgh(168),egh(168),p(90),q(90),cl(12),sl(12)
      dimension agh1(70),agh2(98)
c
      data agh1/-30103.6,-2016.5, 5682.6,-1906.7, 3009.9,-2064.7,1633.0,
     1 	     -58.1, 1278.2,-2142.0, -329.8, 1254.7,  265.9,  831.0,
     2 	    -227.0,  946.9,  792.5,  193.4,  443.8, -265.8, -403.9,
     3 	      53.0,  212.5, -285.2, -220.6,  351.4,	24.5,  262.3,
     4 	     148.4,  -63.8, -161.3, -157.5,  -83.4,  -40.2,   92.3,
     5 	      44.1,   69.9,  -11.2,   27.7,  100.4, -194.3,   77.6,
     6 	      -0.9,  -40.3,    3.8,   -7.9, -108.7,	15.6,   71.5,
     7 	     -53.3,  -76.6,    2.3,  -24.7,   13.4,	-4.5,   -6.4,
     8 	       7.0,    3.2,   24.5,   17.0,  -21.8,	-5.9,  -12.9,
     9 	      11.0,    5.1,    4.9,   -2.6,  -13.9,  -12.6,    5.0/
      data agh2/  -13.8,  -18.0,   -0.1,    5.7,   -2.4,14.5,   12.3,
     1 	     -11.1,    4.9,  -16.7,    9.3,   10.0,  -19.6,    1.6,
     2 	      15.7,  -11.4,    4.9,   10.6,   -3.1,	 0.6,   -4.2,
     3 	      -0.2,    9.7,    0.6,   12.2,    0.5,	-0.2,    0.5,
     4 	       0.3,   -5.0,   -3.3,    1.3,    2.4,	 2.0,   -6.0,
     5 	       2.6,   -1.4,    2.8,    6.6,   -3.6,	 4.6,   -0.2,
     6 	       1.2,    0.3,   -1.8,    3.2,    3.4,	 3.0,   -1.0,
     7 	      -3.4,    2.8,   -1.9,    0.5,   -4.5,	 0.7,    2.9,
     8 	      -0.9,   -1.2,   -1.5,    1.3,    0.3,	-0.8,    0.6,
     9 	       1.9,   -2.1,    3.4,    0.9,   -1.6,	-2.5,    1.7,
     * 	      -0.7,    2.5,    0.3,   -0.5,    0.8,	 0.1,   -1.3,
     1 	      -0.5,    0.1,    0.4,   -0.6,    0.0,	 0.0,    0.8,
     2 	      -1.8,   -0.1,   -1.6,   -0.2,   -0.7,	-0.5,   -0.9,
     3 	       0.3,    0.4,   -2.0,    0.0,    1.4,	 1.1,   -0.2/
      data dgh/ 26.8, 10.0,-10.1,-25.0,  0.3, -2.8,  5.5,-18.9, -3.8,
     1 	  -10.5,  7.2, -4.7,  2.8, -4.7, -6.4, -0.9, -2.2,  5.4,
     2 	   -4.0,  0.7, -2.1,  2.6, -4.6, -0.7,  0.2, -1.0,  0.9,
     3 	    1.3,  2.6, -2.1, -2.7, -0.6,  1.3,  1.3,	1.1,  0.6,
     4 	    0.9, -0.3,	2.3, -0.2,  3.5,  0.2,  0.0, -1.6,  0.8,
     5 	    0.4, -0.4,	2.0, -0.4, -0.2, -1.2, -0.5, -0.2,  0.3,
     6 	    0.0,  0.8,	0.3,  0.6, -0.6,  0.5,  0.0, -0.8,  1.2,
     7 	    0.4,  0.3, -0.2,  0.0, -0.3,  0.4, -0.3, -0.2, -0.3,
     8 	   -0.4,  0.5,	0.6, -0.5, -0.3, -0.6,  0.0,	0.5,88*0.0/
      data egh/ 0.70, 0.00,-0.49,-0.20, 0.00, 0.68, 0.16, 0.00,-0.28,
     1 	   0.00, 0.13,-0.32, 0.00, 0.00, 0.00, 0.00,-0.17, 0.30,
     2 	   0.00, 0.00,-0.14, 0.00,-0.17, 0.16,-0.13,-0.14,-0.14,
     3 	  -0.14, 0.10,-0.10, 0.00, 0.00, 0.00, 0.00, 0.10, 0.07,
     4 	   0.04, 0.00, 0.12, 0.00, 0.11,-0.10, 0.00, 0.00, 0.09,
     5 	   0.00, 0.00, 0.00,120*0.0/
      do 1000 i=1,70
 1000 agh(i)=agh1(i)
      do 2000 j=71,168
 2000 agh(j)=agh2(j-70)
      t     = date - 1975.0
      r     = alt
      one   = colat*0.0174533
      slat  = cos(one)
      clat  = sin(one)
      one   = elong*0.0174533
      cl(1) = cos(one)
      sl(1) = sin(one)
      x     = 0.0
      y     = 0.0
      z=0.
      cd    = 1.0
      sd    = 0.0
      l     = 1
      m     = 1
      n     = 0
      go to (1,2),itype
    1 a2    = 40680925.
      b2    = 40408585.
      one=a2*clat*clat
      two=b2*slat*slat
      three=one+two
      four=sqrt(three)
      r     = sqrt(alt*(alt + 2.0*four) + (a2*one + b2*two)/three)
      cd    = (alt + four)/r
      sd    = (a2 - b2)/four*slat*clat/r
      one   = slat
      slat  = slat*cd - clat*sd
      clat=clat*cd+one*sd
    2 ratio=6371.2/r
      p(1)  = 2.0*slat
      p(2)=2.0*clat
      p(3)  = 4.5*slat*slat - 1.5
      p(4)  = 5.1961524*clat*slat
      q(1)  = -clat
      q(2)  =  slat
      q(3)  = -3.0*clat*slat
      q(4)=1.7320508*(slat*slat-clat*clat)
      do 15 k=1,90
      if(n-m)3,4,4
    3 m     = 0
      n=n+1
      rr    = ratio**(n + 2)
      fn    = n
    4 fm    = m
      if (k-5) 8,5,5
    5 if (m-n) 7,6,7
   6  one   = sqrt(1.0 - 0.5/fm)
      j     = k - n - 1
      p(k)  = (1.0 + 1.0/fm)*one*clat*p(j)
      q(k)  = one*(clat*q(j) + slat/fm*p(j))
      sl(m) = sl(m-1)*cl(1) + cl(m-1)*sl(1)
      cl(m) = cl(m-1)*cl(1) - sl(m-1)*sl(1)
      go to 8
    7 one   = sqrt(fn*fn - fm*fm)
      two   = sqrt((fn - 1.0)**2 - fm*fm)/one
      three = (2.0*fn - 1.0)/one
      i     = k - n
      j     = k - 2*n + 1
      p(k)  = (fn + 1.0)*(three*slat/fn*p(i) - two/(fn - 1.0)*p(j))
      q(k)  = three*(slat*q(i) - clat/fn*p(i)) - two*q(j)
    8 one   = (agh(l)+(dgh(l)+0.5*egh(l)*t)*t)*rr
      if (m) 10,9,10
    9 x     = x + one*q(k)
      z     = z - one*p(k)
      l     = l + 1
      go to 14
   10 two   =(agh(l+1)+(dgh(l+1)+0.5*egh(l+1)*t)*t)*rr
      three = one*cl(m) + two*sl(m)
      x     = x + three*q(k)
      z     = z - three*p(k)
      if (clat) 12,12,11
   11 y     = y + (one*sl(m) - two*cl(m))*fm*p(k)/((fn + 1.0)*clat)
      go to 13
   12 y     = y + (one*sl(m) - two*cl(m))*q(k)*slat
   13 l     = l + 2
   14 m     = m + 1
   15 continue
      one   = x
      x     = x*cd +   z*sd
      z     = z*cd - one*sd
      f     = sqrt(x*x + y*y + z*z)
c
      return
      end
c
c-------------------------------------------------------------------
c
 	subroutine intlin (nin,xin,yin,nout,xout,yout)
c
c 	Subroutine d'interpolation lineaire.
c 	No restriction on extrema. Extrapolation is performed if 
c       necessary.
c 	L'ordre croissant des donnees d'entree est assure par pgmme
c
 	dimension xin(nin),yin(nin)
        dimension xout(nout),yout(nout)
c
 	dimension xxin(2048),yyin(2048)
        dimension xxout(2048),yyout(2048)
c
  	if(nin.gt.2048)then
 	  write(6,*)'Pgmme intlin'
	  write(6,*)'Taille du tableau d''entree trop grande'
 	  write(6,*)'Taille max = 2048. Taille actuelle =',nin
 	  write(6,*)'Programme arrete. Modifier 2048, recompiler'
 	  write(6,*)'et redemarer'
 	  stop
 	endif
  	if(nout.gt.2048)then
 	  write(6,*)'Pgmme intlin'
	  write(6,*)'Taille du tableau de sortie trop grande'
 	  write(6,*)'Taille max = 2048. Taille actuelle =',nout
 	  write(6,*)'Programme arrete. Modifier 2048, recompiler'
 	  write(6,*)'et redemarer'
 	  stop
 	endif

c 	Teste l'ordre des donnees
        iordre=0
c
 	if (xin(1).gt.xin(nin))then
	  do 10 i=1,nin
 	    xxin(i) = xin(nin+1-i)
 	    yyin(i) = yin(nin+1-i)
10 	  continue
 	else
	  do 20 i=1,nin
 	    xxin(i) = xin(i)
 	    yyin(i) = yin(i)
20 	  continue
 	endif
c
 	if (xout(1).gt.xout(nout))then
	  do 30 i=1,nout
 	    xxout(i) = xout(nout+1-i)
	    iordre=1
30 	  continue
 	else
	  do 40 i=1,nout
 	    xxout(i) = xout(i)
40 	  continue
 	endif
c
c 	interpolation lineaire.
c
 	do 60 iout = 1,nout
 	  do 50 iin =1,nin
 	    if(xxout(iout).lt.xxin(1))then
 	      yyout(iout)=(yyin(2)-
     .                   (xxin(2)-xxout(iout))*(yyin(2)-yyin(2-1))
     .                  /(xxin(2)-xxin(2-1)))
 	      go to 60
 	    elseif(xxout(iout).gt.xxin(nin))then
 	      yyout(iout)=(yyin(nin)-
     .                   (xxin(nin)-xxout(iout))*(yyin(nin)-yyin(nin-1))
     .                  /(xxin(nin)-xxin(nin-1)))
 	      go to 60
 	    elseif (xxout(iout).eq.xxin(iin))then
 	      yyout(iout)=(yyin(iin))
 	      go to 60
 	    elseif (xxout(iout).lt.xxin(iin))then
 	      yyout(iout)=(yyin(iin)-
     .                   (xxin(iin)-xxout(iout))*(yyin(iin)-yyin(iin-1))
     .                  /(xxin(iin)-xxin(iin-1)))
 	      go to 60
 	     endif
 50 	  continue
 60 	continue
c
	do 70 iout=1,nout
 	  if(iordre.eq.1)then
	    yout(iout) = yyout(nout+1-iout)
 	  else
	    yout(iout) = yyout(iout)
 	  endif
70 	continue
c
 	return
 	end
c
c------------------------- gaussint -------------------------------
c 
  	subroutine gaussint (ngau,nx,x,fx,sum)
c
c    	Computes |integrale fx(i) dx(i)| using gaussian integration.
c 	Each interval is splitted in ngauss points on which the 
c 	integration is performed. To speed up the program, the
c 	absissas and weight factors for ngauss=6, 12 or 24 are stored.
c 	Otherwise they are computed.
c 	After Handbook of Mathematical functions,
c 	Abramovitz and Stegun, 1970. jl 1990.
c
      	dimension x(nx),fx(nx),absc(100),wt(100)
1000 	format(10('-'),' Warning! Gaussian integration asked on a too',
     &    /,20x,'small (',i3,') number of points. Will be ',
     &    /,20x,'performed on 6 points.')
1010 	format(10('-'),' Warning! Gaussian integration asked on an odd',
     &    /,20x,'(',i3,') number of points. Will be performed ',
     &    /,20x,'on an even (',i3,') number of points.')
c
 	if(ngau.lt.6)then
	  ngauss=6
	  ngauss2=3
	  write(6,1000)ngau
	elseif(float(ngau/2)-float(ngau)/2. .ne. 0.)then
c 	  The integration is performed on an even number of points
c 	  to speed up and simplifie the program.
	  ngauss=ngau+1
	  ngauss2=ngauss/2
	  write(6,1010)ngau,ngauss
 	else
c 	  The gauss parameter being symetical, half of them is enough.
 	  ngauss=ngau
 	  ngauss2 = ngauss/2
 	endif
	do 5 i=1,ngauss2
	  absc(i)=0.
	  wt(i)  =0.
5 	continue
c
 	if(ngauss .eq. 6)then
 	  absc( 1) = 0.112701654434204101562
 	  absc( 2) = 0.500000000000000000000
 	  absc( 3) = 0.887298345565795898438
 	  wt  ( 1) = 0.277777791023254394531
 	  wt  ( 2) = 0.444444447755813598633
 	  wt  ( 3) = 0.277777791023254394531
 	elseif(ngauss .eq. 12)then
 	  absc( 1) = 0.033765256404876708984
 	  absc( 2) = 0.169395297765731811523
 	  absc( 3) = 0.380690395832061767578
 	  absc( 4) = 0.619309604167938232422
 	  absc( 5) = 0.830604672431945800781
 	  absc( 6) = 0.966234743595123291016
 	  wt  ( 1) = 0.085662245750427246094
 	  wt  ( 2) = 0.180380791425704956055
 	  wt  ( 3) = 0.233956962823867797852
 	  wt  ( 4) = 0.233956962823867797852
 	  wt  ( 5) = 0.180380791425704956055
 	  wt  ( 6) = 0.085662245750427246094
 	elseif(ngauss .eq. 24)then
 	  absc( 1) = 0.009219676256179809570
 	  absc( 2) = 0.047941386699676513672
 	  absc( 3) = 0.115048676729202270508
 	  absc( 4) = 0.206341028213500976562
 	  absc( 5) = 0.316084265708923339844
 	  absc( 6) = 0.437383294105529785156
 	  absc( 7) = 0.562616705894470214844
 	  absc( 8) = 0.683915734291076660156
 	  absc( 9) = 0.793658971786499023438
 	  absc(10) = 0.884951353073120117188
 	  absc(11) = 0.952058613300323486328
 	  absc(12) = 0.990780353546142578125
 	  wt  ( 1) = 0.023587668314576148987
 	  wt  ( 2) = 0.053469661623239517212
 	  wt  ( 3) = 0.080039165914058685303
 	  wt  ( 4) = 0.101583711802959442139
 	  wt  ( 5) = 0.116746269166469573975
 	  wt  ( 6) = 0.124573521316051483154
 	  wt  ( 7) = 0.124573521316051483154
 	  wt  ( 8) = 0.116746269166469573975
 	  wt  ( 9) = 0.101583711802959442139
 	  wt  (10) = 0.080039165914058685303
 	  wt  (11) = 0.053469661623239517212
 	  wt  (12) = 0.023587668314576148987
   	else
	  do i=1,ngauss2
	    absc(i)=0.
	    wt(i)=0.
	  enddo
	  call qgauss( ngauss2, absc, wt )
   	endif
c
 	sum = 0.
	do 20 ix=2,nx
	  A = (x(ix)-x(ix-1))/2.
	  B = (x(ix)+x(ix-1))/2.
	  fyy=(fx(ix)-fx(ix-1))/(x(ix)-x(ix-1))
	  do 10 i=1,ngauss2
	    y = A * absc(i) + B
	    fy=fyy*(y-x(ix))+fx(ix)
	    sum = sum + A*fy*wt(i)
	    
	    y = A * (-absc(i)) + B
	    fy=fyy*(y-x(ix))+fx(ix)
	    sum = sum + A*fy*wt(i)
10 	  continue
20 	continue
c
	  if(x(1).gt.x(nx))sum=-1*sum
c
      return
      end
c 
c----------------------- hint ------------------------------
c 
      subroutine hint(nin,z,fction,sum)
c
c    integrate using trapezoid rule
c 	sum=integral(fction(z).dz)
c
      dimension fction(*),z(*)
c
      ninm1 = nin - 1
      sum = 0.0
      do 1 k = 1,ninm1
        j = nin - k
        dz = z(j) - z(j+1)
        den = ( fction(j+1) + fction(j) ) * 0.5
        sum = sum + den * dz
    1 continue
      if (sum.lt.0)sum=-sum
c
      return
      end
c
c------------------------- mnmx --------------------------------------
c
      	subroutine mnmx(tab,ntab,tmin,tmax,linlog)
c
      	dimension tab(ntab)
c
 	if(linlog.eq.0)then
          tmin=tab(1)
          tmax=tab(1)
          do 1 i=2,ntab
	    tmax=max(tmax,tab(i))
	    tmin=min(tmin,tab(i))
1         continue
 	else
c 	  finds the first min non equal to zero, and the max of tab
          tmax=tab(1)
          do 20 i=2,ntab
	    tmax=max(tmax,tab(i))
20         continue
 	  if(tmax.le.0.)then
 	    write(6,*)'Max <= 0, dessin log impossible'
 	    go to  60
   	  endif
 	  do 30 i=1,ntab
 	    if (tab(i).gt.0.)then
 	      tmin=tab(i)
 	      go to 40
 	    endif
 30 	  continue
 40 	  continue
	  do 50 i=1,ntab
 	    if (tab(i).lt.tmin.and.tab(i).gt.0.) tmin=tab(i)
 50 	  continue
 60 	  continue
 	endif
c
          return
          end
c
c-------------------------- mnmxplt -----------------------------------
c
	subroutine mnmxplt(tmin,tmax,linlog)
c
c       find best min and max for nice plot.
c 	linlog = 0 if linear axis
c 	       = 1 if logarithmic axis
c
	if(linlog.eq.1)then
c 	  Quand on est en axe log et qu'il y a moins d'une decade entre
c 	  le min et le max, GREG ne marque pas les unites sur l'axe >
c 	  On se premuni de cela.
 	  if(tmin.ne.0.)then
 	    ttmin = log10(tmin)
 	    ttmax = log10(tmax)
 	    delta = ttmax-ttmin
 	    if(delta .lt. 1.)then
 	      ttminint = float(ifix(ttmin))
 	      ttmaxint = float(ifix(ttmax)+1)
c 	      Quel est le plus proche d'un tick?
 	      deltmin = ttmin-ttminint
 	      deltmax = ttmaxint-ttmax
 	      if(deltmin.le.deltmax)then
c 		C'est tmin!
 	        tmin = 10**ttminint
 	        tmax = tmax*1.2
 	      else
c 		C'est tmax!
 	 	tmax = 10**ttmaxint
 		tmin = tmin/1.2
 	      endif
 	      return
    	    endif
	  endif
c
	  tmin=tmin/1.2
 	  if(tmin.eq.0.)tmin=1.e-05
	  tmax=tmax*1.2
 	else
	  if (tmin.lt.200.) then
	    tmin = float(ifix(tmin)/10 - 1) * 10.
	  elseif (tmin.lt.1000.) then
	    tmin = float(ifix(tmin)/10 - 5) * 10.
	  else
	    tmin = float(ifix(tmin)/10 -10) * 10.
 	  endif
	  if (tmax.lt.200.) then
	    tmax = float(ifix(tmax)/10 + 2) * 10.
	  elseif (tmax.lt.1000.) then
	    tmax = float(ifix(tmax)/10 + 5) * 10.
	  else
	    tmax = float(ifix(tmax)/10 +10) * 10.
 	  endif
 	endif
c
	return
	end
c
c ----------------------- mnmxi --------------------------------
c
	  subroutine mnmxi(tab,ntab,imax,imin)

c 	renvoie les indices du max et min du tableau tab.
      	dimension tab(ntab)
c
          tmin=tab(1)
          tmax=tab(1)
          do i=1,ntab
	    tmax=max(tmax,tab(i))
	    tmin=min(tmin,tab(i))
 	  enddo
          do i=1,ntab
	    if(tmax.eq.tab(i))imax=i
	    if(tmin.eq.tab(i))imin=i
 	  enddo
c
          return
          end
c
c ------------------------- qgauss ------------------------------
c
      subroutine qgauss(m,gmu,gwt)
C
C       COMPUTE WEIGHTS AND ABSCISSAE FOR ORDINARY GAUSSIAN QUADRATURE
C       (NO WEIGHT FUNCTION INSIDE INTEGRAL) ON THE INTERVAL (0,1)
C
C       REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
C                   Integration, Academic Press, New York, pp. 87, 1975.
C
C          METHOD:  Compute the abscissae as roots of the Legendre
C                   Polynomial P-SUB-N using a cubically convergent
C                   refinement of Newton's method.  Compute the
C                   weights from EQ. 2.7.3.8 of Davis/Rabinowitz.
C
C        ACCURACY:  at least 13 significant digits
C
C
C  I N P U T :    M       ORDER OF QUADRATURE RULE
C
C  O U T P U T :  GMU(I)  I = 1 TO M,    ARRAY OF ABSCISSAE
C                 GWT(I)  I = 1 TO M,    ARRAY OF WEIGHTS
C
C  I N T E R N A L    V A R I A B L E S:
C
C    PM2,PM1,P : 3 SUCCESSIVE LEGENDRE POLYNOMIALS
C    PPR       : DERIVATIVE OF LEGENDRE POLYNOMIAL
C    P2PRI     : 2ND DERIVATIVE OF LEGENDRE POLYNOMIAL
C    TOL       : CONVERGENCE CRITERION FOR LEGENDRE POLY ROOT ITERATION
C    X,XI      : SUCCESSIVE ITERATES IN CUBICALLY-
C                CONVERGENT VERSION OF NEWTON'S METHOD
C                ( SEEKING ROOTS OF LEGENDRE POLYNOMIAL )
C+---------------------------------------------------------------------+
      REAL     CONA, GMU(*), GWT(*), PI, T
      INTEGER  LIM, M, NP1
c     DOUBLE   PRECISION  D1MACH
      DOUBLE   PRECISION  EN, NNP1, P, PM1, PM2, PPR, P2PRI, PROD,
     $                    TMP, TOL, X, XI
c     DATA     PI / 0.0 /
      PI = 0.0
C
C
      IF ( PI.EQ.0.0 )  THEN
         PI = 2. * ASIN(1.0)
c 	 La precision indiquee ci apres est heuristique (au dessus
c 	 de m=10, le pgmme plante parfois pour des prec. < 1.e-15)
         TOL = 1.d-30
	 if (M.gt.1) tol=1.d-15
c 	 d1mach est une fonction cray qui donne la precision.
c        TOL = 10. * D1MACH(3)
      END IF
C
      IF ( M.LE.1 )  THEN
         M = 1
         GMU( 1 ) = 0.5
         GWT( 1 ) = 1.0
         RETURN
      END IF
C
      EN   = M
      NP1  = M + 1
      NNP1 = M * NP1
      CONA = FLOAT( M-1 ) / ( 8 * M**3 )
C     INITIAL GUESS FOR K-TH ROOT OF LEGENDRE POLYNOMIAL, FROM
C     DAVIS/RABINOWITZ (2.7.3.3A)
      LIM  = M / 2
      DO 30  K = 1, LIM
         T = ( 4*K - 1 ) * PI / ( 4*M + 2 )
         X = COS ( T + CONA / TAN( T ) )
C        RECURSION RELATION FOR LEGENDRE POLYNOMIALS
10       PM2 = 1.D0
         PM1 = X
         DO 20 NN = 2, M
            P   = ( ( 2*NN - 1 ) * X * PM1 - ( NN-1 ) * PM2 ) / NN
            PM2 = PM1
            PM1 = P
20       CONTINUE
C
         TMP   = 1.D0 / ( 1.D0 - X**2 )
         PPR   = EN * ( PM2 - X * P ) * TMP
         P2PRI = ( 2.D0 * X * PPR - NNP1 * P ) * TMP
         XI    = X - ( P / PPR ) * ( 1.D0 +
     .               ( P / PPR ) * P2PRI / ( 2.D0 * PPR ) )
C
C        CHECK FOR CONVERGENCE
         IF ( DABS(XI-X) .GT. TOL ) THEN
            X = XI
            GO TO 10
         END IF
C        ITERATION FINISHED--CALC. WEIGHTS, ABSCISSAE FOR (-1,1)
         GMU( K ) = - X
         GWT( K ) = 2.D0 / ( TMP * ( EN * PM2 )**2 )
         GMU( NP1 - K ) = - GMU( K )
         GWT( NP1 - K ) =   GWT( K )
30    CONTINUE
C     SET MIDDLE ABSCISSA AND WEIGHT FOR RULES OF ODD ORDER
      IF ( MOD( M,2 ) .NE. 0 )  THEN
         GMU( LIM + 1 ) = 0.0
         PROD = 1.D0
         DO 40 K = 3, M, 2
            PROD = PROD * K / ( K-1 )
40       CONTINUE
         GWT( LIM + 1 ) = 2.D0 / PROD**2
      END IF
C     CONVERT FROM (-1,1) TO (0,1)
      DO 50  K = 1, M
         GMU( K ) = 0.5 * GMU( K ) + 0.5
         GWT( K ) = 0.5 * GWT( K )
50    CONTINUE
C
      RETURN
      END
c
c------------------------------------------------------------------
c
      subroutine xline(nbline,ijfile)                   
      implicit none
      integer,intent(in) :: nbline,ijfile  

      character nc*30
      integer i

      do i=1,nbline        
        read(ijfile,1000) nc    
      end do

 1000 format(a1)              
              
      end subroutine xline                   
c
c---------------------------------------------------------------
c
      SUBROUTINE NUOION (XNUTOT,XNUO,XNUN2,XNUO2, TI,TN,DENO,DENN2,
     &                   DENO2,burnside)
      implicit None

      Real,Intent(In) :: TI, TN,deno,denn2,deno2,burnside
      Real,Intent(out):: XNUTOT,XNUO,XNUN2,XNUO2

      Real TM,XLOG
C
C     ROUTINE BY VINCENT WICKWAR, SRI, MAY 1975.
C
C     ROUTINE TO FIND THE ION-NEUTRAL COLLISION FREQUENCIES
C     FOR O+, WITH O, N2, AND O2.
C     THESE ARE THE *REDUCED* FREQUENCIES adapted from SCHUNK
C     AND WALKER (1973), SCHUNK AND NAGY (1980)
C
C     OUTPUTS:
C		   ION-NEUTRAL COLLISION FREQUENCIES FOR O+
C	  XNUTOT TOTAL (COLLISIONS/SEC)
C	  XNUO   WITH O (COLLISIONS/SEC)
C	  XNUN2  WITH N2 (COLLISIONS/SEC)
C	  XNUO2  WITH O2 (COLLISIONS/SEC)

C     INPUTS:
C	  TI	   ION TEMPERATURE (DEG)
C	  TN	   NEUTRAL TEMPERATURE (DEG)
C	  DENO   DENSITY OF ATOMIC OXYGEN (NB/CM**3)
C	  DENN2  DENSITY OF MOLECULAR NITROGEN (NB/CM**3)
C	  DENO2  DENSITY OF MOLECULAR OXYGEN (NB/CM**3)
C
!      ALOG10(X) = 0.43429448 * ALOG(X)
C
      TM = (TI + TN) / 2.0
      XLOG = LOG10(TM)

c       Formule de Schunk et Walker (1973), abandonnee pour Salah.
c     XNUO = 3.67E-11 * DENO * SQRT(TM) *
c    *   (1.0 - 0.064 * XLOG)**2
c	on remplace cette valeur par celle recommande par le CEDAR
c	qui inclut le facteur de Burnside:1.7 (Chantal, 1993)
c 	On s'autorise quand meme de pouvoir jouer sur ce facteur 
c 	"burnside" devient un parametre libre (jl, 1994)
c	equation 5, Salah, GRL 93, vol 20, p1543-1546

      xnuo =2.3e-11 * deno * tm**.5 *burnside
      XNUN2 = 6.82E-10 * DENN2
      XNUO2 = 6.64E-10 * DENO2
      XNUTOT = XNUO + XNUN2 + XNUO2

      END SUBROUTINE NUOION
C    
c----------------------------------------------------------------------
c
      	SUBROUTINE NOION (XNUTOT,XNUO, XNUN2,XNUO2,DENO,DENN2,
     & 			 DENO2)
      implicit None

      Real,Intent(In) :: deno,denn2,deno2
      Real,Intent(out):: XNUTOT,XNUO,XNUN2,XNUO2
C
C     	ROUTINE BY VINCENT WICKWAR, SRI, MAY 1975.
C
C     	ROUTINE TO FIND THE ION-NEUTRAL COLLISION FREQUENCIES
C     	FOR NO+ WITH O, N2, AND O2.
C     	THESE ARE THE *REDUCED* FREQUENCIES GIVEN BY SCHUNK AND
C     	WALKER (1973), PSS 21, PAGE 1896.
C
C     OUTPUTS:
C		   ION-NEUTRAL COLLISION FREQUENCIES FOR NO+
C	  XNUTOT TOTAL (COLLISIONS/SEC)
C	  XNUO   WITH O (COLLISIONS/SEC)
C	  XNUN2  WITH N2 (COLLISIONS/SEC)
C	  XNUO2  WITH O2 (COLLISIONS/SEC)

C     INPUTS:
C	  DENO   DENSITY OF ATOMIC OXYGEN (NB/CM**3)
C	  DENN2  DENSITY OF MOLECULAR NITROGEN (NB/CM**3)
C	  DENO2  DENSITY OF MOLECULAR OXYGEN (NB/CM**3)
C
      XNUO = 2.44E-10 * DENO
      XNUN2 = 4.34E-10 * DENN2
      XNUO2 = 4.27E-10 * DENO2
      XNUTOT = XNUO + XNUN2 + XNUO2

      END SUBROUTINE NOION
C    
c----------------------------------------------------------------------
c
      	SUBROUTINE O2ION (XNUTOT,XNUO,XNUN2,XNUO2,TI,TN,DENO,DENN2,
     . 			 DENO2)
C
C     	ROUTINE BY VINCENT WICKWAR, SRI, MAY 1975.
C
C     	ROUTINE TO FIND THE ION-NEUTRAL COLLISION FREQUENCIES
C     	FOR O2 WITH O, N2 AND O2.
C     	THESE ARE THE *REDUCED* FREQUENCIES GIVEN BY SCHUNK AND
C     	WALKER (1973), PSS 21, PAGE 1896.
C
C     OUTPUTS:
C		   ION-NEUTRAL COLLISION FREQUENCIES FOR O2
C	  XNUTOT TOTAL (COLLISIONS/SEC)
C	  XNUO   WITH O (COLLISIONS/SEC)
C	  XNUN2  WITH N2 (COLLISIONS/SEC)
C	  XNUO2  WITH O2 (COLLISIONS/SEC)
C     INPUTS:
C	  TI	   ION TEMPERATURE (DEG)
C	  TN	   NEUTRAL TEMPERATURE (DEG)
C	  DENO   DENSITY OF ATOMIC OXYGEN (NB/CM**3)
C	  DENN2  DENSITY OF MOLECULAR NITROGEN (NB/CM**3)
C	  DENO2  DENSITY OF MOLECULAR OXYGEN (NB/CM**3)
C
!      ALOG10(X) = 0.43429448 * ALOG(X)
C
      TM = (TI + TN) / 2.0
      XNUO = 2.31E-10 * DENO
      XNUN2 = 4.13E-10 * DENN2
      IF ( TM - 800.0 ) 1, 1, 3
    1 XNUO2 = 4.08E-10 * DENO2
      GO TO 5
    3 XNUO2 = 2.59E-11 * DENO2 * SQRT(TM) * 
     *   (1.0 - 0.073 * XLOG)**2
    5 XNUTOT = XNUO + XNUN2 + XNUO2

      END SUBROUTINE O2ION
c
c---------------------------- compos -------------------------------
c
      subroutine compos(timeloc,nalt,alt,comp,z50)
c
c 	This subroutine computes the O+ rate [%] (i.e [O+]/Ne)
c 	using the "all season" coefficients fitted by Chantal
c 	Lathuillere (Personnal communication, 1991).
c 	These formulas have been calculated for a time which is
c 	local time + 1
c
      integer nalt
      real*4 comp(nalt),alt(nalt)
      real*4 timeloc,z50
c
      a = 218.3
      b = 8.2
      c = -117.1
      d = 108.9
      e = -26.2
      phitrom=19.21
c
c 	Computes the altitude at which [O+]=50%
c
        hr = timeloc-phitrom/15.
        x = hr/12.
        x2= x*x
        x3= x2*x
        x4= x3*x
        z50 = a + b*x + c*x2+ d*x3 + e*x4
c
c 	computes the standard deviation
c
        deltaz = 0.72 * z50 - 104.
c
c 	computes [O+] rate.
c
        do ialt = 1,nalt
          comp(ialt) = 0.5*(1.+ tanh((alt(ialt)-z50)/deltaz))
         if(comp(ialt).lt.0.) comp(ialt)=0.
        end do
c

        end subroutine compos
c
c----------------------------------------------------------------
c
      subroutine intquad(nin,xin,yin,nout,xout,yout)
c
c    perform quadratic interpolations on array yin(i) vs xin(i)
c    result is yout(i) at the x points given by xout(i)
c
      dimension yin(*),xin(*),xout(*),yout(*)
c
      do i=1,nout
        nl=2
        nu=nin-1
        x=xout(i)
        if((x-xin(nl))*(xin(nu)-x))1,1,2
    2   if(nu.le.nl+1)goto 1
        fint=(nl+1.5)*(xin(nu)-x)+(nu-.5)*(x-xin(nl))
        nm=fint/(xin(nu)-xin(nl))
        if((x-xin(nl))*(xin(nm)-x))4,3,3
    3   nu=nm
        go to 2
    4   nl=nm
        go to 2
    1   nm=nu
        if((x-xin(nl))**2.lt.(x-xin(nu))**2)nm=nl
        yout(i)=yin(nm)*(x-xin(nm-1))*(x-xin(nm+1))/
     .          ((xin(nm)-xin(nm-1))*(xin(nm)-xin(nm+1)))+
     .          (yin(nm-1)*(x-xin(nm+1))/(xin(nm-1)-xin(nm))+
     .          yin(nm+1)*(x-xin(nm-1))/ (xin(nm)-xin(nm+1)))*
     .          (x-xin(nm))/(xin(nm-1)-xin(nm+1))
      end do

      end subroutine intquad
c
c----------------------------------------------------------------------
c
      function lenc(c)
c
c Returns the location of the last non-blank character in a string.
c Arguments :
c	C	C*(*)	Character string			Input
c
      CHARACTER*(*) C
      INTEGER LENC,I
*
      DO I=LEN(C),1,-1
         IF (ICHAR(C(I:I)).GT.32) THEN
            LENC = I
            RETURN
         ENDIF
      ENDDO
      LENC = 0
      END function lenc
c
c----------------------------------------------------------------------
c
 	subroutine clasdesc(tab,ntab,index)
c
c 	Classe un tableau en ordre descendant. 
c 	tab : en entree, tableau a classe
c 	      en sortie, tableau classe (et donc original ecrase)
c 	ntab : nbre de points a classer
c 	index : tableau d'entiers contenant en sortie l'ordre de 
c 	        classement.
 	implicit none
c
 	integer ntab,itab
 	real tab(ntab),trav(3000),tmax,tmin
 	integer index(ntab)
 	integer imax,imin,jmin
c
 	if(ntab.gt.3000) then
 	  write(6,*)'La dimension max. du tableau a classer est 3000'
 	  write(6,*)'dans le sous programme clasdesc'
 	  write(6,*)'La dimension actuelle du tableau est',ntab
 	  write(6,*)'Augmenter la capacite de clasdesc, puis'
 	  write(6,*)'recompiler et relancer'
 	  stop
 	endif
c
 	do itab = 1,ntab
 	  trav(itab) = tab(itab)
 	enddo
c
 	call mnmx(trav,ntab,tmin,tmax,0)
 	do itab = 1,ntab
 	  call mnmxi(trav,ntab,imax,imin)
 	  index(itab) = imax
 	  trav(imax) = tmin - 1.
 	enddo
c
 	do itab = 1,ntab
 	  trav(itab) = tab(itab)
 	enddo
 	do itab = 1,ntab
 	  tab(itab) = trav(index(itab))
 	enddo
c

 	end subroutine clasdesc
c
c----------------------- gridcst ---------------------------------
c 
        subroutine gridcst (ntab,tabmin,tabmax,tab,dtab)
c
c       determine grid with equally spaced points
c 
        real tab(*),dtab(*)
c
        dzo = (tabmax-tabmin)/float(ntab-1)
        tab(1)=tabmax
        dtab(1) = dzo
        do i = 2,ntab-1
          tab(i)= tab(i-1) - dzo
          dtab(i) = dzo
        enddo
        tab(ntab)=tabmin
        dtab(ntab) = dzo
c
        return
        end
c
c----------------------- gridexp ---------------------------------
c 
      	subroutine gridexp (ntab,tabmin,tabmax,tab)
c
c       determine grid with exponentially spaced points
c
      	dimension tab(ntab)
c
      eps=1.0e-04
c
c     dzo = grid spacing between lowest points.
      dzo = (tabmax-tabmin)/(2.*float(ntab))
      dzo = min(1.,dzo)
      ntabm1 = ntab - 1
      flntabm1 = float(ntabm1)
      hsave = 0.0
      itr = 0
   10 itr = itr + 1
      hnu =flntabm1*dzo/log((tabmax + hsave)/(tabmin + hsave)) - tabmin
      resid = abs(hsave - hnu)
      scal = min(abs(hsave) , abs(hnu))
      hsave = hnu
      if(scal .eq. 0.0)go to 15
      resid = resid/scal
   15 if(resid .gt. eps)go to 10
      ratio = (tabmax + hsave)/(tabmin + hsave)
      alrat = log(ratio)
      do 20 i = 1,ntab
        tab( ntab - i + 1) = 
     .         (tabmin + hsave)*ratio**(float(i-1)/flntabm1) - hsave
   20 continue
      tab(1)=tabmax
      tab(ntab)=tabmin
c
      return
      end
c
c------------------------- gridpolo -----------------------------------
c
        subroutine gridpolo(ntab,tmin,tmax,tab,widthtab,spfac)
c
c       computes an array of ntab points between tmin and tmax, with
c       power law spacing
c       tab outputs are in increasing order
c       the growth factor is spfac

        integer,intent(in) :: ntab
        real,intent(out) :: tab(ntab),widthtab(ntab),spfac
        real,intent(in) :: tmin,tmax

c
        tab(1)=tmin
        tab(ntab) = tmax
        flntab =float(ntab)
c
c       Premiere estimation de l'ordre
        a = (tab(ntab)-tab(1))/(2.*tab(1))
        x = 2./((flntab-2.)*(flntab-1.))
        x = x*(a-flntab+1)
        spfac = 1.+x
c
c       Precision sur tmax de 0.01%
        eps = 0.0001
c
c       Ajustement
        iflag = 1
        niter = 0
        do while (iflag.ne.0)
          niter =niter +1
c         write(6,*)'niter ',niter,'[A'
          esup = entab(ntab,spfac,tmin)
          if(esup.gt.tmax*(1.+eps))then
c           Il faut reduire spfac
            x = x/1.2
            spfac = 1.+x
          elseif (esup.lt.tmax*(1.-eps))then
c           Il faut augmenter spfac
            x = x*1.15
            spfac = 1.+x
          else
            iflag = 0
          endif
        enddo
c
*       set up the energy-grid
        de=2.0*tab(1)
        do ien=2,ntab-1
          tab(ien)=tab(ien-1)+de
          de=de*spfac
        enddo
        tab(ntab) = tmax
c
c       Eventually computes the energy width.
        dd = tab(1)
        widthtab(1) = 2.*dd
        do i = 2,ntab
          ener = tab(i-1)+dd
          dd = tab(i)-ener
          widthtab(i) = 2.*dd
        enddo

        end subroutine gridpolo
c
c----------------------------------------------------------------------
c
        pure real function entab(ntab,spfac,tmin)
        implicit none
        real,intent(in) :: spfac,tmin
        integer,intent(in) :: ntab
        real x,summ
        integer i
c
        x = 1.
        summ = 0.
        do i =1,ntab-1
          summ = summ + x
          x = x*spfac
        enddo

        entab = tmin*(1.+2.*summ)

        end function entab
