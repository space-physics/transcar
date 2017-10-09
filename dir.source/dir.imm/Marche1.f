

c=====================================================================
c=====================================================================

      SUBROUTINE MARCHE

      integer*4 par1,par2,par3,dpar1,dpar2,dpar3
      integer*4 dpar11,dpar22,dpar33
      integer*4 dn1,dn2,dnt,dndl
      parameter (par1=1240,par2=2400,par3=2481)
      parameter (dpar1=(25*80),dpar2=(24*160),dpar3=(25*160+1))
      parameter (dpar11=(49*80),dpar22=(48*160),dpar33=(49*160+1))

      real*4 b1(par1),b2(par1),b3(par1)
     %,b4(par1)
     %,b5(par1),b6(par1),xu(par2),xv(par2),
     %PHIi(par1),ui(par1)
     
     %,dxdt(par2),dydt(par2),Tic(par1)
     %,Rnic(par1),hh(par1),
     %aire1(par1),RNNN(par1),TV(par1),ppi(par1)
     %,dppia(par1),dppib(par1),
     %hh1(par1),hh2(par1)     
     %,betaai(par1),Rni1(par1),Ti1(par1)
     %,u4(par1),Rni11(par1),Ti11(par1),
     %utemp(dpar1),FEtemp(dpar1),Etemp(dpar1),alitem(dpar1)

      real*4 RNNNtemp(dpar1),TVtemp(dpar1)
      integer*4 inum(40)
      integer*4 contro(6)
      real*4 TV1(par1),RNNN1(par1)
      real*4 Nmax,mass,u0(par1),T0(par1),rN0(par1)
      common/sol/ener2(par1),ane2(par1)
      dimension ane3(par1),ener3(par1)   
      dimension u(par1)
      dimension E(par1),F(par1),FE(par1)
      dimension result(6,168)
      common/donnees/el,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,
     $hB,n1,n2,Nmax
      common/triang1/n21,n22,n23,h21,h22,h23,Aint1,Aint2,raison
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      dimension temp(10),xn(10),yn(10)
      dimension ntria(10)
      real*4 dTtemp(25)
      real*4 dTtem2(49)

      real*4 alifit(par1),rafou(par1,0:5),iafou(par1,0:5)
      real*4 alii(par1)
      integer*4 par11

	data dTtemp/71.97,69.668,67.279,64.797,62.221,59.548,
     %  56.783,53.932,51.004,48.013,44.977,41.916,38.854,35.816,
     %  32.828,29.916,27.104,24.414,21.864,19.471,17.246,15.197,
     %  13.329,11.645,10.14/

	data dTtem2/90.,88.12,86.24,84.33,82.40,80.43,78.41,
     %  76.33,74.19,71.97,69.668,67.279,64.797,62.221,59.548,
     %  56.783,53.932,51.004,48.013,44.977,41.916,38.854,35.816,
     %  32.828,29.916,27.104,24.414,21.864,19.471,17.246,15.197,
     %  13.329,11.645,10.14,8.82,7.672,6.688,5.859,5.167,4.594,
     %  4.117,3.712,3.351,3.013,2.671,2.304,1.878,1.327,0./

       common/debut/ndl,nt,dn1,dn2,dndl,dnt,ie,istati,
     %iinti,ibrupt,fheure,Bfm,Bfs,rnmax,rtmax,
     %A(par1),B(par1),T(par1),m1(3,par2),m2(par1),
     %Atemp(dpar1),Btemp(dpar1),Ttemp(dpar1),m1temp(3,dpar2),
     %m2temp(dpar1),
     %tabTM(par1,4),tabMT(dpar1,4)

      common/transm/utem2(dpar11),FEtem3(dpar11),Etem3(dpar11),
     %alite2(dpar11)

      common/pot/FEtem2(dpar11),Etem2(dpar11)

      common/resul1/ane1(par1),ener1(par1),Ti(par1),Rni(par1),
     %ali(par1),u3(par1),uo(par1)

      common/debu1/phiopt,muns(dpar3),sp0temp(dpar1),
     %sh0temp(dpar1),ndimat,al1(9),al2(9),al3(9),al4(9),
     %Rnitemp(dpar1),Titemp(dpar1)

cHAO97

	common/fluxioni/FEion(dpar11),Fion(dpar11)
	
cHAO97

cesr

	common/dea/ideaele,ideaion,ideapot

cesr

c     ========================================
c
c     ENTREE DE LA TRIANGULATION : LECTURE SUR LE FICHIER 01 
c             ndl  : nombre de degres de libertes
c             nt   : nombre de triangles
c             n1   : nombre de subdivisions suivant la longitude
c             n2   : nombre de subdivisions suivant la colatitude
c             T    : coordonnees suivant la colatitude:T=teta
c             B    : coordonnees suivant la longitude:B=phi
c             A    : coordonnees suivant la colatitude:A=sin(T)**2
c             hA   : 1/n1
c             hB   : 1/n2
c             Bmax : valeur maximale de la longitude
c             Bmin : valeur minimale de la longitude
c             Amax : valeur maximale de la colatitude
c             Amin : valeur minimale de la colatitude
c             m1   : triangulation du domaine
c             m2   : pointeur des points frontieres
c            


c     ========================================      
c     ndl1=2481
c     n1 DESIGNE LE NOMBRE DE SUBDIVISIONS DE L AXE DES PHI:[0,2*pi]
c     n2 DESIGNE LE NOMBRE DE SUBDIVISIONS DE L AXE DES
c     TETA:[TETA0,pi/2-eps]
c         subdivisions uniformes en sin(t)**2
c     ========================================
c      n1=80
c      n2=30

C ************************************************************ DEBUT
	
	if (ie.eq.0)						then

C ************************************************* D'INITIALISATION

c!!!	write(*,*)'INITIALISATION IN MARCHE'

	ndimat=400000
	n1=40
	n2=30
	iinti=0
	ibrupt=1

CC	C'EST LA OU ON INDIQUE L'HEURE OU S'ARRETE LE RUN.

cesr	fheure=15000.

	Bfm=6.
	Bfs=18.

CC	ON FIXE LA DENSITE DE LA SOURCE

cesr	rnmax=0.4e+06

CC	ON FIXE LA TEMPERATURE DE LA SOURCE DES IONS (rtmax)
CC	ET DES ELECTRONS (emax)

cesr	rtmax=5.
cesr	emax=2.

CC	ON FIXE LA VALEUR DE LA DDP EN PRENANT LA MOITIE

cesr	phiopt=15.

CC	ON COMMENCE LES CALCULS

	rpmax=rtmax*rnmax
	
	heuref=3600.

      dn1=80
      dn2=24
      dnt=2*dn1*dn2
      dndl=dn1*(dn2+1)
      
      pi=4.*atan(1.)
      Bfm=Bfm*2.*pi/24.
      Bfs=Bfs*2.*pi/24.
      ndl=n1*(n2+1)
      nt=2*n1*n2
      nf=2*n1
c     ========================================
c     CONSTRUCTION D UNE TRIANGULATION UNIFORME SUR LE DISQUE
c     ========================================
      call disc1(dn2,dn1,m1temp,m2temp,dnt,0)
      call disc1(n2,n1,m1,m2,nt,0)

c     write(*,*) m1(1,639),m1(2,639),m1(3,639)
c!!!      print 1000,n1,n2,ndl,nt,nf
1000   format ('nombre de subdivisions de l axe des phi:',i3,/,
     &        'nombre de subdivisions de l axe des Teta:',i3,/,
     &        'nombre de degres de liberte',i4,/,
     &        'nombre de triangles:',i4,/,
     &        'nombre de noeuds sur la frontiere:',i3)
     
c!!!      write(*,*) ('ATTENTION: COS(X)**(1/2)')
      
c     ========================================
c     ========================================
      ioptec=0
      if(ioptec.lt.1) goto 502
     
      write(4,1001) ((m1(i,j),j=1,3),i=1,nt)
1001   format ('TRIANGULATION DU DISQUE',/,5(3i5,1x))
      write(4,1002) (m2(i),i=1,ndl)
1002   format ('CONDITIONS AUX LIMITES',/,2480(i1,1x))

c     ========================================
c     CALCUL DES COORDONNEES
c     ========================================
	do 9876 i=1,ndl
	u0(i)=0.
	T0(i)=0.
	rN0(i)=0.
9876	continue	
502     r=6.37e+06

      b0=3.08e-5
      Bmin=0.
      Bmax=2*pi
      Tmin=(90.-71.97)*pi/180.
      teta0=Tmin

      almax=1./cos(10.14*pi/180.)
      almax=almax*almax

      Amax=1./almax
      Tmax=asin(sqrt(Amax))
c!!!      print *,'Tmax=',Tmax
      hB=(Bmax-Bmin)/n1
      Tint1=30./((90.-10.14)-(90.-71.97))*(Tmax-Tmin)
      Tint2=50./((90.-10.14)-(90.-71.97))*(Tmax-Tmin)
      Amin=sin(Tmin)**2
      Aint1=Amax
      Aint2=Amax
      n21=n2
      n22=0
      n23=0
c     triangulation uniforme
      raison=0
c     on impose une triangulation uniforme et donc on a raison=0
c     tout le temps, ce qui rend inutile la lecture de raison.
      h21=(Aint1-Amin)/n21
      
c 	On effectue une triangulation reguliere en L.

c 	Calcul du L minimal.

      almin=(1./(sin((90.-71.97)*pi/180.)**2.))

c 	Calcul du pas en L.

      rlpas=(-almax+almin)/n21
      
c 	Calcul du pas en sin(TETA)**2. pour la frontiere polaire.

      hpas1=almin-rlpas
      hpas1=-(1./almin)+(1./hpas1)
      
c 	Calcul du pas en sin(TETA)**2. pour la frontiere equatoriale.

	hpas2=almax+rlpas
	hpas2=(1./almax)-(1./hpas2)
	
      h22=0.
      hA=h21
      do 561 j=1,n2+1
      if(j.lt.(n21+2)) then
      		rlcour=rlcour-rlpas
		Acour=(1./rlcour)
			endif
      if(j.eq.1) then
      		rlcour=almin
      		Acour=(1./rlcour)
		  endif
      do 566 k=1,n1
      i=(j-1)*n1+k
      A(i)=Acour
566     B(i)=(k-1)*hB
561     continue
      do 55 i=1,ndl
      T(i)=asin(sqrt(A(i)))
55     continue
      Amin=A(1)
      Amax=A(ndl)

      dBmin=0.
      dBmax=2*pi
      dhB=(dBmax-dBmin)/dn1

      do 1561 j=1,dn2+1
      Tcour=(90.-dTtemp(j))*pi/180.
      do 1566 k=1,dn1
      i=(j-1)*dn1+k
      Ttemp(i)=Tcour
1566  Btemp(i)=(k-1)*dhB
1561  continue

      do 1055 i=1,dndl
      Atemp(i)=(sin(Ttemp(i)))**(2.)
1055  continue
      dAmin=Atemp(1)
      dAmax=Atemp(dndl)

C CALCUL DES ELEMENTS DE PASSAGE ENTRE LA GRILLE DU TIEGCM
C ET CELLE DE NOTRE MODELE:
C
C	-) tabTM contient les elements de passage de la grille
C	-) du Tiegcm VERS celle de Notre Modele
C	-) tabMT contient les elements de passage de la grille
C	-) de Notre Modele VERS celle du Tiegcm

c!!!	write(*,*)'calcul des elts de passage de T vers M'

	call TMT(Atemp,Btemp,Ttemp,dndl,dnt,dn1,dn2,
     %m1temp,A,B,T,ndl,nt,n1,n2,m1,tabTM)

c!!!	write(*,*)'calcul des elts de passage de M vers T'

	call TMT(A,B,T,ndl,nt,n1,n2,m1,Atemp,Btemp,
     %Ttemp,dndl,dnt,dn1,dn2,m1temp,tabMT)

c!!!	write(*,*)'fin de calcul des elts de passage'

c!!!	open(4,file='tabTM.lis',status='UNKNOWN')
c!!!	write(4,*)(tabTM(i,1),i=1,ndl)
c!!!	write(4,*)(tabTM(i,2),tabTM(i,3),tabTM(i,4),i=1,ndl)
c!!!	close(4)

c!!!	open(4,file='tabMT.lis',status='UNKNOWN')
c!!!    write(4,*)(tabMT(i,1),i=1,dndl)
c!!!	write(4,*)(tabMT(i,2),tabMT(i,3),tabMT(i,4),i=1,dndl)
c!!!	close(4)
      
c     ========================================
c     ENTREE DES DIVERSES CONSTANTES PHYSIQUES
c      r:RAYON TERRESTRE en metres
c      b0:INDUCTION MAGNETIQUE DANS LE PLAN EQUATORIAL en tesla
c      mass: masse d'un electron en Kg
c      el  : constante d'homogeneite
c      q   : charge d'un electron en Coulomb
c      Nmax:LE NOMBRE MAXIMAL DE PARTICULES DANS DE DOMAINE
c      emax: valeur initiale optimale de l'energie en Kev
c
c     ON CALCULERA LES CONSTANTES:
c      em=sqrt(8*e/pi*mass)
c      ab0=a*a*a*a*a*b0*b0*b0
c     ========================================      
       pi=4.*atan(1.)
       r=6.37e+06
       b0=3.08e-05
       mass=9.e-31
       el=1.6e-16
       q=-1.6e-19
       Nmax=rnmax*r/b0
       Nmax=Nmax/((A(1))**(4.))
c       emax=0.4

	em=sqrt(8*el/(pi*mass))

cHAO97	On multiplie par un coefficient pour tenir compte d'une
cHAO97	non diffusion forte en angle d'attaque pour la precipitation
cHAO97	d'electrons (voir Schumaker et al. (1989), JGR, 10061).

cHAO97	Ici le coefficient est mis a 0.4

	em=em*1.

cHA097	

       ab0=r*r*r*r*r*b0*b0*b0

c      1.  PAS DE TEMPS              : deltat
c      2.  Nombre maximale de triangles parcouru pendant un pas
c          de temps:nemax
c      3.  Les conditions aux limites dependent du temps : iint=1
c                                      sinon                iint=0
c      4.  Choix de la methode :  metod
c                      =3   DECOMPOSITION SM EXACT
c                      =4   DECOMPOSITION SM APPROCHE MASS-LUMPING
c      5.   Valeur extremale du potentiel electostatique sur le bord
c                    phiopt
c
c

      nemax=10
      iint=0
      
c 	Le parametre iinti effectue les meme choses que iint mais
c 	pour les ions.

c 	Le parametre ibrupt permet d'avoir un bord raide ou doux 
c 	(sinusoidal) en azimut:
c 				quand ibrupt=0, pas de bord raide
c 				quand ibrupt=1, bord raide
	
      metod=4

c     ========================================
c     INITIALISATION EN TEMPS
c     ========================================      
c
c             DONNEES DE L ENERGIE ET DU NOMBRE DE PARTICULES AU
c             TEMPS t=0
c
c           ET DE LA PARTIE DU POTENTIEL INDEPENDANTE DU TEMPS: dans
c           uo
c
      heure=0.
      
c 	Initialisation des electrons

      do 1 i=1,ndl
      s=(B(i)-pi/2.)*(3*pi/2.-B(i))
      if((m2(i).eq.1).and.(s.le.0))goto 6
      ane1(i)=0.
      ener1(i)=0.
      goto 1
 6    ane1(i)=1.
      ener1(i)=emax/emax
ccc      s=(B(i)-pi/3.)*(5.*pi/3.-B(i))
ccc      if(s.le.0) goto 7
ccc      ener1(i)=0.5*(1+cos(6*B(i)))*emax
ccc      goto 1
ccc7     ener1(i)=emax
CCC	if (m2(i).eq.1) then
CCC		ener1(i)=emax/emax
CCC		ane1(i)=1.
CCC			else
CCC		ener1(i)=0.
CCC		ane1(i)=0.
CCC			endif
 1     continue

c 	Initialisation des ions 
	
	do 6661 i=1,ndl
	Ti(i)=0.
	Rni(i)=0.
	RNNN(i)=0.
	TV(i)=0.
6661	continue

	do 6651 i=1,dndl
	Titemp(i)=0.
	Rnitemp(i)=0.
	RNNNtemp(i)=0.
	TVtemp(i)=0.
6651	continue


	do 6662 i=1,(1*n1)

cHAO97

cHAO97	call limitIII(Ti(i),i,B(i),ibrupt,iinti,n1,Bfm,Bfs,heure,
cHAO97     %delta2,heuref,1.)
     
cHAO97    	call limitIII(Rni(i),i,B(i),ibrupt,iinti,n1,Bfm,Bfs,heure,
cHAO97     %delta2,heuref,1.)

cHAO97

	Rni(i)=1.
      	Ti(i)=rtmax/rtmax

	s=(B(i)-pi/2.)*(3*pi/2.-B(i))

	if (s.ge.0)	then

        Rni(i)=0.
      	Ti(i)=0.

			endif
      
6662	continue

	do 6652 i=1,(1*dn1)

cHAO97

cHAO97	call limitIII(Titemp(i),i,Btemp(i),ibrupt,iinti,dn1,Bfm,Bfs,
cHAO97     %heure,delta2,heuref,1.)
     
cHAO97  call limitIII(Rnitemp(i),i,Btemp(i),ibrupt,iinti,dn1,Bfm,Bfs,
cHAO97     %heure,delta2,heuref,1.)

cHAO97

	Rnitemp(i)=1.
      	Titemp(i)=rtmax/rtmax

	s=(B(i)-pi/2.)*(3*pi/2.-B(i))

	if (s.ge.0)	then

        Rnitemp(i)=0.
      	Titemp(i)=0.

			endif
     	     	
6652	continue


      do 10 i=1,ndl
      uo(i)=emax*ener1(i)*1.5/Amin**(8./3.)
      uo(i)=uo(i)*((A(i))**(8./3.))
10    continue

	if (istati.eq.0)	goto 2993

c!!!	write(*,*)'Lecture de densii....'

cHAO97

	open(96,file='densii.lis',status='UNKNOWN')
	open(97,file='tempii.lis',status='UNKNOWN')
	open(98,file='tempe1i.lis',status='UNKNOWN')
	open(99,file='nbre1i.lis',status='UNKNOWN')


	do 6789 j=1,62

	read(96,*)(Rni(i),i=1,ndl)

	read(97,*)(Ti(i),i=1,ndl)
	
	read(98,*)(ener1(i),i=1,ndl)
		
	read(99,*)(ane1(i),i=1,ndl)

 6789	continue
	
	close(96)
	close(97)
	close(98)
	close(99)

cHAO97

	open(4,file='uoi.lis',status='UNKNOWN')
	read(4,*)(uo(i),i=1,ndl)
	close(4)

c!	goto 2993
	
	do 2992 i=1,ndl

	ane1(i)=ane1(i)/Nmax
	ener1(i)=ener1(i)*(Amin**(8./3.))/(A(i)**(8./3.))
	ener1(i)=ener1(i)/emax
	Rni(i)=Rni(i)/(rnmax*(A(i)/A(1))**(4.))
	Ti(i)=Ti(i)/(rtmax*(A(i)/A(1))**(8./3.))

 2992	continue
 2993	continue

	write(*,*)(Rni(i),i=1,80)
	write(*,*)(ane1(i),i=1,80)


C *************************************************************** FIN

								endif

C *************************************************  D'INITIALISATION

	write(*,*)(Rni(i),i=1,80)
	write(*,*)(ane1(i),i=1,80)


c	*** INTERPOLATION DU POTENTIEL DE CI SUR MA GRILLE  ***
c	*** ON CONVERTIT DES **temp aux VARIABLES SANS temp ***

c!!!	write(*,*)'J''interpole le potentiel'

	icomp=9
	do 9000 i=1,25
	icomp=icomp+1
		do 9010 j=1,80
		utemp((i-1)*80+j)=utem2((icomp-1)*80+j)
 9010		continue
 9000	continue
			
	call interpolation(Atemp,Btemp,Ttemp,dndl,dnt,dn1,
     %dn2,m1temp,A,B,T,ndl,nt,n1,n2,m1,
     %utemp,utemp,utemp,utemp,utemp,u,u,u,u,u,tabTM)

c!!!	write(*,*)'interpolation du potentiel finie'

c	*** ON AVANCE ***
     
       nnp=0
c     if(nnp.eq.0) write (02) (u(i),i=1,ndl)
       do 2 i=1,ndl
       cor=r*r*b0*2.*pi*(1.e-3)/86400
       u3(i)=u(i)
       u(i)=u(i)-cor*A(i)
       u4(i)=u(i)
       ui(i)=u(i)
       u(i)=-(u(i)*1.e+03
     %+uo(i)*el/q*((A(i)/Amin)**(8./3.)))/((r**2)*b0)
       if (Rni1(i).lt.(1.e-09)) goto 2
       RNNN1(i)=Rni(i)/Rni(1)
       TV1(i)=Ti(i)/Ti(1)
 2     continue

 1965	delta2=deltat

	nemax=10
	iint=0
	metod=4
	
csss	write(*,*)'deltat ',deltat

csss	if (Rni1(i).lt.(1.e-09)) goto 9999
csss	open(4,file='statpot.lis',access='append')
csss	write(*,*)('*',i=1,70)
csss	write(*,*)'heure ',heure
csss	call L2(u0,u3,m1,A,B,T,ndl,nt,n1,n2,res,diff2,relat)
csss	write(*,*)'diff2 ',diff2,'res ',res,'relat ',relat
csss	close(4)
csss	open(4,file='stattem.lis',access='append')
csss	write(*,*)('*',i=1,70)
csss	write(*,*)'heure ',heure
csss	call L2(T0,TV1,m1,A,B,T,ndl,nt,n1,n2,res,diff2,relat)
csss	write(*,*)'diff2 ',diff2,'res ',res,'relat ',relat
csss	close(4)
csss	open(4,file='statnbr.lis',access='append')
csss	write(*,*)('*',i=1,70)
csss	write(*,*)'heure ',heure
csss	call L2(rN0,RNNN1,m1,A,B,T,ndl,nt,n1,n2,res,diff2,relat)
csss	write(*,*)'diff2 ',diff2,'res ',res,'relat ',relat
csss	close(4)
 9999	continue
	do 8765 i=1,ndl
	u0(i)=u3(i)
	T0(i)=TV1(i)
	rN0(i)=RNNN1(i)
 8765	continue

c     ========================================
c     DEBUT DU PROGRAMME EVOLUTIF CONSTITUE DE DEUX ETAPES:
c          -une etape de transport par caracteristiques (caract)
c          -une derniere etape prenant compte de la non linearite
c           (ivol1)
c     ========================================
c     ========================================
c     pour t=temp=0 ,on a entre les valeurs initiales de
c         ane dans ane1:NOMBRE D ELECTRONS 
c         ener dans ener1:ENERGIE DES PARTICULES
c     ========================================
c     ========================================
c     PREMIERE ETAPE EVOLUTIVE ENTRE T ET T+deltat
c     ========================================      
      do 8 i=1,ndl
      ener2(i)=ener1(i)
8     ane2(i)=ane1(i)
c     ========================================
c     RECHERCHE DES TRAJECTOIRES :CARACTERISTIQUES
c         tmin:le temps initial
c         tmax:le temps final
c         itria:le nombre de triangles parcourus
c         isort:0 si la carac passant par i est situee dans le domaine
c                :1 si la carac atteint la frontiere
c         vitess:la vitesse donnee comme le rot de vitess constante
c         par triangle
c         eps:precision
c         xn et yn:les coordonnees des points de contact
c         temp:le temps mis pour effectue chaque etape
c         ntria:le numero des triangles traverses
c         nemax:nombre maximal de triangles traverses
c     ========================================

      heure=0.     
      tmin=heure
      tmax=heure+deltat
      eps=1.e-09
       
c     print 106
c     goto 300
	
c!!!	write(*,*)'Je suis en CARACT'

c *********************************************************
c *** on indique si on bouge les electrons ou les ions: ***
c ***		- ilect=1:	on bouge les electrons  ***
c ***		- ilect=0:	on bouge les ions	***
c *********************************************************

	ilect=1

      call caract(ilect,tmin,tmax,itria,isort,u,eps,xn,
     %yn,temp,ntria,nemax
     %,1,ndl,nt,B,A,m1,m2,xu,xv,ane3,ener3,ane1,ener1,
     %T,metod,iint,uo,b1,result,contro)

c      write (50, 120)
c120   format(1x,'COMPARAISON ENTRE DEUX ITERATIONS:ENERGIE
c      CINETIQUE')
c     call compare(b1,uo,ndl)
      do 56 i=1,ndl
      uo(i)=b1(i)
56     b1(i)=0.

c	***	On moyenne uo sur les triangles 	***

	do 8100 i=1,ndl
	hh(i)=0.
	aire1(i)=0.
 8100 	continue

 	do 8020 i=1,nt 
	x1=B(m1(1,i))
	x2=B(m1(2,i))
	x3=B(m1(3,i))
	y1=A(m1(1,i))
	y2=A(m1(2,i))
	y3=A(m1(3,i))
	nqi=i-(i/(2*n1))*2*n1
	if (nqi.eq.0) x1=2*pi
	if (nqi.eq.0) x2=2*pi
	if (nqi.eq.(2*n1-1)) x2=2*pi
	delta=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)

	do 8030 j=1,3
	noeud=m1(j,i)
	aire1(noeud)=aire1(noeud)+delta
 8030	continue

 8020	continue

	call intmas(uo,hh,ndl,nt,m1,A,B)
 
	do 8040 noeud=n1+1,ndl
	uo(noeud)=6*hh(noeud)/aire1(noeud)
 8040	continue


	ifact=1
csss	write(*,*)'ifact vaut: ',ifact

	if ((ie-(ie/ifact)*ifact).eq.0)	then

c!!!	write(*,*)'Je suis en IVOL1'

      call ivol1(ane3,ener3,ane1,ener1,u,metod,ndl,nt,m1,m2,A,B
     %  ,b1,b2,b3,b4,b5,b6,ifact)

					else

	do 1600 i=n1+1,ndl
	ane1(i)=ane3(i)
	ener1(i)=ener3(i)
	if (abs(ane1(i)).lt.((1.e-05)*(ane1(1)))) then 
						ane1(i)=0.
						ener1(i)=0.
				 	     	endif
 1600	continue

					endif

c!!!  	write(*,*)'COUCOU JE SUIS EN CARACTI'   


ccccc	open(4,file='u.lis',status='UNKNOWN')
ccccc	write(4,*)(u3(i),i=1,ndl)
ccccc	close(4)
	
ccccc	open(4,file='Rni.lis',status='UNKNOWN')
ccccc	write(4,*)(Rni(i),i=1,ndl)
ccccc	close(4)
	
ccccc	open(4,file='Ti.lis',status='UNKNOWN')
ccccc	write(4,*)(Ti(i),i=1,ndl)
ccccc	close(4)
	
	ilect=0
	
	 call caract(ilect,tmin,tmax,itria,isort,ui,eps,xn,
     % yn,temp,ntria,nemax
     %,1,ndl,nt,B,A,m1,m2,dxdt,dydt,Rnic,Tic,Rni,Ti,
     %T,metod,iint,uo,b1,result,contro)
     
cccc     open(4,file='Ntrans.lis',status='UNKNOWN')
cccc	write(4,*)(Rnic(i),i=1,ndl)
cccc	close(4)
	
cccc	open(4,file='Ttrans.lis',status='UNKNOWN')
cccc	write(4,*)(Tic(i),i=1,ndl)
cccc	close(4)

		
cccc	open(4,file='Rni1.lis',status='UNKNOWN')
cccc	write(4,*)(Rni(i),i=1,ndl)
cccc	close(4)
	
cccc	open(4,file='Ti1.lis',status='UNKNOWN')
cccc	write(4,*)(Ti(i),i=1,ndl)
cccc	close(4)
	
c!!!	write(*,*)'Je vais en Ivol5'
     			
        call ivol5(Tic,Rnic,nt,Ti,Rni,ndl,hh,hh1,hh2,
     %aire1,RNNN,TV,Ti1,Rni1,A,B,T,m1,m2,iinti,ibrupt,heure,
     %Bfm,Bfs,delta2,heuref,rtmax,rnmax,b5,b6,znmax,
     %ztmax,inum,heure,ifact,ie)

cccc	open(4,file='Rnic1.lis',status='UNKNOWN')
cccc	write(4,*)(Rnic(i),i=1,ndl)
cccc	close(4)
	
cccc	open(4,file='Tic1.lis',status='UNKNOWN')
cccc	write(4,*)(Tic(i),i=1,ndl)
cccc	close(4)


cccc	open(4,file='Rni2.lis',status='UNKNOWN')
cccc	write(4,*)(Rni(i),i=1,ndl)
cccc	close(4)

cccc	open(4,file='Ti2.lis',status='UNKNOWN')
cccc	write(4,*)(Ti(i),i=1,ndl)
cccc	close(4)	

       
c    	 open(4,file='etat.lis',status='UNKNOWN')
c    	 write(4,*)('noeud :',i,Ti(i),Rni(i),PHIi(i),i=1,ndl)
c     	close(4)

c     ========================================
c      CALCUL DU FLUX F ET DE L'ENERGIE DES ELECTRONS PRECIPITES
c
c        FE : FLUX D'ENERGIE PRECIPITES
c     ========================================
    
      do 4 i=1,ndl
      teme=emax*(ener1(i))*(A(i)**(8./3.))/(Amin**(8./3.))
      E(i)=teme
      F(i)=(((r*b0*A(i))**4)*em/(4*ab0))*Nmax*ane1(i)*
     %   sqrt(teme)
      FE(i)=2*F(i)*teme*1.e+03*el
 4    continue

c***	CALCUL DE LA PRESSION IONIQUE ***

	do 801 i=1,ndl
	if (Rni(i).gt.((1.e-06)*Rni(1))) then
		ppi(i)=Ti(i)*Rni(i)*rtmax*rnmax*(A(i)/A(1))**(20./3.)
		PHIi(i)=ppi(i)
		ppi(i)=1.6e-16*ppi(i)
				else
		ppi(i)=0.
				endif
 801	continue

c***	CALCUL DES GRADIENTS DE PRESSION IONIQUE ***

	call GRADP(aire1,B,A,PHIi,nt,ndl,m1,m2,dppib,dppia,betaai)
	
c***	CALCUL DES COURANTS ALIGNES ***

	do 2994 i=1,ndl

	ali(i)=(betaai(i)*1.6e-16)*2.*
     %(sqrt(4.-3.*A(i)))/(b0*((A(i))**(5.)))
					
 2994	continue


c	*** INTERPOLATION DE Rni, Ti, FE et E DE MA GRILLE SUR  ***
c 	*** CI. ON CONVERTIT DES VARIABLES SANS temp aux **temp ***
c	*** ON INTERPOLERA LE PRODUIT DE Rni et de Ti POUR LE   ***
c	*** CALCUL DU POTENTIEL ET NON Rni et Ti SEPAREMENT 	***
c 	*** POUR EVITER LES PB DE PRECISION			***

c 	*** Rni contiendra le produit et Ti sera egal a 1 	***
c	*** partout pour conformite pour la subroutine ellconi  ***

	do 6640 i=1,ndl
	Rni11(i)=Rni(i)*Ti(i)
	Ti11(i)=1.
 6640	continue

cHAO97

 	goto 7777

cHAO97

c!	On fait un fit a l'ordre 5 en Fourrier de la distribution des
c!	courants alignes

	par11=par1

c!	write(*,*)par11

	call Fourrier(ali,alifit,alii,rafou,iafou,
     & B,par11,n1,n2)

c!	write(*,*)(ali(i),i=1,ndl)

c!!!	write(*,*)'JE SUIS EN INTERPOLATION'

cHAO97

 7777	continue

cHAO97

cHAO97	call interpolation(A,B,T,ndl,nt,n1,
cHA097     %n2,m1,Atemp,Btemp,Ttemp,dndl,dnt,dn1,dn2,m1temp,
cHA097     %Rni11,Ti11,FE,E,ali,Rnitemp,Titemp,FEtemp,Etemp,alitem,
cHAO97     %tabMT)

	call interpolation(A,B,T,ndl,nt,n1,
     %n2,m1,Atemp,Btemp,Ttemp,dndl,dnt,dn1,dn2,m1temp,
     %Rni,Ti,FE,E,ali,Rnitemp,Titemp,FEtemp,Etemp,alitem,
     %tabMT)

c	On calcule les flux ioniques

	
	do 44 i=1,49
		do 444 j=1,80
	
	kk=(i-1)*80+j
	FEion(kk)=0.
	Fion(kk)=0.

	kk1=(i-10)*80+j

	if	((i.ge.10).and.(i.le.34))	then

      temei=rtmax*(Titemp(kk1))*(Atemp(kk1)**(8./3.))/
     %(Atemp(1)**(8./3.))
      Fion(kk)=(((r*b0*Atemp(kk1))**4)*em/(4*ab0))*
     %Nmax*Rnitemp(kk1)*
     %   sqrt(temei)
     %/sqrt(1831.)
     %*(1./(0.4/0.4))
     %*(0.2/0.2)

c    Le facteur (1./sqrt(1831)) prend en compte le rapport de masse
c    des ions (ici protons) et des electrons.
c    Le facteur (1/0.4) est mis pour compenser le facteur des electrons
c    mis a 0.4.
c    Le facteur 0.2 est mis pour reduire la diffusion forte en angle
c    d'ataque des ions.

      FEion(kk)=2*Fion(kk)*temei*1.e+03*el

	Fion(kk)=temei

				endif

 444		continue
 44    continue

cHAO97


	do 7000 i=1,49
		do 7010 j=1,80
		Etem2((i-1)*80+j)=0.
		FEtem2((i-1)*80+j)=0.
		alite2((i-1)*80+j)=0.
 7010	continue
 7000	continue

	icomp=9
	do 8000 i=1,25
	icomp=icomp+1
		do 8010 j=1,80
		Etem2((icomp-1)*80+j)=Etemp((i-1)*80+j)
		FEtem2((icomp-1)*80+j)=FEtemp((i-1)*80+j)

cHAO97 On supprime les J// sur la frontiere car ils n'ont
cHAO97 aucun sens.

		if (icomp.eq.10)	then

		alite2((icomp-1)*80+j)=0.

					else	
			
		alite2((icomp-1)*80+j)=alitem((i-1)*80+j)

					endif

cHAO97

 8010		continue
 8000	continue

c!!!	write(*,*)'Ecriture des grandeurs non adiabatiquees'

c!!!	open(4,file='Enon.lis',status='unknown')
c!!!	write(4,*)(Etem2(i),i=1,dpar11)
c!!!	close(4)

c!!!	open(4,file='FEnon.lis',status='unknown')
c!!!	write(4,*)(FEtem2(i),i=1,dpar11)
c!!!	close(4)

	do 6000	j=1,80

		imax=10

		do 6010 i=imax,imax

cc		if (imax.gt.0)	goto 6010

cc		if (FEtem2((i-1)*80+j).gt.(1.e-06))	then
		imax=i
		dAimax=(90.-dTtem2(imax))*pi/180.
		dAimax=(dAimax**(2.))
		FEimax=FEtem2((imax-1)*80+j)
cc							endif
		
 6010		continue

		do 6020	i=2,(imax-1)

		dAima=(90.-dTtem2(i))*pi/180.
		dAima=(dAima**(2.))
		FEtem2((i-1)*80+j)=FEimax*((dAima/dAimax)**(8.))		
		
 6020		continue

		FEtem2(1)=FEtem2(80+j)+FEtem2(1)

		imax=10

		do 6030 i=imax,imax

cc		if (imax.gt.0)	goto 6030

cc		if (Etem2((i-1)*80+j).gt.(1.e-06))	then
		imax=i
		dAimax=(90.-dTtem2(imax))*pi/180.
		dAimax=(dAimax**(2.))
		Eimax=Etem2((imax-1)*80+j)
cc							endif
		
 6030		continue

		do 6040	i=2,(imax-1)

		dAima=(90.-dTtem2(i))*pi/180.
		dAima=(dAima**(2.))
		Etem2((i-1)*80+j)=Eimax*((dAima/dAimax)**(8./3.))		
		
 6040		continue

		Etem2(1)=Etem2(80+j)+Etem2(1)

 6000	continue

	do 6050	j=2,80
	FEtem2(j)=FEtem2(1)/80.
	Etem2(j)=Etem2(1)/80.
 6050	continue

	FEtem2(1)=FEtem2(2)
	Etem2(1)=Etem2(2)		

c!!!	write(*,*)'INTERPOLATION FINIE'

c***	ON FICHE LES GRANDEURS INTERPOLEES

	femod=0.

	do 1100 i=1,49

	do 1110 j=1,80

	if (FEtem2((i-1)*80+j).gt.femod)	then

	femod=FEtem2((i-1)*80+j)

						endif

 1110	continue

 1100	continue

	do 1120 i=1,49

	do 1130 j=1,80

	FEtem3((i-1)*80+j)=FEtem2((i-1)*80+j)
	Etem3((i-1)*80+j)=Etem2((i-1)*80+j)

	if (Etem3((i-1)*80+j).lt.Etem2(9*80+1))	then

	Etem3((i-1)*80+j)=Etem2(9*80+1)

						endif

	if (FEtem2((i-1)*80+j).lt.(femod/30.))	then

	FEtem3((i-1)*80+j)=1.e-16
	Etem3((i-1)*80+j)=Etem2(9*80+1)

						endif

1130	continue

1120	continue

	ie=ie+1

	return
	end

c==================================================================
c==================================================================

       subroutine disc1(n,l,mm1,mm2,nnt,iopt)
       dimension mm1(3,nnt),mm2(1)
c*****m1:contient pour chaque triangle la liste des noeuds
c*****numerotes dans le sens positif : les dl sont les sommets
c*****m2:egale 0 a l'interieur
c*****         1 sur le bord
       l1=l*2
       l2=l+1
       do 1 i=1,l1
          mod1=i/2
          mod=i-2*mod1
          if(mod.eq.0)goto2
          mm1(1,i)=mod1+1
          mm1(2,i)=mm1(1,i)+1
          mm1(3,i)=mm1(1,i)+l
          goto1
 2        mm1(1,i)=mod1+1
          mm1(2,i)=mm1(1,i)+l
          mm1(3,i)=mm1(2,i)-1
 1     continue
       i=l1-1
       mm1(2,i)=1
       i=l1
       mm1(1,i)=1
       mm1(2,i)=l+1
       l3=l1+1
       do 3 i=l3,nnt
          nq=i/l1
          nr=i-nq*l1
          if (nr.eq.0) goto4
          do 5 j=1,3
             mm1(j,i)=mm1(j,nr)+l*nq
 5        continue
          goto 3
 4        do 6 j=1,3
             mm1(j,i)=mm1(j,l1)+l*(nq-1)
 6        continue
 3     continue
       n3=(n+1)*l
       do 7 i=1,n3
          mm2(i)=0
 7     continue
       do 8 i=1,l
          mm2(i)=mm2(i)+1
 8     continue
       n4=n3-l+1
       do 11 i=n4,n3
          mm2(i)=mm2(i)+2
 11    continue
       if (iopt.lt.1) return
       print 102
       print 100,((mm1(j,i),j=1,3),i=1,nnt)
       print 103
       print 101,(mm2(i),i=1,n3)
 102   format (1x,'TRIANGULATION DU DISQUE','NUMEROTATION')
 100   format(4(3(1x,i4),1x,':'))
 103   format(1x,'CONDITIONS AUX LIMITES')
 101   format(1x,20(i3))
       return
       end

c ===================================================================
c ===================================================================

	subroutine limitIII(Ti,i,B,ibrupt,iinti,nn1,Bfm,Bfs,heure,
     %delta2,heuref,rnorm)
     
     	real*4    Ti,B
	integer*4 ibrupt,iinti,nn1,i
	real*4    Bfm,Bfs,heure,delta2,heuref,rnorm
	
	Ti=rnorm

	if (iinti.eq.1) 			then
			if (heure.le.heuref)	then
						ss=heure/heuref
						Ti=Ti*ss
						endif
			if (heure.lt.(delta2/2.)) then
						Ti=0.
						endif
						endif

	return
	end

c ===============================================================
c ===============================================================

     	subroutine Coeff(x1,x2,x3,y1,y2,y3,N1,N2,N3,N10,N01,N00)

	real*4 x1,x2,x3,y1,y2,y3,N1,N2,N3,N10,N01,N00
	aire=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
	N10=N1*(y2-y3)+N2*(y3-y1)+N3*(y1-y2)
c	write(*,*)'aire :',aire
c	write(*,*)'x1... :',x1,x2,x3,'y1...',y1,y2,y3
	N10=N10/aire
	N01=N1*(x3-x2)+N2*(x1-x3)+N3*(x2-x1)
	N01=N01/aire
	N00=N1-N10*x1-N01*y1
	return
	end
	

c ===================================================================
c ===================================================================

	subroutine L2(psi1,psi2,m1,cooy,coox,T,ndl,nt,n1,n2,res,
     %diff2,relat)
	dimension psi1(1),psi2(1),m1(3,nt),cooy(1),coox(1),T(1)
	real*4 diff2
	diff2=0.
	res=0.
	pi=4.*atan(1.0)
      do 8 i=1,nt
      x1=coox(m1(1,i))
      x2=coox(m1(2,i))
      x3=coox(m1(3,i))
      y1=cooy(m1(1,i))
      y2=cooy(m1(2,i))
      y3=cooy(m1(3,i))
      cost1=2.*cos(T(m1(1,i)))
      cost2=2.*cos(T(m1(2,i)))
      cost3=2.*cos(T(m1(3,i)))
      nqi=i-(i/(2*n1))*2*n1
      if(nqi.eq.0) x1=2*pi
      if(nqi.eq.0) x2=2*pi
      if(nqi.eq.2*n1-1) x2=2*pi
      ps1=psi1(m1(1,i))-psi2(m1(1,i))
      ps1=(ps1**2.)/cost1
      ps2=psi1(m1(2,i))-psi2(m1(2,i))
      ps2=(ps2**2.)/cost2
      ps3=psi1(m1(3,i))-psi2(m1(3,i))
      ps3=(ps3**2.)/cost3
      aire=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
      diff2=diff2+(ps1+ps2+ps3)*aire/6.
      ps1=psi2(m1(1,i))
      ps1=(ps1**2.)/cost1
      ps2=psi2(m1(2,i))
      ps2=(ps2**2.)/cost2
      ps3=psi2(m1(3,i))
      ps3=(ps3**2.)/cost3
      aire=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
      res=res+(ps1+ps2+ps3)*aire/6. 
8      continue
	relat=diff2/res
	return
	end

c==================================================================
c==================================================================

      subroutine caract(ilect,tmin,tmax,itria,isort,psi,eps,
     %xn,yn,temp,
     %     ntria,nemax,iopt,ndl,nt,coox,cooy,m1,m2,xu,xv,ane3,ener3,
     %     ane1,ener1
     %     ,T,metod,iint,uo,u1,result,contro)  
c
c             CALCUL DES CARACTERISTIQUES ENTRE tmin ET tmax;La
c             vitesse
c                est supposee constante en temps dans ce intervalle
c
	real*4 Nmax
      common/donnees/e,q,em,ab0,deltat,r,b0,emax,ymax,ymin,h2,h1,n1,n2
     %,Nmax
      common/triang1/n21,n22,n23,h21,h22,h23,yint1,yint2,raison
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      dimension iv(6),temp(10),xn(10),yn(10),psi(1),result(6,168)
      dimension ntria(10),ane3(1),ener3(1),
     %    xu(nt),xv(nt),coox(1),cooy(1),m1(3,nt),m2(1),T(1)
     %    ,ener1(1),ane1(1),uo(1),u1(1)

      integer*4 compt,compt1,contro(6),H
      real*4 NN1,NN2,NN3,N10,N01,N00
c
c            CALCUL DE LA VITESSE PAR TRIANGLE
c
      compt=0
      compt1=0
      pi=4*atan(1.)

ccc	if (ilect.eq.0)	then
ccc
ccc	do 88 i=1,ndl
ccc	psi(i)=-psi(i)*1000./(r*r*b0)
ccc 88	continue
ccc
ccc			endif

      do 8 i=1,nt

      x1=coox(m1(1,i))
      x2=coox(m1(2,i))
      x3=coox(m1(3,i))
      y1=cooy(m1(1,i))
      y2=cooy(m1(2,i))
      y3=cooy(m1(3,i))
      nqi=i-(i/(2*n1))*2*n1
      if(nqi.eq.0) x1=2*pi
      if(nqi.eq.0) x2=2*pi
      if(nqi.eq.2*n1-1) x2=2*pi
      ps1=psi(m1(1,i))
      ps2=psi(m1(2,i))
      ps3=psi(m1(3,i))
      aire=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
      xv(i)=-((ps1-ps2)*(y1-y3)-(ps1-ps3)*(y1-y2))/(aire)
      xu(i)=-((ps1-ps2)*(x1-x3)-(ps1-ps3)*(x1-x2))/(aire)
      
	if (ilect.eq.0) then
	
c *** 	REMISE A 1 de L'INDICATEUR DE DENSITE NULLE OU PAS ***

	flag1=1
	
c	if (i.lt.161) then
c	write(*,*)'i :',i
c			endif
	phi1=psi(m1(1,i))*1000.
	phi2=psi(m1(2,i))*1000.
	phi3=psi(m1(3,i))*1000.
	T1=ener1(m1(1,i))
	T2=ener1(m1(2,i))
	T3=ener1(m1(3,i))
	NN1=ane1(m1(1,i))
	NN2=ane1(m1(2,i))
	NN3=ane1(m1(3,i))
	x1=coox(m1(1,i))
	x2=coox(m1(2,i))
	x3=coox(m1(3,i))
	y1=cooy(m1(1,i))
	y2=cooy(m1(2,i))
	y3=cooy(m1(3,i))
	nqi=i-(i/(2*n1))*2*n1
	if (nqi.eq.0) x1=2*pi
	if (nqi.eq.0) x2=2*pi
	if (nqi.eq.(2*n1-1) ) x2=2*pi
c	if (i.lt.161) then
c	write(*,*)'x1...y1...',x1,x2,x3,y1,y2,y3
c			endif
	call Coeff(x1,x2,x3,y1,y2,y3,NN1,NN2,NN3,N10,N01,N00)
	call Coeff(x1,x2,x3,y1,y2,y3,T1,T2,T3,T10,T01,T00)
	call Coeff(x1,x2,x3,y1,y2,y3,phi1,phi2,phi3,phi10,phi01,phi00)
	xc=(x1+x2+x3)/3.
	yc=(y1+y2+y3)/3.
	flag1=0
	z=0.
	xv(i)=phi10+flag1*T10*1000.+z*N10
	a=6.37e+06
	B0=3.08e-05
	xv(i)=xv(i)/(a*a*B0)
	xu(i)=phi01+flag1*T01*1000.+z*N01
	xu(i)=-xu(i)/(a*a*B0)
c	open(4,file='derivee.lis',access='append')
c	write(4,*)'T10...',T10,T01,T00,N10,N01,N00
c	close(4)

			endif
 
 8      continue


c     ========================================
c     BOUCLE SUR LES NOEUDS
c     ITRIANG:NOMBRE DE TRIANGLES TRAVERSES
c     ISORTIE:0 SI LE PARCOURS DE I EST SITUE DANS LE DOMAINE
c            :1 SI I SORT DU DOMAINE
c     IS:LE NOEUD COURANT ETUDIE      
c     ========================================      
	
	nmax1=0
	itria=0
	depas=0.
	
      do 1 i=1,ndl
      if (itria.gt.nmax1) nmax1=itria
      if (itria.gt.1) depas=depas+1.
      itria=0
c     iji=0
      isort=0
      do 9 j=1,nemax
      xn(j)=0.
      yn(j)=0.
      temp(j)=0.
      ntria(j)=0.
9     continue
      is=i
      ttol=tmax
c     ========================================
c     DETERMINATION DES VOISINS DE I
c       iv(5)---iv(4)
c         | .     |
c         |   .   |   .
c         |     . |     .
c       iv(6)-----I-----iv(3)
c           .     | .     |
c             .   |   .   |
c               . |     . |
c               iv(1)---iv(2)
c     ========================================
30    nsc=6
            iv(1)=is-n1
      iv(2)=is-n1+1
      iv(6)=is-1
      iv(3)=is+1
      iv(5)=is+n1-1
      iv(4)=is+n1
      if((is-(is/n1)*n1).eq.0) goto 2
      if((is-(is/n1)*n1).eq.1)goto 3
      goto 4
2     continue
      iv(2)=is-2*n1+1
      iv(3)=is-n1+1
      goto4
3     continue
      iv(6)=is+n1-1
      iv(5)=is+2*n1-1
4     continue
      nms=1
      if(is.le.n1) nms=3
      if(is.le.n1) nsc=5
c     ========================================
c     ETUDE DU TRIANGLE TRAVERSE :RECHERCHE DE CE TRIANGLE
c     ========================================      
      x=coox(is)
      y=cooy(is)
      ps=psi(i)
	if (ilect.eq.0)	ps=-psi(i)*1000./(r*r*b0)
      x1=x
      y1=y
      ps1=ps
      do 5 j=nms,nsc
      H=j
      if((is.gt.(ndl-n1)).and.(j.ge.3)) H=nsc
      j1=H+1
      if(H.eq.6) j1=1
      ps2=psi(iv(H))
      ps3=psi(iv(j1))
	if (ilect.eq.0) then
			ps2=-psi(iv(H))*1000./(r*r*b0)
      			ps3=-psi(iv(j1))*1000./(r*r*b0)
			endif
c106   format(1x,'ps1,ps2,ps3',3e12.5)
       s=(ps2-ps)*(ps-ps3)
      if ((abs(ps-ps2).lt.eps).or.(abs(ps-ps3).lt.eps)) goto 789
      if(s.lt.(-eps)) goto 5
 789	continue      
      i2=iv(H)
      i3=iv(j1)
      x2=coox(iv(H))
      y2=cooy(iv(H))
      x3=coox(iv(j1))      
      y3=cooy(iv(j1))
c     ========================================
c     DETERMINATION DES DEUX VOISINS SUR LE TRIANGLE
c     ========================================      
      cb1=ps2-ps1
      cb2=ps1-ps3
      if ((abs(cb1).lt.eps).or.(abs(cb2).lt.eps)) goto 788
      if((cb1.lt.eps).or.(cb2.lt.eps)) goto5
 788	continue
c     ========================================
c     si ps2<ps<ps3 sens des temps positif
c     si ps3<ps<ps2 sens des temps negatif
c     ============================================
c      determinons le triangle traverse
c     ============================================
      if((is-(is/n1)*n1).eq.1) goto21
      if((is-(is/n1)*n1).eq.0) goto22
      goto23
21     if(H.lt.4) goto23
      x1=2*pi
      x=2*pi
       if(H.eq.4) x2=2*pi
       if(H.eq.6) x3=2*pi
      goto 23
22        if(H.ge.4) goto23
      if(H.le.2) x3=2*pi
      if(H.ge.2) x2=2*pi
      goto23
23    xc=1/3.*(x1+x2+x3)
      yc=1/3.*(y1+y2+y3)
c     nota=0
      call dtcar(xc,yc,x,y,ymax,ymin,n1,n2,h1,h2,not,not2)
      
c     if (i.eq.639) write(*,*) i,not,not1,not2,nota
c     if (i.eq.722) write(*,*) i,not,not1,not2,nota
c     if (i.eq.762) write(*,*) i,not,not1,not2,nota
c     ========================================
c     le point de depart est (x,y)
c     le triangle est not
c     dont les sommets sont i1,i2,i3
c     on cherche maintenant le point de contact xc,yc
c     puis le temp mis pour le parcourir
c     ========================================
15    continue
       if(abs(ps2-ps3).gt.eps) goto 92
        xc=x
        yc=y
c	write(*,*) eps,ps2,ps3
      goto 7
92	continue
c 	if (iji.eq.(1)) then
c	write(*,*) itria
c	endif
      xc=x2+(ps-ps2)/(ps3-ps2)*(x3-x2)
      yc=y2+(ps-ps2)/(ps3-ps2)*(y3-y2)
c      if (itria.eq.0) then
c      w2=x2
c      z2=y2
c                      else
c      w2=xn(itria)
c      z2=yn(itria)
c                      endif
c      if((abs(xc-w2).lt.(1.57e-03)).and.(abs(yc-z2).lt.(8.94e-03)))
c     % goto 7
      if ((yc.le.(ymax+eps)).and.(yc.ge.(ymin-eps))) goto 91
      isort=1
      goto 7
91    d=(x-xc)**2+(y-yc)**2
      u=(xu(not)**2)+(xv(not)**2)
      dt=sqrt(d/u)
      ttol=ttol-dt
      if(ttol.le.tmin) goto20
      itria=itria+1
c!!!  if(itria.gt.nemax) then
c      iji=1
c      write(*,*) ('COUCOU C EST MOI'),itria,j,not,i,xc,yc
c      write(*,*) ntria(3),ntria(4),ntria(1),ntria(2)
c!!!  print 114,nemax
c!!!  endif
      temp(itria)=dt
      xn(itria)=xc
      yn(itria)=yc
      ntria(itria)=not
c!!!      if (itria.ge.nemax) then
c!!!      open (4,file='probleme.lis',status='UNKNOWN')
c!!!      write(4,*) ('noeud de depart:'),i
c!!!      write(4,*) ('heure :'),tmin,tmax
c!!!      	do 1987 k=1,nemax
c!!!      write(4,*) ('triangle traverse:'),ntria(k)
c!!!      write(4,*)('Coordonnees du point d,arrivee sur le
c!!!     % triangle'),xn(k),yn(k)
c!!!1987 	continue      
c!!!      close(4)
c     open(4,file='debile.lis',status='UNKNOWN')
c     write(4,*)
c    % H,iv(H),j1,iv(j1),coox(iv(H)),cooy(iv(H)),coox(iv(j1))
c    % ,cooy(iv(j1)) 
c     close(4)
c!!!	STOP
c!!!      endif
      x=xc
      y=yc   
c     TEST
c         SAVOIR SI ON ARRIVE SUR UN NOEUD 
c          SI PS=PS2 =>is=i2
c          SI PS=PS3 =>is=i3
c        
c     ========================================	      
       if(abs(ps-ps2).lt.eps) goto 31
      if(abs(ps-ps3).lt.eps) goto 32
      goto 33
31	continue	    
c	if (iji.eq.(1)) then
c      write(*,*)('NOEUD 2')
c      endif
        is=i2
      goto 30
32	continue
c	if (iji.eq.(1)) then
c	write(*,*)('NOEUD 3')
c	endif
        is=i3
       goto 30
33      continue
c     ===============================================
c      ON VEUT TESTER SI ON EST ARRIVE SUR LA FRONTIERE
c     AVEC UNE VITESSE NORMALE NON NULLE
c     =================================================
ccc      if ((y.gt.ymin).and.(y.lt.ymax)) goto 34

      if ((abs(y-ymin).gt.eps).and.(abs(y-ymax).gt.eps)) goto 34

      if (abs(xv(not)).lt.eps) goto 34
      nu=-1
      if(y.ge.ymax) nu=1
      if (xv(not)*nu.lt.0) then
      isort=1
      endif
      goto 7
c     ========================================
c     recherche du nouveau triangle traverse
c     ========================================
c34    not2=-1
 34    continue

	if (abs(x2-x3).lt.(5.e-04)) then

****** Je suis sur une verticale ******

	ymoy=(y2+y3)/2.
	xmoy=x2
				    endif

	if (abs(y2-y3).lt.(1.e-05)) then

****** Je suis sur une horizontale ******

	xmoy=(x2+x3)/2.
	ymoy=y2
				    endif

	if ((abs(y2-y3).gt.(1.e-05)).and.
     % (abs(x2-x3).gt.(5.e-04))) then


****** Je suis sur une diagonale ******

	ymoy=(y2+y3)/2.
	xmoy=(x2+x3)/2.
				endif

ccc	write(*,*)xmoy,ymoy

      call dtcar(xmoy,ymoy,xmoy,ymoy,ymax,ymin,
     %n1,n2,h1,h2,not1,not2)

c      if (i.eq.639) write(*,*) i,not1,not2,nota
c      if (i.eq.722) write(*,*) i,not1,not2,nota
c      if (i.eq.762) write(*,*) i,not1,not2,nota
c      nota=not
      if(not.eq.not1) goto 26
      if(not.eq.not2) not=not1
      goto 27
26    not=not2
27    do 10 k=1,3
      k1=k+1
      k2=k+2
      if(k.eq.2) k2=1
      if(k.eq.3) k2=2
      if(k.eq.3) k1=1
      if((m1(k,not).eq.i2).and.(m1(k1,not).eq.i3)) goto11
      if((m1(k,not).eq.i3).and.(m1(k1,not).eq.i2)) goto11
10    continue
11    continue
      i4=m1(k2,not)
c     ========================================
c     recherche du cote sur lequel on arrive 
c     ========================================
      ps4=psi(i4)
	if (ilect.eq.0) ps4=-psi(i4)*1000./(r*r*b0)
      if((ps4-ps)*(ps-ps2).gt.0.) goto12
      if((ps4-ps)*(ps-ps3).ge.0.)goto13
c     ========================================
c     erreur
c     ========================================
c100   format ('erreur',i4,1x,i4,1x,i2)
      stop
c     ========================================
c     CALCUL DES 6 COORDONNEES DES SOMMETS DU TRIANGLE
c     ========================================                
12    continue
      i1=i3
      i3=i4
      x1=x3
      y1=y3
      x3=coox(i4)
      y3=cooy(i4)
      ps1=ps3
      ps3=ps4
      npt=not-(not/(2*n1))*2*n1
      if(npt.eq.0) goto40
      if(npt.eq.(2*n1-1)) goto50
      if(npt.eq.1) goto55
      goto14
40      if(k-2) 41,42,43
41      x2=2*pi
        x1=2*pi
       x=2*pi
       goto 14
42       x1=2*pi
       x3=2*pi
       goto14
43      x2=2*pi
       x3=2*pi
       goto 14
50      if(k.lt.3) goto 14
       x3=2*pi
       goto 14
55      if(k.lt.3) goto14
       x=0.
      x1=0.
      x2=0.
      goto14
13    continue
      i1=i2
      i2=i4
      x1=x2
      y1=y2
      x2=coox(i4)
      y2=cooy(i4)
      ps1=ps2
      ps2=ps4
      npt=not-(not/(2*n1))*2*n1
      if(npt.eq.0) goto60   
      if(npt.eq.(2*n1-1)) goto70
      if(npt.eq.1) goto75
      goto14
60      if(k-2) 61,62,63
61      x3=2*pi
       x1=2*pi
       x=2*pi
       goto14
62      x2=2*pi
      x3=2*pi
      goto 14
63      x1=2*pi
        x2=2*pi
       goto 14
70      if(k.lt.3) goto14
      x2=2*pi
      goto14
75     if(k.lt.3) goto14
      x=0.
      x1=0.
      x3=0.
14    continue
      goto15
5      continue
      if((abs(y-ymin).lt.eps).or.(abs(y-ymax).lt.eps)) then
      isort=1
       endif
       goto7
20     continue
c     ==========================================================
c     on s arrete au milieu d un triangle (cause: temps depasse)
c     =========================================================
      xc=x+(xc-x)*(ttol+dt-tmin)/dt
      yc=y+(yc-y)*(ttol+dt-tmin)/dt
      itria=itria+1      
c!!!      if(itria.gt.nemax) then
c!!!      print 114,nemax
	
c!!!      open (4,file='probleme.lis',status='UNKNOWN')
c!!!      write(4,*) ('noeud de depart:'),i
c!!!      write(4,*) ('heure :'),tmin,tmax
c!!!      	do 1988 k=1,nemax
c!!!      write(4,*) ('triangle traverse:'),ntria(k)
c!!!      write(4,*)('Coordonnees du point d,arrivee sur
c!!!     % le triangle'),xn(k),yn(k)
c!!!1988 	continue      
c!!!      close(4)
c!!!114   format ('LE NOMBRE D ETAPES EST SUPERIEUR AU MAX:',i4)
c!!!      STOP
c!!!      endif
      dt=ttol+dt-tmin
      xn(itria)=xc
      yn(itria)=yc
      temp(itria)=dt
      ntria(itria)=not
7     call intcar(ilect,i,itria,ntria,isort,xn,yn,ane3,ener3,
     %temp,tmax,
     %   iint,psi,metod,ndl,
     %   nt,m1,m2,cooy,coox,T,ane1,ener1,uo,u1)

1     continue
110   format ('NOEUD NUMERO:',i4,/,'NOMBRE TOTAL DE TRIANGLES
     % TRAVERSES:
     %',i2,/,'INDICE DE SORTIE (=1 SI SORTIE DU DOMAINE):',i2)
111   format (16x,'POUR CHAQUE TRIANGLE:SON No,LES COORDONNEES
     % DU POINT 
     %D ARRIVEE')
112   format (4x,'TRIANGLE',4x,':',5x,'XN',5x,':',5x,'YN',5x,
     %':',4x,'TEMP
     %S',3x,':')
113   format (6x,i4,6x,':',3(e12.5,':'))  
csss	write(*,*) 'Le nbre max de tri et depas :',nmax1,depas
      return
      end

c==================================================================
c==================================================================

      subroutine dtcar(xm,ym,xp,yp,ymax,ymin,n1,n2,h1,h2,not1,not2)
      common/triang1/n21,n22,n23,h21,h22,h23,yint1,yint2,raison
      common/triang2/almax,almin,rlpas,hpas1,hpas2
c     ========================================
c     DETERMINATION DU TRIANGLE CONTENANT LE POINT DE COORDONNEES
c     XM ET YM ;SI CE POINT N'APPARTIEND PAS AU DOMAINE ,ON
c     DETERMINERA
c     D'ABORD LE POINT DE LA FRONTIERE SITUE SUR LA DROITE MP
c     ET L 'INSTANT 
c     DE PASSAGE SUR CETTE FRONTIERE
c     ========================================
c	indic=0
      eps=1.e-05
      eps1=5.e-04
      eps2=1.e-05
      if(ym.lt.ymin) goto1
      if(ym.gt.ymax) goto2
      goto3
1     continue
c     ========================================
c     on a coupe la frontiere en un point ou ym=0
c     ========================================
      xq=xp-(xm-xp)/(ym-yp)*(yp-ymin)
      ym=ymin
      xm=xq      
      goto3
2     continue
c     ========================================
c     on a coupe la frontiere en un point ou ym=ymax
c     ========================================
      xq=xp+(xm-xp)/(ym-yp)*(ymax-yp)
      ym=ymax
      xm=xq
3     continue
      pi=4.*atan(1.)
      xm=amod(xm,2*pi)
      if(xm.lt.(0.0)) xm=2*pi+xm
      l1=int((xm/h1))
      call couche(ym,l2,ymin1,ymax1,ymin,ymax)
c      if (i.eq.(639)) write(*,*) i,l1,xm
c      if(l1.eq.n1) l1=l1-1
c      if (i.eq.(639)) write(*,*) i,l1,xm
      if (l1.eq.n1) then
      			l1=l1-1
		    endif
      s=-ym+ymin1+(xm-(l1+1)*h1)*(ymin1-ymax1)/h1
      if(s.le.+eps) goto 4
      not1=(2*n1)*l2+2*l1+1
      if(l1.eq.0) goto 5
      if(abs(xm-l1*h1).lt.(eps1))then
      not2=(2*n1)*l2+2*l1
c      indic=1
      endif
c      if ((not2.eq.nota)) goto 6
c      if (indic.eq.1)     return
      goto 6
5      if(abs(xm-l1*h1).lt.(eps1))then
        not2=(2*n1)*(l2+1)
c 	indic=1
	endif
c        if (not2.eq.nota) goto 6
c       if (indic.eq.1)   return
6      if(l2.eq.0) return
      if(abs(ym-ymin1).lt.eps2) not2=not1-2*n1+1
c      if (i.eq.(762)) write(*,*) i,l1,l2
      return
4     continue 
      if (abs(s).lt.eps)then
     	not2=(2*n1)*l2+2*l1+1
     			endif
      not1=(2*n1)*l2+2*(l1+1)
      if(abs(ym-ymax1).lt.eps2)then
       not2=not1+2*n1-1
c       indic=1
      endif
c       if(not2.eq.nota) goto 18
c       if (indic.eq.1)    return
c 18    if (abs(s).lt.eps) then 
c       not2=(2*n1)*l2+2*l1+1
c                          endif      
      if(abs(xm-(l1+1)*h1).lt.(eps1)) then
c      not3=(2*n1)*l2+2*l1+3
	not2=(2*n1)*l2+2*l1+3
c       if (not3.ne.nota) not2=not3
                                   endif
c      open(4,file='result.lis',status='UNKNOWN')
c      if (i.eq.(639)) write(4,*) l2,l1,not1,not2,xm,ym,eps
c      close(4)
c      if(i.eq.(762)) write(*,*)i,l1,l2
      return
      end

c===================================================================
c===================================================================

      subroutine couche(y,l2,ymin1,ymax1,ymin,ymax)
      common/triang1/n21,n22,n23,h21,h22,h23,yint1,yint2,raison
      common/triang2/almax,almin,rlpas,hpas1,hpas2
c      ============================================================
c      l2    :  Numero de la couche
c                  ymin |
c                         0
c               ymin+h21|
c                         1
c            ymin+2*h21 |
c                     ETC......
c      ymin1 :  la valeur inferieur de y sur la couche
c      ymax1 :  la valeur superiuere de y sur la couche
c      =============================================================
      l2=int(((-1/y)+(1/ymin))/rlpas)
      if(l2.eq.n21) l2=l2-1
      ymin1=almin-l2*rlpas
      ymax1=ymin1-rlpas
      ymin1=1./ymin1
      ymax1=1./ymax1
      return
      end

      				
c ==================================================================
c ==================================================================

      subroutine intcar(ilect,i,itria,ntria,isort,xn,yn,
     %  ane3,ener3,temp,tmax,iint,u
     %  ,metod,ndl,nt,m1,m2,cooy,coox,T,ane1,ener1,uo,u1)
	real*4 Nmax
      common/donnees/e,q,em,ab0,deltat,r,b0,emax,ymax,ymin,h2,h1,n1,n2
     %,Nmax
      common/triang1/n21,n22,n23,h21,h22,h23,yint1,yint2,raison
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      dimension ntria(10),xn(10),uo(1),u1(1)
     %   ,yn(10),ener3(1),ane3(1),temp(10),u(1),ener1(1),ane1(1)
      dimension m1(3,nt),m2(1),cooy(1),coox(1),T(1)

      pi=atan(1.)*4.
c     ========================================
c     integration de l equation totale le long des caracteristiques
c     ========================================
c1004   format(1x,'traitement du noeud',i4)
      heure=tmax
      ener3(i)=0.
      ane3(i)=0.
      if(itria.gt.0) goto 84
      itria=1
      xn(1)=coox(i)
      yn(1)=cooy(i)
c      nota=0
      call dtcar(xn(1),yn(1),xn(1),yn(1),ymax,ymin,n1,n2
     %,h1,h2,not,not2)
      j=1
      temp(1)=0.
      x1=coox(m1(1,not))
      x2=coox(m1(2,not))
      x3=coox(m1(3,not))
      y1=cooy(m1(1,not))
      y2=cooy(m1(2,not))
      y3=cooy(m1(3,not))
      npt=not-(not/(2*n1))*2*n1
      if(npt.eq.0) x1=2*pi
      if(npt.eq.0) x2=2*pi
      if(npt.eq.(2*n1-1)) x2=2*pi
      aire=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
      cb2=((xn(j)-x1)*(y3-y1)-(yn(j)-y1)*(x3-x1))/(aire)
      cb3=((yn(j)-y1)*(x2-x1)-(xn(j)-x1)*(y2-y1))/(aire)
      cb1=1-cb2-cb3
      goto 86
84      do 83 j=1,itria
      heure=heure-temp(j)
      not=ntria(j)
      x1=coox(m1(1,not))
      x2=coox(m1(2,not))
      x3=coox(m1(3,not))
      y1=cooy(m1(1,not))
      y2=cooy(m1(2,not))
      y3=cooy(m1(3,not))
      npt=not-(not/(2*n1))*2*n1
      if(npt.eq.0) x1=2*pi
      if(npt.eq.0) x2=2*pi
      if(npt.eq.(2*n1-1)) x2=2*pi
      aire=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
      cb2=((xn(j)-x1)*(y3-y1)-(yn(j)-y1)*(x3-x1))/(aire)
      cb3=((yn(j)-y1)*(x2-x1)-(xn(j)-x1)*(y2-y1))/(aire)
      cb1=1-cb2-cb3
c      
c PRISE EN COMPTE DES CONDITIONS AUX LIMITES
c
83    continue
86    j=itria
	if (ilect.eq.1)	then
      u1(i)=ainter(nt,cb1,cb2,cb3,not,uo,m1)
			endif
      if(isort.eq.0) goto 70
      if(isort.eq.1) goto 71
70    ener3(i)=ener3(i)+ainter(nt,cb1,cb2,cb3,not,ener1,m1)
      ane3(i)=ane3(i)+ainter(nt,cb1,cb2,cb3,not,ane1,m1)
      goto 82
71    continue
      ener3(i)=ener3(i)+ainte1(nt,cb1,cb2,cb3,not
     &    ,m1,cooy,coox,T,ener1)
      ane3(i)=ane3(i)+ainte1(nt,cb1,cb2,cb3,not
     %    ,m1,cooy,coox,T,ane1)
82    continue
       return
      end
      

c===============================================================
c===============================================================

      function ainter(nt,cb1,cb2,cb3,i,f,m1)
      dimension f(1),m1(3,nt)
	real*4 Nmax
      common/donnees/e,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,hB,n1,n2
     %,Nmax
      common/triang2/almax,almin,rlpas,hpas1,hpas2
csss	open(4,file='Rni1111.lis',status='UNKNOWN')
csss	write(4,*)(f(k),k=1,80)
csss	close(4)
      ainter=cb1*f(m1(1,i))+cb2*f(m1(2,i))+cb3*f(m1(3,i))
      return
      end
      
c===================================================================
c===================================================================

      function ainte1(nt,c1,c2,c3,i,m1,A,B,T,ener2)
	real*4 Nmax
      common/donnees/e,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,hB,n1,n2
     %,Nmax
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      dimension A(1),B(1),T(1),m1(3,nt),ener2(1)
      dimension c(3)

      c(1)=c1
      c(2)=c2
      c(3)=c3
      eps=1.e-06
      s=0.
      do 1 j=1,3
      if(c(j).lt.eps) goto1
      m=m1(j,i)
      s=s+c(j)*clndep0(m,A(m),ener2(m),B(m),T(m))
1     continue
      ainte1=s
      return
      end
      
c=====================================================================
c=====================================================================
      
      function clndep0(i,A,e1,B,T)
	real*4 Nmax
      common/donnees/e,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,hB,n1,n2
     %,Nmax
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      
      eps=1.e-06
      if(abs(A-Amin).le.hpas1+eps) goto 5
      if(abs(A-Amax).le.hpas2+eps) goto 6 
      goto 1
6     clndep0=0.
      return
5     clndep0=e1
      return

c!!!1     print 100,T,i
1	continue
100   format(1x,'ERREUR DANS LES CONDITIONS AUX LIMITES DE
     & PARTICULES0:',/,
     &    'la valeur de teta=',e12.5,'noeud',i4)
      clndep0=0.*rnorm
      return
      end
      
c==================================================================
c==================================================================

      subroutine ivol1(ane1,ener1,ane2,ener2,u,metod,ndl,
     %nt,m1,m2,A,phi
     %   ,b1,b2,b3,b4,b5,b6,ifact)
c     ========================================
c     ETAPE EVOLUTIVE PASSAGE DE t A t+deltat
c     EN TENANT COMPTE DES NON LINEARITES
c     ======================================== 
	real*4 Nmax     
      common/donnees/e,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,hB,n1,n2
     %,Nmax
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      dimension b1(1),b2(1),b3(1),b4(1),b5(1),b6(1)
      dimension ener1(1),ener2(1),ane1(1),ane2(1),u(1)
      dimension m1(3,nt),m2(1),A(1),phi(1)
	integer*4 ifact
c     ==============================================
c     metod : methode de resolution
c             = 3 :DECOMPOSITION  INTEGRATION EXACTE DU SECOND MEMBRE
c             = 4 :       "       INTEGRATION APPROCHEE      "
c     ================================================
      do 1 i=1,ndl
      b1(i)=0.
      b2(i)=0.
      b3(i)=0.
      b4(i)=0.
1     continue
c     ========================================
c     CALCUL DU SECOND MEMBRE DE L ENERGIE
c       1. dans b1   K*(a**4)/sqrt(C+a)*u.Vv/v
c       2. dans b2   la partie en E**3/2
c       3. dans b3   q(n)-u.Vq(n)*deltat
c       4. dans b4   la surface du support de wi
c     ========================================      
      pi=4.*atan(1.)
      do 2 i=1,nt
      u12=u(m1(1,i))-u(m1(2,i))
      u13=u(m1(1,i))-u(m1(3,i))
      A12=A(m1(1,i))-A(m1(2,i))
      A13=A(m1(1,i))-A(m1(3,i))
      phi12=phi(m1(1,i))-phi(m1(2,i))
      phi13=phi(m1(1,i))-phi(m1(3,i))
      if((i-2*n1*(i/(2*n1))).eq.(2*n1-1)) phi12=phi(m1(1,i))-2*pi
      if((i-2*n1*(i/(2*n1))).eq.0) phi13=2*pi-phi(m1(3,i))
      aire=phi12*A13-phi13*A12
      do 3 j=1,3
      noeud=m1(j,i)
c     b1(noeud)=b1(noeud)-A12*u13+A13*u12
      b4(noeud)=b4(noeud)+aire
 3     continue
 2     continue
      call intmas(ener1,b5,ndl,nt,m1,A,phi)
      call intmas(ane1,b6,ndl,nt,m1,A,phi)

      do 4 i=n1+1,ndl

ccc      b1(i)=+b1(i)*8./3.*deltat/(A(i)*b4(i))

       b1(i)=0.
       b2(i)=deltat*em/ab0*((r*b0*A(i))**4)/(6*b0*sqrt(4.-3.*A(i)))
       b2(i)=ifact*b2(i)
       b2(i)=sqrt(emax)*b2(i)*(A(i)/Amin)**(4./3.)

      if(metod.eq.4) goto 14
      if(metod-2) 10,11,12
c      ==================
c     metod=4
c     ===================
14    b3(i)=ener1(i)
      ane2(i)=ane1(i)
      goto 13
c     ===================
c     metod=1
c     ===================
10    b3(i)=6*(b5(i)-ifact*deltat*b3(i))/b4(i)
      ane2(i)=6*(b6(i)-ifact*deltat*ane2(i))/b4(i)
      goto 13
c     ====================
c     metod=2
c     =====================
11    b3(i)=ener1(i)-6*ifact*deltat*b3(i)/b4(i)
      ane2(i)=ane1(i)-6*ifact*deltat*ane2(i)/b4(i)
      goto 13
c     ======================
c     metod=3
c     ======================
12    b3(i)=6*b5(i)/b4(i)
      ane2(i)=6*b6(i)/b4(i)
13    ener2(i)=zero(b1(i),b2(i),b3(i),i)
      if (abs(ane2(i)).lt.((1.e-05)*(ane2(1)))) then 
						ane2(i)=0.
						ener2(i)=0.
				 	     	endif
      ane2(i)=ane2(i)/(1+em*sqrt(emax)*((A(i)/Amin)**(4./3.))
     % /(2*ab0)*((r*b0*A(i))
     %    **4)/(b0*sqrt(4.-3.*A(i)))*ifact*deltat*sqrt(ener2(i)))
      
4     continue
      return
      end
 
c==================================================================
c==================================================================

      subroutine intmas(u,v,ndl,nt,m1,A,B)
	real*4 Nmax
      common/donnees/e,q,em,ab0,deltat,r,b0,emax,ymax,ymin,hA,hB,n1
     %,n2,Nmax
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      dimension u(1),v(1),A(1),B(1)
      dimension m1(3,nt)
      do 1 i=1,ndl
1     v(i)=0.
      pi=4*atan(1.)
      do 2 i=1,nt
      A12=A(m1(1,i))-A(m1(2,i))
      A13=A(m1(1,i))-A(m1(3,i))
      B12=B(m1(1,i))-B(m1(2,i))
      B13=B(m1(1,i))-B(m1(3,i))
      if((i-2*n1*(i/(2*n1))).eq.(2*n1-1)) B12=B(m1(1,i))-2*pi
      if((i-2*n1*(i/(2*n1))).eq.0) B13=2*pi-B(m1(3,i))
      aire=B12*A13-B13*A12
      do 3 j=1,3
      noeud=m1(j,i)
      j1=j-1
      if(j1.lt.1)goto4
      do 5 k=1,j1
      v(noeud)=v(noeud)+aire/24.*u(m1(k,i))
      v(m1(k,i))=v(m1(k,i))+aire/24.*u(noeud)
5     continue
4     v(noeud)=v(noeud)+2./24.*u(noeud)*aire
      v(noeud)=v(noeud)
3     continue
2     continue
      return
      end

c==================================================================
c==================================================================

      function zero(a,b,c,i)
c     ========================================
c     Recherche du zero de la function 
c            f(x)=x+a*x+b*(x**3./2.)-c
c     Par une methode de newton
c     ========================================
c     ========================================
c     Test c>0
c     ========================================
       eps=1.e-5
      eps1=1.e-06
      if(c.lt.eps1) goto 1
        h=-c
        hp=1+a

       if(hp.lt.0) goto5
        l=0
        x=0.
3     continue
      l=l+1
      x=x-h/hp
      h=(a+1)*x+b*(x**(3./2.))-c
      if(abs(h).lt.eps)goto 4
      hp=(1+a)+b*3./2.*(x**0.5)
      if(abs(hp).gt.eps) goto 3
5     if(c.lt.eps) goto4
      if(l.gt.10) goto7
      if(l.ne.0) goto 6
      x=(-hp/b)**2
      h=(1+a)*x+b*(x**(3./2.))-c
      hp=1+a+3/2.*b*(x**0.5)
      goto 3
c!!!6     print 103,hp
c!!!103   format (1x,'ERREUR :LA DERIVEE S ANNULE',e12.5)
c!!!      stop

6	continue
	zero=0
	return

4     continue
c100   format ('LA SOLUTION EST OBTENUE EN ',i3,'ITERATIONS ET VAUT
c      & :',e
c     &  12.5,'PRECISION:',e12.5)
      zero=x
      return
c!!!7     print 104,h,hp,x
c!!!104   format(1x,'h,hp,x',3e12.5)

7	continue

      return
1     continue
      eps2=-eps1
      if(c.gt.eps2) goto25
c!!!      print 102,c,i
c!!!102   format ('ERREUR LE COEFFICIENT C EST NEGATIF :',e12.5,i4)
25     continue
       zero=0
      return
      end     

c===================================================================
c===================================================================

	subroutine ivol5(Tic,Rnic,nt,Ti,Rni,ndl,hh,hh1,hh2,
     %  aire1,RNNN,TV,Ti1,Rni1,A,B,T,m1,m2,iinti,ibrupt,heur1,
     %  Bfm,Bfs,delta2,heuref,rtmax,rnmax,b5,b6,znmax,
     %  ztmax,inum,heure,ifact,ie)
     	real*4 Nmax
     	common/donnees/e,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,hB,n1,
     %n2,Nmax
	common/triang2/almax,almin,rlpas,hpas1,hpas2
	real*4	NN1,NN2,NN3,N01,N10,N00
	real*4 NNc,TVc
	dimension Tic(1),Rnic(1),Ti(1),Rni(1),hh(1),hh1(1)
	dimension Ti1(1),Rni1(1),A(1),B(1),T(1),m1(3,nt),m2(1),hh2(1)
	dimension RNNN(1),TV(1),aire1(1),b5(1),b6(1),aire2(1)
	integer*4 inum(40)
	integer*4 iflag1,ifact,ie
	pi=4.*atan(1.)
	eps=1.e-02

csss	write(*,*)'Je suis en Ivol5'
	
cccc	REMARQUE FONDAMENTALE: LES VARIABLES SONT NORMALISEES

c 	CALCUL DE QUELQUES CONSTANTES DE NORMALISATION DU PB

	Volmax=(1./(A(1)))**(4.)
	Rlmax=1/A(1)
	Fmax=1./(2.*sqrt(4-3.*A(1))*r/((A(1))**4.))
	Fmax=Fmax*em/sqrt(1836.)
	Fmax=Fmax*sqrt(rtmax)*(Rlmax**(16./3.))

***
****
cHAO97	On multiplie la diffusion forte en angle d'attaque
cHAO97  par 0.2 pour avoir des flux en energie du bon ordre
cHAO97  de grandeur.

	Fmax=Fmax*1.

cHAO97
****
***
	Femax=Fmax*4./3.
	C1max=rtmax*4000.*(Rlmax**(8./3.))/(r*r*b0)
	C2max=5000.*rtmax*(Rlmax**(8./3.))/(3.*r*r*b0)
	
c 	INITIALISATION DE ZNMAX ET ZTMAX QUI SONT LES VALEURS MAX
c 	DE Z POUR N ET T.

	znmax=0.
	ztmax=0.
	
	iflag2=0
	do 4000 i=ndl,1,-1
	if (Rnic(i).gt.(Rni(1)/1000.))	then
		if (iflag2.eq.0)	then
			icom1=i
			iflag2=1
					endif
					endif
 4000	continue 
 	
	if (icom1.ge.(ndl-n1+1))	then
		icom1=ndl
					else
		nqi=icom1-(icom1/n1)*n1
		if (nqi.eq.0)	then
			icom1=icom1+1
				else
			icom1=((icom1/n1)+1)*n1+1
				endif
					endif
	do 4010 i=1,(1*n1)
	Tic(i)=Ti(i)
	Rnic(i)=Rni(i)
 4010 	continue
 
	do 4020 i=icom1,ndl

	if (Rnic(i).lt.((1.e-03)*Rni(1))) then
	Rnic(i)=0.
	Tic(i)=0.
					  endif
 4020	continue
	
c ***	On calcule le T chapeau moyen et le N chapeau moyen comme 
c 	Laure pour les electrons ***

	call intmas(Rnic,b5,ndl,nt,m1,A,B)
	call intmas(Tic,b6,ndl,nt,m1,A,B)
  	
	do 6663 i=1,ndl
	Ti1(i)=Ti(i)
	Rni1(i)=Rni(i)
 6663	continue
  
c ***	Initialisation des parametres ***

	do 100 i=1,ndl
	hh(i)=0.
	hh1(i)=0.
	hh2(i)=0.
	aire1(i)=0.
 100 	continue

 	do 2 i=1,nt 
	T1=Ti(m1(1,i))
	T2=Ti(m1(2,i))
	T3=Ti(m1(3,i))
	NN1=Rni(m1(1,i))
	NN2=Rni(m1(2,i))
	NN3=Rni(m1(3,i))
	x1=B(m1(1,i))
	x2=B(m1(2,i))
	x3=B(m1(3,i))
	y1=A(m1(1,i))
	y2=A(m1(2,i))
	y3=A(m1(3,i))
	nqi=i-(i/(2*n1))*2*n1
	if (nqi.eq.0) x1=2*pi
	if (nqi.eq.0) x2=2*pi
	if (nqi.eq.(2*n1-1)) x2=2*pi
	delta=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
	call Coeff(x1,x2,x3,y1,y2,y3,NN1,NN2,NN3,N10,N01,N00)
	call Coeff(x1,x2,x3,y1,y2,y3,T1,T2,T3,T10,T01,T00)
	xc=(x1+x2+x3)/3.
	yc=(y1+y2+y3)/3.
	z=N10*T01-N01*T10
	z1=T10
	z2=N10
	
c ***	On fait une boucle pour calculer une moyenne sur le noeud
c 	a partir des triangles environnant ***

	do 3 j=1,3
	noeud=m1(j,i)
		
	hh(noeud)=hh(noeud)+z*delta
	hh1(noeud)=hh1(noeud)+z1*delta
	hh2(noeud)=hh2(noeud)+z2*delta
	aire1(noeud)=aire1(noeud)+delta
 3	continue
 2	continue
 
c ***	On fait la moyenne sur h **

	do 4 noeud=1,ndl
	hh(noeud)=hh(noeud)/aire1(noeud)
	hh1(noeud)=hh1(noeud)/aire1(noeud)
	hh2(noeud)=hh2(noeud)/aire1(noeud)
 4	continue

c **** On rajoute cela pour que les particules ne franchissent pas
c ****
c **** la derniere equipotentielle fermee ****


	do 4030 i=1,icom1
		iflag2=0
		do 4040 j=1,n1
c%%%
		inum(j)=0
c%%%
		if (i.eq.inum(j))	then
		iflag2=1
					endif
 4040	continue
		
		if (abs(Rnic(i)).lt.((1.e-06)*Rni(1))) 	then
		Rnic(i)=0.
		Tic(i)=0.
						endif
		if (iflag2.eq.(1)) goto 4030

		if (abs(Rnic(i)).ge.((1.e-06)*Rni(1)))	then
		Rnic(i)=6.*b5(i)/aire1(i)
		Tic(i)=6.*b6(i)/aire1(i)
						endif

 
 4030	continue

	if ((ie-(ie/ifact)*ifact).eq.0)	goto 1993

	do 1992 i=n1+1,ndl
	Rni(i)=Rnic(i)
	Ti(i)=Tic(i)
	if (abs(Rni(i)).lt.((1.e-03)*Rni(1)))		then
						Ti(i)=0.
						Rni(i)=0.
							endif
 1992	continue
 
	goto 1994

 1993	continue
  
c ***	On resoud maintenant les equations ***

	do 5 i=n1+1,ndl
	
c ***	Arret quand probleme ***

	if ((Ti(i)).lt.(0.))	then
		open(4,file='Neg.lis',access='sequential',position='append')
		write(*,*)'temp :',Ti(i)
		write(*,*)'noeud :',i
		close(4)
		STOP
				endif
				
	if ((Rni(i)).lt.(0.))	then
		open(4,file='Neg.lis',access='sequential',position='append')
		write(*,*)'dens :',Rni(i)
		write(*,*)'noeud :',i
		close(4)
		STOP
				endif
	
c ***	Calcul de N(n+1) ***

c *** 	Calcul de f(T(n)) ***

c***	On calcule ici les flux de precipitation

	z=-Fmax
	z=z*sqrt(Ti(i))
	z=z*(A(i))**(16./3.)
	

c ***	On tient compte ici des flux orthogonaux ***

	if (Rni(i).ge.(Rni(1)*(1.e-06))) then
	w1=(A(i))**(5./3.)
	w1=w1*(hh1(i)*Rni(i)+hh2(i)*Ti(i))
	w1=w1*C1max
c	w1=0.
	z1=w1/Rni(i)
	z1=abs(z1)
	z2=z
	z=z+w1/Rni(i)
					   else
	z=0.
	z1=0.
	z2=0.
					   endif
	
	if ((abs(z)).gt.znmax) znmax=abs(z)
					   
	
c ***	Calcul de N chapeau ***

c	NNc=6.*Rnic(i)/aire1(i)
	NNc=Rnic(i)

c ***	Calcul de NN  l'instant suivant ***

c ***	On calcule NN en tenant compte du fait que les gradients sont
c ***	bien calcules la ou il faut ***

c ***	Ceci est la formule initiale du calcul de RNNN(i) dite implicite

c	RNNN(i)=NNc/(1-ifact*deltat*z)
	
c ***	Ceci est la formule testee du calcul de RNNN(i) dite explicite

	RNNN(i)=Rni(i)*z*ifact*deltat + NNc
					   				
c ***	Calcul de T*V**(2/3)(n+1) ***

c ***	Calcul de g(T(n)) ***
	
	if (Rni(i).lt.((1.e-06)*Rni(1))) 		then
		z=0.
		z1=0.
		z2=0.
							else
							
c ***	On calcule les flux de precipitation

		z=-(A(i))**(16./3.)
		z=z*(Ti(i))**(1./2)
		z=z*Rni(i)
		z=z/4.
		z=z*Femax
	
c ***	On tient compte ici des flux orthogonaux ***

		w1=-hh(i)
		w1=w1*C2max
		w1=w1*((A(i))**(8./3.))
		w1=w1
c		w1=0.
		w2=28.*hh1(i)*Rni(i)/5.
c		w2=0.
		w2=w2-16.*hh2(i)*Ti(i)/15.
		w2=w2*C2max
		w2=w2
		w2=w2*((A(i))**(5./3.))
c		w2=0.
		z1=(w1+w2)/Rni(i)
		z1=abs(z1)
		z2=z/Rni(i)
		z=z+w1+w2
		z=z/Rni(i)
							endif		
c ***	Calcul de kTV*2/3 chapeau ***

c	TVc=6.*Tic(i)/aire1(i)
	TVc=Tic(i)
	
c ***	Calcul de kTV*2/3 a l'instant t(n+1) ***
c ***	On fait de meme attention a l'endroit ou il faut calculer ***
c ***	gradients ***

c ***	Ceci est la formule initiale du calcul de TV(i) dite implicite

c	TV(i)=TVc/(1-ifact*deltat*z)
	
c ***	Ceci est la formule testee du calcul de TV(i) dite explicite

	TV(i)=Ti(i)*z*ifact*deltat + TVc

c ***	On coupe tout a 1.e-03 de Nmax ***

	if (abs(RNNN(i)).lt.((1.e-03)*Rni(1)))		then
							TV(i)=0.
							RNNN(i)=0.
							endif
	
	Rni1(i)=RNNN(i)	
	Ti1(i)=TV(i)
 300	continue
 
c ***	Prise en compte des conditions aux limites

 5	continue
 
	do 8 i=1,ndl
	Ti(i)=Ti1(i)
	Rni(i)=Rni1(i)
 8	enddo

 1994	continue
  	
 	return
	end

c===================================================================
c===================================================================

	subroutine GRADP(aire1,B,A,pp,nt,ndl,m1,m2,dbeta,dalpha,betaa)
	
c***	CETTE SUBROUTINE CALCULE LES GRADIENTS DE PRESSION DES ION
c       OU DES
c 	ELECTRONS EN CHAQUE NOEUD. ELLE CALCULE LES GRADIENTS P/R A
c       r=L*a 
c 	ET beta ***

	real*4 Nmax

	common/donnees/el,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,hB,
     %n1,n2,Nmax
	common/triang2/almax,almin,rlpas,hpas1,hpas2
	
	dimension aire1(1),B(1),A(1),pp(1),m2(1),m1(3,nt),
     %dbeta(1),dalpha(1)
	dimension betaa(1)
	pi=4.* atan(1.)
	do 100 i=1,ndl
	dalpha(i)=0.
	dbeta(i)=0.
	aire1(i)=0.
 100	continue
 	do 2 i=1,nt
	x1=B(m1(1,i))
	x2=B(m1(2,i))
	x3=B(m1(3,i))
	y1=A(m1(1,i))
	y2=A(m1(2,i))
	y3=A(m1(3,i))
	nqi=i-(i/(2*n1))*(2*n1)
	if (nqi.eq.0) x1=2*pi
	if (nqi.eq.0) x2=2*pi
	if (nqi.eq.(2*n1-1)) x2=2*pi
	call Coeff(x1,x2,x3,y1,y2,y3,pp(m1(1,i)),
     %pp(m1(2,i)),pp(m1(3,i)),
     %pp10,pp01,pp00)
     	delta=(x2-x1)*(y3-y1)-(x1-x3)*(y1-y2)
     	do 3 j=1,3
	noeud=m1(j,i)
	dalpha(noeud)=pp01*delta + dalpha(noeud)
	dbeta(noeud)=pp10*delta + dbeta(noeud)
	aire1(noeud)=aire1(noeud)+delta
 3	continue
 2	continue
 	do 4 noeud=1,ndl
	dalpha(noeud)=dalpha(noeud)/aire1(noeud)
	dalpha(noeud)=-dalpha(noeud)/r
	dalpha(noeud)=dalpha(noeud)*((A(noeud))**(2.))
	dbeta(noeud)=dbeta(noeud)/aire1(noeud)
	dbeta(noeud)=dbeta(noeud)/r
	betaa(noeud)=dbeta(noeud)
 4	continue
 	return
	end

c===================================================================
c===================================================================

 	subroutine TMT(Atemp,Btemp,Ttemp,
     %ddndl,ddnt,ddn1,ddn2,m1temp,
     %A,B,T,nndl,nnt,nn1,nn2,m1,tabTM)

	real*4 Atemp(1),Btemp(1),Ttemp(1)
	integer*4 ddndl,ddnt,ddn1,ddn2
	integer*4 m1temp(3,ddnt)
	real*4 A(1),B(1),T(1)
	real*4 tabTM(nndl,4)
	integer*4 m1(3,nnt)

	eps=1.e-09
	pi=4.*atan(1.)

	do 1510	i=nndl,1,-1

		iflag=0

		do 1520 j=ddnt,1,-1

		if (iflag.eq.1)	goto 1520

		x=B(i)
		y=A(i)

		not1=m1temp(1,j)
		not2=m1temp(2,j)
		not3=m1temp(3,j)
		y1=Atemp(not1)
		y2=Atemp(not2)
		y3=Atemp(not3)
		
		irest=j-(j/2)*2

		if (irest.eq.1)	then

c  *****	Triangle inferieur ***

		x1=Btemp(not1)
		jrest=j-(j/(2*ddn1))*(2*ddn1)

		if (jrest.eq.(2*ddn1-1))	then
		x2=2.*pi
						else
		x2=Btemp(not2)	
						endif

		x3=Btemp(not3)
		
		ymin=y1
		ymax=y3
		xmin=x1
		xmax=x2

				else


c  *****	Triangle superieur ***

		jrest=j-(j/(2*ddn1))*(2*ddn1)

		if (jrest.eq.0)			then
		x1=2.*pi
		x2=2.*pi
						else
		x1=Btemp(not1)
		x2=Btemp(not2)	
						endif

		x3=Btemp(not3)
		
		ymin=y1
		ymax=y3
		xmin=x3
		xmax=x2

				endif

		if (abs(y-ymin).lt.(1.e-05)) 	then
						y=ymin+eps
						endif

		if (abs(y-ymax).lt.(1.e-05)) 	then
						y=ymax-eps
						endif

		if (abs(x-xmin).lt.(5.e-04)) 	then
						x=xmin+eps
						endif

		if (abs(x-xmax).lt.(5.e-04)) 	then
						x=xmax-eps
						endif

		if (y.lt.ymin)	goto 1520
		if (y.gt.ymax)	goto 1520
		if (x.lt.xmin)	goto 1520
		if (x.gt.xmax)	goto 1520

*** Est-on sur le triangle inferieur ou le triangle superieur ?***

		if (irest.eq.1)	then

*** Triangle inferieur ***

		adiag=(y3-y2)/(x3-x2)
		bdiag=y3-adiag*x3

		sdiag=(y-adiag*x-bdiag)
		if (sdiag.le.(1.e-05)) iflag=1

				else

*** Triangle superieur ***

		adiag=(y3-y1)/(x3-x1)
		bdiag=y3-adiag*x3

		sdiag=(y-adiag*x-bdiag)
		if (sdiag.ge.(-1.e-05)) iflag=1

				endif

		if (iflag.eq.0) goto 1520

		aire=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
		cb2=((x-x1)*(y3-y1)-(y-y1)*(x3-x1))/aire
		cb2=abs(cb2)
		cb3=((y-y1)*(x2-x1)-(x-x1)*(y2-y1))/aire
		cb3=abs(cb3)
		cb1=1-cb2-cb3
		cb1=abs(cb1)

		tabTM(i,1)=j
		tabTM(i,2)=cb1
		tabTM(i,3)=cb2
		tabTM(i,4)=cb3

 1520 		continue

 1510	continue

	end

c===================================================================
c===================================================================

 	subroutine interpolation(Atemp,Btemp,Ttemp,ddndl,ddnt,ddn1,
     %ddn2,m1temp,A,B,T,nndl,nnt,nn1,nn2,m1,
     %utemp,u1temp,u2temp,u3temp,u4temp,u,u1,u2,u3,u4,tabTM)

	real*4 Atemp(1),Btemp(1),Ttemp(1)
	integer*4 ddndl,ddnt,ddn1,ddn2
	integer*4 m1temp(3,ddnt)
	real*4 A(1),B(1),T(1)
	real*4 u(1),u1(1),u2(1),u3(1),u4(1)
	real*4 utemp(1),u1temp(1),u2temp(1),u3temp(1),u4temp(1)
	real*4 tabTM(nndl,4)
	integer*4 m1(3,nnt)

	eps=1.e-09
	pi=4.*atan(1.)

	rnmax=abs(utemp(1))
	rnmax1=abs(u1temp(1))
	rnmax2=abs(u2temp(1))
	rnmax3=abs(u3temp(1))
	rnmax4=abs(u4temp(1))

	do 1 i=1,ddndl

	if (abs(utemp(i)).ge.rnmax)	rnmax=abs(utemp(i))
	if (abs(u1temp(i)).ge.rnmax1)	rnmax1=abs(u1temp(i))
	if (abs(u2temp(i)).ge.rnmax2)	rnmax2=abs(u2temp(i))
	if (abs(u3temp(i)).ge.rnmax3)	rnmax3=abs(u3temp(i))
	if (abs(u4temp(i)).ge.rnmax4)	rnmax4=abs(u4temp(i))

 1	continue
	
	do 2 i=1,ddndl

	if (abs(utemp(i)).le.(rnmax/1000.))	utemp(i)=0.
	if (abs(u1temp(i)).le.(rnmax1/1000.))	u1temp(i)=0.
	if (abs(u2temp(i)).le.(rnmax2/1000.))	u2temp(i)=0.
	if (abs(u3temp(i)).le.(rnmax3/1000.))	u3temp(i)=0.
	if (abs(u4temp(i)).le.(rnmax4/1000.))	u4temp(i)=0.

 2	continue

	do 10 i=1,nndl
	u(i)=0.
	u1(i)=0.
	u2(i)=0.
	u3(i)=0.
 10	continue

	do 1510	i=1,nndl
	
		j=tabTM(i,1)
		not1=m1temp(1,j)
		not2=m1temp(2,j)
		not3=m1temp(3,j)

		cb1=tabTM(i,2)
		cb2=tabTM(i,3)
		cb3=tabTM(i,4)
			
		u(i)=cb1*utemp(not1)
		u(i)=u(i)+cb2*utemp(not2)
		u(i)=u(i)+cb3*utemp(not3)	

		u1(i)=cb1*u1temp(not1)
		u1(i)=u1(i)+cb2*u1temp(not2)
		u1(i)=u1(i)+cb3*u1temp(not3)

		u2(i)=cb1*u2temp(not1)
		u2(i)=u2(i)+cb2*u2temp(not2)
		u2(i)=u2(i)+cb3*u2temp(not3)	

		u3(i)=cb1*u3temp(not1)
		u3(i)=u3(i)+cb2*u3temp(not2)
		u3(i)=u3(i)+cb3*u3temp(not3)

		u4(i)=cb1*u4temp(not1)
		u4(i)=u4(i)+cb2*u4temp(not2)
		u4(i)=u4(i)+cb3*u4temp(not3)

 1510	continue

	end

c===================================================================
c===================================================================

	subroutine Fourrier(ali,alifit,alii,rafou,iafou,
     & B,par1,n1,n2)

c!	On effectue une transformee de Fourrier en B.
c!	ali est la distribution que l'on veut Fourrieriser.
c!	Le resultat est mis dans alifit puis dans ali.
c!	alii est une variable intermediaire de calcul.
c!	rafou(i,k) et iafou(i,k) contiennent les coefficients de Fourrier
c!	respectivement reels et imaginaires de la transformee
c!	de Fourrier pour la ligne i a l'ordre k.

	integer*4 par33,par44,par55,par1
	real*4 ali(par1)
	real*4 alifit(par1),rafou(par1,0:5),iafou(par1,0:5)
     	real*4 alii(par1)
	real*4 B(par1)
	real*4 rmoy,imoy

	dpi=2.*4.*atan(1.)

	iordi=0
	iordf=5

	par33=n1
	par44=n2+1
	par55=n1+1

	do 5 kk=iordi,iordf
	
	if (kk.eq.0) 	then
 	pas=dpi/(par33*10.)
			else
	pas=dpi/(par33*10.*kk)
			endif

	do 100 i=1,par44
           rafou(i,kk)=0
           iafou(i,kk)=0
 100     continue

        do 200 i=1,par44

	   if (kk.eq.0)	then
	   ifin=par33*10
			else
	   ifin=par33*10*kk
 			endif

	do 250 i1=1,par33
	k=(i-1)*par33+i1
	alii(i1)=ali(k)
 250	continue

	alii(par55)=ali(k-(par33-1))

	do 300 j=1,ifin

***	extrapolation en longitude pour le calcul de la Fourier

		abs1=(j-1)*pas
		abs11=abs1+pas

		do 350	jj=1,par33
		
		eps=1.e-09
		rinf=(jj-1)*dpi/par33-eps
		rsup=rinf+dpi/par33+eps

		if ((abs1.ge.rinf).and.(abs1.le.rsup))	then		
		icompt=jj
							endif

 350		continue

		x1=(icompt-1)*dpi/par33
		x2=x1+dpi/par33
		rpot1=alii(icompt)
		rpot2=alii(icompt+1)
		ra=(rpot1-rpot2)/(x1-x2)
		rb=rpot1-ra*x1
		pot1=ra*abs1+rb
		pot11=ra*abs11+rb

              rmoy=(cos(-kk*abs1)*pot1+cos(-kk*abs11)*pot11)*pas/2.
              rafou(i,kk)=rafou(i,kk)+rmoy
              imoy=(sin(-kk*abs1)*pot1+sin(-kk*abs11)*pot11)*pas/2.
              iafou(i,kk)=iafou(i,kk)+imoy

ccc		if (i.eq.1)	then
ccc		write(*,*)'j',j
ccc		write(*,*)'abs1',abs1,' abs11',abs11
ccc		write(*,*)'pot1',pot1,' pot11',pot11
ccc		write(*,*)'rfou',rfou(i)*2./dpi
ccc		write(*,*)'ifou',ifou(i)*2./dpi
ccc				endif

 300         continue

              rafou(i,kk)=2.*rafou(i,kk)/dpi
              iafou(i,kk)=2.*iafou(i,kk)/dpi

	      if (kk.eq.0)	then
	      			rafou(i,kk)=rafou(i,kk)/2.
				iafou(i,kk)=iafou(i,kk)/2.
				endif

c!	write(*,*)kk,i,rafou(i,kk),iafou(i,kk)
				
 200     continue
       	
 5	continue

	do 6900 i=1,par1

	alifit(i)=0.

 6900	continue

c!	write(*,*)(ali(k),k=1,1240)
	
	do 7000 i=1,par44

		do 7100 j=1,par33

		k=(i-1)*par33+j

		x=B(k)
	
			do 7200 kk=iordi,iordf

			alifit(k)=alifit(k)+rafou(i,kk)*cos(kk*x)
			alifit(k)=alifit(k)-iafou(i,kk)*sin(kk*x)

 7200			continue

		ali(k)=alifit(k)

 7100		continue

 7000	continue

c!	write(*,*)(ali(k),k=1,1240)

	return

	end

