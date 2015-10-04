

c===================================================================
c===================================================================

      SUBROUTINE POTENT

      integer*4 par1,par2,par3,dpar1,dpar2,dpar3
      integer*4	dpar11,dpar22,dpar33
      integer*4 dn1,dn2,dnt,dndl
      parameter (par1=1240,par2=2400,par3=2481)
      parameter (dpar1=(25*80),dpar2=(24*160),dpar3=(25*160+1))
      parameter (dpar11=(49*80),dpar22=(48*160),dpar33=(49*160+1))

      real*4 a11(dpar1),a12(dpar1),
     %a21(dpar1),a22(dpar1),ao(dpar1),ff(dpar1),
     %alp1(dpar1),alp2(dpar1)
     %,a1(400000)
     %,TT(3),DPHI(3),Rnn(3),
     % utemp(dpar1),FEtemp(dpar1),Etemp(dpar1),alitem(dpar1)

      integer*4 ndimat
      real*4 aa(6)
      real*4 Nmax,mass
cHAO97      dimension sptemp(dpar1),shtemp(dpar1)
      common/donnees/el,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,
     $hB,n1,n2,Nmax
      common/triang1/n21,n22,n23,h21,h22,h23,Aint1,Aint2,raison
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      real*4 dTtemp(25)

	data dTtemp/71.97,69.668,67.279,64.797,62.221,59.548,
     %  56.783,53.932,51.004,48.013,44.977,41.916,38.854,35.816,
     %  32.828,29.916,27.104,24.414,21.864,19.471,17.246,15.197,
     %  13.329,11.645,10.14/

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

	common/conductances/sptemp(dpar1),shtemp(dpar1)

cHAO97

cesr

	common/dea/ideaele,ideaion,ideapot

cesr

c!!!
c!!!	Permet Lecture des conductivites de Barbara
c!!!
	real*4 Barped(dpar11),Barhal(dpar11)

	
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
c             muns : pointeur des coefficients diagonaux de la
c		     matrice de
c                    resolution du probleme elliptique
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

c!!!	write(*,*)'INITIALISATION'

	open(4,file='sdm0.lis',status='unknown')
	close(4)

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

c     ========================================
c     CONSTRUCTION D UNE TRIANGULATION UNIFORME SUR LE DISQUE
c     ========================================
      call disc1(dn2,dn1,m1temp,m2temp,dnt,0)

      r=6.37e+06
      b0=3.08e-5
      Bmin=0.
      Bmax=2*pi
      Tmin=(90.-71.97)*pi/180.
      teta0=Tmin
    
      almax=1./cos(10.14*pi/180.)
      almax=almax*almax
      Amax=1./almax
      Tmax=asin(sqrt(Amax))
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
      
c     ========================================
c     CALCUL DES COEFFICIENTS DIAGONAUX
c     ========================================
      call detmns(dndl,dnt,m1temp,muns,3)
c!!!  open(4,file='pointeur.lis',status='UNKNOWN')
c!!!  write(4,*) (muns(i),i=1,(dndl+1))
c!!!  close(4)
      
      dndl1=dndl+1
c     print 100,n1,n2
ccc	open(4,file='disque.lis',status='UNKNOWN')
ccc     write (4,*) ndl,nt,n1,n2,(T(i),i=1,ndl),(B(i),i=1,ndl)
ccc     &       ,hA,hB,Bmax,Bmin,Amin,Amax,((m1(j,i),j=1,3),i=1,nt)
ccc     &       ,(m2(i),i=1,ndl),(A(i),i=1,ndl),(muns(i),i=1,ndl1)
ccc     %   ,Aint1,Aint2,n21,n22,n23,h21,h22,h23,raison
ccc  	close(4)
c     ========================================
c     ENTREE DES MATRICES ELEMENTAIRES DE RIGIDITE
c       al1 : dwi/dx*dwj/dx
c       al2 : dwi/dy*dwj/dy
c       al3 : dwi/dy*dwj/dx
c       al4 : dwi/dx*dwj/dy
c     ========================================      
      al1(1)=1.
      al1(2)=-1.
      al1(3)=0.
      al1(4)=-1.
      al1(5)=1.
      al1(6)=0.
      al1(7)=0.
      al1(8)=0.
      al1(9)=0.

      al2(1)=1.
      al2(2)=0.
      al2(3)=-1.
      al2(4)=0.
      al2(5)=0.
      al2(6)=0.
      al2(7)=-1.
      al2(8)=0.
      al2(9)=1.
      
      al3(1)=1.
      al3(2)=0.
      al3(3)=-1.
      al3(4)=-1.
      al3(5)=0.
      al3(6)=1.
      al3(7)=0.
      al3(8)=0.
      al3(9)=0.
      
      al4(1)=1.
      al4(2)=-1.
      al4(3)=0.
      al4(4)=0.
      al4(5)=0.
      al4(6)=0.
      al4(7)=-1.
      al4(8)=1.
      al4(9)=0.
      
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
ccc     Nmax=rnmax*r/b0
ccc     Nmax=Nmax/((A(1))**(4.))
c!!!	Nmax=9999.
c       emax=0.4

	em=sqrt(8*el/(pi*mass))

cHAO97	On multiplie par un coefficient pour tenir compte d'une
cHAO97	non diffusion forte en angle d'attaque pour la precipitation
cHAO97	d'electrons (voir Schumaker et al. (1989), JGR, 10061).

cHAO97	Ici le coefficient est mis a 0.4

	em=em*1.

cHAO97	

       ab0=r*r*r*r*r*b0*b0*b0

c      1.  Valeur optimale de Sp,le jour:Spjour
c      2.  rapport entre Sh et Sp       : rapp=sh/sp
c      3.  Rapport entre Spjor et Spnuit: rapp1=spjour/spnuit
c                    mis egal a 30 
c      4.  choix d'interpolation pour la loi des conductivites
c                  iopsp=1   lineaire
c                  iopsp=2   quadratique
c      5.  Option d'ecriture du modele de conductivite
c               iopec=0,pas de sortie
c                iopec=1,sortie des valeurs de sp et sh
c      6.  Option d'interpolation       : iopint
c               ioptint =0 , pas d'interpolation
c               ioptint =1 , interpolation sur trois points
c               ioptint =2 , interpolation sur cinq points
c

      spjour=26
      rapp=1.7
      rapp1=30.
      iopsp=1
      iopec=0
      ioptint=1
     
       call tencon1(spjour,rapp,iopec,sp0temp,sh0temp,
     %ioptint,rapp1,dndl,Btemp,Ttemp,dn1)

       do 3 i=1,dndl
       sptemp(i)=sp0temp(i)/1.
 3     shtemp(i)=sh0temp(i)/1.

c!!!	write(*,*)'FEtem2 Etem2 au debut',(FEtem2(i),
c!!!     %Etem2(i),i=1,dpar11)

C ******************************************************* FIN

							endif

C ****************************************** D'INITIALISATION

c     ========================================
c      CALCUL DU FLUX F ET DE L'ENERGIE DES ELECTRONS PRECIPITES
c
c        FE : FLUX D'ENERGIE PRECIPITES
c     ========================================

	icomp=9
	do 9000 i=1,25
	icomp=icomp+1
		do 9010 j=1,80
		FEtemp((i-1)*80+j)=FEtem2((icomp-1)*80+j)
		Etemp((i-1)*80+j)=Etem2((icomp-1)*80+j)
 9010		continue
 9000	continue

      do 1540 i=1,dndl

      sptemp(i)=5.2*(FEtemp(i)**0.5)

cesr
      
      sptemp(i)=ideaele*sptemp(i)
      
cesr
      
      shtemp(i)=0.55*sptemp(i)*((3.*Etemp(i)/2.)**0.6)
            
cesr
      
      shtemp(i)=ideaele*shtemp(i)
      
cesr      

      if(iopsp.eq.1) goto 1550
c     
c     INTERPOLATION QUADRATIQUE DE LA LOI POUR SP ET SH
c     
      sptemp(i)=sqrt(sp0temp(i)**2+sptemp(i)**2)
      shtemp(i)=sqrt(sh0temp(i)**2+shtemp(i)**2)

      goto 1540

 1550    continue
c     
c     INTERPOLATION LINEAIRE
c     
      sptemp(i)=sptemp(i)+sp0temp(i)/1.
      shtemp(i)=shtemp(i)+sh0temp(i)/1.

 1540     continue

c!!!	write(*,*)'COUCOU JE SUIS EN ELLCONI'
csss	write(*,*)(Ti(i),i=1,ndl)

c!!!	
c!!!	On lit les conductivites analytiques
c!!!

c!!!	open(4,file='Pedersen.lis',status='unknown')
c!!!	read(4,*)(Barped(i),i=1,49*80)
c!!!	close(4)

c!!!	open(4,file='Hall.lis',status='unknown')
c!!!	read(4,*)(Barhal(i),i=1,49*80)
c!!!	close(4)

c!!!	icompt=9

c!!!	do 8500 i=1,25

c!!!	icompt=icompt+1

c!!!		do 8510 j=1,80

c!!!		sptemp((i-1)*80+j)=Barped((icompt-1)*80+j)
c!!!		shtemp((i-1)*80+j)=Barhal((icompt-1)*80+j)

c!!!8510		continue

c!!!8500	continue

c!!!	open(4,file='Ped3.lis',status='unknown')
c!!!	write(4,*)(sptemp(i),i=1,dndl)
c!!!	close(4)

c!!!	open(4,file='Hal3.lis',status='unknown')
c!!!	write(4,*)(shtemp(i),i=1,dndl)
c!!!	close(4)

c!!!
c!!!	On lit les courants alignes analytiques
c!!!

c!!!	open(4,file='AligneAnal.lis',status='unknown')
c!!!	read(4,*)(alite2(i),i=1,49*80)
c!!!	close(4)



c!!!
c!!!	On met ds Rnitemp unepartie de alite2 afin de calculer
c!!!	le potentiel directement a partir des J// et non de
c!!!	la moyenne de pression
c!!!

	icompt=9

	do 8600 i=1,25

	icompt=icompt+1

		do 8610 j=1,80

	Rnitemp((i-1)*80+j)=alite2((icompt-1)*80+j)
		
cesr

	Rnitemp((i-1)*80+j)=ideaion*Rnitemp((i-1)*80+j)

cesr		
		
8610		continue

8600	continue

	open(4,file='AligneModele.lis',status='unknown')
	write(4,*)(Rnitemp(i),i=1,25*80)
	close(4)

       call ellconi(sptemp,shtemp,al1,al2,al3,al4,phiopt,utemp,
     %dndl,dnt,
     %Ttemp,Atemp,Btemp,m1temp,m2temp,muns,ndimat,a11,a12,a21,a22,a1,
     %ff,alp1,alp2,ao,Titemp,Rnitemp,Rnn,TT,aa
     %,DPHI,rtmax,rnmax,dn1,dn2,utem2)

c	do 8000 i=1,49
c		do 8010 j=1,80
c		utem2((i-1)*80+j)=0.
c8010		continue
c8000	continue

c!!!	open(4,file='peypot.2kev',status='Unknown')
c!!!	read(4,*)(utem2(i),i=1,49*80)
c!!!	close(4)


	icomp=9
	do 7000 i=1,25
	icomp=icomp+1
		do 7010 j=1,80
		utem2((icomp-1)*80+j)=utemp((i-1)*80+j)
 7010		continue
 7000	continue

c!!!	write(*,*)'utem2',(utem2(9*80+j),j=1,80)

	return
	end

c==================================================================
c==================================================================

       subroutine detmns(nndl,nnt,m,mu,nndlt)
       dimension mu(1),m(3,nnt)
       n1=nndl+1
       do 1 i=1,n1
 1     mu(i)=0
       do 2 i=1,nnt
       nmin=min0(m(1,i),m(2,i),m(3,i))
       do 3 j=1,nndlt
       k=m(j,i)
       k1=k+1    
       k2=(k-nmin)*2+1
       mu(k1)=max0(mu(k1),k2)
 3     continue
 2     continue
       do 4 i=1,n1
       mu(i+1)=mu(i+1)+mu(i)
 4     continue
       return
       end 

c===================================================================
c===================================================================
      
      subroutine tencon1(spjour,rapp,iopec,sp,sh,ioptint,rapp1            
     %  ,nndl,phi,T,nn1)
c
c        CE PROGRAMME CALCULE LE TENSEUR DE CONDUCTIVITE
c
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      dimension sh(1),sp(1)
      dimension phi(1),T(1)
c     ========================================
c     ENTREE DES VALEURS DE LA MATRICE 2x2
c     ========================================
      pi=4.*atan(1.)
      hB=(2.*pi)/nn1
      spnuit=spjour/rapp1
      shnuit=rapp*spnuit
      shjour=rapp*spjour
      pi=4.*atan(1.)
      eps1=hB/5.
      pii=pi/2.
      piii=3*pi/2.
      do 3 i=1,nndl
      ph=phi(i)
      ph1=abs(ph-pii)
      ph2=abs(ph-piii)
      if(ioptint.eq.0) goto1
      i1=i+1
      if(ioptint.eq.2) i1=i+2
      if((ph1.lt.(hB+eps1)).or.(ph2.lt.(hB+eps1)))goto 30
1     if((ph.lt.(piii-eps1)).and.(ph.gt.(pii+eps1))) goto 2
c     ========================================
c     le noeud i est dans l'obscurite,et si iopt >ou= 1 ph app
c     [pi/2-hB,3pi/2+hB]
c     ========================================
      sp(i)=spnuit
      sh(i)=shnuit
      goto 4
2     continue
c     ========================================
c     le noeud i est eclaire
c     ========================================
c attention: cos(x)**(1/2)

      sp(i)=spjour*sqrt(abs(-sin(T(i))*cos(phi(i))))
      sh(i)=shjour*sqrt(abs(-sin(T(i))*cos(phi(i))))
      goto 4
30    continue
c     ========================================
c     iopt est> ou = 1,et (^phi-pi/2^<=hB+eps ou ^phi-3pi/2^<=hB+eps
c     ========================================
      if((ph1.lt.eps1).or.(ph2.lt.eps1)) goto32
      if(ioptint.eq.1) goto 33
      if(((ph-pii).lt.0).or.((ph-piii).gt.0)) goto 34
c     ========================================
c     ici ph-pi/2=hB=>ph=PI/2+Hb
c     ou  ph-2pi/3=-hB=>ph=3pi/2-hB
c     ========================================
c attention: cos(x)**(1/2)

      i1=i+1
      if((ph-piii+hB).gt.-eps1) i1=i+3
      sp(i)=spnuit+(spjour*sqrt(abs(sin(T(i))
     %*abs(cos(phi(i1)))))-spnuit)*9./16.
      sh(i)=sp(i)*rapp
      goto 4
34    continue
c     ========================================
c     ici ph-pi/2=-hB=>ph=pi/2-hB
c     ou  ph-3pi/2=hB=>ph=3pi/2+hB
c     ========================================
c attention: cos(x)**(1/2)

      i1=i+3
      if((ph-piii).gt.0) i1=i+1
      sp(i)=spnuit+(spjour*sqrt(abs(sin(T(i))
     %*abs(cos(phi(i1)))))-spnuit)*1./16.
      sh(i)=sp(i)*rapp
      goto 4
33    continue
c     ========================================
c     ioptint=1
c     ========================================
c attention: cos(x)**(1/2)

      if(((ph-pii).lt.0).or.((ph-piii).gt.0))goto 35
      sp(i)=spjour*sqrt(abs(sin(T(i))*abs(cos(phi(i)))))
      sh(i)=sp(i)*rapp
      goto4
35    sp(i)=spnuit
      sh(i)=shnuit
      goto4
32    sp(i)=spnuit+(spjour*sqrt(abs(+sin(T(i))
     %*abs(cos(phi(i1)))))-spnuit)*1/4.
      sh(i)=sp(i)*rapp
4     continue
       amag=1/sqrt(1.+3*(cos(T(i))**2))
      sh(i)=sh(i)*amag
      sp(i)=sp(i)*amag
3     continue
csss      print 117,spnuit,shnuit,spjour,shjour
117   format ('MODELE DE CONDUCTIVITE:VALEURS OPTIMALES',/,
     &  'LA NUIT:Sp=',e12.5,'et Sh=',e12.5,/,
     &  'LE JOUR:Sp=',e12.5,'et Sh=',e12.5)
      if(iopec.lt.1) return
      print 118,(sp(i),i=1,nndl),(sh(i),i=1,nndl)
118   format ('VERIFICATION ',/,6(e12.5))
      return
      end
      
c 	=========================================================
c 	=========================================================

      subroutine ellconi(sp,sh,al1,al2,al3,al4,phiopt,u,nndl,nnt,
     &  T,A,phi,m1,m2,muns,ndimat,a11,a12,a21,a22,a1,
     &	f,alp1,alp2,ao,Ti,Rni,Rnn,TT,aa
     &  ,DPHI,rtmax,rnmax,nn1,nn2,utem2)
c     ========================================
c     PROGRAMME DE CONSTRUCTION DES MATRICES DE RESOLUTION
c     DU PROBLEME ELLIPTIQUE:APPROX P1
c     LE SYSTEME DE COORDONNEES UTILISE SE DEDUIT DES COORDONNEES
c     SHERIQUES PAR: B=phi=longitude ou le temps local
c                    A=sin(teta)**2=sinus carre de la colatitude
c     ========================================
	real*4 Nmax
      common/donnees/e,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,hB,
     %n1,n2,Nmax
      common/triang2/almax,almin,rlpas,hpas1,hpas2
      dimension a11(1),a12(1),a21(1),a22(1)
      dimension sh(1),sp(1)
      dimension al1(1),al2(1),al3(1),al4(1),u(1)
      dimension alp1(1),alp2(1),no(3),ao(1)
      dimension a1(ndimat)
      dimension at1(9),at(9)
      dimension f(1)
      dimension T(1),A(1),phi(1),m1(3,nnt),m2(1),muns(1)
      dimension Ti(1),TT(1),DPHI(1),Rni(1),Rnn(1),utem2(1)
      real*4 aa(6),x1,x2,x3,y1,y2,y3,Heel(49*80)
      
c ****	INITIALISATION DU SECOND MENBRE DE LA MATRICE ***

	do 900 i=1,nndl
	f(i)=0.

	u(i)=0

 900	continue

      do 16 i=1,nndl

c      sinI=2*cos(T(i))/sqrt(4-3*(sin(T(i))**2))

      sinI=1.

      a11(i)=sp(i)
      a12(i)=-sh(i)/sinI      
      a22(i)=sp(i)/(sinI**2)
      a21(i)=sh(i)/sinI
16    continue
	
c     ========================================
c     ENTREE DE LA FONCTION DE BORD
c     ========================================
c
c           S : CONDUCTANCEDE LA CEINTURE EQUATORIALE
c           r : RAYON TERRESTRE
c
c attention: cos(x)**(1/2)

      S=3.e+08
c     r=6.36e+06
csss      print 119,S,r
      pi2=2.*atan(1.)
            pi3=3*pi2
      eps=1.e-04
      do 4 i=1,nndl
       ao(i)=S/r      
       if(((phi(i)-pi3).gt.eps).or.((phi(i)-pi2).lt.eps))  goto15
       if((phi(i).gt.pi2).and.(phi(i).lt.pi3)) goto 3
 15    ao(i)=ao(i)/30.
       goto 4
 3     ao(i)=ao(i)*sqrt(abs(-sin(T(i))*cos(phi(i))))
c	ao(i)=sin(80./90.*pi2)*97.5
4     continue
c     ========================================
c     ENTREE DES VALEURS DU JACOBIEN DU CHANGEMENT DE VARIABLE
c     ========================================
      do 5 i=1,nndl
      alp1(i)=1./sin(T(i))
      alp2(i)=sin(2*T(i))
5     continue
c     ========================================
c     BALAYAGE DES TRIANGLES ET ASSEMBLAGE DE LA MATRICE
c     ========================================
      munsdl=muns(nndl+1)
      do 1 i=1,munsdl
1     a1(i)=0.
      do 6 i=1,nnt
      do 7 j=1,3
      no(j)=m1(j,i)
7     continue
      x21=phi(no(2))-phi(no(1))
      x31=phi(no(3))-phi(no(1))
      y21=A(no(2))-A(no(1))
      y31=A(no(3))-A(no(1))
      if((i-2*nn1*(i/(2*nn1))).eq.(2*nn1-1)) x21=4*pi2-phi(no(1))
      if((i-2*nn1*(i/(2*nn1))).eq.0) x31=phi(no(3))-4*pi2
      delta=x21*y31-x31*y21
      call maelri(at1,x21,x31,y21,y31,delta,no,a11,a12,a21,a22,
     &  alp1,alp2,al1,al2,al3,al4)
      m=0
      j1=0
      j2=0
      do 8 j=1,3
      if(m2(no(j))-1)8,8,9
9     continue
      m=m+1
      j3=j2
      j2=j1
      j1=j
8     continue
      if(m.lt.2) goto10
      dlt=sqrt((phi(no(j1))-phi(no(j2)))**2+(A(no(j1))-A(no(j2)))**2)
      if((i-2*nn1*(i/(2*nn1))).eq.2*nn1-1) dlt=sqrt(x21**2+y21**2)
      if((i-2*nn1*(i/(2*nn1))).eq.0)
     % dlt=sqrt((x31-x21)**2+(y31-y21)**2)
       call maelbo(at,y21,y31,delta,no,al1,al2,al3,al4,ao,dlt,j1,j2)
      do 11 j=1,9
11     at1(j)=at1(j)+at(j)
10     continue
c	if (i.gt.2320) then
c	open(4,file='frontiere.lis',access='append')
c	write(4,*)'triangle :',i
c	write(4,*)(at(j),j=1,9)
c	close(4)
c			endif
      call assmns(at1,m1,muns,a1,i,nnt,3)
      
c ***	CALCUL DU SECOND MENBRE DE L'EQUATION ELLIPTIQUE ***

	x1=(phi(no(1)))
	x2=(phi(no(2)))
	x3=(phi(no(3)))
	y1=(A(no(1)))
	y2=(A(no(2)))
	y3=(A(no(3)))
	pi=4.*atan(1.0)
	nqi=i-(i/(2*nn1))*2*nn1
	if (nqi.eq.0) x1=nn1*(2*pi/nn1)
	if (nqi.eq.0) x2=nn1*(2*pi/nn1)
	if (nqi.eq.(2*nn1-1)) x2=nn1*(2*pi/nn1)
	
	do 902 k=1,3 
 	f(no(k))=delta/6.+f(no(k))
 902	continue
 6	continue

	do 903 j=1,nndl

	sinI=2*cos(T(j))/sqrt(4-3*(sin(T(j))**2))

c!!!	sinI=1.

	f(j)=f(j)*r*r*sin(T(j))*sinI
	f(j)=f(j)*(Rni(j))
	f(j)=f(j)/sin(2*T(j))

 903	continue
 			
c	open(4,file='sdm0.lis',access='sequential',position='append')
c	write(4,*)'f'
c	write(4,*)(f(i),i=1,nndl)
c	write(4,*)'Rni'
c	write(4,*)(Rni(i),i=1,nndl)
c	write(4,*)('*',i=1,70)
c	write(*,*) (muns(i),i=1,nndl+1)
c	write(*,*)('*',i=1,70)
c	write(*,*)(a1(i),i=2,ndimat)
c	close(4)

c
c             FACTORISATION DE LA MATRICE PAR UNE METHODE DE GAUSS
c
       eps1=1.e-08
       tgv=1.e+18
      do 21 i=1,nndl
       if(m2(i)-1) 21,12,21
12      m=muns(i+1)
       a1(m)=tgv
21     continue
c!!!	write(*,*) 'COUCOU JE SUIS EN CRMC3R'
      call crmc3r(muns,a1,nndl,eps1,0,1,a1,nretou)
c     ========================================
c     CALCUL DES DIFFERENTS SECOND MEMBRE
c       ET PRISE EN COMPTE DES CONDITIONS AUX LIMITES
c     ========================================

	j=0
	
      do 14 i=1,nndl

c	f(i)=0.
      if(m2(i)-1) 14,17,14

17    continue

      f(i)=tgv*phiopt*sin(phi(i))*1000.	

14    continue

c!!!	write(*,*)'limite',(f(i),i=1,80)

c!!!	write(*,*)'COUCOU JE SUIS EN DRCR3R'
c!!!	open(4,file='sdm1.lis',status='UNKNOWN')
c!!!	write(4,*)(f(i),i=1,nndl)
c!!!	close(4)

8100	continue

      call drcr3r(muns,a1,f,nndl,1,2,u)
      
c ***	ON REPASSE EN kV ***

      do 9876 i=1,nndl
      u(i)=u(i)/1000.
 9876 continue
 
      return
119   format ('CONDUCTANCE DE LA CEINTURE EQUATORIALE',e12.5,/,
     &  'RAYON DE LA TERRE',e12.5)
c107   format ('traitement du',i3,'ieme triangle')
c840   format ('VALEUR DES COEFFICIENTS',/,5(e12.5))
c110   format(1x,'nombre de noeuds sur la frontiere',i2)
c112   format('frontiere',4i3)
c111   format('longueur du bout de frontiere',e12.5)
c121   format ('SOLUTION',/,6(e12.5))      
      end

c=================================================================
c=================================================================

      subroutine maelri(at1,x21,x31,y21,y31,delta,no,a11,a12,a21,a22,
     &  alp1,alp2,al1,al2,al3,al4)
c     ========================================
c     LA MATRICE ELEMENTAIRE at1 EST RANGEE DANS L ORDRE
c            j\i    1  4  7
c                   2  5  8
c                   3  6  9
c     ===============================================
      dimension at1(9),no(3),a11(1),a12(1),a21(1),a22(1)
      dimension alp1(1),alp2(1),al1(9),al2(9),al3(9),al4(9)
      dimension f1(3),f2(3),f3(3),f4(3)
      do 1 i=1,9
      at1(i)=0.
1     continue
      do 2 i=1,3
      f1(i)=a11(no(i))*(alp1(no(i)))/alp2(no(i))
      f2(i)=a12(no(i))
      f3(i)=a21(no(i))
      f4(i)=a22(no(i))*alp2(no(i))/alp1(no(i))
2     continue
c
c       INTEGRAtION NUMERIQUE A TROIS POINTS
c
      s1=som1(f1)
      s2=som1(f2)
      s3=som1(f3)
      s4=som1(f4)
      do 3 i=1,9
      at1(i)=at1(i)+s1*(al1(i)*(y31**2)+al2(i)*(y21**2)-(al3(i)+
     %al4(i))*y21*y31)
      at1(i)=at1(i)+s2*(-y31*x31*al1(i)-y21*x21*al2(i)+y31*x21*al3(i)
     &  +y21*x31*al4(i))
      at1(i)=at1(i)+s3*(-x31*y31*al1(i)-y21*x21*al2(i)+x31*y21*al3(i)
     &  +y31*x21*al4(i))      
      at1(i)=at1(i)+s4*(al1(i)*(x31**2)+al2(i)*(x21**2)-(al3(i)
     &+al4(i))*x21*x31)
      at1(i)=at1(i)/delta
3     continue
      return
      end
      
c==================================================================
c==================================================================

      function som1(f)
      dimension f(3)
      s=0.
      do1 i=1,3
      s=s+f(i)
1     continue
      som1=s/6.
      return
      end
      
c=====================================================================
c=====================================================================

      subroutine maelbo(at,y21,y31,delta,no,al1,al2,al3,al4,ao
     &  ,dlt,j1,j2)
c
c           MATRICE ELEMENTAIRE DU BORD
c           INTEGRATION A DEUX POINTS
c
      dimension f(2)
      dimension no(3),al1(9),al2(9),al3(9),al4(9),ao(1),at(9)
      f(1)=ao(no(j1))
      f(2)=ao(no(j2))
      s=som2(f)
      do 3 i=1,9
      at(i)=al1(i)*(y31**2)+al2(i)*(y21**2)-(al3(i)+al4(i))*y31*y21
      at(i)=at(i)*s*dlt/(delta**2)
3     continue
      if(j1.eq.2) goto 1
      at(1)=0.
      at(2)=0.
      at(3)=0.
      at(4)=0.
      at(7)=0.
      return
1     continue
      at(3)=0.
      at(6)=0.
      at(7)=0.
      at(8)=0.
      at(9)=0.
      return
      end
      
c===================================================================
c===================================================================
       
      function som2(f)
      dimension f(2)
      som2=f(1)+f(2)
      som2=som2/2.
      return
      end
c=====================================================================
c=====================================================================

       subroutine assmns(at,m1,mu,a,i,nt,ndlt)
       dimension at(1),mu(1),m1(ndlt,nt),a(1)
       ka=0
       do 1 j=1,ndlt
       do 2 k=1,ndlt
              mj=m1(j,i)
              mk=m1(k,i)
       if (mj-mk)3,4,5
 3     continue
c     =========================================
c     mj<mk
c     =========================================
       l=mu(mk+1)-mk+mj
       goto6
 4     continue
c     ========================================
c     mj=mk
c     ========================================
       l=mu(mj+1)
       goto6
 5     continue
c     =======================================
c     mj>mk
c     =======================================
       l=mu(mj+1)-(mu(mj+1)-mu(mj))/2-mj+mk
 6     continue
       ka=ka+1
       a(l)=a(l)+at(ka)
 2     continue  
 1     continue
       return 
       end
       
c===================================================================
c===================================================================

       SUBROUTINE CRMC3R(MUDL,AO,NTDL,EPS,NENTRE,IMPRE,A,NRETOU)
c....................................................................
       DOUBLE PRECISION S,S1,SAU,SA
       DIMENSION MUDL(1),AO(1),A(1)
 15    FORMAT(' RESULTATS A VERIFIER: LE',I7,G15.7,G15.7,G15.7)
c
       NRETOU=0
       MUDLI=0
            DO 1 I=1,NTDL
	    I1=I-1
c
c      TRAITEMENT DE LA I-EME LIGNE DE L ET COLONNE DE U  (A=L*U)
c
       MUDLI1=MUDL(I+1)
       IH=(MUDLI1-MUDLI)/2
c      IH=NOMBRE DE COEFFICIENTS DE LA I-EME LIGNE (DIAGONALE EXCLUE)
       SA=0.D0
       IF (IH) 2,2,3
 3     IMI=I-IH
c      IMI=NO DE 1-ERE COLONNE NON NULLE DE LA I-EME LIGNE
       IAO=MUDLI
       IA1=MUDLI+IH
       MUDLJ=MUDL(IMI)
       IMI1=IMI-1
c           BOUCLE SUR LES COLONNES DE LA I-EME LIGNE DE L
c           BOUCLE SUR LES LIGNES DE LA I-EME COLONNE DE U
            DO 4 JJ=1,IH
            J=IMI1+JJ
            J1=J-1
            MUDLJ1=MUDL(J+1)
            JH=(MUDLJ1-MUDLJ)/2
            JMI=J-JH
            JAO=MUDLJ
            JA1=MUDLJ+JH
c           JAO POINTE EN TETE DE LA J-EME LIGNE
c           JA1 POINTE EN TETE DE LA J-EME COLONNE
            IS=IMI-JMI
            IF (IS)5,5,6
c
 5          MA=JMI
            IA=IAO-IS
            JA=JAO
            IAA=IA1-IS
            JAA=JA1
            GOTO 7
c
 6          MA=IMI
            IA=IAO
            JA=JAO+IS
            IAA=IA1
            JAA=JA1+IS
c
 7          S=0.D0
            S1=0.D0
c
            K1=J-MA
            IF (K1) 9,9,8
c
 8                  DO 10 K=1,K1
c                   L(I,J)
                    SAU=A(IA+K)
                    S=S+SAU*A(JAA+K)
c                   U(J,I)
                    SAU=A(JA+K)
                    S1=S1+SAU*A(IAA+K)
 10                 CONTINUE
c
 9          IA=IA+K1+1
            IAA=IAA+K1+1
            A(IA)=(AO(IA)-S)/A(MUDLJ1)
            A(IAA)=AO(IAA)-S1
            MUDLJ=MUDLJ1
 4          CONTINUE
c
c     TRAITEMENT DU I-EME COEFFICIENT DIAGONA
c

            DO 11 K=1,IH
            SAU=A(IAO+K)
            SA=SA+SAU*A(IA1+K)
 11         CONTINUE
c
 2    AA=AO(MUDLI1)
      X=AA-SA
c     TEST SUR LA PRECISION DU PIVOT
      IF (ABS(X)-EPS*ABS(AA)) 13,13,14
c     AU DESSOUS DE LA PRECISION
 13   PRINT 15,I,EPS,X,AA
      NRETOU=1
      IF (NENTRE) 14,16,14
c
 14   A(MUDLI1)=X
      MUDLI=MUDLI1
 1    CONTINUE
c
 16   NL=NTDL
      IF (IABS(IMPRE)-6) 17,18,19
 18   NL=10
 19   NL=10
c
 17   RETURN
      END
                    
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                        SP DRCR3R
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  BUT : DESCENTE ET/OU REMONTEE D'UN SYSTEME DONT LA MATRICE EST
c  ---   FACTORISEE SOUS LA FORME A = L * U DITE DE GAUSS
c        L EST UNE MATRICE TRIANGULAIRE INFERIEURE A DIAGONALE UNITE
c        U EST UNE MATRICE TRIANGULAIRE SUPERIEURE
c        L ET U SONT RANGEES DANS A SOUS FORME PROFIL
c        A (I,J) NON NUL => A (J,I) NON NUL (PROFIL SYMETRIQUE)
c        A NON SYMETRIQUE MUDL POINTE SUR LES COEFFICIENTS DIAGONAUX
c        LE SP CRMC3R A DU ETRE EXECUTE AUPARAVANT
c        VERSION REELLE SIMPLE PRECISION
c
c  PARAMETRES D'ENTREE:
c  --------------------
c  MUDL  : MUDL(1)=0,MUDL(I+1)= ADRESSE DU I-EME COEF DIAGONAL
c  A     : MATRICE FACTORISEE A = L * U
c  BO    : LES NDSM SECONDS MEMBRES BO(NDSM,NDTC)
c  NDTL  : ORDRE DE LA MATRICE A
c  NDSM  : NBRE DE SECONDS MEMBRES
c  NIVEAU: 0 U*B=BO        REMONTEES SEULES EXACTES SEULEMENT SI B =
c          BO
c          1 L*B=BO        DESCENTES SEULES 
c          2 L*U*B=BO      DESCENTES,REMONTEES
c
c  PARAMETRES DE SORTIE:
c  ---------------------
c  B     : TABLEAU DES RESULTATS B (NDSM,NTDL)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  PROGRAMMEUR : A.PERRONNET LAN 189 ET IRIA PARIS OCTOBRE 77
c....................................................................
     
c===================================================================
c===================================================================

       SUBROUTINE DRCR3R(MUDL,A,BO,NTDL,NDSM,NIVEAU,B)
c.......................................................................
c
       DOUBLE PRECISION S,S1,SAU,SA
       DIMENSION MUDL(1),A(1),BO(1),B(1)
 2     FORMAT (' ERREUR: DANS DRCR3R CAR A NON SYMETRIQUE.
     %NIVEAU',I12)
c
c      COHERENCE DU PARAMETRE NIVEAU
c      -----------------------------
c
       IF(NIVEAU.GE.0  .AND.  NIVEAU.LE.2)  GOTO 1
       PRINT 2,NIVEAU
       STOP
c      LES DESCENTES
c      -------------
c
 1     IF(NIVEAU.EQ.0) GOTO 10
       IA=0
       MUDLI=0
          DO  500 I=1,NTDL
	  I1=I-1
	  MUDLI1=MUDL(I+1)
	  IH=(MUDLI1-MUDLI)/2
	  IMI=I-IH
	  IAO=IA
	  IBO=(IMI-1)*NDSM
c	  
             DO 3 NC=1,NDSM
	     SA=0.D0
             IA=IAO
	     IB=IBO + NC
	     IF (IH)4,4,5
c
 5              DO 6 K=1,IH
	        SAU=A(IA + K)
	        SA=SA+SAU*B(IB)
	        IB=IB+NDSM
 6              CONTINUE
 4           B(IB)=BO(IB)-SA
 3           CONTINUE
          IA=IA+MUDLI1-MUDLI
	  MUDLI=MUDLI1
 500      CONTINUE
          IF (NIVEAU.EQ.1) GOTO 100
c	  
c         LES REMONTEES
c         -------------
c
 10       NTDL1=NTDL + 1
          IB = NTDL1
	  MUDLI1=MUDL(NTDL1)
	  IA=MUDLI1 + 1
	     DO 11 II=1,NTDL
	     I=NTDL1-II
	     MUDLI=MUDL(I)
	     IH=(MUDLI1-MUDLI)/2
	     IA=IA-1
	     IB=IB-1
	     IBN=IB*NDSM
	     S=A(IA)
c
                 DO 12 NC=1,NDSM
		 B(IBN)=B(IBN)/S
		 IBN=IBN-1
 12              CONTINUE
c
             IF (IH) 13,13,14
 14          IBO=IBN+NDSM
c
                 DO 15 K=1,IH
		 IA=IA-1
		 IBB=IBO
		 S=A(IA)
c
                    DO 16 NC=1,NDSM
		    B(IBN)=B(IBN)-B(IBB)*S
		    IBB=IBB-1
		    IBN=IBN-1
 16                 CONTINUE
c
 15              CONTINUE
c
             IA=IA-IH
 13          MUDLI1=MUDLI
 11          CONTINUE
c
 100   RETURN
       END
