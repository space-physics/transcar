	subroutine Drive11(iecompteur,continu)
c access write(*,*) tencon1 ellconi ivol1 ivol5 Rni1111 statio
c**********************************************************
c*	PROGRAMME TESTANT LA METHODE DE COUPLAGE ENTRE LE *
c*	TIEGCM ET LE MIM				  *
c**********************************************************

CC	PROGRAMME PRINCIPAL

CC	LE PROGRAMME FAIT APPEL A LA ROUTINE ELLCONI QUI
CC	CALCULE LE POTENTIEL. ENSUITE IL FAIT APPEL A LA
CC	ROUTINE MARCHE QUI FAIT AVANCER LES ELECTRONS TOUT
CC	D'ABORD ET LES IONS ENSUITE.
CC	LA ROUTINE ELLCONI CALCULE LE POTENTIEL utemp DANS
CC	LA GRILLE DE CICELEY. CE POTENTIEL EST ENSUITE INTERPOLE
CC	SUR LA GRILLE A PARTIR DE LAQUELLE EST CALCULEE  LE MOU
CC	VEMENT DES PARTICULES.
CC	LA ROUTINE MARCHE CALCULE LE NOMBRE DE PARTICULES CONTE
CC	NUES DANS UN TUBE DE FLUX MAGNETIQUE UNITE ane1 POUR LES
CC	ELECTRONS ET Rni POUR LES IONS. ELLE CALCULE AUSSI LES
CC	INVARIANTS ADIABATIQUES T V**(2/3) ener1 POUR LES ELEC
CC	TRONS ET Ti POUR LES IONS.
CC	A PARTIR DE ane1, ener1, Rni et Ti SONT CALCULES LES
CC	DENSITES ELECTRONIQUES ET IONIQUES, LES TEMPERATURES
CC	ELECTRONIQUES ET IONIQUES, LES FLUX DE PRECIPITATION
CC	DES ELECTRONS ET LES COURANTS ALIGNES. CES PARAMETRES
CC	SONT INTERPOLES SUR LA GRILLE DE CICELEY ET TRANSMIS
CC	A LA SUBROUTINE ELLCONI QUI CALCULE ALORS LE NOUVEAU
CC	POTENTIEL utemp.
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC	ATTENTION: LE COMMENTAIRE QUI SUIT N'EST PAS TOUT FAIT C
CC	VRAI:                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC	"IL FAUT NOTER QUE LES COURANTS ALIGNES N'ONT PAS BESOIN
CC	D'ETRE TRANSMIS A LA ROUTINE ELLCONI QUI N'UTILISE EN
CC	FAIT QUE Ti ET Rni DANS LE SECOND MEMBRE DE L'EQUATION
CC	ELLIPTIQUE DU POTENTIEL POUR PRENDRE EN COMPTE L'EFFET
CC	DES COURANTS ALIGNES".
		integer*4 iecompteur	 	
		logical	continu	
      integer*4 par1,par2,par3,dpar1,dpar2,dpar3
      integer*4	dpar11,dpar22,dpar33
      integer*4 dn1,dn2,dnt,dndl
      parameter (par1=1240,par2=2400,par3=2481)
      parameter (dpar1=(25*80),dpar2=(24*160),dpar3=(25*160+1))
      parameter (dpar11=(49*80),dpar22=(48*160),dpar33=(49*160+1))
      real*4 Nmax

	character*70 mot2,mot3,mot5,mot8

       common/debut/ndl,nt,dn1,dn2,dndl,dnt,ie,istati,
     %iinti,ibrupt,fheure,Bfm,Bfs,rnmax,rtmax,
     %A(par1),B(par1),T(par1),m1(3,par2),m2(par1),
     %Atemp(dpar1),Btemp(dpar1),Ttemp(dpar1),m1temp(3,dpar2),
     %m2temp(dpar1),
     %tabTM(par1,4),tabMT(dpar1,4)

      common/donnees/el,q,em,ab0,deltat,r,b0,emax,Amax,Amin,hA,
     $hB,n1,n2,Nmax
      common/triang1/n21,n22,n23,h21,h22,h23,Aint1,Aint2,raison
      common/triang2/almax,almin,rlpas,hpas1,hpas2

      common/transm/utem2(dpar11),FEtem3(dpar11),Etem3(dpar11),
     %alite2(dpar11)

      common/pot/FEtem2(dpar11),Etem2(dpar11)

      common/resul1/ane1(par1),ener1(par1),Ti(par1),Rni(par1),
     %ali(par1),u3(par1),uo(par1)

      common/debu1/phiopt,muns(dpar3),sp0temp(dpar1),
     %sh0temp(dpar1),ndimat,al1(9),al2(9),al3(9),al4(9),
     %Rnitemp(dpar1),Titemp(dpar1)

      common/conductances/sptemp(dpar1),shtemp(dpar1)
 
cesr
      
      common/dea/ideaele,ideaion,ideapot

cesr

cHAO97	
	common/fluxioni/FEion(dpar11),Fion(dpar11)
cHAO97

	real*4 transp(dpar11),transh(dpar11)


CC	INITIALISATION:

c	open(4,file='dir.source/dir.imm/dir.in_out/Neg.lis',status='unknown')
c	close(4)

c	open(4,file='dir.source/dir.imm/dir.in_out/Teg.lis',status='unknown')
c	close(4)
	ie=iecompteur	
	if (ie.eq.0)then

cc		open(4,file='dir.source/dir.imm/dir.in_out/donnees.dat',status='unknown')		
		open(4,file='dir.in_out/donnees.dat',status='unknown')
cc		On lit l'heure finale ou s'arrete le run
		read(4,*)fheure						

cc		On lit les temperatures et densites du plasma
cc		dans la source. On lit aussi la valeur de la moitie
cc		de la ddp a travers la calotte polaire.

cc		Pour les electrons, la densite rnmax et la temperature
cc		emax sont:

		read(4,*)rnmax						
		read(4,*)emax						


cc		Pour les ions, la densite est la meme que celle des
cc		electrons et leur temperature vaut:

		read(4,*)rtmax						

cc		La moitie de la ddp vaut:

		read(4,*)phiopt						

cc		Si ideaele=1, on prend en compte l'effet des precipitations
cc		sinon si ideaele=0 on ne le prend pas en compte:

		read(4,*)ideaele					

cc		Si ideaion=1, on prend en compte l'effet des courants alignes
cc		sinon si ideaion=0 on ne le prend pas en compte:

		read(4,*)ideaion					


cc		Si ideapot=1, on fait varier au cours du temps le potentiel
cc		sinon si ideapot =0 on ne le prend pas en compte

		read(4,*)ideapot					


cc		On lit la valeur du compteur d'ecriture ni. On ecrit 1 fois
cc		sur ni

		read(4,*)ni					     	


		close(4)						
							   	
		heure=0.

CC		LECTURE D'UN FICHIER INITIAL POUR L'INITIALISATION
CC		DES PARTICULES DANS MARCHE (POUR COUPLAGE AVEC LE
CC		TIEGCM SEULEMENT)
CC			-)	istati=0:	on ne lit pas de
CC						fichier initial
CC			-)	istati=1:	on lit le
CC						fichier initial

		istati=0

CC		PAS EN TEMPS EGAL A 240s POUR SIMULER CE QUE L'ON FERA AVEC
CC		LE TIEGCM:
CC			-)	UN PAS EN TEMPS DE 240s POUR LE TIEGCM
CC			-)	UN PAS EN TEMPS DE 120s POUR LA MAGNETO
CC			-) SPHERE.
CC		ON FERA DONC TOURNER LE MODELE DE MAGNETOSPHERE PENDANT 2
CC		ITERATIONS SUCCESSIVES DE 120s CHACUNE ALORS QUE LE TIEGCM
CC		TOURNERA LUI UNE SEULE ITERATION DE 240s.

		deltat=120.

		delta2=1.*deltat

CC		FIN DES INITIALISATIONS

CC		ITERATION EN TEMPS

cc		On ecrira tous les resultats dans le fichier de sortie
cc		tgcm.lis

cc		open(9,file='dir.source/dir.imm/dir.in_out/tgcm.lis',status='unknown')
		open(9,file='dir.in_out/tgcm.lis',status='unknown')

	end if
 
	write(*,*)'iteration: ',ie

CC	CALCUL DU POTENTIEL

	write(*,*)'POTENT'
	
cesr

	if (ie.ge.1)	then
	
	if (ideapot.eq.1)	then

	call POTENT
	
				endif
				
			else
			
	call POTENT
	
			
			endif

cesr

CC	DEPLACEMENT DES PARTICULES

	write(*,*)'MARCHE'

	call MARCHE

       heure=heure+delta2

CC	ECRITURE DES RESULTATS

CC	LE COMPTEUR ni PERMET D'ECRIRE TOUTES LES ni
CC	ITERATIONS.

CC	ni=5

	if ((ie-(ie/ni)*ni).eq.(0))	then
     
  	write(*,*) 'heure',heure,' ecriture sur fichier'

					endif

CC	ECRITURE DU POTENTIEL SANS COROTATION DANS LA GRILLE
CC	DU TIEGCM


	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*) (utem2(i),i=1,dpar11)
					endif
					
CC	ECRITURE DES FLUX D'ENERGIE ELECTRONIQUES
CC	INTERPOLES DANS LA GRILLE DU TIEGCM

	if ((ie-(ie/ni)*ni).eq.(0))	then
 	write(9,*) (FEtem2(i),i=1,dpar11)
					endif
					
CC	ECRITURE DE LA TEMPERATURE ELECTRONIQUE
CC	DANS LA GRILLE DU TIEGCM


	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*) (Etem2(i),i=1,dpar11)
					endif
					
c***	ECRITURE DES COURANTS ALIGNES DANS LA GRILLE DU TIEGCM***

	
	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*) (Alite2(i),i=1,dpar11)
					endif
					
CC	ECRITURE DU NOMBRE D'IONS
		
	
	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*) (Rni(i),i=1,ndl)
					endif
					
CC	ECRITURE DE LA TEMPERATURE DES IONS

	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*) (Ti(i),i=1,ndl)
					endif
					
CC	ECRITURE DE LA TEMPERATURE DES ELECTRONS

	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*) (ener1(i),i=1,ndl)
					endif
					
CC	ECRITURE DU NOMBRE D'ELECTRONS

	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*) (ane1(i),i=1,ndl)
					endif
					
CC	ECRITURE DE L'INVARIANT ADIABATIQUE DES ELECTRONS

	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*)(uo(i),i=1,ndl)
					endif
	
c***	ON CALCULE LES CONDUCTANCES ***

	icomp=9
	do 1000 i=1,25
	icomp=icomp+1
		do 1010 j=1,80
		transp((icomp-1)*80+j)=sptemp((i-1)*80+j)
		transh((icomp-1)*80+j)=shtemp((i-1)*80+j)
 1010		continue
 1000	continue

CC	ECRITURE DES CONDUCTIVITES DE PEDERSEN DANS LA GRILLE
CC	DU TIEGCM

	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*) (transp(i),i=1,dpar11)
					endif

CC	ECRITURE DES CONDUCTIVITES DE HALL DANS LA GRILLE
CC	DU TIEGCM

	
	if ((ie-(ie/ni)*ni).eq.(0))	then
	write(9,*) (transh(i),i=1,dpar11)
					endif

	
	if (ni.eq.1)			then

		if (ie.ge.1)			then

	close(9)

						endif

					endif


	if (fheure.lt.heure)	continu=.false.

c	open(4,file='compteur.lis',status='unknown')
c	write(4,*)ie
c	write(4,*)ni
c	close(4)

	iecompteur=ie
	return
	end subroutine Drive11
      	
