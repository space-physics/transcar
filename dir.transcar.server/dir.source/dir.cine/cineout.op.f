c
 	subroutine cineout(nalt,Chaufelec,denelc,prodiontot,
     .          Ne_supra,courant_supra,Te_supra,Chaleur_supra,
     .		Ne,npt,indlim,nx,zlim,zlim_1,z,Heat,Ph,Po,Po2,Pn2,Pn,
     .          Ne_sup,courant_sup,Te_sup,Chaleur_sup)
c
        implicit none
c
        integer npt
c
        include 'TRANSPORT.INC'
c
c 	INPUTS
 	real Chaufelec(nbralt),denelc(nbralt),prodiontot(nbrsp*2,nbralt)
        real Ne_supra(nbralt),courant_supra(nbralt),Te_supra(nbralt),
     .          Chaleur_supra(nbralt)
        real z(npt),Ne(npt)
 	integer nalt,nx
 	integer indlim
        real zlim,zlim_1,coef
c
c 	INPUTS/OUTPUTS
        real Ph(npt),Po(npt),Po2(npt),Pn2(npt),Pn(npt)
c
c 	OUTPUTS
        real Heat(npt)
        real Ne_sup(npt),courant_sup(npt),Te_sup(npt),Chaleur_sup(npt)
c
c 	INTERNAL
 	integer ialt
	real xo_1,xo,xo2_1,xo2,xn2_1,xn2,xn_1,xn,xh_1,xh
 	real xqe_1,xqe
        real hauteur,x1,x2,y1,y2
        real HQe_1,HPo_1,HPo2_1,HPn2_1,HPn_1,HPh_1
 	real dz
	real xnes,xnes_1,xjes,xjes_1,xtes,xtes_1,xqes,xqes_1
	real Hnes_1,Hjes_1,Htes_1,Hqes_1
        integer indo2
c
c ====================================================================
c       definition de la fonction de hauteur d echelle
c	definition of the function of scale height
c ====================================================================
c
        hauteur(y1,y2,x1,x2)=(x1-x2)/(log(y2)-log(y1))
c
c ====================================================================
c
c
	do ialt=1,indlim
c	  Heat (indlim+1-ialt)=Chaufelec(ialt)/denelc(ialt)*1.6e-12
          Heat (indlim+1-ialt)=Chaufelec(ialt)*1.6e-12
	  Pn2(indlim+1-ialt)= prodiontot(1,ialt)
	  Po2(indlim+1-ialt)= prodiontot(2,ialt)
	  Po (indlim+1-ialt)= prodiontot(3,ialt)
 	  Ph (indlim+1-ialt)=prodiontot(4,ialt)
	  Pn (indlim+1-ialt)=prodiontot(6,ialt)
	  Ne_sup (indlim+1-ialt)=Ne_supra(ialt)
	  courant_sup (indlim+1-ialt)=courant_supra(ialt)
	  Te_sup (indlim+1-ialt)=Te_supra(ialt)
	  Chaleur_sup (indlim+1-ialt)=1.6e-12*Chaleur_supra(ialt)
	enddo
c
c ---- PIERRE LOUIS : REGARDE SI CELA TE SUFFIT. SI TU MODIFIES CETTE
C      ROUTINE, RENVOIE LA MOI MODIFIEE.
c	IF IT LOOKS JUST TE. IF YOU CHANGED THIS ROUTINE RETURNS THE AMENDMENTS TO ME.

        xo_1 =Po(indlim-1)
        xo   =Po(indlim)
        xo2_1=Po2(indlim-1)
        xo2  =Po2(indlim)
        xn2_1=Pn2(indlim-1)
        xn2  =Pn2(indlim)
        xn_1 =Pn(indlim-1)
        xn   =Pn(indlim)
	xqe_1=Chaufelec (2)
	xqe  =Chaufelec (1)
	xnes_1=Ne_supra (2)
	xnes  =Ne_supra (1)
	xjes_1=Courant_supra (2)
	xjes  =Courant_supra (1)
	xtes_1=Te_supra (2)
	xtes  =Te_supra (1)
	xqes_1=Chaleur_supra (2)
	xqes  =Chaleur_supra (1)
c
	if (xqe_1*xqe.ne.0.) then
	  HQe_1=1./abs(hauteur(xqe_1,xqe,zlim_1,zlim))
	else
	  HQe_1=0.
	endif
	if (xo_1*xo.ne.0.) then
	  HPo_1=1./abs(hauteur(xo_1 ,xo ,zlim_1,zlim))
	else
	  HPo_1=0.
	endif
	if (xo2_1*xo2.ne.0.) then
	  HPo2_1=1./abs(hauteur(xo2_1,xo2,zlim_1,zlim))
	else
	  HPo2_1=0.
	endif
	if (xn2_1*xn2.ne.0.) then
	  HPn2_1=1./abs(hauteur(xn2_1,xn2,zlim_1,zlim))
	else
	  HPn2_1=0.
	endif
	if (xn_1*xn.ne.0.) then
	  HPn_1=1./abs(hauteur(xn_1,xn,zlim_1,zlim))
	else
	  HPn_1=0.
	endif
	if (xh_1*xh.ne.0.) then
c	  HPh_1=1./abs(hauteur(xh_1 ,xh ,zlim_1,zlim))
	  HPh_1=0.
	else
	  HPh_1=0.
	endif
	if (xnes_1*xnes.ne.0.) then
	  Hnes_1=1./abs(hauteur(xnes_1,xnes,zlim_1,zlim))
	else
	  Hnes_1=0.
	endif
	if (xjes_1*xjes.ne.0.) then
	  Hjes_1=1./abs(hauteur(xjes_1,xjes,zlim_1,zlim))
	else
	  Hjes_1=0.
	endif
	if (xtes_1*xtes.ne.0.) then
	  Htes_1=1./abs(hauteur(xtes_1,xtes,zlim_1,zlim))
	else
	  Htes_1=0.
	endif
	if (xqes_1*xqes.ne.0.) then
	  Hqes_1=1./abs(hauteur(xqes_1,xqes,zlim_1,zlim))
	else
	  Hqes_1=0.
	endif

	do ialt=indlim+1,nx
	  dz=zlim-z(ialt)
c	  Heat (ialt)=Heat (indlim)*exp(dz*HQe_1)*Ne(indlim)/Ne(ialt)
          Heat (ialt)=Heat (indlim)*exp(dz*HQe_1)
          Po (ialt)=Po (indlim)*exp(dz*HPo_1)
          Po2(ialt)=Po2(indlim)*exp(dz*HPo2_1)
          Pn2(ialt)=Pn2(indlim)*exp(dz*HPn2_1)
          Ph (ialt)=Ph (indlim)*exp(dz*HPh_1)
          Pn (ialt)=Pn (indlim)*exp(dz*HPn_1)
          coef=((6378+zlim)/(6378+z(ialt)))**3
          Ne_sup(ialt)=Ne_sup(indlim)*coef
          Courant_sup(ialt)=Courant_sup(indlim)*coef
          Te_sup(ialt)=Te_sup(indlim)
          Chaleur_sup(ialt)=Chaleur_sup(indlim)*coef
	enddo

c
	return
 	end
