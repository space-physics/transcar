c
 	subroutine cineout(nalt,Chaufelec,denelc,prodiontot,
     .          Ne_supra,courant_supra,Te_supra,Chaleur_supra,
     .		Ne,npt,indlim,nb_alt,z,Heat,Ph,Po,Po2,Pn2,Pn,
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
 	integer nalt,nb_alt
 	integer indlim
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

	do ialt=indlim+1,nb_alt
	  Po (ialt) = Po (ialt-1)*Po (ialt-2)/Po (ialt-3)
	  Pn2(ialt) = Pn2(ialt-1)*Pn2(ialt-2)/Pn2(ialt-3)
	  Po2(ialt) = Po2(ialt-1)*Po2(ialt-2)/Po2(ialt-3)
	  Pn (ialt) = Pn (ialt-1)*Pn (ialt-2)/Pn (ialt-3)
	  Ph (ialt) = Ph (ialt-1)*Ph (ialt-2)/Ph (ialt-3)
	
	  Heat 	     (ialt) = Heat (ialt-1)
     &			     *Heat (ialt-2)/Heat (ialt-3)
	  Ne_sup     (ialt) = Ne_sup(ialt-1)
     &			     *Ne_sup(ialt-2)/Ne_sup(ialt-3)
	  Te_sup     (ialt) = Te_sup(ialt-1)
     &			     *Te_sup(ialt-2)/Te_sup(ialt-3)
	  courant_sup(ialt) = courant_sup(ialt-1)
     &			     *courant_sup(ialt-2)/courant_sup(ialt-3)
	  Chaleur_sup(ialt) = Chaleur_sup(ialt-1)
     &			     *Chaleur_sup(ialt-2)/Chaleur_sup(ialt-3)
	enddo
		
	return
 	end
