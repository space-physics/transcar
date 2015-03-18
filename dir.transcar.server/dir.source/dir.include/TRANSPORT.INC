!
! Computational parameters :
! --------------------------
! nbralt   =	IPM:(MMAX) length of altitude grid (number of layers +1)
! nbrsp    =	IPJ: number of species (1=N2 , 2=O2 , 3=O , 4=H , 5=He)
! nbrango2 =	IPK : (KSTR2) one half of the maximum number of streams
! nbrang   =	IPK : (KSTR2) maximum number of streams
! nbrexc   =	IPJS: maximum number of excitation states
! nbrionst =	IPJSP: number of ionized states
! nbren    =	IPN: maximum length of the energy grid
!
!
         integer nbralt,nbrsp,nbrango2,nbrang,nbrexc,nbrionst,nbren
         PARAMETER (nbralt=201,nbrsp=5,nbrango2 = 16,nbrang=32,
     .               nbrexc=15,nbrionst=5,nbren=400)
!
!       Nom de la planete d'application, et distance au soleil (UA)
        character*10 planete
        real rayonUA
!
! 	For each specie:
! 	  . Name.
! 	  . Atomic mass (number).
!
 	real atomas(8)
!
	character*5 specie(8)
!
 	
 	integer unfic_in_transcar,unfic_out_transcar
 	integer fic_datdeg,fic_degout,fic_transout,fic_dattrans
 	integer fic_crsin_degrad
 	integer unfic_crsout_degrad,unfic_rdtout_degrad
 	integer icrsin,irdtin,icrsphot
 	integer ineutr,ielec,iindat,transcar_dat
        integer ifeldat,ifelprt,ifeltrans
        integer fic_euvcoeff,fic_euvindex
!
        include '../dir.include/CHEMIN.INC'


        data planete/'Terre'/
        data rayonUA /1./
 	
 	data atomas /28.0134,31.9988,15.9994,1.00797,4.0026,14.0067,
     .		     39.948,30.0061/
	data specie/'[N2]','[O2]','[O]','[H]','[He]','[N]','[A]','[NO]'/
 	
 	data ineutr/10/,ielec/12/,iindat /14/,transcar_dat/15/
 	data unfic_in_transcar/16/,unfic_out_transcar/18/
 	data fic_datdeg/20/,fic_degout/22/,fic_transout/24/
 	data fic_crsin_degrad/26/,fic_dattrans/28/
 	data unfic_crsout_degrad/30/,unfic_rdtout_degrad/32/
 	data icrsin/34/,irdtin/36/,icrsphot/38/
        data ifeldat/40/,ifelprt/50/,ifeltrans/52/
        data fic_euvcoeff/54/,fic_euvindex/56/
