!
! Computational parameters :
! --------------------------
! nbralt   =    IPM:(MMAX) length of altitude grid (number of layers +1)
! nbrsp    =    IPJ: number of species (1=N2 , 2=O2 , 3=O , 4=H , 5=He)
! nbrango2 =    IPK : (KSTR2) one half of the maximum number of streams
! nbrang   =    IPK : (KSTR2) maximum number of streams
! nbrexc   =    IPJS: maximum number of excitation states
! nbrionst =    IPJSP: number of ionized states
! nbren    =    IPN: maximum length of the energy grid
!
!
      integer,parameter :: nbralt=201,nbrsp=5,nbrango2 = 16,nbrang=32,     
     &                    nbrexc=15,nbrionst=5,nbren=400
!
!       Nom de la planete d'application, et distance au soleil (UA)

!
!     For each specie:
!       . Name.
!       . Atomic mass (number).
!
      real,parameter :: atomas(8) = [28.0134,31.9988,15.9994,1.00797,
     &   4.0026,14.0067, 39.948,30.0061]       

        include 'CHEMIN.INC'

      real, parameter :: rayonUA=1.,N_0=1.e4,T_0=1000.
     
       character(*), parameter :: specie(*)=['[N2]','[O2]','[O] ',
     &  '[H] ', '[He]','[N] ','[A] ','[NO]']
     
      integer,parameter :: ineutr=10,ielec=12,iindat=14,transcar_dat=15,
     &  unfic_in_transcar=16,unfic_out_transcar=18,
     &  fic_datdeg=20,fic_degout=22,fic_transout=24,
     &  fic_crsin_degrad=26,fic_dattrans=28,
     &  unfic_crsout_degrad=30,unfic_rdtout_degrad=32,
     &  icrsin=34,irdtin=36,icrsphot=38,
     &  ifeldat=40,ifelprt=50,ifeltrans=52,
     & fic_euvcoeff=54,fic_euvindex=56
