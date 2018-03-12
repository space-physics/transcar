      program transconvec_13
       
      use, intrinsic:: ieee_arithmetic, only: ieee_is_nan  
      include 'comm.f'
      include 'comm_sp.f'
        

!     Ce programme demarre le couple de programmes de transports.
!     Anciennement : eiscat.f
      character(80) split,gridfn

      integer,parameter :: npt=500,xcoeffno=1.,ncol0=50,intemps0=300,
     &                     nb_ion=6,
     &  nb_position_max=100     ! Modif DA 02/02 2001
        logical multi_position                  ! Modif DA 02/02 2001
        real longeo_position(nb_position_max)   ! Modif DA 02/02 2001
        real latgeo_position(nb_position_max)   ! Modif DA 02/02 2001
        character lecture_lat_lon*80            ! Modif DA 02/02 2001


!       indice des ions 
!       ---------------
!            (le premier indice faisant toujours refence a l espece)
!       
!       1 O+     2 H+   3 N+    4 N2+    5 NO+   6 O2+
!

        real dt_max
        common/CFL/dt_max



!     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(npt),       LN(npt),       AH (npt)
          Real     RLN(npt),      LH (npt),      RLH(npt)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(npt)
          Common  /FCT_GRID/ LO, LN, AH, RLN, LH, RLH, ROH, RNH, ADUGTH

          Real     SOURCE(npt)
          Common  /FCT_test/ SOURCE
        Integer Ipos1,Iposn,Iposnp,alpha,np,fid_NaN
        logical flag,flagatmos,flgconv,flgpot

        real*8 temps,tempsini,tempsconv,tempsconv_1,tempsdeb

        common /timeinfo/ temps
        
        real*8 tempsfin,tempslim,tempsint,tempsort,tempstube
        real*8 step,postinto,duree,sortie,dto,dt,postint
        real*8 expnu,dexpnu,xnu,mat(nb_ion,nb_ion)
        real*8 chim(nb_ion,nb_ion),y0(nb_ion),y1(nb_ion)
        real N2new(npt),N1new(npt),Nenew(npt)
        real U2new(npt),U1new(npt),Uenew(npt)
        real T2new(npt),T1new(npt),Tenew(npt)
        real q2new(npt),q1new(npt),qenew(npt)
        real N4new(npt),N6new(npt),N5new(npt),Nmnew(npt)
        real N3new(npt),N3old(npt)
        real xn1(npt),xn2(npt),xn3(npt),xne(npt)
        real xn4(npt),xn5(npt),xn6(npt),xnm(npt),Nliminf
        real xn1_1(npt),xn2_1(npt),xn3_1(npt),xne_1(npt)
        real xn4_1(npt),xn5_1(npt),xn6_1(npt),xnm_1(npt)
        real N1tot,N2tot,N3tot,N4tot,N5tot,N6tot
        real N1totold,N2totold,N3totold,N4totold,N5totold,N6totold

        real T1pnew(npt),T1tnew(npt),T2pnew(npt),T2tnew(npt)
        real T3pnew(npt),T3tnew(npt),Tmpnew(npt),Tmtnew(npt)
        real Tepnew(npt),Tetnew(npt)
        real T1pold(npt),T1told(npt),T2pold(npt),T2told(npt)
        real T3pold(npt),T3told(npt),Tmpold(npt),Tmtold(npt)
        real Tepold(npt),Tetold(npt)

        real facteur0(npt),facteur1(npt),facteur2(npt),Tmoy,Umn2
        real latgeo,longeo,latgeo_ini,longeo_ini

!-------Optical calculation variables
        real Po1d(npt),phdisso2(npt),No1dnew(npt),Uo1dnew(npt)
        real xno1d(npt)
        real nuo1di(npt),nuo1dj(npt),nuo1dk(npt),nuo1dl(npt)
        real nuo1dm(npt),nuo1dn(npt),nuo1de(npt)
       
        real Po1s(npt),Pn2a(npt),No1snew(npt),Nn2anew(npt)
        real Po1sdisso2(npt)
 
        real No1dold(npt),Uo1dold(npt),velo1dc(npt),velo1dm(npt)
        real nuo1dN2(npt),nuo1dO2(npt),nuo1dO(npt)

        !Allowed transitions
        real P1NG(npt),P2PG(npt),P1PG(npt),PO3p3P(npt),PMein(npt)
        real Po3p5P(npt)

        !Metastable
        real Noii2p(npt),Poii2P(npt)

        common/photodissociation/phdisso2,Po1sdisso2

        !Production rate common blocks for optical calculations
        real prate(10,16,202)
        integer jsg(11)
        real fhemd(400,201)
        real fhemu(400,201)
        real e(400)
        integer nen,nspec,nalt
        common/pexc/prate,fhemd,fhemu,e
        common/excinds/jsg
        common/kininds/nalt,nspec,nen

        logical exval
!-------MZ

        real Te1new(npt),qe1new(npt)
        real C2qa(npt),D2qa(npt),C2qb(npt),D2qb(npt),D3q(npt),D7q(npt)
        real C2ea(npt),D2ea(npt),C2eb(npt),D2eb(npt),D3e(npt),D7e(npt)
        real T1_15(npt),T2_15(npt),T3_15(npt),Tm_15(npt),Te_15(npt)

        real*4 buffer((ncol0+1)*npt)

        real Pri,Loi,Prj,Loj,Pri1,Pri2,Prj1,Prj2,Prj3
        real Pri0(npt),Prj0(npt),Loi0(npt),Loj0(npt)

        real Prk,Prl,Prl1,Prl2,Prm,Prm1,Prm2,Prm3,Lok,Lol,Lom
        real Prl3,Prm4,Prn,Lon,Prm5
        real Nplus
        real Nplus1(npt)

        real Tn(npt),Nh(npt),No(npt),Nn2(npt),No2(npt)
        real Un(npt),Vn(npt),Wn(npt)

        real NOHot(npt),TnOHot(npt),q_NOHot(npt)                !MZ
        real TrHot                                !MZ

        real N_0,T_0,P_0,Ci0,Cj0,Ce0,Ck0,Cl0,Cm0,Qi0,Qj0,Qe0


        real t0,R0,G0

        real Ph(npt),Po(npt),Po2(npt),Pn2(npt),Pn(npt)

        real Vm(npt),Wm(npt),omega(npt),Vm_2(npt)

        real nu_omega,coef_conv
        real nuji(npt),nujj(npt),nuje(npt),nujO(npt)
        real nujk(npt),nujl(npt),nujm(npt)
        real nujO2(npt),nujN2(npt)
        real nuin(npt),nujn(npt),nuen(npt)

        real nuki(npt),nukj (npt),nuke (npt),nukk(npt),nukl(npt)
        real nukO (npt),nukO2(npt),nukN2(npt),nukH(npt),nukm(npt)
        real nukn(npt),nuln(npt),numn(npt)

        real nuli(npt),nulj (npt),nule (npt),nulk(npt),null(npt)
        real nulO (npt),nulO2(npt),nulN2(npt),nulH(npt),nulm(npt)

        real numi(npt),numj (npt),nume (npt),numk(npt),numl(npt)
        real numO (npt),numO2(npt),numN2(npt),numH(npt),numm(npt)

        real nuni(npt),nunj (npt),nune (npt),nunk(npt),nunl(npt)
        real nunO (npt),nunO2(npt),nunN2(npt),nunH(npt),nunm(npt)
        real nunn(npt)  

        !2 cell convection
        real tmotion,omegaf,lonstart
        !-MZ

        real nuiOHot(npt),nujOHot(npt),nukOHot(npt),nulOHot(npt)   !MZ
        real numOHot(npt),nunOHot(npt),nueOHot(npt)                !MZ
        
        real amb,kb,ak1,ak2,ak3,ak4,ak5,kjN2,kjO2,TjN2,TjO2
        real Tperp(npt),Ter_1(npt)
        real Terh,Ter,dTen,dTen_1
        real akk1,akk2,akk3,akk4,akk5
        real akl1,akl2,akl3,akl4,akl5
        real akm1,akm2,akm3,akm4,akm5
        real nuN2O2,nuN2N2,nuN2O,TN2O,TN2N2,kN2O,Tjr,Tr,Trh,lTr,amN2
        real Tmr,Tr1
        real klO2,kmN2,TlO2,TmN2,TONO,kONO,TN2O2
        real cofterp

        real nuii(npt),nuij(npt),nuie(npt),nuiH(npt)
        real nuik(npt),nuil(npt),nuim(npt)

        real nuiO(npt),nuiO2(npt),nuiN2(npt)

        real nueH(npt),nueO(npt),nueO2(npt),nueN2(npt)
        real nuei(npt),nuej(npt),nuee(npt)
        real nuek(npt),nuel(npt),nuem(npt)

        real Lenrot,LeN2vib,LeO2vib,LeOexc,dd,f,gg,h,Len(npt)
        real Heat(npt),LeOfin
        real Aq(3),Bq(3),Cq(3),Dx(3),Eq(3),epsq(3),Sq(3)


        real q_Nh(npt),q_No(npt),q_No2(npt),q_Nn2(npt),q_Nn(npt)

        real Nn(npt)
        real Nnoold(npt),Nnonew(npt)
        real Unoold(npt),Unonew(npt)
        real Prno(npt),Lono(npt)
        real factno
        real factionm
        real Velnoc(npt),Velnom(npt)


        real Niqi0,Njqj0,Neqe0,Nkqk0,Nlql0,Nmqm0

        real qen

        character(len=80) filetemps
        logical flgne,flgnorm
        real N2old(npt),N1old(npt),Neold(npt)
        real N4old(npt),N6old(npt),N5old(npt)
        real U2old(npt),U1old(npt),Ueold(npt)
        real T2old(npt)
        real T1old(npt)
        real q2old(npt),q1old(npt)
        real kstep(npt)
        real alt_geo (npt),alt_geo_1(npt),G(npt),alt(npt)
        real Radn(npt)
        real mi,mj,me,dr,Re,z0,deltat,deltat_2,deltat_4,zflu,Rflu

        real mk,ml,mm,mz

        real Cz0,Qz0
        real Cze,Czi,Czj,Ciz,Cjz,Cez,Cnz,Czn
        real Umnew(npt),Umold(npt)
        real Tmnew(npt),Tmold(npt)
        real Velmm(npt),Velme(npt),Velmc(npt)
       real Prk0(npt),Lok0(npt),Prl0(npt),Lol0(npt),Prm0(npt),Lom0(npt)
        real D3(npt),D7(npt)

!CCCC
        real Tn2(npt)
!CCCC


        real Cn0,Qn0,Nnqn0
        real Cne,Cni,Cnj,Cin,Cjn,Cen
        real U3new(npt),U3old(npt)
        real T3new(npt),T3old(npt)
            real Teold(npt),qeold(npt)
        real q3new(npt),q3old(npt)
        real Velnm(npt),Velne(npt),Velnc(npt)
        real Velnq(npt)
        real Prn0(npt),Lon0(npt)
        real mn,Qntop
        real nuni0,nunj0,nune0,nunk0,nunl0,nunm0,nunn0

        real coefg(npt)
        real zeddy,Eddy,nueddy,nun,Dmol,Ntot,rhotot,Diff,mred,mmoy

        real tex0,dte

        real coefqk,coefql,coefqm

        real D2a(npt),D2b(npt),D2c(npt),C2a(npt),C2b(npt),C2c(npt)
        real D2d(npt),C2d(npt)
        real Sc(npt)

        real rbc,lbc,Vr,Vl,D2ar,D2al,D2br,D2bl,D2cr,D2cl,D2dr,D2dl
        real zero(npt),fact
        real coefelec,coefini
        real Nibot,Njbot,Uibot,Ujbot,Tibot,Tjbot,ones(npt)
        real Tebot,Qibot,Qjbot,Qebot,Qitop,Qjtop,Qetop
        real f107(3),ap(7),d(8),t(2),zhplus,zinf,zsup,stl,sec,JJ(npt)
        real vartemp,prodtemp
        real ai(6),bij(6,5),ci(6),cis(6),af(6,npt)

        real J0,B0,sinI0,cosI0

        real nu_0,kp,kptime

        real Velic(npt),Velim(npt),Velie(npt),Veliq(npt)
        real Veljc(npt),Veljm(npt),Velje(npt),Veljq(npt)
        real Velee(npt),Veleq(npt)
        real Velec(npt)

        real extra(5),x2b
        real burn
        integer lenc

        real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0,chi0
        real Bmag,dipangle,Enord,Eest,vperpnord,vperpest,vhorizon,vpara
        real ddp,Jtop
        integer ikp
        real*8 dlongeo,dlatgeo,pot,dlonmag,dlatmag,dlonref
        real dTinf,dUinf,dtz,gradTn
        real Nes(npt),Jes(npt),Tes(npt),qes(npt)
        real el_dens(npt),O_plus(npt),N2_plus(npt),NO_plus(npt)
        real O2_plus(npt),temp_o(npt),temp_m(npt),temp_e(npt)

        real rx0
        real*8 dx0
        integer i1
        real thermodiff,thermacc(npt)

        integer fid_temp
        integer kiappel,file_cond
!
        integer itype

        !added for conductivity calculations
        real sigmaP(npt),sigmaH(npt),sigmaC(npt)
        !-----MZ

        !Coefficients for molecular species
        real cofn2,cofo2        

        real or1

        !Store the convection electric field magnitude, f10.7 and ap indices
        real Econvmag, f107dat, f107Adat, apdat
        real incE

        common/neutral/Nh,No,No2,Nn2,Nn,Tn,Un,Vn,Wn
        common/param/    N_0,T_0,P_0,Ci0,Cj0,Ck0,Cl0,Cm0,Ce0,                
     &            Qi0,Qj0,Qe0
        common/adim/R0,t0,G0
        common/prodion/Ph,Po,Po2,Pn2,Pn
        common/champ/Vm,Wm,omega,Vm_2
        common/heate/Heat
        common/fluneut/q_Nh,q_No,q_No2,q_Nn2,q_Nn
        common/limite/indlim,indlim_1,zlim,zlim_1

        common/OHot/ NOHot,TnOHot,q_NOHot                    !MZ

        common/J_aligne/JJ
        common/fluxtop/Qetop


        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi0,         
     &                  Fe0,Ee0,Fi0,Ei0,                                
     &                  Bmag,dipangle,Enord,Eest,                       
     &                  vperpnord,vperpest,vhorizon,vpara,ddp,Jtop

!-------Control of MSIS90 dens.
        common/buff2/cofn2,cofo2
!-------MZ

        common/ventexo/dUinf
        common/exo/dTinf
        common/grad/dtz

        common /supra/Nes,Jes,Tes,qes

        common /noel/el_dens,O_plus,N2_plus,NO_plus,                    
     &               O2_plus,temp_o,temp_m,temp_e

        common /avril/vparaB

        common  /atm/flagatmos
!
        common/fl/flag

        real zNOHotRef,pctOHot,NOHotRef                        !MZ
        real zTnOHotRef,TnOHotRef,TnOHotDecay,TnOHotInf                !MZ
        logical chkOHot,flagOHot                        !MZ
        common /flagsOHot/ chkOHot,flagOHot                    !MZ
        common /OHotParams/ zNOHotRef,pctOHot,NOHotRef,zTnOHotRef,          
     &        TnOHotRef,TnOHotDecay,TnOHotInf

        real Eprec,Fprec
        common /RBR/ Eprec,Fprec
        real tstartprec,tstopprec
        common /prectimes/ tstartprec,tstopprec
        real Emin,Emax,Flux0,Eprev,Enext
        common /einterval/ Flux0,Emin,Emax,Eprev,Enext

!-------precipitating electrons
        real timeser(1024), edist(1024)
        real fluxdist(1024,1024)
        integer ntimeser,nfluxdist,precint,precext
        common /precdist/ nfluxdist,ntimeser,timeser,                   
     &                  edist,fluxdist,precint,precext
        integer ioerr
        character*80 prec_fname 
        real fluxstat,timestat
        character*80 tempchar
!-------MZ

        include 'TRANSPORT.INC'
    

        data rx0,dx0,i1/0.,0.d0,1/
        
        data amu/1.667e-24/
        data Aq/7.883e-6,9.466e-6,1.037e-7/
        data Bq/1.021,.8458,1.633/
        data Cq/1.009,.9444,1.466/
        data epsq/.02,.028,.008/
        data Eq/227.,326.6,98.9/
        data Sq/2.98e-23,1.91e-23,1.32e-27/

        data zeddy/100./

        data tex0,dte/20000.,0.005/

        data fid_temp/73/
        
        data B0,sinI0,cosI0/.542,.9848,.1736/
! 


!       Ignore the hot O stuff for now
!      chkOHot=.true.
!      flagOHot=.false.

!-------Read TRANSCAR input file and grab initial conditions information

        itype=4                         !machine IBM

        open(transcar_dat, file='dir.input/DATCAR', status='old')
            rewind(transcar_dat)
            read(transcar_dat,*)kiappel
            read(transcar_dat,1010) gridfn
1010        format(a)
            gridfn = split(gridfn,' ')

!-------MZ
        print *, 'reading parameters from file: ',gridfn
        open(unfic_in_transcar,file='dir.input/'//gridfn,
     &        form='unformatted',          
     &                access='direct',status='old',recl=4*2*ncol0)


        flgne=.false.

        read(unfic_in_transcar,rec=1) (buffer(i),i=1,2*ncol0)

        nx      =buffer( 1)
        ncol    =buffer( 2)
        iannee  =buffer( 3)
        imois   =buffer( 4)
        ijour   =buffer( 5)
        iheure  =buffer( 6)
        iminute =buffer( 7)
        seconde =buffer( 8)
        intpas  =buffer( 9)
        longeo  =buffer(10)
        latgeo  =buffer(11)
        lonmag  =buffer(12)
        latmag  =buffer(13)
        tmag    =buffer(14)
        f107(2) =buffer(15)
        f107(3) =buffer(16)
        ap(2)   =buffer(17)
        ikp     =buffer(18)
        dTinf   =buffer(19)
        dUinf   =buffer(20)
        cofo    =buffer(21)
        cofh    =buffer(22)
        cofn    =buffer(23)
!        Enord   =buffer(31)            Get this from the DATCAR input file
!        Eest    =buffer(32)            Get this from DATCAR
        iapprox =buffer(37)

        ipos_z=1
        ipos_n1=2
        ipos_n2=3
        ipos_n3=4
        ipos_n4=5
        ipos_n5=6
        ipos_n6=7
        ipos_u1=8
        ipos_u2=9
        ipos_u3=10
        ipos_um=11
        ipos_ue=12
        if (iapprox.eq.13) then
          ipos_t1p=13
          ipos_t1t=14
          ipos_t2p=15
          ipos_t2t=16
          ipos_t3p=17
          ipos_t3t=18
          ipos_tmp=19
          ipos_tmt=20
          ipos_tep=21
          ipos_tet=22
          ipos_q1=23
          ipos_q2=24
          ipos_q3=25
          ipos_qe=26
          ipos_nno=27
          ipos_uno=28
      ipos_po=29
      ipos_ph=30
      ipos_pn=31
      ipos_pn2=32
      ipos_po2=33
      ipos_heat=34
      ipos_no=35
      ipos_nh=36
      ipos_nn=37
      ipos_nn2=38
      ipos_no2=39
      ipos_tn=40
      ipos_un=41
      ipos_vn=42
      ipos_wn=43
      ipos_nes=44
      ipos_jes=45
      ipos_tes=46
      ipos_qes=47
      else
          ipos_t1p=13
          ipos_t1t=13
          ipos_t2p=14
          ipos_t2t=14
          ipos_t3p=15
          ipos_t3t=15
          ipos_tmp=16
          ipos_tmt=16
          ipos_tep=17
          ipos_tet=17
          ipos_q1=18
          ipos_q2=19
          ipos_q3=20
          ipos_qe=21
          ipos_nno=22
          ipos_uno=23
      ipos_po=24
      ipos_ph=25
      ipos_pn=26
      ipos_pn2=27
      ipos_po2=28
      ipos_heat=29
      ipos_no=30
      ipos_nh=31
      ipos_nn=32
      ipos_nn2=33
      ipos_no2=34
      ipos_tn=35
      ipos_un=36
      ipos_vn=37
      ipos_wn=38
      ipos_nes=39
      ipos_jes=40
      ipos_tes=41
      ipos_qes=42
      endif

        if(ncol>60) then
          ipos_po1d=61
          ipos_no1d=62
          ipos_uo1d=63
        endif

      close(unfic_in_transcar)


        nligne=nx+2
        longbuf=nligne*ncol
        longrec=itype*longbuf

        open(unfic_in_transcar,file='dir.input/'//gridfn,
     &        form='unformatted',              
     &                  access='direct',status='old',recl=longrec)

        N_0=1.e4
        T_0=1000.
!cc
![      Different altitude levels used in the program

        kb=   1.38e-16
        Re  =637800000.
        Rekm=Re*1.e-5
        z0  = 100000000.
        zflu= 600.
        zinf= 20000000.
        zsup=862200000.
        zhplus  = 400.
        cofterp=3.*amu/kb
        Nliminf=1.e-1/N_0
        T_min=100./T_0
        r_min=1.e-33

!
![      Altitude and Time step initialisation
!
        dr=(zsup-zinf)/(nx-1.)
        ielecini=10
        coefini=1./float(ielecini)
        coefelec=coefini
!]      
![      Physical constant initialisation

        amu=1.66e-24
        mi=   1.66e-24
        mj=16.*1.66e-24
        me=   9.11e-28

        mk=28.*1.66e-24
        ml=32.*1.66e-24
        mm=30.*1.66e-24
        mn=14.*1.66e-24
        mz=3./(1./mk+1./ml+1./mm)


        omi=mi/amu
        omj=mj/amu
        omn=mn/amu
        omz=mz/amu
        ome=me/amu

        tex0=tex0+dte

!]      

![      Normalization coefficients initialisation

        P_0=kb*N_0*T_0

        G0 =980.665
        Ci0=sqrt(kb*T_0/mi)
        Cj0=sqrt(kb*T_0/mj)
        Ck0=sqrt(kb*T_0/mk)
        Cl0=sqrt(kb*T_0/ml)
        Cm0=sqrt(kb*T_0/mm)
        Cn0=sqrt(kb*T_0/mn)
        Cz0=sqrt(kb*T_0/mz)
        Ce0=sqrt(kb*T_0/me)
        t0=Ci0/G0
        R0=Ci0*t0




        Qitop=0.
        Qjtop=0.
        Qntop=0.

        Cji=Cj0/Ci0
        Cki=Ck0/Ci0
        Cli=Cl0/Ci0
        Cmi=Cm0/Ci0
        Cni=Cn0/Ci0
        Czi=Cz0/Ci0
        Cei=Ce0/Ci0

        Cij=Ci0/Cj0
        Ckj=Ck0/Cj0
        Clj=Cl0/Cj0
        Cmj=Cm0/Cj0
        Cnj=Cn0/Cj0
        Czj=Cz0/Cj0
        Cej=Ce0/Cj0

        Cik=Ci0/Ck0
        Cjk=Cj0/Ck0
        Clk=Cl0/Ck0
        Cmk=Cm0/Ck0
        Cnk=Cn0/Ck0
        Czk=Cz0/Ck0
        Cek=Ce0/Ck0

        Cil=Ci0/Cl0
        Cjl=Cj0/Cl0
        Ckl=Ck0/Cl0
        Cml=Cm0/Cl0
        Cnl=Cn0/Cl0
        Czl=Cz0/Cl0
        Cel=Ce0/Cl0

        Cim=Ci0/Cm0
        Cjm=Cj0/Cm0
        Ckm=Ck0/Cm0
        Clm=Cl0/Cm0
        Cnm=Cn0/Cm0
        Czm=Cz0/Cm0
        Cem=Ce0/Cm0

        Cin=Ci0/Cn0
        Cjn=Cj0/Cn0
        Ckn=Ck0/Cn0
        Cln=Cl0/Cn0
        Cmn=Cm0/Cn0
        Czn=Cz0/Cn0
        Cen=Ce0/Cn0

        Ciz=Ci0/Cz0
        Cjz=Cj0/Cz0
        Ckz=Ck0/Cz0
        Clz=Cl0/Cz0
        Cmz=Cm0/Cz0
        Cnz=Cn0/Cz0
        Cez=Ce0/Cz0

        Cie=Ci0/Ce0
        Cje=Cj0/Ce0
        Cke=Ck0/Ce0
        Cle=Cl0/Ce0
        Cme=Cm0/Ce0
        Cne=Cn0/Ce0
        Cze=Cz0/Ce0

        Qi0=P_0*Ci0
        Qj0=P_0*Cj0
        Qk0=P_0*Ck0
        Ql0=P_0*Cl0
        Qm0=P_0*Cm0
        Qz0=P_0*Cz0
        Qn0=P_0*Cn0
        Qe0=P_0*Ce0

        Niqi0 =N_0/Qi0
        Njqj0 =N_0/Qj0
        Neqe0 =N_0/Qe0
        Nkqk0 =N_0/Qz0
        Nlql0 =N_0/Qz0
        Nmqm0 =N_0/Qz0
        Nnqn0 =N_0/Qn0

        nu_0=T_0**1.5/N_0/t0


        Rcira=(z0+Re)/R0
        Rflu=(zflu+Re)/R0

        thermodiff=.6*54.5/kb*me/T_0**2.5*t0*Qe_0




!------Continue reading input file
        read(transcar_dat,*)dto         ! pas d'integration numerique
        read(transcar_dat,*)sortie      ! intervalle de temps entre deux sorties 
        read(transcar_dat,*)iyd_ini     ! date de la periode simulee
        read(transcar_dat,*)tempsini    ! UT de debut (en secondes)
        read(transcar_dat,*)tempslim    ! UT limite (en secondes)
        read(transcar_dat,*)jpreci
        read(transcar_dat,*)latgeo_ini,longeo_ini

        tempsconv_1=0.d0
        read(transcar_dat,*)tempsconv_1    ! duree de la convection en amont en secondes (<= 0 si pas de convection)
        tempsconv=0.d0
        read(transcar_dat,*)tempsconv    ! duree de la convection en aval en secondes (<= 0 si pas de convection)
        read(transcar_dat,*)step        ! intervalle de temps entre deux tubes
        read(transcar_dat,*)postinto    ! intervalle de temps entre deux appel a transelec
        vparaB=0.
        read(transcar_dat,*)vparaB      ! transport induit le long de la ligne de champ

        if (latgeo_ini.gt.90.) then                                                ! Modif DA 0202 2001 - debut
           nb_position=0                                                ! Modif DA 0202 2001
           multi_position=.true.                                        ! Modif DA 0202 2001
           read(transcar_dat,*)latgeo_ini,longeo_ini                    ! Modif DA 0202 2001
           do while (latgeo_ini.le.90.)                                 ! Modif DA 0202 2001
              nb_position=nb_position+1                                 ! Modif DA 0202 2001
              latgeo_position(nb_position)=latgeo_ini                   ! Modif DA 0202 2001
              longeo_position(nb_position)=longeo_ini                   ! Modif DA 0202 2001
              lecture_lat_lon=' '                                       ! Modif DA 0202 2001
              read(transcar_dat,*)latgeo_ini,longeo_ini                 ! Modif DA 0202 2001
           enddo                                                        ! Modif DA 0202 2001
       else                                                             ! Modif DA 0202 2001
           nb_position=1                                                ! Modif DA 0202 2001
           multi_position=.false.                                       ! Modif DA 0202 2001
       endif                                                            ! Modif DA 0202 2001 - fin
!------MZ

!-------Read in additional adjustable simulation parameters
        read(transcar_dat,*) f107dat
        read(transcar_dat,*) f107Adat
        read(transcar_dat,*) apdat
        read(transcar_dat,*) Econvmag           !not using convection code so store field in Enord for frictional heating calculations
        incE=1.
        read(transcar_dat,*) cofo
        read(transcar_dat,*) cofn2
        read(transcar_dat,*) cofo2
        read(transcar_dat,*) cofn
        read(transcar_dat,*) cofh

        read(transcar_dat,*) Qetop
        Qetop=Qetop/Qe0                         !normalize topside heat flux

        read(transcar_dat,1010) prec_fname
        prec_fname = split(prec_fname,' ')
        
        read(transcar_dat,*) precint
        read(transcar_dat,*) precext
!        read(transcar_dat,*) tstartprec
!        read(transcar_dat,*) tstopprec
!-------MZ

!-------Close input file
        close(transcar_dat)
!-------MZ

!-------Read file containing electron precipitation information, this is messy
!and should be rewritten at some point.
        prec_fname='dir.input/'//prec_fname

       print *, 'reading PRECIPITATION parameters from file: ',
     &   prec_fname
        open(313,file=prec_fname,status='old')
 
        ioerr=1; ntimeser=1; timestat=1;
        do while(timestat .ge. 0 .and. ioerr>-1 .and.                  
     &            ntimeser .le. 1024)                         !this loops over the different distributions in the time series
          read(313,*,iostat=ioerr) timeser(ntimeser)
          print*,timeser(ntimeser)

          m=1; fluxstat=1.e0
          do while(fluxstat .ge. 0. .and. m .le. 1024)        !this loops over the different energy bins in each distribution
            read(313,*,iostat=ioerr) timestat,                          
     &                  fluxdist(m,ntimeser) 
            print*,timestat,fluxdist(m,ntimeser)

            if (timestat .ge. 0.) then
              edist(m)=timestat
              nfluxdist=m
            endif
            fluxstat=fluxdist(m,ntimeser)

            m=m+1
          enddo
          print*,'Number of energy bins for precip.:  ',nfluxdist

          ntimeser=ntimeser+1
        enddo
        ntimeser=ntimeser-1
        print*,'Number of time series in precip.:  ',ntimeser
        close(313)
!        print*,'Interpolation:  ',precint
!        print*,'Extrapolation:  ',precext
        tstartprec=timeser(1)
        tstopprec=timeser(ntimeser)
!-------MZ



        flgconv=.false.
        if (tempsconv_1.gt.0.) then
          flgconv=.true.
        endif
        
        if (tempsconv.gt.0.) then
          flgconv=.true.
        endif
        
        if (flgconv) then
          itube=-1
          ntubmax=int(tempslim/step)+1
           ncol1=47
              longbuf1=nligne*ncol1
          longrec1=itype*longbuf

          open(unfic_out_transcar,file='dir.output/transcar_output',
     &         form='unformatted',access='direct',recl=longrec1,
     &           status='new')


        else
          itube=-2
          ntubmax=-2
        endif

        nrec_lec=0
        nrec_ecr=0
        ierr=0
        do while (itube.le.ntubmax)    ! debut de la boucle sur les lignes de champs "beginning of the loop on the field lines"

        itube=itube+1
          do i_position=1,nb_position             ! Modif DA 0202 2001

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
![      Densities and Velocities profiles initialisation
!ccccccccccccccccccccccccccccccccccccccccccccccccccc

!      if (ierr.ge.0) nrec_lec=nrec_lec+1
        nrec_lec=1
        read(unfic_in_transcar,rec=nrec_lec,iostat=ierr)
     &            (buffer(i),i=1,longbuf)


!        do i=1,128
!        j=min(max(i-22,1),nx)
!             ipos=(j+1)*ncol
!        alt(i)=90.*(3000./90.)**((i-1)/127.)
!        if (i.eq.128) alt(i)=3000.
        do i=1,nx
             ipos=(i+1)*ncol
             alt(i)      =buffer(ipos+ipos_z)

             N1new(i)    =buffer(ipos+ipos_n1)/N_0*1.e-6
             N2new(i)    =buffer(ipos+ipos_n2)/N_0*1.e-6
             N3new(i)    =buffer(ipos+ipos_n3)/N_0*1.e-6
             N4new(i)    =buffer(ipos+ipos_n4)/N_0*1.e-6
             N5new(i)    =buffer(ipos+ipos_n5)/N_0*1.e-6
             N6new(i)    =buffer(ipos+ipos_n6)/N_0*1.e-6
             U1new(i)    =buffer(ipos+ipos_u1)/Cj0*1.e2
             U2new(i)    =buffer(ipos+ipos_u2)/Ci0*1.e2
             U3new(i)    =buffer(ipos+ipos_u3)/Cn0*1.e2
             Umnew(i)    =buffer(ipos+ipos_um)/Cz0*1.e2
             Uenew(i)    =buffer(ipos+ipos_ue)/Ce0*1.e2
             T1pnew(i)    =buffer(ipos+ipos_t1p)/T_0
             T1tnew(i)    =buffer(ipos+ipos_t1t)/T_0
             T2pnew(i)    =buffer(ipos+ipos_t2p)/T_0
             T2tnew(i)    =buffer(ipos+ipos_t2t)/T_0
             T3pnew(i)    =buffer(ipos+ipos_t3p)/T_0
             T3tnew(i)    =buffer(ipos+ipos_t3t)/T_0
             Tmpnew(i)    =buffer(ipos+ipos_tmp)/T_0
             Tmtnew(i)    =buffer(ipos+ipos_tmt)/T_0
             Tepnew(i)    =buffer(ipos+ipos_tep)/T_0
             Tetnew(i)    =buffer(ipos+ipos_tet)/T_0
             q1new(i)    =buffer(ipos+ipos_q1)/Qj0*1.e3
             q2new(i)    =buffer(ipos+ipos_q2)/Qi0*1.e3
             q3new(i)    =buffer(ipos+ipos_q3)/Qn0*1.e3
             qenew(i)    =buffer(ipos+ipos_qe)/Qe0*1.e3
             Nnonew(i)   =buffer(ipos+ipos_nno)/N_0*1.e-6
             Unonew(i)   =buffer(ipos+ipos_uno)/Cm0*1.e2
         Po(i)     =buffer(ipos+ipos_po)*1.e-6
         Ph(i)     =buffer(ipos+ipos_ph)*1.e-6
!         Pn(i)     =buffer(ipos+ipos_pn)*1.e-6
         Pn(i)     =0.
         Pn2(i)     =buffer(ipos+ipos_pn2)*1.e-6
         Po2(i)     =buffer(ipos+ipos_po2)*1.e-6
         Heat(i)     =buffer(ipos+ipos_heat)

             if(ncol>60) then
                Po1d(i)=buffer(ipos+ipos_po1d)*1.e-6
                No1dnew(i)=buffer(ipos+ipos_no1d)/N_0*1.e-6
                Uo1dnew(i)=buffer(ipos+ipos_uo1d)/Cj0*1.e2
             else
                Po1d(i)=0.
                No1dnew(i)=N1new(i)!buffer(ipos+ipos_no)*1.e-6*1.e-4
                Uo1dnew(i)=U1new(i)!0.
             endif

cccccccccccccccccccccccccccccccccccccccccccccccccccc

                Nenew(i)=N2new(i)+N1new(i)+N4new(i)+N6new(i)
     &           +N5new(i)+N3new(i)

                Nmnew(i)=N4new(i)+N5new(i)+N6new(i)

         T1new(i)=(T1pnew(i)+2.*T1tnew(i))/3.
         T2new(i)=(T2pnew(i)+2.*T2tnew(i))/3.
         T3new(i)=(T3pnew(i)+2.*T3tnew(i))/3.
         Tmnew(i)=(Tmpnew(i)+2.*Tmtnew(i))/3.
         Tenew(i)=(Tepnew(i)+2.*Tetnew(i))/3.

!         Heat(i)= Heat(i)/Nenew(i)/N_0*10.
         Heat(i)= Heat(i)*10.
!]      
![           Time independant parameters initialisation

         alt_geo(i)=(Re+100000.*alt(i)/sinI0)/R0
         alt_geo_1(i)=.5*(1.+tanh((alt(i)-300.)/20.))/alt_geo(i)
         alt_geo_1(i)=1./alt_geo(i)
         G(i)=(Re/R0/alt_geo(i))**2*sinI0
         zero(i)=0.
         ones(i)=1.

         if (alt(i).le.500.) indno=i
!]      

        enddo
    
!        nx=128
!        nligne=nx+2
        ncol=63
        longbuf=nligne*ncol
        longrec=itype*longbuf

      iapprox=13

        ipos_t1p=13
        ipos_t1t=14
        ipos_t2p=15
        ipos_t2t=16
        ipos_t3p=17
        ipos_t3t=18
        ipos_tmp=19
        ipos_tmt=20
        ipos_tep=21
        ipos_tet=22
        ipos_q1=23
        ipos_q2=24
        ipos_q3=25
        ipos_qe=26
        ipos_nno=27
        ipos_uno=28
        ipos_po=29
        ipos_ph=30
        ipos_pn=31
        ipos_pn2=32
        ipos_po2=33
        ipos_heat=34
        ipos_no=35
        ipos_nh=36
        ipos_nn=37
        ipos_nn2=38
        ipos_no2=39
        ipos_tn=40
        ipos_un=41
        ipos_vn=42
        ipos_wn=43
        ipos_nes=44
        ipos_jes=45
        ipos_tes=46
        ipos_qes=47
        ipos_NOHot=48
        ipos_TnOHot=49

        ipos_po1d=61
        ipos_no1d=62
        ipos_uo1d=63

![      Grid initialisation


    
        np=nx+1
        Ipos1=1
        Iposn=nx
        Iposnp=np
        alpha=1
        
        
        T1pnew(np)=T1pnew(nx)
        T1tnew(np)=T1tnew(nx)
        T1new(np)=T1new(nx)
        T2pnew(np)=T2pnew(nx)
        T2tnew(np)=T2tnew(nx)
        T2new(np)=T2new(nx)
        T3pnew(np)=T3pnew(nx)
        T3tnew(np)=T3tnew(nx)
        T3new(np)=T3new(nx)
        Tmpnew(np)=Tmpnew(nx)
        Tmtnew(np)=Tmtnew(nx)
        Tmnew(np)=Tmnew(nx)
        Tepnew(np)=Tepnew(nx)
        Tetnew(np)=Tetnew(nx)
        Tenew(np)=Tenew(nx)
          Tenew(np)=2.*Tenew(nx)-Tenew(nx-1)
        q1new(np)=q1new(nx)
        q2new(np)=q2new(nx)
        q3new(np)=q3new(nx)
        qenew(np)=qenew(nx)
    
    
    
        Radn(1)=(3.*alt_geo(1)-alt_geo(2))/2.
        Radn(np)=(3.*alt_geo(nx)-alt_geo(nx-1))/2.
        
            do i=2,nx
             Radn(i)=.5*(alt_geo(i)+alt_geo(i-1))
            enddo

        coef_flnp=(alt_geo(nx)/Radn(np))**3
        coef_fln =(Radn(nx)/Radn(np))**3

        call makegrid(Radn,Radn,Ipos1,Iposnp,alpha)


        extra(1)=alt_geo(nx-2)
        extra(2)=alt_geo(nx-1)
        extra(3)=alt_geo(nx)
        extra(4)=(alt_geo(nx-2)+alt_geo(nx-1)+alt_geo(nx))/3.
        x2b=(alt_geo(nx-2)**2+alt_geo(nx-1)**2+alt_geo(nx)**2)/3.
        extra(5)=x2b-extra(4)**2


!]


! determination du point de depart de la convection

         if (multi_position) then                       ! Modif DA 0202 2001 - debut
            dlongeo=longeo_position(i_position)         ! Modif DA 0202 2001
            dlatgeo=latgeo_position(i_position)         ! Modif DA 0202 2001
            longeo=longeo_position(i_position)          ! Modif DA 0202 2001
            latgeo=latgeo_position(i_position)          ! Modif DA 0202 2001
         else                                           ! Modif DA 0202 2001
             dlongeo=longeo_ini
             dlatgeo=latgeo_ini
             longeo=longeo_ini
             latgeo=latgeo_ini
         endif                                          ! Modif DA 0202 2001 - fin
    
         call geo2mag(dlatgeo,dlongeo,dlatmag,dlonmag,dlonref)
        
         iyd=iyd_ini
         iydfin=iyd_ini
         iydtube=iydfin

        if (itube.ge.0) then

      tempstube=tempsini+itube*step         ! temps de resolution du tube
      do while (tempstube.ge. 86400._dp)        !
        tempstube=tempstube-86400._dp            ! on ramene tempstube entre 0 et 24 heures
        iydtube=iydtube+1                !
      enddo                        !

      tempsfin=tempstube+tempsconv             ! temps desire en fin de convection (version mesure radar)
      iydfin=iydtube
      do while (tempsfin.ge. 86400._dp)            !
        tempsfin=tempsfin-86400._dp            ! on ramene tempsfin entre 0 et 24 heures
        iydfin=iydfin+1                !
      enddo                        !

      tempsdeb=tempstube-tempsconv_1        ! temps desire en debut de convection (version mesure radar)
      iyddeb=iydtube
      do while (tempsdeb.lt. 0._dp)            !
        tempsdeb=tempsdeb+86400._dp            ! on ramene tempsdeb entre 0 et 24 heures
        iyddeb=iyddeb-1                !
      enddo                        !

      if (tempsconv_1.gt.0) then
        call conv_ar(itube,dlongeo,dlatgeo,
     &             dlonmag,dlatmag,dlonref,
     &                 iyddeb,tempsdeb,
     &                 iydtube,tempstube,
     &                 iyd   ,temps   ,dto,postinto)
            longeo=dlongeo
            latgeo=dlatgeo
            iyd=iyddeb
          endif
          temps=tempsdeb

      else
      temps=tempsini
      tempsfin=temps+tempslim
      do while (tempsfin.ge.86400.d0)        !
        tempsfin=tempsfin-86400.d0        ! on ramene tempsfin entre 0 et 24 heures
        iydfin=iydfin+1            !
      enddo                    !
      iydtube=iydfin
      tempstube=tempsfin
      endif

!------------
! initialisation du temps et ouverture du fichier "initialization time and opening the file"
      if (debug) write(stdout,*),"initialization time, opening the file"

      tempsint=0.d0
      tempsort=0.d0
      duree=0.d0
      itmp=int(tempstube)

      flgpot=.true.
    
      stl=mod(temps/3600.+longeo/15.,24.d0)
      iannee=iyd/1000
      call jour_mois(iyd,
     &                   ijour,imois)
      iheure=int(temps/3600.)
      iminute=int(temps/60.-iheure*60.)
      seconde=mod(temps,60.d0)
          if (seconde.ge.60.) then
            iminute=iminute+1
            seconde=seconde-60.
          endif
      tu=iheure+iminute/60.
      call lec_indices(iannee,imois,ijour,tu,ap,f107)
      kp=ap2kp(ap(2))
      ikp=kp*3.

      if (debug) then
          call cpu_time(tic)
          write(stdout,*),tic,'call convec'
      endif

      call convec(iyd,temps,kp,dlongeo,dlatgeo,
     &              dlonmag,dlatmag,dlonref,dx0,pot,flgpot)
          longeo=dlongeo
          latgeo=dlatgeo
      
      if (debug) then
      call cpu_time(tic)
      write(stdout,*),tic,'transconvec: done convec, entering coskhi'
     & !,' longeo,latgeo',longeo,latgeo
      endif

      chi=acos(coskhi(latgeo,longeo,temps/3600.,imois,ijour,i1))
      chideg=chi*rad2deg

      sinI=sin(dipangle*deg2rad)
      cosI=cos(dipangle*deg2rad)


      if (flgconv) then

        filetemps='dir.output/YYYYJJJ_00000_000.mat'
        write(filetemps(12:18),'(i7)')iydtube

        if (itmp.lt.10) then
          write(filetemps(24:24),'(i1)')itmp
        else if (itmp.lt.100) then
          write(filetemps(23:24),'(i2)')itmp
        else if (itmp.lt.1000) then
          write(filetemps(22:24),'(i3)')itmp
        else if (itmp.lt.10000) then
          write(filetemps(21:24),'(i4)')itmp
        else
          write(filetemps(20:24),'(i5)')itmp
        endif

            if (multi_position) then                        ! Modif DA 0202 2001 - debut
               iw_tube=i_position                           ! Modif DA 0202 2001
            else                                            ! Modif DA 0202 2001
               iw_tube=itube                                ! Modif DA 0202 2001
            endif                                           ! Modif DA 0202 2001
            if (iw_tube.lt.10) then                         ! Modif DA 0202 2001
               write(filetemps(28:28),'(i1)')iw_tube        ! Modif DA 0202 2001
            else if (iw_tube.lt.100) then                   ! Modif DA 0202 2001
               write(filetemps(27:28),'(i2)')iw_tube        ! Modif DA 0202 2001
            else if (iw_tube.lt.1000) then                  ! Modif DA 0202 2001
               write(filetemps(26:28),'(i3)')iw_tube        ! Modif DA 0202 2001 - fin
            endif
    
      else

        filetemps='dir.output/transcar_output'

      endif



!
! definitions des parametres initiaux et sauvegarde initiale "definitions of initial parameters and initial backup"
       if (debug) then
       write(stdout,*),"definitions of initial parameters",
     & " and initial backup"
       endif
       
        flagatmos=.true.

        do i=1,nx
          el_dens(i)=Nenew(i)*N_0
          O_plus(i)=N1new(i)*N_0
          N2_plus(i)=N4new(i)*N_0
          NO_plus(i)=N5new(i)*N_0
          O2_plus(i)=N6new(i)*N_0
          temp_e(i)=Tenew(i)*T_0
          temp_o(i)=T1new(i)*T_0
          temp_m(i)=Tmnew(i)*T_0
        enddo

        if ((any(ieee_is_nan(nenew)))
     &  .or.(any(ieee_is_nan(Tenew)))
     &  .or.(any(ieee_is_nan(T1new)))) then
          write(stderr,*),'problem before calling atmos'
          goto 246
        endif
    
        if (debug) then
         call cpu_time(tic)
         write(stdout,*),tic,' transconvec: call atmos'!  latgeo=',latgeo
        endif
        
        call atmos(iyd,real(temps,sp),stl,alt,latgeo,longeo,jpreci,f107,
     &           ap,Nenew,Tenew,T1new,nx,kiappel,file_cond)

        nrectemps=1
        do i=1,longbuf
          buffer(i)=0.
        enddo
        buffer( 1)=nx
        buffer( 2)=ncol
        buffer( 3)=iannee
        buffer( 4)=imois
        buffer( 5)=ijour
        buffer( 6)=iheure
        buffer( 7)=iminute
        buffer( 8)=seconde
        buffer( 9)=intpas
        buffer(10)=longeo
        buffer(11)=latgeo
        buffer(12)=lonmag
        buffer(13)=latmag
        buffer(14)=tmag
        buffer(15)=f107(2)
        buffer(16)=f107(3)
        buffer(17)=ap(2)
        buffer(18)=kptime
        buffer(19)=dTinf
        buffer(20)=dUinf
        buffer(21)=cofo
        buffer(22)=cofh
        buffer(23)=cofn
        buffer(24)=chi0
        buffer(25)=Fe0
        buffer(26)=Ee0
        buffer(27)=Fi0
        buffer(28)=Ei0
        buffer(29)=Bmag
        buffer(30)=dipangle
        buffer(31)=Enord
        buffer(32)=Eest
        buffer(33)=vperpnord
        buffer(34)=vperpest
        buffer(35)=vhorizon
        buffer(36)=vpara
        buffer(37)=iapprox
        buffer(38)=ddp
        buffer(39)=Jtop
        buffer(60)=1
        do i=1,nx
          ipos=(i+1)*ncol
          buffer(ipos+ipos_z)=alt(i)
          buffer(ipos+ipos_n1)=N_0*N1new(i)*1.e6
          buffer(ipos+ipos_n2)=N_0*N2new(i)*1.e6
          buffer(ipos+ipos_n3)=N_0*N3new(i)*1.e6
          buffer(ipos+ipos_n4)=N_0*N4new(i)*1.e6
          buffer(ipos+ipos_n5)=N_0*N5new(i)*1.e6
          buffer(ipos+ipos_n6)=N_0*N6new(i)*1.e6
          buffer(ipos+ipos_u1) =Cj0*U1new(i)/1.e2
          buffer(ipos+ipos_u2) =Ci0*U2new(i)/1.e2
          buffer(ipos+ipos_u3)=Cn0*U3new(i)/1.e2
          buffer(ipos+ipos_um)=Cz0*Umnew(i)/1.e2
          buffer(ipos+ipos_ue)=Ce0*Uenew(i)/1.e2
          buffer(ipos+ipos_t1p)=T_0*T1pnew(i)
          buffer(ipos+ipos_t1t)=T_0*T1tnew(i)
          buffer(ipos+ipos_t2p)=T_0*T2pnew(i)
          buffer(ipos+ipos_t2t)=T_0*T2tnew(i)
          buffer(ipos+ipos_t3p)=T_0*T3pnew(i)
          buffer(ipos+ipos_t3t)=T_0*T3tnew(i)
          buffer(ipos+ipos_tmp)=T_0*Tmpnew(i)
          buffer(ipos+ipos_tmt)=T_0*Tmtnew(i)
          buffer(ipos+ipos_tep)=T_0*Tepnew(i)
          buffer(ipos+ipos_tet)=T_0*Tetnew(i)
          buffer(ipos+ipos_q1)=Qj0*q1new(i)/1.e3
          buffer(ipos+ipos_q2)=Qi0*q2new(i)/1.e3
          buffer(ipos+ipos_q3)=Qn0*q3new(i)/1.e3
          buffer(ipos+ipos_qe)=Qe0*qenew(i)/1.e3
          buffer(ipos+ipos_nno)=N_0*Nnonew(i)*1.e6
          buffer(ipos+ipos_uno)=Cm0*Unonew(i)/1.e2
          buffer(ipos+ipos_po)=Po(i)*1.e6
          buffer(ipos+ipos_ph)=Ph(i)*1.e6
          buffer(ipos+ipos_pn)=Pn(i)*1.e6
          buffer(ipos+ipos_pn2)=Pn2(i)*1.e6
          buffer(ipos+ipos_po2)=Po2(i)*1.e6
          buffer(ipos+ipos_heat)=Heat(i)/10.
          buffer(ipos+ipos_no)=No(i)*1.e6
          buffer(ipos+ipos_nh)=Nh(i)*1.e6
          buffer(ipos+ipos_nn)=Nn(i)*1.e6
          buffer(ipos+ipos_nn2)=Nn2(i)*1.e6
          buffer(ipos+ipos_no2)=No2(i)*1.e6
          buffer(ipos+ipos_tn)=Tn(i)
          buffer(ipos+ipos_un)=Un(i)/1.e2
          buffer(ipos+ipos_vn)=Vn(i)/1.e2
          buffer(ipos+ipos_wn)=Wn(i)/1.e2
          buffer(ipos+ipos_nes)=Nes(i)*1.e6
          buffer(ipos+ipos_jes)=Jes(i)*1.e4
          buffer(ipos+ipos_tes)=Tes(i)
          buffer(ipos+ipos_qes)=Qes(i)*1.e-7
          buffer(ipos+ipos_NOHot)=NOHot(i)*1.e6
          buffer(ipos+ipos_TnOHot)=phdisso2(i)*1.e6
          
          buffer(ipos+ipos_po1d)=Po1d(i)*1.e6
          buffer(ipos+ipos_no1d)=N_0*No1dnew(i)*1.e6
          buffer(ipos+ipos_uo1d)=Cj0*Uo1dnew(i)/1.e2
       enddo
       open(fid_temp,file=filetemps,
     &           form='unformatted',access='direct',recl=longrec,
     &           status='unknown')
       write(fid_temp,rec=nrectemps)(buffer(i),i=1,longbuf)
       close(fid_temp)
       flagatmos=.true.

1001       format(' temps de depart:  ',i8,' tube numero: ',i3)
           write(*,1001) int(temps),itube

        dt_max=0.
        do i=1,nx
          vmax=max(abs(U1new(i)*Ci0),abs(U2new(i)*Cj0))
          vmax=max(vmax,abs(U3new(i)*Cn0))
          vmax=max(vmax,abs(Umnew(i)*Cz0))
          dt_max=max(dt_max,vmax*rln(i))
        enddo
        dt_max=5.*dt_max/R0
    
       if (debug)  write(stdout,*),'call pas_de_temps'
       call pas_de_temps(iyd,temps,dt,postint,dto,postinto)
       tempsint=postint
       tempsort=sortie
       duree=0.

       do while (temps.le.tempsfin.or.iyd.ne.iydfin)    ! debut de la boucle temporelle sur la ligne de champ

         dt_max=0.
         do i=1,nx
           vmax=max(abs(U1new(i)*Ci0),abs(U2new(i)*Cj0))
           vmax=max(vmax,abs(U3new(i)*Cn0))
           vmax=max(vmax,abs(Umnew(i)*Cz0))
           dt_max=max(dt_max,vmax*rln(i))
         enddo
         dt_max=5.*dt_max/R0
    
        call pas_de_temps(iyd,temps,dt,postint,dto,postinto)

      if (duree.ge.tempsort) then
       nrectemps=nrectemps+1
           do i=1,longbuf
             buffer(i)=0.
           enddo
           buffer( 1)=nx
           buffer( 2)=ncol
           buffer( 3)=iannee
           buffer( 4)=imois
           buffer( 5)=ijour
           buffer( 6)=iheure
           buffer( 7)=iminute
           buffer( 8)=seconde
           buffer( 9)=intpas
           buffer(10)=longeo
           buffer(11)=latgeo
           buffer(12)=lonmag
           buffer(13)=latmag
           buffer(14)=tmag
           buffer(15)=f107(2)
           buffer(16)=f107(3)
           buffer(17)=ap(2)
           buffer(18)=kptime
           buffer(19)=dTinf
           buffer(20)=dUinf
           buffer(21)=cofo
           buffer(22)=cofh
           buffer(23)=cofn
           buffer(24)=chi0
           buffer(25)=Fe0
           buffer(26)=Ee0
           buffer(27)=Fi0
           buffer(28)=Ei0
           buffer(29)=Bmag
           buffer(30)=dipangle
           buffer(31)=Enord
           buffer(32)=Eest
           buffer(33)=vperpnord
           buffer(34)=vperpest
           buffer(35)=vhorizon
           buffer(36)=vpara
       buffer(37)=iapprox
       buffer(38)=ddp
       buffer(39)=Jtop
       buffer(60)=1
           do i=1,nx
             ipos=(i+1)*ncol
             buffer(ipos+ipos_z)=alt(i)
             buffer(ipos+ipos_n1)=N_0*N1new(i)*1.e6
             buffer(ipos+ipos_n2)=N_0*N2new(i)*1.e6
             buffer(ipos+ipos_n3)=N_0*N3new(i)*1.e6
             buffer(ipos+ipos_n4)=N_0*N4new(i)*1.e6
             buffer(ipos+ipos_n5)=N_0*N5new(i)*1.e6
             buffer(ipos+ipos_n6)=N_0*N6new(i)*1.e6
             buffer(ipos+ipos_u1) =Cj0*U1new(i)/1.e2
             buffer(ipos+ipos_u2) =Ci0*U2new(i)/1.e2
             buffer(ipos+ipos_u3)=Cn0*U3new(i)/1.e2
             buffer(ipos+ipos_um)=Cz0*Umnew(i)/1.e2
             buffer(ipos+ipos_ue)=Ce0*Uenew(i)/1.e2
             buffer(ipos+ipos_t1p)=T_0*T1pnew(i)
             buffer(ipos+ipos_t1t)=T_0*T1tnew(i)
             buffer(ipos+ipos_t2p)=T_0*T2pnew(i)
             buffer(ipos+ipos_t2t)=T_0*T2tnew(i)
             buffer(ipos+ipos_t3p)=T_0*T3pnew(i)
             buffer(ipos+ipos_t3t)=T_0*T3tnew(i)
             buffer(ipos+ipos_tmp)=T_0*Tmpnew(i)
             buffer(ipos+ipos_tmt)=T_0*Tmtnew(i)
             buffer(ipos+ipos_tep)=T_0*Tepnew(i)
             buffer(ipos+ipos_tet)=T_0*Tetnew(i)
             buffer(ipos+ipos_q1)=Qj0*q1new(i)/1.e3
             buffer(ipos+ipos_q2)=Qi0*q2new(i)/1.e3
             buffer(ipos+ipos_q3)=Qn0*q3new(i)/1.e3
             buffer(ipos+ipos_qe)=Qe0*qenew(i)/1.e3
             buffer(ipos+ipos_nno)=N_0*Nnonew(i)*1.e6
             buffer(ipos+ipos_uno)=Cm0*Unonew(i)/1.e2
             buffer(ipos+ipos_po)=Po(i)*1.e6
             buffer(ipos+ipos_ph)=Ph(i)*1.e6
             buffer(ipos+ipos_pn)=Pn(i)*1.e6
         buffer(ipos+ipos_pn2)=Pn2(i)*1.e6
         buffer(ipos+ipos_po2)=Po2(i)*1.e6
         buffer(ipos+ipos_heat)=Heat(i)/10.
         buffer(ipos+ipos_no)=No(i)*1.e6
         buffer(ipos+ipos_nh)=Nh(i)*1.e6
         buffer(ipos+ipos_nn)=Nn(i)*1.e6
         buffer(ipos+ipos_nn2)=Nn2(i)*1.e6
         buffer(ipos+ipos_no2)=No2(i)*1.e6
         buffer(ipos+ipos_tn)=Tn(i)
         buffer(ipos+ipos_un)=Un(i)/1.e2
         buffer(ipos+ipos_vn)=Vn(i)/1.e2
         buffer(ipos+ipos_wn)=Wn(i)/1.e2
         buffer(ipos+ipos_nes)=Nes(i)*1.e6
         buffer(ipos+ipos_jes)=Jes(i)*1.e4
         buffer(ipos+ipos_tes)=Tes(i)
         buffer(ipos+ipos_qes)=Qes(i)*1.e-7
         buffer(ipos+ipos_NOHot)=NOHot(i)*1.e6
         buffer(ipos+ipos_TnOHot)=phdisso2(i)*1.e6

             buffer(ipos+ipos_po1d)=Po1d(i)*1.e6
             buffer(ipos+ipos_no1d)=N_0*No1dnew(i)*1.e6
             buffer(ipos+ipos_uo1d)=Cj0*Uo1dnew(i)/1.e2
           enddo
           open(fid_temp,file=filetemps,
     &           form='unformatted',access='direct',recl=longrec,
     &           status='unknown')
           write(fid_temp,rec=nrectemps)(buffer(i),i=1,longbuf)
           close(fid_temp)
1000       format(' utc second:  ',i8,' tube number: ',i3)
           write(*,1000) int(temps),itube
       tempsort=tempsort+sortie

!------Output code for optical emissions
!-------Prepare file for output
        inquire(file='dir.output/emissions.dat',exist=exval)
        if(exval) then
          open(43,file='dir.output/emissions.dat'
     &     ,status='old',position='append')
        else
          open(43,file='dir.output/emissions.dat',status='new')
        endif
!-------MZ
          
!-------Output the values for the lines/bands 
        write(43,*) iyd,temps,dipangle,indlim,nen               !day, time, flux tube dip angle in degrees, and indices for alt. and en.
        write(43,*) (alt(m),                                    !altitude  
     &         No1dnew(m)*N_0,No1snew(m),Noii2P(m),Nn2Anew(m),  !metastable concentrations
     &         Po3p3P(m),Po3p5P(m),P1NG(m),PMein(m),P2PG(m),
     &         P1PG(m),m=1,indlim)

        write(43,*) (e(ien),fhemd(ien,1),ien=1,nen)             !energies for precipitation and topside diff. num. flux
        close(43)
!-------MZ

!-------Output upward and downward suprathermal fluxes during aurora
        if(temps.gt.tstartprec .and. temps.le.tstopprec) then
          !Create/open output file
          inquire(file='dir.output/ediffnumflux.dat',exist=exval)
          if(exval) then
            open(44,file='dir.output/ediffnumflux.dat'
     &       ,status='old',position='append')
          else
            open(44,file='dir.output/ediffnumflux.dat',status='new')
          endif

          !Loop over altitudes of interest
          !inc=floor(real(indlim,dp)/10)
       inc=1
          write(44,*) temps,ceiling(real(indlim,dp)/inc),nen
          do ien=1,nen
            write(44,*) e(ien)
          enddo
          do m=1,indlim,inc
            write(44,*) alt(indlim+1-m)
            do ien=1,nen
              write(44,*) fhemd(ien,m)/2./3.1416,                       
     &            fhemu(ien,m)/2./3.1416
            enddo
          enddo

          !Close the diff num flux file
          close(44);
        endif
!-------MZ

        endif

         if (duree.ge.tempsint) then
               tempsint=tempsint+postint
               flagatmos=.true.
         endif

         deltat=dt/t0
         deltat_2=.5d0*deltat
         deltat_4=.25d0*deltat

         duree=duree+dt
         temps=temps+dt
         if (temps.ge.86400.d0) then
           iyd=iyd+1
           temps=temps-86400.d0
         endif

         stl=mod(temps/3600.+longeo/15.,24.d0)
         iannee=iyd/1000
         call jour_mois(iyd,ijour,imois)
         iheure=int(temps/3600.)
        iminute=int(temps/60.-iheure*60.)
         seconde=temps-iheure*3600.-iminute*60.
         if (seconde.ge.60.) then
           iminute=iminute+1
           seconde=seconde-60.
         endif
         tu=iheure+iminute/60.
         call lec_indices(iannee,imois,ijour,tu,ap,f107)
         kp=ap2kp(ap(2))
         ikp=kp*3
    

         if (flgconv) then
           call convec(iyd,temps,kp,dlongeo,dlatgeo,
     &               dlonmag,dlatmag,dlonref,dt,pot,flgpot)
           longeo=dlongeo
           latgeo=dlatgeo
           vtrans=vpara*100./Ci0
         else
           call convec(iyd,temps,kp,dlongeo,dlatgeo,
     &                 dlonmag,dlatmag,dlonref,dx0,pot,flgpot)
           longeo=dlongeo
           latgeo=dlatgeo
           vtrans=0.
         endif

!-------a hard-coded dusk cell, 2-cell convection pattern
!  500 m/s
        tmotion=22.21*3600.
        omegaf=0.0045
        lonstart=75.

!  750 m/s
!        tmotion=22.92*3600.
!        omegaf=0.0067

!  1500 m/s
!        tmotion=24.1030*3600.
!        omegaf=0.0135

!  1000 m/s
!       tmotion=23.6294*3600.
!       omegaf=0.0090
!
!        if (temps .ge. tmotion .or. 
!     &          temps .lt. 1.05*3600.) then          !convect across polar cap
!                if (temps .lt. tmotion) then
!                  latgeo=67.+omegaf*(temps+86400.-tmotion)
!                else
!                  latgeo=67.+omegaf*(temps-tmotion)
!                endif
!                longeo=
!     &            mod(75.-(0.0041667-0.2054)*(temps-tempsini),360.)
!     &            mod(75.-(0.0041667)*(temps-tempsini),360.)
!     &            mod(lonstart-(0.0041667-0.0034)*(temps-tempsini),360.)
!
!                if (latgeo .gt. 90.) then
!                  latgeo=180.-latgeo
!                  longeo=longeo+180.
!                endif
!                if (longeo .gt. 360. .or. longeo .lt. 0.) then
!                  longeo=mod(longeo+360.,360.)
!                endif
!
!                if (temps .gt. 3690. .and. temps .lt. 3810.) then
!                  Enord=67.5
!                else
!                  if(Enord .lt. 50.) then
!                    Enord=10.+0.0444*incE
!                    incE=incE+1.
!                  else
!                    Enord=50.
!                  endif
!                  Enord=25.
!                endif
!        elseif (temps .lt. tmotion .and. 
!     &          temps .gt. tempsini) then            !stay at noon SLT
!                longeo=
!     &            mod(lonstart-0.0041667*(temps-tempsini),360.)
!
!                if (longeo .gt. 360. .or. longeo .lt. 0.) then
!                  longeo=mod(longeo+360.,360.)
!                endif
!                Enord=10.
!        else                                            !corotate
!                latgeo=67.
!                longeo=309.38
!                Enord=67.5
!!                cofo=1.
!        endif
!
!        call magfild(latgeo,longeo,300.,2001,Bmag,dipangle,or1)
!-------MZ


!------Overwrite the convection field assigned by convec.f and force it to DATCAR value
        Enord=Econvmag
        Eest=0.                                 !Only magnitude matters for frictional heating
!------MZ

!------Overwrite the activity indices set by lec_indices.f
        f107(1)=f107dat
        f107(2)=f107dat
        f107(3)=f107Adat
        ap=apdat
        kp=ap2kp(ap(2))
        ikp=kp*3
!------MZ


        if ((any(ieee_is_nan(nenew)))
     &  .or.(any(ieee_is_nan(Tenew)))
     &  .or.(any(ieee_is_nan(T1new)))) then
          write(stderr,*),'problem before calling atmos'
          goto 246
        endif
        call atmos(iyd,real(temps,sp),stl,alt,latgeo,longeo,jpreci,f107,
     &            ap,Nenew,Tenew,T1new,nx,kiappel,file_cond)

        if (vparaB.ne.0.) vtrans=vparaB*100./Ci0

         flagatmos=.false.

        if (flgne) stop'concentration negative'

!]]
![[     former timestep results memorization

      do i=1,nx

        N1old(i)=N1new(i)
        N2old(i)=N2new(i)
        N3old(i)=N3new(i)
        N4old(i)=N4new(i)
        N6old(i)=N6new(i)
        N5old(i)=N5new(i)
        Neold(i)=Nenew(i)
        Nnoold(i)=Nnonew(i)
            No1dold(i)=No1dnew(i)

        U1old(i)=U1new(i)
        U2old(i)=U2new(i)
        U3old(i)=U3new(i)
        Umold(i)=Umnew(i)
        Ueold(i)=Uenew(i)
        Unoold(i)=Unonew(i)
            Uo1dold(i)=Uo1dnew(i)

        T1old(i)=T1new(i)
        T1pold(i)=T1pnew(i)
        T1told(i)=T1tnew(i)
        T2old(i)=T2new(i)
        T2pold(i)=T2pnew(i)
        T2told(i)=T2tnew(i)
        T3old(i)=T3new(i)
        T3pold(i)=T3pnew(i)
        T3told(i)=T3tnew(i)
        Tmold(i)=Tmnew(i)
        Tmpold(i)=Tmpnew(i)
        Tmtold(i)=Tmtnew(i)

            Teold(i)=Tenew(i)
            Tepold(i)=Tepnew(i)
            Tetold(i)=Tetnew(i)
            qeold(i)=qenew(i)

        q1old(i)=q1new(i)
        q2old(i)=q2new(i)
        q3old(i)=q3new(i)


        Ter_1(i)=300./Tenew(i)/T_0

          enddo


      do i=2,nx
            Velnoc(i)=Cmi*.5*(Unonew(i)+Unonew(i-1))
        Velic(i)=.5*(U2new(i)+U2new(i-1))+vtrans
        Veljc(i)=Cji*.5*(U1new(i)+U1new(i-1))+vtrans
        Velec(i)=Cei*(Uenew(i)+Uenew(i-1))+vtrans
        Velmc(i)=Czi*.5*(Umnew(i)+Umnew(i-1))+vtrans
        Velnc(i)=Cni*.5*(U3new(i)+U3new(i-1))+vtrans
            velo1dc(i)=Cji*.5*(Uo1dnew(i)+Uo1dnew(i-1))
          enddo

          Velnoc(1)=Velnoc(2)
          rbc=min(1.,Nnonew(nx-1)/Nnonew(nx-2))
          rbc_2=sqrt(rbc)
          Nnonew(np)=rbc*Nnonew(nx)
          Unonew(np)=.5*(Unonew(nx)+Unonew(nx-1))
          if (Unonew(np).ge.0.) then
            Unonew(np)=max(coef_flnp/rbc_2*Unonew(nx),
     &               coef_fln/rbc*Unonew(np))
           else
            if (Unonew(nx).ge.0.) then
              Unonew(np)=coef_flnp/rbc_2*Unonew(nx)
            else
              Unonew(np)=max(.5*coef_flnp/rbc_2*Unonew(nx),
     &                 coef_fln/rbc*Unonew(np))
             endif
           endif
      Velnoc(np)=Cmi*(Unonew(np)+vtrans)
      Unonew(np)=2.*Unonew(np)-Unonew(nx)

          Velic(1)=Velic(2)
          rbc=min(1.,N2new(nx-1)/N2new(nx-2))
          rbc_2=sqrt(rbc)
          N2new(np)=rbc*N2new(nx)

          U2new(np)=.5*(U2new(nx)+U2new(nx-1))
          if (U2new(np).ge.0.) then
            U2new(np)=max(coef_flnp/rbc_2*U2new(nx),
     &               coef_fln/rbc*U2new(np))
           else
            if (U2new(nx).ge.0.) then
              U2new(np)=coef_flnp/rbc_2*U2new(nx)
            else
              U2new(np)=max(.5*coef_flnp/rbc_2*U2new(nx),
     &                 coef_fln/rbc*U2new(np))
             endif
           endif
!      Velic(np)=(U2new(np)+vtrans)
      Velic(np)=vtrans
      U2new(np)=2.*U2new(np)-U2new(nx)

          Veljc(1)=Veljc(2)
          rbc=min(1.,N1new(nx-1)/N1new(nx-2))
          rbc_2=sqrt(rbc)
          N1new(np)=rbc*N1new(nx)
          U1new(np)=.5*(U1new(nx)+U1new(nx-1))
          if (U1new(np).ge.0.) then
            U1new(np)=max(coef_flnp/rbc_2*U1new(nx),
     &               coef_fln/rbc*U1new(np))
           else
            if (U1new(nx).ge.0.) then
              U1new(np)=coef_flnp/rbc_2*U1new(nx)
            else
              U1new(np)=max(.5*coef_flnp/rbc_2*U1new(nx),
     &                 coef_fln/rbc*U1new(np))
             endif
           endif
      Veljc(np)=Cji*(U1new(np)+vtrans)
      U1new(np)=2.*U1new(np)-U1new(nx)

          Velec(1)=Velec(2)
          rbc=min(1.,Nenew(nx-1)/Nenew(nx-2))
          rbc_2=sqrt(rbc)
          Nenew(np)=rbc*Nenew(nx)
          Uenew(np)=.5*(Uenew(nx)+Uenew(nx-1))
          if (Uenew(np).ge.0.) then
            Uenew(np)=max(coef_flnp/rbc_2*Uenew(nx),
     &               coef_fln/rbc*Uenew(np))
           else
            if (Uenew(nx).ge.0.) then
              Uenew(np)=coef_flnp/rbc_2*Uenew(nx)
            else
              Uenew(np)=max(.5*coef_flnp/rbc_2*Uenew(nx),
     &                 coef_fln/rbc*Uenew(np))
             endif
           endif
      Velec(np)=Cei*(Uenew(np)+vtrans)
      Uenew(np)=2.*Uenew(np)-Uenew(nx)

          Velnc(1)=Velnc(2)
          rbc=min(1.,N3new(nx-1)/N3new(nx-2))
          rbc_2=sqrt(rbc)
          N3new(np)=rbc*N3new(nx)
          U3new(np)=.5*(U3new(nx)+U3new(nx-1))
          if (U3new(np).ge.0.) then
            U3new(np)=max(coef_flnp/rbc_2*U3new(nx),
     &               coef_fln/rbc*U3new(np))
           else
            if (U3new(nx).ge.0.) then
              U3new(np)=coef_flnp/rbc_2*U3new(nx)
            else
              U3new(np)=max(.5*coef_flnp/rbc_2*U3new(nx),
     &                 coef_fln/rbc*U3new(np))
             endif
           endif
      Velnc(np)=Cni*(U3new(np)+vtrans)
      U3new(np)=2.*U3new(np)-U3new(nx)

          Velmc(1)=Velmc(2)
          rbc=min(1.,Nmnew(nx-1)/Nmnew(nx-2))
          rbc_2=sqrt(rbc)
          Nmnew(np)=rbc*Nmnew(nx)
          Umnew(np)=.5*(Umnew(nx)+Umnew(nx-1))
          if (Umnew(np).ge.0.) then
            Umnew(np)=max(coef_flnp/rbc_2*Umnew(nx),
     &               coef_fln/rbc*Umnew(np))
           else
            if (Umnew(nx).ge.0.) then
              Umnew(np)=coef_flnp/rbc_2*Umnew(nx)
            else
              Umnew(np)=max(.5*coef_flnp/rbc_2*Umnew(nx),
     &                 coef_fln/rbc*Umnew(np))
             endif
           endif
      Velmc(np)=Cmi*(Umnew(np)+vtrans)
      Umnew(np)=2.*Umnew(np)-Umnew(nx)

          Velo1dc(1)=Velo1dc(2)
          rbc=min(1.,No1dnew(nx-1)/No1dnew(nx-2))
          rbc_2=sqrt(rbc)
          No1dnew(np)=rbc*No1dnew(nx)
          Uo1dnew(np)=.5*(Uo1dnew(nx)+Uo1dnew(nx-1))
          if (Uo1dnew(np).ge.0.) then
            Uo1dnew(np)=max(coef_flnp/rbc_2*Uo1dnew(nx),
     &               coef_fln/rbc*Uo1dnew(np))
           else
            if (Uo1dnew(nx).ge.0.) then
              Uo1dnew(np)=coef_flnp/rbc_2*Uo1dnew(nx)
            else
              Uo1dnew(np)=max(.5*coef_flnp/rbc_2*Uo1dnew(nx),
     &                 coef_fln/rbc*Uo1dnew(np))
             endif
           endif
      Velo1dc(np)=Cji*(Uo1dnew(np))
      Uo1dnew(np)=2.*Uo1dnew(np)-Uo1dnew(nx)

          do i=1,np
        Velnom(i)=.5*Velnoc(i)
            velo1dm(i)=.5*velo1dc(i)
        Velim(i)=.5*(Velic(i)+vtrans)
        Veljm(i)=.5*(Veljc(i)+vtrans)
        Velnm(i)=.5*(Velnc(i)+vtrans)
        Velmm(i)=.5*(Velmc(i)+vtrans)
    
        velie(i)=Velic(i)
        velje(i)=Veljc(i)
        velme(i)=Velmc(i)
        velne(i)=Velnc(i)
        velee(i)=Velec(i)
        Veliq(i)=Velie(i)
        Veljq(i)=Velje(i)
        Velnq(i)=Velne(i)
        Veleq(i)=Velee(i)
          enddo

!          velnoc(np)=max(velnoc(np),0.)
!          velic(np)=max(velic(np),0.)
!          veljc(np)=max(veljc(np),0.)
!          velec(np)=max(velec(np),0.)
!          velmc(np)=max(velmc(np),0.)
!          velnc(np)=max(velnc(np),0.)
!          velnom(np)=max(velnoc(np),0.)
!          velim(np)=max(velim(np),0.)
!          veljm(np)=max(veljm(np),0.)
!          velmm(np)=max(velmm(np),0.)
!          velnm(np)=max(velnm(np),0.)
!          velie(np)=max(velie(np),0.)
!          velje(np)=max(velje(np),0.)
!          velee(np)=max(velee(np),0.)
!          velme(np)=max(velme(np),0.)
!          velne(np)=max(velne(np),0.)
!          veliq(np)=max(veliq(np),0.)
!          veljq(np)=max(veljq(np),0.)
!          veleq(np)=max(veleq(np),0.)
!          velnq(np)=max(velnq(np),0.)

![[[[


!c      NO continuity equation

      call velocity(Velnoc,Ipos1,indno,xcoeffno*deltat)
    
       do i=1,indno
    
          Tr=(T_0*T1new(i)+Tn(i))/2.
      TrHot=(T_0*T1new(i)+TnOHot(i))/2.                    !MZ

          nujN2(i)=6.82e-10*Nn2(i)                              *t0
          nujO2(i)=6.64e-10*No2(i)                              *t0
          nujO (i)=3.67e-11*No(i)*sqrt(Tr)
     &             *(1.-.064*alog10(Tr))**2                     *t0
      nujOHot(i)=3.67e-11*NOHot(i)*sqrt(TrHot)                !MZ
     &             *(1.-.064*alog10(TrHot))**2                  *t0        !MZ

      nu_omega=(nujN2(i)+nujO2(i)+nujO(i))/omega(i)*omj
      coef_conv=1./(1.+nu_omega**2)
      Tperp(i)=cofterp*coef_conv*Vm_2(i)

      Tjr=T_0*T1new(i)
      TONO=(Tjr*30.+Tn(i)*16.)/46.+(16.*30./46.*Tperp(i))
      kONO=TONO/300.

      if ( TONO .le. 8200. .and. TONO .ge. 100.) then
        ak4 = 5.622e-13 - 6.094e-14*kONO + 5.74e-14*kONO**2
     &        - 1.399e-15*kONO**3 + 1.84e-17*kONO**4
      elseif( TONO .gt. 8200. .and. TONO .le. 30000.) then
        ak4 =-2.22e-11 + 1.64e-12*kONO - 6.7e-17*kONO**2
          elseif( TONO .gt. 30000.) then
            ak4 =-2.22e-11 + 1.64e-12*kONO - 6.7e-17*kONO**2
      endif

CCCCC  This is an attempt to fudge the neutral temperature so that NO
CCCCC  will form a simple layer, and the D-region will be preserved.

                Tn2(i)=Tn(nx)-(Tn(nx)+350)*exp(-0.0175*(alt(i)-alt(1)))

                if(Tn2(i).lt.0) Tn2(i)=0.

CCCCC  This is done simply so I can see Tn2 in the output file
                TnOHot(i)=Tn2(i)

CCCCC                                                                           !MZ


      Prno(i)=(5.e-16*N6new(i)*Nn2(i)+
     &    3.7e-11*No2(i)*N3new(i)
     &      +1.5e-11*Nn(i)*No2(i)/N_0*exp(-3600./Tn2(i))
     &           )                  *t0

      Lono(i)=(4.5e-10*N6new(i)*N_0+
     &    3.3e-10*N4new(i)*N_0+
     &    2.e-11*N3new(i)*N_0+
     &      ak4*N1new(i)*N_0
     &    +(8.3e-6+6.e-7)+1.5e-12*(Tn2(i)**.5)*Nn(i)
     &        )      *t0

        D3(i)= Prno(i)
        D7(i)=-Lono(i)-3.*Cmi*Unonew(i)*alt_geo_1(i)

        Nnonew(i)=Prno(i)/Lono(i)
          enddo
         do i=indno+1,nx
         factno=(Nh(i)/Nh(indno))
         Nnonew(i)=Nnonew(indno)*factno**30.
         enddo
        goto 147

    
!      call sources(Ipos1,indno,xcoeffno*deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,indno,xcoeffno*deltat,7,zero,D7,0.,0.)
          lbc=sqrt(Nnoold(2)/Nnoold(3))
          rbc=(Nh(indno+1)/Nh(indno))**30.
      lbc=1.
      call lcpfct(Nnoold,Nnonew,Ipos1,indno,
     &              lbc,0.,0.,Nnonew(np),.false.,1)
          do i=1,nx
            dexpnu=D7(i)*xcoeffno*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.d0+dexpnu
              dexpnu=1.d0
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Nnonew(i)=dexpnu*(Nnonew(i)-Nnoold(i)
     &             +D3(i)*xcoeffno*deltat)
     &             +expnu*Nnoold(i)
          enddo


          
         do i=1,indno
         Nnonew(i)=max(Nnonew(i),r_min)
         enddo

!]]]]
         do i=indno+1,nx
             factno=(Nh(i)/Nh(indno))
             Nnonew(i)=Nnonew(indno)*factno**30.
         enddo
         Nnonew(np)=min(1.,Nnonew(nx)/Nnonew(nx-1))*Nnonew(nx)

          if (any(ieee_is_nan(Nnonew))) then
            call cpu_time(tic)
        write(stderr,*),tic,'problem when calculating Tenew in loop 1'
            goto 246
          endif

!]]]

![[[
!    NO velocity equation

    
      call velocity(Velnom,Ipos1,indno,xcoeffno*deltat)

       do i=1,indno
    
!cc   calcul de nueddy
             Eddy=exp(-.05*(alt(i)-zeddy)**2)
             if (alt(i).le.zeddy) then
               Eddy=.5*Eddy+.5*exp(.07*(alt(i)-zeddy))
             endif
             Eddy=1.5e6*Eddy
!cc   fin de calcul de nueddy
        Ntot=No2(i)+N_0*Nnonew(i)+Nn2(i)+No(i)+Nn(i)+Nh(i)
             rhotot=32.*No2(i)+30.*N_0*Nnonew(i)+28.*Nn2(i)
     &            +16.*No(i)+14.*Nn(i)+Nh(i)
         mmoy=rhotot/Ntot
         Dmol=1.52e18*sqrt(1./mmoy+1./30.)*sqrt(Tn(i))/Ntot
         Diff=Dmol+Eddy
         nun=kb*Tn(i)/mm/Diff*t0
         mred=(Dmol+mmoy*amu/mm*Eddy)/Diff
    
         C2a(i)=-Cmi*Tn(i)/T_0
         D2a(i)=log(Nnonew(i)*Tn(i))
         D3(i)=-Cim*mred*G(i)+nun*Un(i)/Cm0

         D7(i)=-nun

       enddo
    
          D2al=(D2a(1)+D2a(2))/2.
          D2ar=.5*(log(Nnonew(indno+1)*Tn(indno+1))
     &           +log(Nnonew(indno)*Tn(indno)))
      call sources(Ipos1,indno,xcoeffno*deltat,2,C2a,DD2a,D2al,D2ar)
!      call sources(Ipos1,indno,xcoeffno*deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,indno,xcoeffno*deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=1.
      call lcpfct(Unoold,Unonew,Ipos1,indno,
     &              lbc,0.,0.,Unonew(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*xcoeffno*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Unonew(i)=dexpnu*(Unonew(i)-Unoold(i)
     &             +D3(i)*xcoeffno*deltat)
     &             +expnu*Unoold(i)
          enddo
          if (any(ieee_is_nan(Unonew))) then
            write(stderr,*) 'probleme lors du calcul de Unonew'
            goto 246
          endif

      do i=indno+1,nx
         Unonew(i)=Unonew(indno)
      enddo

!]]]

147    continue
![[       First-half step

!-------Load up production arrays for excited states of interest (and flip them around)
        do m=1,indlim
          !metastable
          Po1d(m)=prate(3,2,indlim+1-m)
          Po1s(m)=prate(3,3,indlim+1-m)
          Pn2a(m)=prate(1,3,indlim+1-m)+prate(1,4,indlim+1-m)
     &            +prate(1,5,indlim+1-m)+0.5*prate(1,7,indlim+1-m)
          Poii2P(m)=0.2*prate(3,jsg(3),indlim+1-m)

          !allowed
          P1NG(m)=0.11*prate(1,jsg(1),indlim+1-m)                   !BR from Rees 1989
          PMein(m)=0.39*prate(1,jsg(1),indlim+1-m)                      !BR from Rees 1989
          P2PG(m)=0.5*prate(1,7,indlim+1-m)
          P1PG(m)=0.5*prate(1,7,indlim+1-m)+prate(1,4,indlim+1-m)
     &                  +prate(1,5,indlim+1-m)                                  !1PG is complicated, this is approximate
          Po3p3P(m)=0.75*prate(3,5,indlim+1-m)
     &                  +1.e-2*prate(2,jsg(2),indlim+1-m)                             !O and O2
          Po3p5P(m)=prate(3,4,indlim+1-m)                                               !e on O
     &     +4.2/2.2*1.e-2*prate(2,jsg(2),indlim+1-m)                       !e on O2
     &     +N1new(m)*N_0*Nenew(m)*N_0*3.7e-12*(250./Tenew(m)/T_0)**.7         !rad. recomb.
        enddo
!-------MZ


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!O(1D) continuity equation
!        call velocity(velo1dc,Ipos1,Iposnp,deltat_2)

        do i=1,indlim
          D3(i)=Po1d(i)/N_0*t0+phdisso2(i)/N_0*t0+
     &      1.95e-7*(300./(Tenew(i)*T_0))**.7*N6new(i)*Nenew(i)*N_0*t0+
     &      .569*(Tenew(i)*T_0)**.5*(9329.+Tenew(i)*T_0)
     &      /(51183.+Tenew(i)*T_0)**3.*
     &      exp(-22756./(Tenew(i)*T_0))*Nenew(i)*No(i)*t0+
     &      1.215*No1snew(i)/N_0*t0+
     &      0.7*2.e-10*N3new(i)*No2(i)*t0
          D7(i)=-(
     &         5.85e-3+
     &         7.7e-3+
     &         2.0e-11*exp(107.8/Tn(i))*Nn2(i)+
     &         8.1e-10*(Tenew(i)*T_0/300.)**.5*Nenew(i)*N_0+
     &         2.9e-11*exp(67.5/Tn(i))*No2(i)+
     &         (3.73+.11965*sqrt(Tn(i))*(-6.5898e-4)*Tn(i))*10.**(-12.)
     &         )*t0
     
        enddo

        do i=1,indlim
          dexpnu=D7(i)*deltat_2
          if (abs(dexpnu).lt.1.e-7) then
            expnu=1.+dexpnu
            dexpnu=1.
          else
            expnu=exp(dexpnu)
            dexpnu=(expnu-1.)/dexpnu
          endif

          No1dnew(i)=dexpnu*(No1dnew(i)-No1dold(i)+D3(i)*deltat_2)
     &             +expnu*No1dold(i)
        enddo

!        do i=1,nx
!          No1dnew(i)=max(No1dnew(i),rmin)
!        enddo
!        No1dnew(np)=min(1.,No1dnew(nx-1)/No1dnew(nx-2))*No1dnew(nx)
!        No1dnew(nx)=min(1.,No1dnew(nx-2)/No1dnew(nx-3))*No1dnew(nx-1)
        
!        do i=1,nx
!                xno1d(i)=No1dnew(i)+Nliminf
!        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!O(1S) continuity equation

        do i=1,indlim
          !Calculate N2(A3 Sigu+) density via chemical equil.
          Nn2anew(i)=Pn2a(i)/(0.38+2.8e-11*No(i))

          !Production
          D3(i)=Po1s(i)+Po1sdisso2(i)+
!     &    0.12833*1.6e-7*(300./(Tenew(i)*T_0))**.5*
!     &    N6new(i)*N_0*Nenew(i)*N_0*+
     &    8.1e-10*(1150./(Tenew(i)*T_0))**(-1.47)*
     &    N6new(i)*N_0*Nenew(i)*N_0+
     &    2.5e-11*Nn(i)*N6new(i)*N_0+
!     &    0.36*2.8e-11*Nn2anew(i)*No(i)+               !this braching is recommended by dayglow literature
     &    0.2*2.8e-11*Nn2anew(i)*No(i)+               
     &    4.7e-33*(300./Tn(i))**2.*
     &    No(i)**3.*(Nn2(i)+No2(i))/(12.75*No2(i)+205.75*No(i))

          !Loss
          D7(i)=-(
     &          1.215+
     &          0.076+
     &          4.0e-12*exp(-865./Tn(i))*No2(i)+
     &          2e-14*No(i)
     &          )

          !Chemical equilibrium
          No1snew(i)=D3(i)/(-D7(i))
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!O+(2P) continuity equation
        do i=1,indlim
          D3(i)=Poii2P(i)
          D7(i)=-(4.8e-10*Nn2(i)+5.2e-11*No(i)+
     &            1.89e-7*sqrt(tenew(i)*T_0/300.)*Nenew(i)*N_0+
     &            0.173+0.047)

          Noii2P(i)=D3(i)/(-D7(i))
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      

!-------------------------------------!
!       equations de continuite       !
!-------------------------------------!

!    ion O+

      call velocity(Veljc,Ipos1,Iposnp,deltat_2)
          lbc=sqrt(N1old(2)/N1old(3))
      lbc=1.
      rbc=(N1old(nx)/N1old(nx-1))
      call lcpfct(N1old,N1new,Ipos1,Iposn,
     &              lbc,0.,0.,N1new(np),.false.,1)
    
!    ion H+

      call velocity(Velic,Ipos1,Iposnp,deltat_2)
          lbc=sqrt(N2old(2)/N2old(3))
      lbc=1.
      rbc=sqrt(N2old(nx-1)/N2old(nx-2))
      call lcpfct(N2old,N2new,Ipos1,Iposn,
     &              lbc,0.,0.,N2new(np),.false.,1)
    
!    ion N+

      call velocity(Velnc,Ipos1,Iposnp,deltat_2)
          lbc=sqrt(N3old(2)/N3old(3))
      lbc=1.
      rbc=(N3old(nx)/N3old(nx-1))
      call lcpfct(N3old,N3new,Ipos1,Iposn,
     &              lbc,0.,0.,N3new(np),.false.,1)
    
!    ion N2+

      call velocity(Velmc,Ipos1,Iposnp,deltat_2)
          lbc=sqrt(N4old(2)/N4old(3))
      lbc=1.
      rbc=(N4old(nx)/N4old(nx-1))
      call lcpfct(N4old,N4new,Ipos1,Iposn,
     &              lbc,0.,0.,N4new(np),.false.,1)
    
!    ion NO+

          lbc=sqrt(N5old(2)/N5old(3))
!      lbc=1.
      rbc=(N5old(nx)/N5old(nx-1))
      call lcpfct(N5old,N5new,Ipos1,Iposn,
     &              lbc,0.,0.,N5new(np),.false.,1)
    
!    ion O2+

          lbc=sqrt(N6old(2)/N6old(3))
      lbc=1.
      rbc=(N6old(nx)/N6old(nx-1))
      call lcpfct(N6old,N6new,Ipos1,Iposn,
     &              lbc,0.,0.,N6new(np),.false.,1)
    
      do i=1,nx

    
        Tjr=T_0*T1new(i)
        Tmr=T_0*Tmnew(i)
        Tr=(Tn(i)+Tmr)/2.
        Trh=sqrt(Tr)
        lTr=log10(Tr)

            nukO2(i)=0.449e-9*No2(i)                            *t0
            nukO (i)=0.258e-9*No (i)                            *t0
            nukN2 (i)=5.14e-11*Nn2 (i)
     &             *(1.-.069*alog10(Tr))**2*sqrt(Tr)            *t0
            nulN2(i)=0.413e-9*Nn2(i)                            *t0
            nulO (i)=0.231e-9*No (i)                            *t0
            nulO2 (i)=2.59e-11*No2 (i)
     &             *(1.-.073*alog10(Tr))**2*sqrt(Tr)            *t0
            numN2(i)=0.434e-9*Nn2(i)                            *t0
            numO (i)=0.244e-9*No (i)                            *t0
            numO2(i)=0.427e-9*No2 (i)                           *t0

        nu_omega=(nukN2(i)+nukO2(i)+nukO(i)
     &               +nulN2(i)+nulO2(i)+nulO(i)
     &               +numN2(i)+numO2(i)+numO(i))/omega(i)*omz
        coef_conv=1./(1.+nu_omega**2)
        Tperp(i)=cofterp*coef_conv*Vm_2(i)

    
        TjO2=(Tjr*32.+Tn(i)*16.)/48.+(16.*32./48.*Tperp(i))
        kjO2=TjO2/300.
        TjN2=(Tjr*28.+Tn(i)*16.)/44.+(16.*28./44.*Tperp(i))
        kjN2=TjN2/300.


        TN2O=(Tmnew(i)*T_0*16.+Tn(i)*28.)/44.+(16.*28./44.*Tperp(i))
        kN2O=TN2O/300.

        TONO=(Tjr*30.+Tn(i)*16.)/46.+(16.*30./46.*Tperp(i))
        kONO=TONO/300.

        TlO2=(32.*Tmr+16.*Tn(i))/48.+(16.*32./48.*Tperp(i))
        klO2=TlO2/300.

        TN2O2=(Tmr*32.+Tn(i)*28.)/60.+(32.*28./60.*Tperp(i))

!****************************************************************
! O+ + N2 -> NO+ + N  in an O buffer: *
!**************************************
      if (TjN2 .ge. 100. .and. TjN2 .le. 6200.) then
      ak1 = 1.248e-12 - 1.751e-13*kN2O
     &           - 5.101e-14*kN2O**2 + 1.345e-14*kN2O**3
     &           - 2.922e-16*kN2O**4
        elseif ( TjN2 .gt. 6200. .and. TjN2 .le. 22000.) then
      ak1 = -9.626e-11 + 6.994e-12*kN2O
     &           - 2.315e-14*kN2O**2
        elseif ( TjN2 .gt. 22000.) then
      ak1 = -9.626e-11 + 6.994e-12*kN2O
     &           - 2.315e-14*kN2O**2
        endif
!*****************************************
! O+ + N2 -> NO+ + N  in an N2 buffer:   *
!*****************************************
      if (TjN2 .ge. 100. .and. TjN2 .le. 5500.) then
      ak1 = 1.417e-12 - 3.223e-13*kN2O - 2.362e-14*kN2O**2
     &         + 1.247e-14*kN2O**3 - 3.030e-16*kN2O**4    
      elseif(TjN2 .gt. 5500. .and. TjN2 .le. 29000.) then
      ak1 =-4.74e-11 + 3.74e-12*kN2O + 2.80e-14*kN2O**2
     &         - 1.88e-16*kN2O**3
        elseif ( TjN2 .gt. 29000.) then
      ak1 =-4.74e-11 + 3.74e-12*kN2O + 2.80e-14*kN2O**2
     &         - 1.88e-16*kN2O**3
        endif


!****************************************
! O+ + O2 -> O2+ + O in an O buffer:    *
!****************************************
      if (TjO2 .ge. 100. .and. TjO2 .le. 6400.) then
       ak2 = 2.836e-11 - 7.521e-12*kjO2 + 1.039e-12*kjO2**2
     &            - 4.981e-14*kjO2**3 + 9.087e-16*kjO2**4
      elseif (TjO2 .gt. 6400. .and. TjO2 .le. 22000.) then
       ak2 =-3.42e-11 + 4.08e-12*kjO2 - 1.70e-14*kjO2**2
        elseif (TjO2 .gt. 22000.) then
           ak2 =-3.42e-11 + 4.08e-12*kjO2 - 1.70e-14*kjO2**2
        endif
!****************************************
! O+ + O2 -> O2+ + O in an N2 buffer:    *
!****************************************
      if (TjO2 .ge. 100. .and. TjO2 .le. 8400. ) then
       ak2 = 2.763e-11 - 6.733e-12*kjO2 + 8.383e-13*kjO2**2
     &            - 3.317e-14*kjO2**3 + 4.805e-16*kjO2**4
      elseif (TjO2 .gt. 8400. .and. TjO2 .le. 31000.) then
       ak2 =-2.57e-11 + 3.48e-12*kjO2 - 1.01e-14*kjO2**2
        elseif (TjO2 .gt. 31000.) then
           ak2 =-2.57e-11 + 3.48e-12*kjO2 - 1.01e-14*kjO2**2
      endif


!*************************
!  N2+ + O -> O+ + N2    *
!*************************
        if (TN2O.le.1500.) then
          ak3=1.e-11*kN2O**(-.23)
          ak5=1.4e-10*kN2O**(-.44)
        else
          ak3=3.6e-12*kN2O**.41
          ak5=5.2e-11*kN2O**.2
        endif


!****************************************
! O+ + NO -> NO+ + O in an  O buffer:   *
!****************************************
      if (TONO .gt. 100. .and. TONO .le. 6300. ) then
      ak4 = 5.974e-13 - 9.422e-14*kONO + 6.583e-14*kONO**2
     &           - 2.156e-15*kONO**3 + 3.957e-17*kONO**4
        elseif( TONO .gt. 6300. .and. TONO .le. 22000.) then
      ak4 =-1.557e-11 + 1.397e-12*kONO + 2.461e-15*kONO**2
        elseif( TONO .gt. 22000.) then
          ak4 =-1.557e-11 + 1.397e-12*kONO + 2.461e-15*kONO**2
      endif
!****************************************
! O+ + NO -> NO+ + O in an  N2 buffer:   *
!****************************************
      if (TONO .gt. 100. .and. TONO .le. 8200. ) then
      ak4 = 5.622e-13 - 6.094e-14*kONO + 5.74e-14*kONO**2
     &           - 1.399e-15*kONO**3 + 1.84e-17*kONO**4
        elseif( TONO .gt. 8200. .and. TONO .le. 30000.) then
      ak4 =-2.22e-11 + 1.64e-12*kONO - 6.7e-17*kONO**2
        elseif( TONO .gt. 30000.) then
          ak4 =-2.22e-11 + 1.64e-12*kONO - 6.7e-17*kONO**2
      endif
          
!****************************************
!  O+ + O2 -> O2+ + O in an N2 buffer:  *
!****************************************
        if( TlO2 .ge. 100. .and. TlO2 .le. 8400. ) then
          akl2 = 2.763e-11 - 6.733e-12*klO2 + 8.383e-13*klO2**2
     &               - 3.317e-14*klO2**3 + 4.805e-16*klO2**4
        elseif( TlO2 .gt. 8400. .and. TlO2 .le. 31000. ) then
          akl2 =-2.57e-11 + 3.48e-12*klO2 - 1.01e-14*klO2**2
            elseif( TlO2 .gt. 31000.) then
              akl2 =-2.57e-11 + 3.48e-12*klO2 - 1.01e-14*klO2**2
        endif

          
              y0  (1  )=N1old(i)
              y1  (1  )=(N1new(i)-N1old(i))/deltat_2+Po(i)/N_0*t0
              chim(1,1)=-(ak1*Nn2(i)+ak2*No2(i)+ak4*N_0*Nnonew(i)
     &                +2.5e-11*sqrt(Tn(i)+T1new(i)*T_0/16.
     &                +1.2e-8*(U1new(i)*Cj0)**2)*Nh(i))        *t0
     &                -3.*Cji*U1new(i)*alt_geo_1(i)
              chim(1,2)=2.2e-11*sqrt(T_0*T2new(i)+Tn(i)/16.
     &             +1.2e-8*(U2new(i)*Ci0)**2)*No(i)*t0
              chim(1,3)=(3.7e-11*No2(i)+5.e-13*No(i))*t0
              chim(1,4)=ak3*No(i)*t0
              chim(1,5)=0.d0
              chim(1,6)=0.d0
          
              y0  (2  )=N2old(i)
              y1  (2  )=(N2new(i)-N2old(i))/deltat_2+Ph(i)/N_0*t0
              chim(2,1)=2.5d-11*sqrt(Tn(i)+T1new(i)*T_0/16.+1.2e-8*
     &               (U1new(i)*Cj0)**2)*Nh(i)*t0
              chim(2,2)=-2.2d-11*sqrt(T_0*T2new(i)+Tn(i)/16.
     &              +1.2d-8*(U2new(i)*Ci0)**2)*No(i)        *t0
     &              -3.*U2new(i)*alt_geo_1(i)
              chim(2,3)=3.6d-12*Nh(i)*t0
              chim(2,4)=0.d0
              chim(2,5)=0.d0
              chim(2,6)=0.d0
          
              y0  (3  )=N3old(i)
              y1  (3  )=(N3new(i)-N3old(i))/deltat_2+0.21*Pn2(i)*t0/N_0
              chim(3,1)=0.d0
              chim(3,2)=0.d0
              chim(3,3)=-(2.6e-10*No2(i)+3.1e-10*No2(i)+3.7e-11*No2(i)
     &                +2.e-11*N_0*Nnonew(i)+3.6e-12*Nh(i)
     &                +5.e-13*No(i))*t0
     &              -3*Cni*U3new(i)*alt_geo_1(i)
              chim(3,4)=0.d0
              chim(3,5)=0.d0
              chim(3,6)=0.d0
          
              y0  (4  )=N4old(i)
              y1  (4  )=(N4new(i)-N4old(i))/deltat_2+Pn2(i)*0.79*t0/N_0          
        chim(4,1)=0.d0
              chim(4,2)=0.d0
              chim(4,3)=0.d0
             chim(4,4)=-(5.e-11*(300./TN2O2)*No2(i)+ak3*No(i)+ak5*No(i)
     &                +3.3e-10*N_0*Nnonew(i)
     &                +1.8e-7*Ter_1(i)**.39*Nenew(i)*N_0)*t0
     &              -3.*Czi*Umnew(i)*alt_geo_1(i)
              chim(4,5)=0.d0
              chim(4,6)=0.d0
          
              y0  (5  )=N5old(i)
              y1  (5  )=(N5new(i)-N5old(i))/deltat_2+6.e-7*Nnonew(i)*t0          
              chim(5,1)=(ak4*N_0*Nnonew(i)+ak1*Nn2(i))*t0
              chim(5,2)=0.d0
              chim(5,3)=(No2(i)*2.6e-10+2.e-11*N_0*Nnonew(i))*t0
              chim(5,4)=ak5*No(i)*t0
              chim(5,5)=-4.2e-7*Ter_1(i)**.85*Nenew(i)*N_0*t0
     &              -3.*Czi*Umnew(i)*alt_geo_1(i)
              chim(5,6)=(5.e-16*Nn2(i)+1.2e-10*Nn(i)
     &                 +4.5e-10*N_0*Nnonew(i))*t0
           
              y0  (6  )=N6old(i)
              y1  (6  )=(N6new(i)-N6old(i))/deltat_2+Po2(i)/N_0*t0
              chim(6,1)=akl2*No2(i)*t0
              chim(6,2)=0.d0
              chim(6,3)=3.1e-10*No2(i)*t0
              chim(6,4)=5.e-11*(300./TN2O2)*No2(i)*t0
              chim(6,5)=0.d0
              chim(6,6)=-(5.e-16*Nn2(i)+1.2e-10*Nn(i)
     &              +4.5e-10*N_0*Nnonew(i)
     &                +1.6e-7*Ter_1(i)**.55*Nenew(i)*N_0)*t0
     &              -3.*Czi*Umnew(i)*alt_geo_1(i)

        call solve_ODE(nb_ion,y0,y1,chim,mat,deltat_2)
    
        N1new(i)=max(y0(1),r_min)
        N2new(i)=max(y0(2),r_min)
        N3new(i)=max(y0(3),r_min)
        N4new(i)=max(y0(4),r_min)
        N5new(i)=max(y0(5),r_min)
        N6new(i)=max(y0(6),r_min)

      enddo
    
          N1new(np)=min(1.,N1new(nx-1)/N1new(nx-2))*N1new(nx)
          N2new(np)=min(1.,N2new(nx-1)/N2new(nx-2))*N2new(nx)
          N3new(np)=min(1.,N3new(nx-1)/N3new(nx-2))*N3new(nx)
          N4new(np)=min(1.,N4new(nx-1)/N4new(nx-2))*N4new(nx)
          N5new(np)=min(1.,N5new(nx-1)/N5new(nx-2))*N5new(nx)
          N6new(np)=min(1.,N6new(nx-1)/N6new(nx-2))*N6new(nx)

      if (any(ieee_is_nan(N1new))) then
        write(stderr,*) 'problem N1new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N2new))) then
        write(stderr,*) 'problem N2new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N3new))) then
        write(stderr,*) 'problem N3new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N4new))) then
        write(stderr,*) 'problem N4new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N5new))) then
        write(stderr,*) 'problem N5new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N6new))) then
        write(stderr,*) 'problem N6new dans la boucle 1'
        goto 246
      endif
    
    
      do i=1,nx
        Nenew(i)=N2new(i)+N1new(i)+N4new(i)
     &              +N6new(i)+N5new(i)+N3new(i)
             if (alt(i).le.500.) then
          xne(i)=Nenew(i)+Nliminf
          xn1(i)=N1new(i)+Nliminf
          xn2(i)=N2new(i)+Nliminf
          xn3(i)=N3new(i)+Nliminf
          xn4(i)=N4new(i)+Nliminf
          xn5(i)=N5new(i)+Nliminf
          xn6(i)=N6new(i)+Nliminf
          xnm(i)=Nmnew(i)+Nliminf
    
          xne_1(i)=Nenew(i)/xne(i)/xne(i)
          xn1_1(i)=N1new(i)/xn1(i)/xn1(i)
          xn2_1(i)=N2new(i)/xn2(i)/xn2(i)
          xn3_1(i)=N3new(i)/xn3(i)/xn3(i)
          xn4_1(i)=N4new(i)/xn4(i)/xn4(i)
          xn5_1(i)=N5new(i)/xn5(i)/xn5(i)
          xn6_1(i)=N6new(i)/xn6(i)/xn6(i)
          xnm_1(i)=Nmnew(i)/xnm(i)/xnm(i)
        else
          xne(i)=Nenew(i)
          xn1(i)=N1new(i)
          xn2(i)=N2new(i)
          xn3(i)=N3new(i)
          xn4(i)=N4new(i)
          xn5(i)=N5new(i)
          xn6(i)=N6new(i)
          xnm(i)=Nmnew(i)
    
          xne_1(i)=1./Nenew(i)
          xn1_1(i)=1./N1new(i)
          xn2_1(i)=1./N2new(i)
          xn3_1(i)=1./N3new(i)
          xn4_1(i)=1./N4new(i)
          xn5_1(i)=1./N5new(i)
          xn6_1(i)=1./N6new(i)
          xnm_1(i)=1./Nmnew(i)
        endif
      enddo
      xne(np)=Nenew(np)
      xn1(np)=N1new(np)
      xn2(np)=N2new(np)
      xn3(np)=N3new(np)
      xn4(np)=N4new(np)
      xn5(np)=N5new(np)
      xn6(np)=N6new(np)
      xnm(np)=Nmnew(np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!O(1D) momentum equation
        call velocity(velo1dm,Ipos1,Iposnp,deltat_2)

        do i=1,nx
          diffo1d=1./(No(i)/5.37e17/sqrt(Tn(i))+
     &          Nn2(i)/4.76e17/sqrt(Tn(i))+
     &          No2(i)/6.58e17/sqrt(Tn(i)))
          nuo1d=1.38e-17*Tn(i)/(16.*1.66e-27)/diffo1d*t0

          Tr=(Tn(i)+T_0*T1new(i))/2.
          nuo1dj(i)=N_0*N1new(i)*3.67e-11*sqrt(Tr)
     &               *(1.-0.064*alog10(Tr))**2                 *t0
     
          nuo1di(i)=N_0*N2new(i)*1/16*6.61e-11*
     &       sqrt(T_0*T2new(i))*(1.-0.047*
     &       alog10(T_0*T2new(i)))**2                          *t0
          
          nuo1dn(i)=N_0*N3new(i)*14./16.*4.42e-10              *t0
          
          nuo1dk(i)=N_0*N4new(i)*28./16.*2.58e-10              *t0
          
          nuo1dm(i)=N_0*N5new(i)*30./26.*2.44e-10              *t0
          
          nuo1dl(i)=N_0*N6new(i)*2.*2.31e-10                   *t0
          
          nuo1de(i)=N_0*Nenew(i)*8.9e-11*3.428e-5
     &          *(1.+5.7e-4*T_0*Tenew(i))
     &          *sqrt(T_0*Tenew(i))                            *t0


          D3(i)=-Cij*G(i)+
     &           nuo1d*Un(i)/Cj0+
     &           nuo1dj(i)*U1new(i)+nuo1di(i)*Cij*U2new(i)+
     &           (nuo1dk(i)+nuo1dl(i)+nuo1dm(i))*Czj*Umnew(i)+
     &          nuo1dn(i)*Cnj*U3new(i)+nuo1de(i)*Cej*Uenew(i)
          D7(i)=-(nuo1d+
     &          nuo1dj(i)+nuo1di(i)+nuo1dk(i)+nuo1dl(i)+
     &          nuo1dm(i)+nuo1dn(i)+nuo1de(i))
          
          C2a(i)=-Cji*Tn(i)/T_0*No1dnew(i)/xno1d(i)
          D2a(i)=log(xno1d(i)*Tn(i)/T_0)
        enddo
        
        D2al=(D2a(1)+D2a(2))/2.
        D2ar=(xno1d(np)*Tn(np)/T_0*xno1d(nx)*Tn(i)/T_0)
        call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
        call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
        
        lbc=1.
        call lcpfct(Uo1dold,Uo1dnew,ipos1,iposnp,
     &           lbc,0.,0.,Uo1dnew(np),.false.,0)

        do i=1,nx
          dexpnu=D7(i)*deltat_2
          if (abs(dexpnu).lt.1.e-7) then
            expnu=1.+dexpnu
            dexpnu=1.
          else
            expnu=exp(dexpnu)
            dexpnu=(expnu-1.)/dexpnu
          endif

          Uo1dnew(i)=dexpnu*(Uo1dnew(i)-Uo1dold(i)+D3(i)*deltat_2)
     &             +expnu*Uo1dold(i)
        Uo1dnew(i)=0.
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

![[[    H+ momentum equation resolution

      call velocity(Velim,Ipos1,Iposnp,deltat_2)
      
      call cpu_time(tic)
      write(stdout,*),tic,'H+ momentum equation resolution'
      
      do i=1,nx

          Tr1=T_0*T2new(i)
          Tr=(Tn(i)+Tr1)/2.

      T1_15(i)=T1new(i)**(-1.5)/nu_0
      T2_15(i)=T2new(i)**(-1.5)/nu_0
      T3_15(i)=T3new(i)**(-1.5)/nu_0
      Tm_15(i)=Tmnew(i)**(-1.5)/nu_0
      Te_15(i)=Tenew(i)**(-1.5)/nu_0
    
      thermacc(i)=thermodiff*qenew(i)/Tenew(i)**2.5

          nuiN2(i)=3.36e-9*Nn2(i)                               *t0

          nuiO2(i)=3.20e-9*No2(i)                               *t0
          nuiO (i)=6.61e-11*No (i)
     &             *(1.-.047*alog10(Tr1))**2*sqrt(Tr1)          *t0
          nuiOHot(i)=6.61e-11*NOHot(i)                        !MZ
     &             *(1.-.047*alog10(Tr1))**2*sqrt(Tr1)          *t0        !MZ
 
          nuiH (i)=2.65e-10*Nh (i)
     &             *(1.-.083*alog10(Tr))**2*sqrt(Tr)            *t0
          nuij (i) = 1.23*N1new(i)
     &              /((16.*T2new(i)+T1new(i))/17.)**1.5        /nu_0
          nuii (i) = .9*N2new(i)
     &              /(T2new(i))**1.5                 /nu_0
          nuin (i) = 1.23*N3new(i)
     &              /((14.*T2new(i)+T3new(i))/15.)**1.5     /nu_0
          nuik (i) = 1.25*N4new(i)
     &              /((28.*T2new(i)+Tmnew(i))/27.)**1.5        /nu_0
          nuim (i) = 1.25*N5new(i)
     &              /((30.*T2new(i)+Tmnew(i))/31.)**1.5      /nu_0
          nuil (i) = 1.25*N6new(i)
     &              /((32.*T2new(i)+Tmnew(i))/33.)**1.5        /nu_0
          nuie (i)=0.03    *Nenew(i)*Te_15(i)



        C2a(i)=-Tepnew(i)
        D2a(i)=log(Nenew(i)*Tepnew(i))
            C2b(i)=-T2pnew(i)*N2new(i)/xn2(i)
        D2b(i)=log(xn2(i)*T2pnew(i))
        D3(i)=-G(i)+thermacc(i)/mi/Ci0
     &            -3.*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*(T2pnew(i)-T2tnew(i))*alt_geo_1(i)
     &          +nuij(i)*Cji*U1new(i)
     &          +(nuik(i)*Czi+nuil(i)*Czi+nuim(i)*Czi)*Umnew(i)
     &          +nuin(i)*Cni*U3new(i)

        D3(i)=D3(i)+.6*q2new(i)*xn2_1(i)
     &               *(28.*nuiN2(i)/(Tn(i)/T_0+28.*T2new(i))
     &                  +32.*nuiO2(i)/(Tn(i)/T_0+32.*T2new(i))
     &                  +16.*nuiO(i)/(Tn(i)/T_0+16.*T2new(i))
     &                  +16.*nuiOHot(i)/(TnOHot(i)/T_0+16.*T2new(i))        !MZ
     &                  +    nuiH(i)/(Tn(i)/T_0+T2new(i))
     &                  +16.*nuij(i)/(T1new(i)+16.*T2new(i))
     &                  +28.*nuik(i)/(Tmnew(i)+28.*T2new(i))
     &                  +32.*nuil(i)/(Tmnew(i)+32.*T2new(i))
     &                  +30.*nuim(i)/(Tmnew(i)+30.*T2new(i))
     &                  +14.*nuin(i)/(T3new(i)+14.*T2new(i))
     &                  +5.46e-4*nuie(i)/(Tenew(i)+5.46e-4*T2new(i)))

        D3(i)=D3(i)-.6*q1new(i)*xn1_1(i)*Cji
     &              *nuij(i)/(T1new(i)+16.*T2new(i))

        D3(i)=D3(i)-.6*qenew(i)*xne_1(i)*Cei
     &              *nuie(i)/(Tenew(i)+5.46e-4*T2new(i))

        D3(i)=D3(i)-.6*q3new(i)*xn3_1(i)*Cni
     &              *nuin(i)/(T3new(i)+14.*T2new(i))
    
        D3(i)=D3(i)+(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+nuiH(i))        !MZ
     &                    *(Un(i)/Ci0)
     &                 -.6*nuiN2(i)/(Tn(i)/T_0+28.*T2new(i))
     &                    *(q_Nn2(i)*Niqi0)
     &                 -.6*nuiO2(i)/(Tn(i)/T_0+32.*T2new(i))
     &                    *(q_No2(i)*Niqi0)
     &                 -.6*nuiO(i)/(Tn(i)/T_0+16.*T2new(i))
     &                    *(q_No(i)*Niqi0)
     &                 -.6*nuiOHot(i)/(TnOHot(i)/T_0+16.*T2new(i))        !MZ
     &                    *(q_NOHot(i)*Niqi0)                    !MZ
     &                 -.6*nuiH(i)/(Tn(i)/T_0+T2new(i))
     &                    *(q_Nh(i)*Niqi0)
     &                 +nuie(i)*Cei*Uenew(i)

        D7(i)=-(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+nuiH(i)        !MZ
     &              +nuij(i)+nuik(i)+nuil(i)+nuim(i)
     &              +nuin(i)+nuie(i))

      enddo
          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
      D2br=.5*log(xn2(np)*T2pnew(np)*xn2(nx)*T2pnew(nx))

      D2br=log(xn2(nx)*T2pnew(nx))                        !MZ

      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
!      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=1.
          
      call cpu_time(tic)
      write(stdout,*),tic,'call lcpfct'
      call lcpfct(U2old,U2new,Ipos1,Iposn,
     &              lbc,0.,0.,U2new(np),.false.,0)
     
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif

            U2new(i)=dexpnu*(U2new(i)-U2old(i)+D3(i)*deltat_2)
     &             +expnu*U2old(i)
          enddo

          if (any(ieee_is_nan(U2new))) then
            call cpu_time(tic)
            write(stderr,*) tic,'problem calc U2new loop  1'
            goto 246
          endif


!]]]

![[[    O+ momentum equation resolution
        call cpu_time(tic)
      write(stdout,*),tic,'O+ momentum equation resolution'

      call velocity(Veljm,Ipos1,Iposnp,deltat_2)
      
    
      do i=1,nx

          Tr=(T_0*T1new(i)+Tn(i))/2.
      TrHot=(T_0*T1new(i)+TnOHot(i))/2.                    !MZ

          nujN2(i)=6.82e-10*Nn2(i)                              *t0
          nujO2(i)=6.64e-10*No2(i)                              *t0
          nujO (i)=3.67e-11*No(i)*sqrt(Tr)
     &             *(1.-.064*alog10(Tr))**2                     *t0
          nujOHot(i)=3.67e-11*NOHot(i)*sqrt(TrHot)
     &             *(1.-.064*alog10(TrHot))**2                  *t0        !MZ
          nujj(i)=.22*N1new(i)
     &           /(T1new(i))**1.5                /nu_0
          nuji(i)=0.077*N2new(i)
     &           /((T1new(i)+16.*T2new(i))/17.)**1.5        /nu_0
          nujn(i)=0.22*N3new(i)
     &           /((14.*T1new(i)+16.*T3new(i))/30.)**1.5     /nu_0
          nujk(i)=0.25*N4new(i)
     &           /((28.*T1new(i)+16.*Tmnew(i))/44.)**1.5    /nu_0
          nujm(i)=0.26*N5new(i)
     &           /((30.*T1new(i)+16.*Tmnew(i))/46.)**1.5    /nu_0
          nujl(i)=0.26*N6new(i)
     &             /((32.*T1new(i)+16.*Tmnew(i))/48.)**1.5    /nu_0
          nuje(i)=1.87e-3*Nenew(i)*Te_15(i)


        C2a(i)=-Cji*Tepnew(i)
        D2a(i)=log(Nenew(i)*Tepnew(i))
            C2b(i)=-Cji*T1pnew(i)*N1new(i)/xn1(i)
            D2b(i)=log(xn1(i)*T1pnew(i))
        D3(i)=-Cij*G(i)+thermacc(i)/mj/Cj0
     &            -3.*Cji*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*Cji*(T1pnew(i)-T1tnew(i))*alt_geo_1(i)
     &          +nuji(i)*Cij*U2new(i)

        D3(i)=D3(i)+(nujk(i)+nujl(i)+nujm(i))*Czj*Umnew(i)
     &               +nujn(i)*Cnj*U3new(i)

        D3(i)=D3(i)-.6*q2new(i)*xn2_1(i)*Cij*16.
     &            *nuji(i)/(16.*T2new(i)+T1new(i))

        D3(i)=D3(i)+.6*q1new(i)*xn1_1(i)
     &            *(28.*nujN2(i)/(16.*Tn(i)/T_0+28.*T1new(i))
     &                   +32.*nujO2(i)/(16.*Tn(i)/T_0+32.*T1new(i))
     &                   +16.*nujO(i)/(16.*Tn(i)/T_0+16.*T1new(i))
     &             +16.*nujOHot(i)/(16.*TnOHot(i)/T_0+16.*T1new(i))        !MZ
     &                   +    nuji(i)/(16.*T2new(i)+T1new(i))
     &                   +28.*nujk(i)/(16.*Tmnew(i)+28.*T1new(i))
     &                   +32.*nujl(i)/(16.*Tmnew(i)+32.*T1new(i))
     &                   +30.*nujm(i)/(16.*Tmnew(i)+30.*T1new(i))
     &                   +14.*nujn(i)/(16.*T3new(i)+14.*T1new(i))
     &                   +    nuje(i)/(16.*Tenew(i)+5.46e-4*T1new(i))
     &                        *5.46e-4)

        D3(i)=D3(i)-.6*qenew(i)*xne_1(i)*Cej*16.
     &            *nuje(i)/(16.*Tenew(i)+5.46e-4*T1new(i))

        D3(i)=D3(i)-.6*q3new(i)*xn3_1(i)*Cnj*16.
     &            *nujn(i)/(16.*T3new(i)+14.*T1new(i))

        D3(i)=D3(i)+(nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i))*(Un(i)/Cj0)    !MZ
     &                 -.6*nujN2(i)/(16.*Tn(i)/T_0+28.*T1new(i))
     &                  *(16.*q_Nn2(i)*Njqj0)
     &                 -.6*nujO2(i)/(16.*Tn(i)/T_0+32.*T1new(i))
     &                  *(16.*q_No2(i)*Njqj0)
     &                 -.6*nujO(i)/(16.*Tn(i)/T_0+16.*T1new(i))
     &                  *(16.*q_No(i)*Njqj0)
     &                 -.6*nujOHot(i)/(16.*TnOHot(i)/T_0+16.*T1new(i))        !MZ
     &                  *(16.*q_NOHot(i)*Njqj0)                    !MZ
     &                 +nuje(i)*Cje*Uenew(i)

        D7(i)=-(nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i)            !MZ
     &              +nuji(i)+nujn(i)
     &              +nuje(i)+nujk(i)+nujl(i)+nujm(i))

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
      D2br=.5*log(xn1(np)*T1pnew(np)*xn1(nx)*T1pnew(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
!      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=1.
      
      call cpu_time(tic)
      write(stdout,*),tic,'H+ momentum LCPFCT'
          
      call lcpfct(U1old,U1new,Ipos1,Iposn,
     &              lbc,0.,0.,U1new(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            U1new(i)=dexpnu*(U1new(i)-U1old(i)+D3(i)*deltat_2)
     &             +expnu*U1old(i)

          enddo
!      U1new(nx)=max(U1new(nx),-1000./Cj0)
!      U1new(nx)=max(U1new(nx),0.)

          if (any(ieee_is_nan(U1new))) then
          call cpu_time(tic)
         write(stderr,*) tic,'problem calc U1new loop 1'
            goto 246
          endif




!]]]

![[[    heavy ions momentum equation resolution
      call cpu_time(tic)
      write(stdout,*),tic,'Heavy Ions equation resolution'
      call velocity(Velmm,Ipos1,Iposnp,deltat_2)
    
      do i=1,nx

          Tr=(T_0*Tmnew(i)+Tn(i))/2.

          nukH(i)=0.074e-9*Nh(i)                                *t0
          nukO2(i)=0.449e-9*No2(i)                              *t0
          nukO (i)=0.258e-9*No (i)                              *t0
          nukOHot(i)=0.258e-9*NOHot(i)                          *t0        !MZ
          nukN2 (i)=5.14e-11*Nn2 (i)
     &             *(1.-.069*alog10(Tr))**2*sqrt(Tr)            *t0

          nulH(i)=0.065e-9*Nh(i)                                *t0
          nulN2(i)=0.413e-9*Nn2(i)                              *t0
          nulO (i)=0.231e-9*No (i)                              *t0
          nulOHot(i)=0.231e-9*NOHot(i)                          *t0        !MZ
          nulO2 (i)=2.59e-11*No2 (i)
     &             *(1.-.073*alog10(Tr))**2*sqrt(Tr)            *t0

          numH(i)=0.069e-9*Nh(i)                                *t0
          numN2(i)=0.434e-9*Nn2(i)                              *t0
          numO (i)=0.244e-9*No (i)                              *t0
          numOHot(i)=0.244e-9*NOHot(i)                          *t0        !MZ
          numO2(i)=0.427e-9*No2 (i)                             *t0

          nukj(i)=0.15*N1new(i)
     &            / ((16.*Tmnew(i)+28.*T1new(i))/44.)**1.5    /nu_0
          nuki(i)=0.045*N2new(i)
     &            / ((Tmnew(i)+28.*T2new(i))/29.)**1.5        /nu_0
          nukn(i)=0.14*N3new(i)
     &            / ((14.*Tmnew(i)+28.*T3new(i))/42.)**1.5    /nu_0
          nukk(i)=0.17*N4new(i)
     &            /(Tmnew(i))**1.5                /nu_0
          nukm(i)=0.17*N5new(i)
     &            /(Tmnew(i))**1.5                /nu_0
          nukl(i)=0.18*N6new(i)
     &            /(Tmnew(i))**1.5                /nu_0
          nuke(i)=1.07e-3*Nenew(i)/(Tenew(i))**1.5        /nu_0

          nulj(i)=0.13*N1new(i)
     &            / ((16.*Tmnew(i)+32.*T1new(i))/48.)**1.5    /nu_0
          nuli(i)=0.039*N2new(i)
     &            / ((Tmnew(i)+32.*T2new(i))/33.)**1.5        /nu_0
          nuln(i)=0.12*N3new(i)
     &            / ((14.*Tmnew(i)+32.*T3new(i))/46)**1.5    /nu_0
          nulk(i)=0.15*N4new(i)/(Tmnew(i))**1.5        /nu_0
          nulm(i)=0.16*N5new(i)/(Tmnew(i))**1.5        /nu_0
          null(i)=0.16*N6new(i)/(Tmnew(i))**1.5        /nu_0
          nule(i)=9.347e-4*Nenew(i)/(Tenew(i))**1.5    /nu_0
          numj(i)=0.14*N1new(i)
     &            / ((16.*Tmnew(i)+30.*T1new(i))/46.)**1.5    /nu_0
          numi(i)=0.042*N2new(i)
     &            /((Tmnew(i)+30.*T2new(i))/31.)**1.5        /nu_0
          numn(i)=0.13*N3new(i)
     &            /((14.*Tmnew(i)+30.*T3new(i))/44.)**1.5    /nu_0
          numk(i)=0.16*N4new(i)/(Tmnew(i))**1.5             /nu_0
          numm(i)=0.16*N5new(i)/(Tmnew(i))**1.5            /nu_0
          numl(i)=0.17*N6new(i)/(Tmnew(i))**1.5            /nu_0
          nume(i)=9.97e-4*Nenew(i)/(Tenew(i))**1.5        /nu_0

      C2a(i)=-Czi*Tepnew(i)
      D2a(i)=log(Nenew(i)*Tepnew(i))

          C2b(i)=-Czi*Tmpnew(i)*Nmnew(i)/xnm(i)
          D2b(i)= log(xnm(i)*Tmpnew(i))

          D3(i)=-Ciz*G(i)+thermacc(i)/mz/Cz0
     &            -3.*Czi*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*Czi*(Tmpnew(i)-Tmtnew(i))*alt_geo_1(i)
     &         +(nukj(i)+nulj(i)+numj(i))*Cjz/3.*U1new(i)
     &         +(nuki(i)+nuli(i)+numi(i))*Ciz/3.*U2new(i)
     &         +(nukn(i)+nuln(i)+numn(i))*Cnz/3.*U3new(i)

      D3(i)=D3(i)-.2*q2new(i)*xn2_1(i)*Ciz
     &              *(28.*nuki(i)/(28.*T2new(i)+Tmnew(i))
     &                 +32.*nuli(i)/(32.*T2new(i)+Tmnew(i))
     &                 +30.*numi(i)/(30.*T2new(i)+Tmnew(i)))

      D3(i)=D3(i)-.2*q1new(i)*xn1_1(i)*Cjz
     &              *(28.*nukj(i)/(28.*T1new(i)+16.*Tmnew(i))
     &                 +32.*nulj(i)/(32.*T1new(i)+16.*Tmnew(i))
     &                 +30.*numj(i)/(30.*T1new(i)+16.*Tmnew(i)))

      D3(i)=D3(i)-.2*qenew(i)*xne_1(i)*Cez
     &              *(28.*nuke(i)/(28.*Tenew(i)+5.46e-4*Tmnew(i))
     &                 +32.*nule(i)/(32.*Tenew(i)+5.46e-4*Tmnew(i))
     &                 +30.*nume(i)/(30.*Tenew(i)+5.46e-4*Tmnew(i)))

      D3(i)=D3(i)-.2*q3new(i)*xn3_1(i)*Cnz
     &              *(28.*nukn(i)/(28.*T3new(i)+14.*Tmnew(i))
     &                 +32.*nuln(i)/(32.*T3new(i)+14.*Tmnew(i))
     &                 +30.*numn(i)/(30.*T3new(i)+14.*Tmnew(i)))

      D3(i)=D3(i)+(nukN2(i)+nukO2(i)+nukO(i)+nukOHot(i)            !MZ
     &               +nulN2(i)+nulO2(i)+nulO(i)+nulOHot(i)            !MZ
     &        +numN2(i)+numO2(i)+numO(i)+numOHot(i))*Un(i)/Cz0/3.        !MZ
     &           -.2*
     &               (nukN2(i)/(28.*Tn(i)/T_0+28.*Tmnew(i))
     &           *(28.*q_Nn2(i)*Nkqk0)
     &               +nukO2(i)/(28.*Tn(i)/T_0+32.*Tmnew(i))
     &           *(28.*q_No2(i)*Nkqk0)
     &               +nukO(i)/(28.*Tn(i)/T_0+16.*Tmnew(i))
     &           *(28.*q_No(i)*Nkqk0)
     &        +nukOHot(i)/(28.*TnOHot(i)/T_0+16.*Tmnew(i))            !MZ
     &        *(28.*q_NOHot(i)*Nkqk0))                    !MZ
     &           -.2*
     &               (nulN2(i)/(32.*Tn(i)/T_0+28.*Tmnew(i))
     &           *(32.*q_Nn2(i)*Nlql0)
     &               +nulO2(i)/(32.*Tn(i)/T_0+32.*Tmnew(i))
     &           *(32.*q_No2(i)*Nlql0)
     &               +nulO(i)/(32.*Tn(i)/T_0+16.*Tmnew(i))
     &           *(32.*q_No(i)*Nlql0)
     &        +nulOHot(i)/(32.*TnOHot(i)/T_0+16.*Tmnew(i))            !MZ
     &        *(32.*q_NOHot(i)*Nlql0))                    !MZ
     &           -.2*
     &               (numN2(i)/(30.*Tn(i)/T_0+28.*Tmnew(i))
     &           *(30.*q_Nn2(i)*Nmqm0)
     &                +numO2(i)/(30.*Tn(i)/T_0+32.*Tmnew(i))
     &           *(30.*q_No2(i)*Nmqm0)
     &                +numO(i)/(30.*Tn(i)/T_0+16.*Tmnew(i))
     &           *(30.*q_No(i)*Nmqm0)
     &        +numOHot(i)/(30.*TnOHot(i)/T_0+16.*Tmnew(i))        !MZ
     &           *(30.*q_NOHot(i)*Nmqm0))                    !MZ
     &           +Uenew(i)*Cez*(nuke(i)+nule(i)+nume(i))/3.

          D7(i)=-(nukN2(i)+nukO2(i)+nukO(i)+nukOHot(i)                !MZ
     &            +nuki(i)+nuke(i)+nukj(i)+nukn(i)
     &            +nulN2(i)+nulO2(i)+nulO(i)+nulOHot(i)                !MZ
     &            +nuli(i)+nule(i)+nulj(i)+nuln(i)
     &            +numN2(i)+numO2(i)+numO(i)+numOHot(i)                !MZ
     &            +numi(i)+nume(i)+numj(i)+numn(i))/3.
    

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
      D2br=.5*log(xnm(np)*Tmpnew(np)*xnm(nx)*Tmpnew(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
!      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U1old,extra,nx)
      rbc=1.
      call lcpfct(Umold,Umnew,Ipos1,Iposn,
     &              lbc,0.,0.,Umnew(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Umnew(i)=dexpnu*(Umnew(i)-Umold(i)+D3(i)*deltat_2)
     &             +expnu*Umold(i)
          enddo

          if (any(ieee_is_nan(Umnew))) then
          call cpu_time(tic)
       write(stderr,*) tic,'Heavy Ions problem calc  Umnew loop 1'
            goto 246
          endif


![[[[   Boundaries conditions

!]]]

![[[    N+ momentum equation resolution
      call cpu_time(tic)
      write(stdout,*),tic,'N+ momentum equation resolution'
    
      call velocity(Velnm,Ipos1,Iposnp,deltat_2)
      do i=1,nx

          nunH (i)=0.145e-9*Nh(i)                               *t0
          nunN2(i)=0.747e-9*Nn2(i)                              *t0
          nunO (i)=0.442e-9*No (i)                              *t0
          nunOHot(i)=0.442e-9*NOHot(i)                          *t0        !MZ
          nunO2(i)=0.725e-9*No2 (i)                             *t0
          nunj(i)=0.25*N1new(i)
     &           / ((16.*T3new(i)+14.*T1new(i))/30.)**1.5    /nu_0
          nuni(i)=0.088*N2new(i)
     &            / (( T3new(i)+14.*T2new(i))/15.)**1.5        /nu_0
          nunn(i)=0.24*N3new(i)/ (T3new(i))**1.5        /nu_0
          nunk(i)=0.28*N4new(i)
     &            / ((28.*T3new(i)+14.*Tmnew(i))/42.)**1.5    /nu_0
          nunm(i)=0.28*N5new(i)
     &            / ((30.*T3new(i)+14.*Tmnew(i))/44.)**1.5    /nu_0
          nunl(i)=0.28*N6new(i)
     &            / ((32.*T3new(i)+14.*Tmnew(i))/46.)**1.5    /nu_0
          nune(i)=2.136e-3*Nenew(i)*Te_15(i)


        C2a(i)=-Cni*Tepnew(i)
        D2a(i)=log(Nenew(i)*Tepnew(i))
        C2b(i)=-Cni*T3pnew(i)*N3new(i)/xn3(i)
            D2b(i)=log(xn3(i)*T3pnew(i))


        D3(i)=-Cin*G(i)+thermacc(i)/mn/Cn0
     &            -3.*Cni*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*Cni*(T3pnew(i)-T3tnew(i))*alt_geo_1(i)
     &          +nuni(i)*Cin*U2new(i)+nunj(i)*Cjn*U1new(i)
     &          +(nunk(i)+nunl(i)+nunm(i))*Czn*Umnew(i)
    
      D3(i)=D3(i)-.6*q2new(i)*xn2_1(i)*Cin*14.
     &            *nuni(i)/(14.*T2new(i)+T3new(i))


      D3(i)=D3(i)-.6*q1new(i)*xn1_1(i)*Cjn*14.
     &            *nunj(i)/(14.*T1new(i)+16.*T3new(i))


      D3(i)=D3(i)-.6*qenew(i)*xne_1(i)*Cen*14.
     &            *nune(i)/(14.*Tenew(i)+5.46e-4*T3new(i))


      D3(i)=D3(i)+.6*q3new(i)*xn3_1(i)
     &            *(28.*nunN2(i)/(14.*Tn(i)/T_0+28.*T3new(i))
     &                  +32.*nunO2(i)/(14.*Tn(i)/T_0+32.*T3new(i))
     &                  +16.*nunO(i)/(14.*Tn(i)/T_0+16.*T3new(i))
     &                +16.*nunOHot(i)/(14.*TnOHot(i)/T_0+16.*T3new(i))    !MZ
     &                  +    nuni(i)/(14.*T2new(i)+T3new(i))
     &                  +28.*nunk(i)/(14.*Tmnew(i)+28.*T3new(i))
     &                  +32.*nunl(i)/(14.*Tmnew(i)+32.*T3new(i))
     &                  +30.*nunm(i)/(14.*Tmnew(i)+30.*T3new(i))
     &                  +    nune(i)/(14.*Tenew(i)+5.46e-4*T3new(i))
     &                       *5.46e-4)

      D3(i)=D3(i)+(nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i))*Un(i)/Cn0        !MZ
     &           -.6*nunN2(i)/(14.*Tn(i)/T_0+28.*T3new(i))
     &           *(14.*q_Nn2(i)*Nnqn0)
     &           -.6*nunO2(i)/(14.*Tn(i)/T_0+32.*T3new(i))
     &           *(14.*q_No2(i)*Nnqn0)
     &           -.6*nunO(i)/(14.*Tn(i)/T_0+16.*T3new(i))
     &           *(14.*q_No(i)*Nnqn0)
     &           -.6*nunOHot(i)/(14.*TnOHot(i)/T_0+16.*T3new(i))        !MZ
     &           *(14.*q_NOHot(i)*Nnqn0)                    !MZ
     &           +nune(i)*Cen*Uenew(i)

        D7(i)=-(nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i)            !MZ
     &             +nuni(i)+nunj(i)+nune(i)
     &             +nunk(i)+nunl(i)+nunm(i))

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
      D2br=.5*log(xn3(np)*T3pnew(np)*xn3(nx)*T3pnew(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
!      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U3old,extra,nx)
      rbc=1.
      call lcpfct(U3old,U3new,Ipos1,Iposn,
     &              lbc,0.,0.,U3new(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            U3new(i)=dexpnu*(U3new(i)-U3old(i)+D3(i)*deltat_2)
     &             +expnu*U3old(i)
          enddo
!      U3new(nx)=max(U3new(nx),0./Cn0)

!    do i=1,nx
!      U3new(i)=U1new(i)
!    enddo

          if (any(ieee_is_nan(U3new))) then
          call cpu_time(tic)
       write(stderr,*) tic,'N+ mom. problem calc U3new in loop 1'
            goto 246
          endif


!]]]


![[[    Velocities corrections

      do i=1,nx

      Uenew(i)=(Ci0*N2new(i)*U2new(i)
     &           +Cj0*N1new(i)*U1new(i)
     &           +Cn0*N3new(i)*U3new(i)
     &           +Cz0*(N4new(i)+N6new(i)+N5new(i))*Umnew(i)
     &            -JJ(i)/N_0-0.*Jes(i)/N_0)
     &               /Nenew(i)/Ce0

      enddo




      call velocity(Veliq,Ipos1,Iposnp,deltat_2)



![[[    H+ heat flow equation resolution
      call cpu_time(tic)
      write(stdout,*),tic,'H+ heatflow equation resolution'
      do i=1,nx

        C2a(i)=-2.2*q2new(i)
        D2a(i)=U2new(i)
        C2b(i)=-(11./18.*T2pnew(i)+8./9.*T2tnew(i))*N2new(i)
        D2b(i)=T2pnew(i)
        C2c(i)=-(17./9.*T2pnew(i)-8./9.*T2tnew(i))*N2new(i)
        D2c(i)=T2tnew(i)
        C2d(i)=4./9.*(T2pnew(i)-T2tnew(i))**2
        D2d(i)=N2new(i)
      D3(i)=alt_geo_1(i)*N2new(i)*(T2pnew(i)-T2tnew(i))
     &                           *(T2pnew(i)-4.*T2tnew(i))
     &             -N2new(i)*T2new(i)*U2new(i)
     &                     *(nuij(i)*18.5/17.
     &                      +nuik(i)*30.5/29.
     &                      +nuil(i)*34.5/33.
     &                      +nuim(i)*32.5/31.
     &                      +nuin(i)*16.5/15.
     &                      +nuie(i)*2.5
     &                      +nuiN2(i)*30.5/29.
     &                      +nuiO2(i)*34.5/33.
     &                      +nuiO(i)*18.5/17.
     &                      +nuiOHot(i)*18.5/17.                !MZ
     &                      +nuiH(i)*3.5/2.)
     &            + 2.5*N2new(i)*T2new(i)*U2new(i)*(
     &            nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+nuiH(i)        !MZ
     &                      +nuij(i)
     &                      +nuik(i)+nuil(i)+nuim(i)
     &                      +nuin(i)
     &                      +nuie(i))

            D3(i)=D3(i)+N2new(i)*T2new(i)*U1new(i)*Cji*
     &                           (nuij(i)*18.5/17.
     &                           -2.5*nuij(i))

            D3(i)=D3(i)+N2new(i)*T2new(i)*Umnew(i)*Czi*
     &                           (nuik(i)*30.5/29.
     &                           +nuil(i)*34.5/33.+nuim(i)*32.5/31.
     &                           -2.5*(nuik(i)+nuil(i)+nuim(i))     )

            D3(i)=D3(i)+N2new(i)*T2new(i)*U3new(i)*Cni*
     &                           (nuin(i)*16.5/15.
     &                           -2.5*nuin(i)         )

            D3(i)=D3(i)+N2new(i)*q1new(i)*Cji*xn1_1(i)*
     &            nuij(i)*(.0612+1.5*T2new(i)/(T1new(i)+16.*T2new(i)))

            D3(i)=D3(i)+N2new(i)*qenew(i)*Cei*xne_1(i)*
     &            nuie(i)*1.5*(T2new(i)/(Tenew(i)+5.46e-4*T2new(i))-1)

            D3(i)=D3(i)+N2new(i)*q3new(i)*Cni*xn3_1(i)*
     &            nuin(i)*(.068+1.5*T2new(i)/(T3new(i)+14.*T2new(i)))

            D3(i)=D3(i)+N2new(i)*(
     &                            nuiN2(i)*(1.07*q_Nn2(i)*Niqi0
     &                                     +1.052*T2new(i)*Un(i)/Ci0)
     &                           +nuiO2(i)*(1.084*q_No2(i)*Niqi0
     &                                     +1.045*T2new(i)*Un(i)/Ci0)
     &                           +nuiO(i) *(.98*q_No(i)*Niqi0
     &                                     +1.088*T2new(i)*Un(i)/Ci0)
     &                           +nuiOHot(i) *(.98*q_NOHot(i)*Niqi0        !MZ
     &                                     +1.088*T2new(i)*Un(i)/Ci0)        !MZ
     &                           +nuiH(i) *(-.075*q_Nh(i)*Niqi0
     &                                     +1.75*T2new(i)*Un(i)/Ci0 )
     &                           )
     &          -2.5*T2new(i)*N2new(i)*        (
     &                   (nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+nuiH(i))    !MZ
     &                   *(Un(i)/Ci0)
     &                   -.6*(nuiN2(i)/(Tn(i)/T_0+28.*T2new(i))
     &                       *(q_Nn2(i)*Niqi0)
     &                       +nuiO2(i)/(Tn(i)/T_0+32.*T2new(i))
     &                       *(q_No2(i)*Niqi0)
     &                       +nuiO (i)/(Tn(i)/T_0+16.*T2new(i))
     &                       *(q_No(i)*Niqi0)
     &                       +nuiOHot(i)/(TnOHot(i)/T_0+16.*T2new(i))        !MZ
     &                       *(q_NOHot(i)*Niqi0)                !MZ
     &                       +nuiH (i)/(Tn(i)/T_0+    T2new(i))
     &                       *(q_Nh(i)*Niqi0))  )

            D7(i)=
     &            (
     &              nuiN2(i)*.1795+nuiO2(i)*.1824
     &              +nuiO(i)*.1612+nuiOHot(i)*.1612-nuiH(i)*.725        !MZ
     &              -nuii(i)*.8+nuij(i)*.1612+nuik(i)*.1795
     &              +nuil(i)*.1824+nuim(i)*.1811+nuin(i)*.1547
     &              -nuie(i)*3.
     &            )
     &            -1.5*T2new(i)*(
     &                           nuiN2(i)*28./(Tn(i)/T_0+28.*T2new(i))
     &                          +nuiO2(i)*32./(Tn(i)/T_0+32.*T2new(i))
     &                          +nuiO (i)*16./(Tn(i)/T_0+16.*T2new(i))
     &                +nuiOHot(i)*16./(TnOHot(i)/T_0+16.*T2new(i))    !MZ
     &                          +nuiH (i)/(Tn(i)/T_0+T2new(i))
     &                          +nuij(i)*16./(T1new(i)+16.*T2new(i))
     &                          +nuik(i)*28./(Tmnew(i)+28.*T2new(i))
     &                          +nuil(i)*32./(Tmnew(i)+32.*T2new(i))
     &                          +nuim(i)*30./(Tmnew(i)+30.*T2new(i))
     &                          +nuin(i)*14./(T3new(i)+14.*T2new(i))
     &                          +nuie(i)*5.46e-4/
     &                               (Tenew(i)+5.46e-4*T2new(i))     )
     &            -4.2*U2new(i)*alt_geo_1(i)

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
          D2cl=(D2c(1)+D2c(2))/2.
          D2cr=ylimd(Radn(np),D2c,extra,nx)
          D2dl=(D2d(1)+D2d(2))/2.
          D2dr=ylimd(Radn(np),D2d,extra,nx)
    
      D2ar=.5*(U2new(np)+U2new(nx))
      D2br=.5*(T2pnew(np)+T2pnew(nx))
      D2cr=.5*(T2tnew(np)+T2tnew(nx))
      D2dr=sqrt(N2new(np)*N2new(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
      call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
      call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
!      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),q2old,extra,nx)
      rbc=1.
      call lcpfct(q2old,q2new,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            q2new(i)=dexpnu*(q2new(i)-q2old(i)+D3(i)*deltat_2)
     &             +expnu*q2old(i)
          enddo
          q2new(nx)=max(0.,q2new(nx))
          q2new(np)=q2new(nx)

          if (any(ieee_is_nan(q2new))) then
          call cpu_time(tic)
            write(stderr,*) 'H+ heatflow problem calc q2new loop 1'
            goto 246
          endif


!]]]


      call velocity(Veljq,Ipos1,Iposnp,deltat_2)

![[[    O+ heat flow equation resolution
       call cpu_time(tic)
       write(stdout,*),tic,'O+ heat flow equation resolution'
      do i=1,nx

        C2a(i)=-2.2*Cji*q1new(i)
        D2a(i)=U1new(i)
        C2b(i)=-Cji*(11./18.*T1pnew(i)+8./9.*T1tnew(i))*N1new(i)
        D2b(i)=T1pnew(i)
        C2c(i)=-Cji*(17./9.*T1pnew(i)-8./9.*T1tnew(i))*N1new(i)
        D2c(i)=T1tnew(i)
        C2d(i)=4./9.*Cji*(T1pnew(i)-T1tnew(i))**2
        D2d(i)=N1new(i)
      D3(i)=Cji*alt_geo_1(i)*N1new(i)*(T1pnew(i)-T1tnew(i))
     &                           *(T1pnew(i)-4.*T1tnew(i))
     &             +N1new(i)*T1new(i)*U2new(i)*Cij*
     &                   nuji(i)*(2.412-2.5)

          D3(i)=D3(i)-N1new(i)*T1new(i)*U1new(i)*          (
     &                         nujN2(i)*1.545+nujO2(i)*1.5
     &                  +nujO(i)*1.75+nujOHot(i)*1.75+nuji(i)*2.412    !MZ
     &                        +nujn(i)*1.8   +nujk(i)*1.545
     &                        +nujl(i)*1.5   +nujm(i)*1.522
     &                        +nuje(i)*2.5                 )
     &               + 2.5*N1new(i)*T1new(i)*U1new(i)*(nujN2(i)
     &                              +nujO2(i)+nujO(i)+nujOHot(i)        !MZ
     &                              +nuji(i)+nujn(i)
     &                              +nujk(i)+nujl(i)+nujm(i)+nuje(i))

          D3(i)=D3(i)+N1new(i)*T1new(i)*Umnew(i)*Czj*                (
     &                         nujk(i)*1.545+nujl(i)*1.5+nujm(i)*1.522
     &                        - 2.5*(nujk(i)+nujl(i)+nujm(i))        )

          D3(i)=D3(i)-N1new(i)*T1new(i)*U3new(i)*Cnj*.7*nujn(i)

          D3(i)=D3(i)-N1new(i)*q2new(i)*Cij*xn2_1(i)*nuji(i)*
     &                   (1.262- 24.*T1new(i)/(16.*T2new(i)+T1new(i)))

           D3(i)=D3(i)-N1new(i)*qenew(i)*Cej*xne_1(i)*nuje(i)*
     &             (1.5-24.*T1new(i)/(16.*Tenew(i)+5.46e-4*T1new(i)))

           D3(i)=D3(i)-N1new(i)*q3new(i)*Cnj*xn3_1(i)*nujn(i)*
     &             (.128-24.*T1new(i)/(16.*T3new(i)+14.*T1new(i)))

           D3(i)=D3(i)+N1new(i)*(
     &                           nujN2(i)*(.139*q_Nn2(i)*Njqj0
     &                                    +1.545*T1new(i)*Un(i)/Cj0)
     &                          +nujO2(i)*(.2*q_No2(i)*Njqj0
     &                                    +1.5*T1new(i)*Un(i)/Cj0)
     &                          +nujO(i)* (-.075*q_No(i)*Njqj0
     &                                    +1.75*T1new(i)*Un(i)/Cj0)
     &                          +nujOHot(i)* (-.075*q_NOHot(i)*Njqj0        !MZ
     &                                    +1.75*T1new(i)*Un(i)/Cj0)        !MZ
     &                          )
     &           -2.5*T1new(i)*N1new(i)*(
     &            (nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i))*(Un(i)/Cj0)    !MZ
     &               -9.6*Njqj0*
     &                  (q_Nn2(i)*nujN2(i)/(16.*Tn(i)/T_0+28.*T1new(i))
     &                  +q_No2(i)*nujO2(i)/(16.*Tn(i)/T_0+32.*T1new(i))
     &                  +q_No(i)*nujO (i)/(16.*Tn(i)/T_0+16.*T1new(i))
     &          +q_NOHot(i)*nujOHot(i)/(16.*TnOHot(i)/T_0+16.*T1new(i))    !MZ
     &                  )               )

         D7(i)=-(nujN2(i)*.34+nujO2(i)*.27+nujO(i)*.725+nujOHot(i)*.725    !MZ
     &           +nuji(i)*2.66+nujk(i)*.34+nujl(i)*.27
     &           +nujm(i)*.3+nujn(i)*.835+nujj(i)*.8+nuje(i)*3.)
     &          -1.5*T1new(i)*(
     &                        nujN2(i)*28./(16.*Tn(i)/T_0+28.*T1new(i))
     &                       +nujO2(i)*32./(16.*Tn(i)/T_0+32.*T1new(i))
     &                       +nujO (i)*16./(16.*Tn(i)/T_0+16.*T1new(i))
     &            +nujOHot(i)*16./(16.*TnOHot(i)/T_0+16.*T1new(i))    !MZ
     &                       +nuji(i)/(16.*T2new(i)+T1new(i))
     &                       +nujk(i)*28./(16.*Tmnew(i)+28.*T1new(i))
     &                       +nujl(i)*32./(16.*Tmnew(i)+32.*T1new(i))
     &                       +nujm(i)*30./(16.*Tmnew(i)+30.*T1new(i))
     &                       +nujn(i)*14./(16.*T3new(i)+14.*T1new(i))
     &                       +nuje(i)*5.46e-4
     &                               /(16.*Tenew(i)+5.46e-4*T1new(i)))
     &          -4.2*Cji*U1new(i)*alt_geo_1(i)

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
          D2cl=(D2c(1)+D2c(2))/2.
          D2cr=ylimd(Radn(np),D2c,extra,nx)
          D2dl=(D2d(1)+D2d(2))/2.
          D2dr=ylimd(Radn(np),D2d,extra,nx)
    
      D2ar=.5*(U1new(np)+U1new(nx))
      D2br=.5*(T1pnew(np)+T1pnew(nx))
      D2cr=.5*(T1tnew(np)+T1tnew(nx))
      D2dr=sqrt(N1new(np)*N1new(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
      call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
      call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
!      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),q1old,extra,nx)
      rbc=1.
      call lcpfct(q1old,q1new,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            q1new(i)=dexpnu*(q1new(i)-q1old(i)+D3(i)*deltat_2)
     &             +expnu*q1old(i)
          enddo
          q1new(nx)=max(0.,q1new(nx))
          q1new(np)=q1new(nx)

      if (any(ieee_is_nan(q1new))) then
          call cpu_time(tic)
          write(stderr,*) tic,'L3557: problem calc q1new in loop 1'
            goto 246
      endif


!]]]
       call cpu_time(tic)
       write(stdout,*),tic,'N+ heat flow'
      call velocity(Velnq,Ipos1,Iposnp,deltat_2)

![[[    N+ heat flow equation resolution

      do i=1,nx


       C2a(i)=-2.2*Cni*q3new(i)
       D2a(i)=U3new(i)
       C2b(i)=-Cni*(11./18.*T3pnew(i)+8./9.*T3tnew(i))*N3new(i)
       D2b(i)=T3pnew(i)
       C2c(i)=-Cni*(17./9.*T3pnew(i)-8./9.*T3tnew(i))*N3new(i)
       D2c(i)=T3tnew(i)
       C2d(i)=4./9.*Cni*(T3pnew(i)-T3tnew(i))**2
       D2d(i)=N3new(i)

      D3(i)=Cni*alt_geo_1(i)*N3new(i)*(T3pnew(i)-T3tnew(i))
     &                           *(T3pnew(i)-4.*T3tnew(i))
     &            -.1*N3new(i)*U2new(i)*Cin*T3new(i)*nuni(i)
     &          -.8*N3new(i)*U1new(i)*Cjn*T3new(i)*nunj(i)
     &          -   N3new(i)*Umnew(i)*Czn*T3new(i)
     &                        *(nunk(i)+1.04*nunl(i)+1.02*nunm(i))

          D3(i)=D3(i)+N3new(i)*T3new(i)*U3new(i)*               (
     &                   nunN2(i)+nunO2(i)*1.04+nunO(i)*.8
     &            +nunOHot(i)*.8                        !MZ
     &                        +nuni(i)*.1+nunj(i)*.8+nunk(i)
     &                        +nunl(i)*1.04+nunm(i)*1.02        )

          D3(i)=D3(i)+N3new(i)*q2new(i)*Cin*xn2_1(i)*nuni(i)*
     &               (21.*T3new(i)/(14.*T2new(i)+T3new(i))-1.232)

          D3(i)=D3(i)+N3new(i)*q1new(i)*Cjn*xn1_1(i)*nunj(i)*
     &               (21.*T3new(i)/(14.*T1new(i)+16.*T3new(i))-.028)

          D3(i)=D3(i)+N3new(i)*qenew(i)*Cen*xne_1(i)*nune(i)*
     &               (21.*T3new(i)/(14.*Tenew(i)+5.46e-4*T3new(i))-1.5)

          D3(i)=D3(i)+N3new(i)*  (
     &            nunN2(i)* (.1*q_Nn2(i)*Nnqn0
     &                      +1.5*T3new(i)*Un(i)/Cn0)
     &            +nunO2(i)*(.115*q_No2(i)*Nnqn0
     &                      +1.456*T3new(i)*Un(i)/Cn0)
     &            +nunO(i)* (-.028*q_No(i)*Nnqn0
     &                      +1.7*T3new(i)*Un(i)/Cn0)
     &            +nunOHot(i)* (-.028*q_NOHot(i)*Nnqn0                !MZ
     &                      +1.7*T3new(i)*Un(i)/Cn0)                !MZ
     &                           )
     &          -2.5*T3new(i)*N3new(i)*(
     &            (nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i))*(Un(i)/Cn0)    !MZ
     &           -8.4*Nnqn0*(
     &                  nunN2(i)*q_Nn2(i)/(14.*Tn(i)/T_0+28.*T3new(i))
     &                 +nunO2(i)*q_No2(i)/(14.*Tn(i)/T_0+32.*T3new(i))
     &                 +nunO (i)*q_No(i)/(14.*Tn(i)/T_0 +16.*T3new(i))
     &        +nunOHot(i)*q_NOHot(i)/(14.*TnOHot(i)/T_0 +16.*T3new(i)))    !MZ
     &                                 )


          D7(i)=-( nunN2(i)*.27+nunO2(i)*.2+nunO(i)*.62
     &            +nunOHot(i)*.62+nuni(i)*2.62                !MZ
     &            +nunk(i)*.27+nunl(i)*.2 +nunm(i)*.23
     &            +nunj(i)*.62+nunn(i)*.8 +nune(i)*3.)
     &            -1.5*T3new(i)*(
     &                 nunN2(i)*28./(14.*Tn(i)/T_0+28.*T3new(i))
     &                +nunO2(i)*32./(14.*Tn(i)/T_0+32.*T3new(i))
     &                +nunO (i)*16./(14.*Tn(i)/T_0+16.*T3new(i))
     &                +nunOHot(i)*16./(14.*TnOHot(i)/T_0+16.*T3new(i))        !MZ
     &                +nuni(i)/(14.*T2new(i)+T3new(i))
     &                +nunk(i)*28./(14.*Tmnew(i)+28.*T3new(i))
     &                +nunl(i)*32./(14.*Tmnew(i)+32.*T3new(i))
     &                +nunm(i)*30./(14.*Tmnew(i)+30.*T3new(i))
     &                +nunj(i)*16./(14.*T1new(i)+16.*T3new(i))
     &                +nune(i)*5.46e-4/(14.*Tenew(i)+5.46e-4*T3new(i)))
     &          -4.2*Cni*U3new(i)*alt_geo_1(i)

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
          D2cl=(D2c(1)+D2c(2))/2.
          D2cr=ylimd(Radn(np),D2c,extra,nx)
          D2dl=(D2d(1)+D2d(2))/2.
          D2dr=ylimd(Radn(np),D2d,extra,nx)
    
      D2ar=.5*(U3new(np)+U3new(nx))
      D2br=.5*(T3pnew(np)+T3pnew(nx))
      D2cr=.5*(T3tnew(np)+T3tnew(nx))
      D2dr=sqrt(N3new(np)*N3new(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
      call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
      call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
!      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),q3old,extra,nx)
      rbc=1.
      call lcpfct(q3old,q3new,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            q3new(i)=dexpnu*(q3new(i)-q3old(i)+D3(i)*deltat_2)
     &             +expnu*q3old(i)
          enddo

![[[   Boundaries conditions

          q3new(nx)=max(0.,q3new(nx))
          q3new(np)=q3new(nx)

      if (any(ieee_is_nan(q3new))) then
        write(stderr,*) 'problem when calculating q3new in loop 1'
        goto 246
      endif


!]]]

    
    

![[[    Electron energy and heat flow equation resolution (1)
      call cpu_time(tic)
      write(stdout,*),tic,'e- energy & heat flow eqn resolution (1)'
    
      do i=1,nx

        Ter=T_0*Tenew(i)
        Terh=sqrt(Ter)

        nueN2(i)=2.33e-11*Nn2(i)
     &               *exp(-1.21e-4*Ter)*Ter            *t0
        nueO2(i)=1.82e-10*No2(i)
     &             *(1.+3.6e-2*Terh)*Terh            *t0
        nueO (i)=8.90e-11*No (i)
     &             *(1.+5.7e-4*Ter)*Terh            *t0
        nueOHot(i)=8.90e-11*NOHot(i)                    !MZ
     &             *(1.+5.7e-4*Ter)*Terh            *t0        !MZ
        nueH (i)=4.5e-9*Nh(i)
     &             *exp(-1.35e-4*Ter)*Terh            *t0
        nuei (i)=54.5*N2new(i)*Te_15(i)
        nuej (i)=54.5*N1new(i)*Te_15(i)
        nuee (i)=54.5*Nenew(i)*Te_15(i)/sqrt(2.)
            nuek (i)=54.5*N4new(i)*Te_15(i)
            nuel (i)=54.5*N6new(i)*Te_15(i)
            nuem (i)=54.5*N5new(i)*Te_15(i)
            nuen (i)=54.5*N3new(i)*Te_15(i)

            C2ea(i)=Cei*Tenew(i)/3.
            D2ea(i)=Uenew(i)
            C2eb(i)=-.667*Cei*xne_1(i)
            D2eb(i)=qenew(i)

        !compression term
            C2qa(i)=-11./5.*Cei*qenew(i)
            D2qa(i)=Uenew(i)
            
        !Temperature gradient (thermal conduction)
        C2qb(i)=-2.5*Cei*Nenew(i)*Tenew(i)
            D2qb(i)=Tenew(i)

          qen=Tenew(i)*(Un(i)/Ce0-Uenew(i))

          D3q(i)=-1.5*Nenew(i)*Tenew(i)*(nuei(i)*Cie*U2new(i)
     &                                 +nuej(i)*Cje*U1new(i)
     &                                 +(nuek(i)+nuel(i)+nuem(i))
     &                                 *Cze*Umnew(i)
     &                                 +nuen(i)*Cne*U3new(i))

          D3q(i)=D3q(i)+Nenew(i)*q2new(i)*Cie*xn2_1(i)*nuei(i)*
     &                         (6.55e-4+8.19e-4*Tenew(i)
     &                         /(5.46e-4*T2new(i)+Tenew(i)))

          D3q(i)=D3q(i)+Nenew(i)*q1new(i)*Cje*xn1_1(i)*nuej(i)*
     &                         (4.09e-5+8.19e-4*Tenew(i)
     &                         /(5.46e-4*T1new(i)+16.*Tenew(i)))

          D3q(i)=D3q(i)+Nenew(i)*q3new(i)*Cne*xn3_1(i)*nuen(i)*
     &                         (4.68e-5+8.19e-4*Tenew(i)
     &                         /(5.46e-4*T3new(i)+14.*Tenew(i)))

          D3q(i)=D3q(i)+Nenew(i)* (
     &            nueN2(i)*(1.2*q_Nn2(i)*Neqe0+qen)
     &           +nueO2(i)*(1.2*q_No2(i)*Neqe0+qen)
     &           +nueO(i)* (1.2*q_No(i)*Neqe0 +qen)
     &           +nueOHot(i)* (1.2*q_NOHot(i)*Neqe0 +qen)            !MZ
     &           +nueH(i)* (1.2*q_Nh(i)*Neqe0 +qen)
     &          +1.5*(nuei(i)+nuej(i)+nuek(i)+nuel(i)+nuem(i)+nuen(i))
     &              *Tenew(i)*Uenew(i)
     &                          )
     &          -2.5*Nenew(i)*Tenew(i)*(
     &            (nueN2(i)+nueO2(i)+nueO(i)+nueOHot(i)+nueH(i))        !MZ
     &              *(Un(i)/Ce0-Uenew(i))
     &              -.6*nueN2(i)/(5.46e-4*Tn(i)/T_0+28.*Tenew(i))
     &                 *(5.46e-4*q_Nn2(i)*Neqe0)
     &              -.6*nueO2(i)/(5.46e-4*Tn(i)/T_0+32.*Tenew(i))
     &                 *(5.46e-4*q_No2(i)*Neqe0)
     &              -.6*nueO (i)/(5.46e-4*Tn(i)/T_0+16.*Tenew(i))
     &                 *(5.46e-4*q_No(i)*Neqe0)
     &              -.6*nueOHot(i)/(5.46e-4*TnOHot(i)/T_0+16.*Tenew(i))        !MZ
     &                 *(5.46e-4*q_NOHot(i)*Neqe0)                !MZ
     &              - .6*nueH (i)/(5.46e-4*Tn(i)/T_0+Tenew(i))
     &                 *(5.46e-4*q_Nh(i)*Neqe0)
     &                                 )

          D7q(i)=.2*(nueN2(i)+nueO2(i)+nueO(i)+nueOHot(i)+nueH(i)        !MZ
     &             +nuei(i) +nuej(i) +nuek(i)+nuel(i)
     &             +nuem(i) +nuen(i) -nuee(i)*4       )
     &          -1.5*Tenew(i)*(
     &                 nueN2(i)*28./(5.46e-4*Tn(i)/T_0+28.*Tenew(i))
     &                +nueO2(i)*32./(5.46e-4*Tn(i)/T_0+32.*Tenew(i))
     &                +nueO (i)*16./(5.46e-4*Tn(i)/T_0+16.*Tenew(i))
     &            +nueOHot(i)*16./(5.46e-4*TnOHot(i)/T_0+16.*Tenew(i))    !MZ
     &                +nueH (i)/(5.46e-4*Tn(i)/T_0+Tenew(i))
     &                +nuei(i)/(5.46e-4*T2new(i)+Tenew(i))
     &                +nuej(i)*16./(5.46e-4*T1new(i)+16.*Tenew(i))
     &                +nuek(i)*28./(5.46e-4*Tmnew(i)+28.*Tenew(i))
     &                +nuel(i)*32./(5.46e-4*Tmnew(i)+32.*Tenew(i))
     &                +nuem(i)*30./(5.46e-4*Tmnew(i)+30.*Tenew(i))
     &                +nuen(i)*14./(5.46e-4*T3new(i)+14.*Tenew(i)))
     &          -4.2*Cei*Uenew(i)*alt_geo_1(i)                    !geometric term

      if (abs(Ter-Tn(i)).lt.1.e-2) then
         dTen=1.e-2
         if (Ter.ge.Tn(i)) dTen=-dTen
         Ter=Tn(i)-dTen
      else
         dTen=Tn(i)-Ter
      endif

      dTen_1=dTen/Ter/Tn(i)

          Lenrot=(4.6e-26*Nn2(i)+1.1e-25*No2(i))
     &            *dTen/sqrt(Ter)

          if (Ter.lt.10000.) then
            f=exp(2.2e-3*(Ter-1800.))
            f=1.06e4+7.51e3*(f-1.)/(f+1.)
          else
            f=18110.
          endif

          gg=3300+(1.233-2.056e-4*(Ter-4000.))*(Ter-1000.)


            LeN2vib=4.78e-24*Nn2(i)
     &              *exp(f*(Ter-2000.)/2000./Ter)
     &           *(exp(gg*dTen_1)-1.)

          h=3300.-839.*sin(1.91e-4*(Ter-2700.))

          LeO2vib=8.3e-25*No2(i)
     &            *exp(h*(Ter-700.)/700./Ter)
     &         *(exp(2770.*dTen_1)-1.)

      dd=2.4e4+(.3-1.947e-5*(Ter-4000))*(Ter-1500.)

          LeOexc=2.5e-24*No(i)
     &      *exp(dd*(Ter-3000.)/Ter/3000.)
     &      *(exp(22713.*dTen_1)-1.)


          LeOfin=0.
          do j=1,3
            if (j.le.2) then
              Ex=exp(-Eq(j)/Ter)
              Dx(j)=exp(-Eq(j)/Tn(i))
            else
              Ex=exp(-Eq(3)/Ter-Eq(1)/Tn(i))
              Dx(j)=exp(-Eq(2)/Tn(i))
            endif
            F=epsq(j)*(Dx(j)-Ex)-5.91e-9*dTen
     &          *((1.+Bq(j))*Dx(j)+Ex*(Eq(j)/Ter+1.+Bq(j)))
            LeOfin=LeOfin+Aq(j)*Cq(j)*F*(Ter)**(Bq(j)-.5)
          enddo
    
          Z=5.+3.*Dx(1)+Dx(2)
          LeOfin=1.38e-17*No(i)/Z*LeOfin

          LeOfin=Sq(3)*Ter**.6*exp(Eq(3)/Tn(i))
     &              *(exp(Eq(3)*dTen_1)-1.)
    
          do j=1,2
        LeOfin=LeOfin+Sq(j)*(exp(Eq(j)*dTen_1)-1.)
      enddo
    
          Z=5.+3.*Dx(1)+Dx(2)
          LeOfin=No(i)/Z*LeOfin

      Len(i)=(LeOfin+Lenrot+LeN2vib+LeO2vib+LeOexc)*N_0*t0/P_0
      Len(i)=Len(i)/dTen

      D7e(i)=-.667*Len(i)*T_0
!      Len(i)=Heat(i)*N_0*t0/P_0+Len(i)*Tn(i)
      D3e(i)=.667*(Heat(i)/xne(i)*t0/P_0+Len(i)*Tn(i))

!    fin de modif

          D3e(i)= D3e(i)+T2new(i)*1.09e-3*nuei(i)
     &          +T1new(i)*6.8e-5*nuej(i)
     &          +Tmnew(i)*1.09e-3*(nuek(i)/28.+nuel(i)/32.+nuem(i)/30.)
     &          +T3new(i)*7.8e-5*nuen(i)


        nu_omega=(nueN2(i)+nueO2(i)+nueO(i)+nueOHot(i)+
     &            nueH(i))/omega(i)*ome
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Ce0**2
    
          D3e(i)=D3e(i)+1.09e-3*(nueN2(i)/28.+nueO2(i)/32.+nueO(i)/16.
     &                 +nueH(i))*Tn(i)/T_0
     &            +1.09e-3*nueOHot(i)/16.*TnOHot(i)/T_0            !MZ
     &               + .667*(nueN2(i)*((Uenew(i)-Un(i)/Ce0)**2
     &               + terme_Joule)
     &               + nueO2(i)*((Uenew(i)-Un(i)/Ce0)**2
     &               + terme_Joule)
     &               + nueO(i)*((Uenew(i)-Un(i)/Ce0)**2
     &               + terme_Joule)
     &               + nueOHot(i)*((Uenew(i)-Un(i)/Ce0)**2            !MZ
     &               + terme_Joule)                        !MZ
     &               + nueH(i)*((Uenew(i)-Un(i)/Ce0)**2
     &               + terme_Joule)
     &               + nuei(i)*(Cie*U2new(i)-Uenew(i))**2
     &               + nuej(i)*(Cje*U1new(i)-Uenew(i))**2
     &               + nuek(i)*(Cze*Umnew(i)-Uenew(i))**2
     &               + nuel(i)*(Cze*Umnew(i)-Uenew(i))**2
     &               + nuem(i)*(Cze*Umnew(i)-Uenew(i))**2
     &               + nuen(i)*(Cne*U3new(i)-Uenew(i))**2
     &                        )
     &               - 2.*(Cei*qenew(i)*alt_geo_1(i)*xne_1(i))

      D7e(i)=D7e(i)
     &          - 1.09e-3*(nueN2(i)/28.+ nueO2(i)/32. + nueO(i)/16.
     &          +nueOHot(i)/16.+ nueH(i) + nuei(i) + nuej(i)/16.        !MZ
     &          + nuek(i)/28. + nuel(i)/32. + nuem(i)/30.+nuen(i)/14.)
!     &          - 2.*(Cei*Uenew(i)*alt_geo_1(i))

      enddo

!    premiere demi-boucle

          D2qal=(D2qa(1)+D2qa(2))/2.
          D2qbl=(D2qb(1)+D2qb(2))/2.

          D2qar=.5*(Uenew(np)+Uenew(nx))
          D2qbr=.5*(Tenew(np)+Tenew(nx))
    
      call velocity(Veleq,Ipos1,Iposnp,deltat_4)
      call sources(Ipos1,Iposn,deltat_4,2,C2qa,D2qa,D2qal,D2qar)
      call sources(Ipos1,Iposn,deltat_4,2,C2qb,D2qb,D2qbl,D2qbr)
!      call sources(Ipos1,Iposn,deltat_4,3,zero,D3q,0.,0.)
      call sources(Ipos1,Iposn,deltat_4,7,zero,D7q,0.,0.)
          lbc=1.
          rbc=1.
      call lcpfct(qeold,qenew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          qenew(nx)=qetop

     
          do i=1,nx
            dexpnu=D7q(i)*deltat_4
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            qenew(i)=dexpnu*(qenew(i)-qeold(i)+deltat_4*D3q(i))
     &            +expnu*qeold(i)
          enddo
          qenew(nx)=qetop
          qenew(np)=qenew(nx)

      if (any(ieee_is_nan(qenew))) then
        call cpu_time(tic)
        write(stderr,*) tic,'problem when calculating qenew in loop 1'
          goto 246
      endif



          D2eal=(D2ea(1)+D2ea(2))/2.
          D2ebl=(D2eb(1)+D2eb(2))/2.

          D2ear=.5*(Uenew(np)+Uenew(nx))
          D2ebr=.5*(qenew(np)+qenew(nx))
    
      call sources(Ipos1,Iposn,deltat_4,2,C2ea,D2ea,D2eal,D2ear)
      call sources(Ipos1,Iposn,deltat_4,2,C2eb,D2eb,D2ebl,D2ebr)
!      call sources(Ipos1,Iposn,deltat_4,3,zero,D3e,0.,0.)
      call sources(Ipos1,Iposn,deltat_4,7,zero,D7e,0.,0.)
          lbc=1.
      rbc=1.
       call lcpfct(Teold,Tenew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)

          do i=1,nx
            dexpnu=D7e(i)*deltat_4
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Tenew(i)=dexpnu*(Tenew(i)-Teold(i)+deltat_4*D3e(i))
     &             +expnu*Teold(i)
          enddo
          tenew(np)=tenew(nx)
          Tenew(np)=2.*Tenew(nx)-Tenew(nx-1)

          if (any(ieee_is_nan(Tenew))) then
            call cpu_time(tic)
        write(stderr,*) tic,'problem when calculating Tenew in loop 1'
            goto 246
          endif


      call stabenerg(nx,xne,Teold,Tenew,qeold,qenew,
     &            D7e,D7q,Cei,deltat_4)


          if (any(ieee_is_nan(qenew))) then
            call cpu_time(tic)
        write(stderr,*) tic,'problem in stabenerg with qe in loop 1'
            goto 246
          endif
          if (any(ieee_is_nan(Tenew))) then
            call cpu_time(tic)
        write(stderr,*) tic,'problem in stabenerg with Te in loop 1'
            goto 246
          endif

! deuxieme demi-boucle

          do i=1,np
            Tenew(i)=max(Tenew(i),T_min)
        Tepnew(i)=Tenew(i)
        Tetnew(i)=Tenew(i)
            qeold(i)=qenew(i)
            Teold(i)=Tenew(i)
          enddo

          do i=1,nx
            C2ea(i)=Cei*Tenew(i)/3.
            D2eb(i)=qenew(i)

            C2qa(i)=-11./5.*Cei*qenew(i)

            C2qb(i)=-2.5*Cei*Nenew(i)*Tenew(i)
            D2qb(i)=Tenew(i)
      enddo

          D2qbl=(D2qb(1)+D2qb(2))/2.
          D2qbr=.5*(Tenew(np)+Tenew(nx))
    
      call sources(Ipos1,Iposn,deltat_4,2,C2qa,D2qa,D2qal,D2qar)
      call sources(Ipos1,Iposn,deltat_4,2,C2qb,D2qb,D2qbl,D2qbr)
!      call sources(Ipos1,Iposn,deltat_4,3,zero,D3q,0.,0.)
      call sources(Ipos1,Iposn,deltat_4,7,zero,D7q,0.,0.)
          lbc=1.
          rbc=1.
      call lcpfct(qeold,qenew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          qenew(nx)=qetop

     
          do i=1,nx
            dexpnu=D7q(i)*deltat_4
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            qenew(i)=dexpnu*(qenew(i)-qeold(i)+deltat_4*D3q(i))
     &            +expnu*qeold(i)
          enddo
          qenew(nx)=qetop
          qenew(np)=qenew(nx)

          if (any(ieee_is_nan(qenew))) then
            write(stderr,*) 'problem qenew dans la boucle 1'
            goto 246
          endif



          D2ebl=(D2eb(1)+D2eb(2))/2.
          D2ebr=.5*(qenew(np)+qenew(nx))
    
      call sources(Ipos1,Iposn,deltat_4,2,C2ea,D2ea,D2eal,D2ear)
      call sources(Ipos1,Iposn,deltat_4,2,C2eb,D2eb,D2ebl,D2ebr)
c      call sources(Ipos1,Iposn,deltat_4,3,zero,D3e,0.,0.)
      call sources(Ipos1,Iposn,deltat_4,7,zero,D7e,0.,0.)
          lbc=1.
      rbc=1.
       call lcpfct(Teold,Tenew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)

          do i=1,nx
            dexpnu=D7e(i)*deltat_4
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Tenew(i)=dexpnu*(Tenew(i)-Teold(i)+deltat_4*D3e(i))
     &             +expnu*Teold(i)
          enddo
          tenew(np)=tenew(nx)
          Tenew(np)=2.*Tenew(nx)-Tenew(nx-1)

          if (any(ieee_is_nan(Tenew))) then
            call cpu_time(tic)
            write(stderr,*) tic,'problem when calc Tenew in loop 1'
            goto 246
          endif


      call stabenerg(nx,xne,Teold,Tenew,qeold,qenew,
     &            D7e,D7q,Cei,deltat_4)


          if (any(ieee_is_nan(qenew))) then
            call cpu_time(tic)
         write(stderr,*) tic,'problem in stabenerg with qe in loop 1'
            goto 246
          endif
          if (any(ieee_is_nan(Tenew))) then
            write(stderr,*) 'problem stabenerg Te dans la boucle 1'
            goto 246
          endif
              


          do i=1,np
            Tenew(i)=max(Tenew(i),T_min)
        Tepnew(i)=Tenew(i)
        Tetnew(i)=Tenew(i)
            qeold(i)=qenew(i)
            Teold(i)=Tenew(i)
          enddo

C ]]]

C[[[    H+ energy equation resolution

      call velocity(Velie,Ipos1,Iposnp,deltat_2)

c    temperature parallele

      do i=1,nx

      C2a(i)=-T2pnew(i)
      D2a(i)=U2new(i)
      C2b(i)=-1.2*xn2_1(i)
      D2b(i)=q2new(i)

        Q_cross2=max(.3,min(.45,.45-.15*(T2new(i)-2.5)))
        Q_cross2=.3
          D7(i)=
     &      -(.9862*nuiN2(i)+.9818*nuiO2(i)
     &        +(.1176+.9412*Q_cross2)*nuiO(i)+(1.+.5*Q_cross2)*nuiH(i)
     &        +(.1176+.9412*Q_cross2)*nuiOHot(i)                    !MZ
     &       +nuij(i)*.8706+nuie(i)*1.9993
     &       +nuik(i)*.8414+nuil(i)*.8364+nuim(i)*.8387+nuin(i)*.88)
     &        -.8*nuii(i)

          D3(i)=T2tnew(i)*(.7529*nuij(i)+.7467*nuin(i)+.7724*nuik(i)
     &                    +.7758*nuil(i) +.7742*nuim(i)
     &              +4.37e-4*nuie(i)
     &                    +.9172*nuiN2(i)+.9212*nuiO2(i)
     &                    +.9412*Q_cross2*nuiO(i)+.5*Q_cross2*nuiH(i)
     &              +.9412*Q_cross2*nuiOHot(i))                !MZ
     &           +T1pnew(i)*.0706*nuij(i)
     &         +Tmpnew(i)*(.0414*nuik(i)+.0364*nuil(i) +.0387*nuim(i))
     &         +Tepnew(i)*  1.1993*nuie(i)
     &         +T3pnew(i)* .0800*nuin(i)
     &         +T1tnew(i)*.0471*nuij(i)
     &         +Tmtnew(i)*(.0276*nuik(i)+.0242*nuil(i) +.0258*nuim(i))
     &         +Tetnew(i)*  .7996*nuie(i)
     &         +T3tnew(i)* .0533*nuin(i)
     &           +.8*nuii(i)*T2tnew(i)

        nu_omega=(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+
     &            nuiH(i))/omega(i)*omi    !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Ci0**2
        dUn2=(U2new(i)-Un(i)/Ci0)**2

          D3(i)=D3(i)
     &           +Tn(i)/T_0*(.069*nuiN2(i)+.0606*nuiO2(i)
     &                    +nuiH(i)+.1176*nuiO(i))
     &            +TnOHot(i)/T_0*.1176*nuiOhot(i)                !MZ
     &         +nuiN2(i)*(1.0138*dUn2+.9172*terme_joule)
     &         +nuiO2(i)*(1.0182*dUn2+.9212*terme_joule)
     &         +nuiO(i)*.9412*(2.*dUn2+Q_cross2*(terme_joule-dUn2))
     &         +nuiOHot(i)*.9412*(2.*dUn2+Q_cross2*(terme_joule-dUn2))        !MZ
     &         +nuiH(i)*.5*(2.*dUn2+Q_cross2*(terme_joule-dUn2))
     &         +.667*(
     &         + nuij(i)*.94*(Cji*U1new(i)-U2new(i))**2
     &         + (nuik(i)+ nuil(i)+ nuim(i))
     &           *.97*(Czi*Umnew(i)-U2new(i))*(Czi*Umnew(i)-U2new(i))
     &         + nuin(i)*.93*(Cni*U3new(i)-U2new(i))**2
     &         + nuie(i)*5.46e-4*(Cei*Uenew(i)-U2new(i))**2)
     &         - 1.2*(q2new(i)*alt_geo_1(i)*xn2_1(i))

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
          D2ar=.5*(U2new(np)+U2new(nx))
          D2br=.5*(q2new(np)+q2new(nx))
        
          call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
          call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
          call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T2pold,extra,nx)
      rbc=1.
       call lcpfct(T2pold,T2pnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T2pnew(i)=dexpnu*(T2pnew(i)-T2pold(i)+D3(i)*deltat_2)
     &             +expnu*T2pold(i)
          enddo
          t2pnew(np)=t2pnew(nx)

          do  i=1,np
            T2pnew(i)=max(T2pnew(i),T_min)
      enddo

      if (any(ieee_is_nan(T2pnew))) then
        write(stderr,*) 'problem T2pnew dans la boucle 1'
        goto 246
      endif



C]]]

c    temperature perpendiculaire

      do i=1,nx

      C2a(i)=T2tnew(i)
      D2a(i)=U2new(i)
      C2b(i)=-.4*xn2_1(i)
      D2b(i)=q2new(i)

        Q_cross2=max(.3,min(.45,.45-.15*(T2new(i)-2.5)))/2.
        Q_cross2=.3
          D7(i)=
     &      -(.5276*nuiN2(i)+.5212*nuiO2(i)
     &        +(.1176+.9412*Q_cross2)*nuiO(i)+(1.+.5*Q_cross2)*nuiH(i)
     &        +(.1176+.9412*Q_cross2)*nuiOHot(i)                    !MZ
     &       +nuij(i)*.4941+nuie(i)*1.9991
     &       +nuik(i)*.4552+nuil(i)*.4485+nuim(i)*.4516+nuin(i)*.5067)
     &        -.4*nuii(i)

          D3(i)=T2pnew(i)*(.3765*nuij(i)+.3733*nuin(i)+.3862*nuik(i)
     &                    +.3879*nuil(i) +.3871*nuim(i)
     &              +2.1828e-4*nuie(i)
     &                    +.4586*nuiN2(i)+.4606*nuiO2(i)
     &                    +.9412*Q_cross2*nuiO(i)+.5*Q_cross2*nuiH(i)
     &              +.9412*Q_cross2*nuiOHot(i))                !MZ
     &           +T1tnew(i)*.0941*nuij(i)
     &         +Tmtnew(i)*(.0552*nuik(i)+.0485*nuil(i) +.0516*nuim(i))
     &         +Tetnew(i)*  1.5991*nuie(i)
     &         +T3tnew(i)* .1067*nuin(i)
     &         +T1pnew(i)*.0235*nuij(i)
     &         +Tmpnew(i)*(.0138*nuik(i)+.0121*nuil(i) +.0129*nuim(i))
     &         +Tepnew(i)* .3998*nuie(i)
     &         +T3pnew(i)* .0267*nuin(i)
     &           +.4*nuii(i)*T2pnew(i)

        nu_omega=(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+ 
     &            nuiH(i))/omega(i)*omi    !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Ci0**2
        dUn2=(U2new(i)-Un(i)/Ci0)**2

          D3(i)=D3(i)
     &           +Tn(i)/T_0*(.069*nuiN2(i)+.0606*nuiO2(i)
     &        +nuiH(i)+.1176*nuiO(i))
     &        +TnOHot(i)/T_0*.1176*nuiOHot(i)                    !MZ
     &         +nuiN2(i)*(.4586*dUn2+1.4724*terme_joule)
     &         +nuiO2(i)*(.4606*dUn2+1.4788*terme_joule)
     &         +nuiO(i)*.9412*((2.-Q_cross2)*terme_joule+Q_cross2*dUn2)
     &      +nuiOHot(i)*.9412*((2.-Q_cross2)*terme_joule+Q_cross2*dUn2)    !MZ
     &         +nuiH(i)*.5*((2.-Q_cross2)*terme_joule+Q_cross2*dUn2)
     &         +.667*(
     &         + nuij(i)*.94*(Cji*U1new(i)-U2new(i))**2
     &         + (nuik(i)+ nuil(i)+ nuim(i))
     &           *.97*(Czi*Umnew(i)-U2new(i))*(Czi*Umnew(i)-U2new(i))
     &         + nuin(i)*.93*(Cni*U3new(i)-U2new(i))**2
     &         + nuie(i)*5.46e-4*(Cei*Uenew(i)-U2new(i))**2)
     &         - 2.4*(q2new(i)*alt_geo_1(i)*xn2_1(i))

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
          D2ar=.5*(U2new(np)+U2new(nx))
          D2br=.5*(q2new(np)+q2new(nx))
        
           call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
          call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
          call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T2told,extra,nx)
          rbc=1.
          call lcpfct(T2told,T2tnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T2tnew(i)=dexpnu*(T2tnew(i)-T2told(i)+D3(i)*deltat_2)
     &             +expnu*T2told(i)
          enddo
          T2tnew(np)=T2tnew(nx)

          do i=1,np
            T2tnew(i)=max(T2tnew(i),T_min)
            T2new(i)=(T2pnew(i)+2.*T2tnew(i))/3.
          enddo

      if (any(ieee_is_nan(T2tnew))) then
        write(stderr,*) 'problem T2tnew dans la boucle 1'
        goto 246
      endif


C]]]

       call velocity(Velje,Ipos1,Iposnp,deltat_2)

C[[[    O+ energy equation resolution


c    temperature parallele

      do i=1,nx

      C2a(i)=-Cji*T1pnew(i)
      D2a(i)=U1new(i)
          C2b(i)=-1.2*Cji*xn1_1(i)
      D2b(i)=q1new(i)

        Q_cross2=max(.3,min(.45,.45-.15*(T1new(i)-2.5)))/2.
        Q_cross2=.3
          D7(i)=
     &      -(1.3318*nujN2(i)+1.3*nujO2(i)
     &        +(1.+Q_cross2)*nujO(i)
     &        +(1.+Q_cross2)*nujOHot(i)                        !MZ
     &       +nuji(i)*1.9294+nuje(i)*2.
     &       +nujk(i)*1.2364+nujl(i)*1.2+nujm(i)*1.2174+nujn(i)*1.44)
     &        -.8*nujj(i)

          D3(i)=T1tnew(i)*(.0471*nuji(i)+.3733*nujn(i)+.5091*nujk(i)
     &                    +.5333*nujl(i) +.5217*nujm(i)
     &              +2.7299e-5*nuje(i)
     &                    +.6045*nujN2(i)+.6333*nujO2(i)
     &              +Q_cross2*nujO(i)
     &              +Q_cross2*nujOHot(i))                    !MZ
     &           +T2pnew(i)*1.1294*nuji(i)
     &         +Tmpnew(i)*(.4364*nujk(i)+.4*nujl(i) +.4174*nujm(i))
     &         +Tepnew(i)*  1.2*nuje(i)
     &         +T3pnew(i)* .64*nujn(i)
     &         +T2tnew(i)*.7529*nuji(i)
     &         +Tmtnew(i)*(.2909*nujk(i)+.2667*nujl(i) +.2783*nujm(i))
     &         +Tetnew(i)*  .8*nuje(i)
     &         +T3tnew(i)* .4267*nujn(i)
     &           +.8*nujj(i)*T1tnew(i)

        nu_omega=(nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i))/omega(i)*omj        !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cj0**2
        dUn2=(U1new(i)-Un(i)/Cj0)**2

         D3(i)=D3(i)
     &           +Tn(i)/T_0*(.7273*nujN2(i)+.6667*nujO2(i)
     &                    +nujO(i))
     &            +TnOHot(i)/T_0*nujOHot(i)                !MZ
     &         +nujN2(i)*(.6682*dUn2+.6045*terme_joule)
     &         +nujO2(i)*(.7*dUn2+.6333*terme_joule)
     &         +nujO(i)*(dUn2+Q_cross2*(terme_joule-dUn2))
     &         +nujOHot(i)*(dUn2+Q_cross2*(terme_joule-dUn2))            !MZ
     &      + .667*(nuji(i)/17.*(Cij*U2new(i)-U1new(i))**2
     &      + (nujk(i)*.64+ nujl(i)*.667+ nujm(i)*.65)
     &                *(Czj*Umnew(i)-U1new(i))*(Czj*Umnew(i)-U1new(i))
     &      + nujn(i)*.467*(Cnj*U3new(i)-U1new(i))**2
     &      + nuje(i)*5.46e-4/16.*(Cej*Uenew(i)-U1new(i))**2)
     &      - 1.2*(Cji*q1new(i)*alt_geo_1(i)*xn1_1(i))


        enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*(U1new(np)+U1new(nx))
      D2br=.5*(q1new(np)+q1new(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T1pold,extra,nx)
      rbc=1.
       call lcpfct(T1pold,T1pnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T1pnew(i)=dexpnu*(T1pnew(i)-T1pold(i)+D3(i)*deltat_2)
     &             +expnu*T1pold(i)
          enddo
          T1pnew(np)=T1pnew(nx)

      do i=1,np
        T1pnew(i)=max(T1pnew(i),T_min)
      enddo

      if (any(ieee_is_nan(T1pnew))) then
        write(stderr,*) 'problem T1pnew dans la boucle 1'
        goto 246
      endif



c    temperature perpendiculaire

      do i=1,nx

      C2a(i)=Cji*T1tnew(i)
      D2a(i)=U1new(i)
          C2b(i)=-.4*Cji*xn1_1(i)
      D2b(i)=q1new(i)

        Q_cross2=max(.3,min(.45,.45-.15*(T1new(i)-2.5)))/4.
        Q_cross2=.3
          D7(i)=
     &      -(1.0295*nujN2(i)+.9833*nujO2(i)
     &        +(1.+Q_cross2)*nujO(i)
     &        +(1.+Q_cross2)*nujOHot(i)                        !MZ
     &       +nuji(i)*1.9059+nuje(i)*1.9999
     &       +nujk(i)*.9818+nujl(i)*.9333+nujm(i)*.9565+nujn(i)*1.2533)
     &        -.4*nujj(i)

         D3(i)=T1pnew(i)*(.0235*nuji(i)+.1867*nujn(i)+.2545*nujk(i)
     &                    +.2667*nujl(i) +.2609*nujm(i)
     &              +1.365e-5*nuje(i)
     &                    +.3023*nujN2(i)+.3167*nujO2(i)
     &                    +Q_cross2*nujO(i)
     &                    +Q_cross2*nujOHot(i))                    !MZ
     &           +T2tnew(i)*1.5059*nuji(i)
     &         +Tmtnew(i)*(.5818*nujk(i)+.5333*nujl(i) +.5565*nujm(i))
     &         +Tetnew(i)*  1.5999*nuje(i)
     &         +T3tnew(i)* .8533*nujn(i)
     &         +T2pnew(i)*.3765*nuji(i)
     &         +Tmpnew(i)*(.1455*nujk(i)+.1333*nujl(i) +.1391*nujm(i))
     &         +Tepnew(i)*  .4*nuje(i)
     &         +T3pnew(i)* .2133*nujn(i)
     &           +.4*nujj(i)*T1pnew(i)


        nu_omega=(nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i))/omega(i)*omj        !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cj0**2
        dUn2=(U1new(i)-Un(i)/Cj0)**2

        D3(i)=D3(i)
     &           +Tn(i)/T_0*(.7273*nujN2(i)+.6667*nujO2(i)
     &                    +nujO(i))
     &            +TnOHot(i)/T_0*nujOHot(i)                !MZ
     &         +nujN2(i)*(.3023*dUn2+.9705*terme_joule)
     &         +nujO2(i)*(.3167*dUn2+1.0167*terme_joule)
     &         +nujO(i)*(Q_cross2*dUn2+(1.-Q_cross2)*terme_joule)
     &         +nujOHot(i)*(Q_cross2*dUn2+(1.-Q_cross2)*terme_joule)        !MZ
     &      + .667*(nuji(i)/17.*(Cij*U2new(i)-U1new(i))**2
     &      + (nujk(i)*.64+ nujl(i)*.667+ nujm(i)*.65)
     &                *(Czj*Umnew(i)-U1new(i))**2
     &      + nujn(i)*.467*(Cnj*U3new(i)-U1new(i))**2
     &      + nuje(i)*5.46e-4/16.*(Cej*Uenew(i)-U1new(i))**2)
     &      - 2.4*(Cji*q1new(i)*alt_geo_1(i)*xn1_1(i))

       enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
          D2ar=.5*(U1new(np)+U1new(nx))
          D2br=.5*(q1new(np)+q1new(nx))
        
          call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
          call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
          call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T1told,extra,nx)
        rbc=1.
        call lcpfct(T1told,T1tnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T1tnew(i)=dexpnu*(T1tnew(i)-T1told(i)+D3(i)*deltat_2)
     &             +expnu*T1told(i)
          enddo
          T1tnew(np)=T1tnew(nx)

          do  i=1,np
            T1tnew(i)=max(T1tnew(i),T_min)
        T1new(i)=(T1pnew(i)+2.*T1tnew(i))/3.
          enddo


      if (any(ieee_is_nan(T1tnew))) then
        write(stderr,*) 'problem T1tnew dans la boucle 1'
        goto 246
      endif


C]]]



       call velocity(Velme,Ipos1,Iposnp,deltat_2)

C[[[    heavy ions energy equation resolution

c    temperature parallele

      do i=1,nx

      C2a(i)=-Czi*Tmpnew(i)
      D2a(i)=Umnew(i)

      Q_cross2=max(.3,min(.45,.45-.15*(Tmnew(i)-2.5)))/2.
        Q_cross2=.3

          D7(i)=
     &      -((1.+Q_cross2)*nukN2(i)+1.44*nukO2(i)+1.6182*nukO(i)
     &        +1.6182*nukOHot(i)                            !MZ
     &       +nuki(i)*1.9586+nukj(i)*1.5636+nuke(i)*2.
     &       +nukn(i)*1.6)
     &        -(1.51*nulN2(i)+(1.+Q_cross2)*nulO2(i)+1.65*nulO(i)
     &        +1.65*nulOHot(i)                        !MZ
     &       +nuli(i)*1.9636+nulj(i)*1.6+nule(i)*2.
     &       +nuln(i)*1.6348)
     &      -(1.4931*numN2(i)+1.4581*numO2(i)+1.6348*numO(i)
     &        +1.6348*numOHot(i)                        !MZ
     &       +numi(i)*1.9613+numj(i)*1.5826+nume(i)*2.
     &       +numn(i)*1.44)
     &         -.8*(nukk(i)+nukl(i)+nukm(i)
     &           +nulk(i)+null(i)+nulm(i)
     &           +numk(i)+numl(i)+numm(i))

        D3(i)=Tmtnew(i)*(.0276*nuki(i)+.2909*nukj(i)+.2667*nukn(i)
     &                 +1.56e-5*nuke(i)
     &                 +.0242*nuli(i)+.2667*nulj(i)+.2435*nuln(i)
     &                 + 1.365e-5*nule(i)
     &                 +.0258*numi(i)+.2783*numj(i)+.2545*numn(i)
     &                 +1.456e-5*nume(i)
     &                 +Q_cross2*nukN2(i)+.5067*nukO2(i)+.3455*nukO(i)
     &               +.3455*nukOHot(i)                    !MZ
     &                 +.4433*nulN2(i)+Q_cross2*nulO2(i)+.3167*nulO(i)
     &               +.3167*nulOHot(i)                    !MZ
     &                 +.4586*numN2(i)+.4903*numO2(i)+.3304*numO(i)
     &               +.3304*numOHot(i))                    !MZ
     &           +T2pnew(i)*(1.1586*nuki(i)
     &                    +1.1636*nuli(i)
     &                    +1.1613*numi(i))
     &         +T1pnew(i)*(.7636*nukj(i)
     &                    +.8*nulj(i)
     &                    +.7826*numj(i))
     &         +Tepnew(i)*( 1.2*nuke(i)
     &                    +1.2*nule(i)
     &                    +1.2*nume(i))
     &         +T3pnew(i)*(.8*nukn(i)
     &                    +.8348*nuln(i)
     &                    +.8182*numn(i))
     &         +T2tnew(i)*(.7724*nuki(i)
     &                    +.7758*nuli(i)
     &                    +.7742*numi(i))
     &         +T1tnew(i)*(.5091*nukj(i)
     &                    +.5333*nulj(i)
     &                    +.5217*numj(i))
     &         +Tetnew(i)*( .8*nuke(i)
     &                    +.8*nule(i)
     &                    +.8*nume(i))
     &         +T3tnew(i)*(.5333*nukn(i)
     &                    +.5565*nuln(i)
     &                    +.5455*numn(i))
     &           +.8*(nukk(i)+nukl(i)+nukm(i)
     &               +nulk(i)+null(i)+nulm(i)
     &               +numk(i)+numl(i)+numm(i))*Tmtnew(i)


        nu_omega=(nukN2(i)+nukO2(i)+nukO(i)+nukOHot(i)                !MZ
     &               +nulN2(i)+nulO2(i)+nulO(i)+nulOHot(i)            !MZ
     &               +numN2(i)+numO2(i)+numO(i)+numOHot(i))/omega(i)*omz    !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cz0**2
        dUn2=(Umnew(i)-Un(i)/Cz0)**2

        D3(i)=D3(i)
     &           +Tn(i)/T_0
     &        *(nukN2(i)+.9333*nukO2(i)+1.2727*nukO(i)
     &         +1.0667*nulN2(i)+nulO2(i)+1.3333*nulO(i)
     &         +1.0345*numN2(i)+0.9677*numO2(i)+1.3043*numO(i))
     &        +TnOHot(i)/T_0*(1.2727*nukOHot(i)+1.3333*nulOHot(i)        !MZ
     &         +1.3043*numOHot(i))                        !MZ
     &         +nukN2(i)*(dUn2+Q_cross2*(terme_joule-dUn2))
     &         +nukO2(i)*(.5242*dUn2+.4743*terme_joule)
     &         +nukO(i)*(.3574*dUn2+.3234*terme_joule)
     &         +nukOHot(i)*(.3574*dUn2+.3234*terme_joule)            !MZ
     &         +nulN2(i)*(.5242*dUn2+.4743*terme_joule)
     &         +nulO2(i)*(dUn2+Q_cross2*(terme_joule-dUn2))
     &         +nulO(i)*(.3744*dUn2+.3388*terme_joule)
     &         +nulOHot(i)*(.3744*dUn2+.3388*terme_joule)            !MZ
     &         +numN2(i)*(.5084*dUn2+.46*terme_joule)
     &         +numO2(i)*(.5435*dUn2+.4918*terme_joule)
     &         +numO(i)*(.3663*dUn2+.3314*terme_joule)
     &         +numOHot(i)*(.3663*dUn2+.3314*terme_joule)            !MZ
     &      + .22*(        (nuki(i)/29.+nuli(i)/33.+numi(i)/31.)
     &                       *(Ciz*U2new(i)-Umnew(i))**2
     &            +    16.*(nukj(i)/44.+nulj(i)/48.+numj(i)/46.)
     &                   *(Cjz*U1new(i)-Umnew(i))**2
     &            +    14.*(nukn(i)/42.+nuln(i)/46.+numn(i)/44.)
     &                   *(Cnz*U3new(i)-Umnew(i))**2
     &            +5.46e-4*(nuke(i)/28.+nule(i)/32.+nume(i)/30.)
     &               *(Cez*Uenew(i)-Umnew(i))**2
     &          )

 
        enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2ar=(D2a(nx)+D2a(nx-1))/2.
    
          D2ar=.5*(Umnew(np)+Umnew(nx))
        
          call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
c      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
        call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
      rbc=bclimd(Radn(nx+1),Tmpold,extra,nx)
      rbc=1.
       call lcpfct(Tmpold,Tmpnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif

            Tmpnew(i)=dexpnu*(Tmpnew(i)-Tmpold(i)+D3(i)*deltat_2)
     &             +expnu*Tmpold(i)
          enddo
          Tmpnew(np)=Tmpnew(nx)

      do i=1,np
         Tmpnew(i)=max(Tmpnew(i),T_min)
          enddo
c    flag=.false.

      if (any(ieee_is_nan(Tmpnew))) then
        write(stderr,*) 'problem Tmpnew dans la boucle 1'
        goto 246
      endif


c    temperature perpendiculaire

        do i=1,nx

      C2a(i)=Czi*Tmtnew(i)
      D2a(i)=Umnew(i)

        Q_cross2=max(.3,min(.45,.45-.15*(Tmnew(i)-2.5)))/4.
        Q_cross2=.3

          D7(i)=
     &      -((1.+Q_cross2)*nukN2(i)+1.1867*nukO2(i)+1.4455*nukO(i)
     &        +1.4455*nukOHot(i)                            !MZ
     &       +nuki(i)*1.9448+nukj(i)*1.4182+nuke(i)*2.
     &       +nukn(i)*1.4667)
     &      -(1.2883*nulN2(i)+(1.+Q_cross2)*nulO2(i)+1.4917*nulO(i)
     &        +1.4917*nulOHot(i)                            !MZ
     &       +nuli(i)*1.9515+nulj(i)*1.4667+nule(i)*2.
     &       +nuln(i)*1.513)
     &      -(1.2638*numN2(i)+1.2129*numO2(i)+1.4696*numO(i)
     &        +1.4696*numOHot(i)                            !MZ
     &       +numi(i)*1.9484+numj(i)*1.4435+nume(i)*2.
     &       +numn(i)*1.4909)
     &      -.4*(nukk(i)+nukl(i)+nukm(i)
     &          +nulk(i)+null(i)+nulm(i)
     &          +numk(i)+numl(i)+numm(i))


        D3(i)=Tmpnew(i)*(.0138*nuki(i)+.1455*nukj(i)+.1333*nukn(i)
     &                 +7.7998e-6*nuke(i)
     &                 +.0121*nuli(i)+.1333*nulj(i)+.1217*nuln(i)
     &                 +6.8249e-6*nule(i)
     &                 +.0129*numi(i)+.1391*numj(i)+.1273*numn(i)
     &                 +7.2799e-6*nume(i)
     &                 +Q_cross2*nukN2(i)+.2533*nukO2(i)+.1727*nukO(i)
     &            +.1727*nukOHot(i)                    !MZ
     &                 +.2217*nulN2(i)+Q_cross2*nulO2(i)+.1583*nulO(i)
     &            +.1583*nulOHot(i)                    !MZ
     &                 +.2293*numN2(i)+.2452*numO2(i)+.1652*numO(i)
     &            +.1652*numOHot(i))                    !MZ
     &        +T2tnew(i)*(1.5448*nuki(i)+1.5515*nuli(i)+1.5484*numi(i))
     &        +T1tnew(i)*(1.0182*nukj(i)+1.0667*nulj(i)+1.0435*numj(i))
     &         +Tetnew(i)*( 1.6*nuke(i)+1.6*nule(i)+1.6*nume(i))
     &        +T3tnew(i)*(1.0667*nukn(i)+1.1130*nuln(i)+1.0909*numn(i))
     &         +T2pnew(i)*(.3862*nuki(i)+.3879*nuli(i)+.3871*numi(i))
     &         +T1pnew(i)*(.2545*nukj(i)+.2667*nulj(i)+.2609*numj(i))
     &         +Tepnew(i)*(.4*nuke(i)+.4*nule(i)+.4*nume(i))
     &         +T3pnew(i)*(.2667*nukn(i)+.2783*nuln(i)+.2727*numn(i))
     &         +.4*(nukk(i)+nukl(i)+nukm(i)
     &             +nulk(i)+null(i)+nulm(i)
     &             +numk(i)+numl(i)+numm(i))*Tmpnew(i)

        nu_omega=(nukN2(i)+nukO2(i)+nukO(i)+nukOHot(i)                !MZ
     &               +nulN2(i)+nulO2(i)+nulO(i)+nulOHot(i)            !MZ
     &               +numN2(i)+numO2(i)+numO(i)+numOHot(i))/omega(i)*omz    !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cz0**2
        dUn2=(Umnew(i)-Un(i)/Cz0)**2

        D3(i)=D3(i)
     &           +Tn(i)/T_0
     &        *(nukN2(i)+.9333*nukO2(i)+1.2727*nukO(i)
     &         +1.0667*nulN2(i)+nulO2(i)+1.3333*nulO(i)
     &         +1.0345*numN2(i)+0.9677*numO2(i)+1.3043*numO(i))
     &        +TnOHot(i)/T_0*(1.2727*nukOHot(i)+1.3333*nulOHot(i)        !MZ
     &         +1.3043*numOHot(i))                        !MZ
     &         +nukN2(i)*(Q_cross2*dUn2+(1.-Q_cross2)*terme_joule)
     &         +nukO2(i)*(.2371*dUn2+.7614*terme_joule)
     &         +nukO(i)*(.1617*dUn2+.5191*terme_joule)
     &         +nukOHot(i)*(.1617*dUn2+.5191*terme_joule)            !MZ
     &         +nulN2(i)*(.2371*dUn2+.7614*terme_joule)
     &         +nulO2(i)*(Q_cross2*dUn2+(1.-Q_cross2)*terme_joule)
     &         +nulO(i)*(.1694*dUn2+.5438*terme_joule)
     &         +nulOHot(i)*(.1694*dUn2+.5438*terme_joule)            !MZ
     &         +numN2(i)*(.23*dUn2+.7384*terme_joule)
     &         +numO2(i)*(.1657*dUn2+.7894*terme_joule)
     &         +numO(i)*(.3663*dUn2+.532*terme_joule)
     &         +numOHot(i)*(.3663*dUn2+.532*terme_joule)            !MZ
     &      + .22*(        (nuki(i)/29.+nuli(i)/33.+numi(i)/31.)
     &                       *(Ciz*U2new(i)-Umnew(i))**2
     &            +    16.*(nukj(i)/44.+nulj(i)/48.+numj(i)/46.)
     &                   *(Cjz*U1new(i)-Umnew(i))**2
     &            +    14.*(nukn(i)/42.+nuln(i)/46.+numn(i)/44.)
     &                   *(Cnz*U3new(i)-Umnew(i))**2
     &            +5.46e-4*(nuke(i)/28.+nule(i)/32.+nume(i)/30.)
     &               *(Cez*Uenew(i)-Umnew(i))**2
     &          )

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
    
      D2ar=.5*(Umnew(np)+Umnew(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
c      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
      rbc=bclimd(Radn(nx+1),Tmtold,extra,nx)
      rbc=1.
       call lcpfct(Tmtold,Tmtnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Tmtnew(i)=dexpnu*(Tmtnew(i)-Tmtold(i)+D3(i)*deltat_2)
     &             +expnu*Tmtold(i)
          enddo
          Tmtnew(np)=Tmtnew(nx)


      do i=1,np
         Tmtnew(i)=max(Tmtnew(i),T_min)
        Tmnew(i)=(Tmpnew(i)+2.*Tmtnew(i))/3.
          enddo

      if (any(ieee_is_nan(Tmtnew))) then
        write(stderr,*) 'problem Tmtnew dans la boucle 1'
        goto 246
      endif


C]]]

        call velocity(Velne,Ipos1,Iposnp,deltat_2)

C[[[    N+ energy equation resolution

c    temperature parallele

      do i=1,nx

      C2a(i)=-Cni*T3pnew(i)
      D2a(i)=U3new(i)
      C2b(i)=-1.2*Cni*xn3_1(i)
      D2b(i)=q3new(i)

          D7(i)=
     &      -(1.3*nunN2(i)+1.2696*nunO2(i)+1.44*nunO(i)+1.44*nunOHot(i)        !MZ
     &       +nuni(i)*1.92+nunj(i)*1.36+nune(i)*2.
     &       +nunk(i)*1.2+nunl(i)*1.1652+nunm(i)*1.1818)
     &        -.8*nunn(i)

          D3(i)=T3tnew(i)*(.0533*nuni(i)+.4267*nunj(i)+.5333*nunk(i)
     &                    +.5565*nunl(i) +.5455*nunm(i)
     &              +3.1199e-5*nune(i)
     &                    +.6333*nunN2(i)+.6609*nunO2(i)+.5067*nunO(i)
     &              +.5067*nunOHot(i))                    !MZ
     &           +T1pnew(i)*.56*nunj(i)
     &         +T2pnew(i)* 1.12*nuni(i)
     &         +Tmpnew(i)*(.4*nunk(i)+.3652*nunl(i) +.3818*nunm(i))
     &         +Tepnew(i)*  1.2*nune(i)
     &         +T1tnew(i)* .3733*nunj(i)
     &         +T2tnew(i)*.7467*nuni(i)
     &         +Tmtnew(i)*(.2667*nunk(i)+.2435*nunl(i) +.2545*nunm(i))
     &         +Tetnew(i)*  .8*nune(i)
     &           +.8*nunn(i)*T3tnew(i)

        nu_omega=(nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i))/omega(i)*omn        !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cn0**2
        dUn2=(U3new(i)-Un(i)/Cn0)**2

         D3(i)=D3(i)
     &           +Tn(i)/T_0*(.6667*nunN2(i)+.6087*nunO2(i)
     &                    +.9333*nunO(i))
     &        +TnOHot(i)/T_0*.9333*nunOHot(i)                    !MZ
     &         +nunN2(i)*(.7*dUn2+.6333*terme_joule)
     &         +nunO2(i)*(.7304*dUn2+.6609*terme_joule)
     &         +nunO(i)*(.56*dUn2+.5067*terme_joule)
     &         +nunOHot(i)*(.56*dUn2+.5067*terme_joule)                !MZ
     &      + .667*(nuni(i)*.067*(Cin*U2new(i)-U3new(i))**2
     &               + nunk(i)*.667*(Czn*Umnew(i)-U3new(i))**2
     &               + nunl(i)*.7*(Czn*Umnew(i)-U3new(i))**2
     &               + nunm(i)*.68*(Czn*Umnew(i)-U3new(i))**2
     &               + nunj(i)*.53*(Cjn*U1new(i)-U3new(i))**2
     &               + nune(i)*3.9e-5*(Cen*Uenew(i)-U3new(i))**2)
     &               - 1.2*(Cni*q3new(i)*alt_geo_1(i)*xn3_1(i))

       enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*(U3new(np)+U3new(nx))
      D2br=.5*(q3new(np)+q3new(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T3pold,extra,nx)
      rbc=1.
       call lcpfct(T3pold,T3pnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T3pnew(i)=dexpnu*(T3pnew(i)-T3pold(i)+D3(i)*deltat_2)
     &             +expnu*T3pold(i)
          enddo
          T3pnew(np)=T3pnew(nx)

      do i=1,np
            T3pnew(i)=max(T3pnew(i),T_min)
c            T3pnew(i)=T1pnew(i)
      enddo

      if (any(ieee_is_nan(T3pnew))) then
        write(stderr,*) 'problem T3pnew dans la boucle 1'
        goto 246
      endif


c    temperature perpendiculaire

      do i=1,nx

      C2a(i)=Cni*T3tnew(i)
      D2a(i)=U3new(i)
      C2b(i)=-.4*Cni*xn3_1(i)
      D2b(i)=q3new(i)

          D7(i)=
     &      -(.9833*nunN2(i)+.9391*nunO2(i)+1.1867*nunO(i)
     &        +1.1867*nunOHot(i)                        !MZ
     &       +nuni(i)*1.8933+nunj(i)*1.1467+nune(i)*1.9999
     &       +nunk(i)*.9333+nunl(i)*.887+nunm(i)*.9091)
     &        -.4*nunn(i)

         D3(i)=T3pnew(i)*(.0267*nuni(i)+.2133*nunj(i)+.2667*nunk(i)
     &                    +.2783*nunl(i) +.2727*nunm(i)
     &              +1.5599e-5*nune(i)
     &                    +.3167*nunN2(i)+.3304*nunO2(i)+.2533*nunO(i)
     &              +.2533*nunOHot(i))                    !MZ
     &           +T1tnew(i)*.7467*nunj(i)
     &         +T2tnew(i)* 1.4933*nuni(i)
     &         +Tmtnew(i)*(.5333*nunk(i)+.487*nunl(i) +.5091*nunm(i))
     &         +Tetnew(i)*  1.5999*nune(i)
     &         +T1tnew(i)*.1867*nunj(i)
     &         +T2pnew(i)*.3733*nuni(i)
     &         +Tmpnew(i)*(.1333*nunk(i)+.1217*nunl(i) +.1273*nunm(i))
     &         +Tepnew(i)*  .4*nune(i)
     &           +.4*nunn(i)*T3pnew(i)

        nu_omega=(nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i))/omega(i)*omn        !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cn0**2
        dUn2=(U3new(i)-Un(i)/Cn0)**2

         D3(i)=D3(i)
     &           +Tn(i)/T_0*(.6667*nunN2(i)+.6087*nunO2(i)
     &                    +.9333*nunO(i))
     &        +TnOHot(i)/T_0*.9333*nunOHot(i)                    !MZ
     &         +nunN2(i)*(.3167*dUn2+1.0167*terme_joule)
     &         +nunO2(i)*(.3304*dUn2+1.0609*terme_joule)
     &         +nunO(i)*(.2533*dUn2+.8133*terme_joule)
     &         +nunOHot(i)*(.2533*dUn2+.8133*terme_joule)            !MZ
     &      + .667*(nuni(i)*.067*(Cin*U2new(i)-U3new(i))**2
     &               + nunk(i)*.667*(Czn*Umnew(i)-U3new(i))**2
     &               + nunl(i)*.7*(Czn*Umnew(i)-U3new(i))**2
     &               + nunm(i)*.68*(Czn*Umnew(i)-U3new(i))**2
     &               + nunj(i)*.53*(Cjn*U1new(i)-U3new(i))**2
     &               + nune(i)*3.9e-5*(Cen*Uenew(i)-U3new(i))**2)
     &               - 2.4*(Cni*q3new(i)*alt_geo_1(i)*xn3_1(i))

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=U3new(np)
      D2br=q3new(np)
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T3told,extra,nx)
      rbc=1.
       call lcpfct(T3told,T3tnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat_2
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T3tnew(i)=dexpnu*(T3tnew(i)-T3told(i)+D3(i)*deltat_2)
     &             +expnu*T3told(i)
          enddo
          T3tnew(np)=T3tnew(nx)

      do i=1,np
            T3tnew(i)=max(T3tnew(i),T_min)
c            T3tnew(i)=T1tnew(i)
        T3new(i)=(T3pnew(i)+2.*T3tnew(i))/3.
      enddo

      if (any(ieee_is_nan(T3tnew))) then
        write(stderr,*) 'problem T3tnew dans la boucle 1'
        goto 246
      endif

123    continue

C]]]

C[[     second half-step

      do i=2,nx
        Velic(i)=.5*(U2new(i)+U2new(i-1))+vtrans
        Veljc(i)=Cji*.5*(U1new(i)+U1new(i-1))+vtrans
        Velec(i)=Cei*.5*(Uenew(i)+Uenew(i-1))+vtrans
        Velmc(i)=Czi*.5*(Umnew(i)+Umnew(i-1))+vtrans
        Velnc(i)=Cni*.5*(U3new(i)+U3new(i-1))+vtrans
          enddo

          Velnoc(1)=Velnoc(2)
          rbc=min(1.,Nnonew(nx-1)/Nnonew(nx-2))
          rbc_2=sqrt(rbc)
          Nnonew(np)=rbc*Nnonew(nx)
          Unonew(np)=.5*(Unonew(nx)+Unonew(nx-1))
          if (Unonew(np).ge.0.) then
            Unonew(np)=max(coef_flnp/rbc_2*Unonew(nx),
     &               coef_fln/rbc*Unonew(np))
           else
            if (Unonew(nx).ge.0.) then
              Unonew(np)=coef_flnp/rbc_2*Unonew(nx)
            else
              Unonew(np)=max(.5*coef_flnp/rbc_2*Unonew(nx),
     &                 coef_fln/rbc*Unonew(np))
             endif
           endif
      Velnoc(np)=Cmi*(Unonew(np)+vtrans)
      Unonew(np)=2.*Unonew(np)-Unonew(nx)

          Velic(1)=Velic(2)
          rbc=min(1.,N2new(nx-1)/N2new(nx-2))
          rbc_2=sqrt(rbc)
          N2new(np)=rbc*N2new(nx)
          U2new(np)=.5*(U2new(nx)+U2new(nx-1))
          if (U2new(np).ge.0.) then
            U2new(np)=max(coef_flnp/rbc_2*U2new(nx),
     &               coef_fln/rbc*U2new(np))
           else
            if (U2new(nx).ge.0.) then
              U2new(np)=coef_flnp/rbc_2*U2new(nx)
            else
              U2new(np)=max(.5*coef_flnp/rbc_2*U2new(nx),
     &                 coef_fln/rbc*U2new(np))
             endif
           endif
C      Velic(np)=(U2new(np)+vtrans)
      Velic(np)=vtrans
      U2new(np)=2.*U2new(np)-U2new(nx)

          Veljc(1)=Veljc(2)
          rbc=min(1.,N1new(nx-1)/N1new(nx-2))
          rbc_2=sqrt(rbc)
          N1new(np)=rbc*N1new(nx)
          U1new(np)=.5*(U1new(nx)+U1new(nx-1))
          if (U1new(np).ge.0.) then
            U1new(np)=max(coef_flnp/rbc_2*U1new(nx),
     &               coef_fln/rbc*U1new(np))
           else
            if (U1new(nx).ge.0.) then
              U1new(np)=coef_flnp/rbc_2*U1new(nx)
            else
              U1new(np)=max(.5*coef_flnp/rbc_2*U1new(nx),
     &                 coef_fln/rbc*U1new(np))
             endif
           endif
      Veljc(np)=Cji*(U1new(np)+vtrans)
      U1new(np)=2.*U1new(np)-U1new(nx)

          Velec(1)=Velec(2)
          rbc=min(1.,Nenew(nx-1)/Nenew(nx-2))
          rbc_2=sqrt(rbc)
          Nenew(np)=rbc*Nenew(nx)
          Uenew(np)=.5*(Uenew(nx)+Uenew(nx-1))
          if (Uenew(np).ge.0.) then
            Uenew(np)=max(coef_flnp/rbc_2*Uenew(nx),
     &               coef_fln/rbc*Uenew(np))
           else
            if (Uenew(nx).ge.0.) then
              Uenew(np)=coef_flnp/rbc_2*Uenew(nx)
            else
              Uenew(np)=max(.5*coef_flnp/rbc_2*Uenew(nx),
     &                 coef_fln/rbc*Uenew(np))
             endif
           endif
      Velec(np)=Cei*(Uenew(np)+vtrans)
      Uenew(np)=2.*Uenew(np)-Uenew(nx)

          Velnc(1)=Velnc(2)
          rbc=min(1.,N3new(nx-1)/N3new(nx-2))
          rbc_2=sqrt(rbc)
          N3new(np)=rbc*N3new(nx)
          U3new(np)=.5*(U3new(nx)+U3new(nx-1))
          if (U3new(np).ge.0.) then
            U3new(np)=max(coef_flnp/rbc_2*U3new(nx),
     &               coef_fln/rbc*U3new(np))
           else
            if (U3new(nx).ge.0.) then
              U3new(np)=coef_flnp/rbc_2*U3new(nx)
            else
              U3new(np)=max(.5*coef_flnp/rbc_2*U3new(nx),
     &                 coef_fln/rbc*U3new(np))
             endif
           endif
      Velnc(np)=Cni*(U3new(np)+vtrans)
      U3new(np)=2.*U3new(np)-U3new(nx)

          Velmc(1)=Velmc(2)
          rbc=min(1.,Nmnew(nx-1)/Nmnew(nx-2))
          rbc_2=sqrt(rbc)
          Nmnew(np)=rbc*Nmnew(nx)
          Umnew(np)=.5*(Umnew(nx)+Umnew(nx-1))
          if (Umnew(np).ge.0.) then
            Umnew(np)=max(coef_flnp/rbc_2*Umnew(nx),
     &               coef_fln/rbc*Umnew(np))
           else
            if (Umnew(nx).ge.0.) then
              Umnew(np)=coef_flnp/rbc_2*Umnew(nx)
            else
              Umnew(np)=max(.5*coef_flnp/rbc_2*Umnew(nx),
     &                 coef_fln/rbc*Umnew(np))
             endif
           endif
      Velmc(np)=Cmi*(Umnew(np)+vtrans)
      Umnew(np)=2.*Umnew(np)-Umnew(nx)

          do i=1,np
        Velnom(i)=.5*Velnoc(i)
        Velim(i)=.5*(Velic(i)+vtrans)
        Veljm(i)=.5*(Veljc(i)+vtrans)
        Velnm(i)=.5*(Velnc(i)+vtrans)
        Velmm(i)=.5*(Velmc(i)+vtrans)
    
        velie(i)=Velic(i)
        velje(i)=Veljc(i)
        velme(i)=Velmc(i)
        velne(i)=Velnc(i)
        velee(i)=Velec(i)
        Veliq(i)=Velie(i)
        Veljq(i)=Velje(i)
        Velnq(i)=Velne(i)
        Veleq(i)=Velee(i)
          enddo

c          velnoc(np)=max(velnoc(np),0.)
c          velic(np)=max(velic(np),0.)
c          veljc(np)=max(veljc(np),0.)
c          velec(np)=max(velec(np),0.)
c          velmc(np)=max(velmc(np),0.)
c          velnc(np)=max(velnc(np),0.)
c          velnom(np)=max(velnoc(np),0.)
c          velim(np)=max(velim(np),0.)
c          veljm(np)=max(veljm(np),0.)
c          velmm(np)=max(velmm(np),0.)
c          velnm(np)=max(velnm(np),0.)
c          velie(np)=max(velie(np),0.)
c          velje(np)=max(velje(np),0.)
c          velee(np)=max(velee(np),0.)
c          velme(np)=max(velme(np),0.)
c          velne(np)=max(velne(np),0.)
c          veliq(np)=max(veliq(np),0.)
c          veljq(np)=max(veljq(np),0.)
c          veleq(np)=max(veleq(np),0.)
c          velnq(np)=max(velnq(np),0.)


          do i=1,nx

            Ter_1(i)=300./Tenew(i)/T_0

          enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!O(1D) continuity equation
!        call velocity(velo1dc,Ipos1,Iposnp,deltat)

        do i=1,indlim
          D3(i)=Po1d(i)/N_0*t0+phdisso2(i)/N_0*t0+
     &      1.95e-7*(300./(Tenew(i)*T_0))**.7*N6new(i)*Nenew(i)*N_0*t0+
     &      .569*(Tenew(i)*T_0)**.5*(9329.+Tenew(i)*T_0)
     &      /(51183.+Tenew(i)*T_0)**3.*
     &      exp(-22756./(Tenew(i)*T_0))*Nenew(i)*No(i)*t0
          D7(i)=-(
     &          5.85e-3+
     &          7.7e-3+
     &          2.0e-11*exp(107.8/Tn(i))*Nn2(i)+
     &          8.1e-10*(Tenew(i)*T_0/300.)**.5*Nenew(i)*N_0+
     &          2.9e-11*exp(67.5/Tn(i))*No2(i)
     &          )*t0
        enddo

        do i=1,indlim
          dexpnu=D7(i)*deltat
          if (abs(dexpnu).lt.1.e-7) then
            expnu=1.+dexpnu
            dexpnu=1.
          else
            expnu=exp(dexpnu)
            dexpnu=(expnu-1.)/dexpnu
          endif

          No1dnew(i)=dexpnu*(No1dnew(i)-No1dold(i)+D3(i)*deltat)
     &             +expnu*No1dold(i)
        enddo

!        do i=1,nx
!          No1dnew(i)=max(No1dnew(i),rmin)
!        enddo
!        No1dnew(np)=min(1.,No1dnew(nx-1)/No1dnew(nx-2))*No1dnew(nx)
!        No1dnew(nx)=min(1.,No1dnew(nx-2)/No1dnew(nx-3))*No1dnew(nx-1)
!        
!        do i=1,nx
!                xno1d(i)=No1dnew(i)+Nliminf
!        enddo
      

!-------------------------------------!
!       equations de continuite       !
!-------------------------------------!

!    ion O+

      call velocity(Veljc,Ipos1,Iposnp,deltat)
          lbc=sqrt(N1old(2)/N1old(3))
      lbc=1.
      rbc=(N1old(nx)/N1old(nx-1))
      call lcpfct(N1old,N1new,Ipos1,Iposn,
     &              lbc,0.,0.,N1new(np),.false.,1)
    
!    ion H+

      call velocity(Velic,Ipos1,Iposnp,deltat)
          lbc=sqrt(N2old(2)/N2old(3))
      lbc=1.
      rbc=sqrt(N2old(nx-1)/N2old(nx-2))
      call lcpfct(N2old,N2new,Ipos1,Iposn,
     &              lbc,0.,0.,N2new(np),.false.,1)
    
!    ion N+

      call velocity(Velnc,Ipos1,Iposnp,deltat)
          lbc=sqrt(N3old(2)/N3old(3))
      lbc=1.
      rbc=(N3old(nx)/N3old(nx-1))
      call lcpfct(N3old,N3new,Ipos1,Iposn,
     &              lbc,0.,0.,N3new(np),.false.,1)
    
!    ion N2+

      call velocity(Velmc,Ipos1,Iposnp,deltat)
          lbc=sqrt(N4old(2)/N4old(3))
      lbc=1.
      rbc=(N4old(nx)/N4old(nx-1))
      call lcpfct(N4old,N4new,Ipos1,Iposn,
     &              lbc,0.,0.,N4new(np),.false.,1)
    
!    ion NO+

          lbc=sqrt(N5old(2)/N5old(3))
      lbc=1.
      rbc=(N5old(nx)/N5old(nx-1))
      call lcpfct(N5old,N5new,Ipos1,Iposn,
     &              lbc,0.,0.,N5new(np),.false.,1)
    
!    ion O2+

          lbc=sqrt(N6old(2)/N6old(3))
      lbc=1.
      rbc=(N6old(nx)/N6old(nx-1))
      call lcpfct(N6old,N6new,Ipos1,Iposn,
     &              lbc,0.,0.,N6new(np),.false.,1)
    
      do i=1,nx

    
        Tjr=T_0*T1new(i)
        Tmr=T_0*Tmnew(i)
        Tr=(Tn(i)+Tmr)/2.
        Trh=sqrt(Tr)
        lTr=log10(Tr)

            nukO2(i)=0.449e-9*No2(i)                            *t0
            nukO (i)=0.258e-9*No (i)                            *t0
            nukN2 (i)=5.14e-11*Nn2 (i)
     &             *(1.-.069*alog10(Tr))**2*sqrt(Tr)            *t0
            nulN2(i)=0.413e-9*Nn2(i)                            *t0
            nulO (i)=0.231e-9*No (i)                            *t0
            nulO2 (i)=2.59e-11*No2 (i)
     &             *(1.-.073*alog10(Tr))**2*sqrt(Tr)            *t0
            numN2(i)=0.434e-9*Nn2(i)                            *t0
            numO (i)=0.244e-9*No (i)                            *t0
            numO2(i)=0.427e-9*No2 (i)                           *t0

        nu_omega=(nukN2(i)+nukO2(i)+nukO(i)
     &               +nulN2(i)+nulO2(i)+nulO(i)
     &               +numN2(i)+numO2(i)+numO(i))/omega(i)*omz
        coef_conv=1./(1.+nu_omega**2)
        Tperp(i)=cofterp*coef_conv*Vm_2(i)

    
        TjO2=(Tjr*32.+Tn(i)*16.)/48.+(16.*32./48.*Tperp(i))
        kjO2=TjO2/300.
        TjN2=(Tjr*28.+Tn(i)*16.)/44.+(16.*28./44.*Tperp(i))
        kjN2=TjN2/300.


        TN2O=(Tmnew(i)*T_0*16.+Tn(i)*28.)/44.+(16.*28./44.*Tperp(i))
        kN2O=TN2O/300.

        TONO=(Tjr*30.+Tn(i)*16.)/46.+(16.*30./46.*Tperp(i))
        kONO=TONO/300.

        TlO2=(32.*Tmr+16.*Tn(i))/48.+(16.*32./48.*Tperp(i))
        klO2=TlO2/300.

        TN2O2=(Tmr*32.+Tn(i)*28.)/60.+(32.*28./60.*Tperp(i))

c****************************************************************
c O+ + N2 -> NO+ + N  in an O buffer: *
c**************************************
      if (TjN2 .ge. 100. .and. TjN2 .le. 6200.) then
      ak1 = 1.248e-12 - 1.751e-13*kN2O
     &           - 5.101e-14*kN2O**2 + 1.345e-14*kN2O**3
     &           - 2.922e-16*kN2O**4
        elseif ( TjN2 .gt. 6200. .and. TjN2 .le. 22000.) then
      ak1 = -9.626e-11 + 6.994e-12*kN2O
     &           - 2.315e-14*kN2O**2
        elseif ( TjN2 .gt. 22000.) then
      ak1 = -9.626e-11 + 6.994e-12*kN2O
     &           - 2.315e-14*kN2O**2
        endif
c*****************************************
c O+ + N2 -> NO+ + N  in an N2 buffer:   *
c*****************************************
      if (TjN2 .ge. 100. .and. TjN2 .le. 5500.) then
      ak1 = 1.417e-12 - 3.223e-13*kN2O - 2.362e-14*kN2O**2
     &         + 1.247e-14*kN2O**3 - 3.030e-16*kN2O**4    
      elseif(TjN2 .gt. 5500. .and. TjN2 .le. 29000.) then
      ak1 =-4.74e-11 + 3.74e-12*kN2O + 2.80e-14*kN2O**2
     &         - 1.88e-16*kN2O**3
        elseif ( TjN2 .gt. 29000.) then
      ak1 =-4.74e-11 + 3.74e-12*kN2O + 2.80e-14*kN2O**2
     &         - 1.88e-16*kN2O**3
        endif


c****************************************
c O+ + O2 -> O2+ + O in an O buffer:    *
c****************************************
      if (TjO2 .ge. 100. .and. TjO2 .le. 6400.) then
       ak2 = 2.836e-11 - 7.521e-12*kjO2 + 1.039e-12*kjO2**2
     &            - 4.981e-14*kjO2**3 + 9.087e-16*kjO2**4
      elseif (TjO2 .gt. 6400. .and. TjO2 .le. 22000.) then
       ak2 =-3.42e-11 + 4.08e-12*kjO2 - 1.70e-14*kjO2**2
        elseif (TjO2 .gt. 22000.) then
           ak2 =-3.42e-11 + 4.08e-12*kjO2 - 1.70e-14*kjO2**2
        endif
c****************************************
c O+ + O2 -> O2+ + O in an N2 buffer:    *
c****************************************
      if (TjO2 .ge. 100. .and. TjO2 .le. 8400. ) then
       ak2 = 2.763e-11 - 6.733e-12*kjO2 + 8.383e-13*kjO2**2
     &            - 3.317e-14*kjO2**3 + 4.805e-16*kjO2**4
      elseif (TjO2 .gt. 8400. .and. TjO2 .le. 31000.) then
       ak2 =-2.57e-11 + 3.48e-12*kjO2 - 1.01e-14*kjO2**2
        elseif (TjO2 .gt. 31000.) then
           ak2 =-2.57e-11 + 3.48e-12*kjO2 - 1.01e-14*kjO2**2
      endif


c*************************
c  N2+ + O -> O+ + N2    *
c*************************
        if (TN2O.le.1500.) then
          ak3=1.e-11*kN2O**(-.23)
          ak5=1.4e-10*kN2O**(-.44)
        else
          ak3=3.6e-12*kN2O**.41
          ak5=5.2e-11*kN2O**.2
        endif


c****************************************
c O+ + NO -> NO+ + O in an  O buffer:   *
c****************************************
      if (TONO .gt. 100. .and. TONO .le. 6300. ) then
      ak4 = 5.974e-13 - 9.422e-14*kONO + 6.583e-14*kONO**2
     &           - 2.156e-15*kONO**3 + 3.957e-17*kONO**4
        elseif( TONO .gt. 6300. .and. TONO .le. 22000.) then
      ak4 =-1.557e-11 + 1.397e-12*kONO + 2.461e-15*kONO**2
        elseif( TONO .gt. 22000.) then
          ak4 =-1.557e-11 + 1.397e-12*kONO + 2.461e-15*kONO**2
      endif
c****************************************
c O+ + NO -> NO+ + O in an  N2 buffer:   *
c****************************************
      if (TONO .gt. 100. .and. TONO .le. 8200. ) then
      ak4 = 5.622e-13 - 6.094e-14*kONO + 5.74e-14*kONO**2
     &           - 1.399e-15*kONO**3 + 1.84e-17*kONO**4
        elseif( TONO .gt. 8200. .and. TONO .le. 30000.) then
      ak4 =-2.22e-11 + 1.64e-12*kONO - 6.7e-17*kONO**2
        elseif( TONO .gt. 30000.) then
          ak4 =-2.22e-11 + 1.64e-12*kONO - 6.7e-17*kONO**2
      endif
          
c****************************************
c  O+ + O2 -> O2+ + O in an N2 buffer:  *
c****************************************
        if( TlO2 .ge. 100. .and. TlO2 .le. 8400. ) then
          akl2 = 2.763e-11 - 6.733e-12*klO2 + 8.383e-13*klO2**2
     &               - 3.317e-14*klO2**3 + 4.805e-16*klO2**4
        elseif( TlO2 .gt. 8400. .and. TlO2 .le. 31000. ) then
          akl2 =-2.57e-11 + 3.48e-12*klO2 - 1.01e-14*klO2**2
            elseif( TlO2 .gt. 31000.) then
              akl2 =-2.57e-11 + 3.48e-12*klO2 - 1.01e-14*klO2**2
        endif

          
              y0  (1  )=N1old(i)
              y1  (1  )=(N1new(i)-N1old(i))/deltat+Po(i)/N_0*t0
              chim(1,1)=-(ak1*Nn2(i)+ak2*No2(i)+ak4*N_0*Nnonew(i)
     &                +2.5e-11*sqrt(Tn(i)+T1new(i)*T_0/16.
     &                +1.2e-8*(U1new(i)*Cj0)**2)*Nh(i))        *t0
     &                -3.*Cji*U1new(i)*alt_geo_1(i)
              chim(1,2)=2.2e-11*sqrt(T_0*T2new(i)+Tn(i)/16.
     &             +1.2e-8*(U2new(i)*Ci0)**2)*No(i)*t0
              chim(1,3)=(3.7e-11*No2(i)+5.e-13*No(i))*t0
              chim(1,4)=ak3*No(i)*t0
              chim(1,5)=0.d0
              chim(1,6)=0.d0
          
              y0  (2  )=N2old(i)
              y1  (2  )=(N2new(i)-N2old(i))/deltat+Ph(i)/N_0*t0
              chim(2,1)=2.5d-11*sqrt(Tn(i)+T1new(i)*T_0/16.+1.2e-8*
     &               (U1new(i)*Cj0)**2)*Nh(i)*t0
              chim(2,2)=-2.2d-11*sqrt(T_0*T2new(i)+Tn(i)/16.
     &              +1.2d-8*(U2new(i)*Ci0)**2)*No(i)        *t0
     &              -3.*U2new(i)*alt_geo_1(i)
              chim(2,3)=3.6d-12*Nh(i)*t0
              chim(2,4)=0.d0
              chim(2,5)=0.d0
              chim(2,6)=0.d0
          
              y0  (3  )=N3old(i)
              y1  (3  )=(N3new(i)-N3old(i))/deltat+0.21*Pn2(i)*t0/N_0
              chim(3,1)=0.d0
              chim(3,2)=0.d0
              chim(3,3)=-(2.6e-10*No2(i)+3.1e-10*No2(i)+3.7e-11*No2(i)
     &                +2.e-11*N_0*Nnonew(i)+3.6e-12*Nh(i)
     &                +5.e-13*No(i))*t0
     &              -3*Cni*U3new(i)*alt_geo_1(i)
              chim(3,4)=0.d0
              chim(3,5)=0.d0
              chim(3,6)=0.d0
          
              y0  (4  )=N4old(i)
              y1  (4  )=(N4new(i)-N4old(i))/deltat+Pn2(i)*0.79*t0/N_0          
        chim(4,1)=0.d0
              chim(4,2)=0.d0
              chim(4,3)=0.d0
             chim(4,4)=-(5.e-11*(300./TN2O2)*No2(i)+ak3*No(i)+ak5*No(i)
     &                +3.3e-10*N_0*Nnonew(i)
     &                +1.8e-7*Ter_1(i)**.39*Nenew(i)*N_0)*t0
     &              -3.*Czi*Umnew(i)*alt_geo_1(i)
              chim(4,5)=0.d0
              chim(4,6)=0.d0
          
              y0  (5  )=N5old(i)
              y1  (5  )=(N5new(i)-N5old(i))/deltat+6.e-7*Nnonew(i)*t0          
              chim(5,1)=(ak4*N_0*Nnonew(i)+ak1*Nn2(i))*t0
              chim(5,2)=0.d0
              chim(5,3)=(No2(i)*2.6e-10+2.e-11*N_0*Nnonew(i))*t0
              chim(5,4)=ak5*No(i)*t0
              chim(5,5)=-4.2e-7*Ter_1(i)**.85*Nenew(i)*N_0*t0
     &              -3.*Czi*Umnew(i)*alt_geo_1(i)
              chim(5,6)=(5.e-16*Nn2(i)+1.2e-10*Nn(i)
     &                 +4.5e-10*N_0*Nnonew(i))*t0
          
              y0  (6  )=N6old(i)
              y1  (6  )=(N6new(i)-N6old(i))/deltat+Po2(i)/N_0*t0
              chim(6,1)=akl2*No2(i)*t0
              chim(6,2)=0.d0
              chim(6,3)=3.1e-10*No2(i)*t0
              chim(6,4)=5.e-11*(300./TN2O2)*No2(i)*t0
              chim(6,5)=0.d0
              chim(6,6)=-(5.e-16*Nn2(i)+1.2e-10*Nn(i)
     &              +4.5e-10*N_0*Nnonew(i)
     &                +1.6e-7*Ter_1(i)**.55*Nenew(i)*N_0)*t0
     &              -3.*Czi*Umnew(i)*alt_geo_1(i)

        call solve_ODE(nb_ion,y0,y1,chim,mat,deltat)
    
        N1new(i)=max(y0(1),r_min)
        N2new(i)=max(y0(2),r_min)
        N3new(i)=max(y0(3),r_min)
        N4new(i)=max(y0(4),r_min)
        N5new(i)=max(y0(5),r_min)
        N6new(i)=max(y0(6),r_min)

       enddo
    
          N1new(np)=min(1.,N1new(nx-1)/N1new(nx-2))*N1new(nx)
          N2new(np)=min(1.,N2new(nx-1)/N2new(nx-2))*N2new(nx)
          N3new(np)=min(1.,N3new(nx-1)/N3new(nx-2))*N3new(nx)
          N4new(np)=min(1.,N4new(nx-1)/N4new(nx-2))*N4new(nx)
          N5new(np)=min(1.,N5new(nx-1)/N5new(nx-2))*N5new(nx)
          N6new(np)=min(1.,N6new(nx-1)/N6new(nx-2))*N6new(nx)

      if (any(ieee_is_nan(N1new))) then
        write(stderr,*) 'problem N1new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N2new))) then
        write(stderr,*) 'problem  N2new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N3new))) then
        write(stderr,*) 'problem N3new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N4new))) then
        write(stderr,*) 'problem  N4new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N5new))) then
        write(stderr,*) 'problem  N5new dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(N6new))) then
        write(stderr,*) 'problem  N6new dans la boucle 1'
        goto 246
      endif

          do  i=1,nx
            Nenew(i)=N2new(i)+N1new(i)+N4new(i)
     &              +N6new(i)+N5new(i)+N3new(i)
             if (alt(i).le.500.) then
          xne(i)=Nenew(i)+Nliminf
          xn1(i)=N1new(i)+Nliminf
          xn2(i)=N2new(i)+Nliminf
          xn3(i)=N3new(i)+Nliminf
          xn4(i)=N4new(i)+Nliminf
          xn5(i)=N5new(i)+Nliminf
          xn6(i)=N6new(i)+Nliminf
          xnm(i)=Nmnew(i)+Nliminf
    
          xne_1(i)=Nenew(i)/xne(i)/xne(i)
          xn1_1(i)=N1new(i)/xn1(i)/xn1(i)
          xn2_1(i)=N2new(i)/xn2(i)/xn2(i)
          xn3_1(i)=N3new(i)/xn3(i)/xn3(i)
          xn4_1(i)=N4new(i)/xn4(i)/xn4(i)
          xn5_1(i)=N5new(i)/xn5(i)/xn5(i)
          xn6_1(i)=N6new(i)/xn6(i)/xn6(i)
          xnm_1(i)=Nmnew(i)/xnm(i)/xnm(i)
        else
          xne(i)=Nenew(i)
          xn1(i)=N1new(i)
          xn2(i)=N2new(i)
          xn3(i)=N3new(i)
          xn4(i)=N4new(i)
          xn5(i)=N5new(i)
          xn6(i)=N6new(i)
          xnm(i)=Nmnew(i)
    
          xne_1(i)=1./Nenew(i)
          xn1_1(i)=1./N1new(i)
          xn2_1(i)=1./N2new(i)
          xn3_1(i)=1./N3new(i)
          xn4_1(i)=1./N4new(i)
          xn5_1(i)=1./N5new(i)
          xn6_1(i)=1./N6new(i)
          xnm_1(i)=1./Nmnew(i)
        endif
           enddo
      xne(np)=Nenew(np)
      xn1(np)=N1new(np)
      xn2(np)=N2new(np)
      xn3(np)=N3new(np)
      xn4(np)=N4new(np)
      xn5(np)=N5new(np)
      xn6(np)=N6new(np)
      xnm(np)=Nmnew(np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!O(1D) momentum equation
        call velocity(velo1dm,Ipos1,Iposnp,deltat)

        do i=1,nx
          diffo1d=1./(No(i)/5.37e17/sqrt(Tn(i))+
     &          Nn2(i)/4.76e17/sqrt(Tn(i))+
     &          No2(i)/6.58e17/sqrt(Tn(i)))
          nuo1d=1.38e-17*Tn(i)/(16.*1.66e-27)/diffo1d*t0

          Tr=(Tn(i)+T_0*T1new(i))/2.
          nuo1dj(i)=N_0*N1new(i)*3.67e-11*sqrt(Tr)
     &               *(1.-0.064*alog10(Tr))**2                 *t0
     
          nuo1di(i)=N_0*N2new(i)*1/16*6.61e-11*
     &       sqrt(T_0*T2new(i))*(1.-0.047*
     &       alog10(T_0*T2new(i)))**2                          *t0
          
          nuo1dn(i)=N_0*N3new(i)*14./16.*4.42e-10              *t0
          
          nuo1dk(i)=N_0*N4new(i)*28./16.*2.58e-10              *t0
          
          nuo1dm(i)=N_0*N5new(i)*30./26.*2.44e-10              *t0
          
          nuo1dl(i)=N_0*N6new(i)*2.*2.31e-10                   *t0
          
          nuo1de(i)=N_0*Nenew(i)*8.9e-11*3.428e-5
     &          *(1.+5.7e-4*T_0*Tenew(i))
     &          *sqrt(T_0*Tenew(i))                            *t0

          thermacc(i)=thermodiff*qenew(i)/Tenew(i)**2.5

          D3(i)=-Cij*G(i)+
     &           nuo1d*Un(i)/Cj0+
     &           nuo1dj(i)*U1new(i)+nuo1di(i)*Cij*U2new(i)+
     &           (nuo1dk(i)+nuo1dl(i)+nuo1dm(i))*Czj*Umnew(i)+
     &          nuo1dn(i)*Cnj*U3new(i)+nuo1de(i)*Cej*Uenew(i)
          D7(i)=-(nuo1d+
     &          nuo1dj(i)+nuo1di(i)+nuo1dk(i)+nuo1dl(i)+
     &          nuo1dm(i)+nuo1dn(i)+nuo1de(i))
        
          C2a(i)=-Cji*Tn(i)/T_0*No1dnew(i)/xno1d(i)
          D2a(i)=log(xno1d(i)*Tn(i)/T_0)
        enddo
        
        D2al=(D2a(1)+D2a(2))/2.
        D2ar=(xno1d(np)*Tn(np)/T_0*xno1d(nx)*Tn(np)/T_0)
        call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
        call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
        
        lbc=1.
        call lcpfct(Uo1dold,Uo1dnew,ipos1,iposnp,
     &           lbc,0.,0.,Uo1dnew(np),.false.,0)

        do i=1,nx
          dexpnu=D7(i)*deltat
          if (abs(dexpnu).lt.1.e-7) then
            expnu=1.+dexpnu
            dexpnu=1.
          else
            expnu=exp(dexpnu)
            dexpnu=(expnu-1.)/dexpnu
          endif

          Uo1dnew(i)=dexpnu*(Uo1dnew(i)-Uo1dold(i)+D3(i)*deltat)
     &             +expnu*Uo1dold(i)
        Uo1dnew(i)=0.
        enddo
        Uo1dnew(nx)=0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C[[[    H+ momentum equation resolution

      call velocity(Velim,Ipos1,Iposnp,deltat)

          do i=1,nx

          Tr1=T_0*T2new(i)
          Tr=(Tn(i)+Tr1)/2.

      T1_15(i)=T1new(i)**(-1.5)/nu_0
      T2_15(i)=T1new(i)**(-1.5)/nu_0
      T3_15(i)=T3new(i)**(-1.5)/nu_0
      Tm_15(i)=Tmnew(i)**(-1.5)/nu_0
      Te_15(i)=Tenew(i)**(-1.5)/nu_0


      thermacc(i)=thermodiff*qenew(i)/Tenew(i)**2.5

          nuiN2(i)=3.36e-9*Nn2(i)                               *t0

          nuiO2(i)=3.20e-9*No2(i)                               *t0
          nuiO (i)=6.61e-11*No (i)
     &             *(1.-.047*alog10(Tr1))**2*sqrt(Tr1)          *t0
          nuiOHot(i)=6.61e-11*NOHot(i)                        !MZ
     &             *(1.-.047*alog10(Tr1))**2*sqrt(Tr1)          *t0        !MZ
          nuiH (i)=2.65e-10*Nh (i)
     &             *(1.-.083*alog10(Tr))**2*sqrt(Tr)            *t0
          nuij (i) = 1.23*N1new(i)
     &              /((16.*T2new(i)+T1new(i))/17.)**1.5        /nu_0
          nuii (i) = .9*N2new(i)
     &              /(T2new(i))**1.5                 /nu_0
          nuin (i) = 1.23*N3new(i)
     &              /((14.*T2new(i)+T3new(i))/15.)**1.5     /nu_0
          nuik (i) = 1.25*N4new(i)
     &              /((28.*T2new(i)+Tmnew(i))/27.)**1.5        /nu_0
          nuim (i) = 1.25*N5new(i)
     &              /((30.*T2new(i)+Tmnew(i))/31.)**1.5      /nu_0
          nuil (i) = 1.25*N6new(i)
     &              /((32.*T2new(i)+Tmnew(i))/33.)**1.5        /nu_0
          nuie (i)=0.03    *Nenew(i)*Te_15(i)


            C2a(i)=-Tepnew(i)
        D2a(i)=log(Nenew(i)*Tepnew(i))
            C2b(i)=-T2pnew(i)*N2new(i)/xn2(i)
            D2b(i)=log(xn2(i)*T2pnew(i))
            D3(i)=-G(i)+thermacc(i)/mi/Ci0
     &            -3.*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*(T2pnew(i)-T2tnew(i))*alt_geo_1(i)
     &            +nuij(i)*Cji*U1new(i)
     &            +(nuik(i)*Czi+nuil(i)*Czi+nuim(i)*Czi)*Umnew(i)
     &            +nuin(i)*Cni*U3new(i)

            D3(i)=D3(i)+.6*q2new(i)*xn2_1(i)
     &                 *(28.*nuiN2(i)/(Tn(i)/T_0+28.*T2new(i))
     &                  +32.*nuiO2(i)/(Tn(i)/T_0+32.*T2new(i))
     &                  +16.*nuiO(i)/(Tn(i)/T_0+16.*T2new(i))
     &                  +16.*nuiOHot(i)/(TnOHot(i)/T_0+16.*T2new(i))        !MZ
     &                  +    nuiH(i)/(Tn(i)/T_0+T2new(i))
     &                  +16.*nuij(i)/(T1new(i)+16.*T2new(i))
     &                  +28.*nuik(i)/(Tmnew(i)+28.*T2new(i))
     &                  +32.*nuil(i)/(Tmnew(i)+32.*T2new(i))
     &                  +30.*nuim(i)/(Tmnew(i)+30.*T2new(i))
     &                  +14.*nuin(i)/(T3new(i)+14.*T2new(i))
     &                  +5.46e-4*nuie(i)/(Tenew(i)+5.46e-4*T2new(i)))

            D3(i)=D3(i)-.6*q1new(i)*xn1_1(i)*Cji
     &                    *nuij(i)/(T1new(i)+16.*T2new(i))

            D3(i)=D3(i)-.6*qenew(i)*xne_1(i)*Cei
     &                    *nuie(i)/(Tenew(i)+5.46e-4*T2new(i))

            D3(i)=D3(i)-.6*q3new(i)*xn3_1(i)*Cni
     &                    *nuin(i)/(T3new(i)+14.*T2new(i))

            D3(i)=D3(i)+(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+nuiH(i))        !MZ
     &                    *(Un(i)/Ci0)
     &                 -.6*nuiN2(i)/(Tn(i)/T_0+28.*T2new(i))
     &                    *(q_Nn2(i)*Niqi0)
     &                 -.6*nuiO2(i)/(Tn(i)/T_0+32.*T2new(i))
     &                    *(q_No2(i)*Niqi0)
     &                 -.6*nuiO(i)/(Tn(i)/T_0+16.*T2new(i))
     &                    *(q_No(i)*Niqi0)
     &                 -.6*nuiOHot(i)/(TnOHot(i)/T_0+16.*T2new(i))        !MZ
     &                    *(q_NOHot(i)*Niqi0)                    !MZ
     &                 -.6*nuiH(i)/(Tn(i)/T_0+T2new(i))
     &                    *(q_Nh(i)*Niqi0)
     &                 +nuie(i)*Cei*Uenew(i)

            D7(i)=-(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+nuiH(i)        !MZ
     &              +nuij(i)+nuik(i)+nuil(i)+nuim(i)
     &              +nuin(i)+nuie(i))

          enddo
          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
      D2br=.5*log(xn2(np)*T2pnew(np)*xn2(nx)*T2pnew(nx))

      D2br=log(xn2(nx)*T2pnew(nx))                        !MZ
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U2old,extra,nx)
      rbc=1.
      call lcpfct(U2old,U2new,Ipos1,Iposn,
     &              lbc,0.,0.,U2new(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            U2new(i)=dexpnu*(U2new(i)-U2old(i)+D3(i)*deltat)
     &             +expnu*U2old(i)
          enddo

      if (any(ieee_is_nan(U2new))) then
        write(stderr,*) 'problem U2new dans la boucle 2'
        goto 246
      endif


C]]]

C[[[    O+ momentum equation resolution

      call velocity(Veljm,Ipos1,Iposnp,deltat)

          do i=1,nx

          Tr=(T_0*T1new(i)+Tn(i))/2.
      TrHot=(T_0*T1new(i)+TnOHot(i))/2.                    !MZ

          nujN2(i)=6.82e-10*Nn2(i)                              *t0
          nujO2(i)=6.64e-10*No2(i)                              *t0
          nujO (i)=3.67e-11*No(i)*sqrt(Tr)
     &             *(1.-.064*alog10(Tr))**2                     *t0
          nujOHot(i)=3.67e-11*NOHot(i)*sqrt(TrHot)                !MZ
     &             *(1.-.064*alog10(TrHot))**2                  *t0        !MZ
          nujj(i)=.22*N1new(i)
     &           /(T1new(i))**1.5                /nu_0
          nuji(i)=0.077*N2new(i)
     &           /((T1new(i)+16.*T2new(i))/17.)**1.5        /nu_0
          nujn(i)=0.22*N3new(i)
     &           /((14.*T1new(i)+16.*T3new(i))/30.)**1.5     /nu_0
          nujk(i)=0.25*N4new(i)
     &           /((28.*T1new(i)+16.*Tmnew(i))/44.)**1.5    /nu_0
          nujm(i)=0.26*N5new(i)
     &           /((30.*T1new(i)+16.*Tmnew(i))/46.)**1.5    /nu_0
          nujl(i)=0.26*N6new(i)
     &             /((32.*T1new(i)+16.*Tmnew(i))/48.)**1.5    /nu_0
          nuje(i)=1.87e-3*Nenew(i)*Te_15(i)


            C2a(i)=-Cji*Tepnew(i)
        D2a(i)=log(Nenew(i)*Tepnew(i))
            C2b(i)=-Cji*T1pnew(i)*N1new(i)/xn1(i)
            D2b(i)=log(xn1(i)*T1pnew(i))
            D3(i)=-Cij*G(i)+thermacc(i)/mj/Cj0
     &            -3.*Cji*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*Cji*(T1pnew(i)-T1tnew(i))*alt_geo_1(i)
     &            +nuji(i)*Cij*U2new(i)

            D3(i)=D3(i)+(nujk(i)+nujl(i)+nujm(i))*Czj*Umnew(i)
     &                 +nujn(i)*Cnj*U3new(i)

            D3(i)=D3(i)-.6*q2new(i)*xn2_1(i)*Cij*16.
     &                  *nuji(i)/(16.*T2new(i)+T1new(i))

            D3(i)=D3(i)+.6*q1new(i)*xn1_1(i)
     &                  *(28.*nujN2(i)/(16.*Tn(i)/T_0+28.*T1new(i))
     &                   +32.*nujO2(i)/(16.*Tn(i)/T_0+32.*T1new(i))
     &                   +16.*nujO(i)/(16.*Tn(i)/T_0+16.*T1new(i))
     &             +16.*nujOHot(i)/(16.*TnOHot(i)/T_0+16.*T1new(i))        !MZ
     &                   +    nuji(i)/(16.*T2new(i)+T1new(i))
     &                   +28.*nujk(i)/(16.*Tmnew(i)+28.*T1new(i))
     &                   +32.*nujl(i)/(16.*Tmnew(i)+32.*T1new(i))
     &                   +30.*nujm(i)/(16.*Tmnew(i)+30.*T1new(i))
     &                   +14.*nujn(i)/(16.*T3new(i)+14.*T1new(i))
     &                   +    nuje(i)/(16.*Tenew(i)+5.46e-4*T1new(i))
     &                        *5.46e-4)


            D3(i)=D3(i)-.6*qenew(i)*xne_1(i)*Cej*16.
     &                  *nuje(i)/(16.*Tenew(i)+5.46e-4*T1new(i))

            D3(i)=D3(i)-.6*q3new(i)*xn3_1(i)*Cnj*16.
     &                  *nujn(i)/(16.*T3new(i)+14.*T1new(i))

            D3(i)=D3(i)+(nujN2(i)+nujO2(i)+nujO(i)
     &            +nujOHot(i))*(Un(i)/Cj0)                !MZ
     &                 -.6*nujN2(i)/(16.*Tn(i)/T_0+28.*T1new(i))
     &                  *(16.*q_Nn2(i)*Njqj0)
     &                 -.6*nujO2(i)/(16.*Tn(i)/T_0+32.*T1new(i))
     &                  *(16.*q_No2(i)*Njqj0)
     &                 -.6*nujO(i)/(16.*Tn(i)/T_0+16.*T1new(i))
     &                  *(16.*q_No(i)*Njqj0)
     &                 -.6*nujOHot(i)/(16.*TnOHot(i)/T_0+16.*T1new(i))        !MZ
     &                  *(16.*q_NOHot(i)*Njqj0)                    !MZ
     &                 +nuje(i)*Cje*Uenew(i)

            D7(i)=-(nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i)            !MZ
     &              +nuji(i)+nujn(i)
     &              +nuje(i)+nujk(i)+nujl(i)+nujm(i))

          enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
      D2br=.5*log(xn1(np)*T1pnew(np)*xn1(nx)*T1pnew(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U1old,extra,nx)
      rbc=1.
      call lcpfct(U1old,U1new,Ipos1,Iposn,
     &              lbc,0.,0.,U1new(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            U1new(i)=dexpnu*(U1new(i)-U1old(i)+D3(i)*deltat)
     &             +expnu*U1old(i)
          enddo

      if (any(ieee_is_nan(U1new))) then
        write(stderr,*) 'problem U1new dans la boucle 2'
        goto 246
      endif


C]]]

C[[[    heavy ions momentum equation resolution

      call velocity(Velmm,Ipos1,Iposnp,deltat)

          do i=1,nx

          Tr=(T_0*Tmnew(i)+Tn(i))/2.

          nukH(i)=0.074e-9*Nh(i)                                *t0
          nukO2(i)=0.449e-9*No2(i)                              *t0
          nukO (i)=0.258e-9*No (i)                              *t0
          nukOHot(i)=0.258e-9*NOHot(i)                          *t0        !MZ
          nukN2 (i)=5.14e-11*Nn2 (i)
     &             *(1.-.069*alog10(Tr))**2*sqrt(Tr)            *t0

          nulH(i)=0.065e-9*Nh(i)                                *t0
          nulN2(i)=0.413e-9*Nn2(i)                              *t0
          nulO (i)=0.231e-9*No (i)                              *t0
          nulOHot(i)=0.231e-9*NOHot(i)                          *t0        !MZ
          nulO2 (i)=2.59e-11*No2 (i)
     &             *(1.-.073*alog10(Tr))**2*sqrt(Tr)            *t0

          numH(i)=0.069e-9*Nh(i)                                *t0
          numN2(i)=0.434e-9*Nn2(i)                              *t0
          numO (i)=0.244e-9*No (i)                              *t0
          numOHot(i)=0.244e-9*NOHot(i)                          *t0        !MZ
          numO2(i)=0.427e-9*No2 (i)                             *t0

          nukj(i)=0.15*N1new(i)
     &            / ((16.*Tmnew(i)+28.*T1new(i))/44.)**1.5    /nu_0
          nuki(i)=0.045*N2new(i)
     &            / ((Tmnew(i)+28.*T2new(i))/29.)**1.5        /nu_0
          nukn(i)=0.14*N3new(i)
     &            / ((14.*Tmnew(i)+28.*T3new(i))/42.)**1.5    /nu_0
          nukk(i)=0.17*N4new(i)
     &            /(Tmnew(i))**1.5                /nu_0
          nukm(i)=0.17*N5new(i)
     &            /(Tmnew(i))**1.5                /nu_0
          nukl(i)=0.18*N6new(i)
     &            /(Tmnew(i))**1.5                /nu_0
          nuke(i)=1.07e-3*Nenew(i)/(Tenew(i))**1.5        /nu_0

          nulj(i)=0.13*N1new(i)
     &            / ((16.*Tmnew(i)+32.*T1new(i))/48.)**1.5    /nu_0
          nuli(i)=0.039*N2new(i)
     &            / ((Tmnew(i)+32.*T2new(i))/33.)**1.5        /nu_0
          nuln(i)=0.12*N3new(i)
     &            / ((14.*Tmnew(i)+32.*T3new(i))/46)**1.5    /nu_0
          nulk(i)=0.15*N4new(i)/(Tmnew(i))**1.5        /nu_0
          nulm(i)=0.16*N5new(i)/(Tmnew(i))**1.5        /nu_0
          null(i)=0.16*N6new(i)/(Tmnew(i))**1.5        /nu_0
          nule(i)=9.347e-4*Nenew(i)/(Tenew(i))**1.5    /nu_0
          numj(i)=0.14*N1new(i)
     &            / ((16.*Tmnew(i)+30.*T1new(i))/46.)**1.5    /nu_0
          numi(i)=0.042*N2new(i)
     &            /((Tmnew(i)+30.*T2new(i))/31.)**1.5        /nu_0
          numn(i)=0.13*N3new(i)
     &            /((14.*Tmnew(i)+30.*T3new(i))/44.)**1.5    /nu_0
          numk(i)=0.16*N4new(i)/(Tmnew(i))**1.5             /nu_0
          numm(i)=0.16*N5new(i)/(Tmnew(i))**1.5            /nu_0
          numl(i)=0.17*N6new(i)/(Tmnew(i))**1.5            /nu_0
          nume(i)=9.97e-4*Nenew(i)/(Tenew(i))**1.5        /nu_0


          C2a(i)=-Czi*Tepnew(i)
          D2a(i)=log(Nenew(i)*Tepnew(i))

          C2b(i)=-Czi*Tmpnew(i)*Nmnew(i)/xnm(i)
          D2b(i)= log(xnm(i)*Tmpnew(i))

          D3(i)=-Ciz*G(i)+thermacc(i)/mz/Cz0
     &            -3.*Czi*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*Czi*(Tmpnew(i)-Tmtnew(i))*alt_geo_1(i)
     &           +(nukj(i)+nulj(i)+numj(i))*Cjz/3.*U1new(i)
     &           +(nuki(i)+nuli(i)+numi(i))*Ciz/3.*U2new(i)
     &           +(nukn(i)+nuln(i)+numn(i))*Cnz/3.*U3new(i)

          D3(i)=D3(i)-.2*q2new(i)*xn2_1(i)*Ciz
     &                *(28.*nuki(i)/(28.*T2new(i)+Tmnew(i))
     &                 +32.*nuli(i)/(32.*T2new(i)+Tmnew(i))
     &                 +30.*numi(i)/(30.*T2new(i)+Tmnew(i)))


          D3(i)=D3(i)-.2*q1new(i)*xn1_1(i)*Cjz
     &                *(28.*nukj(i)/(28.*T1new(i)+16.*Tmnew(i))
     &                 +32.*nulj(i)/(32.*T1new(i)+16.*Tmnew(i))
     &                 +30.*numj(i)/(30.*T1new(i)+16.*Tmnew(i)))


          D3(i)=D3(i)-.2*qenew(i)*xne_1(i)*Cez
     &                *(28.*nuke(i)/(28.*Tenew(i)+5.46e-4*Tmnew(i))
     &                 +32.*nule(i)/(32.*Tenew(i)+5.46e-4*Tmnew(i))
     &                 +30.*nume(i)/(30.*Tenew(i)+5.46e-4*Tmnew(i)))

          D3(i)=D3(i)-.2*q3new(i)*xn3_1(i)*Cnz
     &                *(28.*nukn(i)/(28.*T3new(i)+14.*Tmnew(i))
     &                 +32.*nuln(i)/(32.*T3new(i)+14.*Tmnew(i))
     &                 +30.*numn(i)/(30.*T3new(i)+14.*Tmnew(i)))


          D3(i)=D3(i)+(nukN2(i)+nukO2(i)+nukO(i)+nukOHot(i)            !MZ
     &               +nulN2(i)+nulO2(i)+nulO(i)+nulOHot(i)            !MZ
     &               +numN2(i)+numO2(i)+numO(i)+numOHot(i))*Un(i)/Cz0/3.    !MZ
     &           -.2*
     &               (nukN2(i)/(28.*Tn(i)/T_0+28.*Tmnew(i))
     &           *(28.*q_Nn2(i)*Nkqk0)
     &               +nukO2(i)/(28.*Tn(i)/T_0+32.*Tmnew(i))
     &           *(28.*q_No2(i)*Nkqk0)
     &               +nukO(i)/(28.*Tn(i)/T_0+16.*Tmnew(i))
     &           *(28.*q_No(i)*Nkqk0)
     &               +nukOHot(i)/(28.*TnOHot(i)/T_0+16.*Tmnew(i))        !MZ
     &           *(28.*q_NOHot(i)*Nkqk0))                    !MZ
     &           -.2*
     &               (nulN2(i)/(32.*Tn(i)/T_0+28.*Tmnew(i))
     &           *(32.*q_Nn2(i)*Nlql0)
     &               +nulO2(i)/(32.*Tn(i)/T_0+32.*Tmnew(i))
     &           *(32.*q_No2(i)*Nlql0)
     &               +nulO(i)/(32.*Tn(i)/T_0+16.*Tmnew(i))
     &           *(32.*q_No(i)*Nlql0)
     &               +nulOHot(i)/(32.*TnOhot(i)/T_0+16.*Tmnew(i))        !MZ
     &           *(32.*q_NOHot(i)*Nlql0))                    !MZ
     &           -.2*
     &               (numN2(i)/(30.*Tn(i)/T_0+28.*Tmnew(i))
     &           *(30.*q_Nn2(i)*Nmqm0)
     &                +numO2(i)/(30.*Tn(i)/T_0+32.*Tmnew(i))
     &           *(30.*q_No2(i)*Nmqm0)
     &                +numO(i)/(30.*Tn(i)/T_0+16.*Tmnew(i))
     &           *(30.*q_No(i)*Nmqm0)
     &                +numOHot(i)/(30.*TnOHot(i)/T_0+16.*Tmnew(i))        !MZ
     &           *(30.*q_NOHot(i)*Nmqm0))                    !MZ
     &           +Uenew(i)*Cez*(nuke(i)+nule(i)+nume(i))/3.

          D7(i)=-(nukN2(i)+nukO2(i)+nukO(i)+nukOHot(i)                !MZ
     &            +nuki(i)+nuke(i)+nukj(i)+nukn(i)
     &            +nulN2(i)+nulO2(i)+nulO(i)+nulOHot(i)                !MZ
     &            +nuli(i)+nule(i)+nulj(i)+nuln(i)
     &            +numN2(i)+numO2(i)+numO(i)+numOHot(i)                !MZ
     &            +numi(i)+nume(i)+numj(i)+numn(i))/3.

          enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
      D2br=.5*log(xnm(np)*Tmpnew(np)*xnm(nx)*Tmpnew(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U1old,extra,nx)
      rbc=1.
      call lcpfct(Umold,Umnew,Ipos1,Iposn,
     &              lbc,0.,0.,Umnew(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Umnew(i)=dexpnu*(Umnew(i)-Umold(i)+D3(i)*deltat)
     &             +expnu*Umold(i)
          enddo

      if (any(ieee_is_nan(Umnew))) then
        write(stderr,*) 'problem  Umnew dans la boucle 2'
        goto 246
      endif


C[[[[   Boundaries conditions

C]]]

C[[[    N+ momentum equation resolution
    
      call velocity(Velnm,Ipos1,Iposnp,deltat)

          do i=1,nx

          nunH (i)=0.145e-9*Nh(i)                               *t0
          nunN2(i)=0.747e-9*Nn2(i)                              *t0
          nunO (i)=0.442e-9*No (i)                              *t0
          nunOHot(i)=0.442e-9*NOHot(i)                          *t0        !MZ
          nunO2(i)=0.725e-9*No2 (i)                             *t0
          nunj(i)=0.25*N1new(i)
     &           / ((16.*T3new(i)+14.*T1new(i))/30.)**1.5    /nu_0
          nuni(i)=0.088*N2new(i)
     &            / (( T3new(i)+14.*T2new(i))/15.)**1.5        /nu_0
          nunn(i)=0.24*N3new(i)/ (T3new(i))**1.5        /nu_0
          nunk(i)=0.28*N4new(i)
     &            / ((28.*T3new(i)+14.*Tmnew(i))/42.)**1.5    /nu_0
          nunm(i)=0.28*N5new(i)
     &            / ((30.*T3new(i)+14.*Tmnew(i))/44.)**1.5    /nu_0
          nunl(i)=0.28*N6new(i)
     &            / ((32.*T3new(i)+14.*Tmnew(i))/46.)**1.5    /nu_0
          nune(i)=2.136e-3*Nenew(i)*Te_15(i)


            C2a(i)=-Cni*Tepnew(i)
        D2a(i)=log(Nenew(i)*Tepnew(i))
            C2b(i)=-Cni*T3pnew(i)*N3new(i)/xn3(i)
            D2b(i)=log(xn3(i)*T3pnew(i))

            D3(i)=-Cin*G(i)+thermacc(i)/mn/Cn0
     &            -3.*Cni*(Tepnew(i)-Tetnew(i))*alt_geo_1(i)
     &            -3.*Cni*(T3pnew(i)-T3tnew(i))*alt_geo_1(i)
     &            +nuni(i)*Cin*U2new(i)+nunj(i)*Cjn*U1new(i)
     &            +(nunk(i)+nunl(i)+nunm(i))*Czn*Umnew(i)

            D7(i)=-(nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i)            !MZ
     &             +nuni(i)+nunj(i)+nune(i)
     &             +nunk(i)+nunl(i)+nunm(i))


          D3(i)=D3(i)-.6*q2new(i)*xn2_1(i)*Cin*14.
     &                  *nuni(i)/(14.*T2new(i)+T3new(i))

          D3(i)=D3(i)-.6*q1new(i)*xn1_1(i)*Cjn*14.
     &                  *nunj(i)/(14.*T1new(i)+16.*T3new(i))

          D3(i)=D3(i)-.6*qenew(i)*xne_1(i)*Cen*14.
     &                  *nune(i)/(14.*Tenew(i)+5.46e-4*T3new(i))

          D3(i)=D3(i)+.6*q3new(i)*xn3_1(i)
     &                  *(28.*nunN2(i)/(14.*Tn(i)/T_0+28.*T3new(i))
     &                  +32.*nunO2(i)/(14.*Tn(i)/T_0+32.*T3new(i))
     &                  +16.*nunO(i)/(14.*Tn(i)/T_0+16.*T3new(i))
     &                  +16.*nunOHot(i)/(14.*TnOHot(i)/T_0+16.*T3new(i))    !MZ
     &                  +    nuni(i)/(14.*T2new(i)+T3new(i))
     &                  +28.*nunk(i)/(14.*Tmnew(i)+28.*T3new(i))
     &                  +32.*nunl(i)/(14.*Tmnew(i)+32.*T3new(i))
     &                  +30.*nunm(i)/(14.*Tmnew(i)+30.*T3new(i))
     &                  +    nune(i)/(14.*Tenew(i)+5.46e-4*T3new(i))
     &                       *5.46e-4)


          D3(i)=D3(i)+(nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i))*Un(i)/Cn0        !MZ
     &           -.6*nunN2(i)/(14.*Tn(i)/T_0+28.*T3new(i))
     &           *(14.*q_Nn2(i)*Nnqn0)
     &           -.6*nunO2(i)/(14.*Tn(i)/T_0+32.*T3new(i))
     &           *(14.*q_No2(i)*Nnqn0)
     &           -.6*nunO(i)/(14.*Tn(i)/T_0+16.*T3new(i))
     &           *(14.*q_No(i)*Nnqn0)
     &           -.6*nunOHot(i)/(14.*TnOHot(i)/T_0+16.*T3new(i))        !MZ
     &           *(14.*q_NOHot(i)*Nnqn0)                    !MZ
     &           +nune(i)*Cen*Uenew(i)

          enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*log(Nenew(np)*Tepnew(np)*Nenew(nx)*Tepnew(nx))
      D2br=.5*log(xn3(np)*T3pnew(np)*xn3(nx)*T3pnew(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimvd(Radn(nx+1),U3old,extra,nx)
      rbc=1.
      call lcpfct(U3old,U3new,Ipos1,Iposn,
     &              lbc,0.,0.,U3new(np),.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            U3new(i)=dexpnu*(U3new(i)-U3old(i)+D3(i)*deltat)
     &             +expnu*U3old(i)
          enddo

      if (any(ieee_is_nan(U3new))) then
        write(stderr,*) 'problem U3new dans la boucle 2'
        goto 246
      endif


C]]]


C[[[    Velocities corrections

          do 303 i=1,nx

         Uenew(i)=(Ci0*N2new(i)*U2new(i)
     &           +Cj0*N1new(i)*U1new(i)
     &           +Cn0*N3new(i)*U3new(i)
     &           +Cz0*(N4new(i)+N6new(i)+N5new(i))*Umnew(i)
     &            -JJ(i)/N_0-0.*Jes(i)/N_0)
     &               /Nenew(i)/Ce0

303        continue

      call velocity(Veliq,Ipos1,Iposnp,deltat)


C[[[    H+ heat flow equation resolution

          do i=1,nx

       C2a(i)=-2.2*q2new(i)
       D2a(i)=U2new(i)
       C2b(i)=-(11./18.*T2pnew(i)+8./9.*T2tnew(i))*N2new(i)
       D2b(i)=T2pnew(i)
       C2c(i)=-(17./9.*T2pnew(i)-8./9.*T2tnew(i))*N2new(i)
       D2c(i)=T2tnew(i)
       C2d(i)=4./9.*(T2pnew(i)-T2tnew(i))**2
       D2d(i)=N2new(i)

      D3(i)=alt_geo_1(i)*N2new(i)*(T2pnew(i)-T2tnew(i))
     &                           *(T2pnew(i)-4.*T2tnew(i))
     &            -N2new(i)*T2new(i)*U2new(i)
     &                     *(nuij(i)*18.5/17.
     &                      +nuik(i)*30.5/29.
     &                      +nuil(i)*34.5/33.
     &                      +nuim(i)*32.5/31.
     &                      +nuin(i)*16.5/15.
     &                      +nuie(i)*2.5
     &                      +nuiN2(i)*30.5/29.
     &                      +nuiO2(i)*34.5/33.
     &                      +nuiO(i)*18.5/17.
     &                      +nuiOHot(i)*18.5/17.                !MZ
     &                      +nuiH(i)*3.5/2.)
     &            + 2.5*N2new(i)*T2new(i)*U2new(i)*(
     &                nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+nuiH(i)    !MZ
     &                      +nuij(i)
     &                      +nuik(i)+nuil(i)+nuim(i)
     &                      +nuin(i)
     &                      +nuie(i))

            D3(i)=D3(i)+N2new(i)*T2new(i)*U1new(i)*Cji*
     &                           (nuij(i)*18.5/17.
     &                           -2.5*nuij(i))

            D3(i)=D3(i)+N2new(i)*T2new(i)*Umnew(i)*Czi*
     &                           (nuik(i)*30.5/29.
     &                           +nuil(i)*34.5/33.+nuim(i)*32.5/31.
     &                           -2.5*(nuik(i)+nuil(i)+nuim(i))     )

            D3(i)=D3(i)+N2new(i)*T2new(i)*U3new(i)*Cni*
     &                           (nuin(i)*16.5/15.
     &                           -2.5*nuin(i)         )

            D3(i)=D3(i)+N2new(i)*q1new(i)*Cji*xn1_1(i)*
     &            nuij(i)*(.0612+1.5*T2new(i)/(T1new(i)+16.*T2new(i)))

            D3(i)=D3(i)+N2new(i)*qenew(i)*Cei*xne_1(i)*
     &            nuie(i)*1.5*(T2new(i)/(Tenew(i)+5.46e-4*T2new(i))-1)

            D3(i)=D3(i)+N2new(i)*q3new(i)*Cni*xn3_1(i)*
     &            nuin(i)*(.068+1.5*T2new(i)/(T3new(i)+14.*T2new(i)))

            D3(i)=D3(i)+N2new(i)*(
     &                            nuiN2(i)*(1.07*q_Nn2(i)*Niqi0
     &                                     +1.052*T2new(i)*Un(i)/Ci0)
     &                           +nuiO2(i)*(1.084*q_No2(i)*Niqi0
     &                                     +1.045*T2new(i)*Un(i)/Ci0)
     &                           +nuiO(i) *(.98*q_No(i)*Niqi0
     &                                     +1.088*T2new(i)*Un(i)/Ci0)
     &                           +nuiOHot(i) *(.98*q_NOHot(i)*Niqi0        !MZ
     &                                     +1.088*T2new(i)*Un(i)/Ci0)        !MZ
     &                           +nuiH(i) *(-.075*q_Nh(i)*Niqi0
     &                                     +1.75*T2new(i)*Un(i)/Ci0 )
     &                           )
     &          -2.5*T2new(i)*N2new(i)*        (
     &                    (nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+nuiH(i))    !MZ
     &                   *(Un(i)/Ci0)
     &                   -.6*(nuiN2(i)/(Tn(i)/T_0+28.*T2new(i))
     &                       *(q_Nn2(i)*Niqi0)
     &                       +nuiO2(i)/(Tn(i)/T_0+32.*T2new(i))
     &                       *(q_No2(i)*Niqi0)
     &                       +nuiO (i)/(Tn(i)/T_0+16.*T2new(i))
     &                       *(q_No(i)*Niqi0)
     &                       +nuiOHot(i)/(TnOHot(i)/T_0+16.*T2new(i))        !MZ
     &                       *(q_NOHot(i)*Niqi0)                !MZ
     &                       +nuiH (i)/(Tn(i)/T_0+    T2new(i))
     &                       *(q_Nh(i)*Niqi0))  )

            D7(i)=
     &            (
     &              nuiN2(i)*.1795+nuiO2(i)*.1824
     &              +nuiO(i)*.1612+nuiOHot(i)*.1612-nuiH(i)*.725        !MZ
     &              -nuii(i)*.8+nuij(i)*.1612+nuik(i)*.1795
     &              +nuil(i)*.1824+nuim(i)*.1811+nuin(i)*.1547
     &              -nuie(i)*3.
     &            )
     &            -1.5*T2new(i)*(
     &                           nuiN2(i)*28./(Tn(i)/T_0+28.*T2new(i))
     &                          +nuiO2(i)*32./(Tn(i)/T_0+32.*T2new(i))
     &                          +nuiO (i)*16./(Tn(i)/T_0+16.*T2new(i))
     &                +nuiOHot(i)*16./(TnOHot(i)/T_0+16.*T2new(i))    !MZ
     &                          +nuiH (i)/(Tn(i)/T_0+T2new(i))
     &                          +nuij(i)*16./(T1new(i)+16.*T2new(i))
     &                          +nuik(i)*28./(Tmnew(i)+28.*T2new(i))
     &                          +nuil(i)*32./(Tmnew(i)+32.*T2new(i))
     &                          +nuim(i)*30./(Tmnew(i)+30.*T2new(i))
     &                          +nuin(i)*14./(T3new(i)+14.*T2new(i))
     &                          +nuie(i)*5.46e-4/
     &                               (Tenew(i)+5.46e-4*T2new(i))     )
     &            -4.2*U2new(i)*alt_geo_1(i)

          enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
          D2cl=(D2c(1)+D2c(2))/2.
          D2cr=ylimd(Radn(np),D2c,extra,nx)
          D2dl=(D2d(1)+D2d(2))/2.
          D2dr=ylimd(Radn(np),D2d,extra,nx)
    
      D2ar=.5*(U2new(np)+U2new(nx))
      D2br=.5*(T2pnew(np)+T2pnew(nx))
      D2cr=.5*(T2tnew(np)+T2tnew(nx))
      D2dr=sqrt(N2new(np)*N2new(nx))
    
      call sources(Ipos1,Iposn,deltat_2,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat_2,2,C2b,D2b,D2bl,D2br)
      call sources(Ipos1,Iposn,deltat_2,2,C2c,D2c,D2cl,D2cr)
      call sources(Ipos1,Iposn,deltat_2,2,C2d,D2d,D2dl,D2dr)
c      call sources(Ipos1,Iposn,deltat_2,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat_2,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),q2old,extra,nx)
      rbc=1.
      call lcpfct(q2old,q2new,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            q2new(i)=dexpnu*(q2new(i)-q2old(i)+D3(i)*deltat)
     &             +expnu*q2old(i)
          enddo
          q2new(nx)=max(0.,q2new(nx))
          q2new(np)=q2new(nx)

      if (any(ieee_is_nan(q2new))) then
        write(stderr,*) 'problem q2new dans la boucle 2'
        goto 246
      endif


C]]]



      call velocity(Veljq,Ipos1,Iposnp,deltat)

C[[[    O+ heat flow equation resolution

         do i=1,nx

       C2a(i)=-2.2*Cji*q1new(i)
       D2a(i)=U1new(i)
       C2b(i)=-Cji*(11./18.*T1pnew(i)+8./9.*T1tnew(i))*N1new(i)
       D2b(i)=T1pnew(i)
       C2c(i)=-Cji*(17./9.*T1pnew(i)-8./9.*T1tnew(i))*N1new(i)
       D2c(i)=T1tnew(i)
       C2d(i)=4./9.*Cji*(T1pnew(i)-T1tnew(i))**2
       D2d(i)=N1new(i)

      D3(i)=Cji*alt_geo_1(i)*N1new(i)*(T1pnew(i)-T1tnew(i))
     &                           *(T1pnew(i)-4.*T1tnew(i))
     &            +N1new(i)*T1new(i)*U2new(i)*Cij*
     &                   nuji(i)*(2.412-2.5)

          D3(i)=D3(i)-N1new(i)*T1new(i)*U1new(i)*                    (
     &                         nujN2(i)*1.545+nujO2(i)*1.5
     &                        +nujO(i)*1.75  +nuji(i)*2.412
     &                        +nujOHot(i)*1.75                    !MZ
     &                        +nujn(i)*1.8   +nujk(i)*1.545
     &                        +nujl(i)*1.5   +nujm(i)*1.522
     &                        +nuje(i)*2.5
     &                - 2.5*(nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i)    !MZ
     &                              +nuji(i)+nujn(i)
     &                              +nujk(i)+nujl(i)+nujm(i)+nuje(i)))

          D3(i)=D3(i)+N1new(i)*T1new(i)*Umnew(i)*Czj*                (
     &                         nujk(i)*1.545+nujl(i)*1.5+nujm(i)*1.522
     &                        - 2.5*(nujk(i)+nujl(i)+nujm(i))        )

          D3(i)=D3(i)-N1new(i)*T1new(i)*U3new(i)*Cnj*.7*nujn(i)

          D3(i)=D3(i)-N1new(i)*q2new(i)*Cij*xn2_1(i)*nuji(i)*
     &                   (1.262- 24.*T1new(i)/(16.*T2new(i)+T1new(i)))

           D3(i)=D3(i)-N1new(i)*qenew(i)*Cej*xne_1(i)*nuje(i)*
     &             (1.5-24.*T1new(i)/(16.*Tenew(i)+5.46e-4*T1new(i)))

           D3(i)=D3(i)-N1new(i)*q3new(i)*Cnj*xn3_1(i)*nujn(i)*
     &             (.128-24.*T1new(i)/(16.*T3new(i)+14.*T1new(i)))

           D3(i)=D3(i)+N1new(i)*(
     &                           nujN2(i)*(.139*q_Nn2(i)*Njqj0
     &                                    +1.545*T1new(i)*Un(i)/Cj0)
     &                          +nujO2(i)*(.2*q_No2(i)*Njqj0
     &                                    +1.5*T1new(i)*Un(i)/Cj0)
     &                          +nujO(i)* (-.075*q_No(i)*Njqj0
     &                                    +1.75*T1new(i)*Un(i)/Cj0)
     &                          +nujOHot(i)* (-.075*q_NOHot(i)*Njqj0        !MZ
     &                                    +1.75*T1new(i)*Un(i)/Cj0)        !MZ
     &                          )
     &           -2.5*T1new(i)*N1new(i)*(
     &            (nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i))*(Un(i)/Cj0)    !MZ
     &               -9.6*Njqj0*
     &                  (q_Nn2(i)*nujN2(i)/(16.*Tn(i)/T_0+28.*T1new(i))
     &                  +q_No2(i)*nujO2(i)/(16.*Tn(i)/T_0+32.*T1new(i))
     &                  +q_No(i)*nujO (i)/(16.*Tn(i)/T_0+16.*T1new(i))
     &          +q_NOHot(i)*nujOHot(i)/(16.*TnOHot(i)/T_0+16.*T1new(i))    !MZ
     &                  )               )

         D7(i)=-(nujN2(i)*.34+nujO2(i)*.27+nujO(i)*.725+nujOHot(i)*.725    !MZ
     &           +nuji(i)*2.66+nujk(i)*.34+nujl(i)*.27
     &           +nujm(i)*.3+nujn(i)*.835+nujj(i)*.8+nuje(i)*3.)
     &          -1.5*T1new(i)*(
     &                        nujN2(i)*28./(16.*Tn(i)/T_0+28.*T1new(i))
     &                       +nujO2(i)*32./(16.*Tn(i)/T_0+32.*T1new(i))
     &                       +nujO (i)*16./(16.*Tn(i)/T_0+16.*T1new(i))
     &                 +nujOHot(i)*16./(16.*TnOHot(i)/T_0+16.*T1new(i))    !MZ
     &                       +nuji(i)/(16.*T2new(i)+T1new(i))
     &                       +nujk(i)*28./(16.*Tmnew(i)+28.*T1new(i))
     &                       +nujl(i)*32./(16.*Tmnew(i)+32.*T1new(i))
     &                       +nujm(i)*30./(16.*Tmnew(i)+30.*T1new(i))
     &                       +nujn(i)*14./(16.*T3new(i)+14.*T1new(i))
     &                       +nuje(i)*5.46e-4
     &                               /(16.*Tenew(i)+5.46e-4*T1new(i)))
     &          -4.2*Cji*U1new(i)*alt_geo_1(i)

         enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
          D2cl=(D2c(1)+D2c(2))/2.
          D2cr=ylimd(Radn(np),D2c,extra,nx)
          D2dl=(D2d(1)+D2d(2))/2.
          D2dr=ylimd(Radn(np),D2d,extra,nx)
    
      D2ar=.5*(U1new(np)+U1new(nx))
      D2br=.5*(T1pnew(np)+T1pnew(nx))
      D2cr=.5*(T1tnew(np)+T1tnew(nx))
      D2dr=sqrt(N1new(np)*N1new(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
      call sources(Ipos1,Iposn,deltat,2,C2c,D2c,D2cl,D2cr)
      call sources(Ipos1,Iposn,deltat,2,C2d,D2d,D2dl,D2dr)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),q1old,extra,nx)
      rbc=1.
      call lcpfct(q1old,q1new,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            q1new(i)=dexpnu*(q1new(i)-q1old(i)+D3(i)*deltat)
     &             +expnu*q1old(i)
          enddo
          q1new(nx)=max(0.,q1new(nx))
          q1new(np)=q1new(nx)

      if (any(ieee_is_nan(q1new))) then
        write(stderr,*) 'problem q1new dans la boucle 2'
        goto 246
      endif



C]]]

      call velocity(Velnq,Ipos1,Iposnp,deltat)


C[[[    N+ heat flow equation resolution

         do i=1,nx

       C2a(i)=-2.2*Cni*q3new(i)
       D2a(i)=U3new(i)
       C2b(i)=-Cni*(11./18.*T3pnew(i)+8./9.*T3tnew(i))*N3new(i)
       D2b(i)=T3pnew(i)
       C2c(i)=-Cni*(17./9.*T3pnew(i)-8./9.*T3tnew(i))*N3new(i)
       D2c(i)=T3tnew(i)
       C2d(i)=4./9.*Cni*(T3pnew(i)-T3tnew(i))**2
       D2d(i)=N3new(i)

      D3(i)=Cni*alt_geo_1(i)*N3new(i)*(T3pnew(i)-T3tnew(i))
     &                           *(T3pnew(i)-4.*T3tnew(i))
     &            -.1*N3new(i)*U2new(i)*Cin*T3new(i)*nuni(i)
     &          -.8*N3new(i)*U1new(i)*Cjn*T3new(i)*nunj(i)
     &          -   N3new(i)*Umnew(i)*Czn*T3new(i)
     &                        *(nunk(i)+1.04*nunl(i)+1.02*nunm(i))

          D3(i)=D3(i)+N3new(i)*T3new(i)*U3new(i)*               (
     &                         nunN2(i)+nunO2(i)*1.04+nunO(i)*.8
     &                +nunOHot(i)*.8                    !MZ
     &                        +nuni(i)*.1+nunj(i)*.8+nunk(i)
     &                        +nunl(i)*1.04+nunm(i)*1.02        )

          D3(i)=D3(i)+N3new(i)*q2new(i)*Cin*xn2_1(i)*nuni(i)*
     &               (21.*T3new(i)/(14.*T2new(i)+T3new(i))-1.232)

          D3(i)=D3(i)+N3new(i)*q1new(i)*Cjn*xn1_1(i)*nunj(i)*
     &               (21.*T3new(i)/(14.*T1new(i)+16.*T3new(i))-.028)

          D3(i)=D3(i)+N3new(i)*qenew(i)*Cen*xne_1(i)*nune(i)*
     &               (21.*T3new(i)/(14.*Tenew(i)+5.46e-4*T3new(i))-1.5)

          D3(i)=D3(i)+N3new(i)*  (
     &            nunN2(i)* (.1*q_Nn2(i)*Nnqn0
     &                      +1.5*T3new(i)*Un(i)/Cn0)
     &            +nunO2(i)*(.115*q_No2(i)*Nnqn0
     &                      +1.456*T3new(i)*Un(i)/Cn0)
     &            +nunO(i)* (-.028*q_No(i)*Nnqn0
     &                      +1.7*T3new(i)*Un(i)/Cn0)
     &            +nunOHot(i)* (-.028*q_NOHot(i)*Nnqn0                !MZ
     &                      +1.7*T3new(i)*Un(i)/Cn0)                !MZ
     &                           )
     &          -2.5*T3new(i)*N3new(i)*(
     &            (nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i))*(Un(i)/Cn0)    !MZ
     &           -8.4*Nnqn0*(
     &                  nunN2(i)*q_Nn2(i)/(14.*Tn(i)/T_0+28.*T3new(i))
     &                 +nunO2(i)*q_No2(i)/(14.*Tn(i)/T_0+32.*T3new(i))
     &                 +nunO (i)*q_No(i)/(14.*Tn(i)/T_0 +16.*T3new(i))
     &        +nunOHot(i)*q_NOHot(i)/(14.*TnOHot(i)/T_0 +16.*T3new(i)))    !MZ
     &                                 )

          D7(i)=-( nunN2(i)*.27+nunO2(i)*.2+nunO(i)*.62
     &            +nunOHot(i)*.62+nuni(i)*2.62                !MZ
     &            +nunk(i)*.27+nunl(i)*.2 +nunm(i)*.23
     &            +nunj(i)*.62+nunn(i)*.8 +nune(i)*3.)
     &            -1.5*T3new(i)*(
     &                 nunN2(i)*28./(14.*Tn(i)/T_0+28.*T3new(i))
     &                +nunO2(i)*32./(14.*Tn(i)/T_0+32.*T3new(i))
     &                +nunO (i)*16./(14.*Tn(i)/T_0+16.*T3new(i))
     &            +nunOHot(i)*16./(14.*TnOHot(i)/T_0+16.*T3new(i))    !MZ
     &                +nuni(i)/(14.*T2new(i)+T3new(i))
     &                +nunk(i)*28./(14.*Tmnew(i)+28.*T3new(i))
     &                +nunl(i)*32./(14.*Tmnew(i)+32.*T3new(i))
     &                +nunm(i)*30./(14.*Tmnew(i)+30.*T3new(i))
     &                +nunj(i)*16./(14.*T1new(i)+16.*T3new(i))
     &                +nune(i)*5.46e-4/(14.*Tenew(i)+5.46e-4*T3new(i)))
     &          -4.2*Cni*U3new(i)*alt_geo_1(i)

         enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
          D2cl=(D2c(1)+D2c(2))/2.
          D2cr=ylimd(Radn(np),D2c,extra,nx)
          D2dl=(D2d(1)+D2d(2))/2.
          D2dr=ylimd(Radn(np),D2d,extra,nx)
    
      D2ar=.5*(U3new(np)+U3new(nx))
      D2br=.5*(T3pnew(np)+T3pnew(nx))
      D2cr=.5*(T3tnew(np)+T3tnew(nx))
      D2dr=sqrt(N3new(np)*N3new(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
      call sources(Ipos1,Iposn,deltat,2,C2c,D2c,D2cl,D2cr)
      call sources(Ipos1,Iposn,deltat,2,C2d,D2d,D2dl,D2dr)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),q3old,extra,nx)
      rbc=1.
      call lcpfct(q3old,q3new,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            q3new(i)=dexpnu*(q3new(i)-q3old(i)+D3(i)*deltat)
     &             +expnu*q3old(i)
          enddo

C[[[   Boundaries conditions

          q3new(nx)=max(0.,q3new(nx))
          q3new(np)=q3new(nx)

      if (any(ieee_is_nan(q3new))) then
        write(stderr,*) 'problem when calculating q3new in loop 2'
        goto 246
      endif


C]]]



C[[[

C[[[    Electron energy and heat flow equation resolution (2)

    

      do i=1,nx

        Ter=T_0*Tenew(i)
        Terh=sqrt(Ter)

        nueN2(i)=2.33e-11*Nn2(i)
     &               *exp(-1.21e-4*Ter)*Ter            *t0
        nueO2(i)=1.82e-10*No2(i)
     &             *(1.+3.6e-2*Terh)*Terh            *t0
        nueO (i)=8.90e-11*No (i)
     &             *(1.+5.7e-4*Ter)*Terh            *t0
        nueOHot(i)=8.90e-11*NOHot(i)                    !MZ
     &             *(1.+5.7e-4*Ter)*Terh            *t0        !MZ
        nueH (i)=4.5e-9*Nh(i)
     &             *exp(-1.35e-4*Ter)*Terh            *t0
        nuei (i)=54.5*N2new(i)*Te_15(i)
        nuej (i)=54.5*N1new(i)*Te_15(i)
        nuee (i)=54.5*Nenew(i)*Te_15(i)/sqrt(2.)
            nuek (i)=54.5*N4new(i)*Te_15(i)
            nuel (i)=54.5*N6new(i)*Te_15(i)
            nuem (i)=54.5*N5new(i)*Te_15(i)
            nuen (i)=54.5*N3new(i)*Te_15(i)

            C2ea(i)=Cei*Tenew(i)/3.
            D2ea(i)=Uenew(i)
            C2eb(i)=-.667*Cei*xne_1(i)
            D2eb(i)=qenew(i)
        
        !compression term
            C2qa(i)=-11./5.*Cei*qenew(i)
            D2qa(i)=Uenew(i)

        !Temperature gradient (thermal conduction)
            C2qb(i)=-2.5*Cei*Nenew(i)*Tenew(i)
            D2qb(i)=Tenew(i)


          qen=Tenew(i)*(Un(i)/Ce0-Uenew(i))

          D3q(i)=-1.5*Nenew(i)*Tenew(i)*(nuei(i)*Cie*U2new(i)
     &                                 +nuej(i)*Cje*U1new(i)
     &                                 +(nuek(i)+nuel(i)+nuem(i))
     &                                 *Cze*Umnew(i)
     &                                 +nuen(i)*Cne*U3new(i))

          D3q(i)=D3q(i)+Nenew(i)*q2new(i)*Cie*xn2_1(i)*nuei(i)*
     &                         (6.55e-4+8.19e-4*Tenew(i)
     &                         /(5.46e-4*T2new(i)+Tenew(i)))

          D3q(i)=D3q(i)+Nenew(i)*q1new(i)*Cje*xn1_1(i)*nuej(i)*
     &                         (4.09e-5+8.19e-4*Tenew(i)
     &                         /(5.46e-4*T1new(i)+16.*Tenew(i)))

          D3q(i)=D3q(i)+Nenew(i)*q3new(i)*Cne*xn3_1(i)*nuen(i)*
     &                         (4.68e-5+8.19e-4*Tenew(i)
     &                         /(5.46e-4*T3new(i)+14.*Tenew(i)))

          D3q(i)=D3q(i)+Nenew(i)* (
     &            nueN2(i)*(1.2*q_Nn2(i)*Neqe0+qen)
     &           +nueO2(i)*(1.2*q_No2(i)*Neqe0+qen)
     &           +nueO(i)* (1.2*q_No(i)*Neqe0 +qen)
     &           +nueOHot(i)* (1.2*q_NOHot(i)*Neqe0 +qen)            !MZ
     &           +nueH(i)* (1.2*q_Nh(i)*Neqe0 +qen)
     &          +1.5*(nuei(i)+nuej(i)+nuek(i)+nuel(i)+nuem(i)+nuen(i))
     &              *Tenew(i)*Uenew(i)
     &                          )
     &          -2.5*Nenew(i)*Tenew(i)*(
     &               (nueN2(i)+nueO2(i)+nueO(i)+nueOHot(i)+nueH(i))        !MZ
     &              *(Un(i)/Ce0-Uenew(i))
     &              -.6*nueN2(i)/(5.46e-4*Tn(i)/T_0+28.*Tenew(i))
     &                 *(5.46e-4*q_Nn2(i)*Neqe0)
     &              -.6*nueO2(i)/(5.46e-4*Tn(i)/T_0+32.*Tenew(i))
     &                 *(5.46e-4*q_No2(i)*Neqe0)
     &              -.6*nueO (i)/(5.46e-4*Tn(i)/T_0+16.*Tenew(i))
     &                 *(5.46e-4*q_No(i)*Neqe0)
     &              -.6*nueOHot(i)/(5.46e-4*TnOHot(i)/T_0+16.*Tenew(i))        !MZ
     &                 *(5.46e-4*q_NOHot(i)*Neqe0)                !MZ
     &              - .6*nueH (i)/(5.46e-4*Tn(i)/T_0+Tenew(i))
     &                 *(5.46e-4*q_Nh(i)*Neqe0)
     &                                 )

          D7q(i)=.2*(nueN2(i)+nueO2(i)+nueO(i)+nueOHot(i)+nueH(i)        !MZ
     &             +nuei(i) +nuej(i) +nuek(i)+nuel(i)
     &             +nuem(i) +nuen(i) -nuee(i)*4       )
     &          -1.5*Tenew(i)*(
     &                 nueN2(i)*28./(5.46e-4*Tn(i)/T_0+28.*Tenew(i))
     &                +nueO2(i)*32./(5.46e-4*Tn(i)/T_0+32.*Tenew(i))
     &                +nueO (i)*16./(5.46e-4*Tn(i)/T_0+16.*Tenew(i))
     &            +nueOHot(i)*16./(5.46e-4*TnOHot(i)/T_0+16.*Tenew(i))    !MZ
     &                +nueH (i)/(5.46e-4*Tn(i)/T_0+Tenew(i))
     &                +nuei(i)/(5.46e-4*T2new(i)+Tenew(i))
     &                +nuej(i)*16./(5.46e-4*T1new(i)+16.*Tenew(i))
     &                +nuek(i)*28./(5.46e-4*Tmnew(i)+28.*Tenew(i))
     &                +nuel(i)*32./(5.46e-4*Tmnew(i)+32.*Tenew(i))
     &                +nuem(i)*30./(5.46e-4*Tmnew(i)+30.*Tenew(i))
     &                +nuen(i)*14./(5.46e-4*T3new(i)+14.*Tenew(i)))
     &          -4.2*Cei*Uenew(i)*alt_geo_1(i)                                  !geometric term

      if (abs(Ter-Tn(i)).lt.1.e-2) then
         dTen=1.e-2
         if (Ter.ge.Tn(i)) dTen=-dTen
         Ter=Tn(i)-dTen
      else
         dTen=Tn(i)-Ter
      endif
    
      dTen_1=(Tn(i)-Ter)/Ter/Tn(i)
    
          Lenrot=(4.6e-26*Nn2(i)+1.1e-25*No2(i))
     &            *dTen/sqrt(Ter)

          if (Ter.lt.10000.) then
            f=exp(2.2e-3*(Ter-1800.))
            f=1.06e4+7.51e3*(f-1.)/(f+1.)
          else
            f=18110.
          endif

          gg=3300+(1.233-2.056e-4*(Ter-4000.))*(Ter-1000.)


            LeN2vib=4.78e-24*Nn2(i)
     &              *exp(f*(Ter-2000.)/2000./Ter)
     &           *(exp(gg*dTen_1)-1.)

          h=3300.-839.*sin(1.91e-4*(Ter-2700.))

          LeO2vib=8.3e-25*No2(i)
     &            *exp(h*(Ter-700.)/700./Ter)
     &         *(exp(2770.*dTen_1)-1.)

      dd=2.4e4+(.3-1.947e-5*(Ter-4000))*(Ter-1500.)

          LeOexc=2.5e-24*No(i)
     &      *exp(dd*(Ter-3000.)/Ter/3000.)
     &      *(exp(22713.*dTen_1)-1.)


          LeOfin=0.
          do j=1,3
            if (j.le.2) then
              Ex=exp(-Eq(j)/Ter)
              Dx(j)=exp(-Eq(j)/Tn(i))
            else
              Ex=exp(-Eq(3)/Ter-Eq(1)/Tn(i))
              Dx(j)=exp(-Eq(2)/Tn(i))
            endif
            F=epsq(j)*(Dx(j)-Ex)-5.91e-9*dTen
     &          *((1.+Bq(j))*Dx(j)+Ex*(Eq(j)/Ter+1.+Bq(j)))
            LeOfin=LeOfin+Aq(j)*Cq(j)*F*(Ter)**(Bq(j)-.5)
          enddo
          Z=5.+3.*Dx(1)+Dx(2)
          LeOfin=1.38e-17*No(i)/Z*LeOfin

          LeOfin=Sq(3)*Ter**.6*exp(Eq(3)/Tn(i))
     &              *(exp(Eq(3)*dTen_1)-1.)
    
          do j=1,2
        LeOfin=LeOfin+Sq(j)*(exp(Eq(j)*dTen_1)-1.)
      enddo
    
          Z=5.+3.*Dx(1)+Dx(2)
          LeOfin=No(i)/Z*LeOfin

c    modif cargese
      Len(i)=(LeOfin+Lenrot+LeN2vib+LeO2vib+LeOexc)*N_0*t0/P_0
      Len(i)=Len(i)/dTen

      D7e(i)=-.667*Len(i)*T_0
c      Len(i)=Heat(i)*N_0*t0/P_0+Len(i)*Tn(i)
      D3e(i)=.667*(Heat(i)/xne(i)*t0/P_0+Len(i)*Tn(i))

c    fin de modif


          D3e(i)= D3e(i)+T2new(i)*1.09e-3*nuei(i)
     &          +T1new(i)*6.8e-5*nuej(i)
     &          +Tmnew(i)*1.09e-3*(nuek(i)/28.+nuel(i)/32.+nuem(i)/30.)
     &          +T3new(i)*7.8e-5*nuen(i)


        nu_omega=(nueN2(i)+nueO2(i)+nueO(i)+
     &        nueOHot(i)+nueH(i))/omega(i)*ome                !MZ
         coef_conv=1./(1.+nu_omega**2)
         terme_Joule=.5*coef_conv*Vm_2(i)/Ce0**2
    
          D3e(i)=D3e(i)+1.09e-3*(nueN2(i)/28.+nueO2(i)/32.+nueO(i)/16.
     &                        +nueH(i))*Tn(i)/T_0
     &            +1.09e-3*nueOHot(i)/16.*TnOHot(i)/T_0            !MZ
     &               + .667*(nueN2(i)*((Uenew(i)-Un(i)/Ce0)**2
     &               + terme_Joule)
     &               + nueO2(i)*((Uenew(i)-Un(i)/Ce0)**2
     &               + terme_Joule)
     &               + nueO(i)*((Uenew(i)-Un(i)/Ce0)**2
     &               + terme_Joule)
     &               + nueOHot(i)*((Uenew(i)-Un(i)/Ce0)**2            !MZ
     &               + terme_Joule)
     &               + nueH(i)*((Uenew(i)-Un(i)/Ce0)**2
     &               + terme_Joule)
     &               + nuei(i)*(Cie*U2new(i)-Uenew(i))**2
     &               + nuej(i)*(Cje*U1new(i)-Uenew(i))**2
     &               + nuek(i)*(Cze*Umnew(i)-Uenew(i))**2
     &               + nuel(i)*(Cze*Umnew(i)-Uenew(i))**2
     &               + nuem(i)*(Cze*Umnew(i)-Uenew(i))**2
     &               + nuen(i)*(Cne*U3new(i)-Uenew(i))**2
     &                        )
     &               - 2.*(Cei*qenew(i)*alt_geo_1(i)*xne_1(i))

      D7e(i)=D7e(i)
     &          - 1.09e-3*(nueN2(i)/28.+ nueO2(i)/32. + nueO(i)/16.
     &          + nueOHot(i)/16.+nueH(i) + nuei(i) + nuej(i)/16.        !MZ
     &          + nuek(i)/28. + nuel(i)/32. + nuem(i)/30.+nuen(i)/14.)
c     &          - 2.*(Cei*Uenew(i)*alt_geo_1(i))

       enddo

    
c    premiere demi-boucle

          D2qal=(D2qa(1)+D2qa(2))/2.
          D2qbl=(D2qb(1)+D2qb(2))/2.

          D2qar=.5*(Uenew(np)+Uenew(nx))
          D2qbr=.5*(Tenew(np)+Tenew(nx))
    
          call velocity(Veleq,Ipos1,Iposnp,deltat_4)
          call sources(Ipos1,Iposn,deltat_4,2,C2qa,D2qa,D2qal,D2qar)
          call sources(Ipos1,Iposn,deltat_4,2,C2qb,D2qb,D2qbl,D2qbr)
c      call sources(Ipos1,Iposn,deltat_4,3,zero,D3q,0.,0.)
        call sources(Ipos1,Iposn,deltat_4,7,zero,D7q,0.,0.)
          lbc=1.
          rbc=1.
        call lcpfct(qeold,qenew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          qenew(nx)=qetop

     
          do i=1,nx
            dexpnu=D7q(i)*deltat_4
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            qenew(i)=dexpnu*(qenew(i)-qeold(i)+deltat_4*D3q(i))
     &            +expnu*qeold(i)
          enddo
          qenew(nx)=qetop
          qenew(np)=qenew(nx)

      if (any(ieee_is_nan(qenew))) then
        call cpu_time(tic)
        write(stderr,*) tic,'L6910: problem calc  qenew in loop 1'
        goto 246
      endif

      D2eal=(D2ea(1)+D2ea(2))/2.
      D2ebl=(D2eb(1)+D2eb(2))/2.

      D2ear=.5*(Uenew(np)+Uenew(nx))
      D2ebr=.5*(qenew(np)+qenew(nx))
    
      call sources(Ipos1,Iposn,deltat_4,2,C2ea,D2ea,D2eal,D2ear)
      call sources(Ipos1,Iposn,deltat_4,2,C2eb,D2eb,D2ebl,D2ebr)
c      call sources(Ipos1,Iposn,deltat_4,3,zero,D3e,0.,0.)
      call sources(Ipos1,Iposn,deltat_4,7,zero,D7e,0.,0.)
          lbc=1.
      rbc=1.
       call lcpfct(Teold,Tenew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)

          do i=1,nx
            dexpnu=D7e(i)*deltat_4
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Tenew(i)=dexpnu*(Tenew(i)-Teold(i)+deltat_4*D3e(i))
     &             +expnu*Teold(i)
          enddo
          tenew(np)=tenew(nx)
          Tenew(np)=2.*Tenew(nx)-Tenew(nx-1)

      if (any(ieee_is_nan(Tenew))) then
        call cpu_time(tic)
        write(stderr,*) tic,'L6966: problem when calc Tenew  loop 1'
        goto 246
      endif


      call stabenerg(nx,xne,Teold,Tenew,qeold,qenew,
     &            D7e,D7q,Cei,deltat_4)


      if (any(ieee_is_nan(qenew))) then
        write(stderr,*) 'problem in stabenerg with qe in loop 1'
        goto 246
      endif
      if (any(ieee_is_nan(Tenew))) then
        write(stderr,*) 'problem in stabenerg with Te in loop 1'
        goto 246
      endif

c deuxieme demi-boucle

          do i=1,np
            Tenew(i)=max(Tenew(i),T_min)
        Tepnew(i)=Tenew(i)
        Tetnew(i)=Tenew(i)
            qeold(i)=qenew(i)
            Teold(i)=Tenew(i)
          enddo

          do i=1,nx
            C2ea(i)=Cei*Tenew(i)/3.
            D2eb(i)=qenew(i)

            C2qa(i)=-11./5.*Cei*qenew(i)

            C2qb(i)=-2.5*Cei*Nenew(i)*Tenew(i)
            D2qb(i)=Tenew(i)
      enddo

          D2qbl=(D2qb(1)+D2qb(2))/2.
          D2qbr=.5*(Tenew(np)+Tenew(nx))
    
      call sources(Ipos1,Iposn,deltat_4,2,C2qa,D2qa,D2qal,D2qar)
      call sources(Ipos1,Iposn,deltat_4,2,C2qb,D2qb,D2qbl,D2qbr)
c      call sources(Ipos1,Iposn,deltat_4,3,zero,D3q,0.,0.)
      call sources(Ipos1,Iposn,deltat_4,7,zero,D7q,0.,0.)
          lbc=1.
          rbc=1.
      call lcpfct(qeold,qenew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          qenew(nx)=qetop

     
          do i=1,nx
            dexpnu=D7q(i)*deltat_4
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            qenew(i)=dexpnu*(qenew(i)-qeold(i)+deltat_4*D3q(i))
     &            +expnu*qeold(i)
          enddo
          qenew(nx)=qetop
          qenew(np)=qenew(nx)

      if (any(ieee_is_nan(qenew))) then
        call cpu_time(tic)
        write(stderr,*) tic,'L7014: problem when calc qenew in loop 1'
        goto 246
      endif



          D2ebl=(D2eb(1)+D2eb(2))/2.
          D2ebr=.5*(qenew(np)+qenew(nx))
    
      call sources(Ipos1,Iposn,deltat_4,2,C2ea,D2ea,D2eal,D2ear)
      call sources(Ipos1,Iposn,deltat_4,2,C2eb,D2eb,D2ebl,D2ebr)
c      call sources(Ipos1,Iposn,deltat_4,3,zero,D3e,0.,0.)
      call sources(Ipos1,Iposn,deltat_4,7,zero,D7e,0.,0.)
          lbc=1.
      rbc=1.
       call lcpfct(Teold,Tenew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)

          do i=1,nx
            dexpnu=D7e(i)*deltat_4
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Tenew(i)=dexpnu*(Tenew(i)-Teold(i)+deltat_4*D3e(i))
     &             +expnu*Teold(i)
          enddo
          tenew(np)=tenew(nx)
          Tenew(np)=2.*Tenew(nx)-Tenew(nx-1)

      if (any(ieee_is_nan(Tenew))) then
        call cpu_time(tic)
        write(stderr,*) tic,'L7049: problem when calc Tenew in loop 1'
        goto 246
      endif


      call stabenerg(nx,xne,Teold,Tenew,qeold,qenew,
     &            D7e,D7q,Cei,deltat_4)


      if (any(ieee_is_nan(qenew))) then
        write(stderr,*) 'problem stabenerg avec qe dans la boucle 1'
        goto 246
      endif
      if (any(ieee_is_nan(Tenew))) then
        write(stderr,*) 'problem stabenerg avec Te dans la boucle 1'
        goto 246
      endif
              


          do i=1,np
            Tenew(i)=max(Tenew(i),T_min)
        Tepnew(i)=Tenew(i)
        Tetnew(i)=Tenew(i)
            qeold(i)=qenew(i)
            Teold(i)=Tenew(i)
          enddo
    

C ]]]

C [[[

      call velocity(Velie,Ipos1,Iposnp,deltat)

C[[[    H+ energy equation resolution


c    temperature parallele

      do i=1,nx

      C2a(i)=-T2pnew(i)
      D2a(i)=U2new(i)
      C2b(i)=-1.2*xn2_1(i)
      D2b(i)=q2new(i)

       Q_cross2=max(.3,min(.45,.45-.15*(T2new(i)-2.5)))
        Q_cross2=.3

          D7(i)=
     &      -(.9862*nuiN2(i)+.9818*nuiO2(i)
     &        +(.1176+.9412*Q_cross2)*nuiO(i)+(1.+.5*Q_cross2)*nuiH(i)
     &        +(.1176+.9412*Q_cross2)*nuiOHot(i)                    !MZ
     &       +nuij(i)*.8706+nuie(i)*1.9993
     &       +nuik(i)*.8414+nuil(i)*.8364+nuim(i)*.8387+nuin(i)*.88)
     &        -.8*nuii(i)

          D3(i)=T2tnew(i)*(.7529*nuij(i)+.7467*nuin(i)+.7724*nuik(i)
     &                    +.7758*nuil(i) +.7742*nuim(i)
     &              +4.37e-4*nuie(i)
     &                    +.9172*nuiN2(i)+.9212*nuiO2(i)
     &                    +.9412*Q_cross2*nuiO(i)+.5*Q_cross2*nuiH(i)
     &              +.9412*Q_cross2*nuiOHot(i))                !MZ
     &           +T1pnew(i)*.0706*nuij(i)
     &         +Tmpnew(i)*(.0414*nuik(i)+.0364*nuil(i) +.0387*nuim(i))
     &         +Tepnew(i)*  1.1993*nuie(i)
     &         +T3pnew(i)* .0800*nuin(i)
     &         +T1tnew(i)*.0471*nuij(i)
     &         +Tmtnew(i)*(.0276*nuik(i)+.0242*nuil(i) +.0258*nuim(i))
     &         +Tetnew(i)*  .7996*nuie(i)
     &         +T3tnew(i)* .0533*nuin(i)
     &           +.8*nuii(i)*T2tnew(i)

        nu_omega=(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+
     &               nuiH(i))/omega(i)*omi    !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Ci0**2
        dUn2=(U2new(i)-Un(i)/Ci0)**2

          D3(i)=D3(i)
     &           +Tn(i)/T_0*(.069*nuiN2(i)+.0606*nuiO2(i)
     &                    +nuiH(i)+.1176*nuiO(i))
     &        +TnOHot(i)/T_0*.1176*nuiOHot(i)                    !MZ
     &         +nuiN2(i)*(1.0138*dUn2+.9172*terme_joule)
     &         +nuiO2(i)*(1.0182*dUn2+.9212*terme_joule)
     &         +nuiO(i)*.9412*(2.*dUn2+Q_cross2*(terme_joule-dUn2))
     &         +nuiOHot(i)*.9412*(2.*dUn2+Q_cross2*(terme_joule-dUn2))        !MZ
     &         +nuiH(i)*.5*(2.*dUn2+Q_cross2*(terme_joule-dUn2))
     &         +.667*(
     &         + nuij(i)*.94*(Cji*U1new(i)-U2new(i))**2
     &         + (nuik(i)+ nuil(i)+ nuim(i))
     &           *.97*(Czi*Umnew(i)-U2new(i))*(Czi*Umnew(i)-U2new(i))
     &         + nuin(i)*.93*(Cni*U3new(i)-U2new(i))**2
     &         + nuie(i)*5.46e-4*(Cei*Uenew(i)-U2new(i))**2)
     &         - 1.2*(q2new(i)*alt_geo_1(i)*xn2_1(i))

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*(U2new(np)+U2new(nx))
      D2br=.5*(q2new(np)+q2new(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T2pold,extra,nx)
      rbc=1.
       call lcpfct(T2pold,T2pnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T2pnew(i)=dexpnu*(T2pnew(i)-T2pold(i)+D3(i)*deltat)
     &             +expnu*T2pold(i)
          enddo
          T2pnew(np)=T2pnew(nx)

          do  i=1,np
            T2pnew(i)=max(T2pnew(i),T_min)
      enddo

      if (any(ieee_is_nan(T2pnew))) then
        write(stderr,*) 'problem T2pnew dans la boucle 2'
        goto 246
      endif



C]]]

c    temperature perpendiculaire

      do i=1,nx

      C2a(i)=T2tnew(i)
      D2a(i)=U2new(i)
      C2b(i)=-.4*xn2_1(i)
      D2b(i)=q2new(i)

        Q_cross2=max(.3,min(.45,.45-.15*(T2new(i)-2.5)))/2.
        Q_cross2=.3

          D7(i)=
     &      -(.5276*nuiN2(i)+.5212*nuiO2(i)
     &        +(.1176+.9412*Q_cross2)*nuiO(i)+(1.+.5*Q_cross2)*nuiH(i)
     &        +(.1176+.9412*Q_cross2)*nuiOHot(i)                    !MZ
     &       +nuij(i)*.4941+nuie(i)*1.9991
     &       +nuik(i)*.4552+nuil(i)*.4485+nuim(i)*.4516+nuin(i)*.5067)
     &        -.4*nuii(i)

          D3(i)=T2pnew(i)*(.3765*nuij(i)+.3733*nuin(i)+.3862*nuik(i)
     &                    +.3879*nuil(i) +.3871*nuim(i)
     &              +2.1828e-4*nuie(i)
     &                    +.4586*nuiN2(i)+.4606*nuiO2(i)
     &                    +.9412*Q_cross2*nuiO(i)+.5*Q_cross2*nuiH(i)
     &              +.9412*Q_cross2*nuiOHot(i))                !MZ
     &           +T1tnew(i)*.0941*nuij(i)
     &         +Tmtnew(i)*(.0552*nuik(i)+.0485*nuil(i) +.0516*nuim(i))
     &         +Tetnew(i)*  1.5991*nuie(i)
     &         +T3tnew(i)* .1067*nuin(i)
     &         +T1pnew(i)*.0235*nuij(i)
     &         +Tmpnew(i)*(.0138*nuik(i)+.0121*nuil(i) +.0129*nuim(i))
     &         +Tepnew(i)* .3998*nuie(i)
     &         +T3pnew(i)* .0267*nuin(i)
     &           +.4*nuii(i)*T2pnew(i)

        nu_omega=(nuiN2(i)+nuiO2(i)+nuiO(i)+nuiOHot(i)+
     &              nuiH(i))/omega(i)*omi    !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Ci0**2
        dUn2=(U2new(i)-Un(i)/Ci0)**2

          D3(i)=D3(i)
     &           +Tn(i)/T_0*(.069*nuiN2(i)+.0606*nuiO2(i)
     &        +nuiH(i)+.1176*nuiO(i))
     &        +TnOHot(i)/T_0*.1176*nuiOHot(i)                    !MZ
     &         +nuiN2(i)*(.4586*dUn2+1.4724*terme_joule)
     &         +nuiO2(i)*(.4606*dUn2+1.4788*terme_joule)
     &         +nuiO(i)*.9412*((2.-Q_cross2)*terme_joule+Q_cross2*dUn2)
     &      +nuiOHot(i)*.9412*((2.-Q_cross2)*terme_joule+Q_cross2*dUn2)    !MZ
     &         +nuiH(i)*.5*((2.-Q_cross2)*terme_joule+Q_cross2*dUn2)
     &         +.667*(
     &         + nuij(i)*.94*(Cji*U1new(i)-U2new(i))**2
     &         + (nuik(i)+ nuil(i)+ nuim(i))
     &           *.97*(Czi*Umnew(i)-U2new(i))*(Czi*Umnew(i)-U2new(i))
     &         + nuin(i)*.93*(Cni*U3new(i)-U2new(i))**2
     &         + nuie(i)*5.46e-4*(Cei*Uenew(i)-U2new(i))**2)
     &         - 2.4*(q2new(i)*alt_geo_1(i)*xn2_1(i))

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*(U2new(np)+U2new(nx))
      D2br=.5*(q2new(np)+q2new(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T2told,extra,nx)
      rbc=1.
       call lcpfct(T2told,T2tnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T2tnew(i)=dexpnu*(T2tnew(i)-T2told(i)+D3(i)*deltat)
     &             +expnu*T2told(i)
          enddo
          T2tnew(np)=T2tnew(nx)

          do i=1,np
            T2tnew(i)=max(T2tnew(i),T_min)
        T2new(i)=(T2pnew(i)+2.*T2tnew(i))/3.
      enddo

      if (any(ieee_is_nan(T2tnew))) then
        write(stderr,*) 'problem T2tnew dans la boucle 2'
        goto 246
      endif


C]]]


      call velocity(Velje,Ipos1,Iposnp,deltat)

C[[[    O+ energy equation resolution


c    temperature parallele

      do i=1,nx

      C2a(i)=-Cji*T1pnew(i)
      D2a(i)=U1new(i)
          C2b(i)=-1.2*Cji*xn1_1(i)
      D2b(i)=q1new(i)

        Q_cross2=max(.3,min(.45,.45-.15*(T1new(i)-2.5)))/2.
        Q_cross2=.3

          D7(i)=
     &      -(1.3318*nujN2(i)+1.3*nujO2(i)
     &        +(1.+Q_cross2)*nujO(i)
     &        +(1.+Q_cross2)*nujOHot(i)                        !MZ
     &       +nuji(i)*1.9294+nuje(i)*2.
     &       +nujk(i)*1.2364+nujl(i)*1.2+nujm(i)*1.2174+nujn(i)*1.44)
     &        -.8*nujj(i)

          D3(i)=T1tnew(i)*(.0471*nuji(i)+.3733*nujn(i)+.5091*nujk(i)
     &                    +.5333*nujl(i) +.5217*nujm(i)
     &              +2.7299e-5*nuje(i)
     &                    +.6045*nujN2(i)+.6333*nujO2(i)
     &              +Q_cross2*nujO(i)
     &              +Q_cross2*nujOHot(i))                    !MZ
     &           +T2pnew(i)*1.1294*nuji(i)
     &         +Tmpnew(i)*(.4364*nujk(i)+.4*nujl(i) +.4174*nujm(i))
     &         +Tepnew(i)*  1.2*nuje(i)
     &         +T3pnew(i)* .64*nujn(i)
     &         +T2tnew(i)*.7529*nuji(i)
     &         +Tmtnew(i)*(.2909*nujk(i)+.2667*nujl(i) +.2783*nujm(i))
     &         +Tetnew(i)*  .8*nuje(i)
     &         +T3tnew(i)* .4267*nujn(i)
     &           +.8*nujj(i)*T1tnew(i)



        nu_omega=(nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i))/omega(i)*omj        !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cj0**2
        dUn2=(U1new(i)-Un(i)/Cj0)**2

         D3(i)=D3(i)
     &           +Tn(i)/T_0*(.7273*nujN2(i)+.6667*nujO2(i)
     &        +nujO(i))
     &        +TnOHot(i)/T_0*nujOHot(i)                    !MZ
     &         +nujN2(i)*(.6682*dUn2+.6045*terme_joule)
     &         +nujO2(i)*(.7*dUn2+.6333*terme_joule)
     &         +nujO(i)*(dUn2+Q_cross2*(terme_joule-dUn2))
     &        +nujOHot(i)*(dUn2+Q_cross2*(terme_joule-dUn2))            !MZ
     &      + .667*(nuji(i)/17.*(Cij*U2new(i)-U1new(i))**2
     &      + (nujk(i)*.64+ nujl(i)*.667+ nujm(i)*.65)
     &          *(Czj*Umnew(i)-U1new(i))*(Czj*Umnew(i)-U1new(i))
     &      + nujn(i)*.467*(Cnj*U3new(i)-U1new(i))**2
     &      + nuje(i)*5.46e-4/16.*(Cej*Uenew(i)-U1new(i))**2)
     &      - 1.2*(Cji*q1new(i)*alt_geo_1(i)*xn1_1(i))


      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*(U1new(np)+U1new(nx))
      D2br=.5*(q1new(np)+q1new(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T1pold,extra,nx)
      rbc=1.
       call lcpfct(T1pold,T1pnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T1pnew(i)=dexpnu*(T1pnew(i)-T1pold(i)+D3(i)*deltat)
     &             +expnu*T1pold(i)
          enddo
          T1pnew(np)=T1pnew(nx)

      do i=1,np
        T1pnew(i)=max(T1pnew(i),T_min)
      enddo

      if (any(ieee_is_nan(T1pnew))) then
        write(stderr,*) 'problem T1pnew dans la boucle 2'
        goto 246
      endif



c    temperature perpendiculaire

      do i=1,nx

      C2a(i)=Cji*T1tnew(i)
      D2a(i)=U1new(i)
          C2b(i)=-.4*Cji*xn1_1(i)
      D2b(i)=q1new(i)

        Q_cross2=max(.3,min(.45,.45-.15*(T1new(i)-2.5)))/4.
        Q_cross2=.3

          D7(i)=
     &      -(1.0295*nujN2(i)+.9833*nujO2(i)
     &        +(1.+Q_cross2)*nujO(i)
     &        +(1.+Q_cross2)*nujOHot(i)                        !MZ
     &       +nuji(i)*1.9059+nuje(i)*1.9999
     &       +nujk(i)*.9818+nujl(i)*.9333+nujm(i)*.9565+nujn(i)*1.2533)
     &        -.4*nujj(i)

         D3(i)=T1pnew(i)*(.0235*nuji(i)+.1867*nujn(i)+.2545*nujk(i)
     &                    +.2667*nujl(i) +.2609*nujm(i)
     &              +1.365e-5*nuje(i)
     &                    +.3023*nujN2(i)+.3167*nujO2(i)
     &                    +Q_cross2*nujO(i)
     &                    +Q_cross2*nujOHot(i))                    !MZ
     &           +T2tnew(i)*1.5059*nuji(i)
     &         +Tmtnew(i)*(.5818*nujk(i)+.5333*nujl(i) +.5565*nujm(i))
     &         +Tetnew(i)*  1.5999*nuje(i)
     &         +T3tnew(i)* .8533*nujn(i)
     &         +T2pnew(i)*.3765*nuji(i)
     &         +Tmpnew(i)*(.1455*nujk(i)+.1333*nujl(i) +.1391*nujm(i))
     &         +Tepnew(i)*  .4*nuje(i)
     &         +T3pnew(i)* .2133*nujn(i)
     &           +.4*nujj(i)*T1pnew(i)



        nu_omega=(nujN2(i)+nujO2(i)+nujO(i)+nujOHot(i))/omega(i)*omj        !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cj0**2
        dUn2=(U1new(i)-Un(i)/Cj0)**2

        D3(i)=D3(i)
     &           +Tn(i)/T_0*(.7273*nujN2(i)+.6667*nujO2(i)
     &                    +nujO(i))
     &        +TnOHot(i)/T_0*nujOHot(i)                    !MZ
     &         +nujN2(i)*(.3023*dUn2+.9705*terme_joule)
     &         +nujO2(i)*(.3167*dUn2+1.0167*terme_joule)
     &         +nujO(i)*(Q_cross2*dUn2+(1.-Q_cross2)*terme_joule)
     &         +nujOHot(i)*(Q_cross2*dUn2+(1.-Q_cross2)*terme_joule)        !MZ
     &      + .667*(nuji(i)/17.*(Cij*U2new(i)-U1new(i))**2
     &      + (nujk(i)*.64+ nujl(i)*.667+ nujm(i)*.65)
     &                *(Czj*Umnew(i)-U1new(i))*(Czj*Umnew(i)-U1new(i))
     &      + nujn(i)*.467*(Cnj*U3new(i)-U1new(i))**2
     &      + nuje(i)*5.46e-4/16.*(Cej*Uenew(i)-U1new(i))**2)
     &      - 2.4*(Cji*q1new(i)*alt_geo_1(i)*xn1_1(i))

       enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*(U1new(np)+U1new(nx))
      D2br=.5*(q1new(np)+q1new(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T1pold,extra,nx)
      rbc=1.
       call lcpfct(T1told,T1tnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T1tnew(i)=dexpnu*(T1tnew(i)-T1told(i)+D3(i)*deltat)
     &             +expnu*T1told(i)
          enddo
          T1tnew(np)=T1tnew(nx)

          do  i=1,np
            T1tnew(i)=max(T1tnew(i),T_min)
        T1new(i)=(T1pnew(i)+2.*T1tnew(i))/3.
          enddo

      if (any(ieee_is_nan(T1tnew))) then
        call cpu_time(tic)
        write(stderr,*) tic,'problem when calculating T1tnew in loop 2'
        goto 246
      endif



C]]]


      call velocity(Velme,Ipos1,Iposnp,deltat)

C[[[    heavy ions energy equation resolution

c    temperature parallele

      do i=1,nx

      C2a(i)=-Czi*Tmpnew(i)
      D2a(i)=Umnew(i)

        Q_cross2=max(.3,min(.45,.45-.15*(Tmnew(i)-2.5)))/2.
        Q_cross2=.3


          D7(i)=
     &      -((1.+Q_cross2)*nukN2(i)+1.44*nukO2(i)+1.6182*nukO(i)
     &         +1.6182*nukOHot(i)                            !MZ
     &       +nuki(i)*1.9586+nukj(i)*1.5636+nuke(i)*2.
     &       +nukn(i)*1.6)
     &      -(1.51*nulN2(i)+(1.+Q_cross2)*nulO2(i)+1.65*nulO(i)
     &         +1.65*nulOHot(i)                            !MZ
     &       +nuli(i)*1.9636+nulj(i)*1.6+nule(i)*2.
     &       +nuln(i)*1.6348)
     &      -(1.4931*numN2(i)+1.4581*numO2(i)+1.6348*numO(i)
     &         +1.6348*numOHot(i)                            !MZ
     &       +numi(i)*1.9613+numj(i)*1.5826+nume(i)*2.
     &       +numn(i)*1.44)
     &       -.8*(nukk(i)+nukl(i)+nukm(i)
     &           +nulk(i)+null(i)+nulm(i)
     &           +numk(i)+numl(i)+numm(i))


        D3(i)=Tmtnew(i)*(.0276*nuki(i)+.2909*nukj(i)+.2667*nukn(i)
     &                 +1.56e-5*nuke(i)
     &                 +.0242*nuli(i)+.2667*nulj(i)+.2435*nuln(i)
     &                 + 1.365e-5*nule(i)
     &                 +.0258*numi(i)+.2783*numj(i)+.2545*numn(i)
     &                 +1.456e-5*nume(i)
     &                 +Q_cross2*nukN2(i)+.5067*nukO2(i)+.3455*nukO(i)
     &            +.3455*nukOHot(i)                    !MZ
     &                 +.4433*nulN2(i)+Q_cross2*nulO2(i)+.3167*nulO(i)
     &            +.3167*nulOHot(i)                    !MZ
     &                 +.4586*numN2(i)+.4903*numO2(i)+.3304*numO(i)
     &            +.3304*numOHot(i))                    !MZ
     &           +T2pnew(i)*(1.1586*nuki(i)
     &                    +1.1636*nuli(i)
     &                    +1.1613*numi(i))
     &         +T1pnew(i)*(.7636*nukj(i)
     &                    +.8*nulj(i)
     &                    +.7826*numj(i))
     &         +Tepnew(i)*( 1.2*nuke(i)
     &                    +1.2*nule(i)
     &                    +1.2*nume(i))
     &         +T3pnew(i)*(.8*nukn(i)
     &                    +.8348*nuln(i)
     &                    +.8182*numn(i))
     &         +T2tnew(i)*(.7724*nuki(i)
     &                    +.7758*nuli(i)
     &                    +.7742*numi(i))
     &         +T1tnew(i)*(.5091*nukj(i)
     &                    +.5333*nulj(i)
     &                    +.5217*numj(i))
     &         +Tetnew(i)*( .8*nuke(i)
     &                    +.8*nule(i)
     &                    +.8*nume(i))
     &         +T3tnew(i)*(.5333*nukn(i)
     &                    +.5565*nuln(i)
     &                    +.5455*numn(i))
     &         +.8*(nukk(i)+nukl(i)+nukm(i)
     &             +nulk(i)+null(i)+nulm(i)
     &             +numk(i)+numl(i)+numm(i))*Tmtnew(i)

        nu_omega=(nukN2(i)+nukO2(i)+nukO(i)+nukOHot(i)                !MZ
     &               +nulN2(i)+nulO2(i)+nulO(i)+nulOHot(i)            !MZ
     &              +numN2(i)+numO2(i)+numO(i)+numOHot(i))/omega(i)*omz    !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cz0**2
        dUn2=(Umnew(i)-Un(i)/Cz0)**2

        D3(i)=D3(i)
     &           +Tn(i)/T_0
     &        *(nukN2(i)+.9333*nukO2(i)+1.2727*nukO(i)
     &         +1.0667*nulN2(i)+nulO2(i)+1.3333*nulO(i)
     &         +1.0345*numN2(i)+0.9677*numO2(i)+1.3043*numO(i))
     &        +TnOHot(i)/T_0*(1.2727*nukOHot(i)+1.3333*nulOHot(i)        !MZ
     &            +1.3043*numOHot(i))                    !MZ
     &         +nukN2(i)*(dUn2+Q_cross2*(terme_joule-dUn2))
     &         +nukO2(i)*(.5242*dUn2+.4743*terme_joule)
     &         +nukO(i)*(.3574*dUn2+.3234*terme_joule)
     &         +nukOHot(i)*(.3574*dUn2+.3234*terme_joule)            !MZ
     &         +nulN2(i)*(.5242*dUn2+.4743*terme_joule)
     &         +nulO2(i)*(dUn2+Q_cross2*(terme_joule-dUn2))
     &         +nulO(i)*(.3744*dUn2+.3388*terme_joule)
     &         +nulOHot(i)*(.3744*dUn2+.3388*terme_joule)            !MZ
     &         +numN2(i)*(.5084*dUn2+.46*terme_joule)
     &         +numO2(i)*(.5435*dUn2+.4918*terme_joule)
     &         +numO(i)*(.3663*dUn2+.3314*terme_joule)
     &         +numOHot(i)*(.3663*dUn2+.3314*terme_joule)            !MZ
     &      + .22*(        (nuki(i)/29.+nuli(i)/33.+numi(i)/31.)
     &                       *(Ciz*U2new(i)-Umnew(i))**2
     &            +    16.*(nukj(i)/44.+nulj(i)/48.+numj(i)/46.)
     &                   *(Cjz*U1new(i)-Umnew(i))**2
     &            +    14.*(nukn(i)/42.+nuln(i)/46.+numn(i)/44.)
     &                   *(Cnz*U3new(i)-Umnew(i))**2
     &            +5.46e-4*(nuke(i)/28.+nule(i)/32.+nume(i)/30.)
     &               *(Cez*Uenew(i)-Umnew(i))**2
     &          )


      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2ar=(D2a(nx)+D2a(nx-1))/2.
    
      D2ar=.5*(Umnew(np)+Umnew(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
      rbc=bclimd(Radn(nx+1),Tmpold,extra,nx)
      rbc=1.
       call lcpfct(Tmpold,Tmpnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Tmpnew(i)=dexpnu*(Tmpnew(i)-Tmpold(i)+D3(i)*deltat)
     &             +expnu*Tmpold(i)
          enddo
          Tmpnew(np)=Tmpnew(nx)

      do i=1,np
         Tmpnew(i)=max(Tmpnew(i),T_min)
          enddo
c    flag=.false.

      if (any(ieee_is_nan(Tmpnew))) then
        call cpu_time(tic)
        write(stderr,*) tic,'problem when calculating Tmpnew in loop 2'
        goto 246
      endif

c    temperature perpendiculaire

      do i=1,nx

      C2a(i)=Czi*Tmtnew(i)
      D2a(i)=Umnew(i)

        Q_cross2=max(.3,min(.45,.45-.15*(Tmnew(i)-2.5)))/4.
        Q_cross2=.3


          D7(i)=
     &      -((1.+Q_cross2)*nukN2(i)+1.1867*nukO2(i)+1.4455*nukO(i)
     &        +1.4455*nukOHot(i)                        !MZ
     &       +nuki(i)*1.9448+nukj(i)*1.4182+nuke(i)*2.
     &       +nukn(i)*1.4667)
     &      -(1.2883*nulN2(i)+(1.+Q_cross2)*nulO2(i)+1.4917*nulO(i)
     &        +1.4917*nulOHot(i)                        !MZ
     &       +nuli(i)*1.9515+nulj(i)*1.4667+nule(i)*2.
     &       +nuln(i)*1.513)
     &      -(1.2638*numN2(i)+1.2129*numO2(i)+1.4696*numO(i)
     &        +1.4696*numOHot(i)                        !MZ
     &       +numi(i)*1.9484+numj(i)*1.4435+nume(i)*2.
     &       +numn(i)*1.4909)
     &       -.4*(nukk(i)+nukl(i)+nukm(i)
     &           +nulk(i)+null(i)+nulm(i)
     &           +numk(i)+numl(i)+numm(i))

        D3(i)=Tmpnew(i)*(.0138*nuki(i)+.1455*nukj(i)+.1333*nukn(i)
     &                 +7.7998e-6*nuke(i)
     &                 +.0121*nuli(i)+.1333*nulj(i)+.1217*nuln(i)
     &                 +6.8249e-6*nule(i)
     &                 +.0129*numi(i)+.1391*numj(i)+.1273*numn(i)
     &                 +7.2799e-6*nume(i)
     &               +Q_cross2*nukN2(i)+.2533*nukO2(i)+.1727*nukO(i)
     &            +.1727*nukOHot(i)                    !MZ
     &               +.2217*nulN2(i)+Q_cross2*nulO2(i)+.1583*nulO(i)
     &            +.1583*nulOHot(i)                    !MZ
     &               +.2293*numN2(i)+.2452*numO2(i)+.1652*numO(i)
     &            +.1652*numOHot(i))                    !MZ
     &        +T2tnew(i)*(1.5448*nuki(i)+1.5515*nuli(i)+1.5484*numi(i))
     &        +T1tnew(i)*(1.0182*nukj(i)+1.0667*nulj(i)+1.0435*numj(i))
     &         +Tetnew(i)*( 1.6*nuke(i)+1.6*nule(i)+1.6*nume(i))
     &        +T3tnew(i)*(1.0667*nukn(i)+1.1130*nuln(i)+1.0909*numn(i))
     &         +T2pnew(i)*(.3862*nuki(i)+.3879*nuli(i)+.3871*numi(i))
     &         +T1pnew(i)*(.2545*nukj(i)+.2667*nulj(i)+.2609*numj(i))
     &         +Tepnew(i)*(.4*nuke(i)+.4*nule(i)+.4*nume(i))
     &         +T3pnew(i)*(.2667*nukn(i)+.2783*nuln(i)+.2727*numn(i))
     &         +.4*(nukk(i)+nukl(i)+nukm(i)
     &             +nulk(i)+null(i)+nulm(i)
     &             +numk(i)+numl(i)+numm(i))*Tmpnew(i)

         nu_omega=(nukN2(i)+nukO2(i)+nukO(i)+nukOHot(i)                !MZ
     &               +nulN2(i)+nulO2(i)+nulO(i)+nulOHot(i)            !MZ
     &              +numN2(i)+numO2(i)+numO(i)+numOHot(i))/omega(i)*omz    !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cz0**2
        dUn2=(Umnew(i)-Un(i)/Cz0)**2

        D3(i)=D3(i)
     &           +Tn(i)/T_0
     &        *(nukN2(i)+.9333*nukO2(i)+1.2727*nukO(i)
     &         +1.0667*nulN2(i)+nulO2(i)+1.3333*nulO(i)
     &         +1.0345*numN2(i)+0.9677*numO2(i)+1.3043*numO(i))
     &        +TnOHot(i)/T_0*(1.2727*nukOHot(i)+1.3333*nulOHot(i)        !MZ
     &            +1.3043*numOHot(i))                    !MZ
     &         +nukN2(i)*(Q_cross2*dUn2+(1.-Q_cross2)*terme_joule)
     &         +nukO2(i)*(.2371*dUn2+.7614*terme_joule)
     &         +nukO(i)*(.1617*dUn2+.5191*terme_joule)
     &         +nukOHot(i)*(.1617*dUn2+.5191*terme_joule)            !MZ
     &         +nulN2(i)*(.2371*dUn2+.7614*terme_joule)
     &         +nulO2(i)*(Q_cross2*dUn2+(1.-Q_cross2)*terme_joule)
     &         +nulO(i)*(.1694*dUn2+.5438*terme_joule)
     &         +nulOHot(i)*(.1694*dUn2+.5438*terme_joule)            !MZ
     &         +numN2(i)*(.23*dUn2+.7384*terme_joule)
     &         +numO2(i)*(.1657*dUn2+.7894*terme_joule)
     &         +numO(i)*(.3663*dUn2+.532*terme_joule)
     &         +numOHot(i)*(.3663*dUn2+.532*terme_joule)            !MZ
     &      + .22*(        (nuki(i)/29.+nuli(i)/33.+numi(i)/31.)
     &                       *(Ciz*U2new(i)-Umnew(i))**2
     &            +    16.*(nukj(i)/44.+nulj(i)/48.+numj(i)/46.)
     &                   *(Cjz*U1new(i)-Umnew(i))**2
     &            +    14.*(nukn(i)/42.+nuln(i)/46.+numn(i)/44.)
     &                   *(Cnz*U3new(i)-Umnew(i))**2
     &            +5.46e-4*(nuke(i)/28.+nule(i)/32.+nume(i)/30.)
     &               *(Cez*Uenew(i)-Umnew(i))**2
     &          )

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
    
      D2ar=.5*(Umnew(np)+Umnew(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
      rbc=bclimd(Radn(nx+1),Tmpold,extra,nx)
      rbc=1.
       call lcpfct(Tmtold,Tmtnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            Tmtnew(i)=dexpnu*(Tmtnew(i)-Tmtold(i)+D3(i)*deltat)
     &             +expnu*Tmtold(i)
          enddo
          Tmtnew(np)=Tmtnew(nx)

      do i=1,np
         Tmtnew(i)=max(Tmtnew(i),T_min)
        Tmnew(i)=(Tmpnew(i)+2.*Tmtnew(i))/3.
          enddo

      if (any(ieee_is_nan(Tmtnew))) then
        write(stderr,*) 'problem calc Tmtnew loop 2'
        goto 246
      endif


C]]]

      call velocity(Velne,Ipos1,Iposnp,deltat)

C[[[    N+ energy equation resolution

c    temperature parallele

        do i=1,nx

      C2a(i)=-Cni*T3pnew(i)
      D2a(i)=U3new(i)
      C2b(i)=-1.2*Cni*xn3_1(i)
      D2b(i)=q3new(i)

          D7(i)=
     &      -(1.3*nunN2(i)+1.2696*nunO2(i)+1.44*nunO(i)
     &        +1.44*nunOHot(i)                        !MZ
     &       +nuni(i)*1.92+nunj(i)*1.36+nune(i)*2.
     &       +nunk(i)*1.2+nunl(i)*1.1652+nunm(i)*1.1818)
     &        -.8*nunn(i)

          D3(i)=T3tnew(i)*(.0533*nuni(i)+.4267*nunj(i)+.5333*nunk(i)
     &                    +.5565*nunl(i) +.5455*nunm(i)
     &              +3.1199e-5*nune(i)
     &                    +.6333*nunN2(i)+.6609*nunO2(i)+.5067*nunO(i)
     &              +.5067*nunOHot(i))                    !MZ
     &           +T1pnew(i)*.56*nunj(i)
     &         +T2pnew(i)* 1.12*nuni(i)
     &         +Tmpnew(i)*(.4*nunk(i)+.3652*nunl(i) +.3818*nunm(i))
     &         +Tepnew(i)*  1.2*nune(i)
     &         +T1tnew(i)* .3733*nunj(i)
     &         +T2tnew(i)*.7467*nuni(i)
     &         +Tmtnew(i)*(.2667*nunk(i)+.2435*nunl(i) +.2545*nunm(i))
     &         +Tetnew(i)*  .8*nune(i)
     &           +.8*nunn(i)*T3tnew(i)

        nu_omega=(nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i))/omega(i)*omn        !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cn0**2
        dUn2=(U3new(i)-Un(i)/Cn0)**2

         D3(i)=D3(i)
     &           +Tn(i)/T_0*(.6667*nunN2(i)+.6087*nunO2(i)
     &                    +.9333*nunO(i))
     &        +TnOHot(i)/T_0*.9333*nunOHot(i)                    !MZ
     &         +nunN2(i)*(.7*dUn2+.6333*terme_joule)
     &         +nunO2(i)*(.7304*dUn2+.6609*terme_joule)
     &         +nunO(i)*(.56*dUn2+.5067*terme_joule)
     &         +nunOHot(i)*(.56*dUn2+.5067*terme_joule)                !MZ
     &      + .667*(nuni(i)*.067*(Cin*U2new(i)-U3new(i))**2
     &               + nunk(i)*.667*(Czn*Umnew(i)-U3new(i))**2
     &               + nunl(i)*.7*(Czn*Umnew(i)-U3new(i))**2
     &               + nunm(i)*.68*(Czn*Umnew(i)-U3new(i))**2
     &               + nunj(i)*.53*(Cjn*U1new(i)-U3new(i))**2
     &               + nune(i)*3.9e-5*(Cen*Uenew(i)-U3new(i))**2)
     &               - 1.2*(Cni*q3new(i)*alt_geo_1(i)*xn3_1(i))

       enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*(U3new(np)+U3new(nx))
      D2br=.5*(q3new(np)+q3new(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T3pold,extra,nx)
      rbc=1.
       call lcpfct(T3pold,T3pnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T3pnew(i)=dexpnu*(T3pnew(i)-T3pold(i)+D3(i)*deltat)
     &             +expnu*T3pold(i)
          enddo
          T3pnew(np)=T3pnew(nx)

      do i=1,np
            T3pnew(i)=max(T3pnew(i),T_min)
c            T3pnew(i)=T1pnew(i)
      enddo

      if (any(ieee_is_nan(T3pnew))) then
        call cpu_time(tic)
        write(stderr,*) tic,'problem when calculating T3pnew in loop 2'
        goto 246
      endif


c    temperature perpendiculaire

      do i=1,nx

      C2a(i)=Cni*T3tnew(i)
      D2a(i)=U3new(i)
      C2b(i)=-.4*Cni*xn3_1(i)
      D2b(i)=q3new(i)

          D7(i)=
     &      -(.9833*nunN2(i)+.9391*nunO2(i)+1.1867*nunO(i)
     &        +1.1867*nunOHot(i)                        !MZ
     &       +nuni(i)*1.8933+nunj(i)*1.1467+nune(i)*1.9999
     &       +nunk(i)*.9333+nunl(i)*.887+nunm(i)*.9091)
     &        -.4*nunn(i)

         D3(i)=T3pnew(i)*(.0267*nuni(i)+.2133*nunj(i)+.2667*nunk(i)
     &                    +.2783*nunl(i) +.2727*nunm(i)
     &              +1.5599e-5*nune(i)
     &                    +.3167*nunN2(i)+.3304*nunO2(i)+.2533*nunO(i)
     &              +.2533*nunOHot(i))                    !MZ
     &           +T1tnew(i)*.7467*nunj(i)
     &         +T2tnew(i)* 1.4933*nuni(i)
     &         +Tmtnew(i)*(.5333*nunk(i)+.487*nunl(i) +.5091*nunm(i))
     &         +Tetnew(i)*  1.5999*nune(i)
     &         +T1tnew(i)*.1867*nunj(i)
     &         +T2pnew(i)*.3733*nuni(i)
     &         +Tmpnew(i)*(.1333*nunk(i)+.1217*nunl(i) +.1273*nunm(i))
     &         +Tepnew(i)*  .4*nune(i)
     &           +.4*nunn(i)*T3pnew(i)

        nu_omega=(nunN2(i)+nunO2(i)+nunO(i)+nunOHot(i))/omega(i)*omn        !MZ
        coef_conv=1./(1.+nu_omega**2)
        terme_Joule=.5*coef_conv*Vm_2(i)/Cn0**2
        dUn2=(U3new(i)-Un(i)/Cn0)**2

         D3(i)=D3(i)
     &           +Tn(i)/T_0*(.6667*nunN2(i)+.6087*nunO2(i)
     &                    +.9333*nunO(i))
     &        +TnOHot(i)/T_0*.9333*nunOHot(i)                    !MZ
     &         +nunN2(i)*(.3167*dUn2+1.0167*terme_joule)
     &         +nunO2(i)*(.3304*dUn2+1.0609*terme_joule)
     &         +nunO(i)*(.2533*dUn2+.8133*terme_joule)
     &         +nunOHot(i)*(.2533*dUn2+.8133*terme_joule)            !MZ
     &      + .667*(nuni(i)*.067*(Cin*U2new(i)-U3new(i))**2
     &               + nunk(i)*.667*(Czn*Umnew(i)-U3new(i))**2
     &               + nunl(i)*.7*(Czn*Umnew(i)-U3new(i))**2
     &               + nunm(i)*.68*(Czn*Umnew(i)-U3new(i))**2
     &               + nunj(i)*.53*(Cjn*U1new(i)-U3new(i))**2
     &               + nune(i)*3.9e-5*(Cen*Uenew(i)-U3new(i))**2)
     &               - 2.4*(Cni*q3new(i)*alt_geo_1(i)*xn3_1(i))

      enddo

          D2al=(D2a(1)+D2a(2))/2.
          D2ar=ylimd(Radn(np),D2a,extra,nx)
          D2bl=(D2b(1)+D2b(2))/2.
          D2br=ylimd(Radn(np),D2b,extra,nx)
    
      D2ar=.5*(U3new(np)+U3new(nx))
      D2br=.5*(q3new(np)+q3new(nx))
    
      call sources(Ipos1,Iposn,deltat,2,C2a,D2a,D2al,D2ar)
      call sources(Ipos1,Iposn,deltat,2,C2b,D2b,D2bl,D2br)
c      call sources(Ipos1,Iposn,deltat,3,zero,D3,0.,0.)
      call sources(Ipos1,Iposn,deltat,7,zero,D7,0.,0.)
          lbc=1.
          rbc=bclimd(Radn(nx+1),T3told,extra,nx)
      rbc=1.
       call lcpfct(T3told,T3tnew,Ipos1,Iposn,
     &              lbc,0.,rbc,0.,.false.,0)
          do i=1,nx
            dexpnu=D7(i)*deltat
            if (abs(dexpnu).lt.1.e-7) then
              expnu=1.+dexpnu
              dexpnu=1.
            else
              expnu=exp(dexpnu)
              dexpnu=(expnu-1.)/dexpnu
            endif
            T3tnew(i)=dexpnu*(T3tnew(i)-T3told(i)+D3(i)*deltat)
     &             +expnu*T3told(i)
          enddo
          T3tnew(np)=T3tnew(nx)

      do i=1,np
            T3tnew(i)=max(T3tnew(i),T_min)
c            T3tnew(i)=T1tnew(i)
        T3new(i)=(T3pnew(i)+2.*T3tnew(i))/3.
      enddo

      if (any(ieee_is_nan(T3tnew))) then
        call cpu_time(tic)
        write(stderr,*) tic,'problem when calculating T3tnew in loop 2'
        goto 246
      endif

456   continue

C]]]
        enddo                    ! fin de boucle temporelle
    
      if (flgconv) then

      if (abs(duree-tempsort).le.sortie/2.) then
        nrectemps=nrectemps+1
            do i=1,longbuf
              buffer(i)=0.
            enddo
            buffer( 1)=nx
            buffer( 2)=ncol
            buffer( 3)=iannee
            buffer( 4)=imois
            buffer( 5)=ijour
            buffer( 6)=iheure
            buffer( 7)=iminute
            buffer( 8)=seconde
            buffer( 9)=intpas
            buffer(10)=longeo
            buffer(11)=latgeo
            buffer(12)=lonmag
            buffer(13)=latmag
            buffer(14)=tmag
            buffer(15)=f107(2)
            buffer(16)=f107(3)
            buffer(17)=ap(2)
            buffer(18)=kptime
            buffer(19)=dTinf
            buffer(20)=dUinf
            buffer(21)=cofo
            buffer(22)=cofh
            buffer(23)=cofn
            buffer(24)=chi0
            buffer(25)=Fe0
            buffer(26)=Ee0
            buffer(27)=Fi0
            buffer(28)=Ei0
            buffer(29)=Bmag
            buffer(30)=dipangle
            buffer(31)=Enord
            buffer(32)=Eest
            buffer(33)=vperpnord
            buffer(34)=vperpest
            buffer(35)=vhorizon
            buffer(36)=vpara
        buffer(37)=iapprox
        buffer(38)=ddp
        buffer(39)=Jtop
            do i=1,nx
              ipos=(i+1)*ncol
              buffer(ipos+ipos_z)=alt(i)
              buffer(ipos+ipos_n1)=N_0*N1new(i)*1.e6
              buffer(ipos+ipos_n2)=N_0*N2new(i)*1.e6
              buffer(ipos+ipos_n3)=N_0*N3new(i)*1.e6
              buffer(ipos+ipos_n4)=N_0*N4new(i)*1.e6
              buffer(ipos+ipos_n5)=N_0*N5new(i)*1.e6
              buffer(ipos+ipos_n6)=N_0*N6new(i)*1.e6
              buffer(ipos+ipos_u1) =Cj0*U1new(i)/1.e2
              buffer(ipos+ipos_u2) =Ci0*U2new(i)/1.e2
              buffer(ipos+ipos_u3)=Cn0*U3new(i)/1.e2
              buffer(ipos+ipos_um)=Cz0*Umnew(i)/1.e2
              buffer(ipos+ipos_ue)=Ce0*Uenew(i)/1.e2
              buffer(ipos+ipos_t1p)=T_0*T1pnew(i)
              buffer(ipos+ipos_t1t)=T_0*T1tnew(i)
              buffer(ipos+ipos_t2p)=T_0*T2pnew(i)
              buffer(ipos+ipos_t2t)=T_0*T2tnew(i)
              buffer(ipos+ipos_t3p)=T_0*T3pnew(i)
              buffer(ipos+ipos_t3t)=T_0*T3tnew(i)
              buffer(ipos+ipos_tmp)=T_0*Tmpnew(i)
              buffer(ipos+ipos_tmt)=T_0*Tmtnew(i)
              buffer(ipos+ipos_tep)=T_0*Tepnew(i)
              buffer(ipos+ipos_tet)=T_0*Tetnew(i)
              buffer(ipos+ipos_q1)=Qj0*q1new(i)/1.e3
              buffer(ipos+ipos_q2)=Qi0*q2new(i)/1.e3
              buffer(ipos+ipos_q3)=Qn0*q3new(i)/1.e3
              buffer(ipos+ipos_qe)=Qe0*qenew(i)/1.e3
              buffer(ipos+ipos_nno)=N_0*Nnonew(i)*1.e6
              buffer(ipos+ipos_uno)=Cm0*Unonew(i)/1.e2
              buffer(ipos+ipos_po)=Po(i)*1.e6
              buffer(ipos+ipos_ph)=Ph(i)*1.e6
              buffer(ipos+ipos_pn)=Pn(i)*1.e6
          buffer(ipos+ipos_pn2)=Pn2(i)*1.e6
          buffer(ipos+ipos_po2)=Po2(i)*1.e6
          buffer(ipos+ipos_heat)=Heat(i)/10.
          buffer(ipos+ipos_no)=No(i)*1.e6
          buffer(ipos+ipos_nh)=Nh(i)*1.e6
          buffer(ipos+ipos_nn)=Nn(i)*1.e6
          buffer(ipos+ipos_nn2)=Nn2(i)*1.e6
          buffer(ipos+ipos_no2)=No2(i)*1.e6
          buffer(ipos+ipos_tn)=Tn(i)
          buffer(ipos+ipos_un)=Un(i)/1.e2
          buffer(ipos+ipos_vn)=Vn(i)/1.e2
          buffer(ipos+ipos_wn)=Wn(i)/1.e2
          buffer(ipos+ipos_nes)=Nes(i)*1.e6
          buffer(ipos+ipos_jes)=Jes(i)*1.e4
          buffer(ipos+ipos_tes)=Tes(i)
          buffer(ipos+ipos_qes)=Qes(i)*1.e-7
          buffer(ipos+ipos_NOHot)=NOHot(i)*1.e6
          buffer(ipos+ipos_TnOHot)=phdisso2(i)*1.e6

              buffer(ipos+ipos_po1d)=Po1d(i)*1.e6
              buffer(ipos+ipos_no1d)=N_0*No1dnew(i)*1.e6
              buffer(ipos+ipos_uo1d)=Cj0*Uo1dnew(i)/1.e2
            enddo
            write(fid_temp,rec=nrectemps)(buffer(i),i=1,longbuf)
          endif
       close(fid_temp)
      nrec_ecr=nrec_ecr+1
        write(unfic_out_transcar,rec=nrec_ecr)(buffer(i),i=1,longbuf)
      endif
        enddo           ! Modif DA 0202 2001: fin de boucle sur i_position
    
      enddo        ! fin de boucle sur les lignes de champ

        close(unfic_out_transcar)
        close(unfic_in_transcar)
        close(fid_temp)

      stop 'fin normale'
    

C    on a debranche ici a cause d'un probleme de NaN
246    continue
        print*,'something caused a NaN'
        print*,'pas de temps et params de norm',dt,deltat,R0,Ci0
        close(unfic_out_transcar)
        close(unfic_in_transcar)
        close(fid_temp)
        fid_NaN=unfic_out_transcar
        open(fid_NaN,file='dir.output/transcar.dump',
     &            form='unformatted',
     &                  access='direct',status='replace',recl=longrec)

      do i=1,longbuf
        buffer(i)=0.
      enddo
      buffer( 1)=nx
      buffer( 2)=ncol
      buffer( 3)=iannee
      buffer( 4)=imois
      buffer( 5)=ijour
      buffer( 6)=iheure
      buffer( 7)=iminute
      buffer( 8)=seconde
      buffer( 9)=intpas
      buffer(10)=longeo
      buffer(11)=latgeo
      buffer(12)=lonmag
      buffer(13)=latmag
      buffer(14)=tmag
      buffer(15)=f107(2)
      buffer(16)=f107(3)
      buffer(17)=ap(2)
      buffer(18)=kptime
      buffer(19)=dTinf
      buffer(20)=dUinf
      buffer(21)=cofo
      buffer(22)=cofh
      buffer(23)=cofn
      buffer(24)=chi0
      buffer(25)=Fe0
      buffer(26)=Ee0
      buffer(27)=Fi0
      buffer(28)=Ei0
      buffer(29)=Bmag
      buffer(30)=dipangle
      buffer(31)=Enord
      buffer(32)=Eest
      buffer(33)=vperpnord
      buffer(34)=vperpest
      buffer(35)=vhorizon
      buffer(36)=vpara
      buffer(37)=iapprox
      buffer(38)=ddp
      buffer(39)=Jtop
      do i=1,nx
        ipos=(i+1)*ncol
        buffer(ipos+ipos_z)=alt(i)
        buffer(ipos+ipos_n1)=N_0*N1old(i)*1.e6
        buffer(ipos+ipos_n2)=N_0*N2old(i)*1.e6
        buffer(ipos+ipos_n3)=N_0*N3old(i)*1.e6
        buffer(ipos+ipos_n4)=N_0*N4old(i)*1.e6
        buffer(ipos+ipos_n5)=N_0*N5old(i)*1.e6
        buffer(ipos+ipos_n6)=N_0*N6old(i)*1.e6
        buffer(ipos+ipos_u1) =Cj0*U1old(i)/1.e2
        buffer(ipos+ipos_u2) =Ci0*U2old(i)/1.e2
        buffer(ipos+ipos_u3)=Cn0*U3old(i)/1.e2
        buffer(ipos+ipos_um)=Cz0*Umold(i)/1.e2
        buffer(ipos+ipos_ue)=Ce0*Ueold(i)/1.e2
        buffer(ipos+ipos_t1p)=T_0*T1pold(i)
        buffer(ipos+ipos_t1t)=T_0*T1told(i)
        buffer(ipos+ipos_t2p)=T_0*T2pold(i)
        buffer(ipos+ipos_t2t)=T_0*T2told(i)
        buffer(ipos+ipos_t3p)=T_0*T3pold(i)
        buffer(ipos+ipos_t3t)=T_0*T3told(i)
        buffer(ipos+ipos_tmp)=T_0*Tmpold(i)
        buffer(ipos+ipos_tmt)=T_0*Tmtold(i)
        buffer(ipos+ipos_tep)=T_0*Tepold(i)
        buffer(ipos+ipos_tet)=T_0*Tetold(i)
        buffer(ipos+ipos_q1)=Qj0*q1old(i)/1.e3
        buffer(ipos+ipos_q2)=Qi0*q2old(i)/1.e3
        buffer(ipos+ipos_q3)=Qn0*q3old(i)/1.e3
        buffer(ipos+ipos_qe)=Qe0*qeold(i)/1.e3
        buffer(ipos+ipos_nno)=N_0*Nnoold(i)*1.e6
        buffer(ipos+ipos_uno)=Cm0*Unoold(i)/1.e2
        buffer(ipos+ipos_po)=Po(i)*1.e6
        buffer(ipos+ipos_ph)=Ph(i)*1.e6
        buffer(ipos+ipos_pn)=Pn(i)*1.e6
        buffer(ipos+ipos_pn2)=Pn2(i)*1.e6
        buffer(ipos+ipos_po2)=Po2(i)*1.e6
        buffer(ipos+ipos_heat)=Heat(i)/10.
        buffer(ipos+ipos_no)=No(i)*1.e6
        buffer(ipos+ipos_nh)=Nh(i)*1.e6
        buffer(ipos+ipos_nn)=Nn(i)*1.e6
        buffer(ipos+ipos_nn2)=Nn2(i)*1.e6
        buffer(ipos+ipos_no2)=No2(i)*1.e6
        buffer(ipos+ipos_tn)=Tn(i)
        buffer(ipos+ipos_un)=Un(i)/1.e2
        buffer(ipos+ipos_vn)=Vn(i)/1.e2
        buffer(ipos+ipos_wn)=Wn(i)/1.e2
        buffer(ipos+ipos_nes)=Nes(i)*1.e6
        buffer(ipos+ipos_jes)=Jes(i)*1.e4
        buffer(ipos+ipos_tes)=Tes(i)
        buffer(ipos+ipos_qes)=Qes(i)*1.e-7
      enddo
      write(fid_NaN,rec=1)(buffer(i),i=1,longbuf)
      close(fid_NaN)
      error stop 'NaN detected'
      end program 


      pure real function signe(x)
      implicit none
      real,intent(in) :: x

      if (x >= 0.) then
        signe=1.
      else
        signe=-1.
      endif

      end function signe


      ! split a string into 2 either side of a delimiter token
      character(80) function split(instr,  delm)
        implicit none
        CHARACTER(80),intent(in) :: instr
        character,intent(in) :: delm
        INTEGER :: idx

        idx = scan(instr,delm)
        split = instr(1:idx-1)

      END function split
