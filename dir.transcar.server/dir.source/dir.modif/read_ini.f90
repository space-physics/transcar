subroutine read_ini



tempsconv_1=0.d0
tempsconv=0.d0
vparaB=0.


!	machine SGI
itype=1
!	machine IBM
itype=4
	
call pwd(chemin,lpath,data_path,lpath_data)

!	lecture des caracteristiques de la simulation dans le fichier DATCAR
!	reading the characteristics of the simulation in the DATCAR file
open(transcar_dat, file=chemin(1:lpath)//'dir.input/DATCAR')
rewind(transcar_dat)

read(transcar_dat,*)kiappel			!
read(transcar_dat,'(a)')filein                   !
read(transcar_dat,*)dto     			! pas d'integration numerique. no numerical integration
read(transcar_dat,*)sortie      		! intervalle de temps entre deux sorties
read(transcar_dat,*)iyd_ini     		! date de la periode simulee
read(transcar_dat,*)tempsini			! UT de debut (en secondes)
read(transcar_dat,*)tempslim    		! UT limite (en secondes)
read(transcar_dat,*)jpreci                      !
read(transcar_dat,*)latgeo_ini,longeo_ini       !
read(transcar_dat,*)tempsconv_1			! duree de la convection en amont en secondes (<= 0 si pas de convection)
read(transcar_dat,*)tempsconv			! duree de la convection en aval en secondes (<= 0 si pas de convection)
read(transcar_dat,*)step        		! intervalle de temps entre deux tubes
read(transcar_dat,*)postinto    		! intervalle de temps entre deux appel a transelec
read(transcar_dat,*)vparaB      		! transport induit le long de la ligne de champ

if (latgeo_ini.gt.90.) then
  nb_position=0
  multi_position=.true.
  read(transcar_dat,*)latgeo_ini,longeo_ini
  do while (latgeo_ini.le.90.)
    nb_position=nb_position+1
    latgeo_position(nb_position)=latgeo_ini
    longeo_position(nb_position)=longeo_ini
    lecture_lat_lon=' '
    read(transcar_dat,*)latgeo_ini,longeo_ini
  enddo
else
  nb_position=1
  multi_position=.false.
endif

close(transcar_dat)

!	fin de lecture

!	definition des parametres de simulation

if (filein(1:1).ne.'/') then
  filein=chemin(1:lpath)//'dir.input/'//filein(1:lenc(filein))
endif
filein = filein(1:lenc(filein))

!	definition des indices pour lecture

call simu_input(unfic_in_transcar,filein,1,4*2*ncol0,2*ncol0,buffer)

nb_alt  =buffer( 1)
ncol_i    =buffer( 2)
nligne=nb_alt+2
longbuf_i=nligne*ncol_i
longrec_i=itype*longbuf

iapprox_i =buffer(37)
ipos_i%z=1
ipos_i%n1=2
ipos_i%n2=3
ipos_i%n3=4
ipos_i%n4=5
ipos_i%n5=6
ipos_i%n6=7
ipos_i%u1=8
ipos_i%u2=9
ipos_i%u3=10
ipos_i%um=11
ipos_i%ue=12
if (iapprox.eq.13) then
  ipos_i%t1p=13
  ipos_i%t1t=14
  ipos_i%t2p=15
  ipos_i%t2t=16
  ipos_i%t3p=17
  ipos_i%t3t=18
  ipos_i%tmp=19
  ipos_i%tmt=20
  ipos_i%tep=21
  ipos_i%tet=22
  ipos_i%q1=23
  ipos_i%q2=24
  ipos_i%q3=25
  ipos_i%qe=26
  ipos_i%nno=27
  ipos_i%uno=28
  ipos_i%po=29
  ipos_i%ph=30
  ipos_i%pn=31
  ipos_i%pn2=32
  ipos_i%po2=33
  ipos_i%heat=34
  ipos_i%no=35
  ipos_i%nh=36
  ipos_i%nn=37
  ipos_i%nn2=38
  ipos_i%no2=39
  ipos_i%tn=40
  ipos_i%un=41
  ipos_i%vn=42
  ipos_i%wn=43
  ipos_i%nes=44
  ipos_i%jes=45
  ipos_i%tes=46
  ipos_i%qes=47
else
  ipos_i%t1p=13
  ipos_i%t1t=13
  ipos_i%t2p=14
  ipos_i%t2t=14
  ipos_i%t3p=15
  ipos_i%t3t=15
  ipos_i%tmp=16
  ipos_i%tmt=16
  ipos_i%tep=17
  ipos_i%tet=17
  ipos_i%q1=18
  ipos_i%q2=19
  ipos_i%q3=20
  ipos_i%qe=21
  ipos_i%nno=22
  ipos_i%uno=23
  ipos_i%po=24
  ipos_i%ph=25
  ipos_i%pn=26
  ipos_i%pn2=27
  ipos_i%po2=28
  ipos_i%heat=29
  ipos_i%no=30
  ipos_i%nh=31
  ipos_i%nn=32
  ipos_i%nn2=33
  ipos_i%no2=34
  ipos_i%tn=35
  ipos_i%un=36
  ipos_i%vn=37
  ipos_i%wn=38
  ipos_i%nes=39
  ipos_i%jes=40
  ipos_i%tes=41
  ipos_i%qes=42
endif

ncol_o=47
longbuf_o=nligne*ncol_o
longrec_o=itype*longbuf_o

iapprox_o=13

ipos_o%t1p=13
ipos_o%t1t=14
ipos_o%t2p=15
ipos_o%t2t=16
ipos_o%t3p=17
ipos_o%t3t=18
ipos_o%tmp=19
ipos_o%tmt=20
ipos_o%tep=21
ipos_o%tet=22
ipos_o%q1=23
ipos_o%q2=24
ipos_o%q3=25
ipos_o%qe=26
ipos_o%nno=27
ipos_o%uno=28
ipos_o%po=29
ipos_o%ph=30
ipos_o%pn=31
ipos_o%pn2=32
ipos_o%po2=33
ipos_o%heat=34
ipos_o%no=35
ipos_o%nh=36
ipos_o%nn=37
ipos_o%nn2=38
ipos_o%no2=39
ipos_o%tn=40
ipos_o%un=41
ipos_o%vn=42
ipos_o%wn=43
ipos_o%nes=44
ipos_o%jes=45
ipos_o%tes=46
ipos_o%qes=47