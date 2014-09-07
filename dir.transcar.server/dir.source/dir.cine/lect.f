c
c-------------------------- reed ----------------------------
c
      	subroutine lect (nspec,knm,nen,nalt,zbot,ztop,hrloc,ut,day,year,
     .		jpreci,tempexo,f107,f107a,Apind,chi,chideg,glat,glong,
     .		albedo,altkm,alt,tneutre,densneut,colden,botE,centE,
     .		ddeng,nang,nango2,angzb,gmu,gwt,fluxdown,fluxup,
     .		denelc,temelc,temion,dipang,smgdpa)
c
 	implicit none
c
	include 'TRANSPORT.INC'
c
 	integer nspec,knmneutral,nalt,jpreci,modatmos,neutspe
 	integer isp,ialt,ien,i,iang
 	real bid
 	real zbot,ztop,hrloc,ut,year,tempexo,f107,f107a,Apind,day,
     . 	 	fctemp,fcdens(nbrsp),glat,glong,alt(nbralt),
     .		altkm(nbralt),tneutre(nbralt),densneut(8,nbralt),
     .     	colden(8,nbralt),albedo
c
 	integer nen
 	real botE(nbren),centE(nbren),ddeng(nbren)
 	integer knm,ne,nang,nango2
 	real ezero,trav(nbren)
 	real angzb(2*nbrango2),gmu(2*nbrango2),gwt(2*nbrango2)
 	real fluxdown(nbren,nbrango2),fluxup(nbren,nbrango2)
 	real denelc(nbralt),temelc(nbralt),dipang(nbralt),
     .	       zel(nbralt),smgdpa(nbralt),temion(nbralt),chi,chideg
 	real pi
c
 	pi = 4.*atan(1.)
c 
c-------------------- lecture NEUTRAL ------------------------
c
c   	read neutral atmosphere parameters.
c     	nalt: # altitude grid points.
c     	zbot: minimum altitude.
c     	ztop: maximum altitude.
c
        open (ineutr,file= data_path(1:lpath_data)
     &                           //'dir.cine/NEUTRAL')
 	rewind(ineutr)
      	call xline(2,ineutr)
      	read (ineutr,*) nspec,knmneutral
      	call xline(4,ineutr)
      	read (ineutr,*) nalt
      	read (ineutr,*) zbot, ztop
c
c     	hrloc: hour angle (0 to 24) (local)
c     	glat: latitude (-90 to +90).
c     	day: julian day (1 to 365).
c     	tempexo: exospheric temperature (use default if below 10 deg).
c     	f107a: solar 10.7-cm line output (average).
c     	f107: solar 10.7-cm line output (this day).
c 	chideg : solar zenith angle
c     	fcdens: density scaling factor for n2,o2,o
c
      	call xline(15,ineutr)
      	read (ineutr,*) hrloc,day,year,jpreci,ut
      	read (ineutr,*) tempexo, f107, f107a,Apind,chideg
      	read (ineutr,*) fctemp,(fcdens(isp),isp=1,nspec)
      	read (ineutr,*) glat,glong,modatmos,albedo
 	chi = chideg*pi/180.
c
c 	Read the altitudes to compute the intensity.
      	call xline(1,ineutr)
      	read(ineutr,*)nalt
      	read(ineutr,*)(altkm(ialt),ialt=1,nalt)
 	do ialt =1,nalt
 	  alt(ialt) = altkm(ialt)*1.e+05
 	enddo
c
      	call xline(1,ineutr)
      	read(ineutr,*)(tneutre(ialt),ialt=1,nalt)
      	read(ineutr,*)neutspe
      	do isp=1,nspec
          call xline(1,ineutr)
          read(ineutr,*)(densneut(isp,ialt),ialt=1,nalt)
          call xline(1,ineutr)
          read(ineutr,*)(colden(isp,ialt),ialt=1,nalt)
      	enddo
 	close (ineutr)
c
c----------------------------------------------------------------
c------------------- Lecture ELEC ---------------------
c
     	open (unit=ielec,file =data_path(1:lpath_data)
     &                               //'dir.cine/ELEC')
      	call xline(4,ielec)
      	read (ielec,*) knm,ezero
 	if (knmneutral.ne.knm)then
 	  write(6,*)'+++=== WARNING ===+++'
 	  write(6,*)'Identification number in ELEC    =',knm
 	  write(6,*)'Identification number in NEUTRAL =',knmneutral
 	endif
      	call xline(1,ielec)
      	read (ielec,*)nen
c 	Bottom energies
      	call xline(1,ielec)
      	read (ielec,*) (botE(nen+1-i),i=1,nen)
c 	Center energies
      	call xline(1,ielec)
      	read (ielec,*) (centE(nen+1-i),i=1,nen)
      	call xline(1,ielec)
      	read (ielec,*) (ddeng(nen+1-i),i=1,nen)
c
c 	Lecture des angles. Les angles descendants vont de 0 a 
c 	90 degres et les angles montant de 90 a 180 degres.
c	gwt : gaussian weigths (         1 ...    nango2 : up, 
c				nbrango2+1 ...  2*nango2 : down)
c	gmu : gaussian angles
c
c	The up and downward intensity, flux etc is defined in terms of
c	positive and negative mu. My definition is, however, independent
c	of the direction of the magnetic field and thus a bit arbitrary.
c	On the northern hemisphere the field is pointing down, and
c	thus a pitch angle of 0 degrees (mu=1) would be down and a pitch
c	angle of 180 degrees (mu=-1) would be pointing up.  However, I
c	define the pitch angle 90-180 (mu=[0,-1]) to be down, and the 
c	pitch angle of 0-90 (mu=[0,+1]) up.  This way the transport code
c	doesn't care about the hemisphere anymore.
c 	The ELEC file is built after the Northern hemisphere assumption,
c 	i.e. {0,90} degrees --> downward, and {90,180} degrees --> 
c 	upward. So that the angle must be inverted for DISORT.
c
      	call xline(1,ielec)
      	read (ielec,*)nang
 	nango2 = nang/2
	if(nango2.gt.nbrango2) then
	  print*,' Too many streams: nbrango2 = ', nbrango2,
     .		 '                   nango2   = ', nango2
	  stop 'error'		! call abort
	end if
      	read(ielec,*) (angzb(iang),iang=1,nang)
        call xline(1,ielec)
      	read(ielec,*) (gmu(iang),iang=1,nang)
        call xline(1,ielec)
      	read(ielec,*) (gwt(iang),iang=1,nang)
c
c 	fluxdown (1--->nango2)=flux down 
c 	fluxup   (1--->nango2)=flux up 
c
        call xline(2,ielec)
        do iang=1,nango2
	  read(ielec,*)bid
          read(ielec,*) (fluxdown(nen+1-ien,iang),ien=1,nen)
 	enddo
        call xline(2,ielec)
        do iang=1,nango2
	  read(ielec,*)bid
          read(ielec,*) (fluxup(nen+1-ien,iang),ien=1,nen)
 	enddo
c
c  	read in ambient electron parameters.
c
      	call xline(6,ielec)
      	read(ielec,*) ne
 	if (ne.ne.nalt)then
 	  write(6,*)'nbre d''altitudes ELEC    = ',ne
 	  write(6,*)'nbre d''altitudes NEUTRAL = ',nalt
 	  write(6,*)'Check and restart'
 	  stop
 	endif
      	read(ielec,*) (zel(ialt), ialt = 1,ne)
      	read(ielec,*) (denelc(ialt),ialt = 1,ne)
      	call xline(1,ielec)
      	read(ielec,*) (temelc(ialt), ialt = 1,ne)
c 	Ion temperature
      	call xline(1,ielec)
      	read(ielec,*) (temion(ialt), ialt = 1,ne)
c 	dumps dTe/dz
      	call xline(1,ielec)
      	read(ielec,*) (bid, ialt = 1,ne)
c 	dumps the O+ composition
      	call xline(1,ielec)
      	read(ielec,*) (bid, ialt = 1,ne)
c 	reads the magnetic dip angle
      	call xline(1,ielec)
      	read(ielec,*) (dipang(ialt), ialt = 1,ne)
c 	reads the sin of the magnetic dip angle
      	call xline(1,ielec)
      	read(ielec,*) (smgdpa(ialt), ialt = 1,ne)
 	close (ielec)
c
      end
