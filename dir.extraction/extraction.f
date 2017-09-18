	program extraction

	parameter (npt=500,ncol0=50,nsiz=(npt+2)*ncol0)

	real buffer(nsiz)
	character*80 filein,fileou

	write(*,*)'Nom du fichier temporel ? '
	read(*,'(a)') filein
	write(*,*)'Nom du fichier de sortie ?'
	read(*,'(a)')fileou

	open(10,file=filein,form='unformatted',access='direct',
     &		status='unknown',recl=9)

	read(10,rec=1) (buffer(i),i=1,9)
	nx=buffer(1)
	ncol=buffer(2)
	ian=buffer(3)
	imois=buffer(4)
	ijour=buffer(5)
	iheure=buffer(6)
	imins=buffer(7)
	isec=buffer(8)
	deltat=buffer(9)

100	format('La date de debut du fichier est le : ',i2,'/',i2,'/',
     &  i2,' … ',i2,'h',i2,'m',i2,'s')
	write(*,100)ijour,imois,ian,iheure,imins,isecs

	close(10)

	nligne=nx+2
	longbuf=nligne*ncol
	longrec=longbuf
	
	open(10,file=filein,form='unformatted',access='direct',
     &		status='unknown',recl=longrec)

	open(20,file=fileou,form='unformatted',access='direct',
     &		status='unknown',recl=longrec)

 	write(*,*)'Quelle heure voulez-vous extraire (JJ HH MM SS) ? '
 	read(*,*)ijourf,iheuf,iminsf,isecf

 	xrec=float((86400*(ijourf-ijour)+3600*(iheuf-iheure)
     &		  +60*(iminsf-imin)+isecsf-isecs))/deltat
 	nrec=int(xrec)
 	if (xrec-nrec.ge..5) nrec=nrec+1
c	write(*,*)'numero du record ?'
c	read(*,*) nrec

200	format('La date trouvee est le : ',i2,'/',i2,'/',
     &  i2,' … ',i2,'h',i2,'m',i2,'s')

997	read(10,rec=nrec,err=999)(buffer(i),i=1,longbuf)	
	goto 998
999	backspace(10)
	inquire(10,nextrec=nrec)
	goto 997
998	continue
	nx=buffer(1)
	ncol=buffer(2)
	ian=buffer(3)
	imois=buffer(4)
	ijour=buffer(5)
	iheure=buffer(6)
	imins=buffer(7)
	isec=buffer(8)
	deltat=buffer(9)

	write(*,200)ijour,imois,ian,iheure,imins,isecs

	write(20,rec=1)(buffer(i),i=1,longbuf)

	close(10)
	close(20)
	end
