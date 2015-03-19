	program euvac

c programme calculant le flux selon Richards et al,JGR,99,8981,1994

	include 'CHEMIN.INC'

	dimension flux(39)

	open(10,file=data_path(1:lpath_data)//'EUVAC.dat')
	open(11,file=data_path(1:lpath_data)//'EUVACOUT')
	
c donner f et fbar
	f=60.
	fbar=62.
	p=f+fbar
	p=p/2.


	do i=1,39

	read(10,*)n,wave1,wave2,fbase,a	
	flux(i)=fbase*(1.+a*(p-80))*1e9
	write(11,*)n,wave1,wave2,flux(i)
	enddo

	return
	end
