	subroutine pwd(file_dir,len_file,data_dir,len_data)

	integer len_file,len_data
	character*(*) file_dir,data_dir
	character*5 archi,architecture

	call getcwd(file_dir)
	len_file=lenc(file_dir)
	if (file_dir(len_file:len_file).ne.'/') then
          len_file=len_file+1
          file_dir(len_file:len_file)='/'
	endif
        data_dir=file_dir(1:len_file)//'dir.data/'
        len_data=lenc(data_dir)
        archi=architecture(data_dir(1:len_data))
        data_dir=data_dir(1:len_data)//'dir.'//archi
        len_data=lenc(data_dir)
	data_dir=data_dir(1:len_data)//'/'
	len_data=len_data+1
	
	return
	end

	subroutine getcwd(file_dir)

	integer*4 ilen,ierror
	character*(*) file_dir

	call pxfgetcwd(file_dir,ilen,ierror)
!	file_dir='/home/blelly/dir.spacegrid'
!	file_dir=''
	return
	end
