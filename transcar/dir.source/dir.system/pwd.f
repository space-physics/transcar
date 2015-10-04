	  subroutine pwd(file_dir,len_file,data_dir,len_data)

	  integer,intent(out):: len_file,len_data
	  character(len=*),intent(inout):: file_dir,data_dir
	  character(len=5) archi,architecture

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

	  end subroutine pwd
