function architecture(path)

    character(len=*)  :: path
    character*5:: architecture,type(2)=(/'unix ','linux'/)
    integer*1  ::  iarch
    integer*2  :: itest=1

    open(99,file=path//'/type',form='unformatted',status='unknown')
    write(99) itest
    rewind(99)
    read(99) iarch
    close(99)	
    architecture=type(iarch+1)
end function architecture
