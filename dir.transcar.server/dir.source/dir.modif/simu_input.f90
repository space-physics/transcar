subroutine simu_input(fid,filename,nrec,longrec,longbuf,buffer)

implicit none
integer		:: i,fid,nrec,longrec,longbuf
real		:: buffer(10000)
character*128	:: filename

open(fid,file=filename,form='unformatted',access='direct',recl=longrec,status='unknown')
read(fid,rec=nrec)(buffer(i),i=1,longbuf)
close(fid)
end subroutine simu_input

subroutine simu_output(fid,filename,nrec,longrec,longbuf,buffer)

implicit none
integer		:: i,fid,nrec,longrec,longbuf
real		:: buffer(10000)
character*128	:: filename

open(fid,file=filename,form='unformatted',access='direct',recl=longrec,status='unknown')
write(fid,rec=nrec)(buffer(i),i=1,longbuf)
close(fid)

end subroutine simu_output
