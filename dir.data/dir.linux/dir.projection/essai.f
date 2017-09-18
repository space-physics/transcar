        program essai
        implicit none	

        real buffer(100)
        integer irec,i

        open(10,file='varpot.dat',status='old',
     &		form='unformatted',access='direct',recl=652)

     	irec=1
1     	read(10,rec=irec) (buffer(i),i=1,163)
     	print*,(buffer(i),i=1,10)
        irec=irec+1
        goto 1

        end program essai
