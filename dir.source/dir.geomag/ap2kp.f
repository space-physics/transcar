	real function ap2kp(ap)

	real ap,xap(28)
	data xap/0.,2.,3.,4.,5.,6.,7.,9.,12.,15.,18.,22.,27.,32.,
     &		39.,48.,56.,67.,80.,94.,111.,132.,154.,179.,207.,
     &		236.,300.,400./

	i=1
	do while (ap.gt.xap(i))
	  i=i+1
	enddo
	i=max(i-1,1)
	dkp=(ap-xap(i))/(xap(i+1)-xap(i))
	ap2kp=(i-1+dkp)/3.

	return
	end
