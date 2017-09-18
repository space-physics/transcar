	real function iad(s)

	real s,a,b,c,d,e,f,alpha,beta

	data a,b,c,d,e,f/28.06,-53.63,32.28,-10.25,2.37,1.02/
	data alpha,beta/.31,1./

	if (s.gt.1..or.s.lt.0.) then
	  iad=0.
	else
	  iad=a*s+b
	  iad=iad*s+c
	  iad=iad*s+d
	  iad=iad*s+e
	  iad=iad*s+f
	  iad=s**alpha*(1.-s)**beta*exp(iad)
	endif
	return
	end
