      Subroutine plegendr(n,x,pn,dpn)

	implicit none

	integer n,m,i
	real*8 x,x2,x3,x4,x5,pn,dpn
	real*8 pn_1,pn_2,dpn_1,dpn_2

	x2=x*x
	x4=x2*x2

	if (n.eq.0) then
	  pn=1.
	  dpn=0.
	elseif (n.eq.1) then
	  pn = x
	  dpn = 1.
	elseif (n.eq.2) then
	  x2=x*x
	  pn=1.5*x2-.5
	  dpn=3.*x
	elseif (n.eq.3) then
	  x2=x*x
	  x3=x2*x
	  pn=2.5*x3-1.5*x
	  dpn=7.5*x2-1.5
	elseif (n.eq.4) then
	  x2=x*x
	  x3=x2*x
	  x4=x3*x
	  pn=4.375*x4-3.75*x2+.375
	  dpn=17.5*x3-7.5*x
	elseif (n.eq.5) then
	  x2=x*x
	  x3=x2*x
	  x4=x3*x
	  x5=x4*x
	  pn=7.875*x5 - 8.75*x3+1.875*x
	  dpn=39.375*x4-26.25*x2+1.875
	else
	  x2=x*x
	  x3=x2*x
	  x4=x3*x
	  x5=x4*x
	  pn_2=4.375*x4-3.75*x2+.375
	  dpn_2=17.5*x3-7.5*x
	  pn_1=7.875*x5 - 8.75*x3+1.875*x
	  dpn_1=39.375*x4-26.25*x2+1.875
	  m=5
	  do i=1,n-m
	    m=m+1
	    pn=(2.-1./m)*x*pn_1-(1.-1./m)*pn_2
	    dpn=(2.-1./m)*(pn_1+x*dpn_1)-(1.-1./m)*dpn_2
	    pn_2=pn_1
	    dpn_2=dpn_1
	    pn_1=pn
	    dpn_1=dpn
	  enddo
	endif
	return
	end
