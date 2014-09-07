	subroutine reacrates(Tn,TnHot,TOx,TH,TOHot,kOH,kHO,kOHotH,kHOHot
     &       ,khothot)

C	This subroutine computes the reaction rates for in units of cm^3/s
C
C	1)	O+ + H --> H+ + O
C	2)	H+ + O --> O+ + H
C
C	as indicated by Stancil, et al, AAS, 1999

	real TOx,TH,TOHot,kOH,kHO,KOHotH,KHOHot
	real Tr,b1(2),c1(2),b2(2),c2(2)

C	Reaction 1) above
	data b1 /2.08e-9,1.11e-11/
	data c1 /4.05e-1,-4.58e-1/

C	Reaction 2) above
	data b2 /1.26e-9,4.25e-10/
	data c2 /5.17e-1,6.69e-3/

C	Cold reaction rates
	Tr=(Tn +16.*TOx)/17.
	kOH=b1(1)*(Tr/1.e4)**c1(1) + b1(2)*(Tr/1.e4)**c1(2)
	Tr=(Tn*16.+TH)/17.
	kHO=(b2(1)*(Tr/1.e4)**c2(1) + b2(2)*(Tr/1.e4)**c2(2))*exp(-227./Tr)

C	Hot reaction rates
	Tr=(Tn +16.*TOHot)/17.
	kOHotH=b1(1)*(Tr/1.e4)**c1(1) + b1(2)*(Tr/1.e4)**c1(2)
	Tr=(TnHot*16.+TH)/17.
	kHOHot=(b2(1)*(Tr/1.e4)**c2(1) + b2(2)*(Tr/1.e4)**c2(2))*exp(-227./Tr)

C	Hot-Hot reaction rates
	Tr=(TOHot+TnHot)/2.
	khothot=b1(1)*(Tr/1.e4)**c1(1) + b1(2)*(Tr/1.e4)**c1(2)

	return
	end

