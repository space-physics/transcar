pure subroutine reacrates(Tn,TnHot,TOx,TH,TOHot,kOH,kHO,kOHotH,kHOHot,khothot)
implicit none
!	This subroutine computes the reaction rates for in units of cm^3/s
!
!	1)	O+ + H --> H+ + O
!	2)	H+ + O --> O+ + H
!
!	as indicated by Stancil, et al, AAS, 1999
real, intent(in) :: Tn, TnHot, Tox, TH,TOHot
real, intent(out):: kOH, kHO,kOHotH, KHOHot,khothot

real Tr

!	Reaction 1) above
real,parameter :: b1(2)=[2.08e-9, 1.11e-11]
real,parameter :: c1(2)=[4.05e-1, -4.58e-1]

!	Reaction 2) above
real,parameter :: b2(2)=[1.26e-9,4.25e-10]
real,parameter :: c2(2)=[5.17e-1,6.69e-3]

!	Cold reaction rates
Tr=(Tn +16.*TOx)/17.
kOH=b1(1)*(Tr/1.e4)**c1(1) + b1(2)*(Tr/1.e4)**c1(2)
Tr=(Tn*16.+TH)/17.
kHO=(b2(1)*(Tr/1.e4)**c2(1) + b2(2)*(Tr/1.e4)**c2(2))* exp(-227./Tr)
!	Hot reaction rates
Tr=(Tn +16.*TOHot)/17.
kOHotH=b1(1)*(Tr/1.e4)**c1(1) + b1(2)*(Tr/1.e4)**c1(2)
Tr=(TnHot*16.+TH)/17.
kHOHot=(b2(1)*(Tr/1.e4)**c2(1) + b2(2)*(Tr/1.e4)**c2(2))*exp(-227./Tr)

!	Hot-Hot reaction rates
Tr=(TOHot+TnHot)/2.
khothot=b1(1)*(Tr/1.e4)**c1(1) + b1(2)*(Tr/1.e4)**c1(2)

end subroutine reacrates
