SUBROUTINE IMM_COUR(ndg,mdg,phicour,Lmin,Lmax,latequi)

	IMPLICIT NONE

        integer :: ndg,mdg,ndeg,mdeg
	real*8  :: Lmin,Lmax,latmin,latmax,latequi
	
	integer,parameter :: rangcoef=1000
	real*8 :: phicour(rangcoef)
	real*8 :: coefpot(rangcoef),coefflux(rangcoef),&
			coefte(rangcoef),coefalg(rangcoef)
				
	common/phi/coefpot,coefflux,coefte,coefalg,latmin,latmax,ndeg,mdeg

INTERFACE
    SUBROUTINE IMM_POT(ndg,mdg,phipot,Lmin,Lmax,latequi)

        integer :: ndeg,mdeg,ndg,mdg
	real*8  :: Lmin,Lmax,latmin,latmax,latequi
	
	integer,parameter :: rangcoef=1000
	real*8 :: phipot(rangcoef)
	real*8 :: coefpot(rangcoef),coefflux(rangcoef),&
			coefte(rangcoef),coefalg(rangcoef)
				
	common/phi/coefpot,coefflux,coefte,coefalg,latmin,latmax,ndeg,mdeg

     END SUBROUTINE IMM_POT


     SUBROUTINE IMM_PREC(ndg,mdg,phienerg,phifluxE,Lmin,Lmax,latequi)
        integer:: ndeg,mdeg,ndg,mdg
	real*8 :: Lmin,Lmax,latmin,latmax,latequi

	integer,parameter :: rangcoef=1000
	real*8 :: phifluxE(rangcoef),phienerg(rangcoef)
	
	real*8 :: coefpot(rangcoef),coefflux(rangcoef),&
			coefte(rangcoef),coefalg(rangcoef)
				
	common/phi/coefpot,coefflux,coefte,coefalg,latmin,latmax,ndeg,mdeg

     END SUBROUTINE IMM_PREC

END INTERFACE

Lmin=latmin
Lmax=latmax
latequi=Lmin-5.d0
ndg=ndeg
mdg=mdeg

phicour(1:(ndg+1)*(2*mdg+1))=coefalg(1:(ndg+1)*(2*mdg+1))

END SUBROUTINE IMM_COUR




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE IMM_POT(ndg,mdg,phipot,Lmin,Lmax,latequi)
	
	IMPLICIT NONE

        integer :: ndeg,mdeg,ndg,mdg
	real*8  :: Lmin,Lmax,latmin,latmax,latequi
	
	integer,parameter :: rangcoef=1000
	real*8 :: phipot(rangcoef)
	real*8 :: coefpot(rangcoef),coefflux(rangcoef),&
			coefte(rangcoef),coefalg(rangcoef)
				
common/phi/coefpot,coefflux,coefte,coefalg,latmin,latmax,ndeg,mdeg

Lmin=latmin
Lmax=latmax
latequi=Lmin-5.d0
ndg=ndeg
mdg=mdeg
phipot(1:(ndg+1)*(2*mdg+1))=coefpot(1:(ndg+1)*(2*mdg+1))


END SUBROUTINE IMM_POT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE IMM_PREC(ndg,mdg,phienerg,phifluxE,Lmin,Lmax,latequi)

	IMPLICIT NONE

        integer:: ndeg,mdeg,ndg,mdg
	real*8 :: Lmin,Lmax,latmin,latmax,latequi
	
	integer,parameter :: rangcoef=1000
	real*8 :: phifluxE(rangcoef),phienerg(rangcoef)
	
	real*8 :: coefpot(rangcoef),coefflux(rangcoef),&
			coefte(rangcoef),coefalg(rangcoef)
				
common/phi/coefpot,coefflux,coefte,coefalg,latmin,latmax,ndeg,mdeg

Lmin=latmin
Lmax=latmax
latequi=Lmin-5.d0
ndg=ndeg
mdg=mdeg
phienerg(1:(ndg+1)*(2*mdg+1))=coefte(1:(ndg+1)*(2*mdg+1))
phifluxE(1:(ndg+1)*(2*mdg+1))=coefflux(1:(ndg+1)*(2*mdg+1))

END SUBROUTINE IMM_PREC
