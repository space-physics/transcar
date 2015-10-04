
	SUBROUTINE SOLR(FE,FA,FINF)
C   Calculates array of photon fluxes in 39 bands from EUVAC model
C   of Phil Richards [1994], which uses the F74113 spectrum given in
C   Torr et al. [1979]. Richards multiplies fluxes from 50-150A by factor of
C   3 (bins 3 and 4) and multiplies fluxes from 150-250A (bins 5 and 6) by
C   factor of 2. FINF values multiplied by 1E13 in PHION to give units m-2s-1
C                                                     - MJB  Dec. 30, 1993.
C   Constrain coronal fluxes to be greater than 80% of the modified F74113 
C   reference flux as per EUVAC paper    -MJB  April 25, 1994.
	REAL FE, FA, P, CHK
  	REAL F74(39), AI(39), FINF(39), ICOR(14)
	INTEGER I,J
	DATA F74 /0.012, 0.050,
     1  1.200, 0.450, 4.800, 3.100, 0.460, 0.210, 1.679, 0.800, 6.900,
     2  0.965, 0.650, 0.314, 0.383, 0.290, 0.285, 0.452, 0.720, 1.270,
     3  0.357, 0.530, 1.590, 0.342, 0.230, 0.360, 0.141, 0.170, 0.260,
     4  0.702, 0.758, 1.625, 3.537, 3.000, 4.400, 1.475, 3.500, 2.100,
     5  2.467/
	DATA AI /5.9375E-2, 3.2184E-2,
     *           1.0017E-2, 7.1250E-3, 1.3375E-2, 1.9450E-2, 2.7750E-3,
     1           1.3768E-1, 2.6467E-2, 2.5000E-2, 3.3333E-3, 2.2450E-2,
     2           6.5917E-3, 3.6542E-2, 7.4083E-3, 7.4917E-3, 2.0225E-2,
     3           8.7583E-3, 3.2667E-3, 5.1583E-3, 3.6583E-3, 1.6175E-2,
     4           3.3250E-3, 1.1800E-2, 4.2667E-3, 3.0417E-3, 4.7500E-3,
     5           3.8500E-3, 1.2808E-2, 3.2750E-3, 4.7667E-3, 4.8167E-3,
     6           5.6750E-3, 4.9833E-3, 3.9417E-3, 4.4167E-3, 5.1833E-3,
     7           5.2833E-3, 4.3750E-3/
	DATA ICOR/3,5,6,8,9,10,12,13,14,15,16,17,22,29/
	P = (FE+FA)/2.
	  DO I=1,39
	  FINF(I)=F74(I)*(1. + AI(I)*(P-80.))
	  WRITE(86,"(1X,I3,F10.4)") I,FINF(I)
Constrain coronal fluxes to be greater than 80% of the modified F74113 
C  reference flux.
	  CHK = 0.8 * F74(I)
	    IF (FINF(I) .LT. CHK) THEN
		DO J= 1,14
		  IF(I .EQ. ICOR(J) )THEN
		  FINF(I) = CHK
	          WRITE(86,"(1X,I3,' Flux constrained to',F10.4)") 
     1                                   I,FINF(I)
		  END IF
	  	END DO
	    END IF
	  END DO
	RETURN
	END

