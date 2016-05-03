      subroutine magfild(latgeo,longeo,zref,year,Bmag,dipangle,orient)

      real,intent(in) :: latgeo,longeo,zref,year
      real,intent(out):: bmag, dipangle,orient

      COMMON/GENER/ UMR,ERA,AQUAD,BQUAD
       logical flginit
       data flginit/.true./

       UMR=0.01745329252
       ERA=6371.2
       AQUAD=4.0680925e7
       BQUAD=4.04085884e7

       if (flginit) then
         print*,'call feldcof'
         CALL FELDCOF(YEAR,DIMO)
         flginit=.false.
       endif
        
        print*,'call feldg'
        CALL FELDG(latgeo,longeo,zref,BNORTH,BEAST,BDOWN,Bmag)
c       CALL SHELLG(latgeo,longeo,zref,DIMO,XL,ICODE,BAB1)
         dipangle=90.-ASIN(BDOWN/Bmag)/UMR
       orient=1.
       if (BNORTH.lt.0.) orient=-1.
c       DEC=ASIN(BEAST/SQRT(BEAST*BEAST+BNORTH*BNORTH))/UMR


      END subroutine magfild
