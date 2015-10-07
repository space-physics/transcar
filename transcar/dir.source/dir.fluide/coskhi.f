        real function coskhi(lat,long,h,mois,jour,ioption)
        implicit none
c******  calcul de cos(ki) ****
c**  lat latitude geographic **   degrees decimal
c**  long longitude **  degrees decimal
c**  h heure tu ** decimal hours
c** jour ** day...
c**    ...of the month if ioption = 1
c**    ...of the year if ioption = 2
c** mois ** month integer
c
        integer,intent(in) :: mois,jour,ioption
        real,intent(in) :: lat, long,h

        real decli,pideg,dlat,hloc,temps,tloc,ki
        integer jj

        real,parameter :: pi=3.141592654,deuxpi=6.283185308
        integer,parameter :: moi(12)=(/0,31,59,90,120,151,181,212,243,
     &    273, 304, 334/)

        common /decl/decli

        pideg = pi/180.
        dlat = lat * pideg
        hloc = (h + mod(long+360.,360.)/15.)/24.
        if (ioption.eq.1) then
            jj = jour + moi(mois)
        else
            jj = jour
        endif
        temps = ((jj - 172.625) + hloc)/365.25 * deuxpi
        tloc = hloc  * deuxpi
        decli = 23.45 * pideg * cos(temps)
        coskhi = sin (decli)*sin(dlat)
     &          - cos(decli)  *  cos(dlat) * cos(tloc)
        coskhi=min(coskhi,1.)

        end function coskhi
