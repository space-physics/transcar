      real function coskhi(lat,long,h,month,day,ioption)
       
       use comm, only: deg2rad, pi
        implicit none

        
c******  calcul de cos(ki) ****
c**  lat latitude geographic **   degrees decimal
c**  long longitude **  degrees decimal
c**  h heure tu ** decimal hours
c** day ** day...
c**    ...of the month if ioption = 1
c**    ...of the year if ioption = 2
c** month ** month integer
c
        integer,intent(in) :: month,day,ioption
        real,intent(in) :: lat, long,h

        real decli,pideg,dlat,hloc,temps,tloc,ki
        integer jj

        integer,parameter :: moi(12)=[0,31,59,90,120,151,181,212,243,
     &    273, 304, 334]

        common /decl/decli

        dlat = lat * deg2rad
        hloc = (h + mod(long+360.,360.)/15.)/24.
        if (ioption == 1) then
            jj = day + moi(month)
        else
            jj = day
        endif
        
        temps = ((jj - 172.625) + hloc)/365.25 * (2*pi)
        tloc = hloc  * (2*pi)
        decli = 23.45 * deg2rad * cos(temps)
        coskhi = sin (decli)*sin(dlat)
     &          - cos(decli)  *  cos(dlat) * cos(tloc)
        coskhi=min(coskhi,1.)

      end function coskhi
