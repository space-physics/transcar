	function coskhi(lat,long,h,mois,jour,ioption)
c
c******  calcul de cos(ki) ****
c**  lat latitude geographique **   degres decimaux
c**  long longitude **  degres decimaux
c**  h heure tu ** heures decimales
c** jour ** entier
c**	quantieme du mois si ioption = 1
c**	quantieme de l'annee si ioption = 2
c** mois ** entier
c
	integer moi(12)
	real lat, long,h,ki
	data pi,deuxpi/3.141592654,6.283185308/
	data moi/0,31,59,90,120,151,181,212,243,273,304,334/
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
	return
	end
