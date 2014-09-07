      subroutine lec_indices(ian,imois,ijour,tu,ap,f107)

      implicit none

      integer ian,imois,ijour,i
      integer ian_deb,iday_deb,num,km,kj,k,num_rec
      real tu,ap(7),f107(3)
      integer*2 idat(6),kdat(7,8)

      include 'CHEMIN.INC'

      open(59,file=data_path(1:lpath_data)
     &             //'dir.geomag/data_geom.bin',
     &    access='direct',recl=2*62,status='unknown')

      read(59,rec=1,err=999) idat,kdat
      ian_deb = idat(1)
      iday_deb = idat(3)

      num = num_rec(ian_deb,iday_deb,ian,imois,ijour)

      read(59,rec=num,err=999) idat,kdat
C
C   Les 6 premiers paramŠtres lus sont respectivement :
C	1- l'ann‚e (4 chiffres)
C	2- le mois * 100 + le jour du mois
C	3- le quantiŠme du jour dans l'ann‚e en cours
C	4- le F10.7 du jour
C	5- Le F10.7 de la veille
C	6- le F10.7 moyenn‚ sur 81 jours centr‚s sur le jour consid‚r‚


      km = idat(2)/100
      kj = idat(2) - km*100
      f107(1) = float(idat(4))/10.
      f107(2) = float(idat(5))/10.
      f107(3) = float(idat(6))/10.

      
C
C  dans chaque tranche horaire ( 1 a 8), 
C  on a 7 valeurs de ap qui sont :
C       1- la moyenne journaliere
C	2- la valeur a l'heure 
C	3- la valeur 3 h avant
C	4- la valeur 6 h avant
C	5- la valeur 9 h avant
C	6- la moyenne entre 12 et 33 heures avant
C	7- la moyenne entre 36 et 57 heures avant
C

      k = tu/3 + 1

      do i=1,7
        ap(i)=kdat(i,k)
      enddo

      close(59)

      return

999   continue

      print*,'il y a erreur de lecture de data_geom.bin: ',num
      print*,ian_deb,iday_deb,ian,imois,ijour,tu
      print*,data_path(1:lpath_data)//'dir.geomag/data_geom.bin'
      print*,idat,kdat
      close(59)
      f107(1)=150.
      f107(2)=150.
      f107(3)=150.
      do i=1,7
	ap(i)=4.
      enddo
      return

      end






	integer function num_rec(ian_deb,iday_deb,ian,imois,ijour)

	integer m1(12),m2(12)

	data m1/31,59,90,120,151,181,212,243,273,304,334,365/
	data m2/31,60,91,121,152,182,213,244,274,305,335,366/

	
	nj_an = 1 - iday_deb
	
	do i = ian_deb+1,ian
	   kd = m1(12)
	   if(mod(i-1,4).eq.0) kd = kd+1
	   nj_an = nj_an + kd
	enddo
		

	if (imois .eq. 1) then
	   num_rec = nj_an + ijour
	   return
	else if (imois .eq. 2) then
	   num_rec = nj_an + ijour + 31
	   return
	endif
	
	njour=mod(ian,4)
	if (njour .eq. 0) then
	   num_rec = nj_an + m2(imois-1) + ijour
	else
	   num_rec = nj_an + m1(imois-1) + ijour 

        endif

        num_rec=max(num_rec,1)

        return
        end
