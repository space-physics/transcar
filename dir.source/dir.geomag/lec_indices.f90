subroutine lec_indices(ian,imois,ijour,tu,ap,f107)
! French->English:
! an      year
! mois    month
! jour    day

implicit none
integer,intent(in) :: ian,imois,ijour
real,intent(in)  :: tu
real,intent(out) :: ap(7),f107(3)

integer   :: recsize=2*62,offset
integer i,ihe,u
integer   :: ian_deb,iday_deb,km,kj,k,num_rec

integer*2 :: idat(6),kdat(7,8)
integer*4 :: ftell,num,max_num
integer*8 :: ftelli8
integer*2 ::tempo(124)


print*,'record size=',recsize

open(newunit=u,file='dir.data/dir.linux/dir.geomag/data_geom.bin', &
        access='direct',recl=recsize,form='unformatted',status='old')

call fseek(u,0,2)
max_num=(ftell(u)-8)/recsize



read(u,rec=1) idat
ian_deb = idat(1)
iday_deb = idat(3)


num = num_rec(ian_deb,iday_deb,ian,imois,ijour)
ihe=tu

print*,'num,max_num=',num,max_num
if (num<=max_num) then

  read(u,rec=num) idat,kdat

!
!   Les 6 premiers parametres lus sont respectivement :
!	1- l'annee (4 chiffres)
!	2- le mois * 100 + le jour du mois
!	3- le quantieme du jour dans l'annee en cours
!	4- le F10.7 du jour
!	5- Le F10.7 de la veille
!	6- le F10.7 moyenne sur 81 jours centrs sur le jour considere


  km = idat(2)/100
  kj = idat(2) - km*100
  f107(1) = float(idat(4))/10.
  f107(2) = float(idat(5))/10.
  f107(3) = float(idat(6))/10.


!
!  dans chaque tranche horaire ( 1 a 8),
!  on a 7 valeurs de ap qui sont :
!       1- la moyenne journaliere
!	2- la valeur a l'heure
!	3- la valeur 3 h avant
!	4- la valeur 6 h avant
!	5- la valeur 9 h avant
!	6- la moyenne entre 12 et 33 heures avant
!	7- la moyenne entre 36 et 57 heures avant
!

  k = ihe/3 + 1

  ap=kdat(:,k)

else

  print*,'Error reading file data_geom.bin for record: ',num
  print*,ian,imois,ijour,ihe
  f107(1)=150.
  f107(2)=150.
  f107(3)=150.
  ap=4.

endif

 close(u)

print*,idat
print*,f107,ap


end



pure integer function num_rec(ian_deb,iday_deb,ian,imois,ijour)
    implicit none
    integer,intent(in):: ian_deb,iday_deb,ian,imois,ijour

    !day of year for non-leap and leap years month ends
    integer,parameter :: m1(12)=(/31,59,90,120,151,181,212,243,273,304,334,365/),	&
                         m2(12)=(/31,60,91,121,152,182,213,244,274,305,335,366/)

    integer nj_an,kd,i,njour

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

end function num_rec
