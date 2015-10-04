        real function flux_val(energ,E0,FE)

        implicit none
        real energ,E0,FE
        real timeser(1024), edist(1024)
        real fluxdist(1024,1024)
        integer ntimeser,nfluxdist,precint,precext
        common /precdist/ nfluxdist,ntimeser,timeser,
     &          edist,fluxdist,precint,precext
        integer m,li
        real logslope,fm1,fm,em1,em,th
        logical precflag
        integer prect
        common/precstate/precflag,prect


        !Check whether to precipitate or not
        if(.not. precflag) then
          flux_val=0.
          return
        endif

!        print*,'Using flux dist.:  ',prect,'; size:  ',nfluxdist
!        write(*,*)(edist(m),fluxdist(m,prect),m=1,nfluxdist)
!        print*,'Interp vals:  ',precint,precext

        !Find energy index into flux array
        if(energ .lt. edist(1)) then
          m=0
        else
          m=2
          do while(energ .ge. edist(m) .and. m .le. nfluxdist)
            m=m+1
          enddo
        endif 

        !Compute value of flux
        if(m .eq. 0 .or. m .gt. nfluxdist) then         !perform extrapolation
          if(precext .eq. 0) then                       !extrapolate
            if(m .eq. 0) then
              fm1=fluxdist(1,prect)
              fm=fluxdist(2,prect)
              em1=edist(1)
              em=edist(2)
            else
              fm1=fluxdist(nfluxdist-1,prect)
              fm=fluxdist(nfluxdist-2,prect)
              em1=edist(nfluxdist-1)
              em=edist(nfluxdist-2)
            endif

            logslope=(log10(fm)-log10(fm1))/(log10(em)-log10(em1))
            flux_val=10**(logslope*(log10(energ)-log10(em))+log10(fm))

!-- for case study.  only extrapolate to low energies w/cutoff at 5eV
            if(m .ne. 0 .or. energ .lt. 5.0) then
              flux_val=0
            endif
!---MZ 
          elseif(precext .eq. 1) then                                   !truncate
            flux_val=0.
          else                                                          !hold
            if(m .eq. 0) then
              flux_val=fluxdist(1,prect)
            else
              flux_val=fluxdist(nfluxdist-1,prect)
            endif
          endif
        else                                            !perform interpolation
          if(precint .eq. 0) then                       !log-linear interpolation
            if(m .eq. nfluxdist) then
              li=nfluxdist-1
            else
              li=m
            endif
            fm1=fluxdist(li-1,prect)
            fm=fluxdist(li,prect)
            em1=edist(li-1)
            em=edist(li)
        
            logslope=(log10(fm)-log10(fm1))/(log10(em)-log10(em1))
            flux_val=10**(logslope*(log10(energ)-log10(em))+log10(fm))
          else                                          !sample and hold
            flux_val=fluxdist(m-1,prect)
          endif
        endif

        return
        end




        subroutine flux_integre(nango2,nen,centE,ddeng,
     .                  gmu,gwt,fluxdown,fluxup,E0,F0)
c
c       Normalise le flux d'entree a une valeur en energie Fe donnee.
c
        implicit none
        include 'TRANSPORT.INC'
c
        real centE(nbren),ddeng(nbren),gmu(nbrang),gwt(nbrang)
        real fluxup(nbren,nbrango2),fluxdown(nbren,nbrango2)
        integer nen,isotro,nango2,iang,ien
        real F0,E0,qsum,qtot,xnorm
c
c       Compute input energy in eV/cm2/sec/sr
        F0 = 0.               ! total flux energy in eV/cm2/sec/sr
        E0 = 0.               ! mean energy in eV
        do iang=1,nango2
          do ien=1,nen
            F0=F0+fluxdown(ien,iang)*gwt(iang)*
     .            gmu(iang)*centE(ien)*ddeng(ien)
            E0=E0+fluxdown(ien,iang)*gwt(iang)*
     .            gmu(iang)*ddeng(ien)
          enddo
        enddo
        E0=F0/E0
        return
        end

