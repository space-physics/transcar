c
c--------------------------- prot ---------------------------------
c
        subroutine prot(centE,Ebot,ddeng,altkm,
     &        iprt,nango2,nalt,jpreci,mcount,
     &        prodprotelec,primprotelec,fluxprimprot,qprimprot,
     &        prodionprot,gmu,gwt,densig,densneut, fic_transout)

        include 'comm.f'
        implicit none
        include 'TRANSPORT.INC'
    
c
c     Lus dans input :
c       production des proto-ions
c         proprotion_init(E,sp,z)        : cm-3.s-1.eV-1
c       production des proto-electrons
c        proelec_init(z)            : cm-3.s-1
c       production de protoelectrons (somme des proto-ions)
c         proprotel_init (z,E)        : cm-3.s-1.eV-1
c
c     Outputs
c
c         prodprotelec(z)         : cm-3.s-1
c         primprotelec(z,E)        : cm-3.s-1.eV-1
c         fluxprimprot(z,E)        : cm-2.s-1.eV-1.sr-1
c         qprimprot(E,z,A)        : cm-2.s-1.eV-1.sr-1
c     Outputs pour le moment gardes en interne (en attendant que ca
c     serve) : ion production rate (m-3 s-1)
c    prodionprot (isp,ialt)                  : cm-3.s-1
c     1,2,3,4,5,6 -> N2+,O2+,O+,H+,He+,N+
c
!     INPUTS
        integer, intent(in) :: fic_transout
        integer iprt(40),mcount(5),nspec
        real centE(nbren),Ebot(nbren),ddeng(nbren)
        real alt(nbralt),altkm(nbralt),densneut(8,nbralt)
        real gmu(2*nbrango2),gwt(2*nbrango2)
        integer nang,nalt,nango2
c
c     INTERNAL PARAMETERS
c
c     Parametres du flux de protons selon Hardy :
        real en_moy_keV,en_moy_eV,eflux,correc,eflux_init
        real Echar_protons
        real fluxdown(nbren,nbrango2),fluxup(nbren,nbrango2)
        real lonmag,latmag,tmag,cofo,cofh,cofn,Fe0,Ee0,Fi0,Ei0
        real Flux_part_int,Flux_ener_int,chi0,Eprotelec
c    prodionprot_init (isp,ialt)                  : cm-3.s-1
c     1,2,3,4 -> N2+,O2+,O+,N+
c     Isotropie du flux
        integer isotro
        real Bmag,dipangle,Enord,Eest,vperpnord,vperpest,vhorizon,vpara
        real ddp,Jtop
        integer ikp,iener
        integer,parameter:: lit=0        ! On passe par le code simplifie
        common/buff/lonmag,latmag,tmag,ikp,cofo,cofh,cofn,chi0,
     &                  Fe0,Ee0,Fi0,Ei0,
     &                  Bmag,dipangle,Enord,Eest,
     &                  vperpnord,vperpest,vhorizon,vpara,ddp,Jtop
        real xnN2(nbralt),xnO2(nbralt),xnO1(nbralt)
        integer iz
        real xmin,xmax,ymin,ymax,tab(300)

c     Lus :
        integer nnen, nz,nspecprot
        real proprotion_init(nbren,2,200)
        real eprot(nbren),z(200)
c
c     Calcules
        real proprotel_init(200,nbren),prodionprot_init(4,nbralt)
        real zmin,zmax,eprotmin,eprotmax,alp
c
        real pi,r2pi,densig(nbralt,nbren)
        real finput(nbren),foutput(nbren),finter(nbren,200)
        real proelec_init(200),prodion(200,nbrsp*2)
        integer ialt,ien,isp,jpreci,iang
c
c     OUTPUTS
        real prodprotelec(nbralt),prodionprot(6,nbralt)
        real qti(nbralt),qia(4,nbralt)
        real primprotelec(nbralt,nbren)
        real fluxprimprot(nbralt,nbren)
        real qprimprot(nbren,nbralt,-nbrango2:nbrango2)
c
900      format (5f10.2)
910      format (5(1pe10.2))
920      format (20x,'Altitude:',f10.2,' km')
1000     format (/,'Input protoelectron production [cm-3.s-1.eV-1]')
!1010     format (/,'# of Energies and altitudes for primary ',/,
!     &        'protoelectron fluxes (nnen,nz):',2i10)
1040     format (/,' Prod. (proprotel_init) (cm-3.s-1.eV-1):')
1050     format (/,'# of Energies and altitudes for secondary ',/,
     &        'electon fluxes (nen,nalt):',2i10)
1060     format (/,' Primary protoelectron flux : (cm-2.s-1.eV-1) :')
1070     format (/,' no of e. and alt. for sec. ',
     &        'el. fluxes (nen,nalt):',2i10)
1080     format ('Primary interpolated prod. :(cm-3.s-1.eV-1):')
1090    format('Electron primary production due to protons impact',/,
     &      '(prodprotelec) [cm-3.s-1]')
c
        do ialt = 1,nalt
          prodprotelec(ialt) = 0.
             do ien = 1,nen
              primprotelec(ialt,ien) = 0.
              fluxprimprot(ialt,ien) = 0.
              do iang = -nango2,nango2
!                            qprimprot(ien,ialt,iang) = 0.
              enddo
            enddo
        enddo
        if(jpreci.le.2) return
c
c     Si on est ici, c'est qu'on a des protons.
        pi=atan(1.)*4.
        r2pi = 1./(pi*2.)
c
c     Read protons from a file, lit = 1
c     calculate via code low_proton, lit = 0
c
c
        if (lit.eq.0)then
c
            iz = 0
            do ialt = 1,nalt
              write(stdout,*),altkm(ialt)
              if(altkm(ialt).le.840.)then
                iz = iz+1
                z(iz) = altkm(ialt)
                xnN2(iz) = densneut(1,ialt)
                xnO2(iz) = densneut(2,ialt)
                xnO1(iz) = densneut(3,ialt)
              endif
            enddo
            nz = iz
c
c       Calcul de l'energie caracteristique des ions
c
          call hardion(latmag,tmag,ikp,Flux_part_int,Flux_ener_int)
c
c         Les sorties de hardion sont en log decimal :
c         Flux de particules integre [cm-2.s-1.sr-1]
c
          Flux_part_int=10.**Flux_part_int
c
c         Flux d'energie integre [keV.cm-2.s-1.sr-1]
c
          Flux_ener_int=10.**Flux_ener_int
c
c         On assume que la forme du flux est une maxwellienne, bien que
c         hardy et al,(JGR Mai90, p4229-4248, 1985) ont publie des
c         formes de flux dont ils precisent autant que faire se peut
c         qu'elles ne sont justement pas maxwelliennes. Helas, ils ne
c         donnent pas une loi generale de representation de leurs
c         precipitations.
c         Pour une maxwellienne F = k E exp(-E/Eo), on a
c         Flux de particules = integrale de  k F(E) dE
c         Flux d'energie     = integrale de  k E F(E) dE
c         Le rapport des 2 donne l'energie moyenne
c             Eo = (Flux d'energie/2*Flux de particules)
c
          en_moy_keV = Flux_ener_int/(2.*Flux_part_int)    ! en keV
          en_moy_eV = en_moy_keV*1000.            ! en eV

          Fi0=6.2832*Flux_ener_int*1.6e-12
          Ei0= en_moy_eV
c
c      1 Erg.cm-2.s-1 = 1 mW.m-2
c       Pour avoir le flux d'energie en mW.m-2 a partir de
c       eV.cm-2.s-1.sr-1 il faut multiplier par
c       1e4(cm2/m2)*1.602 e-19(J/eV)*1 (W.s/J)*1e3(mW/W)
c       *4 pi (sr)* eflux (eV.cm-2.s-1.sr-1)
c       Le 4 pi vient de l'isotropie du flux dans Hardy. Est-ce
c       exact ?
c
       eflux_init = Flux_ener_int*20.131e-12     ! mW.m-2
c
c ----------------
cc       Valeurs test :
cc       Attention ! l'energie moyenne est 2 fois l'energie
cc       caracteristique de la maxwellienne.
c       eflux_init  = 1.                 ! mW.m-2 (=erg.cm-2.s-1)
c         Flux_ener_int=eflux_init/20.131e-12    ! keV.cm-2.s-1.sr-1
c       Echar_protons = 5.            ! keV
c       en_moy_keV = 2.* Echar_protons
c       en_moy_eV = 1000.* en_moy_keV
c ----------------
c
c       On veut des productions par electron volt. Le programme de
c       Marina les donne en global. On va donc calculer dans chaque
c       energie de la grille combien on a de particules,
c       puis faire comme si elles se repartissaient de facon
c       maxwellienne dans chaque element de la grille. On appellera
c       lowflux pour chacune de ces valeurs, et on divisera par
c       la largeur du pas de grille.
c
          isotro=2        ! tout va vers le bas
          call inmaxwl(Flux_ener_int,en_moy_eV,nango2,nen,centE,
     .                  isotro,gmu,fluxdown,fluxup)
          call normflux(Flux_ener_int,nango2,nen,centE,ddeng,gmu,gwt,
     .                  fluxdown,fluxup)
c
c       A partir de maintenant, on fait comme si chaque flux etait
c       repartit dans sa case suivant une maxwellienne
c
       do ien = nen,1,-1
c         Calcul de l'energie de l'electron ejecte. Je fais une
c         approximation honteuse : suivant la these de MG, quelque
c         soit la particule incidente de haute energie, l'energie
c         perdue est 1/alpha *sqrt(seuil*centE(ien)/1836)
c         Le seuil vaut 17,15,16 pour N2,O2,O, soit ... 16 (j'ai
c         honte).
c         alpha vaut 0.91
c         L'energie de l'electron est donc sqrt(centE(ien))/10.
         Eprotelec = sqrt(centE(ien))/10.
c
c           search in which box iener to put the created electron.
            call search (nen,Ebot,ddeng,Eprotelec,iener)

c         Il faut calculer le flux d'energie integree dans
c         chaque boite [eV.cm-2.s-1.sr-1]
         eflux = centE(ien)*ddeng(ien)*fluxdown(ien,1)

c         comme on est completement forward, il n'y a pas de
c         multiplication par 2pi, car on est directement en
c         eV.cm-2.s-1
            eflux = eflux*1.602e-12        ! mW.m-2
         alp = centE(ien)*1.e-3        ! keV

        write(stdout,*),'low_proton: nz=',nz
          call low_proton(eflux,alp,nz,z,xnN2,xnO2,xnO1,
     .        qti,qia)
c
c         Calcul des sorties de prot.f
         do ialt = 1,nz
           prodprotelec(ialt) = prodprotelec(ialt)+qti(ialt)
           primprotelec(ialt,iener)= primprotelec(ialt,iener)+
     .            qti(ialt)/ddeng(ien)
          enddo         ! sur les altitudes
       enddo            ! sur les energies
c       Passage de m-3 s-1 a cm-3 s-1
          do ialt = 1,nz
            prodprotelec(ialt) = prodprotelec(ialt)/1.e+06
            primprotelec(ialt,ien)= primprotelec(ialt,ien)/1.e+06
          enddo

c
c         On passe en flux.
       do ialt = 1,nz
         do ien = 1,nen
              fluxprimprot(ialt,ien)=
     .               primprotelec(ialt,ien)/densig(ialt,ien)
              do iang=-nango2,nango2
                qprimprot(ien,ialt,iang)=primprotelec(ialt,ien)/4./pi
              enddo
            enddo
       enddo
c
c       Cette facon de faire a sans doute change la production qu'
c       on aurait eue normalement. On va renormaliser
c
        call low_proton(eflux_init,en_moy_keV,nz,z,
     .        xnN2,xnO2,xnO1,qti,prodionprot_init)
c       Passage en cm-3.s-1
       do ialt = 1,nz
         qti(ialt) = qti(ialt)/1.e+6
         prodionprot(1,ialt) = prodionprot_init(1,ialt)/1.e+6  !N2+
         prodionprot(2,ialt) = prodionprot_init(2,ialt)/1.e+6  !O2+
         prodionprot(3,ialt) = prodionprot_init(3,ialt)/1.e+6  !O+
         prodionprot(4,ialt) = 0.              !H+ thermique
         prodionprot(5,ialt) = 0.                  !He+
         prodionprot(6,ialt) = prodionprot_init(4,ialt)/1.e+6  !N+
       enddo
c
       do ialt = 1,nz
         if(qti(ialt).eq.0. .or. prodprotelec(ialt).eq.0.) then
           correc = 0.
         else
           correc = qti(ialt)/prodprotelec(ialt)
         endif
         prodprotelec(ialt) = qti(ialt)
            do ien = 1,nen
           primprotelec(ialt,ien)= primprotelec(ialt,ien)*correc
              fluxprimprot(ialt,ien)= fluxprimprot(ialt,ien)*correc
              do iang=-nango2,nango2
                qprimprot(ien,ialt,iang)=qprimprot(ien,ialt,iang)*correc
              enddo
            enddo
       enddo
c
c       Il faut maintenant revenir aux altitudes de transelec
c
      else if(lit.eq.1) then
          open (58,file='dir.data/dir.linux/dir.cine/SOURCE_OUT_V',
     &             status='old')
          rewind(58)
c
c  *      read in primary protoelectron flux data
c         nnen   : # energies.
c         nz     : # altitudes.
c         eprot   : energies.
c         z      : altitudes.
c
       call xline(1,58)
          read  (58,*) nnen, nz,nspecprot
       if (nnen.gt.nbren)then
         write(6,*)'Tableau des productions de proto-electrons'
         write(6,*)'de taille en E > nbren. Declarer dans sub. prot'
         write(6,*)'un tableau de bonne taille et recommencer'
         stop
      endif
       if (nz.gt.200)then
         write(6,*)'Tableau des productions de proto-electrons'
         write(6,*)'de taille en z > 200. Declarer dans sub. prot'
         write(6,*)'un tableau de bonne taille et recommencer'
         stop
      endif

       call xline(1,58)
          read  (58,*) (eprot(ien), ien = 1,nnen)
       call xline(1,58)
          read  (58,*) (z(ialt), ialt = 1,nz)
c         production de H+ et H en cm-3s-1eV-1
       call xline(1,58)
          read  (58,*) (((proprotion_init(ien,isp,ialt),ien=1,nnen),
     .              isp=1,2),ialt=1,nz)
          do ien = 1,nnen
         do ialt = 1,nz
              proprotion_init(ien,1,ialt)=
     .            .75*proprotion_init(ien,1,ialt)
              proprotion_init(ien,2,ialt)=
     .            .75*proprotion_init(ien,2,ialt)
         enddo
       enddo
c         La somme vaut la production des protoelec. en cm-3s-1eV-1
          do ien = 1,nnen
            do ialt = 1,nz
              proprotel_init(ialt,ien)=
     .        proprotion_init(ien,1,ialt)+proprotion_init(ien,2,ialt)
            enddo
          enddo
c         Production des protoelectrons en cm-3.s-1  (prion)
       call xline(1,58)
          read  (58,*) (proelec_init(ialt), ialt = 1,nz)
       do ialt = 1,nz
            proelec_init(ialt)=.75*proelec_init(ialt)
       enddo
       close(58)
c
       call mnmx(z,nz,zmin,zmax,0)
       call mnmx(eprot,nnen,eprotmin,eprotmax,0)
c
            if(iprt(22).eq.1)then
              write (fic_transout,1000)
            write (fic_transout,*) 'Altitudes (z)'
            write (fic_transout,900) (z(ialt), ialt = 1,nz)
            write (fic_transout,*) 'Energies (eprot)'
            write (fic_transout,910) (eprot(ien), ien = 1,nnen)
            write (fic_transout,1040)
            do ialt=1,nz
              write(fic_transout,920)  z(ialt)
              write (fic_transout,910)(proprotel_init(ialt,ien),
     .                    ien = 1,nnen)
c             [proprotel_init]=cm-3.s-1.eV-1
         enddo
            endif
c
c       Interpolations en altitude
c       Production totale
       call intlin(nz,z,proelec_init,nalt,altkm,prodprotelec)
       do ialt = 1,nalt
         if(altkm(ialt).lt.zmin .or. altkm(ialt).gt.zmax)
     .        prodprotelec(ialt)=0.
       enddo
c
c       Interpolation en energie:
      do ialt=1,nz
        do ien=1,nnen
          finput(ien)=proprotel_init(ialt,ien)
         enddo
         call intlin(nnen,eprot,finput,nen,centE,foutput)
        do ien=1,nen
          finter(ien,ialt)=foutput(ien)
         enddo
       enddo
c
c       Interpolation en altitude:
       do ien=1,nen
         if(centE(ien).gt.eprotmin .and. centE(ien).lt.eprotmax) then
          do ialt=1,nz
            finput(ialt)=finter(ien,ialt)
           enddo
           call intlin(nz,z,finput,nalt,altkm,foutput)
          do ialt=1,nalt
             if(altkm(ialt).gt.zmin .and. altkm(ialt).lt.zmax)
     .            primprotelec(ialt,ien)=foutput(ialt)
           enddo
         endif
       enddo
c
          if(iprt(22).eq.1)then
c           print primary production rate.
            write (fic_transout,1070) nen, nalt
            write (fic_transout,*) 'Altitudes (alt)'
            write (fic_transout,900) (altkm(ialt), ialt = 1,nalt)
            write (fic_transout,*) 'Energies (e)'
            write (fic_transout,910) (centE(ien), ien = 1,nen)
            write (fic_transout,1080)
            do ialt=1,nalt,mcount(5)
              write(fic_transout,920)  altkm(ialt)
              write (fic_transout,910)(primprotelec(ialt,ien),ien=1,nen)
         enddo
          endif
c
c       On passe en flux.
       do ien = 1,nen
        do ialt=1,nalt
            fluxprimprot(ialt,ien)=
     .        primprotelec(ialt,ien)/densig(ialt,ien)
              do iang=-nango2,nango2
c             La division par 2 pi vient de l'integration polaire
                qprimprot(ien,ialt,iang)= fluxprimprot(ialt,ien)* r2pi
c             On redivise par 2 parce que le flux primaire est
c             isotrope (la moitie vers le haut et la moitie vers le
c             bas)
                qprimprot(ien,ialt,iang)= qprimprot(ien,ialt,iang)/2.
c             On redivise par 2 pour avoir 0,5 erg et pas 1
                qprimprot(ien,ialt,iang)= qprimprot(ien,ialt,iang)/2.
           enddo
         enddo
       enddo
c
          if(iprt(21).eq.1 .or. iprt(22).eq.1)then
c           print normalized source function.
            write (fic_transout,1050) nen, nalt
            write (fic_transout,*) 'Energies (e)'
            write (fic_transout,910) (centE(ien), ien = 1,nen)
            write (fic_transout,*) 'Altitudes (alt)'
            write (fic_transout,900) (altkm(ialt), ialt = 1,nalt)
            write (fic_transout,1090)
            write (fic_transout,900) (prodprotelec(ialt), ialt = 1,nalt)
            write (fic_transout,1060)
            do ialt=1,nalt,mcount(5)
              write(fic_transout,920)  altkm(ialt)
              write(fic_transout,910)(primprotelec(ialt,ien),ien=1,nen)
         enddo
          endif
c
      endif             ! sur calcul ou lecture

        end subroutine prot
