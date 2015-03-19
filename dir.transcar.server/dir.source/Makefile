# first mkdir dir.transcar.server/dir.obj (one time)
# USAGE: 
# make -s fc=gfortran
# make -s fc=ifort
#
# Note that ifort will error out on msis90.f when using warn all

ifeq ($(strip $(fc)),)
fc=gfortran
endif

F77=$(fc)
F90=$(fc)

# -O3 gives SIGSEGV with gfortran, -O2 doesn't.
ifeq ($(fc),ifort)
OPTIONS=-assume byterecl  -O2 -mtune=native -march=native #-warn all #for ifort
else ifeq ($(fc),gfortran)
# we use std=legacy because there are deleted features used
OPTIONS=-frecord-marker=4 -O2 -mtune=native -march=native  -std=legacy #-pedantic #-Wall #for gfortran 
else
$(info you must specify your compiler with fc=mycompiler.  E.g. make fc=gfortran)
endif

COMPINFO=$(shell $(fc) --version | head -n1)
OPSYS=$(shell uname -a)
$(info $(COMPINFO) on $(OPSYS))


FLAG_COMPIL= $(OPTIONS) -I$(INCLUDE) -c
FLAG_LIEN= $(OPTIONS) -I$(INCLUDE) -o

LIB=$(LD_LIBRARY_PATH)/

#Definition des variables d'environnements pour la compilation
TCDIR := $(shell cd ..; pwd)

SOURCE=$(TCDIR)/dir.source
FLUIDE=$(TCDIR)/dir.source/dir.fluide
CINETIC=$(TCDIR)/dir.source/dir.cine
NUMERIC=$(TCDIR)/dir.source/dir.numeric
SYSTEM=$(TCDIR)/dir.source/dir.system
CONVECTION=$(TCDIR)/dir.source/dir.convection
PROJECTION=$(TCDIR)/dir.source/dir.projection
GEOMAG=$(TCDIR)/dir.source/dir.geomag
IMM=$(TCDIR)/dir.source/dir.imm
IRI=$(TCDIR)/dir.source/dir.cine/dir.iri
INCLUDE=$(TCDIR)/dir.source/dir.include
BUT_OBJ=$(TCDIR)/dir.obj
BUT_EXE=$(TCDIR)

#
all: $(BUT_EXE)/transconvec_13.op.out

clean: 
	$(RM) $(BUT_OBJ)/*.o
#
#Definition des modules
#----------------------
#

mod_fluide_op = $(BUT_OBJ)/atmos.op.o $(BUT_OBJ)/vent.o $(BUT_OBJ)/msis90.o \
             $(BUT_OBJ)/coskhi.o $(BUT_OBJ)/fchap.o

mod_numeric=$(BUT_OBJ)/lcpfct.o $(BUT_OBJ)/stabenerg.o $(BUT_OBJ)/pas_de_temps.o \
	    $(BUT_OBJ)/solve_ODE.o $(BUT_OBJ)/balanc.o $(BUT_OBJ)/balbak.o \
            $(BUT_OBJ)/elmhes.o $(BUT_OBJ)/eltran.o $(BUT_OBJ)/hqr2.o  $(BUT_OBJ)/cdiv.o

mod_systeme=$(BUT_OBJ)/architecture.o $(BUT_OBJ)/pwd.o $(BUT_OBJ)/isnan.o \
            $(BUT_OBJ)/jour_mois.o $(BUT_OBJ)/misc.o

mod_proj   =$(BUT_OBJ)/projection.o $(BUT_OBJ)/val_fit.o $(BUT_OBJ)/plegendr.o \
            $(BUT_OBJ)/coef_prec.o $(BUT_OBJ)/coef_pot.o $(BUT_OBJ)/coef_cour.o \
            $(BUT_OBJ)/precipitation.o $(BUT_OBJ)/potentiel.o \
            $(BUT_OBJ)/courant.o $(BUT_OBJ)/hardelec.o

mod_conv   =$(BUT_OBJ)/convec.o $(BUT_OBJ)/convec_1.o $(BUT_OBJ)/conv_ar.o

mod_geomag =$(BUT_OBJ)/geo2mag.o $(BUT_OBJ)/lec_indices.o $(BUT_OBJ)/ap2kp.o \
            $(BUT_OBJ)/magfild.o $(BUT_OBJ)/shellig.o

mod_cinetic_op=$(BUT_OBJ)/transelec.op.o $(BUT_OBJ)/lect.o  $(BUT_OBJ)/iniflu.o \
            $(BUT_OBJ)/degrad.o $(BUT_OBJ)/felin.op.o $(BUT_OBJ)/solflux.o \
            $(BUT_OBJ)/trans.op.o $(BUT_OBJ)/disort.o $(BUT_OBJ)/cineout.op.o \
            $(BUT_OBJ)/linsub.o  \
            $(BUT_OBJ)/quelle_grille.o $(BUT_OBJ)/quel_flux.o \
            $(BUT_OBJ)/flux_val.o $(BUT_OBJ)/hardion.o \
            $(BUT_OBJ)/iris12.o $(BUT_OBJ)/irif12.o \
            $(BUT_OBJ)/cira86_pour_iri.o $(BUT_OBJ)/ecr.o \
            $(BUT_OBJ)/conduc.o $(BUT_OBJ)/prec_time.o \
            $(BUT_OBJ)/dE_dt.o $(BUT_OBJ)/low_proton.o $(BUT_OBJ)/prot.o \

#	$(BUT_OBJ)/transplot.o $(BUT_OBJ)/felplot.o $(BUT_OBJ)/degplot.o $(BUT_OBJ)/atmplot.o

mod_imm    =$(BUT_OBJ)/IMM2.o $(BUT_OBJ)/IMM_routines.o



#
# Fabrication de transcar
# -----------------------
#
	
$(BUT_EXE)/transconvec_13.op.out:   $(BUT_OBJ)/transconvec_13.op.o \
	$(mod_fluide_op) $(mod_numeric) $(mod_systeme) \
	$(mod_proj) $(mod_conv) $(mod_geomag) $(mod_imm) \
	$(mod_cinetic_op) $(BUT_OBJ)/reacrates.o

	$(F90) $(FLAG_LIEN) $(BUT_EXE)/transconvec_13.op.out $(BUT_OBJ)/transconvec_13.op.o \
	$(mod_fluide_op) $(mod_numeric) $(mod_systeme) \
	$(mod_proj) $(mod_conv) $(mod_geomag) $(mod_imm) \
	$(mod_cinetic_op) $(BUT_OBJ)/reacrates.o -L$(LIB)
#

$(BUT_EXE)/tube0:   $(BUT_OBJ)/tube.o \
        $(BUT_OBJ)/atmos.o $(BUT_OBJ)/pas_de_temps.o \
	$(BUT_OBJ)/vent.o $(BUT_OBJ)/msis90.o $(BUT_OBJ)/coskhi.o \
        $(BUT_OBJ)/lcpfct.o $(BUT_OBJ)/stabenerg.o $(BUT_OBJ)/hardelec.o \
	$(BUT_OBJ)/hardion.o $(BUT_OBJ)/jour_mois.o \
        $(BUT_OBJ)/iad.o\
	$(BUT_OBJ)/transelec.o $(BUT_OBJ)/trans.o $(BUT_OBJ)/degrad.o \
        $(BUT_OBJ)/felin.o $(BUT_OBJ)/solflux.o $(BUT_OBJ)/disort.o \
        $(BUT_OBJ)/lect.o $(BUT_OBJ)/misc.o $(BUT_OBJ)/linsub.o \
	$(BUT_OBJ)/iniflu.o $(BUT_OBJ)/cineout.o $(BUT_OBJ)/dE_dt.o\
	$(BUT_OBJ)/iris12.o $(BUT_OBJ)/irif12.o $(BUT_OBJ)/cira86_pour_iri.o \
	$(BUT_OBJ)/quelle_grille.o $(BUT_OBJ)/quel_flux.o \
	$(BUT_OBJ)/ecr.o $(BUT_OBJ)/flux_val.o $(BUT_OBJ)/contenu.o\
	$(BUT_OBJ)/conduc.o $(BUT_OBJ)/fchap.o \
        $(BUT_OBJ)/convec.o $(BUT_OBJ)/convec_1.o $(BUT_OBJ)/geo2mag.o \
	$(BUT_OBJ)/ap2kp.o $(BUT_OBJ)/magfild.o $(BUT_OBJ)/plegendr.o \
	$(BUT_OBJ)/shellig.o $(BUT_OBJ)/lec_indices.o \
	$(BUT_OBJ)/conv_ar.o \
	$(BUT_OBJ)/low_proton.o  $(BUT_OBJ)/prot.o \
	$(BUT_OBJ)/pwd.o $(BUT_OBJ)/architecture.o $(BUT_OBJ)/isnan.o \
	$(BUT_OBJ)/courant.o $(BUT_OBJ)/potentiel.o  $(BUT_OBJ)/precipitation.o \
	$(BUT_OBJ)/coef_cour.o $(BUT_OBJ)/coef_pot.o  $(BUT_OBJ)/coef_prec.o \
	$(BUT_OBJ)/val_fit.o \
	$(BUT_OBJ)/projection.o	$(BUT_OBJ)/IMM_routines.o $(BUT_OBJ)/IMM2.o
#
	$(F90) $(FLAG_LIEN) $(BUT_EXE)/tube $(BUT_OBJ)/tube.o \
        $(BUT_OBJ)/atmos.o $(BUT_OBJ)/pas_de_temps.o \
	$(BUT_OBJ)/vent.o $(BUT_OBJ)/msis90.o $(BUT_OBJ)/coskhi.o \
        $(BUT_OBJ)/lcpfct.o $(BUT_OBJ)/stabenerg.o $(BUT_OBJ)/hardelec.o \
	$(BUT_OBJ)/hardion.o $(BUT_OBJ)/jour_mois.o \
        $(BUT_OBJ)/iad.o\
	$(BUT_OBJ)/transelec.o $(BUT_OBJ)/trans.o $(BUT_OBJ)/degrad.o \
        $(BUT_OBJ)/felin.o $(BUT_OBJ)/solflux.o $(BUT_OBJ)/disort.o \
        $(BUT_OBJ)/lect.o $(BUT_OBJ)/misc.o $(BUT_OBJ)/linsub.o \
	$(BUT_OBJ)/iniflu.o $(BUT_OBJ)/cineout.o $(BUT_OBJ)/dE_dt.o\
	$(BUT_OBJ)/iris12.o $(BUT_OBJ)/irif12.o $(BUT_OBJ)/cira86_pour_iri.o \
	$(BUT_OBJ)/quelle_grille.o $(BUT_OBJ)/quel_flux.o \
	$(BUT_OBJ)/ecr.o $(BUT_OBJ)/flux_val.o $(BUT_OBJ)/contenu.o\
	$(BUT_OBJ)/conduc.o $(BUT_OBJ)/fchap.o \
        $(BUT_OBJ)/convec.o $(BUT_OBJ)/convec_1.o $(BUT_OBJ)/geo2mag.o \
	$(BUT_OBJ)/ap2kp.o $(BUT_OBJ)/magfild.o $(BUT_OBJ)/plegendr.o \
	$(BUT_OBJ)/shellig.o $(BUT_OBJ)/lec_indices.o \
	$(BUT_OBJ)/conv_ar.o \
	$(BUT_OBJ)/low_proton.o  $(BUT_OBJ)/prot.o \
	$(BUT_OBJ)/pwd.o $(BUT_OBJ)/architecture.o $(BUT_OBJ)/isnan.o \
	$(BUT_OBJ)/courant.o $(BUT_OBJ)/potentiel.o  $(BUT_OBJ)/precipitation.o \
	$(BUT_OBJ)/coef_cour.o $(BUT_OBJ)/coef_pot.o  $(BUT_OBJ)/coef_prec.o \
	$(BUT_OBJ)/val_fit.o \
	$(BUT_OBJ)/projection.o	$(BUT_OBJ)/IMM_routines.o $(BUT_OBJ)/IMM2.o	
#
#
#		COMPILATIONS DES PROGRAMMES
#		---------------------------
#
# Programmes principaux
# -------------------
#

$(BUT_OBJ)/transconvec_13.op.o: $(SOURCE)/transconvec_13.op.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(SOURCE)/transconvec_13.op.f
#
$(BUT_OBJ)/reacrates.o: $(SOURCE)/reacrates.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(SOURCE)/reacrates.f
#
$(BUT_OBJ)/tube.o:	$(CONVECTION)/tube.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CONVECTION)/tube.f
#
$(BUT_OBJ)/extraction.o:	$(SOURCE)/extraction.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(SOURCE)/extraction.f


#
#
# Sources systemes
# ----------------
#
$(BUT_OBJ)/pwd.o:	$(SYSTEM)/pwd.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(SYSTEM)/pwd.f
$(BUT_OBJ)/isnan.o:	$(SYSTEM)/isnan.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(SYSTEM)/isnan.f90
$(BUT_OBJ)/architecture.o:	$(SYSTEM)/architecture.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(SYSTEM)/architecture.f90
#
#
# Sources numeriques
# ------------------
#
$(BUT_OBJ)/lcpfct.o:	$(NUMERIC)/lcpfct.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(NUMERIC)/lcpfct.f
$(BUT_OBJ)/stabenerg.o:	$(NUMERIC)/stabenerg.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(NUMERIC)/stabenerg.f90
$(BUT_OBJ)/pas_de_temps.o:	$(NUMERIC)/pas_de_temps.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(NUMERIC)/pas_de_temps.f
$(BUT_OBJ)/solve_ODE.o:		$(NUMERIC)/solve_ODE.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(NUMERIC)/solve_ODE.f90
$(BUT_OBJ)/balanc.o:		$(NUMERIC)/balanc.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(NUMERIC)/balanc.f
$(BUT_OBJ)/balbak.o:	$(NUMERIC)/balbak.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(NUMERIC)/balbak.f
$(BUT_OBJ)/elmhes.o:	$(NUMERIC)/elmhes.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(NUMERIC)/elmhes.f
$(BUT_OBJ)/eltran.o:	$(NUMERIC)/eltran.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(NUMERIC)/eltran.f
$(BUT_OBJ)/hqr2.o:	$(NUMERIC)/hqr2.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(NUMERIC)/hqr2.f
$(BUT_OBJ)/cdiv.o:	$(NUMERIC)/cdiv.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(NUMERIC)/cdiv.f
	
#
#
# Sources fluide
# --------------
#
$(BUT_OBJ)/atmos.op.o:	$(FLUIDE)/atmos.op.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(FLUIDE)/atmos.op.f
$(BUT_OBJ)/vent.o:	$(FLUIDE)/vent.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(FLUIDE)/vent.f
$(BUT_OBJ)/msis90.o:	$(FLUIDE)/msis90.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(FLUIDE)/msis90.f
$(BUT_OBJ)/coskhi.o:	$(FLUIDE)/coskhi.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(FLUIDE)/coskhi.f
$(BUT_OBJ)/jour_mois.o:	$(FLUIDE)/jour_mois.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(FLUIDE)/jour_mois.f
$(BUT_OBJ)/fchap.o:	$(FLUIDE)/fchap.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(FLUIDE)/fchap.f
$(BUT_OBJ)/courant_time.o:	$(FLUIDE)/courant_time.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(FLUIDE)/courant_time.f

#$(BUT_OBJ)/contenu.o:	$(FLUIDE)/contenu.f
#	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(FLUIDE)/contenu.f
#$(BUT_OBJ)/iad.o:	$(FLUIDE)/iad.f
#	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(FLUIDE)/iad.f


#
#
# Sources cinetique
# -----------------
#
# Transport : alt driver cinetique
$(BUT_OBJ)/transelec.op.o:	$(CINETIC)/transelec.op.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/transelec.op.f
# Lecture des entrees du pgme cinetique si appele independamment
$(BUT_OBJ)/lect.o:		$(CINETIC)/lect.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/lect.f
# Init. des param. du pgme cinetique si appele par le programme fluide.
$(BUT_OBJ)/iniflu.o:	$(CINETIC)/iniflu.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/iniflu.f
# Degradation d'energie
$(BUT_OBJ)/degrad.o:	$(CINETIC)/degrad.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/degrad.f
# Photoproduction primaire
$(BUT_OBJ)/felin.o: 	$(CINETIC)/felin.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/felin.f
#
$(BUT_OBJ)/felin.op.o: $(CINETIC)/felin.op.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/felin.op.f
# Photoproduction primaire : Calcul du flux de photons
$(BUT_OBJ)/solflux.o:	$(CINETIC)/solflux.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/solflux.f
# Transport cinetique : programme principal
$(BUT_OBJ)/trans.o:		$(CINETIC)/trans.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/trans.f
# Transport cinetique : alt programme principal
$(BUT_OBJ)/trans.op.o:		$(CINETIC)/trans.op.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/trans.op.f
# Transport : solver de l'equation differentielle
$(BUT_OBJ)/disort.o:	$(CINETIC)/disort.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/disort.f
# Sorties du transport cinetique pour le tranport fluide
$(BUT_OBJ)/cineout.o:	$(CINETIC)/cineout.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/cineout.f
# Sorties du transport cinetique pour le tranport fluide
$(BUT_OBJ)/cineout.op.o:	$(CINETIC)/cineout.op.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/cineout.op.f
# Divers sous programmes
$(BUT_OBJ)/misc.o:		$(CINETIC)/misc.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/misc.f
# Bibliotheque mathematique issue de linpack
$(BUT_OBJ)/linsub.o:	$(CINETIC)/linsub.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/linsub.f
# Dessins sous greg du transport
#$(BUT_OBJ)/transplot.o:	$(CINETIC)/transplot.f
#	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/transplot.f
# Dessins sous greg de la production primaire
#$(BUT_OBJ)/felplot.o:	$(CINETIC)/felplot.f
#	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/felplot.f
# Dessins sous greg du transport
#$(BUT_OBJ)/degplot.o:	$(CINETIC)/degplot.f
#	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/degplot.f
# Dessins sous greg de l'atmosphere
#$(BUT_OBJ)/atmplot.o:	$(CINETIC)/atmplot.f
#	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/atmplot.f
# choix des grilles d'energie et d'angle
$(BUT_OBJ)/quelle_grille.o:	$(CINETIC)/quelle_grille.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/quelle_grille.f
# choix du flux precipite
$(BUT_OBJ)/quel_flux.o:	$(CINETIC)/quel_flux.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/quel_flux.f
# calcul du flux
$(BUT_OBJ)/flux_val.o:	$(CINETIC)/flux_val.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/flux_val.f
# Modele statistique de precipitation des electrons de Hardy.
$(BUT_OBJ)/hardelec.o:	$(CINETIC)/hardelec.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/hardelec.f
# Modele statistique de precipitation des ions de Hardy.
$(BUT_OBJ)/hardion.o:	$(CINETIC)/hardion.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/hardion.f
# Estimation de Ne pour iniflu par iri
$(BUT_OBJ)/iris12.o:	$(IRI)/iris12.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(IRI)/iris12.f
$(BUT_OBJ)/irif12.o:	$(IRI)/irif12.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(IRI)/irif12.f
$(BUT_OBJ)/cira86_pour_iri.o:	$(IRI)/cira86_pour_iri.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(IRI)/cira86_pour_iri.f
# Ecriture eventuelle dans ELEC et NEUTRAL
$(BUT_OBJ)/ecr.o:	$(CINETIC)/ecr.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/ecr.f
# calcul eventuel des conductivites
$(BUT_OBJ)/conduc.o:	$(CINETIC)/conduc.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/conduc.f
#
$(BUT_OBJ)/prec_time.o:	$(CINETIC)/prec_time.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/prec_time.f
#
$(BUT_OBJ)/dE_dt.o:	$(CINETIC)/dE_dt.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/dE_dt.f
#
$(BUT_OBJ)/low_proton.o:	$(CINETIC)/low_proton.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/low_proton.f
#
$(BUT_OBJ)/prot.o:	$(CINETIC)/prot.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CINETIC)/prot.f
#
#
# Sources convection
# ------------------
#
#
$(BUT_OBJ)/convec.o:		$(CONVECTION)/convec.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CONVECTION)/convec.f
#
$(BUT_OBJ)/conv_ar.o:		$(CONVECTION)/conv_ar.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CONVECTION)/conv_ar.f
#
$(BUT_OBJ)/convec_1.o:		$(CONVECTION)/convec_1.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(CONVECTION)/convec_1.f
#
#
# Sources geomag
# ------------------
#
#
$(BUT_OBJ)/geo2mag.o:		$(GEOMAG)/geo2mag.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(GEOMAG)/geo2mag.f
#
$(BUT_OBJ)/lec_indices.o:		$(GEOMAG)/lec_indices.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(GEOMAG)/lec_indices.f90
#
$(BUT_OBJ)/ap2kp.o:		$(GEOMAG)/ap2kp.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(GEOMAG)/ap2kp.f
#
$(BUT_OBJ)/magfild.o:		$(GEOMAG)/magfild.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(GEOMAG)/magfild.f
#
$(BUT_OBJ)/shellig.o:		$(GEOMAG)/shellig.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(GEOMAG)/shellig.f
#

#
#
# Sources projection
# ------------------
#
#
$(BUT_OBJ)/courant.o:		$(PROJECTION)/courant.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(PROJECTION)/courant.f90
#
$(BUT_OBJ)/potentiel.o:		$(PROJECTION)/potentiel.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(PROJECTION)/potentiel.f90
#
$(BUT_OBJ)/precipitation.o:		$(PROJECTION)/precipitation.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(PROJECTION)/precipitation.f90
#
$(BUT_OBJ)/coef_cour.o:		$(PROJECTION)/coef_cour.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(PROJECTION)/coef_cour.f
#
$(BUT_OBJ)/coef_pot.o:		$(PROJECTION)/coef_pot.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(PROJECTION)/coef_pot.f
#
$(BUT_OBJ)/coef_prec.o:		$(PROJECTION)/coef_prec.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(PROJECTION)/coef_prec.f
#
$(BUT_OBJ)/val_fit.o:		$(PROJECTION)/val_fit.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(PROJECTION)/val_fit.f90
#
$(BUT_OBJ)/plegendr.o:		$(PROJECTION)/plegendr.f
	cd $(BUT_OBJ); $(F77) $(FLAG_COMPIL) $(PROJECTION)/plegendr.f
#
$(BUT_OBJ)/projection.o:		$(PROJECTION)/projection.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(PROJECTION)/projection.f90
#
#
#
# Sources IMM
# ------------------
#
#
$(BUT_OBJ)/IMM2.o:		$(IMM)/IMM2.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(IMM)/IMM2.f90
#
$(BUT_OBJ)/IMM_routines.o:		$(IMM)/IMM_routines.f90
	cd $(BUT_OBJ); $(F90) $(FLAG_COMPIL) $(IMM)/IMM_routines.f90
