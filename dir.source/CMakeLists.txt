add_compile_options(${transcar_optim_flag} ${transcar_file_opts})

add_library(common OBJECT dir.include/comm.f90)

add_library(fluid OBJECT dir.fluide/atmos.op.f dir.fluide/vent.f dir.fluide/msis90.f
dir.fluide/coskhi.f dir.fluide/fchap.f
$<TARGET_OBJECTS:common>
)

add_library(numeric OBJECT dir.numeric/lcpfct.f dir.numeric/stabenerg.f90 dir.numeric/pas_de_temps.f
dir.numeric/solve_ODE.f90 dir.numeric/balanc.f dir.numeric/balbak.f dir.numeric/elmhes.f
dir.numeric/eltran.f dir.numeric/hqr2.f dir.numeric/cdiv.f
)

add_library(system OBJECT
dir.fluide/jour_mois.f
dir.cine/misc.f
)

add_library(proj OBJECT
dir.projection/projection.f90
dir.projection/val_fit.f90
dir.projection/plegendr.f
dir.projection/coef_prec.f
dir.projection/coef_pot.f
dir.projection/coef_cour.f
dir.projection/precipitation.f90
dir.projection/potentiel.f90
dir.projection/courant.f90
${transcar_SOURCE_DIR}/dir.source/dir.cine/hardelec.f
$<TARGET_OBJECTS:common>
)

add_library(convec OBJECT
dir.convection/convec.f
dir.convection/convec_1.f
dir.convection/conv_ar.f
$<TARGET_OBJECTS:proj>
)


add_library(geomag OBJECT
dir.geomag/geo2mag.f
dir.geomag/lec_indices.f90
dir.geomag/ap2kp.f
dir.geomag/magfild.f
dir.geomag/shellig.f
)

add_library(transelec OBJECT dir.cine/transelec.op.f dir.cine/lect.f dir.cine/iniflu.f
dir.cine/degrad.f dir.cine/felin.op.f
dir.cine/solflux.f dir.cine/trans.op.f dir.cine/disort.f dir.cine/cineout.op.f
dir.cine/linsub.f dir.cine/quelle_grille.f dir.cine/quel_flux.f
dir.cine/flux_val.f dir.cine/hardion.f
dir.cine/dir.iri/iris12.f dir.cine/dir.iri/irif12.f dir.cine/dir.iri/cira86_pour_iri.f
dir.cine/ecr.f dir.cine/prec_time.f dir.cine/conduc.f
dir.cine/dE_dt.f dir.cine/low_proton.f dir.cine/prot.f
$<TARGET_OBJECTS:common>
$<TARGET_OBJECTS:fluid>
$<TARGET_OBJECTS:system>
)
target_include_directories(transelec PRIVATE dir.include)


add_library(imm OBJECT
dir.imm/IMM2.f90
dir.imm/IMM_routines.f90
)
