add_executable(transconvec transconvec_13.op.f reacrates.f)
target_link_libraries(transconvec PRIVATE fluid numeric system proj
                       convec geomag transelec imm common)
target_include_directories(transconvec PRIVATE ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/dir.source/dir.include)

add_subdirectory(dir.include)

add_library(fluid OBJECT dir.fluide/atmos.op.f dir.fluide/vent.f dir.fluide/msis90.f
dir.fluide/coskhi.f dir.fluide/fchap.f
)
target_link_libraries(fluid PRIVATE common)


add_library(numeric OBJECT dir.numeric/lcpfct.f dir.numeric/stabenerg.f90 dir.numeric/pas_de_temps.f
dir.numeric/solve_ODE.f90 dir.numeric/balanc.f dir.numeric/balbak.f dir.numeric/elmhes.f
dir.numeric/eltran.f dir.numeric/hqr2.f dir.numeric/cdiv.f
)

add_library(system dir.fluide/jour_mois.f dir.cine/misc.f)

add_subdirectory(dir.projection)

add_library(convec OBJECT dir.convection/convec.f dir.convection/convec_1.f dir.convection/conv_ar.f)
target_link_libraries(convec PRIVATE proj)


add_library(geomag OBJECT dir.geomag/geo2mag.f dir.geomag/lec_indices.f90 dir.geomag/ap2kp.f
dir.geomag/magfild.f dir.geomag/shellig.f
)

add_subdirectory(dir.cine)

add_library(imm OBJECT dir.imm/IMM2.f90 dir.imm/IMM_routines.f90)
