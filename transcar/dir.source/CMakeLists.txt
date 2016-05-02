cmake_minimum_required(VERSION 2.8.12)
project(transcar Fortran)

add_compile_options(-march=native)

# compiler dependent flags necessary for disk files used (sorry)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
    add_compile_options(-frecord-marker=4 -std=legacy)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
    add_compile_options(-assume byterecl)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
    #already assumes byterecl
else()
    message(WARNING "${CMAKE_Fortran_COMPILER_ID} unknown compiler type, disk options may be invalid")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../..)

include_directories(dir.include)

add_executable(transconvec_13.op.out ../transconvec_13.op.f ../reacrates.f)
target_link_libraries(transconvec_13.op.out fluide numeric system proj 
                       convec geomag transelec imm)

#---- build libraries ----------
set(c ../dir.fluide)
add_library(fluide ${c}/atmos.op.f ${c}/vent.f ${c}/msis90.f ${c}/coskhi.f ${c}/fchap.f)

set(c ../dir.numeric)
add_library(numeric ${c}/lcpfct.f ${c}/stabenerg.f90 ${c}/pas_de_temps.f 
        ${c}/solve_ODE.f90 ${c}/balanc.f ${c}/balbak.f ${c}/elmhes.f 
        ${c}/eltran.f ${c}/hqr2.f  ${c}/cdiv.f)

set(c ../dir.system)
add_library(system ${c}/architecture.f90 ${c}/isnan.f90 ../dir.fluide/jour_mois.f ../dir.cine/misc.f)

set(c ../dir.projection)
add_library(proj ${c}/projection.f90 ${c}/val_fit.f90 ${c}/plegendr.f 
            ${c}/coef_prec.f ${c}/coef_pot.f ${c}/coef_cour.f
            ${c}/precipitation.f90 ${c}/potentiel.f90
            ${c}/courant.f90 ../dir.cine/hardelec.f)

set(c ../dir.convection)
add_library(convec ${c}/convec.f ${c}/convec_1.f ${c}/conv_ar.f)
target_link_libraries(convec proj)

set(c ../dir.geomag)
add_library(geomag ${c}/geo2mag.f ${c}/lec_indices.f90 ${c}/ap2kp.f
            ${c}/magfild.f ${c}/shellig.f)

set(c ../dir.cine)
add_library(transelec ${c}/transelec.op.f ${c}/lect.f  ${c}/iniflu.f 
            ${c}/degrad.f ${c}/felin.op.f 
            ${c}/solflux.f ${c}/trans.op.f ${c}/disort.f ${c}/cineout.op.f
            ${c}/linsub.f ${c}/quelle_grille.f ${c}/quel_flux.f
            ${c}/flux_val.f ${c}/hardion.f 
            ${c}/dir.iri/iris12.f ${c}/dir.iri/irif12.f ${c}/dir.iri/cira86_pour_iri.f 
            ${c}/ecr.f ${c}/conduc.f ${c}/prec_time.f
            ${c}/dE_dt.f ${c}/low_proton.f ${c}/prot.f)
target_link_libraries(transelec system fluide)

add_library(imm ../dir.imm/IMM2.f90 ../dir.imm/IMM_routines.f90)
