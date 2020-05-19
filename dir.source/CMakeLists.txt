add_executable(transconvec ${CMAKE_CURRENT_SOURCE_DIR}/transconvec_13.op.f
  ${CMAKE_CURRENT_SOURCE_DIR}/reacrates.f)
target_link_libraries(transconvec PRIVATE fluid numeric system proj
                       convec geomag transelec imm common)
target_include_directories(transconvec PRIVATE ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/dir.source/dir.include)

add_subdirectory(dir.include)

add_subdirectory(dir.fluide)

add_subdirectory(dir.numeric)

add_library(system dir.fluide/jour_mois.f dir.cine/misc.f)

add_subdirectory(dir.projection)

add_subdirectory(dir.convection)

add_subdirectory(dir.geomag)

add_subdirectory(dir.cine)

add_subdirectory(dir.imm)
