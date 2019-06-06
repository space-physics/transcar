cmake_minimum_required(VERSION 3.1)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(transcar Fortran)

# compiler dependent flags necessary for disk files used (sorry)
if(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
  add_compile_options(-frecord-marker=4 -std=legacy
          -march=native -Wline-truncation -fmax-errors=3)# -Warray-bounds=2)
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES Intel)
  if(WIN32)
    add_compile_options(/assume:byterecl)
  else()
    add_compile_options(-assume byterecl)
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES PGI)
    #already assumes byterecl
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES Flang)
    #already assumes byterecl
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(transconvec ${CMAKE_CURRENT_SOURCE_DIR}/transconvec_13.op.f
  ${CMAKE_CURRENT_SOURCE_DIR}/reacrates.f)
target_link_libraries(transconvec PRIVATE fluid numeric system proj
                       convec geomag transelec imm common)
target_include_directories(transconvec PRIVATE ${CMAKE_SOURCE_DIR}/dir.include)

add_subdirectory(dir.include)

add_subdirectory(dir.fluide)

add_subdirectory(dir.numeric)

add_library(system dir.fluide/jour_mois.f dir.cine/misc.f)

add_subdirectory(dir.projection)

add_subdirectory(dir.convection)

add_subdirectory(dir.geomag)

add_subdirectory(dir.cine)

add_subdirectory(dir.imm)

