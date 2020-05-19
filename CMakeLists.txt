cmake_minimum_required(VERSION 3.13)
# allow default to debug build due to non-standard Fortran used
# will segfault with -O3 and GCC 7-9 at least
# -O2 is OK for GCC.

project(transcar
  LANGUAGES Fortran)

# compiler-specific flags necessary for Transcar legacy disk file format
if(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
  add_compile_options(-mtune=native -frecord-marker=4 -std=legacy
          -O2 -Wline-truncation -fmax-errors=3)# -Warray-bounds=2)
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES Intel)
  if(WIN32)
    add_compile_options(/O2 /assume:byterecl)
  else()
    add_compile_options(-O2 -assume byterecl)
  endif()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(dir.source)
