cmake_minimum_required(VERSION 3.20)
# due to non-standard Fortran used
# may segfault with -O3 with GCC 7-9 at least
# -O2 is OK for GCC.

project(transcar
  LANGUAGES Fortran)

# compiler-specific flags necessary for Transcar legacy disk file format
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-frecord-marker=4
  -std=legacy
  -Wline-truncation
  -fmax-errors=3)
  # -Warray-bounds=2)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  add_compile_options(-assume:byterecl)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${transcar_BINARY_DIR})

add_executable(transconvec dir.source/transconvec_13.op.f dir.source/reacrates.f)
target_link_libraries(transconvec PRIVATE fluid numeric system proj convec geomag transelec imm common)
target_include_directories(transconvec PRIVATE ${transcar_BINARY_DIR} dir.source/dir.include
)

add_subdirectory(dir.source)
