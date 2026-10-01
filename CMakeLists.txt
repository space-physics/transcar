cmake_minimum_required(VERSION 3.20)
# due to non-standard Fortran used
# may segfault with -O3 with GCC 7-9 at least
# -O2 is OK for GCC.

project(transcar
  LANGUAGES Fortran)

# compiler-specific flags necessary for Transcar legacy disk file format
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options(
  -std=legacy
  -Wline-truncation
  -fmax-errors=3
  )
  # -Warray-bounds=2)

  if(APPLE AND NOT DEFINED transcar_optim_flag)
    set(transcar_optim_flag -O1)
  endif()

  set(transcar_file_opts -frecord-marker=4)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  set(transcar_file_opts -assume:byterecl)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${transcar_BINARY_DIR})

add_executable(transconvec dir.source/transconvec_13.op.f dir.source/reacrates.f
$<TARGET_OBJECTS:fluid>
$<TARGET_OBJECTS:numeric>
$<TARGET_OBJECTS:system>
$<TARGET_OBJECTS:proj>
$<TARGET_OBJECTS:convec>
$<TARGET_OBJECTS:geomag>
$<TARGET_OBJECTS:transelec>
$<TARGET_OBJECTS:imm>
$<TARGET_OBJECTS:common>
)
target_include_directories(transconvec PRIVATE ${transcar_BINARY_DIR} dir.source/dir.include)
target_compile_options(transconvec PRIVATE ${transcar_optim_flag} ${transcar_file_opts})

add_subdirectory(dir.source)
