add_library(imm OBJECT IMM2.f90 IMM_routines.f90)
set_target_properties(imm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})