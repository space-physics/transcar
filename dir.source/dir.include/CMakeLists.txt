add_library(common comm.f90)
set_target_properties(common PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
