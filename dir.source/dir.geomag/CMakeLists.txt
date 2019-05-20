add_library(geomag geo2mag.f lec_indices.f90 ap2kp.f magfild.f shellig.f)
set_target_properties(geomag PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})