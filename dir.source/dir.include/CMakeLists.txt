add_library(common OBJECT comm.f90)
set_target_properties(common PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
