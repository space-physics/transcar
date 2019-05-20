add_library(convec convec.f convec_1.f conv_ar.f)
target_link_libraries(convec PRIVATE proj)
set_target_properties(convec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})