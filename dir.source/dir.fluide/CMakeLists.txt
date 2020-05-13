add_library(fluid OBJECT atmos.op.f vent.f msis90.f coskhi.f fchap.f)
set_target_properties(fluid PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_link_libraries(fluid PRIVATE common)
