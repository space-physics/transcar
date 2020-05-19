add_library(numeric OBJECT lcpfct.f stabenerg.f90 pas_de_temps.f
        solve_ODE.f90 balanc.f balbak.f elmhes.f eltran.f hqr2.f cdiv.f)
set_target_properties(numeric PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})