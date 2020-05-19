add_library(transelec transelec.op.f lect.f  iniflu.f
            degrad.f felin.op.f
            solflux.f trans.op.f disort.f cineout.op.f
            linsub.f quelle_grille.f quel_flux.f
            flux_val.f hardion.f
            dir.iri/iris12.f dir.iri/irif12.f dir.iri/cira86_pour_iri.f
            ecr.f prec_time.f conduc.f
            dE_dt.f low_proton.f prot.f)
target_link_libraries(transelec PRIVATE system fluid common)
set_target_properties(transelec PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
target_include_directories(transelec PRIVATE ${PROJECT_SOURCE_DIR}/dir.source/dir.include)