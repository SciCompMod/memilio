
if(MEMILIO_HAS_SBML)
    add_executable(sbml_to_memilio sbml_to_memilio.cpp)
    target_link_libraries(sbml_to_memilio PRIVATE sbml memilio)
    target_compile_options(sbml_to_memilio PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()

