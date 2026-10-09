set(MEMILIO_CXX_FLAGS_SBML_GENERATED "${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS};")
string(APPEND MEMILIO_CXX_FLAGS_SBML_GENERATED "-Wno-unused-parameter;-Wno-unused-variable;-Wno-unused-but-set-variable;")

if(MEMILIO_ENABLE_SBML)
    add_executable(sbml_to_memilio sbml_to_memilio.cpp)
    target_link_libraries(sbml_to_memilio PRIVATE sbml memilio)
    memilio_target_warnings(sbml_to_memilio)
endif()
