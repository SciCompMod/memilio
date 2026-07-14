add_executable(hospital_model hospital_model)
target_link_libraries(hospital_model PRIVATE memilio abm)
target_compile_options(hospital_model PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
