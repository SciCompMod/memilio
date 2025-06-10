add_executable(output_processing output_processing.cpp)
target_link_libraries(output_processing PRIVATE memilio)

# target_compile_options(output_processing PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
