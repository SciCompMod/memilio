set(JULIA_BUNDLE "${CMAKE_CURRENT_SOURCE_DIR}/julia_lib")

add_executable(hospital_model hospital_model.cpp)
target_include_directories(hospital_model PRIVATE "${JULIA_BUNDLE}/include")
target_link_directories(hospital_model PRIVATE "${JULIA_BUNDLE}/lib")
target_link_libraries(hospital_model PRIVATE memilio abm juliatestpkg)
target_compile_options(hospital_model PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

set_target_properties(hospital_model PROPERTIES
    BUILD_RPATH   "${JULIA_BUNDLE}/lib;${JULIA_BUNDLE}/lib/julia"
    INSTALL_RPATH "${JULIA_BUNDLE}/lib;${JULIA_BUNDLE}/lib/julia")
