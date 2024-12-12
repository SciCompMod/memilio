add_library(dabm
    model.h
    model.cpp
    simulation.h
    simulation.cpp
    parameters.h
)
target_link_libraries(dabm PUBLIC memilio)
target_include_directories(dabm PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(dabm PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
