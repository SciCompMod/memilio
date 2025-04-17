add_library(smm
    parameters.h
    model.h
    model.cpp
    simulation.h
    simulation.cpp
)
target_link_libraries(smm PUBLIC memilio)
target_include_directories(smm PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(smm PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
