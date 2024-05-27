add_library(sde_seir2v
    infection_state.h
    model.h
    model.cpp
    parameters.h
    simulation.h
)
target_link_libraries(sde_seir2v PUBLIC memilio)
target_include_directories(sde_seir2v PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(sde_seir2v PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
