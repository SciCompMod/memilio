add_library(sde_seirvv
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(sde_seirvv PUBLIC memilio)
target_include_directories(sde_seirvv PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(sde_seirvv PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
