add_library(ode_seir
    infection_state.h
    model.h
    model.cpp
    parameters.h
    state_estimators.h
)
target_link_libraries(ode_seir PUBLIC memilio)
target_include_directories(ode_seir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_seir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
