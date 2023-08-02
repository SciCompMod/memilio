add_library(lct_secir
    infection_state.h
    model.h
    model.cpp
    simulation.h
    simulation.cpp    
    parameters.h
    initialization.h
    initialization.cpp
)
target_link_libraries(lct_secir PUBLIC memilio)
target_include_directories(lct_secir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(lct_secir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})