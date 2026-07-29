add_library(lct_sir
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(lct_sir PUBLIC memilio)
target_include_directories(lct_sir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(lct_sir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
