add_library(hybrid
    temporal_hybrid_model.h
    temporal_hybrid_model.cpp
    conversion_functions.h
    conversion_functions.cpp
    infection_state.h
)
target_link_libraries(hybrid PUBLIC memilio)
target_include_directories(hybrid PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(hybrid PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
