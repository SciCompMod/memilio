# The Halle model is only used through the Python bindings (memilio.simulation.halle_abm), which compile
# halle_model.cpp themselves. This library target only keeps it compiled and warning-free in the C++ build.
add_library(abm_halle_pais halle_model.cpp halle_model.h)
target_link_libraries(abm_halle_pais PUBLIC memilio abm)
target_include_directories(abm_halle_pais PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(abm_halle_pais PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
