add_library(lsmm INTERFACE)
target_link_libraries(lsmm INTERFACE memilio)
target_include_directories(lsmm INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
