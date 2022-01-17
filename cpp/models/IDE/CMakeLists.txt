add_library(IDE
    IDE.h
    IDE.cpp
)
target_link_libraries(IDE PUBLIC memilio)
target_include_directories(IDE PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(IDE PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})