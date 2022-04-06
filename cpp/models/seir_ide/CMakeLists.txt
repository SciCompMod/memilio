add_library(seir_ide
    seir_ide.h
    seir_ide.cpp
)
target_link_libraries(seir_ide PUBLIC memilio)
target_include_directories(seir_ide PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(seir_ide PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})