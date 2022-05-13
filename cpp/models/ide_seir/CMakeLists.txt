add_library(ide_seir
    ide_seir.h
    ide_seir.cpp
)
target_link_libraries(ide_seir PUBLIC memilio)
target_include_directories(ide_seir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ide_seir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})