add_library(AD INTERFACE)
add_library(AD::AD ALIAS AD)

target_include_directories(AD INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(AD INTERFACE
    spdlog::spdlog
)