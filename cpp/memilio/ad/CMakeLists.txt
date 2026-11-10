# header-only library for algorithmic differentiation, see README.md
add_library(AD INTERFACE)
add_library(AD::AD ALIAS AD)

target_sources(AD INTERFACE FILE_SET HEADERS BASE_DIRS include FILES include/ad/ad.hpp)
target_link_libraries(AD INTERFACE spdlog::spdlog)
