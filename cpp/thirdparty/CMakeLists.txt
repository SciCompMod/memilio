
include(FetchContent)

option(EPI_USE_BUNDLED_SPDLOG "Use spdlog bundled with epi" ON)
mark_as_advanced(EPI_USE_BUNDLED_SPDLOG)

if(EPI_USE_BUNDLED_SPDLOG)
  # Fetch pybind11
  FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.5.0
  )
  FetchContent_GetProperties(spdlog)
  if(NOT spdlog_POPULATED)
    FetchContent_Populate(spdlog)
    add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
else()
    find_package(spdlog REQUIRED)
endif()
