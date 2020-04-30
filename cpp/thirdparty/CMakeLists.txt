

set(SPDLOG_INSTALL ON)

if(EPI_USE_BUNDLED_SPDLOG)
    if(CMAKE_VERSION VERSION_LESS 3.11)
        set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")

        include(DownloadProject)
        download_project(PROJ                spdlog
                         GIT_REPOSITORY      https://github.com/gabime/spdlog.git
                         GIT_TAG             v1.5.0
                         UPDATE_DISCONNECTED 1
                         QUIET
        )

        add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_SOURCE_DIR} )
    else()
        include(FetchContent)
        FetchContent_Declare(
          spdlog
          GIT_REPOSITORY https://github.com/gabime/spdlog.git
          GIT_TAG v1.5.0
        )
        FetchContent_GetProperties(spdlog)
        if(NOT spdlog_POPULATED)
          FetchContent_Populate(spdlog)
          add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR} )
        endif()
    endif()
else()
    find_package(spdlog REQUIRED)
endif()
