

set(SPDLOG_INSTALL ON)

if(EPI_USE_BUNDLED_SPDLOG)
    message(STATUS "Downloading Spdlog library")
    if(CMAKE_VERSION VERSION_LESS 3.11)
        set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")

        include(DownloadProject)
        download_project(PROJ                spdlog
                         GIT_REPOSITORY      https://github.com/gabime/spdlog.git
                         GIT_TAG             v1.5.0
                         UPDATE_DISCONNECTED 1
                         QUIET
        )

        add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_SOURCE_DIR} EXCLUDE_FROM_ALL)
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
          add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR} EXCLUDE_FROM_ALL)
        endif()
    endif()
else()
    find_package(spdlog REQUIRED)
endif()

if(EPI_USE_BUNDLED_EIGEN)
    message(STATUS "Downloading Eigen library")
    if(CMAKE_VERSION VERSION_LESS 3.11)
        set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")

        include(DownloadProject)
        download_project(PROJ                eigen
                         URL                 http://gitlab.com/libeigen/eigen/-/archive/3.3.0/eigen-3.3.0.tar.gz
                         UPDATE_DISCONNECTED 1
                         QUIET
        )
    else()
        include(FetchContent)
        FetchContent_Declare(
          eigen
          URL http://gitlab.com/libeigen/eigen/-/archive/3.3.0/eigen-3.3.0.tar.gz
        )
        FetchContent_GetProperties(eigen)
        if(NOT eigen_POPULATED)
          FetchContent_Populate(eigen)
        endif()
    endif()
    add_library(eigen INTERFACE)
    target_include_directories(eigen INTERFACE ${eigen_SOURCE_DIR})
    add_library(Eigen3::Eigen ALIAS eigen)
else()
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
endif()

if (EPI_BUILD_EPI_IO)
  # TODO: support also HDF5 config mode
  find_package(HDF5 REQUIRED COMPONENTS CXX)

  find_package(tixi3 REQUIRED CONFIG)
endif(EPI_BUILD_EPI_IO)