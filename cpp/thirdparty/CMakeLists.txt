option(MEMILIO_USE_CONAN "Use Conan package manager to resolve dependencies." ON)

list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})

include(memilio-conan)

find_package(Boost REQUIRED CONFIG COMPONENTS filesystem)

find_package(Eigen3 REQUIRED CONFIG)

find_package(spdlog REQUIRED CONFIG)

find_package(jsoncpp CONFIG)
memilio_target_alternatives(TARGET JsonCpp::JsonCpp ALTERNATIVES jsoncpp::jsoncpp)
if (TARGET JsonCpp::JsonCpp)
    set(MEMILIO_HAS_JSONCPP ON)
else()
    message(WARNING "Memilio: JsonCpp was not found. Memilio will be built without some IO features. 
        Set CMake variable MEMILIO_USE_BUNDLED_JSONCPP to ON or install JsonCpp and set the jsoncpp_DIR cmake variable 
        to the directory containing the jsoncppConfig.cmake file to build with JsonCpp.")
endif()

find_package(HDF5 CONFIG COMPONENTS C)
if (HDF5_FOUND AND NOT TARGET HDF5::C)
    message(STATUS "Memilio: Creating ALIAS Target HDF5::C. Some HDF5 distributions don't use namespace targets.")
    add_library(HDF5::C INTERFACE IMPORTED)
    target_link_libraries(HDF5::C INTERFACE ${HDF5_C_LIBRARIES})
    target_include_directories(HDF5::C INTERFACE ${HDF5_INCLUDE_DIRS})
endif()
if (TARGET HDF5::C)
    set(MEMILIO_HAS_HDF5 ON)
else()
    message(WARNING "Memilio: HDF5 was not found. Memilio will be built without some IO features. Install HDF5 Libraries and set the HDF5_DIR cmake variable to the directory containing the hdf5-config.cmake file to build with HDF5.")
endif()
