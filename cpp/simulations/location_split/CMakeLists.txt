file(GLOB LOCATION_SPLIT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(location_split_simulation main.cpp ${LOCATION_SPLIT_SOURCES})
target_include_directories(location_split_simulation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(location_split_simulation PRIVATE memilio abm)
target_compile_options(location_split_simulation PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

# Independent runs are executed in parallel if OpenMP is available. Without it (e.g. Apple Clang
# without libomp) the runs are executed one after another.
find_package(OpenMP COMPONENTS CXX)
if(OpenMP_CXX_FOUND)
    target_link_libraries(location_split_simulation PRIVATE OpenMP::OpenMP_CXX)
else()
    message(STATUS "location_split_simulation: OpenMP not found, runs are executed sequentially.")
endif()
