cmake_minimum_required(VERSION 3.4)

project(epidemiology-examples)

if (NOT TARGET "epidemiology")
  find_package(epidemiology REQUIRED)
endif()

set(EPI_EXAMPLES euler_test rungekutta_test seir secir secir_ageres graph)

foreach(EX ${EPI_EXAMPLES})
    add_executable(${EX} ${EX}.cpp)
    target_link_libraries(${EX} PRIVATE epidemiology)
endforeach()
