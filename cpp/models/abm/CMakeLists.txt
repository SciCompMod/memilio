add_library(abm
    location.cpp
    location.h
    household.cpp
    household.h
    simulation.cpp
    simulation.h
    person.cpp
    person.h
    world.cpp
    world.h
    state.h
    location_type.h
    parameters.h
    parameters.cpp
    migration_rules.cpp
    migration_rules.h
    lockdown_rules.cpp
    lockdown_rules.h
    testing_scheme.cpp
    testing_scheme.h
)
target_link_libraries(abm PUBLIC memilio)
target_include_directories(abm PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(abm PRIVATE ${_EPI_CXX_FLAGS_ENABLE_WARNING_ERRORS})