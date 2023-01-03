add_library(abm
    location.cpp
    location.h
    household.cpp
    household.h
    simulation.cpp
    simulation.h
    person.cpp
    person.h
    testing_strategy.cpp
    testing_strategy.h
    world.cpp
    world.h
    location_type.h
    parameters.h
    parameters.cpp
    migration_rules.cpp
    migration_rules.h
    trip_list.cpp
    trip_list.h
    lockdown_rules.cpp
    lockdown_rules.h
    infection.cpp
    infection.h
    infection_state.h
    virus_variant.h
    vaccine.h
    mask.h
    mask.cpp
)
target_link_libraries(abm PUBLIC memilio)
target_include_directories(abm PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(abm PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
