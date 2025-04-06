function(add_memilio_model MODEL_NAME)
    if (${MODEL_NAME} STREQUAL "abm")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/location.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/location.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/location_id.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/household.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/household.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/simulation.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/simulation.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/person.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/person.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/person_id.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/personal_rng.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/personal_rng.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/testing_strategy.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/testing_strategy.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/location_type.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/parameters.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/mobility_rules.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/mobility_rules.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/model_functions.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/model_functions.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/trip_list.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/trip_list.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/lockdown_rules.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/lockdown_rules.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/infection.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/infection.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/virus_variant.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/protection_event.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/mask.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/mask.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/abm/common_abm_loggers.h"
        )
    elseif(${MODEL_NAME} STREQUAL "d_abm")
        set(MODEL_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/d_abm/model.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/d_abm/model.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/d_abm/simulation.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/d_abm/simulation.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/d_abm/parameters.h"
        )
    elseif(${MODEL_NAME} STREQUAL "glct_secir")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/glct_secir/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/glct_secir/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/glct_secir/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/glct_secir/parameters.h"
        )
    elseif(${MODEL_NAME} STREQUAL "graph_abm")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/graph_abm/graph_abmodel.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/graph_abm/graph_abmodel.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/graph_abm/graph_abm_mobility.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/graph_abm/graph_abm_mobility.h"
        )
    elseif(${MODEL_NAME} STREQUAL "ide_secir")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_secir/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_secir/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_secir/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_secir/simulation.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_secir/simulation.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_secir/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_secir/parameters_io.h"
        )
    elseif(${MODEL_NAME} STREQUAL "ide_seir")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_seir/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_seir/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_seir/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ide_seir/infection_state.h"
        )
    elseif(${MODEL_NAME} STREQUAL "lct_secir")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/lct_secir/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/lct_secir/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/lct_secir/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/lct_secir/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/lct_secir/initializer_flows.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/lct_secir/parameters_io.h"
        )
    elseif(${MODEL_NAME} STREQUAL "ode_seair")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_seair/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_seair/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_seair/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_seair/parameters.h"
        )
    elseif(${MODEL_NAME} STREQUAL "ode_secir")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/analyze_result.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/analyze_result.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/parameter_space.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/parameter_space.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/parameters_io.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/parameters_io.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secir/model.cpp"
        )
    elseif(${MODEL_NAME} STREQUAL "ode_secirts")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/analyze_result.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/analyze_result.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/parameter_space.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/parameter_space.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/parameters_io.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/parameters_io.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirts/model.cpp"
        )
    elseif(${MODEL_NAME} STREQUAL "ode_secirvvs")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/analyze_result.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/analyze_result.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/parameter_space.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/parameter_space.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/parameters_io.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/parameters_io.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_secirvvs/model.cpp"
        )
    elseif(${MODEL_NAME} STREQUAL "ode_seir")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_seir/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_seir/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_seir/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_seir/parameters.h"
        )
    elseif(${MODEL_NAME} STREQUAL "ode_sir")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_sir/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_sir/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_sir/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/ode_sir/parameters.h"
        )
    elseif(${MODEL_NAME} STREQUAL "sde_seirvv")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_seirvv/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_seirvv/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_seirvv/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_seirvv/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_seirvv/simulation.h"
        )
    elseif(${MODEL_NAME} STREQUAL "sde_sir")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sir/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sir/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sir/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sir/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sir/simulation.h"
        )
    elseif(${MODEL_NAME} STREQUAL "sde_sirs")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sirs/infection_state.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sirs/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sirs/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sirs/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/sde_sirs/simulation.h"
        )
    elseif(${MODEL_NAME} STREQUAL "smm")
        set(MODEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/smm/parameters.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/smm/model.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/smm/model.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/smm/simulation.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/smm/simulation.cpp"
        )
    else()
        message(FATAL_ERROR "Source list not defined for model ${MODEL_NAME}")
    endif()

    # Optional extra sources handling remains the same
    foreach(EXTRA_SOURCE ${ARGN})
        list(APPEND MODEL_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${MODEL_NAME}/${EXTRA_SOURCE}")
    endforeach()

    if(NOT MODEL_SOURCES)
        message(WARNING "No sources specified for model '${MODEL_NAME}'. Skipping.")
        return()
    endif()

    add_library(${MODEL_NAME} ${MODEL_SOURCES}) # Use the explicitly set list

    target_link_libraries(${MODEL_NAME} PRIVATE memilio)

    # Includes, warnings, sanitizers remain the same
    target_include_directories(${MODEL_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> # Keep this for includes like "models/ode_secir/model.h"
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
    memilio_set_warning_flags(${MODEL_NAME})
    memilio_set_sanitizer_flags(${MODEL_NAME})
endfunction()

# --- Add Models ---
add_memilio_model(abm)
add_memilio_model(d_abm)
add_memilio_model(graph_abm)

add_memilio_model(glct_secir)
add_memilio_model(lct_secir)

add_memilio_model(ide_secir)
add_memilio_model(ide_seir)

add_memilio_model(ode_seair)
add_memilio_model(ode_secir)
add_memilio_model(ode_secirts)
add_memilio_model(ode_secirvvs)
add_memilio_model(ode_seir)
add_memilio_model(ode_sir)

add_memilio_model(sde_seirvv)
add_memilio_model(sde_sir)
add_memilio_model(sde_sirs)
add_memilio_model(smm)
