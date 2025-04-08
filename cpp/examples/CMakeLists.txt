#-------------------------------------------------------
# Helper Function for Examples
#-------------------------------------------------------
function(add_memilio_example NAME SOURCES)
    set(EXTRA_DEPS ${ARGN})
    add_executable(${NAME} ${SOURCES})
    target_link_libraries(${NAME} PRIVATE memilio ${EXTRA_DEPS})
    memilio_set_warning_flags(${NAME})
    memilio_set_sanitizer_flags(${NAME})
endfunction()

#-------------------------------------------------------
# Example Executable Definitions
#-------------------------------------------------------
add_memilio_example(memilio_example_euler            euler_test.cpp)
add_memilio_example(memilio_example_adapt_rk         adapt_rk_test.cpp)
add_memilio_example(memilio_example_history          history.cpp)

add_memilio_example(memilio_example_ode_sir                  ode_sir.cpp                  ode_sir)
add_memilio_example(memilio_example_ode_sir_ageres           ode_sir_ageres.cpp           ode_sir)
add_memilio_example(memilio_example_ode_seir                 ode_seir.cpp                 ode_seir)
add_memilio_example(memilio_example_ode_seir_ageres          ode_seir_ageres.cpp          ode_seir)
add_memilio_example(memilio_example_seir_flows               ode_seir_flows.cpp           ode_seir)
add_memilio_example(memilio_example_ode_secir                ode_secir.cpp                ode_secir)
add_memilio_example(memilio_example_ode_secir_ageres         ode_secir_ageres.cpp         ode_secir)
add_memilio_example(memilio_example_ode_secir_contacts       ode_secir_contact_changes.cpp ode_secir)
add_memilio_example(memilio_example_ode_secir_feedback       ode_secir_feedback.cpp       ode_secir)
add_memilio_example(memilio_example_ode_secir_graph          ode_secir_graph.cpp          ode_secir)
add_memilio_example(memilio_example_ode_secirvvs             ode_secirvvs.cpp             ode_secirvvs)
add_memilio_example(memilio_example_ode_secirts              ode_secirts.cpp              ode_secirts)
add_memilio_example(memilio_example_graph_stochastic_mobility graph_stochastic_mobility.cpp ode_secir)
add_memilio_example(memilio_example_ode_secir_param_sampling ode_secir_parameter_sampling.cpp ode_secir)

add_memilio_example(memilio_example_sde_sir    sde_sir.cpp    sde_sir)
add_memilio_example(memilio_example_sde_sirs   sde_sirs.cpp   sde_sirs)
add_memilio_example(memilio_example_sde_seirvv sde_seirvv.cpp sde_seirvv)

add_memilio_example(memilio_example_ad_square  ad_square_example.cpp)
add_memilio_example(memilio_example_ad_odeint  ad_odeint_example.cpp Boost::boost)
add_memilio_example(memilio_example_ode_seair  ode_seair.cpp         ode_seair)

add_memilio_example(memilio_example_ide_seir        ide_seir.cpp        ide_seir)
add_memilio_example(memilio_example_ide_secir       ide_secir.cpp       ide_secir)
add_memilio_example(memilio_example_ide_secir_ageres ide_secir_ageres.cpp ide_secir)

add_memilio_example(memilio_example_lct_secir  lct_secir.cpp  lct_secir)
add_memilio_example(memilio_example_glct_secir glct_secir.cpp glct_secir)

add_memilio_example(memilio_example_abm_minimal  abm_minimal.cpp      abm)
add_memilio_example(memilio_example_abm_history  abm_history_object.cpp abm)
add_memilio_example(memilio_example_graph_abm    graph_abm.cpp        graph_abm abm)
add_memilio_example(memilio_example_dabm         d_abm.cpp            d_abm)

add_memilio_example(memilio_example_smm smm.cpp smm)

if(MEMILIO_ENABLE_IPOPT)
    message(STATUS "Building examples requiring Ipopt.")
    add_memilio_example(memilio_example_ode_seair_optimization ode_seair_optimization.cpp ode_seair Ipopt::Ipopt)
endif()

if(MEMILIO_HAS_JSONCPP)
    message(STATUS "Building examples requiring JsonCpp.")
    add_memilio_example(memilio_example_cli         cli.cpp                   JsonCpp::JsonCpp)
    add_memilio_example(memilio_example_serialize   serialize.cpp             JsonCpp::JsonCpp)
    add_memilio_example(memilio_example_ode_read_graph ode_secir_read_graph.cpp ode_secir JsonCpp::JsonCpp)
    add_memilio_example(memilio_example_ide_init    ide_initialization.cpp    ide_secir JsonCpp::JsonCpp)
endif()

if(MEMILIO_HAS_HDF5)
    message(STATUS "Building examples requiring HDF5.")
    add_memilio_example(memilio_example_ode_save_results ode_secir_save_results.cpp ode_secir HDF5::HDF5)
endif()

if(MEMILIO_HAS_HDF5 AND MEMILIO_HAS_JSONCPP)
    message(STATUS "Building examples requiring HDF5 and JsonCpp.")
    add_memilio_example(memilio_example_ode_param_study      ode_secir_parameter_study.cpp       ode_secir HDF5::HDF5 JsonCpp::JsonCpp)
    add_memilio_example(memilio_example_ode_param_study_graph ode_secir_parameter_study_graph.cpp ode_secir HDF5::HDF5 JsonCpp::JsonCpp)
endif()
