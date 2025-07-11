if(MEMILIO_ENABLE_IPOPT)

    add_executable(ipopt_exact_derivative ipopt_exact_derivative/ipopt_test.cpp ipopt_exact_derivative/problem_setting.cpp)
    target_link_libraries(ipopt_exact_derivative PRIVATE memilio AD::AD Ipopt::Ipopt)

    add_executable(ipopt_forward_ad ipopt_forward_ad/ipopt_test.cpp ipopt_forward_ad/problem_setting.cpp)
    target_link_libraries(ipopt_forward_ad PRIVATE memilio AD::AD Ipopt::Ipopt)

    add_executable(ipopt_reverse_ad ipopt_reverse_ad/ipopt_test.cpp ipopt_reverse_ad/problem_setting.cpp)
    target_link_libraries(ipopt_reverse_ad PRIVATE memilio AD::AD Ipopt::Ipopt)


    add_executable(seair_pc_individual seair_pc_individual/main.cpp seair_pc_individual/seair_ipopt.cpp)
    target_link_libraries(seair_pc_individual PRIVATE memilio ode_seair AD::AD Ipopt::Ipopt)

    add_executable(seair_pc_global_max seair_pc_global_max/main.cpp seair_pc_global_max/seair_ipopt.cpp)
    target_link_libraries(seair_pc_global_max PRIVATE memilio ode_seair AD::AD Ipopt::Ipopt)


    add_executable(secirvvs_pc_testing secirvvs_pc_testing/main.cpp secirvvs_pc_testing/secirvvs_ipopt.cpp)
    target_link_libraries(secirvvs_pc_testing PRIVATE memilio ode_secirvvs AD::AD Ipopt::Ipopt)

    add_executable(secirvvs_ESID_simulate 
        secirvvs_ESID_simulate/main.cpp 
        secirvvs_ESID_simulate/optimization_model/optimization_model.cpp
        secirvvs_ESID_simulate/optimization_settings/secirvvs_optimization.cpp
        secirvvs_ESID_simulate/ipopt_solver/secirvvs_ipopt.cpp
        secirvvs_ESID_simulate/control_parameters/control_parameters.cpp
        secirvvs_ESID_simulate/constraints/constraints.cpp
    )
    target_link_libraries(secirvvs_ESID_simulate PRIVATE memilio ode_secirvvs AD::AD Ipopt::Ipopt)

    add_executable(secirvvs_ESID_simulation 
        secirvvs_ESID_simulation/main.cpp 
        secirvvs_ESID_simulation/optimization_model/optimization_model.cpp
        secirvvs_ESID_simulation/optimization_settings/secirvvs_optimization.cpp
        secirvvs_ESID_simulation/ipopt_solver/secirvvs_ipopt.cpp
        secirvvs_ESID_simulation/control_parameters/control_parameters.cpp
        secirvvs_ESID_simulation/constraints/constraints.cpp
    )
    target_link_libraries(secirvvs_ESID_simulation PRIVATE memilio ode_secirvvs AD::AD Ipopt::Ipopt)

    add_executable(secirvvs_graph 
        secirvvs_graph/main.cpp 
        secirvvs_graph/optimization_model/optimization_model.cpp 
        secirvvs_graph/control_parameters/control_parameters.cpp
        secirvvs_graph/constraints/constraints.cpp
        secirvvs_graph/optimization_settings/secirvvs_optimization.cpp
        secirvvs_graph/ipopt_solver/secirvvs_ipopt.cpp
    )
    target_link_libraries(secirvvs_graph PRIVATE memilio ode_secirvvs AD::AD Ipopt::Ipopt)

endif()