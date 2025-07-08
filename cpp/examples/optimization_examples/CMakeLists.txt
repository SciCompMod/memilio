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

    add_executable(secirvvs_dampings secirvvs_dampings/main.cpp secirvvs_dampings/secirvvs_ipopt.cpp)
    target_link_libraries(secirvvs_dampings PRIVATE memilio ode_secirvvs AD::AD Ipopt::Ipopt)

    add_executable(secirvvs_graph 
        secirvvs_graph/main.cpp 
        secirvvs_graph/secirvvs_optimization.cpp 
        secirvvs_graph/secirvvs_ipopt.cpp
        secirvvs_graph/validate_constraints.cpp
    )
    target_link_libraries(secirvvs_graph PRIVATE memilio ode_secirvvs AD::AD Ipopt::Ipopt)

    # add_executable(secirvvs_ESID 
    #     secirvvs_ESID/main.cpp 
    #     secirvvs_ESID/optimization_model.cpp
    #     secirvvs_ESID/secirvvs_optimization.cpp
    #     secirvvs_ESID/secirvvs_ipopt.cpp
    # )
    # target_link_libraries(secirvvs_ESID PRIVATE memilio ode_secirvvs AD::AD Ipopt::Ipopt)


    add_executable(secirvvs_ESID_2 
        secirvvs_ESID_2/main.cpp 
        secirvvs_ESID_2/optimization_model/optimization_model.cpp
        secirvvs_ESID_2/optimization_settings/secirvvs_optimization.cpp
        secirvvs_ESID_2/ipopt_solver/secirvvs_ipopt.cpp
        secirvvs_ESID_2/control_parameters/control_parameters.cpp
        secirvvs_ESID_2/constraints/constraints.cpp
    )
    target_link_libraries(secirvvs_ESID_2 PRIVATE memilio ode_secirvvs AD::AD Ipopt::Ipopt)



endif()