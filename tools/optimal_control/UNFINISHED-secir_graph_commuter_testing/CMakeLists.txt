    add_executable(secir_graph_commuter_testing 
        secir_graph_commuter_testing/main.cpp 
        secir_graph_commuter_testing/optimization_model/optimization_model.cpp 
        secir_graph_commuter_testing/control_parameters/control_parameters.cpp
        secir_graph_commuter_testing/constraints/constraints.cpp
        secir_graph_commuter_testing/optimization_settings/secir_optimization.cpp
        secir_graph_commuter_testing/ipopt_solver/secir_ipopt.cpp
    )
    target_link_libraries(secir_graph_commuter_testing PRIVATE memilio ode_secir AD::AD Ipopt::Ipopt)
