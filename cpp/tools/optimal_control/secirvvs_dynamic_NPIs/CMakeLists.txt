if(MEMILIO_ENABLE_IPOPT)

    # Executable for testing or running optimal control
    add_executable(optimal_control_dynamic_NPIs
        optimal_control_dynamic_NPIs.cpp
        constraints/constraints.cpp
        control_parameters/control_parameters.cpp
        ipopt_solver/secirvvs_ipopt.cpp
        optimization_model/optimization_model.cpp
        optimization_settings/optimization_settings.cpp
    )

    target_link_libraries(optimal_control_dynamic_NPIs PRIVATE 
        memilio 
        ode_secirvvs 
        Boost::boost 
        Boost::filesystem 
        AD::AD 
        Ipopt::Ipopt
    )

    target_include_directories(optimal_control_dynamic_NPIs PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/constraints>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/control_parameters>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/helpers>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ipopt_solver>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/optimization_model>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/optimization_settings>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

endif()
