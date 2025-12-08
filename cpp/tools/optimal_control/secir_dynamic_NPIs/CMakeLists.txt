if(MEMILIO_ENABLE_IPOPT)

    # Executable for testing or running optimal control
    add_executable(optimal_control_dynamic_NPIs
        optimal_control_dynamic_NPIs.cpp
        optimization_model/optimization_model.cpp
    )

    target_link_libraries(optimal_control_dynamic_NPIs PRIVATE 
        memilio 
        ode_secir 
        Boost::boost 
        Boost::filesystem 
        AD::AD 
        Ipopt::Ipopt
    )

    target_include_directories(optimal_control_dynamic_NPIs PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/constraints>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/helpers>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

endif()
