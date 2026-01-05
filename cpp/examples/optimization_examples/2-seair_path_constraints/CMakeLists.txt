 
if(MEMILIO_ENABLE_IPOPT)
    add_executable(seair_pc_individual seair_pc_individual/main.cpp seair_pc_individual/seair_ipopt.cpp)
    target_link_libraries(seair_pc_individual PRIVATE memilio ode_seair AD::AD Ipopt::Ipopt)

    add_executable(seair_pc_global_max seair_pc_global_max/main.cpp seair_pc_global_max/seair_ipopt.cpp)
    target_link_libraries(seair_pc_global_max PRIVATE memilio ode_seair AD::AD Ipopt::Ipopt)
endif()
