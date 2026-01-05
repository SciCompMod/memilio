 
if(MEMILIO_ENABLE_IPOPT)
    add_executable(secirvvs_dampings_opt main.cpp secirvvs_ipopt.cpp)
    target_link_libraries(secirvvs_dampings_opt PRIVATE memilio ode_secirvvs AD::AD Ipopt::Ipopt)
endif()
