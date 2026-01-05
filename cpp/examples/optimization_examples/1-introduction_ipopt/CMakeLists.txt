if(MEMILIO_ENABLE_IPOPT)
add_executable(ipopt_exact_derivative ipopt_exact_derivative/ipopt_test.cpp ipopt_exact_derivative/problem_setting.cpp)
    target_link_libraries(ipopt_exact_derivative PRIVATE memilio AD::AD Ipopt::Ipopt)

    add_executable(ipopt_forward_ad ipopt_forward_ad/ipopt_test.cpp ipopt_forward_ad/problem_setting.cpp)
    target_link_libraries(ipopt_forward_ad PRIVATE memilio AD::AD Ipopt::Ipopt)

    add_executable(ipopt_reverse_ad ipopt_reverse_ad/ipopt_test.cpp ipopt_reverse_ad/problem_setting.cpp)
    target_link_libraries(ipopt_reverse_ad PRIVATE memilio AD::AD Ipopt::Ipopt)
    endif()
