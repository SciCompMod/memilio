add_executable(ipopt_1 example1/ipopt_test.cpp example1/header.cpp)
target_link_libraries(ipopt_1 PRIVATE memilio AD::AD Ipopt::Ipopt)

add_executable(ipopt_2 example2/ipopt_test.cpp example2/header.cpp)
target_link_libraries(ipopt_2 PRIVATE memilio AD::AD Ipopt::Ipopt)

add_executable(ipopt_3 example3/ipopt_test.cpp example3/header.cpp)
target_link_libraries(ipopt_3 PRIVATE memilio AD::AD Ipopt::Ipopt)

add_executable(seair_1 seair1/main.cpp seair1/seair_ipopt.cpp)
target_link_libraries(seair_1 PRIVATE memilio AD::AD Ipopt::Ipopt)

add_executable(seair_2 seair2/main.cpp seair2/seair_ipopt.cpp)
target_link_libraries(seair_2 PRIVATE memilio AD::AD Ipopt::Ipopt)

add_executable(seair_3 seair3/main.cpp seair3/seair_ipopt.cpp)
target_link_libraries(seair_3 PRIVATE memilio AD::AD Ipopt::Ipopt)

add_executable(testing testing/main.cpp)
target_link_libraries(testing PRIVATE memilio AD::AD Ipopt::Ipopt)