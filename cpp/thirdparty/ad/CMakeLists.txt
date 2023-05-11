cmake_minimum_required(VERSION 3.5)
project(AD CXX)

add_library(AD INTERFACE)
add_library(AD::AD ALIAS AD)

target_include_directories(AD INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

### installation related stuff ####

 install(TARGETS AD
     EXPORT ADConfig
     LIBRARY DESTINATION lib COMPONENT AD
     ARCHIVE DESTINATION lib COMPONENT AD
     RUNTIME DESTINATION bin COMPONENT AD
     INCLUDES DESTINATION include)

 install(EXPORT ADConfig
     FILE ADConfig.cmake
     NAMESPACE AD::
     DESTINATION lib/cmake/ADConfig
     COMPONENT AD)


 install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
     DESTINATION "include"
     COMPONENT AD
 )
 
