CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
# set the project name
PROJECT(ast_builder)

find_package(LLVM REQUIRED)
include_directories(${LLVM_INCLUDE_DIRS})

# add the executable
ADD_EXECUTABLE(ast_builder ${CMAKE_CURRENT_SOURCE_DIR}/ast_builder.cpp)
target_link_libraries(ast_builder ${LLVM_LIBRARIES})