add_executable(panvXabm 
    main_improved.cpp
    src/constants.cpp
    src/custom_loggers.cpp
    src/parameter_setter.cpp
    src/file_utils.cpp
    src/world_creator.cpp
    src/simulation_runner.cpp
    src/config_manager.cpp
    src/city_builder.cpp
    src/event_simulator.cpp
    src/multi_run_simulator.cpp
)

# Add include directories
target_include_directories(panvXabm PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../models
)

# Link required libraries
target_link_libraries(panvXabm PRIVATE 
    memilio 
    abm
    Boost::filesystem
)

# Set compiler options
target_compile_options(panvXabm PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})