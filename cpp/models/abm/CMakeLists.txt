add_library(cuda_abm_interface INTERFACE interface_cuda.h time.h)

# Then define the main ABM library
add_library(abm
    location.cpp
    location.h
    household.cpp
    household.h
    simulation.cpp
    simulation.h
    person.cpp
    person.h
    person_id.h
    personal_rng.cpp
    personal_rng.h
    testing_strategy.cpp
    testing_strategy.h
    model.cpp
    model.h
    location_type.h
    parameters.h
    parameters.cpp
    mobility_rules.cpp
    mobility_rules.h
    model_functions.cpp
    model_functions.h
    trip_list.cpp
    trip_list.h
    lockdown_rules.cpp
    lockdown_rules.h
    infection.cpp
    infection.h
    infection_state.h
    virus_variant.h
    protection_event.h
    mask.h
    mask.cpp
    common_abm_loggers.h
)

# Link abm to memilio
target_link_libraries(abm PUBLIC memilio cuda_abm_interface)
target_include_directories(abm PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(abm PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

# Add CUDA implementation only if CUDA is enabled
if(MEMILIO_ENABLE_CUDA)
    
    add_library(cuda_logger STATIC common_abm_loggers_cuda.cu)
    # Set CUDA properties
    set_target_properties(cuda_logger PROPERTIES 
                         CUDA_SEPARABLE_COMPILATION ON
                         CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    
    # Link dependencies, but don't create a cycle
    target_link_libraries(cuda_logger PUBLIC ${CUDA_LIBRARIES} ${CUDA_CUDART_LIBRARY} cuda_abm_interface)
    target_include_directories(cuda_logger PRIVATE 
                             ${CUDA_INCLUDE_DIRS}
                             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)
    
    # Now link the CUDA logger to the main ABM library
    target_link_libraries(abm PUBLIC cuda_logger
    cuda_abm_interface)

    add_library(cuda_example STATIC perform_mobility_rules.cu)
    # Set CUDA properties
    set_target_properties(cuda_example PROPERTIES 
                         CUDA_SEPARABLE_COMPILATION ON
                         CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    
    # Link dependencies, but don't create a cycle
    target_link_libraries(cuda_example PUBLIC ${CUDA_LIBRARIES} ${CUDA_CUDART_LIBRARY} cuda_abm_interface)
    target_include_directories(cuda_example PRIVATE 
                             ${CUDA_INCLUDE_DIRS}
                             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)
    
    # Now link the CUDA logger to the main ABM library
    target_link_libraries(abm PUBLIC cuda_example
    cuda_abm_interface)
endif()
