option(MEMILIO_ENABLE_IPOPT "Enable numerical optimization with Ipopt, requires a Fortran compiler." ON)

add_library(optimal_control
    constraints/constraints.cpp
    constraints/constraints.h
    constraints/update_constraints.h
    control_parameters/control_activation.h
    control_parameters/control_parameters.cpp
    control_parameters/control_parameters.h
    control_parameters/damping_controls.h
    helpers/ad_type.h
    helpers/integrator_selector.h
    helpers/make_time_grid.h
    ipopt_solver/constraint_feasibility.h
    ipopt_solver/constraint_function.h
    ipopt_solver/objective_function.h
    ipopt_solver/save_solution.h
    ipopt_solver/mio_ipopt.h
    optimization_model/optimization_model.h
    optimization_model/osecirvvs_utils.h
    optimization_settings/optimization_settings.h
)
target_link_libraries(optimal_control PUBLIC memilio ode_secirvvs Boost::boost Boost::filesystem AD::AD Ipopt::Ipopt)
target_include_directories(optimal_control PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(optimal_control PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})