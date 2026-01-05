#pragma once

#include "../optimization_model/optimization_model.h"
#include "../control_parameters/control_parameters.h"
#include "../control_parameters/control_activation.h"
#include "../constraints/constraints.h"
#include "../helpers/integrator_selector.h"
#include "../helpers/ad_type.h"

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstddef>

/* Stores the optimization model and the settings used in optimal control */
class SecirvvsOptimization
{
public:
    SecirvvsOptimization(const OptimizationModel& optimization_model, size_t num_control_intervals,
                         size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                         size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                         std::vector<DampingControlParameter> control_parameters,
                         std::vector<Constraint> path_constraints, std::vector<Constraint> terminal_constraints,
                         ControlActivation activation, size_t num_simulation_runs, unsigned int base_seed);

    const OptimizationModel& optimization_model() const;
    double t0() const;
    double tmax() const;

    size_t num_control_intervals() const; // Number of control intervals
    size_t pc_resolution() const; // Path constraint resolution
    bool random_start() const; // Whether to randomize initial control values
    IntegratorType integrator_type() const; // Type of integrator used
    size_t integrator_resolution() const; // Resolution of the integrator
    ADType ad_eval_f() const; // AD type for function evaluation
    ADType ad_eval_jac() const; // AD type for Jacobian evaluation

    // Control parameters used in the optimization
    std::vector<DampingControlParameter>& control_parameters();
    const std::vector<DampingControlParameter>& control_parameters() const;
    // Path constraints used in the optimization
    std::vector<Constraint>& path_constraints();
    const std::vector<Constraint>& path_constraints() const;
    // Terminal constraints used in the optimization
    std::vector<Constraint>& terminal_constraints();
    const std::vector<Constraint>& terminal_constraints() const;

    // Mapping f:[0,1]->[0,1] for control parameters
    ControlActivationFunction activation_function() const;

    size_t num_intervals() const; // Number of intervals to check for path constraints (each day).
    double dt() const; // Time step size used in the integrator.
    size_t num_control_parameters() const; // Number of control parameters.
    size_t num_path_constraints() const; // Number of path constraints.
    size_t num_terminal_constraints() const; // Number of terminal constraints.

    size_t num_simulation_runs() const; // Number of simulation runs with different random seeds.
    unsigned int base_seed() const; // Base seed for random number generation.

private:
    /* Constructor members */
    const OptimizationModel& m_optimization_model;
    double m_t0;
    double m_tmax;
    size_t m_num_control_intervals;
    size_t m_pc_resolution;
    bool m_random_start;
    IntegratorType m_integrator_type;
    size_t m_integrator_resolution;
    ADType m_ad_eval_f;
    ADType m_ad_eval_jac;
    std::vector<DampingControlParameter> m_control_parameters;
    std::vector<Constraint> m_path_constraints;
    std::vector<Constraint> m_terminal_constraints;
    ControlActivationFunction m_activation_function;
    size_t m_num_simulation_runs;
    unsigned int m_base_seed;
    /* Helper members */
    size_t m_num_intervals;
    double m_dt;
    size_t m_num_control_parameters;
    size_t m_num_path_constraints;
    size_t m_num_terminal_constraints;
};
