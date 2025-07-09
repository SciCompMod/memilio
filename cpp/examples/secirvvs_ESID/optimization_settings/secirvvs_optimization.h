#pragma once

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstddef>

#include "memilio/utils/logging.h"

#include "../optimization_model/optimization_model.h"
#include "../control_parameters/control_parameters.h"
#include "../control_parameters/control_activation.h"
#include "../constraints/constraints.h"
#include "../constraints/infection_state_utils.h"
#include "../optimization_settings/secirvvs_optimization.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"
#include "../helpers/ad_type.h"

class SecirvvsOptimization
{
public:
    SecirvvsOptimization(const OptimizationModel& optimization_model, size_t num_control_intervals,
                         size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                         size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                         std::vector<ControlParameter> control_parameters, std::vector<Constraint> path_constraints,
                         std::vector<Constraint> terminal_constraints, ControlActivation activation);

    const OptimizationModel& optimization_model() const;
    double t0() const;
    double tmax() const;

    size_t num_control_intervals() const;
    size_t pc_resolution() const;
    bool random_start() const;
    IntegratorType integrator_type() const;
    size_t integrator_resolution() const;
    ADType ad_eval_f() const;
    ADType ad_eval_jac() const;

    std::vector<ControlParameter>& control_parameters();
    const std::vector<ControlParameter>& control_parameters() const;
    std::vector<Constraint>& path_constraints();
    const std::vector<Constraint>& path_constraints() const;
    std::vector<Constraint>& terminal_constraints();
    const std::vector<Constraint>& terminal_constraints() const;

    ControlActivationFunction activation_function() const;

    size_t num_intervals() const;
    double dt() const;
    size_t num_control_parameters() const;
    size_t num_path_constraints() const;
    size_t num_terminal_constraints() const;

    void check_constraint_feasability();

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
    std::vector<ControlParameter> m_control_parameters;
    std::vector<Constraint> m_path_constraints;
    std::vector<Constraint> m_terminal_constraints;
    ControlActivationFunction m_activation_function;
    /* Helper members */
    size_t m_num_intervals;
    double m_dt;
    size_t m_num_control_parameters;
    size_t m_num_path_constraints;
    size_t m_num_terminal_constraints;
};
