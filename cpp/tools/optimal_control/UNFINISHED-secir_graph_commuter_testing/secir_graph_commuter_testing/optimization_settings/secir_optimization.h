#pragma once

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstddef>

#include "memilio/utils/logging.h"

#include "../optimization_model/optimization_model.h"
#include "../control_parameters/control_parameters.h"
#include "../constraints/constraints.h"
#include "../constraints/infection_state_utils.h"
#include "../optimization_settings/secir_optimization.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"
#include "../helpers/ad_type.h"

class SecirOptimization
{
public:
    SecirOptimization(const OptimizationModel& optimization_model, size_t simulation_days, bool random_start,
                      IntegratorType integrator_type, size_t integrator_resolution, ADType ad_eval_f,
                      ADType ad_eval_jac, ControlParameter commuter_testing_parameter,
                      std::vector<ControlParameter> dynamic_NPI_parameters, std::vector<Constraint> path_constraints,
                      std::vector<Constraint> terminal_constraints, double available_tests);

    const OptimizationModel& optimization_model() const;

    size_t simulation_days() const;
    bool random_start() const;
    IntegratorType integrator_type() const;
    size_t integrator_resolution() const;
    ADType ad_eval_f() const;
    ADType ad_eval_jac() const;

    ControlParameter& commuter_testing_parameter();
    const ControlParameter& commuter_testing_parameter() const;
    std::vector<ControlParameter>& dynamic_NPI_parameters();
    const std::vector<ControlParameter>& dynamic_NPI_parameters() const;
    std::vector<Constraint>& path_constraints();
    const std::vector<Constraint>& path_constraints() const;
    std::vector<Constraint>& terminal_constraints();
    const std::vector<Constraint>& terminal_constraints() const;

    double available_tests() const;

    double t0() const;
    double tmax() const;
    double dt() const;
    size_t num_dynamic_NPI_parameters() const;
    size_t num_path_constraints() const;
    size_t num_terminal_constraints() const;

    void check_constraint_feasability();

private:
    /* Constructor members */
    const OptimizationModel& m_optimization_model;
    size_t m_simulation_days;
    bool m_random_start;
    IntegratorType m_integrator_type;
    size_t m_integrator_resolution;
    ADType m_ad_eval_f;
    ADType m_ad_eval_jac;
    ControlParameter m_commuter_testing_parameter;
    std::vector<ControlParameter> m_dynamic_NPI_parameters;
    std::vector<Constraint> m_path_constraints;
    std::vector<Constraint> m_terminal_constraints;
    double m_available_tests;
    /* Helper members */
    double m_t0;
    double m_tmax;
    double m_dt;
    size_t m_num_dynamic_NPI_parameters;
    size_t m_num_path_constraints;
    size_t m_num_terminal_constraints;
};
