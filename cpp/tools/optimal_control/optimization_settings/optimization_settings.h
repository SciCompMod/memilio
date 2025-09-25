#pragma once

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstddef>
#include <utility>

#include "memilio/utils/logging.h"

#include "tools/optimal_control/optimization_model/optimization_model.h"
#include "tools/optimal_control/control_parameters/control_parameters.h"
#include "tools/optimal_control/control_parameters/control_activation.h"
#include "tools/optimal_control/constraints/constraints.h"

#include "tools/optimal_control/helpers/integrator_selector.h"
#include "tools/optimal_control/helpers/ad_type.h"

template <template <typename> class Model, template <typename> class Simulation, class OptimizationModel>
class OptimizationSettings
{
public:
    template <typename FP>
    using SimulationTemplate = Simulation<FP>;
    template <typename FP>
    using ModelTemplate = Model<FP>; 
    using InfectionState = typename Model<double>::Compartments;
    template <typename FP>
    using Graph = typename OptimizationModel::template Graph<FP>;
    using OptimizationModelType = OptimizationModel; 

    
    OptimizationSettings(OptimizationModel& optimization_model, size_t num_control_intervals,
                                            size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                                            size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                                            std::vector<ControlParameter> control_parameters,
                                            std::vector<Constraint> path_constraints,
                                            std::vector<Constraint> terminal_constraints, ActivationVariant activation_function,
                                            std::vector<std::pair<std::string, InfectionState>> states_strings)
        : m_optimization_model(optimization_model)
        , m_states_strings(states_strings)
        , m_t0(m_optimization_model.t0())
        , m_tmax(m_optimization_model.tmax())
        , m_integrator_type(integrator_type)
        , m_integrator_resolution(integrator_resolution)
        , m_num_control_intervals(num_control_intervals)
        , m_pc_resolution(pc_resolution)
        , m_random_start(random_start)
        , m_ad_eval_f(ad_eval_f)
        , m_ad_eval_jac(ad_eval_jac)
        , m_control_parameters(control_parameters)
        , m_path_constraints(path_constraints)
        , m_terminal_constraints(terminal_constraints)
        , m_activation_function(activation_function)
    {
        m_num_intervals            = m_num_control_intervals * m_pc_resolution;
        m_dt                       = (m_tmax - m_t0) / (m_num_intervals * m_integrator_resolution);
        m_num_control_parameters   = control_parameters.size();
        m_num_path_constraints     = path_constraints.size();
        m_num_terminal_constraints = terminal_constraints.size();
    }

    OptimizationModel& optimization_model() const
    {
        return m_optimization_model;
    }

    double t0() const
    {
        return m_t0;
    }

    double tmax() const
    {
        return m_tmax;
    }

    size_t num_control_intervals() const
    {
        return m_num_control_intervals;
    }

    size_t pc_resolution() const
    {
        return m_pc_resolution;
    }

    bool random_start() const
    {
        return m_random_start;
    }

    IntegratorType integrator_type() const
    {
        return m_integrator_type;
    }

    size_t integrator_resolution() const
    {
        return m_integrator_resolution;
    }

    ADType ad_eval_f() const
    {
        return m_ad_eval_f;
    }

    ADType ad_eval_jac() const
    {
        return m_ad_eval_jac;
    }

    std::vector<ControlParameter>& control_parameters()
    {
        return m_control_parameters;
    }

    const std::vector<ControlParameter>& control_parameters() const
    {
        return m_control_parameters;
    }

    std::vector<Constraint>& path_constraints()
    {
        return m_path_constraints;
    }

    const std::vector<Constraint>& path_constraints() const
    {
        return m_path_constraints;
    }

    std::vector<Constraint>& terminal_constraints()
    {
        return m_terminal_constraints;
    }

    const std::vector<Constraint>& terminal_constraints() const
    {
        return m_terminal_constraints;
    }

    template <typename FP>
    FP activation_function(FP x) const {
        return std::visit([&](auto&& activation_function) {
            return activation_function(x);   // std::visit for resolving the variant
        }, m_activation_function);
    }

    std::vector<std::pair<std::string, InfectionState>>& states_strings()
    {
        return m_states_strings;
    }

    const std::vector<std::pair<std::string, InfectionState>>& states_strings() const
    {
        return m_states_strings;
    }

    size_t num_intervals() const
    {
        return m_num_intervals;
    }
    
    double dt() const
    {
        return m_dt;
    }

    size_t num_control_parameters() const
    {
        return m_num_control_parameters;
    }

    size_t num_path_constraints() const
    {
        return m_num_path_constraints;
    }

    
    size_t num_terminal_constraints() const
    {
        return m_num_terminal_constraints;
    }


private:
    /* Constructor members */
    /* OptimizationModel */
    OptimizationModel& m_optimization_model;
    std::vector<std::pair<std::string, InfectionState>> m_states_strings;
    double m_t0;
    double m_tmax;
    IntegratorType m_integrator_type;
    size_t m_integrator_resolution;

    /* IPOpt */
    size_t m_num_control_intervals;
    size_t m_pc_resolution;
    bool m_random_start;
    ADType m_ad_eval_f;
    ADType m_ad_eval_jac;
    std::vector<ControlParameter> m_control_parameters;
    std::vector<Constraint> m_path_constraints;
    std::vector<Constraint> m_terminal_constraints;
    ActivationVariant m_activation_function;

    /* Helper members */
    size_t m_num_intervals;
    double m_dt;
    size_t m_num_control_parameters;
    size_t m_num_path_constraints;
    size_t m_num_terminal_constraints;
};