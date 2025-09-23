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

template <template <typename> class Model, template <typename> class Simulation>
class OptimizationSettings
{
public:
    template <typename FP>
    using ModelTemplate = Model<FP>;
    template <typename FP>
    using SimulationTemplate = Simulation<FP>;  
    using InfectionState = typename Model<double>::Compartments;
    template <typename FP>    
    using ContactPatterns = OptimizationModel::Contactpatterns<SP>

    OptimizationSettings(const OptimizationModel& optimization_model, size_t num_control_intervals,
                         size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                         size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                         std::vector<ControlParameter> control_parameters, std::vector<Constraint> path_constraints,
                         std::vector<Constraint> terminal_constraints, ActivationVariant activation_function, std::vector<std::pair<std::string, InfectionState>> states_strings);

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
    std::vector<std::pair<std::string, InfectionState>>& states_strings();
    const std::vector<std::pair<std::string, InfectionState>>& states_strings() const;

    template <typename FP>
    FP OptimizationSettings<Model, Simulation>::activation_function(FP x) const;

    size_t num_intervals() const;
    double dt() const;
    size_t num_control_parameters() const;
    size_t num_path_constraints() const;
    size_t num_terminal_constraints() const;

private:
    /* Constructor members */
    /* OptimizationModel */
    const OptimizationModel& m_optimization_model;
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

template <template <typename> class Model, template <typename> class Simulation>
OptimizationSettings<Model, Simulation>::OptimizationSettings(const OptimizationModel& optimization_model, size_t num_control_intervals,
                                           size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                                           size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                                           std::vector<ControlParameter> control_parameters,
                                           std::vector<Constraint> path_constraints,
                                           std::vector<Constraint> terminal_constraints, ActivationVariant activation_function,
                                           std::vector<std::pair<std::string, InfectionState>> states_strings)
    : m_optimization_model(optimization_model)
    , m_t0(m_optimization_model.t0())
    , m_tmax(m_optimization_model.tmax())
    , m_num_control_intervals(num_control_intervals)
    , m_pc_resolution(pc_resolution)
    , m_random_start(random_start)
    , m_integrator_type(integrator_type)
    , m_integrator_resolution(integrator_resolution)
    , m_ad_eval_f(ad_eval_f)
    , m_ad_eval_jac(ad_eval_jac)
    , m_control_parameters(control_parameters)
    , m_path_constraints(path_constraints)
    , m_terminal_constraints(terminal_constraints)
    , m_activation_function(activation_function)
    , m_states_strings(states_strings)
{
    m_num_intervals            = m_num_control_intervals * m_pc_resolution;
    m_dt                       = (m_tmax - m_t0) / (m_num_intervals * m_integrator_resolution);
    m_num_control_parameters   = control_parameters.size();
    m_num_path_constraints     = path_constraints.size();
    m_num_terminal_constraints = terminal_constraints.size();
}

template <template <typename> class Model, template <typename> class Simulation>
const OptimizationModel& OptimizationSettings<Model, Simulation>::optimization_model() const
{
    return m_optimization_model;
}

template <template <typename> class Model, template <typename> class Simulation>
double OptimizationSettings<Model, Simulation>::t0() const
{
    return m_t0;
}

template <template <typename> class Model, template <typename> class Simulation>
double OptimizationSettings<Model, Simulation>::tmax() const
{
    return m_tmax;
}

template <template <typename> class Model, template <typename> class Simulation>
size_t OptimizationSettings<Model, Simulation>::num_control_intervals() const
{
    return m_num_control_intervals;
}

template <template <typename> class Model, template <typename> class Simulation>
size_t OptimizationSettings<Model, Simulation>::pc_resolution() const
{
    return m_pc_resolution;
}

template <template <typename> class Model, template <typename> class Simulation>
bool OptimizationSettings<Model, Simulation>::random_start() const
{
    return m_random_start;
}

template <template <typename> class Model, template <typename> class Simulation>
IntegratorType OptimizationSettings<Model, Simulation>::integrator_type() const
{
    return m_integrator_type;
}

template <template <typename> class Model, template <typename> class Simulation>
size_t OptimizationSettings<Model, Simulation>::integrator_resolution() const
{
    return m_integrator_resolution;
}

template <template <typename> class Model, template <typename> class Simulation>
ADType OptimizationSettings<Model, Simulation>::ad_eval_f() const
{
    return m_ad_eval_f;
}

template <template <typename> class Model, template <typename> class Simulation>
ADType OptimizationSettings<Model, Simulation>::ad_eval_jac() const
{
    return m_ad_eval_jac;
}

template <template <typename> class Model, template <typename> class Simulation>
std::vector<ControlParameter>& OptimizationSettings<Model, Simulation>::control_parameters()
{
    return m_control_parameters;
}

template <template <typename> class Model, template <typename> class Simulation>
const std::vector<ControlParameter>& OptimizationSettings<Model, Simulation>::control_parameters() const
{
    return m_control_parameters;
}

template <template <typename> class Model, template <typename> class Simulation>
std::vector<Constraint>& OptimizationSettings<Model, Simulation>::path_constraints()
{
    return m_path_constraints;
}

template <template <typename> class Model, template <typename> class Simulation>
const std::vector<Constraint>& OptimizationSettings<Model, Simulation>::path_constraints() const
{
    return m_path_constraints;
}

template <template <typename> class Model, template <typename> class Simulation>
std::vector<Constraint>& OptimizationSettings<Model, Simulation>::terminal_constraints()
{
    return m_terminal_constraints;
}

template <template <typename> class Model, template <typename> class Simulation>
const std::vector<Constraint>& OptimizationSettings<Model, Simulation>::terminal_constraints() const
{
    return m_terminal_constraints;
}


template <template <typename> class Model, template <typename> class Simulation>
template <typename FP>
FP OptimizationSettings<Model, Simulation>::activation_function(FP x) const {
    return std::visit([&](auto&& activation_function) {
        return activation_function(x);   // std::visit for resolving the variant
    }, m_activation_function);
}

template <template <typename> class Model, template <typename> class Simulation>
std::vector<std::pair<std::string, typename OptimizationSettings<Model, Simulation>::InfectionState>>& OptimizationSettings<Model, Simulation>::states_strings()
{
    return m_states_strings;
}

template <template <typename> class Model, template <typename> class Simulation>
const std::vector<std::pair<std::string, typename OptimizationSettings<Model, Simulation>::InfectionState>>& OptimizationSettings<Model, Simulation>::states_strings() const
{
    return m_states_strings;
}

template <template <typename> class Model, template <typename> class Simulation>
size_t OptimizationSettings<Model, Simulation>::num_intervals() const
{
    return m_num_intervals;
}

template <template <typename> class Model, template <typename> class Simulation>
double OptimizationSettings<Model, Simulation>::dt() const
{
    return m_dt;
}

template <template <typename> class Model, template <typename> class Simulation>
size_t OptimizationSettings<Model, Simulation>::num_control_parameters() const
{
    return m_num_control_parameters;
}

template <template <typename> class Model, template <typename> class Simulation>
size_t OptimizationSettings<Model, Simulation>::num_path_constraints() const
{
    return m_num_path_constraints;
}

template <template <typename> class Model, template <typename> class Simulation>
size_t OptimizationSettings<Model, Simulation>::num_terminal_constraints() const
{
    return m_num_terminal_constraints;
}
