#pragma once

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstddef>

#include "memilio/utils/logging.h"

#include "tools/optimal_control/optimization_model/optimization_model.h"
#include "tools/optimal_control/control_parameters/control_parameters.h"
#include "tools/optimal_control/control_parameters/control_activation.h"
#include "tools/optimal_control/control_parameters/damping_controls.h"
#include "tools/optimal_control/constraints/update_constraints.h"
#include "tools/optimal_control/constraints/constraints.h"
#include "tools/optimal_control/constraints/infection_state_utils.h"

#include "tools/optimal_control/helpers/integrator_selector.h"
#include "tools/optimal_control/helpers/make_time_grid.h"
#include "tools/optimal_control/helpers/ad_type.h"

template <template <typename> class Model>
class OptimizationSettings
{
public:
    template <typename FP>
    using ModelTemplate = Model<FP>;  
    using InfectionState = typename Model<double>::Compartments;

    OptimizationSettings(const OptimizationModel& optimization_model, size_t num_control_intervals,
                         size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                         size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                         std::vector<ControlParameter> control_parameters, std::vector<Constraint> path_constraints,
                         std::vector<Constraint> terminal_constraints, ControlActivation activation, std::vector<std::pair<std::string, InfectionState>> states_strings);

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

    ControlActivationFunction activation_function() const;

    size_t num_intervals() const;
    double dt() const;
    size_t num_control_parameters() const;
    size_t num_path_constraints() const;
    size_t num_terminal_constraints() const;

    void check_constraint_feasability();

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
    ControlActivationFunction m_activation_function;

    /* Helper members */
    size_t m_num_intervals;
    double m_dt;
    size_t m_num_control_parameters;
    size_t m_num_path_constraints;
    size_t m_num_terminal_constraints;
};

template <template <typename> class Model>
OptimizationSettings<Model>::OptimizationSettings(const OptimizationModel& optimization_model, size_t num_control_intervals,
                                           size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                                           size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                                           std::vector<ControlParameter> control_parameters,
                                           std::vector<Constraint> path_constraints,
                                           std::vector<Constraint> terminal_constraints, ControlActivation activation,
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
    , m_activation_function(ControlActivationFunction(activation))
    , m_states_strings(states_strings)
{
    m_num_intervals            = m_num_control_intervals * m_pc_resolution;
    m_dt                       = (m_tmax - m_t0) / (m_num_intervals * m_integrator_resolution);
    m_num_control_parameters   = control_parameters.size();
    m_num_path_constraints     = path_constraints.size();
    m_num_terminal_constraints = terminal_constraints.size();

    check_constraint_feasability();
}

template <template <typename> class Model>
const OptimizationModel& OptimizationSettings<Model>::optimization_model() const
{
    return m_optimization_model;
}

template <template <typename> class Model>
double OptimizationSettings<Model>::t0() const
{
    return m_t0;
}

template <template <typename> class Model>
double OptimizationSettings<Model>::tmax() const
{
    return m_tmax;
}

template <template <typename> class Model>
size_t OptimizationSettings<Model>::num_control_intervals() const
{
    return m_num_control_intervals;
}

template <template <typename> class Model>
size_t OptimizationSettings<Model>::pc_resolution() const
{
    return m_pc_resolution;
}

template <template <typename> class Model>
bool OptimizationSettings<Model>::random_start() const
{
    return m_random_start;
}

template <template <typename> class Model>
IntegratorType OptimizationSettings<Model>::integrator_type() const
{
    return m_integrator_type;
}

template <template <typename> class Model>
size_t OptimizationSettings<Model>::integrator_resolution() const
{
    return m_integrator_resolution;
}

template <template <typename> class Model>
ADType OptimizationSettings<Model>::ad_eval_f() const
{
    return m_ad_eval_f;
}

template <template <typename> class Model>
ADType OptimizationSettings<Model>::ad_eval_jac() const
{
    return m_ad_eval_jac;
}

template <template <typename> class Model>
std::vector<ControlParameter>& OptimizationSettings<Model>::control_parameters()
{
    return m_control_parameters;
}

template <template <typename> class Model>
const std::vector<ControlParameter>& OptimizationSettings<Model>::control_parameters() const
{
    return m_control_parameters;
}

template <template <typename> class Model>
std::vector<Constraint>& OptimizationSettings<Model>::path_constraints()
{
    return m_path_constraints;
}

template <template <typename> class Model>
const std::vector<Constraint>& OptimizationSettings<Model>::path_constraints() const
{
    return m_path_constraints;
}

template <template <typename> class Model>
std::vector<Constraint>& OptimizationSettings<Model>::terminal_constraints()
{
    return m_terminal_constraints;
}

template <template <typename> class Model>
const std::vector<Constraint>& OptimizationSettings<Model>::terminal_constraints() const
{
    return m_terminal_constraints;
}

template <template <typename> class Model>
ControlActivationFunction OptimizationSettings<Model>::activation_function() const
{
    return m_activation_function;
}

template <template <typename> class Model>
std::vector<std::pair<std::string, typename OptimizationSettings<Model>::InfectionState>>& OptimizationSettings<Model>::states_strings()
{
    return m_states_strings;
}

template <template <typename> class Model>
const std::vector<std::pair<std::string, typename OptimizationSettings<Model>::InfectionState>>& OptimizationSettings<Model>::states_strings() const
{
    return m_states_strings;
}

template <template <typename> class Model>
size_t OptimizationSettings<Model>::num_intervals() const
{
    return m_num_intervals;
}

template <template <typename> class Model>
double OptimizationSettings<Model>::dt() const
{
    return m_dt;
}

template <template <typename> class Model>
size_t OptimizationSettings<Model>::num_control_parameters() const
{
    return m_num_control_parameters;
}

template <template <typename> class Model>
size_t OptimizationSettings<Model>::num_path_constraints() const
{
    return m_num_path_constraints;
}

template <template <typename> class Model>
size_t OptimizationSettings<Model>::num_terminal_constraints() const
{
    return m_num_terminal_constraints;
}

template <template <typename> class Model>
void OptimizationSettings<Model>::check_constraint_feasability()
{
    std::vector<double> path_constraint_values(num_path_constraints(), 0.0);
    std::vector<double> terminal_constraint_values(num_terminal_constraints(), 0.0);

    Model<double> model = optimization_model().create_model<double>();

    auto set_most_restrictive_control = [&](std::vector<double>& parameters, const std::string& name) {
        size_t control_index = static_cast<size_t>(string_to_control(name));
        for (size_t control_interval = 0; control_interval < num_control_intervals(); control_interval++) {
            parameters[control_index + control_interval * num_control_parameters()] =
                control_parameters()[control_index].max();
        }
    };

    std::vector<double> parameters(num_control_parameters() * num_control_intervals());
    set_most_restrictive_control(parameters, "SchoolClosure");
    set_most_restrictive_control(parameters, "HomeOffice");
    set_most_restrictive_control(parameters, "PhysicalDistancingSchool");
    set_most_restrictive_control(parameters, "PhysicalDistancingWork");
    set_most_restrictive_control(parameters, "PhysicalDistancingOther");

    set_control_dampings<double, std::remove_reference_t<decltype(*this)>>(*this, model, parameters);

    auto integrator = make_integrator<double>(integrator_type(), dt());

    std::vector<double> time_steps = make_time_grid<double>(t0(), tmax(), num_intervals());

    update_path_constraint<double, std::remove_reference_t<decltype(*this)>>(*this, model, path_constraint_values);
    for (size_t interval = 0; interval < num_intervals(); interval++) {

        mio::TimeSeries<double> result = mio::simulate<double, Model<double>>(
            time_steps[interval], time_steps[interval + 1], dt(), model, integrator);
        const auto& final_state = result.get_last_value();

        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (size_t state_index = 0; state_index < num_infection_states<InfectionState>(); state_index++) {
                size_t idx = age_group.get() * num_infection_states<InfectionState>() + state_index;
                model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
            }
        }
        update_path_constraint<double, std::remove_reference_t<decltype(*this)>>(*this, model, path_constraint_values);
    }
    update_terminal_constraint<double, std::remove_reference_t<decltype(*this)>>(*this, model, terminal_constraint_values);

    double epsilon = 1e-3;

    // Validate and clamp path constraint values
    for (size_t i = 0; i < num_path_constraints(); ++i) {
        double path_constraint_value = path_constraint_values[i];

        Constraint& constraint = path_constraints()[i];
        double min_val         = constraint.min();
        double max_val         = constraint.max();

        if (max_val < path_constraint_value) {
            std::cout << "Path Constraint [" << i << "] \"" << constraint.name() << "\": Increasing from " << max_val
                      << " to " << path_constraint_value * (1 + epsilon) << ".\n";

            assert(0.0 <= path_constraint_value);
            constraint.set_range({min_val, path_constraint_value * (1 + epsilon)});
        }
    }

    // Validate and clamp terminal constraint values
    for (size_t i = 0; i < num_terminal_constraints(); ++i) {
        double terminal_constraint_value = terminal_constraint_values[i];

        Constraint& constraint = terminal_constraints()[i];
        double min_val         = constraint.min();
        double max_val         = constraint.max();

        if (max_val < terminal_constraint_value) {
            std::cout << "Terminal Constraint [" << i << "] \"" << constraint.name() << "\": Increasing from "
                      << max_val << " to " << terminal_constraint_value * (1 + epsilon) << ".\n";

            assert(0.0 <= terminal_constraint_value);
            constraint.set_range({min_val, terminal_constraint_value * (1 + epsilon)});
        }
    }
}