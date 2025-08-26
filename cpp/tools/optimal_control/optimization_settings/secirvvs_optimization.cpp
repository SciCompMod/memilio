#include "tools/optimal_control/optimization_settings/secirvvs_optimization.h"

#include "tools/optimal_control/control_parameters/damping_controls.h"
#include "tools/optimal_control/constraints/update_constraints.h"

SecirvvsOptimization::SecirvvsOptimization(const OptimizationModel& optimization_model, size_t num_control_intervals,
                                           size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                                           size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                                           std::vector<ControlParameter> control_parameters,
                                           std::vector<Constraint> path_constraints,
                                           std::vector<Constraint> terminal_constraints, ControlActivation activation)
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
{
    m_num_intervals            = m_num_control_intervals * m_pc_resolution;
    m_dt                       = (m_tmax - m_t0) / (m_num_intervals * m_integrator_resolution);
    m_num_control_parameters   = control_parameters.size();
    m_num_path_constraints     = path_constraints.size();
    m_num_terminal_constraints = terminal_constraints.size();

    check_constraint_feasability();
}

const OptimizationModel& SecirvvsOptimization::optimization_model() const
{
    return m_optimization_model;
}

double SecirvvsOptimization::t0() const
{
    return m_t0;
}

double SecirvvsOptimization::tmax() const
{
    return m_tmax;
}

size_t SecirvvsOptimization::num_control_intervals() const
{
    return m_num_control_intervals;
}

size_t SecirvvsOptimization::pc_resolution() const
{
    return m_pc_resolution;
}

bool SecirvvsOptimization::random_start() const
{
    return m_random_start;
}

IntegratorType SecirvvsOptimization::integrator_type() const
{
    return m_integrator_type;
}

size_t SecirvvsOptimization::integrator_resolution() const
{
    return m_integrator_resolution;
}

ADType SecirvvsOptimization::ad_eval_f() const
{
    return m_ad_eval_f;
}

ADType SecirvvsOptimization::ad_eval_jac() const
{
    return m_ad_eval_jac;
}

std::vector<ControlParameter>& SecirvvsOptimization::control_parameters()
{
    return m_control_parameters;
}

const std::vector<ControlParameter>& SecirvvsOptimization::control_parameters() const
{
    return m_control_parameters;
}

std::vector<Constraint>& SecirvvsOptimization::path_constraints()
{
    return m_path_constraints;
}

const std::vector<Constraint>& SecirvvsOptimization::path_constraints() const
{
    return m_path_constraints;
}

std::vector<Constraint>& SecirvvsOptimization::terminal_constraints()
{
    return m_terminal_constraints;
}

const std::vector<Constraint>& SecirvvsOptimization::terminal_constraints() const
{
    return m_terminal_constraints;
}

ControlActivationFunction SecirvvsOptimization::activation_function() const
{
    return m_activation_function;
}

size_t SecirvvsOptimization::num_intervals() const
{
    return m_num_intervals;
}

double SecirvvsOptimization::dt() const
{
    return m_dt;
}

size_t SecirvvsOptimization::num_control_parameters() const
{
    return m_num_control_parameters;
}

size_t SecirvvsOptimization::num_path_constraints() const
{
    return m_num_path_constraints;
}

size_t SecirvvsOptimization::num_terminal_constraints() const
{
    return m_num_terminal_constraints;
}

void SecirvvsOptimization::check_constraint_feasability()
{
    std::vector<double> path_constraint_values(num_path_constraints(), 0.0);
    std::vector<double> terminal_constraint_values(num_terminal_constraints(), 0.0);

    mio::osecirvvs::Model<double> model = optimization_model().create_model<double>();

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

    set_control_dampings<double>(*this, model, parameters);

    auto integrator = make_integrator<double>(integrator_type(), dt());

    std::vector<double> time_steps = make_time_grid<double>(t0(), tmax(), num_intervals());

    update_path_constraint<double>(*this, model, path_constraint_values);
    for (size_t interval = 0; interval < num_intervals(); interval++) {

        mio::TimeSeries<double> result = mio::simulate<double, mio::osecirvvs::Model<double>>(
            time_steps[interval], time_steps[interval + 1], dt(), model, integrator);
        const auto& final_state = result.get_last_value();

        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (size_t state_index = 0; state_index < num_infection_states(); state_index++) {
                size_t idx = age_group.get() * num_infection_states() + state_index;
                model.populations[{age_group, mio::osecirvvs::InfectionState(state_index)}] = final_state[idx];
            }
        }
        update_path_constraint<double>(*this, model, path_constraint_values);
    }
    update_terminal_constraint<double>(*this, model, terminal_constraint_values);

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
