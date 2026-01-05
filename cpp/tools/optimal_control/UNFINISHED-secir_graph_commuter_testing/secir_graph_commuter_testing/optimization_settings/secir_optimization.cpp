#include "secir_optimization.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include <utility>

SecirOptimization::SecirOptimization(const OptimizationModel& optimization_model, size_t simulation_days,
                                     bool random_start, IntegratorType integrator_type, size_t integrator_resolution,
                                     ADType ad_eval_f, ADType ad_eval_jac, ControlParameter commuter_testing_parameter,
                                     std::vector<ControlParameter> dynamic_NPI_parameters,
                                     std::vector<Constraint> path_constraints,
                                     std::vector<Constraint> terminal_constraints, double available_tests)
    : m_optimization_model(optimization_model)
    , m_simulation_days(simulation_days)
    , m_random_start(random_start)
    , m_integrator_type(integrator_type)
    , m_integrator_resolution(integrator_resolution)
    , m_ad_eval_f(ad_eval_f)
    , m_ad_eval_jac(ad_eval_jac)
    , m_commuter_testing_parameter(commuter_testing_parameter)
    , m_dynamic_NPI_parameters(dynamic_NPI_parameters)
    , m_path_constraints(path_constraints)
    , m_terminal_constraints(terminal_constraints)
    , m_available_tests(available_tests)
{
    m_t0   = 0.0;
    m_tmax = static_cast<double>(simulation_days);
    m_dt   = 1.0 / m_integrator_resolution;

    m_num_dynamic_NPI_parameters = dynamic_NPI_parameters.size();
    m_num_path_constraints       = path_constraints.size();
    m_num_terminal_constraints   = terminal_constraints.size();

    check_constraint_feasability();
}

const OptimizationModel& SecirOptimization::optimization_model() const
{
    return m_optimization_model;
}

size_t SecirOptimization::simulation_days() const
{
    return m_simulation_days;
}

bool SecirOptimization::random_start() const
{
    return m_random_start;
}

IntegratorType SecirOptimization::integrator_type() const
{
    return m_integrator_type;
}

size_t SecirOptimization::integrator_resolution() const
{
    return m_integrator_resolution;
}

ADType SecirOptimization::ad_eval_f() const
{
    return m_ad_eval_f;
}

ADType SecirOptimization::ad_eval_jac() const
{
    return m_ad_eval_jac;
}

ControlParameter& SecirOptimization::commuter_testing_parameter()
{
    return m_commuter_testing_parameter;
}

const ControlParameter& SecirOptimization::commuter_testing_parameter() const
{
    return m_commuter_testing_parameter;
}

std::vector<ControlParameter>& SecirOptimization::dynamic_NPI_parameters()
{
    return m_dynamic_NPI_parameters;
}

const std::vector<ControlParameter>& SecirOptimization::dynamic_NPI_parameters() const
{
    return m_dynamic_NPI_parameters;
}

std::vector<Constraint>& SecirOptimization::path_constraints()
{
    return m_path_constraints;
}

const std::vector<Constraint>& SecirOptimization::path_constraints() const
{
    return m_path_constraints;
}

std::vector<Constraint>& SecirOptimization::terminal_constraints()
{
    return m_terminal_constraints;
}

const std::vector<Constraint>& SecirOptimization::terminal_constraints() const
{
    return m_terminal_constraints;
}

double SecirOptimization::available_tests() const
{
    return m_available_tests;
}

double SecirOptimization::t0() const
{
    return m_t0;
}

double SecirOptimization::tmax() const
{
    return m_tmax;
}

double SecirOptimization::dt() const
{
    return m_dt;
}

size_t SecirOptimization::num_dynamic_NPI_parameters() const
{
    return m_num_dynamic_NPI_parameters;
}

size_t SecirOptimization::num_path_constraints() const
{
    return m_num_path_constraints;
}

size_t SecirOptimization::num_terminal_constraints() const
{
    return m_num_terminal_constraints;
}

void SecirOptimization::check_constraint_feasability()
{
    auto simulation_graph  = optimization_model().get_graph_model<double>();
    size_t num_states      = 400;
    size_t num_graph_nodes = simulation_graph.nodes().size();

    std::vector<double> path_constraint_values(num_path_constraints(), 0.0);
    std::vector<double> terminal_constraint_values(num_terminal_constraints(), 0.0);

    // Dynamic NPI parameters (e.g., school closures, distancing measures)
    std::vector<std::pair<double, double>> dynamic_NPI_values(num_dynamic_NPI_parameters());
    // Commuter testing parameters (e.g., commuter non-detection)
    // std::vector<double> commuter_testing_values(num_states, commuter_testing_parameter().min());
    std::vector<double> commuter_testing_values(num_states, 0.4);

    // Lambda to set the most restrictive control for a given NPI
    auto set_most_restrictive_control = [&](std::vector<std::pair<double, double>>& parameters,
                                            const std::string& control_name) {
        size_t index      = static_cast<size_t>(string_to_control(control_name));
        parameters[index] = {0.0, dynamic_NPI_parameters()[index].max()}; // threshold, strength
    };
    // Apply the most restrictive controls for selected NPIs
    for (const auto& dynamic_NPI : dynamic_NPI_parameters()) {
        set_most_restrictive_control(dynamic_NPI_values, dynamic_NPI.name());
    }

    // set_dynamic_NPIs<double>(*this, simulation_graph, dynamic_NPI_values);
    set_commuter_testing<double>(*this, simulation_graph, commuter_testing_values);

    for (auto& node : simulation_graph.nodes()) {
        node.property.get_simulation().set_integrator_core(std::move(make_integrator<double>(integrator_type(), dt())));
    }

    std::vector<double> time_steps = make_time_grid<double>(t0(), tmax(), simulation_days());
    auto graph_sim_mobility        = mio::make_mobility_sim<double>(t0(), dt(), std::move(simulation_graph));

    update_path_constraint(*this, graph_sim_mobility, path_constraint_values);
    for (size_t day = 0; day < simulation_days(); day++) {
        std::cout << day << std::endl;
        graph_sim_mobility.advance(time_steps[day + 1]);
        update_path_constraint(*this, graph_sim_mobility, path_constraint_values);
    }
    update_terminal_constraint(*this, graph_sim_mobility, terminal_constraint_values);

    double tests_used = get_tests_used(*this, graph_sim_mobility);

    std::cout << "Tests used in feasibility check: " << tests_used << "/" << available_tests() << "\n";

    double epsilon = 1e-3;

    // Validate and clamp path constraint values
    for (size_t i = 0; i < num_path_constraints(); ++i) {
        double path_constraint_value = path_constraint_values[i];

        Constraint& constraint = path_constraints()[i];
        double min_val         = constraint.min();
        double max_val         = constraint.max();

        std::cout << "Path Constraint [" << i << "] \"" << constraint.name() << "\": value " << path_constraint_value
                  << ".\n";

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

        std::cout << "Terminal Constraint [" << i << "] \"" << constraint.name() << "\": value "
                  << terminal_constraint_value << ".\n";

        if (max_val < terminal_constraint_value) {
            std::cout << "Terminal Constraint [" << i << "] \"" << constraint.name() << "\": Increasing from "
                      << max_val << " to " << terminal_constraint_value * (1 + epsilon) << ".\n";

            assert(0.0 <= terminal_constraint_value);
            constraint.set_range({min_val, terminal_constraint_value * (1 + epsilon)});
        }
    }
}
