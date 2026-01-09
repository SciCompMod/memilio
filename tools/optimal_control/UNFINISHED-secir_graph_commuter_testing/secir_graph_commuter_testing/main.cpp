#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <typeindex>
#include <unordered_map>
#include <any>
#include <mutex>

#include "memilio/utils/logging.h"
#include "memilio/io/mobility_io.h"

#include "optimization_model/optimization_model.h"
#include "helpers/integrator_selector.h"
#include "helpers/ad_type.h"
#include "control_parameters/control_parameters.h"
#include "constraints/constraints.h"
#include "optimization_settings/secir_optimization.h"

#include "ipopt_solver/secir_ipopt.h"
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "control_parameters/damping_controls.h"
#include "constraints/update_constraints.h"

int main(int argc, char* argv[])
{
    // ------------------------------ //
    // --- Parse 'data_directory' --- //
    // ------------------------------ //
    boost::filesystem::path data_directory                     = argv[1];
    boost::filesystem::path optimal_control_settings_directory = argv[2];

    mio::set_log_level(mio::LogLevel::warn);
    double INF = 1e19;

    // ---------------------------- //
    // --- Create model problem --- //
    // ---------------------------- //
    int simulation_days = 56;

    // General dynamic NPIs settings
    double npi_duration = 14.0;
    double npi_interval = 4.0;
    double npi_base     = 100'000; // Measured per 100'000 people

    size_t num_age_groups = 1;

    OptimizationModel optimization_model(data_directory, optimal_control_settings_directory, simulation_days,
                                         num_age_groups, npi_duration, npi_interval, npi_base);
    auto simulation_graph = optimization_model.get_graph_model<double>();

    const auto& nodes = simulation_graph.nodes();
    const auto& edges = simulation_graph.edges();
    std::cout << "Graph info: nodes=" << nodes.size() << ", edges=" << edges.size() << ", edge density="
              << (nodes.size() > 1 ? 100.0 * edges.size() / (nodes.size() * (nodes.size() - 1)) : 0.0) << "%\n";

    IntegratorType integrator_type = IntegratorType::ControlledCashKarp54;
    size_t integrator_resolution   = 10;

    double t0   = 0.0;
    double tmax = simulation_days;
    double dt   = 1.0 / integrator_resolution;

    // --------------------------------------- //
    // --- Configure optimization settings --- //
    // --------------------------------------- //
    // Set to true to randomize initial control values
    bool random_start = false;

    // ADType::Reverse can run out of tape size: Segmentation fault
    ADType ad_eval_f   = ADType::Forward;
    ADType ad_eval_jac = ADType::Forward;

    // -----------------------------------------------------------------------------
    // Commuter Testing
    // -----------------------------------------------------------------------------
    // Test detection probability: proportion of infected commuters correctly detected
    double test_detection = 0.75; // 75% detection rate
    //
    // Effective non-detection rate per week (commuter_nondetection):
    // commuter_nondetection = 1.0 - test_detection * test_days / 5
    //
    // Examples (5-day working week):
    //   test_days = 0 → 1.00 (no testing)
    //   test_days = 1 → 0.85 (tested once per week)
    //   test_days = 2 → 0.70 (tested twice per week)
    //   test_days = 3 → 0.55 (tested three times per week)
    //   test_days = 4 → 0.40 (tested four times per week)
    //   test_days = 5 → 0.25 (tested daily)
    //
    // Control parameter for commuter testing:
    //   - Range: [1.0 - test_detection, 1.0]
    //   - Effectiveness: test_detection
    //   - Cost: 1.0
    ControlParameter commuter_testing_parameter = {"CommuterTesting", {1.0 - test_detection, 1.0}, test_detection, 1.0};
    double available_tests                      = 200'000;

    // -----------------------------------------------------------------------------
    // Dynamic NPI Parameters
    // -----------------------------------------------------------------------------
    // General format:
    //   { name, {lower_bound, upper_bound}, effectiveness, cost }
    //
    // Notes:
    //   - Damping values: 0.0 = no intervention, 1.0 = full intervention
    //   - Cost is currently normalized to 1.0 for all interventions
    //
    // Interventions included:
    //   - School closures
    //   - Home office adoption
    //   - Physical distancing (school, work, other settings)
    std::vector<ControlParameter> dynamic_npi_parameters = {
        // {"SchoolClosure", {0.0, 1.0}, 1.00, 1.0},
        // {"HomeOffice", {0.0, 1.0}, 0.25, 1.0},
        // {"PhysicalDistancingSchool", {0.0, 1.0}, 0.25, 1.0},
        // {"PhysicalDistancingWork", {0.0, 1.0}, 0.25, 1.0},
        // {"PhysicalDistancingOther", {0.0, 1.0}, 0.35, 1.0},
    };

    // -----------------------------------------------------------------------------
    // Path Constraints
    // -----------------------------------------------------------------------------
    // Purpose:
    //   Constraints that remain active throughout the entire simulation.
    // Format:
    //   { name, {lower_bound, upper_bound}, node_id (optional) }
    // Notes:
    //   - If node_id is std::nullopt → applies globally
    //   - Otherwise, constraint is specific to the given node
    std::vector<Constraint> path_constraints = {};

    // -----------------------------------------------------------------------------
    // Terminal Constraints
    // -----------------------------------------------------------------------------
    // Purpose:
    //   Constraints that apply only at the final time step of the simulation.
    // Format:
    //   { name, {lower_bound, upper_bound}, node_id (optional) }
    // Notes:
    //   - If node_id is std::nullopt → applies globally
    //   - Otherwise, constraint is specific to the given node
    std::vector<Constraint> terminal_constraints = {{"InfectedSymptoms", {0.0, 1'500'000}, std::nullopt}};

    SecirOptimization settings(optimization_model, simulation_days, random_start, integrator_type,
                               integrator_resolution, ad_eval_f, ad_eval_jac, commuter_testing_parameter,
                               dynamic_npi_parameters, path_constraints, terminal_constraints, available_tests);

    // // ----------------------------- //
    // // --- Create NLP and solver --- //
    // // ----------------------------- //
    // Ipopt::SmartPtr<Ipopt::TNLP> mynlp           = new Secir_NLP(settings);
    // Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // int verbose            = 5;
    // int print_frequency    = 1;
    // bool print_timings     = true;
    // int max_iter           = 25;
    // double tol             = 1e-2;
    // bool use_exact_hessian = false;
    // assert(!use_exact_hessian);
    // std::string hessianApproximation = use_exact_hessian ? "exact" : "limited-memory";
    // std::string linearSolver         = "mumps";
    // std::string muStrategy           = "adaptive";
    // std::string memoryUpdateType     = "bfgs";

    // app->Options()->SetIntegerValue("print_level", verbose);
    // app->Options()->SetIntegerValue("print_frequency_iter", print_frequency);
    // app->Options()->SetIntegerValue("max_iter", max_iter);
    // app->Options()->SetNumericValue("tol", tol);
    // app->Options()->SetStringValue("linear_solver", linearSolver);
    // app->Options()->SetStringValue("hessian_approximation", hessianApproximation);
    // app->Options()->SetStringValue("mu_strategy", muStrategy);
    // app->Options()->SetStringValue("limited_memory_update_type", memoryUpdateType);
    // if (print_timings) {
    //     app->Options()->SetStringValue("timing_statistics", "yes");
    //     app->Options()->SetStringValue("print_timing_statistics", "yes");
    // }

    // // Initialize and solve
    // Ipopt::ApplicationReturnStatus status = app->Initialize();
    // if (status != Ipopt::Solve_Succeeded) {
    //     std::cout << "\n*** Error during initialization!\n";
    //     return (int)status;
    // }

    // status = app->OptimizeTNLP(mynlp);

    // if (status == Ipopt::Solve_Succeeded) {
    //     std::cout << "\n*** The problem was solved successfully!\n";
    // }
    // else {
    //     std::cout << "\n*** The problem failed to solve!\n";
    // }

    // // Run Model

    // {

    //     auto simulation_graph  = settings.optimization_model().get_graph_model<double>();
    //     size_t num_graph_nodes = simulation_graph.nodes().size();

    //     std::vector<double> path_constraint_values(settings.num_path_constraints(), 0.0);
    //     std::vector<double> terminal_constraint_values(settings.num_terminal_constraints(), 0.0);

    //     // Dynamic NPI parameters (e.g., school closures, distancing measures)
    //     std::vector<std::pair<double, double>> dynamic_NPI_values(settings.num_dynamic_NPI_parameters());
    //     // Commuter testing parameters (e.g., commuter non-detection)
    //     // std::vector<double> commuter_testing_values(num_graph_nodes, settings.commuter_testing_parameter().min());

    //     std::vector<double> commuter_testing_values(num_graph_nodes, 0.25);

    //     // 0.25 -> 175.861 used
    //     // 0.4 -> 177.052 used
    //     // 0.55 -> 178.623 used | 242.789 used
    //     // 0.85 -> 139.165 used
    //     // 0.95 -> 69.366.7 used

    //     // Lambda to set the most restrictive control for a given NPI
    //     auto set_most_restrictive_control = [&](std::vector<std::pair<double, double>>& parameters,
    //                                             const std::string& control_name) {
    //         size_t index      = static_cast<size_t>(string_to_control(control_name));
    //         parameters[index] = {0.0, settings.dynamic_NPI_parameters()[index].min()}; // threshold, strength
    //     };
    //     // Apply the most restrictive controls for selected NPIs
    //     for (const auto& dynamic_NPI : settings.dynamic_NPI_parameters()) {
    //         set_most_restrictive_control(dynamic_NPI_values, dynamic_NPI.name());
    //     }

    //     set_dynamic_NPIs<double>(settings, simulation_graph, dynamic_NPI_values);
    //     set_commuter_testing<double>(settings, simulation_graph, commuter_testing_values);

    //     for (auto& node : simulation_graph.nodes()) {
    //         node.property.get_simulation().set_integrator_core(
    //             std::move(make_integrator<double>(settings.integrator_type(), settings.dt())));
    //     }

    //     std::vector<double> time_steps =
    //         make_time_grid<double>(settings.t0(), settings.tmax(), settings.simulation_days());
    //     auto graph_sim_mobility =
    //         mio::make_mobility_sim<double>(settings.t0(), settings.dt(), std::move(simulation_graph));

    //     update_path_constraint(settings, graph_sim_mobility, path_constraint_values);
    //     for (size_t day = 0; day < settings.simulation_days(); day++) {
    //         graph_sim_mobility.advance(time_steps[day + 1]);
    //         update_path_constraint(settings, graph_sim_mobility, path_constraint_values);
    //     }
    //     update_terminal_constraint(settings, graph_sim_mobility, terminal_constraint_values);

    //     double tests_used = get_tests_used(settings, graph_sim_mobility);

    //     std::cout << "Tests used: " << tests_used << " / " << settings.available_tests() << "\n";
    // }

    return 0;
}
