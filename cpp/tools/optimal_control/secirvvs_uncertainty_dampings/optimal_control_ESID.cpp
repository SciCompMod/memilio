#include "boost/filesystem.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"
#include "memilio/epidemiology/damping.h"

#include "create_control_parameters.h"
#include "create_constraints_ESID.h"
#include "check_feasibility.h"

#include "optimization_model/optimization_model.h"
#include "optimization_settings/optimization_settings.h"

#include "helpers/integrator_selector.h"
#include "helpers/ad_type.h"
#include "control_parameters/control_activation.h"
#include "constraints/constraints.h"

#include "ipopt_solver/secirvvs_ipopt.h"
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <set>

#include <chrono>

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " [data directory] [optimal_control_settings directory]\n";
        return 1;
    }

    mio::set_log_level(mio::LogLevel::critical);

    boost::filesystem::path data_directory                     = argv[1];
    boost::filesystem::path optimal_control_settings_directory = argv[2];

    boost::filesystem::path interventions_file = optimal_control_settings_directory / "intervention_list.json";
    boost::filesystem::path constraints_file   = optimal_control_settings_directory / "constraints.json";

    size_t num_age_groups = 6;
    std::vector<DampingControlParameter> control_parameters =
        create_control_parameters(interventions_file, num_age_groups);

    std::vector<Constraint> path_constraints;
    std::vector<Constraint> terminal_constraints;
    load_constraints_from_file(constraints_file, path_constraints, terminal_constraints);

    size_t simulation_weeks = 8;
    size_t days_per_week    = 7;
    size_t simulation_days  = simulation_weeks * days_per_week;

    OptimizationModel optimization_model(data_directory, simulation_days, num_age_groups);

    size_t num_simulation_runs = 50;
    unsigned int base_seed     = 42;

    auto graph_model  = optimization_model.get_graph_model<double>(base_seed);
    const auto& nodes = graph_model.nodes();

    std::cout << "Number of graph nodes: " << nodes.size() << std::endl;

    // --------------------------------------- //
    // --- Configure optimization settings --- //
    // --------------------------------------- //
    size_t num_control_intervals = simulation_weeks; // Number of control intervals (e.g., one for each week)
    size_t pc_resolution         = days_per_week; // Path constraint resolution (e.g., 7 days per control interval)
    bool random_start            = false; // Set to true to randomize initial control values

    IntegratorType integrator_type = IntegratorType::ControlledCashKarp54;
    size_t integrator_resolution   = 10;

    // ADType::Reverse can run out of tape size: Segmentation fault
    ADType ad_eval_f   = ADType::Forward;
    ADType ad_eval_jac = ADType::Forward;

    // Mapping f:[0,1]->[0,1], f(0)=0, f(1/2)=1/2, f(1)=1.
    // ControlActivation::Sigmoid can only be used if controls are in [0,1].
    ControlActivation activation = ControlActivation::Linear;

    SecirvvsOptimization settings(optimization_model, num_control_intervals, pc_resolution, random_start,
                                  integrator_type, integrator_resolution, ad_eval_f, ad_eval_jac, control_parameters,
                                  path_constraints, terminal_constraints, activation, num_simulation_runs, base_seed);

    //  check_constraint_feasibility(settings);

    // ----------------------------- //
    // --- Create NLP and solver --- //
    // ----------------------------- //
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp           = new Secirvvs_NLP(settings);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    int verbose            = 5;
    int print_frequency    = 1;
    bool print_timings     = true;
    int max_iter           = 500;
    double tol             = 1e-3;
    bool use_exact_hessian = false;
    assert(!use_exact_hessian);
    std::string hessianApproximation = use_exact_hessian ? "exact" : "limited-memory";
    std::string linearSolver         = "mumps";
    std::string muStrategy           = "adaptive";
    std::string memoryUpdateType     = "bfgs";

    app->Options()->SetIntegerValue("print_level", verbose);
    app->Options()->SetIntegerValue("print_frequency_iter", print_frequency);
    app->Options()->SetIntegerValue("max_iter", max_iter);
    app->Options()->SetNumericValue("tol", tol);
    app->Options()->SetStringValue("linear_solver", linearSolver);
    app->Options()->SetStringValue("hessian_approximation", hessianApproximation);
    app->Options()->SetStringValue("mu_strategy", muStrategy);
    app->Options()->SetStringValue("limited_memory_update_type", memoryUpdateType);
    if (print_timings) {
        app->Options()->SetStringValue("timing_statistics", "yes");
        app->Options()->SetStringValue("print_timing_statistics", "yes");
    }

    // Initialize and solve
    Ipopt::ApplicationReturnStatus status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cout << "\n*** Error during initialization!\n";
        return (int)status;
    }

    status = app->OptimizeTNLP(mynlp);

    if (status == Ipopt::Solve_Succeeded) {
        std::cout << "\n*** The problem was solved successfully!\n";
    }
    else {
        std::cout << "\n*** The problem failed to solve!\n";
    }

    return 0;
}
