#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <cstddef>

#include "memilio/utils/logging.h"

#include "optimization_model/optimization_model.h"
#include "optimization_settings/secirvvs_optimization.h"
#include "control_parameters/control_parameters.h"
#include "control_parameters/control_activation.h"
#include "constraints/constraints.h"

#include "helpers/integrator_selector.h"
#include "helpers/ad_type.h"

// #include "ipopt_solver/secirvvs_ipopt.h"
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "constraints/infection_state_utils.h"

// Example: ./../build/bin/secirvvs_ESID /home/jli/Memilio-Branches/memilio/data/Germany

// We use strings to gather infection state constraints:
// "Infected+Immunity Critical"
int main(int argc, char* argv[])
{
    std::cout << "Example: \"Infected+Immunity+NoSymptoms Critical\":" << std::endl;
    for (auto state : query_infection_states("Infected+Immunity+NoSymptoms Critical")) {
        std::cout << infection_state_to_string(state) << std::endl;
    }

    // ------------------------------ //
    // --- Parse 'data_directory' --- //
    // ------------------------------ //
    std::filesystem::path data_directory;
    if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [data_directory]\n";
        return 1;
    }
    if (argc == 2) {
        data_directory = argv[1];
    }
    else {
        data_directory = std::filesystem::current_path();
    }

    mio::set_log_level(mio::LogLevel::warn);

    // double INF = 1e19;

    // ---------------------------- //
    // --- Create model problem --- //
    // ---------------------------- //
    double t0   = 0.0;
    double tmax = 56.0; // 8 weeks in days
    OptimizationModel model(data_directory, t0, tmax);

    // --------------------------------------- //
    // --- Configure optimization settings --- //
    // --------------------------------------- //
    size_t num_control_intervals = 8; // Number of control intervals (e.g., one for each week)
    size_t pc_resolution         = 7; // Path constraint resolution (e.g., 7 days per control interval)
    bool random_start            = false; // Set to true to randomize initial control values

    IntegratorType integrator_type = IntegratorType::ControlledFehlberg78;
    size_t integrator_resolution   = 10;

    // ADType::Reverse can run out of tape size: Segmentation fault
    ADType ad_eval_f   = ADType::Forward;
    ADType ad_eval_jac = ADType::Forward;

    // Mapping f:[0,1]->[0,1], f(0)=0, f(1/2)=1/2, f(1)=1.
    // ControlActivation::Sigmoid can only be used if controls are in [0,1].
    ControlActivation activation = ControlActivation::Linear;

    // clang-format off
    // Control parameters: 0: No intervention, 1: Full intervention
    // { name, {lower, upper}, effectiveness, cost }
    std::vector<ControlParameter> control_parameters = {
        {"SchoolClosure",            {0.0, 1.0}, 1.00, 1.0},
        {"HomeOffice",               {0.0, 1.0}, 0.25, 1.0},
        {"PhysicalDistancingSchool", {0.0, 1.0}, 0.25, 1.0},
        {"PhysicalDistancingWork",   {0.0, 1.0}, 0.25, 1.0},
        {"PhysicalDistancingOther",  {0.0, 1.0}, 0.35, 1.0}
    };

    // Path Constraints: Active during entire simulation
    // { name, {lower, upper} }
    std::vector<Constraint> path_constraints = {
        {"Severe", {0.0, 125'000}},
        {"Critical", {0.0, 125}}
    };

    // Terminal constraints: Active at the last time step
    // { name, {lower, upper} }
    std::vector<Constraint> terminal_constraints = {
        {"Dead", {0.0, 224'000}}
    };
    // clang-format on

    SecirvvsOptimization settings(model, num_control_intervals, pc_resolution, random_start, integrator_type,
                                  integrator_resolution, ad_eval_f, ad_eval_jac, control_parameters, path_constraints,
                                  terminal_constraints, activation);

    // ----------------------------- //
    // --- Create NLP and solver --- //
    // ----------------------------- //
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp           = new Secirvvs_NLP(settings);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    int verbose            = 5;
    int print_frequency    = 1;
    bool print_timings     = true;
    int max_iter           = 300;
    double tol             = 1e-4;
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
