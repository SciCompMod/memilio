// clang-format off

#include <iostream>
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "settings.h"
#include "secirvvs_ipopt.h"

#include "memilio/ad/ad.h"

#include <filesystem>
#include <string>
#include "memilio/io/io.h"
#include "memilio/io/mobility_io.h"

int main()
{
    constexpr double INF = 1e19;

    mio::set_log_level(mio::LogLevel::warn);

    // --------------------------------------- //
    // --- Configure optimization settings --- //
    // --------------------------------------- //
    double t0                 = 0.0;
    double tmax               = 56.0; // 8 weeks in days
    int num_control_intervals = 8;    // Number of control intervals (e.g., one for each week)
    int pc_resolution         = 7;    // Path constraint resolution (e.g., 7 days per control interval)
    int integrator_resolution = 15;   // simulation substeps
    PathConstraintMode pc_mode = PathConstraintMode::GlobalMax;
    bool random_start = false;        // Set to true to randomize initial control values

    // Vector of control bounds for various parameters: { name, {lower, upper}, initial value }
    std::vector<std::tuple<std::string, std::pair<double, double>, double>> control_bounds = {
        {"SchoolClosure", {0.0, 1.0}, 0.5},
        {"HomeOffice", {0.0, 1.0}, 0.5},
        {"PhysicalDistancingSchool", {0.0, 1.0}, 0.5}, 
        {"PhysicalDistancingWork", {0.0, 1.0}, 0.5}, 
        {"PhysicalDistancingOther", {0.0, 1.0}, 0.5}
    };
    std::vector<double> control_costs = {1.0, 1.0, 1.0, 1.0, 1.0}; // Penalty weights for each control
    
    // Vector of path constraint bounds: { name, {lower, upper} }
    std::vector<std::pair<std::string, std::pair<double, double>>> path_constraints = {
        {"Severe", {0.0, 125'000}}, 
        {"Critical", {0.0, 125}}
    };

    // Vector of terminal bounds: { name, {lower, upper} }
    std::vector<std::pair<std::string, std::pair<double, double>>> terminal_constraints = {
        {"Dead", {0.0, 224'000}}
    };

    std::random_device rd;
    std::mt19937 gen(rd());
    if (random_start) {
        for (auto& cb : control_bounds) {
            auto [name, bounds, init] = cb;
            std::uniform_real_distribution<> dis(bounds.first, bounds.second);
            init = dis(gen);
            cb = std::make_tuple(name, bounds, init);
        }
    }

    // ---------------------------------------------- //
    // --- Create problem and simulation settings --- //
    // ---------------------------------------------- //
    int num_age_groups = 6;
    std::filesystem::path data_directory = "/home/jli/Memilio-Branches/memilio/data/Germany";
    std::filesystem::path initialization_file = data_directory / "compartment_initialization_2025-05-26" / "initialization_sum.txt";

    ProblemSettings problem_settings(
        num_control_intervals, pc_resolution, t0, tmax, integrator_resolution, pc_mode,
        control_bounds, path_constraints, terminal_constraints
    );
    SimulationSettings simulation_settings(
        num_age_groups, data_directory, initialization_file
    );

    // ----------------------------- //
    // --- Create NLP and solver --- //
    // ----------------------------- //
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new Secirvvs_NLP(problem_settings, simulation_settings);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // Configure solver
    int verbose = 5;
    int print_frequency = 1;
    bool print_timings = true;
    int max_iter = 300;
    double tol = 1e-6;
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
}
