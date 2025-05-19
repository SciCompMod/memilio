#include <iostream>

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "settings.h"
#include "secirvvs_ipopt.h"

#include "ad/ad.hpp"

#include "memilio/utils/logging.h"

// ============================================================================
// --- Main Function ---------------------------------------------------------
// ============================================================================

int main()
{
    // Switch off logging for mio
    mio::set_log_level(mio::LogLevel::off);

    int num_control_intervals = 20;
    int control_interval_resolution = 5;
    double t_start = 0.0;
    double t_end = 100.0;
    int integrator_resolution = 5;

    std::vector<std::tuple<std::string, std::pair<double, double>, double>> controlBounds = {
        {"SocialDistancing", {0.05, 0.5}, 0.2},
        {"Quarantined",      {0.01, 0.3}, 0.2},
        {"TestingRate",      {0.15, 0.3}, 0.2}
    };

    std::vector<std::pair<std::string, std::pair<double, double>>> pathConstraints = {
        {"Infected+Confirmed Severe Critical", {0.0, 1'000'000.0}},
        {"Severe", {0.0, 80'000.0}},
        {"Dead", {0.0, 10'000.0}}
    };

    std::vector<std::pair<std::string, std::pair<double, double>>> terminalConstraints = {
        {"Infected+Confirmed Severe Critical", {0.0, 700'000.0}},
        {"Severe", {0.0, 60'000.0}},
        {"Critical", {0.0, 7'000.0}},
    };

    int num_age_groups = 1;
    int num_graph_nodes = 1;

    ProblemSettings settings(
        num_control_intervals, control_interval_resolution,
        t_start, t_end, integrator_resolution,
        controlBounds, pathConstraints, terminalConstraints,
        num_age_groups, num_graph_nodes
    );

    settings.print(); // optional


    std::vector<double> p(1000 * 1000);
    std::vector<double> grad(1000 * 1000);

    constraint_functions<double>(settings, p.data(), 1, grad.data(), 1);

    save_solution<double>(settings, 1, nullptr, nullptr, nullptr, 1, nullptr,nullptr,0.0);







    // auto one_state = query_infection_states("ExposedNaive"); 


    // auto states = query_infection_states("Dead+t Inf");
    // for (const auto& state : states) {
    //     std::cout << infection_state_to_string(state) << std::endl;
    // }







    // // Set up problem settings
    // ProblemSettings problem;

    // // Create NLP and solver
    // Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new Seair_NLP(problem);
    // Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // // Configure solver
    // int verbose = 5;
    // bool print_timings = true;
    // int max_iter = 500;
    // double tol = 1e-6;
    // bool useHessianApproximation = false; 
    // assert(!useHessianApproximation);
    // std::string hessianApproximation = useHessianApproximation ? "exact" : "limited-memory";
    // std::string linearSolver = "mumps";
    // std::string muStrategy = "adaptive";
    // std::string memoryUpdateType = "bfgs";

    // app->Options()->SetIntegerValue("print_level", verbose);
    // app->Options()->SetIntegerValue("print_frequency_iter", 10);
    // app->Options()->SetIntegerValue("max_iter", max_iter);
    // app->Options()->SetNumericValue("tol", tol);
    // app->Options()->SetStringValue("linear_solver", linearSolver);
    // app->Options()->SetStringValue("hessian_approximation", hessianApproximation);
    // app->Options()->SetStringValue("mu_strategy", muStrategy);
    // app->Options()->SetStringValue("limited_memory_update_type", memoryUpdateType);
    // if(print_timings){
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

    return 0;
}
