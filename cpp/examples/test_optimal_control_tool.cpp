#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <cstddef>

#include "memilio/utils/logging.h"

#include "tools/optimal_control/optimization_model/optimization_model.h"
#include "tools/optimal_control/optimization_model/osecirvvs_utils.h"
#include "tools/optimal_control/optimization_settings/optimization_settings.h"
#include "tools/optimal_control/control_parameters/control_parameters.h"
#include "tools/optimal_control/control_parameters/control_activation.h"
#include "tools/optimal_control/constraints/constraints.h"
#include "tools/optimal_control/ipopt_solver/mio_ipopt.h"
#include "tools/optimal_control/helpers/integrator_selector.h"
#include "tools/optimal_control/helpers/ad_type.h"

#include "boost/filesystem.hpp"

#include "IpIpoptApplication.hpp"

int main(int argc, char* argv[])
{
    // Parse 'data_directory' 
    boost::filesystem::path data_directory;
    if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [data_directory]\n";
        return 1;
    }
    if (argc == 2) {
        data_directory = argv[1];
    }
    else {
        data_directory = boost::filesystem::current_path();
    }

    mio::set_log_level(mio::LogLevel::warn);

    // double INF = 1e19;

    // Create model problem 
    double t0   = 0.0;
    double tmax = 56.0; // 8 weeks in days
    
    enum class ContactLocation
    {
        Home = 0,
        School,
        Work,
        Other,
        Count,
    };

    const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                      {ContactLocation::School, "school_pf_eig"},
                                                                      {ContactLocation::Work, "work"},
                                                                      {ContactLocation::Other, "other"}
                                                                     };

    OptimizationModel<ContactLocation> model(data_directory, t0, tmax, true);

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
    // SigmoidActivation can only be used if controls are in [0,1].
    ActivationVariant activation = LinearActivation{};

    // clang-format off
    // Control parameters: 0: No intervention, 1: Full intervention
    // { name, {lower, upper}, effectiveness, cost }
    std::vector<ControlParameter> control_parameters = {
        {"SchoolClosure",            {0.0, 1.0}, 1.00, 1.0, mio::DampingLevel(0), mio::DampingType(0), {(size_t)ContactLocation::School}, 6},
        {"HomeOffice",               {0.0, 1.0}, 0.25, 1.0, mio::DampingLevel(0), mio::DampingType(1), {(size_t)ContactLocation::Home}, 6},
        {"PhysicalDistancingSchool", {0.0, 1.0}, 0.25, 1.0, mio::DampingLevel(1), mio::DampingType(2), {(size_t)ContactLocation::School}, 6},
        {"PhysicalDistancingWork",   {0.0, 1.0}, 0.25, 1.0, mio::DampingLevel(1), mio::DampingType(2), {(size_t)ContactLocation::Work}, 6},
        {"PhysicalDistancingOther",  {0.0, 1.0}, 0.35, 1.0, mio::DampingLevel(1), mio::DampingType(2), {(size_t)ContactLocation::Other}, 6}
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

    OptimizationSettings<mio::osecirvvs::Model, mio::osecirvvs::Simulation, OptimizationModel<ContactLocation>> settings(model, num_control_intervals, pc_resolution, random_start, integrator_type,
                                  integrator_resolution, ad_eval_f, ad_eval_jac, control_parameters, path_constraints,
                                  terminal_constraints, activation, osecirvvs_states_strings);

    // ----------------------------- //
    // --- Create NLP and solver --- //
    // ----------------------------- //
    // Ipopt::SmartPtr<Ipopt::TNLP> mynlp           = new MIO_NLP(settings);
    // Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // int verbose            = 5;
    // int print_frequency    = 1;
    // bool print_timings     = true;
    // int max_iter           = 300;
    // double tol             = 1e-4;
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

    return 0;
}
