#include <iostream>

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "settings.h"
#include "seair_ipopt.h"

#include "ad/ad.hpp"

int main()
{
    // Switch off logging for mio
    mio::set_log_level(mio::LogLevel::off);

    // Set up problem settings
    ProblemSettings problem;

    // Create NLP and solver
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp           = new Seair_NLP(problem);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // Configure solver
    int verbose                  = 5;
    bool print_timings           = true;
    int max_iter                 = 500;
    double tol                   = 1e-6;
    bool useHessianApproximation = false;
    assert(!useHessianApproximation);
    std::string hessianApproximation = useHessianApproximation ? "exact" : "limited-memory";
    std::string linearSolver         = "mumps";
    std::string muStrategy           = "adaptive";
    std::string memoryUpdateType     = "bfgs";

    app->Options()->SetIntegerValue("print_level", verbose);
    app->Options()->SetIntegerValue("print_frequency_iter", 1);
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
