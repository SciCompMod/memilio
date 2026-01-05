#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <array>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <iomanip>

#include "IpIpoptApplication.hpp"
#include "problem_setting.h"

// IPOPT documentation: https://coin-or.github.io/Ipopt/OPTIONS.html

int main()
{
    // Sets the default verbosity level for console output.
    // The valid range for this integer option is 0 <= print_level <= 12 and its default value is 5.
    int verbose = 0;
    // Switch to print timing statistics.
    bool print_timings = false;
    // Controls whether IPOPT receives second-order derivative information via eval_h.
    // - true: Use the second derivatives provided by the NLP.
    // - false: Second-order derivatives are unavailable; IPOPT uses a limited-memory quasi-Newton approximation.
    // In this example, we provide second-order derivatives, so the flag is set to true.
    // You can experiment with setting it to false. This increases the computation time as more time is spend obtaining the hessian.
    bool use_hessian_approximation = true;
    // Desired convergence tolerance (relative).
    // The valid range for this real option is 0 < tol and its default value is 1e-8.
    double tolerance = 1e-10;

    // Create the NLP (Nonlinear Programming) problem instance.
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new TutorialCpp_NLP(use_hessian_approximation);
    // Create the IPOPT application instance that will solve the NLP.
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // Configure IPOPT solver options.
    app->Options()->SetIntegerValue("print_level", verbose);
    app->Options()->SetNumericValue("tol", tolerance);
    app->Options()->SetStringValue("linear_solver", "mumps");
    if (use_hessian_approximation) {
        app->Options()->SetStringValue("hessian_approximation", "exact");
    }
    else {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }
    if (print_timings) {
        app->Options()->SetStringValue("timing_statistics", "yes");
        app->Options()->SetStringValue("print_timing_statistics", "yes");
        app->Options()->SetIntegerValue("print_level", 4);
    }
    if (app->Initialize() != Ipopt::ApplicationReturnStatus::Solve_Succeeded) {
        std::cerr << "Error during Ipopt initialization.\n";
        return EXIT_FAILURE;
    }

    auto start = std::chrono::high_resolution_clock::now();
    // Run the IPOPT solver on the NLP problem
    auto status = app->OptimizeTNLP(mynlp);
    // Record the end time of the optimization
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Optimization took " << duration.count() << " ms.\n";

    if (status == Ipopt::ApplicationReturnStatus::Solve_Succeeded) {
        std::cout << "*** The problem solved!\n";
    }
    else {
        std::cerr << "*** The problem FAILED!\n";
    }

    return static_cast<int>(status);
}
