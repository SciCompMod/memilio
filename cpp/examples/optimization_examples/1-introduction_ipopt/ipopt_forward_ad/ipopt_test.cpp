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

// Gradient information of the objective function and constraints is provided accurately using
// automatic differentiation (AD). AD is a computational technique that applies the chain rule
// systematically to compute derivatives of functions defined by computer programs.
// Unlike numerical differentiation (finite differences), AD produces exact derivatives up to machine precision
// without symbolic computation.

// Forward mode AD, specifically, propagates derivatives alongside the computation of the function itself.
// For a function f: R^n -> R^m, forward mode efficiently computes the derivative of f with respect to one input variable at a time,
// which is particularly efficient when the number of inputs (variables) is small relative to the number of outputs.
// This allows IPOPT to obtain the exact gradient of the objective and constraints without manually coding the derivatives.
// In addition, since each forward-mode derivative computation for a specific input variable is independent,
// we can compute the gradients in parallel by looping over the n input variables simultaneously.
int main()
{
    int verbose        = 0;
    bool print_timings = false;
    // Controls whether IPOPT receives second-order derivative information via eval_h.
    // - true: Use the second derivatives provided by the NLP.
    // - false: Second-order derivatives are unavailable; IPOPT uses a limited-memory quasi-Newton approximation.
    // In this example, we do NOT provide exact second-order derivatives, so the flag is set to false.
    bool use_hessian_approximation = false;
    double tolerance               = 1e-10;

    Ipopt::SmartPtr<Ipopt::TNLP> mynlp           = new TutorialCpp_NLP(use_hessian_approximation);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

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

    auto start                                         = std::chrono::high_resolution_clock::now();
    auto status                                        = app->OptimizeTNLP(mynlp);
    auto end                                           = std::chrono::high_resolution_clock::now();
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
