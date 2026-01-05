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

// Reverse mode AD works by first recording the sequence of operations in a “tape” while evaluating the function (forward pass),
// and then propagating derivatives backward through the tape (backward pass) to compute gradients efficiently.
// Unlike forward mode, reverse mode is particularly efficient when the number of outputs is small
// and the number of inputs is large (f: R^n -> R^m, n >> m), because it computes the full gradient
// with a single backward pass rather than one forward pass per input variable.

// However, with our current AD library, reverse mode requires allocation of a global tape that stores
// all intermediate variables and operations. Because the tape is global, reverse-mode gradient computation
// is **not thread-safe**, which prevents parallel evaluation of multiple independent function calls
// or simultaneous evaluations across multiple threads. This is an important limitation compared to forward mode,
// which can be parallelized more easily since it does not rely on a shared global tape.
int main()
{
    int verbose  = 0;
    bool timings = false;
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
    if (timings) {
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
