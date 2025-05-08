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
#include "header.h"

int main() {
    int verbose = 0;
    bool timings = false;
    bool use_hessian_approximation = true;
    double tolerance = 1e-10;

    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new TutorialCpp_NLP(use_hessian_approximation);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetIntegerValue("print_level", verbose);
    app->Options()->SetNumericValue("tol", tolerance);
    app->Options()->SetStringValue("linear_solver", "mumps");

    if(use_hessian_approximation) {
        app->Options()->SetStringValue("hessian_approximation", "exact");
    } else {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }
    if(timings) {
        app->Options()->SetStringValue("timing_statistics", "yes");
        app->Options()->SetStringValue("print_timing_statistics", "yes");
        app->Options()->SetIntegerValue("print_level", 4);
    }

    if (app->Initialize() != Ipopt::ApplicationReturnStatus::Solve_Succeeded) {
        std::cerr << "Error during Ipopt initialization.\n";
        return EXIT_FAILURE;
    }

    auto start = std::chrono::high_resolution_clock::now();
    auto status = app->OptimizeTNLP(mynlp);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Optimization took " << duration.count() << " ms.\n";

    if (status == Ipopt::ApplicationReturnStatus::Solve_Succeeded) {
        std::cout << "*** The problem solved!\n";
    } else {
        std::cerr << "*** The problem FAILED!\n";
    }

    return static_cast<int>(status);
}

