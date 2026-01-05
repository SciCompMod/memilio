#include "problem_setting.h"

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

#include "IpIpoptData.hpp"

// ---------------
// TutorialCpp_NLP
// ---------------
// This tutorial demonstrates how to define and solve a nonlinear programming (NLP) problem
// using IPOPT in C++. The problem we solve here is a standard test problem with:
// - 4 decision variables: x0, x1, x2, x3
//     1 <= xi <= 5 for all decision variables
// - 2 constraints:
//     g0(x) = x0 * x1 * x2 * x3 >= 25
//     g1(x) = x0^2 + x1^2 + x2^2 + x3^2 == 40
// - Objective function:
//     f(x) = x0 * x3 * (x0 + x1 + x2) + x2

// This example demonstrates:
// 1. Defining variable and constraint bounds
// 2. Providing starting points
// 3. Computing objective, gradient, constraints, Jacobian, and Hessian
// 4. Using exact second-order derivatives or quasi-Newton approximation
// 5. Retrieving and displaying the solution

constexpr Ipopt::Number INF = 1e19;

// Constructor: set up problem dimensions, bounds, and Hessian option
TutorialCpp_NLP::TutorialCpp_NLP(bool use_hessian_approximation) noexcept
{
    n_                         = 4; // number of decision variables
    m_                         = 2; // number of constraints
    nnz_jac_g_                 = n_ * m_; // full Jacobian (2x4)
    nnz_h_lag_                 = use_hessian_approximation ? n_ * (n_ + 1) / 2 : 0; // lower-triangle of Hessian
    use_hessian_approximation_ = use_hessian_approximation;

    // Variable bounds
    x_l_.assign(n_, 1.0);
    x_u_.assign(n_, 5.0);

    // Constraint bounds
    g_l_.resize(m_);
    g_u_.resize(m_);
    g_l_[0] = 25.0;
    g_u_[0] = INF;
    g_l_[1] = 40.0;
    g_u_[1] = 40.0;
}

// Provide IPOPT with basic problem info (dimensions, sparsity)
bool TutorialCpp_NLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                                   Ipopt::TNLP::IndexStyleEnum& index_style)
{
    n           = n_;
    m           = m_;
    nnz_jac_g   = nnz_jac_g_;
    nnz_h_lag   = nnz_h_lag_;
    index_style = Ipopt::TNLP::C_STYLE;
    return true;
}

// Provide bounds for variables and constraints
bool TutorialCpp_NLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m,
                                      Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    assert(n == n_ && m == m_);
    for (Ipopt::Index i = 0; i < n_; ++i) {
        x_l[i] = x_l_[i];
        x_u[i] = x_u_[i];
    }

    for (Ipopt::Index i = 0; i < m_; ++i) {
        g_l[i] = g_l_[i];
        g_u[i] = g_u_[i];
    }

    return true;
}

// Provide a starting point for the solver
bool TutorialCpp_NLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
                                         Ipopt::Number* /*z_L*/, Ipopt::Number* /*z_U*/, Ipopt::Index m,
                                         bool init_lambda, Ipopt::Number* /*lambda*/
)
{
    assert(init_x && !init_z && !init_lambda);
    assert(n == n_ && m == m_);

    // IPOPT does still work even if the starting point is infeasible.
    for (Ipopt::Index i = 0; i < n; ++i) {
        if (std::isfinite(x_l_[i]) && std::isfinite(x_u_[i])) {
            // Finite lower and upper bounds -> start in the middle
            x[i] = 0.5 * (x_l_[i] + x_u_[i]);
        }
        else if (std::isfinite(x_l_[i])) {
            // Only lower bound -> start a bit above
            x[i] = x_l_[i] + 1.0;
        }
        else if (std::isfinite(x_u_[i])) {
            // Only upper bound -> start a bit below
            x[i] = x_u_[i] - 1.0;
        }
        else {
            // No bounds -> default to zero
            x[i] = 0.0;
        }
    }

    return true;
}

// Compute the objective function value
bool TutorialCpp_NLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& obj_value)
{
    assert(n == n_);
    obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
    return true;
}

// Compute the gradient of the objective
bool TutorialCpp_NLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f)
{
    assert(n == n_);

    grad_f[0] = x[3] * (2 * x[0] + x[1] + x[2]);
    grad_f[1] = x[0] * x[3];
    grad_f[2] = x[0] * x[3] + 1.0;
    grad_f[3] = x[0] * (x[0] + x[1] + x[2]);

    return true;
}

// Compute constraint function values
bool TutorialCpp_NLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number* g)
{
    assert(n == n_ && m == m_);

    g[0] = x[0] * x[1] * x[2] * x[3];
    g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];

    return true;
}

// Compute the Jacobian of constraints
bool TutorialCpp_NLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                                 Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    assert(n == n_ && m == m_ && nele_jac == nnz_jac_g_);

    // If values is nullptr, just provide the sparsity pattern
    if (values == nullptr) {
        // Return the structure of the Jacobian: 2 constraints × 4 vars = 8 non‐zeros
        Ipopt::Index idx = 0;
        // Constraint 0: g0 = x0*x1*x2*x3
        for (Ipopt::Index var = 0; var < n_; ++var) {
            iRow[idx] = 0; // row 0
            jCol[idx] = var; // cols 0,1,2,3
            ++idx;
        }
        // Constraint 1: g1 = x0^2 + x1^2 + x2^2 + x3^2
        for (Ipopt::Index var = 0; var < n_; ++var) {
            iRow[idx] = 1; // row 1
            jCol[idx] = var; // cols 0,1,2,3
            ++idx;
        }
    }
    else {
        // Return the values of the Jacobian
        Ipopt::Index idx = 0;
        // ∂g0/∂x_i = product of all x_j for j≠i
        values[idx++] = x[1] * x[2] * x[3]; // ∂g0/∂x0
        values[idx++] = x[0] * x[2] * x[3]; // ∂g0/∂x1
        values[idx++] = x[0] * x[1] * x[3]; // ∂g0/∂x2
        values[idx++] = x[0] * x[1] * x[2]; // ∂g0/∂x3

        // ∂g1/∂x_i = 2 * x_i
        values[idx++] = 2.0 * x[0]; // ∂g1/∂x0
        values[idx++] = 2.0 * x[1]; // ∂g1/∂x1
        values[idx++] = 2.0 * x[2]; // ∂g1/∂x2
        values[idx++] = 2.0 * x[3]; // ∂g1/∂x3
    }

    return true;
}

// Compute the Hessian of the Lagrangian (exact or approximation)
bool TutorialCpp_NLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number obj_factor,
                             Ipopt::Index m, const Ipopt::Number* lambda, bool /*new_lambda*/, Ipopt::Index nele_hess,
                             Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    assert(n == n_ && m == m_ && nele_hess == nnz_h_lag_);

    // If values is nullptr, just provide the sparsity pattern
    if (values == nullptr) {
        // 1) Return the structure (lower‐triangle) of the Hessian:
        //    rows 0→3, for j >= i
        Ipopt::Index idx = 0;
        // Row 0: (0,0),(0,1),(0,2),(0,3)
        iRow[idx]   = 0;
        jCol[idx++] = 0;
        iRow[idx]   = 0;
        jCol[idx++] = 1;
        iRow[idx]   = 0;
        jCol[idx++] = 2;
        iRow[idx]   = 0;
        jCol[idx++] = 3;
        // Row 1: (1,1),(1,2),(1,3)
        iRow[idx]   = 1;
        jCol[idx++] = 1;
        iRow[idx]   = 1;
        jCol[idx++] = 2;
        iRow[idx]   = 1;
        jCol[idx++] = 3;
        // Row 2: (2,2),(2,3)
        iRow[idx]   = 2;
        jCol[idx++] = 2;
        iRow[idx]   = 2;
        jCol[idx++] = 3;
        // Row 3: (3,3)
        iRow[idx]   = 3;
        jCol[idx++] = 3;
    }
    else {
        // 2) Return the values of each of those entries:
        // Precompute for readability
        Ipopt::Number l0 = lambda[0];
        Ipopt::Number l1 = lambda[1];

        Ipopt::Index idx = 0;

        values[idx++] = obj_factor * 2.0 * x[3] + l1 * 2.0; // (0,0)
        values[idx++] = obj_factor * x[3] + l0 * (x[2] * x[3]); // (0,1)
        values[idx++] = obj_factor * x[3] + l0 * (x[1] * x[3]); // (0,2)
        values[idx++] = obj_factor * (2.0 * x[0] + x[1] + x[2]) + l0 * (x[1] * x[2]); // (0,3)

        values[idx++] = l1 * 2.0; // (1,1)
        values[idx++] = l0 * (x[0] * x[3]); // (1,2)
        values[idx++] = obj_factor * x[0] + l0 * (x[0] * x[2]); // (1,3)

        values[idx++] = l1 * 2.0; // (2,2)
        values[idx++] = obj_factor * x[0] + l0 * (x[0] * x[1]); // (2,3)

        values[idx++] = l1 * 2.0; // (3,3)
    }

    return true;
}

// Finalize the optimization solution: display results and solver information.
// This function can also be adapted to save the final solution, variable values,
// and constraint evaluations to a file for further analysis or record-keeping.
void TutorialCpp_NLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                        const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                        const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                        const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* /*ip_cq*/
)
{
    std::cout << "IPOPT with Exact Derivative\n";

    // Print the final objective value
    std::cout << "Final Objective Value: " << obj_value << "\n";

    // Print the final solution for x
    std::cout << "Optimal Solution (x): \n";
    for (Ipopt::Index i = 0; i < n; ++i) {
        std::cout << "x[" << i << "] = " << x[i] << "\n";
    }

    // Optionally, print the constraints values and multipliers
    std::cout << "Constraints (g): \n";
    for (Ipopt::Index i = 0; i < m; ++i) {
        std::cout << "g[" << i << "] = " << g[i] << "\n";
    }

    if (ip_data) {
        std::cout << "Number of iterations: " << ip_data->iter_count() << "\n";
    }
}
