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

template <typename T>
T objective_function(const T* x, std::size_t n)
{
    assert(n >= 4);
    return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
}

template <typename T>
void constraint_functions(const T* x, std::size_t n, T* g, std::size_t m)
{
    assert(n >= 4);
    assert(m == 2);
    g[0] = x[0] * x[1] * x[2] * x[3];
    g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
}

constexpr Ipopt::Number INF = 1e19;

TutorialCpp_NLP::TutorialCpp_NLP(bool use_hessian_approximation) noexcept
{
    n_                         = 4;
    m_                         = 2;
    nnz_jac_g_                 = n_ * m_;
    nnz_h_lag_                 = use_hessian_approximation ? n_ * (n_ + 1) / 2 : 0;
    use_hessian_approximation_ = use_hessian_approximation;

    x_l_.assign(n_, 1.0);
    x_u_.assign(n_, 5.0);

    g_l_.resize(m_);
    g_u_.resize(m_);
    g_l_[0] = 25.0;
    g_u_[0] = INF;
    g_l_[1] = 40.0;
    g_u_[1] = 40.0;

    if (!ad::ga1s<double>::global_tape) {
        ad::ga1s<double>::global_tape = tape_t::create();
    }
    tape_ = ad::ga1s<double>::global_tape;
}

TutorialCpp_NLP::~TutorialCpp_NLP()
{
    if (tape_) {
        tape_t::remove(tape_);
        tape_ = nullptr;
    }
}

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

bool TutorialCpp_NLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
                                         Ipopt::Number* /*z_L*/, Ipopt::Number* /*z_U*/, Ipopt::Index m,
                                         bool init_lambda, Ipopt::Number* /*lambda*/
)
{
    assert(init_x && !init_z && !init_lambda);
    assert(n == n_ && m == m_);

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

bool TutorialCpp_NLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& obj_value)
{
    assert(n == n_);
    obj_value = objective_function(x, static_cast<std::size_t>(n));
    return true;
}

bool TutorialCpp_NLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f)
{
    assert(n == n_);
    using ad_t = ad::ga1s<double>::type; // AD type for reverse-mode differentiation

    tape_->reset();

    std::vector<ad_t> x_ad(n);

    for (Ipopt::Index i = 0; i < n; ++i) {
        ad::value(x_ad[i])      = x[i];
        ad::derivative(x_ad[i]) = 0.0;
        tape_->register_variable(x_ad[i]);
    }

    ad_t obj_ad = objective_function(x_ad.data(), static_cast<std::size_t>(n));
    tape_->register_output_variable(obj_ad);
    ad::derivative(obj_ad) = 1.0;

    tape_->interpret_adjoint(); // Propagate derivatives backward through the tape

    for (Ipopt::Index i = 0; i < n; ++i) {
        grad_f[i] = ad::derivative(x_ad[i]);
    }

    return true;
}

bool TutorialCpp_NLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number* g)
{
    assert(n == n_ && m == m_);

    constraint_functions(x, n, g, m);

    return true;
}

bool TutorialCpp_NLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                                 Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    assert(n == n_ && m == m_ && nele_jac == nnz_jac_g_);

    // If values is nullptr, just provide the sparsity pattern
    if (values == nullptr) {
        Ipopt::Index idx = 0;
        for (Ipopt::Index row = 0; row < m; ++row) {
            for (Ipopt::Index column = 0; column < n; ++column) {
                iRow[idx] = row;
                jCol[idx] = column;
                ++idx;
            }
        }
    }
    else {
        using ad_t = ad::ga1s<double>::type; // AD type for reverse-mode differentiation

        tape_->reset();

        std::vector<ad_t> x_ad(n);
        for (Ipopt::Index j = 0; j < n; ++j) {
            ad::value(x_ad[j])      = x[j];
            ad::derivative(x_ad[j]) = 0.0;
            tape_->register_variable(x_ad[j]);
        }

        std::vector<ad_t> g_ad(m);
        constraint_functions(x_ad.data(), n, g_ad.data(), m);

        for (Ipopt::Index i = 0; i < m; ++i) {
            tape_->register_output_variable(g_ad[i]);
        }

        for (Ipopt::Index i = 0; i < m; ++i) {
            tape_->zero_adjoints();
            ad::derivative(g_ad[i]) = 1.0;
            tape_->interpret_adjoint(); // Propagate derivatives backward through the tape
            for (Ipopt::Index j = 0; j < n; ++j) {
                values[i * n + j] = ad::derivative(x_ad[j]);
            }
        }
    }

    return true;
}

bool TutorialCpp_NLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number obj_factor,
                             Ipopt::Index m, const Ipopt::Number* lambda, bool /*new_lambda*/, Ipopt::Index nele_hess,
                             Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    assert(!use_hessian_approximation_);
    return false;
}

void TutorialCpp_NLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                        const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                        const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                        const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* /*ip_cq*/
)
{
    std::cout << "IPOPT with Reverse-Mode AD\n";
    std::cout << "(Allocating the tape may need additional time)\n";

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
