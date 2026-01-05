// clang-format off

#include "secirvvs_ipopt.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <chrono>

#include "memilio/utils/mioomp.h"

#include "settings.h"
#include "IpIpoptData.hpp"

Secirvvs_NLP::Secirvvs_NLP(const ProblemSettings& problem_settings, const SimulationSettings& simulation_settings)
    : problem_settings_(problem_settings)
    , simulation_settings_(simulation_settings)
    , model_double_(simulation_settings.numAgeGroups())
    , model_tangent_linear_(simulation_settings.numAgeGroups())
    , model_adjoint_(simulation_settings.numAgeGroups())
{
    int num_control_intervals = problem_settings_.numControlIntervals();
    int pc_resolution = problem_settings_.pcResolution();
    int num_intervals = problem_settings_.numIntervals();
    int num_controls = problem_settings_.numControls();
    int num_path_constraints = problem_settings_.numPathConstraints();
    int num_terminal_constraints = problem_settings_.numTerminalConstraints();
    assert(num_control_intervals * pc_resolution == num_intervals);

    int effective_pc_nodes = 1;
    switch (problem_settings_.pathConstraintMode()) {
        case PathConstraintMode::Individual:
            effective_pc_nodes = num_intervals;
            break;
        case PathConstraintMode::GlobalMax:
            effective_pc_nodes = 1;
            break;
        default:
            throw std::runtime_error("Unsupported path constraint mode");
    }

    n_ = num_control_intervals * num_controls;
    m_ = effective_pc_nodes * num_path_constraints + num_terminal_constraints;

    nnz_jac_g_ = n_ * m_;
    bool use_exact_hessian = false;
    nnz_h_lag_ = use_exact_hessian ? n_ * (n_ + 1) / 2 : 0;

    // Resize vectors
    x_l_.resize(n_); x_u_.resize(n_);
    g_l_.resize(m_); g_u_.resize(m_);

    // Fill x_l_ and x_u_ from control bounds
    for (int controlInterval = 0; controlInterval < num_control_intervals; controlInterval++) {
        for (int controlIndex = 0; controlIndex < num_controls; controlIndex++) {
            int idx = controlInterval * num_controls + controlIndex;
            const auto& bounds = std::get<1>(problem_settings_.controlBounds()[controlIndex]);
            x_l_[idx] = bounds.first;
            x_u_[idx] = bounds.second;
        }
    }

    for (int interval = 0; interval < effective_pc_nodes; interval++) {
        for (int constraintIndex = 0; constraintIndex < num_path_constraints; ++constraintIndex) {
            int idx = interval * num_path_constraints + constraintIndex;
            const auto& bounds = std::get<1>(problem_settings_.pathConstraints()[constraintIndex]);
            g_l_[idx] = bounds.first;
            g_u_[idx] = bounds.second;
        }
    }

    int terminal_offset = effective_pc_nodes * num_path_constraints;
    for (int constraintIndex = 0; constraintIndex < num_terminal_constraints; constraintIndex++) {
        int idx = terminal_offset + constraintIndex;
        const auto& bounds = std::get<1>(problem_settings_.terminalConstraints()[constraintIndex]);
        g_l_[idx] = bounds.first;
        g_u_[idx] = bounds.second;
    }

    if (!ad::ga1s<double>::global_tape) {
        ad::ga1s<double>::global_tape = tape_t::create();
    }
    tape_ = ad::ga1s<double>::global_tape;

    // plain-double model for objective & constraints in native double precision
    set_initial_values(model_double_, problem_settings_, simulation_settings_);
    // forward-mode (tangent-linear) AD model for Jacobian via GT1S
    set_initial_values(model_tangent_linear_, problem_settings_, simulation_settings_);
    // reverse-mode (adjoint) AD model for gradient via GA1S
    set_initial_values(model_adjoint_, problem_settings_, simulation_settings_);
}

Secirvvs_NLP::~Secirvvs_NLP()
{
    if (tape_) {
        tape_t::remove(tape_);
        tape_ = nullptr;
    }
}

bool Secirvvs_NLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                                Ipopt::TNLP::IndexStyleEnum& index_style)
{
    n           = n_;
    m           = m_;
    nnz_jac_g   = nnz_jac_g_;
    nnz_h_lag   = nnz_h_lag_;
    index_style = Ipopt::TNLP::C_STYLE;
    return true;
}

bool Secirvvs_NLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m,
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

bool Secirvvs_NLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
                                      Ipopt::Number* /*z_L*/, Ipopt::Number* /*z_U*/, Ipopt::Index m, bool init_lambda,
                                      Ipopt::Number* /*lambda*/
)
{
    assert(init_x && !init_z && !init_lambda);
    assert(n == n_ && m == m_);

    int num_control_intervals = problem_settings_.numControlIntervals();
    int num_controls = problem_settings_.numControls();

    for (int controlInterval = 0; controlInterval < num_control_intervals; controlInterval++) {
        for (int controlIndex = 0; controlIndex < num_controls; controlIndex++) {
            int idx = controlInterval * num_controls + controlIndex;
            const auto& controlTuple = problem_settings_.controlBounds()[controlIndex];
            x[idx] = std::get<2>(controlTuple);
        }
    }

    return true;
}

bool Secirvvs_NLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& obj_value)
{
    assert(n == n_);

    obj_value = objective_function(model_double_, problem_settings_, x, n);

    return true;
}

bool Secirvvs_NLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number* g)
{
    assert(n == n_ && m == m_);

    constraint_functions(model_double_, problem_settings_, x, n, g, m);

    return true;
}

bool Secirvvs_NLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    bool use_forward_mode = true;

    if (use_forward_mode) {
        // Use forward mode for Jacobian evaluation
        return eval_grad_f_forward(n, x, new_x, grad_f);
    }
    else {
        // Use reverse mode for Jacobian evaluation
        return eval_grad_f_reverse(n, x, new_x, grad_f);
    }
}

bool Secirvvs_NLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                              Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    assert(n == n_ && m == m_ && nele_jac == nnz_jac_g_);

    bool use_forward_mode = true;
    switch (problem_settings_.pathConstraintMode()) {
        case PathConstraintMode::GlobalMax:
            use_forward_mode = true;
            break;
        case PathConstraintMode::Individual:
            use_forward_mode = true;
            break;
        default:
            throw std::runtime_error("Unsupported path constraint mode");
    }

    if (use_forward_mode) {
        // Use forward mode for Jacobian evaluation
        return eval_jac_g_forward(n, x, new_x, m, nele_jac, iRow, jCol, values);
    }
    else {
        // Use reverse mode for Jacobian evaluation
        return eval_jac_g_reverse(n, x, new_x, m, nele_jac, iRow, jCol, values);
    }
}

bool Secirvvs_NLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number obj_factor,
                          Ipopt::Index m, const Ipopt::Number* lambda, bool /*new_lambda*/, Ipopt::Index nele_hess,
                          Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    assert(!use_hessian_approximation_);
    return false;
}

bool Secirvvs_NLP::eval_grad_f_forward(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f)
{
    assert(n == n_);

    using ad_t = ad::gt1s<double>::type;

    PRAGMA_OMP(parallel)
    {
        std::vector<ad_t> x_ad(n);
        for (Ipopt::Index i = 0; i < n; ++i) {
            x_ad[i]                 = x[i];
            ad::derivative(x_ad[i]) = 0.0;
        }

        PRAGMA_OMP(for)
        for (Ipopt::Index column = 0; column < n; ++column) {
            ad::derivative(x_ad[column]) = 1.0;
            ad_t obj_ad                  = objective_function(model_tangent_linear_, problem_settings_, x_ad.data(), n);
            grad_f[column]               = ad::derivative(obj_ad);
            ad::derivative(x_ad[column]) = 0.0;
        }
    }

    return true;
}

bool Secirvvs_NLP::eval_grad_f_reverse(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f)
{
    assert(n == n_);
    using ad_t = ad::ga1s<double>::type;

    tape_->reset();

    std::vector<ad_t> x_ad(n);

    for (Ipopt::Index i = 0; i < n; ++i) {
        ad::value(x_ad[i])      = x[i];
        ad::derivative(x_ad[i]) = 0.0;
        tape_->register_variable(x_ad[i]);
    }

    ad_t obj_ad = objective_function(model_adjoint_, problem_settings_, x_ad.data(), n);
    tape_->register_output_variable(obj_ad);
    ad::derivative(obj_ad) = 1.0;

    tape_->interpret_adjoint();

    for (Ipopt::Index i = 0; i < n; ++i) {
        grad_f[i] = ad::derivative(x_ad[i]);
    }

    return true;
}

bool Secirvvs_NLP::eval_jac_g_forward(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                                      Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol,
                                      Ipopt::Number* values)
{
    assert(n == n_ && m == m_ && nele_jac == nnz_jac_g_);

    if (values == nullptr) {
        Ipopt::Index idx = 0;
        for (Ipopt::Index column = 0; column < n; ++column) {
            for (Ipopt::Index row = 0; row < m; ++row) {
                iRow[idx] = row;
                jCol[idx] = column;
                ++idx;
            }
        }
        return true;
    }

    using ad_t = ad::gt1s<double>::type;

    PRAGMA_OMP(parallel)
    {
        std::vector<ad_t> x_ad(n);
        std::vector<ad_t> g_ad(m);

        for (Ipopt::Index i = 0; i < n; ++i) {
            x_ad[i]                 = x[i];
            ad::derivative(x_ad[i]) = 0.0;
        }

        PRAGMA_OMP(for)
        for (Ipopt::Index column = 0; column < n; ++column) {
            ad::derivative(x_ad[column]) = 1.0;
            constraint_functions(model_tangent_linear_, problem_settings_, x_ad.data(), n, g_ad.data(), m);

            for (Ipopt::Index row = 0; row < m; ++row) {
                values[column * m + row] = ad::derivative(g_ad[row]);
            }

            ad::derivative(x_ad[column]) = 0.0;
        }
    }

    // Critical check: Abort if gradient is too large
    const double gradient_threshold = 1e15;
    for (Ipopt::Index column = 0; column < n; ++column) {
        for (Ipopt::Index row = 0; row < m; ++row) {
            const double grad_val = values[column * m + row];
            if (std::abs(grad_val) > gradient_threshold) {
                mio::log_critical("Gradient explosion detected: |∂g_{}/∂x_{}| = {:e} at x[{}] = {:e}",
                    row, column, grad_val, column, x[column]
                );
                return false;
            }
        }
    }

    return true;
}

bool Secirvvs_NLP::eval_jac_g_reverse(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                                      Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol,
                                      Ipopt::Number* values)
{
    assert(n == n_ && m == m_ && nele_jac == nnz_jac_g_);

    if (values == nullptr) {
        Ipopt::Index idx = 0;
        for (Ipopt::Index row = 0; row < m; ++row) {
            for (Ipopt::Index column = 0; column < n; ++column) {
                iRow[idx] = row;
                jCol[idx] = column;
                ++idx;
            }
        }
        return true;
    }

    using ad_t = ad::ga1s<double>::type;

    tape_->reset();

    std::vector<ad_t> x_ad(n);
    for (Ipopt::Index column = 0; column < n; ++column) {
        ad::value(x_ad[column]) = x[column];
        ad::derivative(x_ad[column]) = 0.0;
        tape_->register_variable(x_ad[column]);
    }

    std::vector<ad_t> g_ad(m);
    constraint_functions(model_adjoint_, problem_settings_, x_ad.data(), n, g_ad.data(), m);

    for (Ipopt::Index row = 0; row < m; ++row) {
        tape_->register_output_variable(g_ad[row]);
    }

    for (Ipopt::Index row = 0; row < m; ++row) {
        tape_->zero_adjoints();
        ad::derivative(g_ad[row]) = 1.0;
        tape_->interpret_adjoint(); // Takes up most of the time
        for (Ipopt::Index column = 0; column < n; ++column) {
            values[row * n + column] = ad::derivative(x_ad[column]);
        }
    }

    // Critical check: Abort if gradient is too large
    const double gradient_threshold = 1e15;
    for (Ipopt::Index row = 0; row < m; ++row) {
        for (Ipopt::Index column = 0; column < n; ++column) {
            const double grad_val = values[row * n + column];
            if (std::abs(grad_val) > gradient_threshold) {
                mio::log_critical("Gradient explosion detected: |∂g_{}/∂x_{}| = {:e} at x[{}] = {:e}",
                    row, column, grad_val, column, x[column]
                );
                return false;
            }
        }
    }

    return true;
}

void Secirvvs_NLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                     const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                     const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                     const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    std::cout << "Final Objective Value: " << obj_value << "\n";

    if (ip_data) {
        std::cout << "Number of iterations: " << ip_data->iter_count() << "\n";
    }

    save_solution(model_double_, problem_settings_, n, x, z_L, z_U, m, g, lambda, obj_value);

    return;
}
