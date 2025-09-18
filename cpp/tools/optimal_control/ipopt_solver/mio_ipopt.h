#pragma once

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
#include <random>
#include <memory>

#include "IpTNLP.hpp"
#include "IpIpoptData.hpp"
#include "ad/ad.hpp"

#include "tools/optimal_control/ipopt_solver/objective_function.h"
#include "tools/optimal_control/ipopt_solver/constraint_function.h"
#include "tools/optimal_control/ipopt_solver/save_solution.h"

#include "models/ode_secirvvs/model.h"
#include "memilio/utils/mioomp.h"

template <class OptimizationSettings>
class MIO_NLP : public Ipopt::TNLP
{
public:
    MIO_NLP(const OptimizationSettings& settings);

    MIO_NLP(const MIO_NLP&)            = delete;
    MIO_NLP(MIO_NLP&&)                 = delete;
    MIO_NLP& operator=(const MIO_NLP&) = delete;
    MIO_NLP& operator=(MIO_NLP&&)      = delete;

    ~MIO_NLP() override;

    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      Ipopt::TNLP::IndexStyleEnum& index_style) override;

    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l,
                         Ipopt::Number* g_u) override;

    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
                            Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;

    bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) override;
    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) override;

    bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) override;
    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                    Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) override;

    bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
                const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                Ipopt::Index* jCol, Ipopt::Number* values) override;

    bool eval_grad_f_forward(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);
    bool eval_grad_f_reverse(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);

    bool eval_jac_g_forward(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                            Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values);
    bool eval_jac_g_reverse(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                            Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values);

    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L,
                           const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g,
                           const Ipopt::Number* lambda, Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
                           Ipopt::IpoptCalculatedQuantities* ip_cq) override;

private:
    const OptimizationSettings& m_settings;

    // plain-double model for objective & constraints in native double precision
    std::unique_ptr<OptimizationSettings::ModelTemplate<double>> m_model_double;
    // forward-mode (tangent-linear) AD model for Jacobian via GT1S
    std::unique_ptr<OptimizationSettings::ModelTemplate<ad::gt1s<double>::type>> m_model_tangent_linear;
    // reverse-mode (adjoint) AD model for gradient via GA1S
    std::unique_ptr<OptimizationSettings::ModelTemplate<ad::ga1s<double>::type>> m_model_adjoint;

    Ipopt::Index m_n;
    Ipopt::Index m_m;
    Ipopt::Index m_nnz_jac_g;
    Ipopt::Index m_nnz_h_lag;

    std::vector<Ipopt::Number> m_x_l;
    std::vector<Ipopt::Number> m_x_u;
    std::vector<Ipopt::Number> m_g_l;
    std::vector<Ipopt::Number> m_g_u;

    using tape_t   = ad::ga1s<double>::tape_t;
    tape_t* m_tape = nullptr;
};

template <class OptimizationSettings>
MIO_NLP<OptimizationSettings>::MIO_NLP(const OptimizationSettings& settings)
    : m_settings(settings)
    // plain-double model for objective & constraints in native double precision
    , m_model_double(std::make_unique<OptimizationSettings::ModelTemplate<double>>(settings.optimization_model().create_model<double>()))
    // forward-mode (tangent-linear) AD model for Jacobian via GT1S
    , m_model_tangent_linear(std::make_unique<OptimizationSettings::ModelTemplate<<ad::gt1s<double>>>(settings.optimization_model().create_model<ad::gt1s<double>::type>()))
    // reverse-mode (adjoint) AD model for gradient via GA1S
    , m_model_adjoint(std::make_unique<OptimizationSettings::ModelTemplate<<ad::ga1s<double>>>(settings.optimization_model().create_model<ad::ga1s<double>::type>()))
{
    size_t num_control_intervals    = m_settings.num_control_intervals();
    size_t pc_resolution            = m_settings.pc_resolution();
    size_t num_intervals            = m_settings.num_intervals();
    size_t num_control_parameters   = m_settings.num_control_parameters();
    size_t num_path_constraints     = m_settings.num_path_constraints();
    size_t num_terminal_constraints = m_settings.num_terminal_constraints();
    assert(num_control_intervals * pc_resolution == num_intervals);

    m_n = num_control_intervals * num_control_parameters;
    m_m = num_path_constraints + num_terminal_constraints;

    m_nnz_jac_g            = m_n * m_m;
    bool use_exact_hessian = false;
    m_nnz_h_lag            = use_exact_hessian ? m_n * (m_n + 1) / 2 : 0;

    // Resize vectors
    m_x_l.resize(m_n);
    m_x_u.resize(m_n);
    m_g_l.resize(m_m);
    m_g_u.resize(m_m);

    // Fill m_x_l and m_x_u from control parameters
    for (size_t control_interval = 0; control_interval < num_control_intervals; control_interval++) {
        for (size_t control_index = 0; control_index < num_control_parameters; control_index++) {
            size_t idx = control_index + control_interval * num_control_parameters;
            m_x_l[idx] = m_settings.control_parameters()[control_index].min();
            m_x_u[idx] = m_settings.control_parameters()[control_index].max();
        }
    }
    // Fill m_g_l and m_g_u from path- and terminal constraints
    size_t constraint_index = 0;
    for (size_t path_constraint_index = 0; path_constraint_index < num_path_constraints; path_constraint_index++) {
        m_g_l[constraint_index] = m_settings.path_constraints()[path_constraint_index].min();
        m_g_u[constraint_index] = m_settings.path_constraints()[path_constraint_index].max();
        constraint_index++;
    }
    for (size_t terminal_constraint_index = 0; terminal_constraint_index < num_terminal_constraints;
         terminal_constraint_index++) {
        m_g_l[constraint_index] = m_settings.terminal_constraints()[terminal_constraint_index].min();
        m_g_u[constraint_index] = m_settings.terminal_constraints()[terminal_constraint_index].max();
        constraint_index++;
    }

    if (m_settings.ad_eval_f() == ADType::Reverse || m_settings.ad_eval_jac() == ADType::Reverse) {
        if (!ad::ga1s<double>::global_tape) {
            ad::ga1s<double>::global_tape = tape_t::create();
        }
        m_tape = ad::ga1s<double>::global_tape;
    }
}

template <class OptimizationSettings>
MIO_NLP<OptimizationSettings>::~MIO_NLP()
{
    if (m_tape) {
        tape_t::remove(m_tape);
        m_tape = nullptr;
    }
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                                Ipopt::TNLP::IndexStyleEnum& index_style)
{
    n           = m_n;
    m           = m_m;
    nnz_jac_g   = m_nnz_jac_g;
    nnz_h_lag   = m_nnz_h_lag;
    index_style = Ipopt::TNLP::C_STYLE;
    return true;
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m,
                                   Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    assert(n == m_n && m == m_m);
    for (Ipopt::Index i = 0; i < n; i++) {
        x_l[i] = m_x_l[i];
        x_u[i] = m_x_u[i];
    }

    for (Ipopt::Index i = 0; i < m; i++) {
        g_l[i] = m_g_l[i];
        g_u[i] = m_g_u[i];
    }

    return true;
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
                                      Ipopt::Number* /*z_L*/, Ipopt::Number* /*z_U*/, Ipopt::Index m, bool init_lambda,
                                      Ipopt::Number* /*lambda*/
)
{
    assert(init_x && !init_z && !init_lambda);
    assert(n == m_n && m == m_m);

    size_t num_control_intervals  = m_settings.num_control_intervals();
    size_t num_control_parameters = m_settings.num_control_parameters();

    std::random_device rd;
    std::mt19937 gen(rd());

    for (size_t control_interval = 0; control_interval < num_control_intervals; control_interval++) {
        for (size_t control_index = 0; control_index < num_control_parameters; control_index++) {
            size_t idx = control_index + control_interval * num_control_parameters;

            double min_val = m_settings.control_parameters()[control_index].min();
            double max_val = m_settings.control_parameters()[control_index].max();

            if (m_settings.random_start()) {
                std::uniform_real_distribution<> dis(min_val, max_val);
                x[idx] = dis(gen);
            }
            else {
                x[idx] = 0.5 * (min_val + max_val);
            }
        }
    }

    return true;
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& obj_value)
{
    assert(n == m_n);

    obj_value = objective_function(m_settings, *m_model_double, x, n);

    return true;
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number* g)
{
    assert(n == m_n && m == m_m);

    constraint_function(m_settings, *m_model_double, x, n, g, m);

    return true;
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    assert(n == m_n);

    if (m_settings.ad_eval_f() == ADType::Forward) {
        // Use forward mode for Jacobian evaluation
        return eval_grad_f_forward(n, x, new_x, grad_f);
    }
    else {
        // Use reverse mode for Jacobian evaluation
        return eval_grad_f_reverse(n, x, new_x, grad_f);
    }
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                              Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    assert(n == m_n && m == m_m && nele_jac == m_nnz_jac_g);

    if (m_settings.ad_eval_jac() == ADType::Forward) {
        // Use forward mode for Jacobian evaluation
        return eval_jac_g_forward(n, x, new_x, m, nele_jac, iRow, jCol, values);
    }
    else {
        // Use reverse mode for Jacobian evaluation
        return eval_jac_g_reverse(n, x, new_x, m, nele_jac, iRow, jCol, values);
    }
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number obj_factor,
                          Ipopt::Index m, const Ipopt::Number* lambda, bool /*new_lambda*/, Ipopt::Index nele_hess,
                          Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    return false;
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::eval_grad_f_forward(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f)
{
    assert(n == m_n);

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
            ad_t obj_ad                  = objective_function(m_settings, *m_model_tangent_linear, x_ad.data(), n);
            grad_f[column]               = ad::derivative(obj_ad);
            ad::derivative(x_ad[column]) = 0.0;
        }
    }

    return true;
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::eval_grad_f_reverse(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f)
{
    assert(n == m_n);
    using ad_t = ad::ga1s<double>::type;

    m_tape->reset();

    std::vector<ad_t> x_ad(n);

    for (Ipopt::Index i = 0; i < n; ++i) {
        ad::value(x_ad[i])      = x[i];
        ad::derivative(x_ad[i]) = 0.0;
        m_tape->register_variable(x_ad[i]);
    }

    ad_t obj_ad = objective_function(m_settings, *m_model_adjoint, x_ad.data(), n);
    m_tape->register_output_variable(obj_ad);
    ad::derivative(obj_ad) = 1.0;

    m_tape->interpret_adjoint();

    for (Ipopt::Index i = 0; i < n; ++i) {
        grad_f[i] = ad::derivative(x_ad[i]);
    }

    return true;
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::eval_jac_g_forward(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                                      Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol,
                                      Ipopt::Number* values)
{
    assert(n == m_n && m == m_m && nele_jac == m_nnz_jac_g);

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
            constraint_function(m_settings, *m_model_tangent_linear, x_ad.data(), n, g_ad.data(), m);

            for (Ipopt::Index row = 0; row < m; ++row) {
                values[column * m + row] = ad::derivative(g_ad[row]);
            }

            ad::derivative(x_ad[column]) = 0.0;
        }
    }

    // ⚠️ Critical check: Abort if gradient is too large
    const double gradient_threshold = 1e15;
    for (Ipopt::Index column = 0; column < n; ++column) {
        for (Ipopt::Index row = 0; row < m; ++row) {
            const double grad_val = values[column * m + row];
            if (std::abs(grad_val) > gradient_threshold) {
                mio::log_critical("Gradient explosion detected: |∂g_{}/∂x_{}| = {:e} at x[{}] = {:e}", row, column,
                                  grad_val, column, x[column]);
                return false;
            }
        }
    }

    return true;
}

template <class OptimizationSettings>
bool MIO_NLP<OptimizationSettings>::eval_jac_g_reverse(Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                                      Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol,
                                      Ipopt::Number* values)
{
    assert(n == m_n && m == m_m && nele_jac == m_nnz_jac_g);

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

    m_tape->reset();

    std::vector<ad_t> x_ad(n);
    for (Ipopt::Index column = 0; column < n; ++column) {
        ad::value(x_ad[column])      = x[column];
        ad::derivative(x_ad[column]) = 0.0;
        m_tape->register_variable(x_ad[column]);
    }

    std::vector<ad_t> g_ad(m);
    constraint_function(m_settings, *m_model_adjoint, x_ad.data(), n, g_ad.data(), m);

    for (Ipopt::Index row = 0; row < m; ++row) {
        m_tape->register_output_variable(g_ad[row]);
    }

    for (Ipopt::Index row = 0; row < m; ++row) {
        m_tape->zero_adjoints();
        ad::derivative(g_ad[row]) = 1.0;
        m_tape->interpret_adjoint(); // Takes up most of the time
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
                mio::log_critical("Gradient explosion detected: |∂g_{}/∂x_{}| = {:e} at x[{}] = {:e}", row, column,
                                  grad_val, column, x[column]);
                return false;
            }
        }
    }

    return true;
}

template <class OptimizationSettings>
void MIO_NLP<OptimizationSettings>::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                     const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                     const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                     const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    std::cout << "Final Objective Value: " << obj_value << "\n";

    if (ip_data) {
        std::cout << "Number of iterations: " << ip_data->iter_count() << "\n";
    }

    save_solution(m_settings, *m_model_double, n, x, z_L, z_U, m, g, lambda, obj_value);

    return;
}
