#pragma once

#include <vector>

#include "IpTNLP.hpp"
#include "ad/ad.hpp"

#include "../optimization_settings/secirvvs_optimization.h"
#include "models/ode_secirvvs/model.h"

class Secirvvs_NLP : public Ipopt::TNLP
{
public:
    Secirvvs_NLP(const SecirvvsOptimization& settings);

    Secirvvs_NLP(const Secirvvs_NLP&)            = delete;
    Secirvvs_NLP(Secirvvs_NLP&&)                 = delete;
    Secirvvs_NLP& operator=(const Secirvvs_NLP&) = delete;
    Secirvvs_NLP& operator=(Secirvvs_NLP&&)      = delete;

    ~Secirvvs_NLP() override;

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
    const SecirvvsOptimization& m_settings;

    // plain-double model for objective & constraints in native double precision
    mio::osecirvvs::Model<double> m_model_double;
    // forward-mode (tangent-linear) AD model for Jacobian via GT1S
    mio::osecirvvs::Model<ad::gt1s<double>::type> m_model_tangent_linear;
    // reverse-mode (adjoint) AD model for gradient via GA1S
    mio::osecirvvs::Model<ad::ga1s<double>::type> m_model_adjoint;

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