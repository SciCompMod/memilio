#pragma once

#include "IpTNLP.hpp"
#include "memilio/ad/ad.h"

#include <vector>

class TutorialCpp_NLP : public Ipopt::TNLP
{
public:
    TutorialCpp_NLP(bool use_hessian_approximation) noexcept;

    TutorialCpp_NLP(const TutorialCpp_NLP&)            = delete;
    TutorialCpp_NLP(TutorialCpp_NLP&&)                 = delete;
    TutorialCpp_NLP& operator=(const TutorialCpp_NLP&) = delete;
    TutorialCpp_NLP& operator=(TutorialCpp_NLP&&)      = delete;

    ~TutorialCpp_NLP() override;

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

    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L,
                           const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g,
                           const Ipopt::Number* lambda, Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
                           Ipopt::IpoptCalculatedQuantities* ip_cq) override;

private:
    // Problem dimensions
    Ipopt::Index n_; // Number of decision variables
    Ipopt::Index m_; // Number of constraints
    Ipopt::Index nnz_jac_g_; // Number of nonzeros in the constraint Jacobian
    Ipopt::Index nnz_h_lag_; // Number of nonzeros in the Hessian of the Lagrangian

    bool use_hessian_approximation_; // False: use limited-memory quasi-Newton approximation

    // Variable and constraint bounds
    std::vector<Ipopt::Number> x_l_; // Lower bounds for variables
    std::vector<Ipopt::Number> x_u_; // Upper bounds for variables
    std::vector<Ipopt::Number> g_l_; // Lower bounds for constraints
    std::vector<Ipopt::Number> g_u_; // Upper bounds for constraints

    // Reverse-mode AD tape for computing derivatives
    using tape_t = ad::ga1s<double>::tape_t; // Type of the AD tape
    tape_t* tape_; // Pointer to the AD tape used for recording operations
};
