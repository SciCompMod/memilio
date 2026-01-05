#pragma once

#include "IpTNLP.hpp"
#include "memilio/ad/ad.h"

#include <vector>

// ---------------
// TutorialCpp_NLP
// ---------------
// Example implementation of an IPOPT TNLP (Nonlinear Programming) problem in C++.
//
// This class defines the structure of a nonlinear optimization problem, including:
//   - Variable bounds
//   - Constraint bounds
//   - Objective function and its gradient
//   - Constraint functions and its Jacobian
//   - Hessian (second-order derivatives) or Hessian approximation
//
// The class overrides the key TNLP virtual functions from Ipopt::TNLP:
//
//   get_nlp_info        -> provide problem dimensions and sparsity structure
//   get_bounds_info     -> provide bounds for variables and constraints
//   get_starting_point  -> provide an initial guess for the solver
//   eval_f              -> compute the objective function
//   eval_grad_f         -> compute the gradient of the objective
//   eval_g              -> compute constraint functions
//   eval_jac_g          -> compute Jacobian of constraints
//   eval_h              -> compute Hessian of the Lagrangian (optional)
//   finalize_solution   -> handle solution output after optimization
//
// The constructor allows toggling whether to use exact second-order derivatives or an approximation via 'use_hessian_approximation'.

class TutorialCpp_NLP : public Ipopt::TNLP
{
public:
    TutorialCpp_NLP(bool use_hessian_approximation) noexcept;

    TutorialCpp_NLP(const TutorialCpp_NLP&)            = delete;
    TutorialCpp_NLP(TutorialCpp_NLP&&)                 = delete;
    TutorialCpp_NLP& operator=(const TutorialCpp_NLP&) = delete;
    TutorialCpp_NLP& operator=(TutorialCpp_NLP&&)      = delete;
    ~TutorialCpp_NLP() override                        = default;

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
};
