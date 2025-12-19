/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Maximilian Betz
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
* 
* The documentation of the Ipopt::TNLP member functions in Secirvvs_NLP
* is extracted from the Ipopt documentation
*/

// This program implements direct single shooting for the optimal control of
// nonpharmazeutical intervation in pandemic ordinary differential equation (ODE) models.
// The socioeconomic costs are to be minmized using an objective functional which is the
// weighted integral of the nonpharmazeutical interventions and testing costs.
// In direct single shooting (also known as control vector parameterization) the
// continuous control variables are paremeterized using a discretization such that
// continuous optimal contorl prblem is transcribed into a finite dimensional nonlinear
// program which is then solved by the interior point method IPOPT. Here, we use
// a piecewise constant discretization.

#include "simulate_optimal_control_model.h"

class Secirvvs_NLP : public Ipopt::TNLP
{
public:

    Secirvvs_NLP(std::string data_dir, std::string temp_dir, std::string output_file_name, double constraint_compartment_upper_bound, 
                 std::vector<std::string> control_names, double tmax = 56);
    Secirvvs_NLP(const Secirvvs_NLP&) = delete;
    Secirvvs_NLP(Secirvvs_NLP&&)      = delete;
    Secirvvs_NLP& operator=(const Secirvvs_NLP&) = delete;
    Secirvvs_NLP& operator=(Secirvvs_NLP&&) = delete;

    /** Method to request the initial information about the problem.
    *
    *  %Ipopt uses this information when allocating the arrays
    *  that it will later ask you to fill with values. Be careful in this
    *  method since incorrect values will cause memory bugs which may be very
    *  difficult to find.
    *
    *  @param n           (out) Storage for the number of variables \f$x\f$
    *  @param m           (out) Storage for the number of constraints \f$g(x)\f$
    *  @param nnz_jac_g   (out) Storage for the number of nonzero entries in the Jacobian
    *  @param nnz_h_lag   (out) Storage for the number of nonzero entries in the Hessian
    *  @param index_style (out) Storage for the index style,
    *                     the numbering style used for row/col entries in the sparse matrix format
    *                     (TNLP::C_STYLE: 0-based, TNLP::FORTRAN_STYLE: 1-based; see also \ref TRIPLET)
    */
    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      IndexStyleEnum& index_style) override;

    /** Method to request bounds on the variables and constraints.
    *
    *  @param n   (in) the number of variables \f$x\f$ in the problem
    *  @param x_l (out) the lower bounds \f$x^L\f$ for the variables \f$x\f$
    *  @param x_u (out) the upper bounds \f$x^U\f$ for the variables \f$x\f$
    *  @param m   (in) the number of constraints \f$g(x)\f$ in the problem
    *  @param g_l (out) the lower bounds \f$g^L\f$ for the constraints \f$g(x)\f$
    *  @param g_u (out) the upper bounds \f$g^U\f$ for the constraints \f$g(x)\f$
    *
    *  @return true if success, false otherwise.
    *
    * The values of `n` and `m` that were specified in TNLP::get_nlp_info are passed
    * here for debug checking. Setting a lower bound to a value less than or
    * equal to the value of the option \ref OPT_nlp_lower_bound_inf "nlp_lower_bound_inf"
    * will cause %Ipopt to assume no lower bound. Likewise, specifying the upper bound above or
    * equal to the value of the option \ref OPT_nlp_upper_bound_inf "nlp_upper_bound_inf"
    * will cause %Ipopt to assume no upper bound. These options are set to -10<sup>19</sup> and
    * 10<sup>19</sup>, respectively, by default, but may be modified by changing these
    * options.
    */
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l,
                         Ipopt::Number* g_u) override;

    /** Method to request the starting point before iterating.
    *
    *  @param n      (in) the number of variables \f$x\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param init_x (in) if true, this method must provide an initial value for \f$x\f$
    *  @param x      (out) the initial values for the primal variables \f$x\f$
    *  @param init_z (in) if true, this method must provide an initial value for the bound multipliers \f$z^L\f$ and \f$z^U\f$
    *  @param z_L    (out) the initial values for the bound multipliers \f$z^L\f$
    *  @param z_U    (out) the initial values for the bound multipliers \f$z^U\f$
    *  @param m      (in) the number of constraints \f$g(x)\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param init_lambda (in) if true, this method must provide an initial value for the constraint multipliers \f$\lambda\f$
    *  @param lambda (out) the initial values for the constraint multipliers, \f$\lambda\f$
    *
    *  @return true if success, false otherwise.
    */
    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
                            Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;

    /** Method to request the value of the objective function.
    *
    *  @param n     (in) the number of variables \f$x\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param x     (in) the values for the primal variables \f$x\f$ at which the objective function \f$f(x)\f$ is to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise.
    *                    This can be helpful when users have efficient implementations that calculate multiple outputs at once.
    *                    %Ipopt internally caches results from the TNLP and generally, this flag can be ignored.
    *  @param obj_value (out) storage for the value of the objective function \f$f(x)\f$
    *
    *  @return true if success, false otherwise.
    */
    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) override;

    /** Method to request the gradient of the objective function.
    *
    *  @param n     (in) the number of variables \f$x\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param x     (in) the values for the primal variables \f$x\f$ at which the gradient \f$\nabla f(x)\f$ is to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise; see also TNLP::eval_f
    *  @param grad_f (out) array to store values of the gradient of the objective function \f$\nabla f(x)\f$.
    *                      The gradient array is in the same order as the \f$x\f$ variables
    *                      (i.e., the gradient of the objective with respect to `x[2]` should be put in `grad_f[2]`).
    *
    *  @return true if success, false otherwise.
    */
    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) override;

    /** Method to request the constraint values.
    *
    *  @param n     (in) the number of variables \f$x\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param x     (in) the values for the primal variables \f$x\f$ at which the constraint functions \f$g(x)\f$ are to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise; see also TNLP::eval_f
    *  @param m     (in) the number of constraints \f$g(x)\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param g     (out) array to store constraint function values \f$g(x)\f$, do not add or subtract the bound values \f$g^L\f$ or \f$g^U\f$.
    *
    *  @return true if success, false otherwise.
    */
    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) override;

    /** Method to request either the sparsity structure or the values of the Jacobian of the constraints.
    *
    * The Jacobian is the matrix of derivatives where the derivative of
    * constraint function \f$g_i\f$ with respect to variable \f$x_j\f$ is placed in row
    * \f$i\f$ and column \f$j\f$.
    * See \ref TRIPLET for a discussion of the sparse matrix format used in this method.
    *
    *  @param n     (in) the number of variables \f$x\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param x     (in) first call: NULL; later calls: the values for the primal variables \f$x\f$ at which the constraint Jacobian \f$\nabla g(x)^T\f$ is to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise; see also TNLP::eval_f
    *  @param m     (in) the number of constraints \f$g(x)\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param nele_jac (in) the number of nonzero elements in the Jacobian; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param iRow  (out) first call: array of length nele_jac to store the row indices of entries in the Jacobian of the constraints; later calls: NULL
    *  @param jCol  (out) first call: array of length nele_jac to store the column indices of entries in the Jacobian of the constraints; later calls: NULL
    *  @param values (out) first call: NULL; later calls: array of length nele_jac to store the values of the entries in the Jacobian of the constraints
    *
    *  @return true if success, false otherwise.
    *  */
    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                            Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) override;

    /** Method to request either the sparsity structure or the values of the Hessian of the Lagrangian.
    *
    * The Hessian matrix that %Ipopt uses is
    * \f[ \sigma_f \nabla^2 f(x_k) + \sum_{i=1}^m\lambda_i\nabla^2 g_i(x_k) \f]
    * for the given values for \f$x\f$, \f$\sigma_f\f$, and \f$\lambda\f$.
    * See \ref TRIPLET for a discussion of the sparse matrix format used in this method.
    *
    *  @param n     (in) the number of variables \f$x\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param x     (in) first call: NULL; later calls: the values for the primal variables \f$x\f$ at which the Hessian is to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise; see also TNLP::eval_f
    *  @param obj_factor (in) factor \f$\sigma_f\f$ in front of the objective term in the Hessian
    *  @param m     (in) the number of constraints \f$g(x)\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param lambda (in) the values for the constraint multipliers \f$\lambda\f$ at which the Hessian is to be evaluated
    *  @param new_lambda (in) false if any evaluation method was previously called with the same values in lambda, true otherwise
    *  @param nele_hess (in) the number of nonzero elements in the Hessian; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param iRow  (out) first call: array of length nele_hess to store the row indices of entries in the Hessian; later calls: NULL
    *  @param jCol  (out) first call: array of length nele_hess to store the column indices of entries in the Hessian; later calls: NULL
    *  @param values (out) first call: NULL; later calls: array of length nele_hess to store the values of the entries in the Hessian
    *
    *  @return true if success, false otherwise.
    */
    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
                        const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                        Ipopt::Index* jCol, Ipopt::Number* values) override;

    /** This method is called when the algorithm has finished (successfully or not) so the TNLP can digest the outcome, e.g., store/write the solution, if any.
    *
    *  @param status @parblock (in) gives the status of the algorithm
    *   - SUCCESS: Algorithm terminated successfully at a locally optimal
    *     point, satisfying the convergence tolerances (can be specified
    *     by options).
    *   - MAXITER_EXCEEDED: Maximum number of iterations exceeded (can be specified by an option).
    *   - CPUTIME_EXCEEDED: Maximum number of CPU seconds exceeded (can be specified by an option).
    *   - STOP_AT_TINY_STEP: Algorithm proceeds with very little progress.
    *   - STOP_AT_ACCEPTABLE_POINT: Algorithm stopped at a point that was converged, not to "desired" tolerances, but to "acceptable" tolerances (see the acceptable-... options).
    *   - LOCAL_INFEASIBILITY: Algorithm converged to a point of local infeasibility. Problem may be infeasible.
    *   - USER_REQUESTED_STOP: The user call-back function TNLP::intermediate_callback returned false, i.e., the user code requested a premature termination of the optimization.
    *   - DIVERGING_ITERATES: It seems that the iterates diverge.
    *   - RESTORATION_FAILURE: Restoration phase failed, algorithm doesn't know how to proceed.
    *   - ERROR_IN_STEP_COMPUTATION: An unrecoverable error occurred while %Ipopt tried to compute the search direction.
    *   - INVALID_NUMBER_DETECTED: Algorithm received an invalid number (such as NaN or Inf) from the NLP; see also option check_derivatives_for_nan_inf).
    *   - INTERNAL_ERROR: An unknown internal error occurred.
    *   @endparblock
    *  @param n     (in) the number of variables \f$x\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param x     (in) the final values for the primal variables
    *  @param z_L   (in) the final values for the lower bound multipliers
    *  @param z_U   (in) the final values for the upper bound multipliers
    *  @param m     (in) the number of constraints \f$g(x)\f$ in the problem; it will have the same value that was specified in TNLP::get_nlp_info
    *  @param g     (in) the final values of the constraint functions
    *  @param lambda (in) the final values of the constraint multipliers
    *  @param obj_value (in) the final value of the objective function
    *  @param ip_data (in) provided for expert users
    *  @param ip_cq (in) provided for expert users
    */
    virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                   const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                   const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                   const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) override;

    /** Intermediate Callback method for the user.
    *
    * This method is called once per iteration (during the convergence check),
    * and can be used to obtain information about the optimization status while
    * %Ipopt solves the problem, and also to request a premature termination.
    *
    * The information provided by the entities in the argument list correspond
    * to what %Ipopt prints in the iteration summary (see also \ref OUTPUT),
    * except for inf_pr, which by default corresponds to the original problem
    * in the log but to the scaled internal problem in this callback.
    * Further information can be obtained from the ip_data and ip_cq objects.
    * The current iterate and violations of feasibility and optimality can be
    * accessed via the methods Ipopt::TNLP::get_curr_iterate() and
    * Ipopt::TNLP::get_curr_violations().
    * These methods translate values for the *internal representation* of
    * the problem from `ip_data` and `ip_cq` objects into the TNLP representation.
    *
    * @return If this method returns false, %Ipopt will terminate with the
    *   User_Requested_Stop status.
    *
    * It is not required to implement (overload) this method.
    * The default implementation always returns true.
    */
    bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
                               Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number DeathRate,
                               Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du,
                               Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
                               Ipopt::IpoptCalculatedQuantities* ip_cq) override;

    /**
     * @brief Secirvvs_NLP::eval_objective_constraints evaluates the objective and the constraints of the NLP
     * @param x optimization variables of the NLP
     * @param constraints are the constraints of the NLP
     * @param objective is the objectie of the NLP
     * @param model is the ode model
     */
    template <class FP>
    void eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective, const EBModel<FP>& model);

public:
    int getN()
    {
        return n_;
    }
    int getM()
    {
        return m_;
    }

private:
    const std::string output_file_path_;
    
    // Scenario definition
    const double t0_;
    const double tmax_;
    const int num_age_groups_;
    const double constraint_compartment_upper_bound_;

    const int controlInterval_;
    const int numControlIntervals_; // number of piecewise constants interval for controls (same for all control variables)
    const int pcresolution_; // the resultion of path constraints is by this factor higher than the control discretization

    const std::vector<std::string> control_names_;
    const int numControls_; // number of control variables
    const int numPathConstraints_; // number of path constraints
    const int numIntervals_ = pcresolution_ * numControlIntervals_; // number of integration intervals
    const int n_            = numControlIntervals_ * numControls_; // number of optimization variables in the NLP
    const int m_            = numPathConstraints_; // number of constraints in the NLP

    std::unique_ptr<EBModel<gt1s_type>> model_ad_gt1s_;
    std::unique_ptr<EBModel<internal_type>> model_double_;
};

Secirvvs_NLP::Secirvvs_NLP(std::string data_dir, std::string temp_dir, std::string output_file_name, double constraint_compartment_upper_bound, 
                           std::vector<std::string> control_names, double tmax)
    : output_file_path_(mio::path_join(temp_dir, output_file_name))
    , t0_(0)
    , tmax_(tmax)
    , num_age_groups_(6)
    , constraint_compartment_upper_bound_(constraint_compartment_upper_bound)
    , controlInterval_(7)
    , numControlIntervals_(((int)tmax_ - 1) / controlInterval_ + 1)
    , pcresolution_(7)
    , control_names_(control_names)
    , numControls_(control_names_.size())
    , numPathConstraints_(1)
{
    double num_age_groups = num_age_groups_;

    mio::osecirvvs::Model<gt1s_type> model_ad_gt1s(num_age_groups);
    mio::osecirvvs::Model<internal_type> model_double(num_age_groups);

    auto out = set_initial_values(model_double, model_ad_gt1s, data_dir, temp_dir, tmax_);

    model_ad_gt1s_ = std::make_unique<mio::osecirvvs::Model<gt1s_type>>(model_ad_gt1s);
    model_double_ = std::make_unique<mio::osecirvvs::Model<internal_type>>(model_double);
}

bool Secirvvs_NLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                             IndexStyleEnum& index_style)
{
    n = n_;
    m = m_;

    // in this example the jacobian is dense and contains 8 nonzeros
    nnz_jac_g = m_ * n_;

    // the Hessian is also dense and has 16 total nonzeros, but we
    // only need the lower left corner (since it is symmetric)
    nnz_h_lag = n_ * n_;

    // use the C style indexing (0-based)
    index_style = Ipopt::TNLP::C_STYLE;

    return true;
}

bool Secirvvs_NLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m,
                                Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    mio::unused(n, m);
    // controls order: 1.-5. Damping variable 
    for (int n_intervals = 0; n_intervals < numControlIntervals_; ++n_intervals) {
        for (int n_controls = 0; n_controls < numControls_; ++n_controls) {
            x_l[n_intervals + n_controls * numControlIntervals_] = 0.0 + mio::Limits<double>::zero_tolerance(); // lower bound of 
            x_u[n_intervals + n_controls * numControlIntervals_] = 1.0 - mio::Limits<double>::zero_tolerance(); // upper bound of 
        } 
    }

    // path constraints
    for (int i = 0; i < m_; ++i) {
        g_l[i] = 0.0 + mio::Limits<double>::zero_tolerance();
        g_u[i] = constraint_compartment_upper_bound_ - mio::Limits<double>::zero_tolerance();
    }
    return true;
}

bool Secirvvs_NLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
                                   Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda)
{
    mio::unused(init_x, init_z, z_L, z_U, m, init_lambda, lambda);
    assert(init_z == false);
    assert(init_lambda == false);

    for (int i = 0; i < n; ++i) {
        x[i] = 0.5;
    }
    return true;
}

bool Secirvvs_NLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
    mio::unused(new_x);
    std::vector<double> xx(getN());
    std::vector<double> constraints(getM());
    for (int i = 0; i < n; ++i)
        xx[i] = x[i];
    eval_objective_constraints(xx, constraints, obj_value, *model_double_);
    return true;
}

bool Secirvvs_NLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    mio::unused(new_x);
    using FP = gt1s_type;
    std::vector<FP> xx(getN());
    std::vector<FP> constraints(getM());
    FP objective;
    for (int i = 0; i < n; ++i)
        ad::value(xx[i]) = x[i];
    for (int i = 0; i < n; ++i) {
        ad::derivative(xx[i]) = 1.0;
        eval_objective_constraints(xx, constraints, objective, *model_ad_gt1s_);
        grad_f[i]             = ad::derivative(objective);
        ad::derivative(xx[i]) = 0.0;
    }
    return true;
}

bool Secirvvs_NLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
    mio::unused(new_x);
    std::vector<double> xx(getN());
    std::vector<double> constraints(getM());
    double obj_value = 0;
    for (int i = 0; i < n; ++i)
        xx[i] = x[i];
    eval_objective_constraints(xx, constraints, obj_value, *model_double_);
    for (int i = 0; i < m; ++i)
        g[i] = constraints[i];
    return true;
}

bool Secirvvs_NLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                           Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    mio::unused(new_x, nele_jac);

    if (values == nullptr) {
        int jac_index = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                iRow[jac_index] = j;
                jCol[jac_index] = i;
                ++jac_index;
            }
        }
    }
    else {
        using FP = gt1s_type;
        std::vector<FP> xx(getN());
        std::vector<FP> constraints(getM());
        FP objective;
        int jac_index = 0;
        for (int i = 0; i < n; ++i)
            ad::value(xx[i]) = x[i];
        for (int i = 0; i < n; ++i) {
            ad::derivative(xx[i]) = 1.0;
            eval_objective_constraints(xx, constraints, objective, *model_ad_gt1s_);
            for (int j = 0; j < m; ++j) {
                values[jac_index] = ad::derivative(constraints[j]);
                ++jac_index;
                // std::cout << values[jac_index] << "\n";
            }
            ad::derivative(xx[i]) = 0.0;
        }
    }
    return true;
}

bool Secirvvs_NLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
                       const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                       Ipopt::Index* jCol, Ipopt::Number* values)
{
    mio::unused(n, x, new_x, obj_factor, m, lambda, new_lambda, new_lambda, nele_hess, iRow, jCol, values);
    return true;
}

template <class FP>
void Secirvvs_NLP::eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective, const EBModel<FP>& model)
{
    using std::max;

    FP t0   = t0_; // initial time
    FP tmax = tmax_; // stop time
    FP dt   = 0.2; // hint for initial step size for the integrator
    auto num_age_groups = num_age_groups_;

    std::vector<FP> grid(numIntervals_ + 1);
    for (int i = 0; i < numIntervals_ + 1; ++i) {
        grid[i] = (tmax / numIntervals_) * i + (t0 / numIntervals_) * (numIntervals_ - i);
    }

    std::vector<FP> xx(n_);
    for (int i = 0; i < n_; ++i) {
        xx[i] = x[i];
    }

    // mio::osecirvvs::Simulation<FP> sim(model, t0, dt);
    mio::Simulation<FP, mio::osecirvvs::Model<FP>> sim(model, t0, dt);
    // sim.set_integrator(std::make_unique<mio::ExplicitStepperWrapper<FP, boost::numeric::odeint::runge_kutta_fehlberg78>>());
    auto out = set_control_values<FP>(sim.get_model(), t0, tmax, xx, numControlIntervals_, control_names_);
    
    int gridindex = 0;
    objective     = 0.0;
    for (int controlIndex = 0; controlIndex < numControlIntervals_; ++controlIndex) {

        FP obj_value = 0;
        // the objective is the weighted integral of the controls, since we use a piecewise constant dicretization
        // we just can add the values
        for (int controlVariableIndex = 0; controlVariableIndex < numControls_; controlVariableIndex++) {
            obj_value += xx[controlIndex + controlVariableIndex * numControlIntervals_];
        }
        objective += pcresolution_ * obj_value;
        
        for (int i = 0; i < pcresolution_; ++i, ++gridindex) {
            
            sim.advance(grid[gridindex + 1]);
            auto result = sim.get_result().get_last_value();
            
            FP current_constraint = 0;
            for (int agegroup = 0; agegroup < num_age_groups; agegroup++)
            {
                auto age_group_offset = agegroup * (int)mio::osecirvvs::InfectionState::Count;

                current_constraint += result[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive + age_group_offset] + result[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity + age_group_offset] + result[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity + age_group_offset] + 
                    result[(int)mio::osecirvvs::InfectionState::InfectedSymptomsNaive + age_group_offset] + result[(int)mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity + age_group_offset] + result[(int)mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity + age_group_offset];
                
            }

            constraints[0] = max(constraints[0], current_constraint);
        }
    }
    
    return;
}

void Secirvvs_NLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                    const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                    const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                    const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    mio::unused(status, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);

    std::vector<double> grid(numIntervals_ + 1);
    for (int i = 0; i < numIntervals_ + 1; ++i) {
        grid[i] = (tmax_ / numIntervals_) * i + (t0_ / numIntervals_) * (numIntervals_ - i);
    }

    std::vector<double> xx(n);
    for (int i = 0; i < n; ++i)
        xx[i] = double(x[i]);
    
    auto round_to_decimal = [](double value, int decimal_places) {
        const double multiplier = std::pow(10.0, decimal_places);
        return std::round(value * multiplier) / multiplier;
    };

    // open file for parameter output
    mio::TimeSeries<double> resultControlVariables = mio::TimeSeries<double>(numControls_);
    Eigen::VectorX<double> currentControlVariables(numControls_); 

    int gridindex = 0;
    for (int controlIndex = 0; controlIndex < numControlIntervals_; ++controlIndex) {

        for (int controlVariableIndex = 0; controlVariableIndex < numControls_; controlVariableIndex++) {
            currentControlVariables[controlVariableIndex] = std::abs(round_to_decimal(xx[controlIndex + controlVariableIndex * numControlIntervals_], 6));
        }
        resultControlVariables.add_time_point(grid[gridindex], currentControlVariables);

        for (int i = 0; i < pcresolution_; ++i) {
            ++gridindex;
        }
    }
    resultControlVariables.add_time_point(grid[gridindex], currentControlVariables);

    auto out = resultControlVariables.export_csv(output_file_path_, control_names_, ' ', 6);

    return;
}

bool Secirvvs_NLP::intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
                                      Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number DeathRate,
                                      Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du,
                                      Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
                                      Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    mio::unused(mode, iter, obj_value, inf_pr, inf_du, DeathRate, d_norm, regularization_size, alpha_du, alpha_pr,
                ls_trials, ip_data, ip_cq);
    return true;
}

// Changes name to the currently used ones in optimal control. 
// Maybe should be same as the ones in Esid in the future.
std::string normalize_control_name(const std::string& name)
{
    if (name == "School closure")
        return "SchoolClosure";
    if (name == "Remote work")
        return "RemoteWork";
    if (name == "Face masks & social distancing School")
        return "PhysicalDistancingSchool";
    if (name == "Face masks & social distancing Work")
        return "PhysicalDistancingWork";
    if (name == "Face masks & social distancing Other")
        return "PhysicalDistancingOther";

    return name; // fallback
}

std::vector<std::string> get_linked_intervention_names(std::string scenario_data_run_file, std::string intervention_list_file)
{
    // Extract control names
    auto scenario_data_run = mio::read_json(scenario_data_run_file).value();
    auto intervention_list = mio::read_json(intervention_list_file).value();
    // Extract linked interventions
    auto linked_interventions = scenario_data_run["linkedInterventions"];

    std::vector<std::string> control_names;

    for (auto linked_intervention : scenario_data_run["linkedInterventions"])
    {
        std::string intervention_id = linked_intervention["interventionId"].asString();
        bool found = false;

        for (auto entry : intervention_list)
        {
            if (entry["id"].asString() == intervention_id)
            {
                control_names.push_back(
                    normalize_control_name(entry["name"].asString()));
                found = true;
                break;
            }
        }

        if (!found)
        {
            std::cerr << "Intervention " << intervention_id << " not found in intervention list\n";
        }
    }

    return control_names;
}

int main(int argc, char** argv)
{
    //switch of logging for mio
    mio::set_log_level(mio::LogLevel::err);

    auto cli_parameters = mio::cli::ParameterSetBuilder()
                            .add<"DataDirectory">(std::string(mio::path_join(mio::base_dir(), "data")))
                            .add<"TempDirectory">(std::string(mio::path_join(mio::base_dir(), "data")))
                            .add<"OutputFileName">(std::string{"OptimalControl_output.txt"})
                            .add<"ConstraintInfectedCases">(1e5)
                            .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, cli_parameters, {"DataDirectory"});
    if (!cli_result) {
        std::cout << cli_result.error().message();  
        return cli_result.error().code().value();  
    }

    // Extract linked interventions to use as controls
    std::vector<std::string> control_names = 
            get_linked_intervention_names(mio::path_join(cli_parameters.get<"TempDirectory">(), "scenario_data_run.json"), 
                                          mio::path_join(cli_parameters.get<"TempDirectory">(), "intervention_list.json"));;
    // std::vector<std::string> control_names = {"SchoolClosure", "RemoteWork", "PhysicalDistancingSchool", "PhysicalDistancingWork", "PhysicalDistancingOther"};

    // Create a new instance of your nlp
    std::cout << "*** Initialize Model!" << std::endl;
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new Secirvvs_NLP(cli_parameters.get<"DataDirectory">(), cli_parameters.get<"TempDirectory">(), cli_parameters.get<"OutputFileName">(), cli_parameters.get<"ConstraintInfectedCases">(), control_names);

    // Create a new instance of IpoptApplication
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // Change some options
    app->Options()->SetNumericValue("tol", 1e-5);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("limited_memory_update_type", "bfgs");
    app->Options()->SetIntegerValue("max_iter", 500);

    // Initialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int)status;
    }

    // Ask Ipopt to solve the problem
    std::cout << "*** Start Optimization!" << std::endl;
    status = app->OptimizeTNLP(mynlp);

    if (status == Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else {
        std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    }
    
    return (int)status;
}
