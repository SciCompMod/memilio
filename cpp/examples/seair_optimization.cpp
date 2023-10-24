/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
*/
#include "ad/ad.hpp"


#include "ode_seair/model.h"
#include "ode_seair/infection_state.h"
#include "ode_seair/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/time_series_to_file.h"
#include "IpTNLP.hpp" // IWYU pragma: keep
#include <fstream>




class Seair_NLP: public Ipopt::TNLP
{
public:
    static constexpr double N = 327167434;
    Seair_NLP() = default;
    Seair_NLP(const Seair_NLP&) = delete;
    Seair_NLP(Seair_NLP&&) = delete;
    Seair_NLP& operator=(const Seair_NLP&) = delete;
    Seair_NLP& operator=(Seair_NLP&&) = delete;
    ~Seair_NLP() = default;
    bool get_nlp_info(
        Ipopt::Index&          n,
        Ipopt::Index&          m,
        Ipopt::Index&          nnz_jac_g,
        Ipopt::Index&          nnz_h_lag,
        IndexStyleEnum& index_style
        ) override;
    bool get_bounds_info(
        Ipopt::Index   n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index   m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u
    ) override;

    bool get_starting_point(
        Ipopt::Index   n,
        bool    init_x,
        Ipopt::Number* x,
        bool    init_z,
        Ipopt::Number* z_L,
        Ipopt::Number* z_U,
        Ipopt::Index   m,
        bool    init_lambda,
        Ipopt::Number* lambda
    ) override;

    /** Method to return the objective value */
    virtual bool eval_f(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool          new_x,
        Ipopt::Number&       obj_value
        ) override;

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool          new_x,
        Ipopt::Number*       grad_f
        ) override;

    /** Method to return the constraint residuals */
    virtual bool eval_g(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool          new_x,
        Ipopt::Index         m,
        Ipopt::Number*       g
        ) override;

    /** Method to return:
    *   1) The structure of the jacobian (if "values" is NULL)
    *   2) The values of the jacobian (if "values" is not NULL)
    */
    virtual bool eval_jac_g(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool          new_x,
        Ipopt::Index         m,
        Ipopt::Index         nele_jac,
        Ipopt::Index*        iRow,
        Ipopt::Index*        jCol,
        Ipopt::Number*       values
        ) override;

    /** Method to return:
    *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
    *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
    */
    virtual bool eval_h(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool          new_x,
        Ipopt::Number        obj_factor,
        Ipopt::Index         m,
        const Ipopt::Number* lambda,
        bool          new_lambda,
        Ipopt::Index         nele_hess,
        Ipopt::Index*        iRow,
        Ipopt::Index*        jCol,
        Ipopt::Number*       values
        ) override;

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(
        Ipopt::SolverReturn               status,
        Ipopt::Index                      n,
        const Ipopt::Number*              x,
        const Ipopt::Number*              z_L,
        const Ipopt::Number*              z_U,
        Ipopt::Index                      m,
        const Ipopt::Number*              g,
        const Ipopt::Number*              lambda,
        Ipopt::Number                     obj_value,
        const Ipopt::IpoptData*           ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq
        ) override;
    //@}

    bool intermediate_callback(
        Ipopt::AlgorithmMode              mode,
        Ipopt::Index                      iter,
        Ipopt::Number                     obj_value,
        Ipopt::Number                     inf_pr,
        Ipopt::Number                     inf_du,
        Ipopt::Number                     mu,
        Ipopt::Number                     d_norm,
        Ipopt::Number                     regularization_size,
        Ipopt::Number                     alpha_du,
        Ipopt::Number                     alpha_pr,
        Ipopt::Index                      ls_trials,
        const Ipopt::IpoptData*           ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq
        ) override;


    template<typename FP=double>
    void eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective);

private:
    const int numControlIntervals_ = 20;
    const int numControls_ = 3;
    const int numPathConstraints_ = 1;
    const int pcresolution_ = 5;
    const int numIntervals_ = pcresolution_ * numControlIntervals_;
    const int n_ = numControlIntervals_ * numControlIntervals_;
    const int m_ = numIntervals_ * numPathConstraints_;


};


template<typename FP>
void set_initial_values(mio::oseair::Model<FP>& model) {
    const double N = 327167434;// total population of the United States

    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] = 0.9977558755803503;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}]   = 0.0003451395725394549;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}]   = 0.00037846880968213874;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}]  = (337072.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] = (17448.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Perished)}]   = (9619.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::ObjectiveFunction)}]   = 0.0;

}

template<typename FP>
void Seair_NLP::eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective){

    FP t0   = 0;
    FP tmax = 100;
    FP dt=0.2;
    std::vector<FP> grid(numIntervals_+1);
    for(int i = 0; i < numIntervals_+1; ++i) {
        grid[i] = (tmax/numIntervals_)*i +(t0/numIntervals_)*(numIntervals_-i);
    }
    mio::oseair::Model<FP> model;
    set_initial_values(model);
    for(int i = 0; i < numIntervals_; ++i) {
        auto seair1 = mio::simulate<mio::oseair::Model<FP>,FP>(grid[i], grid[i+1], dt, model);
    }

}

int main()
{

    return 0;

}

bool Seair_NLP::get_nlp_info(
    Ipopt::Index&          n,
    Ipopt::Index&          m,
    Ipopt::Index&          nnz_jac_g,
    Ipopt::Index&          nnz_h_lag,
    IndexStyleEnum& index_style
    )
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


bool Seair_NLP::get_bounds_info(
    Ipopt::Index   n,
    Ipopt::Number* x_l,
    Ipopt::Number* x_u,
    Ipopt::Index   m,
    Ipopt::Number* g_l,
    Ipopt::Number* g_u
    )
{
    // controls order: 1. alpha_a, 2. alpha_i, 3. kappa
    for(int i=0; i < numControls_; ++i) {
        x_l[i] = 0.05; // lower bound of alpha_a
        x_u[i] = 0.5;  // upper bound of alpha_a
        x_l[i + numControlIntervals_] = 0.01; // lower bound of alpha_i
        x_u[i + numControlIntervals_] = 0.3;  // upper bound of alpha_i
        x_l[i + 2 * numControlIntervals_] = 0.15; // lower bound of kappa
        x_u[i + 2 * numControlIntervals_] = 0.3;  // upper bound of kappa
    }

    // path constraints
    for(int i=0; i < m_; ++i) {
        g_l[i] = 0.0;
        g_u[i] = 1e6/N;
    }
    return true;
}

bool Seair_NLP::get_starting_point(
    Ipopt::Index   n,
    bool    init_x,
    Ipopt::Number* x,
    bool    init_z,
    Ipopt::Number* z_L,
    Ipopt::Number* z_U,
    Ipopt::Index   m,
    bool    init_lambda,
    Ipopt::Number* lambda
    )
{
    assert(init_z == false);
    assert(init_lambda == false);

    for(int i = 0; i < n; ++i) {
        x[i] = 0.2;
    }
    return true;
}

bool Seair_NLP::eval_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number &obj_value)
{
    return true;
}

bool Seair_NLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number *grad_f)
{
    return true;
}

bool Seair_NLP::eval_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Number *g)
{
    return true;
}

bool Seair_NLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
    return true;
}

bool Seair_NLP::eval_h(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number *lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
    return true;
}

void Seair_NLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number *x, const Ipopt::Number *z_L, const Ipopt::Number *z_U, Ipopt::Index m, const Ipopt::Number *g, const Ipopt::Number *lambda, Ipopt::Number obj_value, const Ipopt::IpoptData *ip_data, Ipopt::IpoptCalculatedQuantities *ip_cq)
{
    return;
}

bool Seair_NLP::intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData *ip_data, Ipopt::IpoptCalculatedQuantities *ip_cq)
{
    return true;
}
