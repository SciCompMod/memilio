/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Ralf Hannemann-Tamas
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
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include <fstream>

class Seair_NLP : public Ipopt::TNLP
{
public:
    static constexpr double N              = 327167434;
    Seair_NLP()                            = default;
    Seair_NLP(const Seair_NLP&)            = delete;
    Seair_NLP(Seair_NLP&&)                 = delete;
    Seair_NLP& operator=(const Seair_NLP&) = delete;
    Seair_NLP& operator=(Seair_NLP&&)      = delete;
    ~Seair_NLP()                           = default;
    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      IndexStyleEnum& index_style) override;
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l,
                         Ipopt::Number* g_u) override;

    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
                            Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;

    /** Method to return the objective value */
    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) override;

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) override;

    /** Method to return the constraint residuals */
    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) override;

    /** Method to return:
    *   1) The structure of the jacobian (if "values" is NULL)
    *   2) The values of the jacobian (if "values" is not NULL)
    */
    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                            Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) override;

    /** Method to return:
    *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
    *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
    */
    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
                        const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                        Ipopt::Index* jCol, Ipopt::Number* values) override;

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                   const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                   const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                   const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) override;
    //@}

    bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
                               Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm,
                               Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr,
                               Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
                               Ipopt::IpoptCalculatedQuantities* ip_cq) override;

    template <typename FP = double>
    void eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective);

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
    const int numControlIntervals_ = 20;
    const int numControls_         = 3;
    const int numPathConstraints_  = 1;
    const int pcresolution_        = 5;
    const int numIntervals_        = pcresolution_ * numControlIntervals_;
    const int n_                   = numControlIntervals_ * numControls_;
    const int m_                   = numIntervals_ * numPathConstraints_;
};

template <typename FP>
void set_initial_values(mio::oseair::Model<FP>& model)
{
    const double N = Seair_NLP::N; // total population of the United States

    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] =
        0.9977558755803503;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}] =
        0.0003451395725394549;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] =
        0.00037846880968213874;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}] =
        (337072.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] =
        (17448.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Perished)}] = (9619.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::ObjectiveFunction)}] = 0.0;
}

template <typename FP>
void Seair_NLP::eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective)
{

    FP t0   = 0;
    FP tmax = 100;
    FP dt   = 0.2;
    std::vector<FP> grid(numIntervals_ + 1);
    for (int i = 0; i < numIntervals_ + 1; ++i) {
        grid[i] = (tmax / numIntervals_) * i + (t0 / numIntervals_) * (numIntervals_ - i);
    }
    mio::oseair::Model<FP> model;
    auto& params = model.parameters;

    set_initial_values(model);
    int gridindex = 0;
    for (int controlIndex = 0; controlIndex < numControlIntervals_; ++controlIndex) {
        model.parameters.template get<mio::oseair::AlphaA<FP>>() = x[controlIndex];
        model.parameters.template get<mio::oseair::AlphaI<FP>>() = x[controlIndex + numControlIntervals_];
        model.parameters.template get<mio::oseair::Kappa<FP>>()  = x[controlIndex + 2 * numControlIntervals_];

        for (int i = 0; i < pcresolution_; ++i, ++gridindex) {

            auto result = mio::simulate<FP, mio::oseair::Model<FP>>(grid[gridindex], grid[gridindex + 1], dt, model);

            for (int j = 0; j < (int)mio::oseair::InfectionState::Count; ++j) {
                model.populations[mio::oseair::InfectionState(j)] = result.get_last_value()[j];
            }
            if (gridindex == numIntervals_ - 1) {
                objective = result.get_last_value()[(int)mio::oseair::InfectionState::ObjectiveFunction] - (tmax - t0);
            }
            constraints[gridindex] = result.get_last_value()[(int)mio::oseair::InfectionState::Infected];
        }
    }

    return;
}

int main()
{
    //switch of logging for mio
    mio::set_log_level(mio::LogLevel::off);

    // Create a new instance of your nlp
    //  (use a SmartPtr, not raw)
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new Seair_NLP();

    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // Change some options
    // Note: The following choices are only examples, they might not be
    //       suitable for your optimization problem.
    app->Options()->SetNumericValue("tol", 1e-6);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("limited_memory_update_type", "bfgs");

    // Initialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int)status;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if (status == Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else {
        std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    }

    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.

    return (int)status;
}

int main_old()
{
    Seair_NLP nlp;
    using FP = ad::gt1s<double>::type;
    std::vector<double> x(nlp.getN());
    std::vector<FP> t1_x(nlp.getN());
    std::vector<FP> constraints(nlp.getM());
    FP objective;

    for (size_t i = 0; i < x.size(); ++i) {
        ad::value(t1_x[i]) = x[i];
    }
    ad::derivative(t1_x[0]) = 1;

    nlp.get_starting_point(nlp.getN(), true, x.data(), false, nullptr, nullptr, nlp.getM(), false, nullptr);
    nlp.eval_objective_constraints(t1_x, constraints, objective);
    std::cout << "directional derivitive of objective function is  " << ad::derivative(objective) << std::endl;

    return 0;
}

bool Seair_NLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
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

bool Seair_NLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m,
                                Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    // controls order: 1. alpha_a, 2. alpha_i, 3. kappa
    for (int i = 0; i < numControlIntervals_; ++i) {
        x_l[i]                            = 0.05; // lower bound of alpha_a
        x_u[i]                            = 0.5; // upper bound of alpha_a
        x_l[i + numControlIntervals_]     = 0.01; // lower bound of alpha_i
        x_u[i + numControlIntervals_]     = 0.3; // upper bound of alpha_i
        x_l[i + 2 * numControlIntervals_] = 0.15; // lower bound of kappa
        x_u[i + 2 * numControlIntervals_] = 0.3; // upper bound of kappa
    }

    // path constraints
    for (int i = 0; i < m_; ++i) {
        g_l[i] = 0.0;
        g_u[i] = 1e6 / N;
    }
    return true;
}

bool Seair_NLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
                                   Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda)
{
    assert(init_z == false);
    assert(init_lambda == false);

    for (int i = 0; i < n; ++i) {
        x[i] = 0.2;
    }
    return true;
}

bool Seair_NLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
    std::vector<double> xx(getN());
    std::vector<double> constraints(getM());
    for (int i = 0; i < n; ++i)
        xx[i] = x[i];
    eval_objective_constraints(xx, constraints, obj_value);
    return true;
}

bool Seair_NLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    using FP = ad::gt1s<double>::type;
    std::vector<FP> xx(getN());
    std::vector<FP> constraints(getM());
    FP objective;
    for (int i = 0; i < n; ++i)
        ad::value(xx[i]) = x[i];
    for (int i = 0; i < n; ++i) {
        ad::derivative(xx[i]) = 1.0;
        eval_objective_constraints(xx, constraints, objective);
        grad_f[i]             = ad::derivative(objective);
        ad::derivative(xx[i]) = 0.0;
    }
    return true;
}

bool Seair_NLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
    std::vector<double> xx(getN());
    std::vector<double> constraints(getM());
    double obj_value = 0;
    for (int i = 0; i < n; ++i)
        xx[i] = x[i];
    eval_objective_constraints(xx, constraints, obj_value);
    for (int i = 0; i < m; ++i)
        g[i] = constraints[i];
    return true;
}

bool Seair_NLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                           Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{

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
        using FP = ad::gt1s<double>::type;
        std::vector<FP> xx(getN());
        std::vector<FP> constraints(getM());
        FP objective;
        int jac_index = 0;
        for (int i = 0; i < n; ++i)
            ad::value(xx[i]) = x[i];
        for (int i = 0; i < n; ++i) {
            ad::derivative(xx[i]) = 1.0;
            eval_objective_constraints(xx, constraints, objective);
            for (int j = 0; j < m; ++j) {
                values[jac_index] = ad::derivative(constraints[j]);
                ++jac_index;
            }
            ad::derivative(xx[i]) = 0.0;
        }
    }
    return true;
}

bool Seair_NLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
                       const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                       Ipopt::Index* jCol, Ipopt::Number* values)
{
    return true;
}

void Seair_NLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                  const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                  const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                  const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    std::cout << "optimal solution is\n";
    std::cout << "Writing output to text files" << std::endl;
    using FP  = double;
    using IS  = mio::oseair::InfectionState;
    using Idx = mio::Index<mio::oseair::InfectionState>;

    FP t0   = 0;
    FP tmax = 100;
    FP dt   = 0.2;
    std::vector<FP> grid(numIntervals_ + 1);
    for (int i = 0; i < numIntervals_ + 1; ++i) {
        grid[i] = (tmax / numIntervals_) * i + (t0 / numIntervals_) * (numIntervals_ - i);
    }
    mio::oseair::Model<FP> model;
    auto& params = model.parameters;

    //open files for parameter output
    std::ofstream outFileAlphaA("alpha_a.txt");
    std::ofstream outFileAlphaI("alpha_i.txt");
    std::ofstream outFileKappa("kappa.txt");

    //open files for state output
    std::ofstream outFileS("s.txt");
    std::ofstream outFileE("e.txt");
    std::ofstream outFileA("a.txt");
    std::ofstream outFileI("i.txt");
    std::ofstream outFileR("r.txt");
    std::ofstream outFileP("p.txt");

    set_initial_values(model);
    int gridindex = 0;
    for (int controlIndex = 0; controlIndex < numControlIntervals_; ++controlIndex) {
        model.parameters.template get<mio::oseair::AlphaA<FP>>() = x[controlIndex];
        model.parameters.template get<mio::oseair::AlphaI<FP>>() = x[controlIndex + numControlIntervals_];
        model.parameters.template get<mio::oseair::Kappa<FP>>()  = x[controlIndex + 2 * numControlIntervals_];

        outFileAlphaA << grid[gridindex] << " " << model.parameters.template get<mio::oseair::AlphaA<FP>>() << "\n";
        outFileAlphaI << grid[gridindex] << " " << model.parameters.template get<mio::oseair::AlphaI<FP>>() << "\n";
        outFileKappa << grid[gridindex] << " " << model.parameters.template get<mio::oseair::Kappa<FP>>() << "\n";

        outFileS << grid[gridindex] << " " << model.populations[{Idx(IS::Susceptible)}] * N / 1000.0 << "\n";
        outFileE << grid[gridindex] << " " << model.populations[{Idx(IS::Exposed)}] * N / 1000.0 << "\n";
        outFileA << grid[gridindex] << " " << model.populations[{Idx(IS::Asymptomatic)}] * N / 1000.0 << "\n";
        outFileI << grid[gridindex] << " " << model.populations[{Idx(IS::Infected)}] * N / 1000.0 << "\n";
        outFileR << grid[gridindex] << " " << model.populations[{Idx(IS::Recovered)}] * N / 1000.0 << "\n";
        outFileP << grid[gridindex] << " " << model.populations[{Idx(IS::Perished)}] * N / 1000.0 << "\n";

        for (int i = 0; i < pcresolution_; ++i, ++gridindex) {

            auto result = mio::simulate<FP, mio::oseair::Model<FP>>(grid[gridindex], grid[gridindex + 1], dt, model);

            for (int j = 0; j < (int)mio::oseair::InfectionState::Count; ++j) {
                model.populations[mio::oseair::InfectionState(j)] = result.get_last_value()[j];
            }
            outFileS << grid[gridindex + 1] << " " << model.populations[{Idx(IS::Susceptible)}] * N / 1000.0 << "\n";
            outFileE << grid[gridindex + 1] << " " << model.populations[{Idx(IS::Exposed)}] * N / 1000.0 << "\n";
            outFileA << grid[gridindex + 1] << " " << model.populations[{Idx(IS::Asymptomatic)}] * N / 1000.0 << "\n";
            outFileI << grid[gridindex + 1] << " " << model.populations[{Idx(IS::Infected)}] * N / 1000.0 << "\n";
            outFileR << grid[gridindex + 1] << " " << model.populations[{Idx(IS::Recovered)}] * N / 1000.0 << "\n";
            outFileP << grid[gridindex + 1] << " " << model.populations[{Idx(IS::Perished)}] * N / 1000.0 << "\n";
        }

        outFileAlphaA << grid[gridindex] << " " << model.parameters.template get<mio::oseair::AlphaA<FP>>() << "\n";
        outFileAlphaI << grid[gridindex] << " " << model.parameters.template get<mio::oseair::AlphaI<FP>>() << "\n";
        outFileKappa << grid[gridindex] << " " << model.parameters.template get<mio::oseair::Kappa<FP>>() << "\n";
    }

    //close files
    outFileAlphaA.close();
    outFileAlphaI.close();
    outFileKappa.close();
    outFileS.close();
    outFileE.close();
    outFileA.close();
    outFileI.close();
    outFileR.close();
    outFileP.close();

    return;
}

bool Seair_NLP::intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
                                      Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
                                      Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du,
                                      Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
                                      Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    return true;
}
