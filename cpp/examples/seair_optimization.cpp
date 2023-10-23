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



int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    using FP=typename ad::gt1s<double>::type;



    FP  t0   = 0;
    FP  tmax = 100;
    FP  dt   = 0.2;

    mio::log_info("Simulating SEAIR; t={} ... {} with dt = {}.", ad::value(t0), ad::value(tmax), ad::value(dt));

    mio::oseair::Model<FP> model1;
    set_initial_values(model1);
    FP value = model1.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}];
    ad::derivative(value) = 1.0;
    model1.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] = value;
    model1.check_constraints();
    auto seair1 = mio::simulate<mio::oseair::Model<FP>,FP>(t0, tmax, dt, model1);

    mio::oseair::Model<double> model2;
    set_initial_values(model2);
    double h = 1e-4;
    model2.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] += h;
    model2.check_constraints();
    auto seair2 = mio::simulate<mio::oseair::Model<double>,double>(0, 100, 0.2, model2);

    const std::string file_name = "seair.txt";
    std::cout << "Writing output to " << file_name << std::endl;
    mio::time_series_to_file(seair1, file_name);

    auto last1 = seair1.get_last_value();
    auto last2 = seair2.get_last_value();

    std::cout << "Last value is" << std::endl;
    std::cout <<  last2 << std::endl;

    std::cout << "Compare algorithmic differentiation (AD) with finite differences (FD)" << std::endl;
    for (size_t i=0; i < last1.size(); ++i) {
        std::cout << "AD: " << ad::derivative(last1[i]) << "  FD: " << (1./h)*(last2[i]-ad::value(last1[i])) << std::endl;
    }


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
