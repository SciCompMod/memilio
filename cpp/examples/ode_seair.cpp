/*
* Copyright (C) 2020-2025 MEmilio
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

// This model is a extented SEIR type model of the COVID-19 pandemic in the US
// that als includes asymptomatic and dead people.
// A detailed description of the model can be found in the publication
// Tsay et al. (2020), Modeling, state estimation, and optimal control for the US COVID-19 outbreak.

#include "memilio/ad/ad.hpp"

#include "ode_seair/model.h"
#include "ode_seair/infection_state.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"

/**
 * @brief set_initial_values sets the initial value of the mio::oseair::Model<FP> model according to
 * the publication of Tsay et al. (2020): Modeling, state estimation, and optimal control for the US COVID-19 outbreak.
 * Note that the total population, i.e., the sum of all compartments, is normalized to one.
 * @tparam FP floating point type, e.g., double
 * @param model an instance of the mio::oseair::Model<FP> which is a compartmental model.
 */
template <typename FP>
void set_initial_values(mio::oseair::Model<FP>& model)
{
    const double N = 327167434; // total population of the United States

    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] =
        0.9977558755803503 * N;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}] =
        0.0003451395725394549 * N;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] =
        0.00037846880968213874 * N;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}]  = 337072.0;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] = 17448.0;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Dead)}]      = 9619.0;
}

int main()
{

    mio::set_log_level(mio::LogLevel::debug);

    using FP = typename ad::gt1s<double>::type; // algorithmic differentiation data type: scalar tangent-linear mode

    FP t0   = 0;
    FP tmax = 10;
    FP dt   = 0.3;

    mio::log_info("Simulating SEAIR; t={} ... {} with dt = {}.", ad::value(t0), ad::value(tmax), ad::value(dt));

    // Compute derivative of the final states of the model with respect to the initial value of the suscetible population
    mio::oseair::Model<FP> model1;
    set_initial_values(model1);
    FP value = model1.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}];
    ad::derivative(value)                                                                                   = 1.0;
    model1.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] = value;
    model1.check_constraints();
    auto seair1 = mio::simulate<FP, mio::oseair::Model<FP>>(t0, tmax, dt, model1);

    // We want to compare the derivatives computed ba algorithmic differention with difference quotient.
    // To this end we perturbe the corresponding initial value of the by an increment h and simulate again.
    mio::oseair::Model<double> model2;
    set_initial_values(model2);
    double h = 1e-4;
    model2.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] += h;
    model2.check_constraints();
    mio::TimeSeries<double> seair2 =
        mio::simulate<double, mio::oseair::Model<double>>(ad::value(t0), ad::value(tmax), ad::value(dt), model2);

    const std::string file_name = "seair-compare.csv";
    std::ofstream file(file_name);
    std::cout << "Writing output to " << file_name << std::endl;
    seair1.print_table({}, 21, 10, file);
    file.close();

    auto last1 = seair1.get_last_value();
    auto last2 = seair2.get_last_value();

    std::cout << "Last value is" << std::endl;
    std::cout << last2 << std::endl;

    std::cout << "Compare algorithmic differentiation (AD) with finite differences (FD)" << std::endl;
    for (Eigen::Index i = 0; i < last1.size(); ++i) {
        std::cout << "AD: " << ad::derivative(last1[i]) << "  FD: " << (1. / h) * (last2[i] - ad::value(last1[i]))
                  << std::endl;
    }

    return 0;
}
