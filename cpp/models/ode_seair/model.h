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

/* This model is a extented SEIR type model of the COVID-19 pandemic in the US
 * that als includes asymptomatic and perished people.
 * A detailed description of the model can be found in the publication
 * Tsay et al. (2020), Modeling, state estimation, and optimal control for the US COVID-19 outbreak */

#ifndef ODESEAIR_MODEL_H
#define ODESEAIR_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h" // IWYU pragma: keep
#include "ode_seair/infection_state.h"
#include "ode_seair/parameters.h"

namespace mio
{
namespace oseair
{

/********************
    * define the model *
    ********************/
template <typename FP = double>
class Model : public mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>>
{
    using Base = mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>>;

public:
    Model()
        : Base(mio::Populations<FP, InfectionState>({InfectionState::Count}, 0.), typename Base::ParameterSet())
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::Matrix<FP, Eigen::Dynamic, 1>> /* pop */,
                         Eigen::Ref<const Eigen::Matrix<FP, Eigen::Dynamic, 1>> y, FP /* t */,
                         Eigen::Ref<Eigen::Matrix<FP, Eigen::Dynamic, 1>> dydt) const override
    {
        auto& params = this->parameters;

        auto& alpha_a          = params.template get<AlphaA<FP>>();
        auto& alpha_i          = params.template get<AlphaI<FP>>();
        auto& kappa            = params.template get<Kappa<FP>>();
        auto& beta             = params.template get<Beta<FP>>();
        auto& mu               = params.template get<Mu<FP>>();
        auto& t_latent_inverse = params.template get<TLatentInverse<FP>>();
        auto& rho              = params.template get<Rho<FP>>();
        auto& gamma            = params.template get<Gamma<FP>>();

        const auto& s = y[(size_t)InfectionState::Susceptible];
        const auto& e = y[(size_t)InfectionState::Exposed];
        const auto& a = y[(size_t)InfectionState::Asymptomatic];
        const auto& i = y[(size_t)InfectionState::Infected];
        const auto& r = y[(size_t)InfectionState::Recovered];

        dydt[(size_t)InfectionState::Susceptible]       = -alpha_a * s * a - alpha_i * s * i + gamma * r;
        dydt[(size_t)InfectionState::Exposed]           = alpha_a * s * a + alpha_i * s * i - t_latent_inverse * e;
        dydt[(size_t)InfectionState::Asymptomatic]      = t_latent_inverse * e - kappa * a - rho * a;
        dydt[(size_t)InfectionState::Infected]          = kappa * a - beta * i - mu * i;
        dydt[(size_t)InfectionState::Recovered]         = rho * a + beta * i - gamma * r;
        dydt[(size_t)InfectionState::Perished]          = mu * i;
        dydt[(size_t)InfectionState::ObjectiveFunction] = 1 - alpha_i - alpha_a + 0.1 * kappa;
    }
};

} // namespace oseair
} // namespace mio

#endif // ODESEAIR_MODEL_H
