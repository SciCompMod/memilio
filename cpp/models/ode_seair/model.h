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

/* This model is an extented SEIR type model of the COVID-19 pandemic in the US
 * that also includes asymptomatic and dead people.
 * A detailed description of the model can be found in the publication
 * Tsay et al. (2020), Modeling, state estimation, and optimal control for the US COVID-19 outbreak */

#ifndef ODESEAIR_MODEL_H
#define ODESEAIR_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/config.h"
#include "memilio/epidemiology/populations.h"
#include "ode_seair/infection_state.h"
#include "ode_seair/parameters.h"

namespace mio
{
namespace oseair
{

/********************
 * define the model *
 ********************/
template <typename FP = ScalarType>
class Model : public mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>>
{
    using Base = mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>>;

public:
    Model()
        : Base(mio::Populations<FP, InfectionState>({InfectionState::Count}, 0.), typename Base::ParameterSet())
    {
    }

    void get_derivatives(Eigen::Ref<const Vector<FP>> pop, Eigen::Ref<const Vector<FP>> y, FP /* t */,
                         Eigen::Ref<Vector<FP>> dydt) const override
    {
        auto& params = this->parameters;

        const auto pop_total = pop.sum();

        auto& alpha_a          = params.template get<AlphaA<FP>>();
        auto& alpha_i          = params.template get<AlphaI<FP>>();
        auto& kappa            = params.template get<Kappa<FP>>();
        auto& beta             = params.template get<Beta<FP>>();
        auto& mu               = params.template get<Mu<FP>>();
        auto& t_latent_inverse = params.template get<TLatentInverse<FP>>();
        auto& rho              = params.template get<Rho<FP>>();
        auto& gamma            = params.template get<Gamma<FP>>();

        dydt[(size_t)InfectionState::Susceptible] =
            -alpha_a / pop_total * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Asymptomatic] -
            alpha_i / pop_total * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected] +
            gamma * y[(size_t)InfectionState::Recovered];
        dydt[(size_t)InfectionState::Exposed] =

            alpha_a / pop_total * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Asymptomatic] +
            alpha_i / pop_total * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected] -
            t_latent_inverse * y[(size_t)InfectionState::Exposed];
        dydt[(size_t)InfectionState::Asymptomatic] = t_latent_inverse * y[(size_t)InfectionState::Exposed] -
                                                     (kappa + rho) * y[(size_t)InfectionState::Asymptomatic];
        dydt[(size_t)InfectionState::Infected] =
            kappa * y[(size_t)InfectionState::Asymptomatic] - (beta + mu) * y[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Recovered] = rho * y[(size_t)InfectionState::Asymptomatic] +
                                                  beta * y[(size_t)InfectionState::Infected] -
                                                  gamma * y[(size_t)InfectionState::Recovered];
        dydt[(size_t)InfectionState::Dead] = mu * y[(size_t)InfectionState::Infected];
    }
};

} // namespace oseair
} // namespace mio

#endif // ODESEAIR_MODEL_H
