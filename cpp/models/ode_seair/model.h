/*
* Copyright (C) 2020-2026 MEmilio
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

#ifndef ODESEAIR_MODEL_H
#define ODESEAIR_MODEL_H

#include "memilio/compartments/compartmental_model.h"
#include "memilio/config.h"
#include "memilio/epidemiology/populations.h"
#include "ode_seair/infection_state.h"
#include "ode_seair/parameters.h"

namespace mio
{
namespace oseair
{

template <typename FP>
class Model : public mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>>
{
    using Base = mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;

    Model()
        : Base(Populations({InfectionState::Count}, 0.0), ParameterSet())
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP /* t */,
                         Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
    {
        auto& params         = this->parameters;
        const auto pop_total = pop.sum();

        dydt[(size_t)InfectionState::Susceptible] =
            -params.template get<SocialDistancing<FP>>() / pop_total * y[(size_t)InfectionState::Susceptible] *
                pop[(size_t)InfectionState::Asymptomatic] -
            params.template get<Quarantined<FP>>() / pop_total * y[(size_t)InfectionState::Susceptible] *
                pop[(size_t)InfectionState::Infected] +
            params.template get<TimeRecoveredInv<FP>>() * y[(size_t)InfectionState::Recovered];
        dydt[(size_t)InfectionState::Exposed] =
            params.template get<SocialDistancing<FP>>() / pop_total * y[(size_t)InfectionState::Susceptible] *
                pop[(size_t)InfectionState::Asymptomatic] +
            params.template get<Quarantined<FP>>() / pop_total * y[(size_t)InfectionState::Susceptible] *
                pop[(size_t)InfectionState::Infected] -
            y[(size_t)InfectionState::Exposed] / params.template get<TimeExposed<FP>>();
        dydt[(size_t)InfectionState::Asymptomatic] =
            y[(size_t)InfectionState::Exposed] / params.template get<TimeExposed<FP>>() -
            (params.template get<TestingRate<FP>>() + params.template get<RecoveryRateFromAsymptomatic<FP>>()) *
                y[(size_t)InfectionState::Asymptomatic];
        dydt[(size_t)InfectionState::Infected] =
            params.template get<TestingRate<FP>>() * y[(size_t)InfectionState::Asymptomatic] -
            (params.template get<RecoveryRate<FP>>() + params.template get<DeathRate<FP>>()) *
                y[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Recovered] =
            params.template get<RecoveryRateFromAsymptomatic<FP>>() * y[(size_t)InfectionState::Asymptomatic] +
            params.template get<RecoveryRate<FP>>() * y[(size_t)InfectionState::Infected] -
            params.template get<TimeRecoveredInv<FP>>() * y[(size_t)InfectionState::Recovered];
        dydt[(size_t)InfectionState::Dead] = params.template get<DeathRate<FP>>() * y[(size_t)InfectionState::Infected];
    }
};

} // namespace oseair
} // namespace mio

#endif // ODESEAIR_MODEL_H
