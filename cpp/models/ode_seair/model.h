/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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
#ifndef SEAIR_MODEL_H
#define SEAIR_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "ode_seair/infection_state.h"
#include "ode_seair/parameters.h"

namespace mio
{
namespace oseair
{

/********************
    * define the model *
    ********************/

class Model : public CompartmentalModel<InfectionState, Populations<InfectionState>, Parameters>
{
    using Base = CompartmentalModel<InfectionState, mio::Populations<InfectionState>, Parameters>;

public:
    Model()
        : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        auto& params     = this->parameters;
        double coeffStoE = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                           params.get<TransmissionProbabilityOnContact>() / populations.get_total();

        dydt[(size_t)InfectionState::Susceptible] =
            -coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Exposed] =
            coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected] -
            (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed];
        dydt[(size_t)InfectionState::Infected] =
            (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed] -
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Recovered] =
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected];
    }
};

} // namespace oseir
} // namespace mio

#endif // SEAIR_MODEL_H
