/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#ifndef SEIR_MODEL_H
#define SEIR_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"

namespace mio
{
namespace oseir
{

    /********************
    * define the model *
    ********************/

    class Model : public CompartmentalModel<InfectionState, Populations<InfectionState>, ParametersBase>
    {
        using Base = CompartmentalModel<InfectionState, mio::Populations<InfectionState>, ParametersBase>;

    public:
        Model()
            : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
        {
        }

        void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                             Eigen::Ref<Eigen::VectorXd> dydt) const override
        {
            auto& params     = this->parameters;
            double S2E_coeff = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                               params.get<InfectionProbabilityFromContact>() / populations.get_total();

            dydt[(size_t)InfectionState::Susceptible] =
                -S2E_coeff * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
            dydt[(size_t)InfectionState::Exposed] =
                S2E_coeff * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected] -
                (1.0 / params.get<LatentTime>()) * y[(size_t)InfectionState::Exposed];
            dydt[(size_t)InfectionState::Infected] =
                (1.0 / params.get<LatentTime>()) * y[(size_t)InfectionState::Exposed] -
                (1.0 / params.get<InfectiousTime>()) * y[(size_t)InfectionState::Infected];
            dydt[(size_t)InfectionState::Recovered] =
                (1.0 / params.get<InfectiousTime>()) * y[(size_t)InfectionState::Infected];
        }
    };

} // namespace oseir
} // namespace mio

#endif // SEIR_MODEL_H
