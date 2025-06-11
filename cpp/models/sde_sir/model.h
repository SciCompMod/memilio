/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Nils Wassmuth, Rene Schmieding, Martin J. Kuehn
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

#ifndef MIO_SDE_SIR_MODEL_H
#define MIO_SDE_SIR_MODEL_H

#include "memilio/compartments/stochastic_model.h"
#include "memilio/compartments/flow_model.h"
#include "memilio/epidemiology/populations.h"
#include "sde_sir/infection_state.h"
#include "sde_sir/parameters.h"

namespace mio
{
namespace ssir
{

/********************
 * define the model *
 ********************/

using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Infected>,
                       Flow<InfectionState::Infected, InfectionState::Recovered>>;

class Model : public mio::StochasticModel<ScalarType, InfectionState, mio::Populations<ScalarType, InfectionState>,
                                          Parameters, Flows>
{
public:
    using Base = mio::StochasticModel<ScalarType, InfectionState, mio::Populations<ScalarType, InfectionState>,
                                      Parameters, Flows>;

    Model()
        : Base(typename Base::Populations({InfectionState::Count}, 0.), typename Base::ParameterSet())
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                   ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> flows) const
    {
        auto& params         = Base::parameters;
        ScalarType coeffStoI = params.template get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                               params.template get<TransmissionProbabilityOnContact>() / Base::populations.get_total();

        flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>()] =
            coeffStoI * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
        flows[this->template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>()] =
            (1.0 / params.template get<TimeInfected>()) * y[(size_t)InfectionState::Infected];
    }

    void get_noise(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                   ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> noise) const
    {
        get_flows(pop, y, t, noise);
        noise = noise.array().sqrt();
    }
};

} // namespace ssir
} // namespace mio

#endif // MIO_SDE_SIR_MODEL_H
