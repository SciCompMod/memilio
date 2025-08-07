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

#ifndef MIO_SDE_SEIRVV_MODEL_H
#define MIO_SDE_SEIRVV_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/stochastic_model.h"
#include "memilio/epidemiology/populations.h"
#include "sde_seirvv/infection_state.h"
#include "sde_seirvv/parameters.h"

namespace mio
{
namespace sseirvv
{

/********************
 * define the model *
 ********************/

using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::ExposedV1>,
                       Flow<InfectionState::Susceptible, InfectionState::ExposedV2>,
                       Flow<InfectionState::ExposedV1, InfectionState::InfectedV1>,
                       Flow<InfectionState::ExposedV2, InfectionState::InfectedV2>,
                       Flow<InfectionState::InfectedV1, InfectionState::RecoveredV1>,
                       Flow<InfectionState::InfectedV2, InfectionState::RecoveredV2>,
                       Flow<InfectionState::RecoveredV1, InfectionState::ExposedV1V2>,
                       Flow<InfectionState::ExposedV1V2, InfectionState::InfectedV1V2>,
                       Flow<InfectionState::InfectedV1V2, InfectionState::RecoveredV1V2>>;

class Model : public mio::StochasticModel<ScalarType, InfectionState, mio::Populations<ScalarType, InfectionState>,
                                          Parameters, Flows>
{
    using Base = mio::StochasticModel<ScalarType, InfectionState, mio::Populations<ScalarType, InfectionState>,
                                      Parameters, Flows>;

public:
    Model()
        : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                   ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> flows) const
    {
        const auto& params = this->parameters;
        params.get<ContactPatterns>().get_matrix_at(t)(0, 0);
        ScalarType coeffStoIV1 = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                                 params.get<TransmissionProbabilityOnContactV1>() / populations.get_total();
        ScalarType coeffStoIV2 = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                                 params.get<TransmissionProbabilityOnContactV2>() / populations.get_total();

        flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV1>()] =
            coeffStoIV1 * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::InfectedV1];

        flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV2>()] =
            coeffStoIV2 * y[(size_t)InfectionState::Susceptible] *
            (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2]);

        flows[get_flat_flow_index<InfectionState::ExposedV1, InfectionState::InfectedV1>()] =
            (1.0 / params.get<TimeExposedV1>()) * y[(size_t)InfectionState::ExposedV1];

        flows[get_flat_flow_index<InfectionState::ExposedV2, InfectionState::InfectedV2>()] =
            (1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV2];

        flows[get_flat_flow_index<InfectionState::InfectedV1, InfectionState::RecoveredV1>()] =
            (1.0 / params.get<TimeInfectedV1>()) * y[(size_t)InfectionState::InfectedV1];

        flows[get_flat_flow_index<InfectionState::InfectedV2, InfectionState::RecoveredV2>()] =
            (1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV2];

        flows[get_flat_flow_index<InfectionState::RecoveredV1, InfectionState::ExposedV1V2>()] =
            coeffStoIV2 * y[(size_t)InfectionState::RecoveredV1] *
            (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2]);

        flows[get_flat_flow_index<InfectionState::ExposedV1V2, InfectionState::InfectedV1V2>()] =
            (1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV1V2];

        flows[get_flat_flow_index<InfectionState::InfectedV1V2, InfectionState::RecoveredV1V2>()] =
            (1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV1V2];
    }

    void get_noise(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                   ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> noise) const
    {
        Eigen::VectorX<ScalarType> flows(Flows::size());
        get_flows(pop, y, t, flows);
        flows = flows.array().sqrt() * Base::white_noise(Flows::size()).array();
        get_derivatives(flows, noise);
    }
};

} // namespace sseirvv
} // namespace mio

#endif // MIO_SDE_SEIRVV_MODEL_H
