/*
* Copyright (C) 2020-2026 MEmilio
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

template <typename FP>
class Model
    : public mio::StochasticModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>, Flows>
{
    using Base = mio::StochasticModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>, Flows>;

public:
    Model()
        : Base(typename Base::Populations({InfectionState::Count}, 0.0), typename Base::ParameterSet())
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const
    {
        const auto& params = this->parameters;
        params.template get<ContactPatterns<FP>>().get_matrix_at(SimulationTime<FP>(t))(0, 0);
        FP coeffStoIV1 = params.template get<ContactPatterns<FP>>().get_matrix_at(SimulationTime<FP>(t))(0, 0) *
                         params.template get<TransmissionProbabilityOnContactV1<FP>>() / this->populations.get_total();
        FP coeffStoIV2 = params.template get<ContactPatterns<FP>>().get_matrix_at(SimulationTime<FP>(t))(0, 0) *
                         params.template get<TransmissionProbabilityOnContactV2<FP>>() / this->populations.get_total();

        flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV1>()] =
            coeffStoIV1 * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::InfectedV1];

        flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV2>()] =
            coeffStoIV2 * y[(size_t)InfectionState::Susceptible] *
            (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2]);

        flows[this->template get_flat_flow_index<InfectionState::ExposedV1, InfectionState::InfectedV1>()] =
            (1.0 / params.template get<TimeExposedV1<FP>>()) * y[(size_t)InfectionState::ExposedV1];

        flows[this->template get_flat_flow_index<InfectionState::ExposedV2, InfectionState::InfectedV2>()] =
            (1.0 / params.template get<TimeExposedV2<FP>>()) * y[(size_t)InfectionState::ExposedV2];

        flows[this->template get_flat_flow_index<InfectionState::InfectedV1, InfectionState::RecoveredV1>()] =
            (1.0 / params.template get<TimeInfectedV1<FP>>()) * y[(size_t)InfectionState::InfectedV1];

        flows[this->template get_flat_flow_index<InfectionState::InfectedV2, InfectionState::RecoveredV2>()] =
            (1.0 / params.template get<TimeInfectedV2<FP>>()) * y[(size_t)InfectionState::InfectedV2];

        flows[this->template get_flat_flow_index<InfectionState::RecoveredV1, InfectionState::ExposedV1V2>()] =
            coeffStoIV2 * y[(size_t)InfectionState::RecoveredV1] *
            (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2]);

        flows[this->template get_flat_flow_index<InfectionState::ExposedV1V2, InfectionState::InfectedV1V2>()] =
            (1.0 / params.template get<TimeExposedV2<FP>>()) * y[(size_t)InfectionState::ExposedV1V2];

        flows[this->template get_flat_flow_index<InfectionState::InfectedV1V2, InfectionState::RecoveredV1V2>()] =
            (1.0 / params.template get<TimeInfectedV2<FP>>()) * y[(size_t)InfectionState::InfectedV1V2];
    }

    void get_noise(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> noise) const
    {
        Eigen::VectorX<FP> flows(Flows::size());
        get_flows(pop, y, t, flows);
        flows = flows.array().sqrt() * Base::white_noise(Flows::size()).array();
        this->get_derivatives(flows, noise);
    }
};

} // namespace sseirvv
} // namespace mio

#endif // MIO_SDE_SEIRVV_MODEL_H
