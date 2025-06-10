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
#include "memilio/epidemiology/populations.h"
#include "memilio/utils/random_number_generator.h"
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
class Model : public FlowModel<FP, InfectionState, Populations<FP, InfectionState>, Parameters<FP>, Flows>
{
    using Base = FlowModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>, Flows>;

public:
    Model()
        : Base(typename Base::Populations({InfectionState::Count}, 0.0), typename Base::ParameterSet())
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const
    {
        using std::sqrt;

        auto& params = this->parameters;
        params.template get<ContactPatterns<FP>>().get_matrix_at(SimulationTime<FP>(t))(0, 0);
        FP coeffStoIV1 = params.template get<ContactPatterns<FP>>().get_matrix_at(SimulationTime<FP>(t))(0, 0) *
                         params.template get<TransmissionProbabilityOnContactV1<FP>>() / this->populations.get_total();
        FP coeffStoIV2 = params.template get<ContactPatterns<FP>>().get_matrix_at(SimulationTime<FP>(t))(0, 0) *
                         params.template get<TransmissionProbabilityOnContactV2<FP>>() / this->populations.get_total();

        // Normal distributed values for the stochastic part of the flows, variables are encoded
        // in the following way: x_y is the stochastic part for the flow from x to y. Variant
        // specific compartments also get an addendum v1 or v2 denoting the relevant variant.

        FP s_ev1       = mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        FP s_ev2       = mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        FP ev1_iv1     = mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        FP ev2_iv2     = mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        FP iv1_rv1     = mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        FP iv2_rv2     = mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        FP rv1_ev1v2   = mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        FP ev1v2_iv1v2 = mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        FP iv1v2_rv1v2 = mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);

        // Assuming that no person can change its InfectionState twice in a single time step,
        // take the minimum of the calculated flow and the source compartment, to ensure that
        // no compartment attains negative values.

        // Calculate inv_step_size and inv_sqrt_step_size for optimization.
        FP inv_step_size      = 1.0 / step_size;
        FP inv_sqrt_step_size = 1.0 / sqrt(step_size);

        // Two outgoing flows from S so will clamp their sum to S * inv_step_size to ensure non-negative S.
        const FP outflow1 = std::clamp<FP>(
            coeffStoIV1 * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::InfectedV1] +
                sqrt(coeffStoIV1 * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::InfectedV1]) *
                    inv_sqrt_step_size * s_ev1,
            0.0, y[(size_t)InfectionState::Susceptible] * inv_step_size);

        const FP outflow2 = std::clamp<FP>(
            coeffStoIV1 * y[(size_t)InfectionState::Susceptible] *
                    (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2]) +
                sqrt(coeffStoIV2 * y[(size_t)InfectionState::Susceptible] *
                     (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2])) *
                    inv_sqrt_step_size * s_ev2,
            0.0, y[(size_t)InfectionState::Susceptible] * inv_step_size);

        const FP outflow_sum = outflow1 + outflow2;
        if (outflow_sum > 0) {
            const FP scale =
                std::clamp(outflow_sum, 0.0, y[(size_t)InfectionState::Susceptible] * inv_step_size) / outflow_sum;
            flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV1>()] =
                outflow1 * scale;
            flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV2>()] =
                outflow2 * scale;
        }
        else {
            flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV1>()] = 0;
            flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV2>()] = 0;
        }

        flows[this->template get_flat_flow_index<InfectionState::ExposedV1, InfectionState::InfectedV1>()] =
            std::clamp<FP>(
                (1.0 / params.template get<TimeExposedV1<FP>>()) * y[(size_t)InfectionState::ExposedV1] +
                    sqrt((1.0 / params.template get<TimeExposedV1<FP>>()) * y[(size_t)InfectionState::ExposedV1]) *
                        inv_sqrt_step_size * ev1_iv1,
                0.0, y[(size_t)InfectionState::ExposedV1] * inv_step_size);

        flows[this->template get_flat_flow_index<InfectionState::ExposedV2, InfectionState::InfectedV2>()] =
            std::clamp<FP>(
                (1.0 / params.template get<TimeExposedV2<FP>>()) * y[(size_t)InfectionState::ExposedV2] +
                    sqrt((1.0 / params.template get<TimeExposedV2<FP>>()) * y[(size_t)InfectionState::ExposedV2]) *
                        inv_sqrt_step_size * ev2_iv2,
                0.0, y[(size_t)InfectionState::ExposedV2] * inv_step_size);

        flows[this->template get_flat_flow_index<InfectionState::InfectedV1, InfectionState::RecoveredV1>()] =
            std::clamp<FP>(
                (1.0 / params.template get<TimeInfectedV1<FP>>()) * y[(size_t)InfectionState::InfectedV1] +
                    sqrt((1.0 / params.template get<TimeInfectedV1<FP>>()) * y[(size_t)InfectionState::InfectedV1]) *
                        inv_sqrt_step_size * iv1_rv1,
                0.0, y[(size_t)InfectionState::InfectedV1] * inv_step_size);

        flows[this->template get_flat_flow_index<InfectionState::InfectedV2, InfectionState::RecoveredV2>()] =
            std::clamp<FP>(
                (1.0 / params.template get<TimeInfectedV2<FP>>()) * y[(size_t)InfectionState::InfectedV2] +
                    sqrt((1.0 / params.template get<TimeInfectedV2<FP>>()) * y[(size_t)InfectionState::InfectedV2]) *
                        inv_sqrt_step_size * iv2_rv2,
                0.0, y[(size_t)InfectionState::InfectedV2] * inv_step_size);

        flows[this->template get_flat_flow_index<InfectionState::RecoveredV1, InfectionState::ExposedV1V2>()] =
            std::clamp<FP>(
                coeffStoIV2 * y[(size_t)InfectionState::RecoveredV1] *
                        (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2]) +
                    sqrt(coeffStoIV2 * y[(size_t)InfectionState::RecoveredV1] *
                         (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2])) *
                        inv_sqrt_step_size * rv1_ev1v2,
                0.0, y[(size_t)InfectionState::RecoveredV1] * inv_step_size);

        flows[this->template get_flat_flow_index<InfectionState::ExposedV1V2, InfectionState::InfectedV1V2>()] =
            std::clamp<FP>(
                (1.0 / params.template get<TimeExposedV2<FP>>()) * y[(size_t)InfectionState::ExposedV1V2] +
                    sqrt((1.0 / params.template get<TimeExposedV2<FP>>()) * y[(size_t)InfectionState::ExposedV1V2]) /
                        sqrt(step_size) * ev1v2_iv1v2,
                0.0, y[(size_t)InfectionState::ExposedV1V2] * inv_step_size);

        flows[this->template get_flat_flow_index<InfectionState::InfectedV1V2, InfectionState::RecoveredV1V2>()] =
            std::clamp<FP>(
                (1.0 / params.template get<TimeInfectedV2<FP>>()) * y[(size_t)InfectionState::InfectedV1V2] +
                    sqrt((1.0 / params.template get<TimeInfectedV2<FP>>()) * y[(size_t)InfectionState::InfectedV1V2]) /
                        sqrt(step_size) * iv1v2_rv1v2,
                0.0, y[(size_t)InfectionState::InfectedV1V2] * inv_step_size);
    }

    FP step_size; ///< A step size of the model with which the stochastic process is realized.
    mutable RandomNumberGenerator rng;

private:
};

} // namespace sseirvv
} // namespace mio

#endif // MIO_SDE_SEIRVV_MODEL_H
