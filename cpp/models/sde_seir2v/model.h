/* 
* Copyright (C) 2020-2024 MEmilio
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

#ifndef MIO_SDE_SEIR2V_MODEL_H
#define MIO_SDE_SEIR2V_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/utils/random_number_generator.h"
#include "sde_seir2v/infection_state.h"
#include "sde_seir2v/parameters.h"

namespace mio
{
namespace sseir2v
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

class Model : public FlowModel<InfectionState, Populations<InfectionState>, Parameters, Flows>
{
    using Base = FlowModel<InfectionState, mio::Populations<InfectionState>, Parameters, Flows>;

public:
    Model()
        : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                   Eigen::Ref<Eigen::VectorXd> flows) const
    {
        auto& params     = this->parameters;
        params.get<ContactPatterns>().get_matrix_at(t)(0, 0);
        double coeffStoIV1 = /*params.get<ContactPatterns>().get_matrix_at(t)(0, 0) **/
                           params.get<TransmissionProbabilityOnContactV1>() / populations.get_total();
        double coeffStoIV2 = /*params.get<ContactPatterns>().get_matrix_at(t)(0, 0) **/
                           params.get<TransmissionProbabilityOnContactV2>() / populations.get_total();

        std::initializer_list<uint32_t> seeds = {14159265u, 35897932u};
        rng.seed(seeds);

        double s_ev1 = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double s_ev2 = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double ev1_iv1 = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double ev2_iv2 = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double iv1_rv1 = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double iv2_rv2 = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double rv1_ev1v2 = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double ev1v2_iv1v2 = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double iv1v2_rv1v2 = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);  


        const double inv_sqrt_dt = 1 / sqrt(step_size);

        // Assuming that no person can change its InfectionState twice in a single time step,
        // take the minimum of the calculated flow and the source compartment, to ensure that
        // no compartment attains negative values.

        flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV1>()] = std::max(
            std::min(
                coeffStoIV1 * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::InfectedV1] +
                    sqrt(coeffStoIV1 * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::InfectedV1]) *
                        inv_sqrt_dt * s_ev1,
                y[(size_t)InfectionState::Susceptible] / step_size),
            0.);
         
        flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV2>()] = std::max(
            std::min(
                coeffStoIV1 * y[(size_t)InfectionState::Susceptible] * 
                    (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2]) +
                    sqrt(coeffStoIV2 * y[(size_t)InfectionState::Susceptible] * (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2])) *
                        inv_sqrt_dt * s_ev2,
                y[(size_t)InfectionState::Susceptible] / step_size),
            0.);  

        flows[get_flat_flow_index<InfectionState::ExposedV1, InfectionState::InfectedV1>()] = std::max(
            std::min(
                (1.0 / params.get<TimeExposedV1>()) * y[(size_t)InfectionState::ExposedV1] +
                    sqrt((1.0 / params.get<TimeExposedV1>()) * y[(size_t)InfectionState::ExposedV1]) * inv_sqrt_dt * ev1_iv1,
                y[(size_t)InfectionState::ExposedV1] / step_size),
            0.);            

        flows[get_flat_flow_index<InfectionState::ExposedV2, InfectionState::InfectedV2>()] = std::max(
            std::min(
                (1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV2] +
                    sqrt((1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV2]) * inv_sqrt_dt * ev2_iv2,
                y[(size_t)InfectionState::ExposedV2] / step_size),
            0.);

        flows[get_flat_flow_index<InfectionState::InfectedV1, InfectionState::RecoveredV1>()] = std::max(
            std::min(
                (1.0 / params.get<TimeInfectedV1>()) * y[(size_t)InfectionState::InfectedV1] +
                    sqrt((1.0 / params.get<TimeInfectedV1>()) * y[(size_t)InfectionState::InfectedV1]) * inv_sqrt_dt * iv1_rv1,
                y[(size_t)InfectionState::InfectedV1] / step_size),
            0.);

        flows[get_flat_flow_index<InfectionState::InfectedV2, InfectionState::RecoveredV2>()] = std::max(
            std::min(
                (1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV2] +
                    sqrt((1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV2]) * inv_sqrt_dt * iv2_rv2,
                y[(size_t)InfectionState::InfectedV2] / step_size),
            0.);

        flows[get_flat_flow_index<InfectionState::RecoveredV1, InfectionState::ExposedV1V2>()] = std::max(
            std::min(
                coeffStoIV2 * y[(size_t)InfectionState::RecoveredV1] 
                    * pop[(size_t)InfectionState::InfectedV1V2]+ pop[(size_t)InfectionState::InfectedV2] +
                    sqrt(coeffStoIV2 * y[(size_t)InfectionState::RecoveredV1] * (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2])) *
                        inv_sqrt_dt * rv1_ev1v2,
                y[(size_t)InfectionState::RecoveredV1] / step_size),
            0.);

        flows[get_flat_flow_index<InfectionState::ExposedV1V2, InfectionState::InfectedV1V2>()] = std::max(
            std::min(
                (1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV1V2] +
                    sqrt((1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV1V2]) * inv_sqrt_dt * ev1v2_iv1v2,
                y[(size_t)InfectionState::ExposedV1V2] / step_size),
            0.);

        flows[get_flat_flow_index<InfectionState::InfectedV1V2, InfectionState::RecoveredV1V2>()] = std::max(
            std::min(
                (1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV1V2] +
                    sqrt((1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV1V2]) * inv_sqrt_dt * iv1v2_rv1v2,
                y[(size_t)InfectionState::InfectedV1V2] / step_size),
            0.);
    }

    ScalarType step_size = 1.; ///< A step size of the model with which the stochastic process is realized.
    mutable RandomNumberGenerator rng;

private:
};

} // namespace sseir2v
} // namespace mio

#endif // MIO_SDE_SEIR2V_MODEL_H
