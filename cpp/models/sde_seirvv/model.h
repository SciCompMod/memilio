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

class Model : public FlowModel<ScalarType, InfectionState, Populations<ScalarType, InfectionState>, Parameters, Flows>
{
    using Base = FlowModel<ScalarType, InfectionState, mio::Populations<ScalarType, InfectionState>, Parameters, Flows>;

public:
    Model()
        : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                   ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> flows) const
    {
        std::cout << "calc flows" << std::endl;
        auto& params = this->parameters;

        // Hole den Kontaktmuster-Wert aus der Matrix bei Zeit t
        ScalarType cp_val = params.get<ContactPatterns>().get_matrix_at(t)(0, 0);
        std::cout << "ContactPatterns value: " << cp_val << std::endl;

        // Berechne Koeffizienten für die Transmission
        ScalarType coeffStoIV1 = cp_val * params.get<TransmissionProbabilityOnContactV1>() / populations.get_total();
        std::cout << "coeffStoIV1: " << coeffStoIV1 << std::endl;

        ScalarType coeffStoIV2 = cp_val * params.get<TransmissionProbabilityOnContactV2>() / populations.get_total();
        std::cout << "coeffStoIV2: " << coeffStoIV2 << std::endl;

        // Normalverteilte Zufallszahlen für den stochastischen Anteil der Flows
        ScalarType s_ev1 =
            mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        std::cout << "s_ev1: " << s_ev1 << std::endl;

        ScalarType s_ev2 =
            mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        std::cout << "s_ev2: " << s_ev2 << std::endl;

        ScalarType ev1_iv1 =
            mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        std::cout << "ev1_iv1: " << ev1_iv1 << std::endl;

        ScalarType ev2_iv2 =
            mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        std::cout << "ev2_iv2: " << ev2_iv2 << std::endl;

        ScalarType iv1_rv1 =
            mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        std::cout << "iv1_rv1: " << iv1_rv1 << std::endl;

        ScalarType iv2_rv2 =
            mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        std::cout << "iv2_rv2: " << iv2_rv2 << std::endl;

        ScalarType rv1_ev1v2 =
            mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        std::cout << "rv1_ev1v2: " << rv1_ev1v2 << std::endl;

        ScalarType ev1v2_iv1v2 =
            mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        std::cout << "ev1v2_iv1v2: " << ev1v2_iv1v2 << std::endl;

        ScalarType iv1v2_rv1v2 =
            mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(rng, 0.0, 1.0);
        std::cout << "iv1v2_rv1v2: " << iv1v2_rv1v2 << std::endl;

        // Berechne Optimierungsgrößen
        ScalarType inv_step_size = 1.0 / step_size;
        std::cout << "inv_step_size: " << inv_step_size << std::endl;

        ScalarType inv_sqrt_step_size = 1.0 / sqrt(step_size);
        std::cout << "inv_sqrt_step_size: " << inv_sqrt_step_size << std::endl;

        // Berechne den ersten Outflow aus S
        const ScalarType outflow1 = std::clamp(
            coeffStoIV1 * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::InfectedV1] +
                sqrt(coeffStoIV1 * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::InfectedV1]) *
                    inv_sqrt_step_size * s_ev1,
            0.0, y[(size_t)InfectionState::Susceptible] * inv_step_size);
        std::cout << "outflow1: " << outflow1 << std::endl;

        // Berechne den zweiten Outflow aus S
        const ScalarType outflow2 =
            std::clamp(coeffStoIV1 * y[(size_t)InfectionState::Susceptible] *
                               (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2]) +
                           sqrt(coeffStoIV2 * y[(size_t)InfectionState::Susceptible] *
                                (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2])) *
                               inv_sqrt_step_size * s_ev2,
                       0.0, y[(size_t)InfectionState::Susceptible] * inv_step_size);
        std::cout << "outflow2: " << outflow2 << std::endl;

        // Summe der Outflows
        const ScalarType outflow_sum = outflow1 + outflow2;
        std::cout << "outflow_sum: " << outflow_sum << std::endl;

        if (outflow_sum > 0) {
            const ScalarType scale =
                std::clamp(outflow_sum, 0.0, y[(size_t)InfectionState::Susceptible] * inv_step_size) / outflow_sum;
            std::cout << "scale: " << scale << std::endl;
            flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV1>()] = outflow1 * scale;
            std::cout << "flow Sus -> ExpV1: "
                      << flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV1>()]
                      << std::endl;
            flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV2>()] = outflow2 * scale;
            std::cout << "flow Sus -> ExpV2: "
                      << flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV2>()]
                      << std::endl;
        }
        else {
            flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV1>()] = 0;
            std::cout << "flow Sus -> ExpV1 auf 0 gesetzt" << std::endl;
            flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::ExposedV2>()] = 0;
            std::cout << "flow Sus -> ExpV2 auf 0 gesetzt" << std::endl;
        }

        // Fluss von ExposedV1 zu InfectedV1
        flows[get_flat_flow_index<InfectionState::ExposedV1, InfectionState::InfectedV1>()] =
            std::clamp((1.0 / params.get<TimeExposedV1>()) * y[(size_t)InfectionState::ExposedV1] +
                           sqrt((1.0 / params.get<TimeExposedV1>()) * y[(size_t)InfectionState::ExposedV1]) *
                               inv_sqrt_step_size * ev1_iv1,
                       0.0, y[(size_t)InfectionState::ExposedV1] * inv_step_size);
        std::cout << "flow ExpV1 -> InfV1: "
                  << flows[get_flat_flow_index<InfectionState::ExposedV1, InfectionState::InfectedV1>()] << std::endl;

        // Fluss von ExposedV2 zu InfectedV2
        flows[get_flat_flow_index<InfectionState::ExposedV2, InfectionState::InfectedV2>()] =
            std::clamp((1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV2] +
                           sqrt((1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV2]) *
                               inv_sqrt_step_size * ev2_iv2,
                       0.0, y[(size_t)InfectionState::ExposedV2] * inv_step_size);
        std::cout << "flow ExpV2 -> InfV2: "
                  << flows[get_flat_flow_index<InfectionState::ExposedV2, InfectionState::InfectedV2>()] << std::endl;

        // Fluss von InfectedV1 zu RecoveredV1
        flows[get_flat_flow_index<InfectionState::InfectedV1, InfectionState::RecoveredV1>()] =
            std::clamp((1.0 / params.get<TimeInfectedV1>()) * y[(size_t)InfectionState::InfectedV1] +
                           sqrt((1.0 / params.get<TimeInfectedV1>()) * y[(size_t)InfectionState::InfectedV1]) *
                               inv_sqrt_step_size * iv1_rv1,
                       0.0, y[(size_t)InfectionState::InfectedV1] * inv_step_size);
        std::cout << "flow InfV1 -> RecV1: "
                  << flows[get_flat_flow_index<InfectionState::InfectedV1, InfectionState::RecoveredV1>()] << std::endl;

        // Fluss von InfectedV2 zu RecoveredV2
        flows[get_flat_flow_index<InfectionState::InfectedV2, InfectionState::RecoveredV2>()] =
            std::clamp((1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV2] +
                           sqrt((1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV2]) *
                               inv_sqrt_step_size * iv2_rv2,
                       0.0, y[(size_t)InfectionState::InfectedV2] * inv_step_size);
        std::cout << "flow InfV2 -> RecV2: "
                  << flows[get_flat_flow_index<InfectionState::InfectedV2, InfectionState::RecoveredV2>()] << std::endl;

        // Fluss von RecoveredV1 zu ExposedV1V2
        flows[get_flat_flow_index<InfectionState::RecoveredV1, InfectionState::ExposedV1V2>()] =
            std::clamp(coeffStoIV2 * y[(size_t)InfectionState::RecoveredV1] *
                               (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2]) +
                           sqrt(coeffStoIV2 * y[(size_t)InfectionState::RecoveredV1] *
                                (pop[(size_t)InfectionState::InfectedV1V2] + pop[(size_t)InfectionState::InfectedV2])) *
                               inv_sqrt_step_size * rv1_ev1v2,
                       0.0, y[(size_t)InfectionState::RecoveredV1] * inv_step_size);
        std::cout << "flow RecV1 -> ExpV1V2: "
                  << flows[get_flat_flow_index<InfectionState::RecoveredV1, InfectionState::ExposedV1V2>()]
                  << std::endl;

        // Fluss von ExposedV1V2 zu InfectedV1V2
        flows[get_flat_flow_index<InfectionState::ExposedV1V2, InfectionState::InfectedV1V2>()] =
            std::clamp((1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV1V2] +
                           sqrt((1.0 / params.get<TimeExposedV2>()) * y[(size_t)InfectionState::ExposedV1V2]) /
                               sqrt(step_size) * ev1v2_iv1v2,
                       0.0, y[(size_t)InfectionState::ExposedV1V2] * inv_step_size);
        std::cout << "flow ExpV1V2 -> InfV1V2: "
                  << flows[get_flat_flow_index<InfectionState::ExposedV1V2, InfectionState::InfectedV1V2>()]
                  << std::endl;

        // Fluss von InfectedV1V2 zu RecoveredV1V2
        double time_infected_v2 = params.get<TimeInfectedV2>();
        double infected_v1v2    = y[(size_t)InfectionState::InfectedV1V2];
        double sqrt_term        = sqrt((1.0 / time_infected_v2) * infected_v1v2);
        auto sqrt_term_2        = sqrt(step_size);

        std::cout << "time_infected_v2: " << time_infected_v2 << std::endl;
        std::cout << "infected_v1v2: " << infected_v1v2 << std::endl;
        std::cout << "sqrt_term_1: " << sqrt_term << std::endl;
        std::cout << "sqrt_term_2: " << sqrt_term_2 << std::endl;
        std::cout << "iv1v2_rv1v2: " << iv1v2_rv1v2 << std::endl;
        std::cout << "step_size: " << step_size << std::endl;
        std::cout << "inv_step_size: " << inv_step_size << std::endl;

        double flow_value = std::clamp((1.0 / time_infected_v2) * infected_v1v2 + sqrt_term / sqrt_term_2 * iv1v2_rv1v2,
                                       0.0, infected_v1v2 * inv_step_size);

        std::cout << "flow_value: " << flow_value << std::endl;

        flows[get_flat_flow_index<InfectionState::InfectedV1V2, InfectionState::RecoveredV1V2>()] =
            std::clamp((1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV1V2] +
                           sqrt((1.0 / params.get<TimeInfectedV2>()) * y[(size_t)InfectionState::InfectedV1V2]) /
                               sqrt(step_size) * iv1v2_rv1v2,
                       0.0, y[(size_t)InfectionState::InfectedV1V2] * inv_step_size);
        std::cout << "flow InfV1V2 -> RecV1V2: "
                  << flows[get_flat_flow_index<InfectionState::InfectedV1V2, InfectionState::RecoveredV1V2>()]
                  << std::endl;

        std::cout << "calc flows done" << std::endl;
    }

    ScalarType step_size; ///< A step size of the model with which the stochastic process is realized.
    mutable RandomNumberGenerator rng;

private:
};

} // namespace sseirvv
} // namespace mio

#endif // MIO_SDE_SEIRVV_MODEL_H
