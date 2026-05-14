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

#ifndef MIO_SDE_SIR_MODEL_H
#define MIO_SDE_SIR_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/stochastic_model.h"
#include "memilio/epidemiology/populations.h"
#include "sde_sirs/infection_state.h"
#include "sde_sirs/parameters.h"

namespace mio
{
namespace ssirs
{

/********************
 * define the model *
 ********************/

using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Infected>,
                       Flow<InfectionState::Infected, InfectionState::Recovered>,
                       Flow<InfectionState::Recovered, InfectionState::Susceptible>>;

template <typename FP>
class Model
    : public mio::StochasticModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>, Flows>
{
public:
    using Base = mio::StochasticModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>, Flows>;

    Model()
        : Base(typename Base::Populations({InfectionState::Count}, 0.0), typename Base::ParameterSet())
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const
    {
        auto& params = Base::parameters;
        // effective contact rate by contact rate between groups i and j and damping j
        FP season_val =
            (1 + params.template get<Seasonality<FP>>() *
                     sin(std::numbers::pi_v<ScalarType> * ((params.template get<StartDay<FP>>() + t) / 182.5 + 0.5)));

        FP coeffStoI = season_val *
                       params.template get<ContactPatterns<FP>>().get_matrix_at(SimulationTime<FP>(t))(0, 0) *
                       params.template get<TransmissionProbabilityOnContact<FP>>() / Base::populations.get_total();

        flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>()] =
            coeffStoI * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
        flows[this->template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>()] =
            (1.0 / params.template get<TimeInfected<FP>>()) * y[(size_t)InfectionState::Infected];
        flows[this->template get_flat_flow_index<InfectionState::Recovered, InfectionState::Susceptible>()] =
            (1.0 / params.template get<TimeImmune<FP>>()) * y[(size_t)InfectionState::Recovered];
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

} // namespace ssirs
} // namespace mio

#endif // MIO_SDE_SIRS_MODEL_H
