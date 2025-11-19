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

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/stochastic_model.h"
#include "memilio/epidemiology/populations.h"
#include "sde_sirs/infection_state.h"
#include "sde_sirs/parameters.h"
#include <cmath>

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
        // FP season_val =
        //     (1 + params.template get<Seasonality<FP>>() *
        //              sin(std::numbers::pi_v<ScalarType> * ((params.template get<StartDay<FP>>() + t) / 182.5 + 0.5)));
        FP t_season =
            params.template get<SeasonalityPeak<FP>>() - params.template get<StartDay<FP>>() + int(t / 365.) * 365.;
        FP season_end = t_season + 182.5;
        if (t >= season_end) {
            t_season = params.template get<SeasonalityPeak<FP>>() - params.template get<StartDay<FP>>() +
                       (int(t / 365.) + 1) * 365.;
        }
        FP season_val = params.template get<Seasonality<FP>>() + (1 - params.template get<Seasonality<FP>>()) *
                                                                     (gaussian(t - t_season) / gaussian(0.)) * zeta(t);

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

    double gaussian(FP t) const
    {
        auto& params = Base::parameters;
        return 1. /
               (std::sqrt(2. * std::numbers::pi_v<ScalarType> * params.template get<SeasonalitySigma<FP>>() *
                          params.template get<SeasonalitySigma<FP>>())) *
               std::exp(-0.5 * std::pow(t / params.template get<SeasonalitySigma<FP>>(), 2));
    }

    double zeta(FP t) const
    {
        auto& params     = Base::parameters;
        double tau_minus = params.template get<SeasonalityPeak<FP>>() - params.template get<StartDay<FP>>() +
                           365. * int(t / 365.) - 182.5;
        double tau_plus = params.template get<SeasonalityPeak<FP>>() - params.template get<StartDay<FP>>() +
                          365. * int(t / 365.) - 182.5;
        if (t <= params.template get<SeasonalityPeak<FP>>() - params.template get<StartDay<FP>>() - 182.5 + 30) {
            return (t - tau_minus + 30) / 60.;
        }
        else if ((tau_plus - 30 <= t) && (t <= tau_plus)) {
            return (-t + tau_plus + 30) / 60.;
        }
        else if ((tau_plus <= t) && (t <= tau_plus + 30)) {
            return (t - tau_plus + 30) / 60.;
        }
        else {
            return 1.;
        }
    }
};

} // namespace ssirs
} // namespace mio

#endif // MIO_SDE_SIRS_MODEL_H
