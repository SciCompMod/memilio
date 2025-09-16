/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef ODE_MSEIRS4_MODEL_H
#define ODE_MSEIRS4_MODEL_H

#include "memilio/compartments/compartmental_model.h"
#include "memilio/epidemiology/populations.h"
#include "ode_mseirs4/infection_state.h"
#include "ode_mseirs4/parameters.h"

namespace mio
{
namespace omseirs4
{

template <typename FP = ScalarType>
class Model : public mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>>
{
    using Base = mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;

    Model()
        : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                         Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
    {
        auto& params   = this->parameters;
        const FP b0    = params.template get<BaseTransmissionRate<FP>>();
        const FP b1    = params.template get<SeasonalAmplitude<FP>>();
        const FP phi   = params.template get<SeasonalPhase<FP>>();
        const FP mu    = params.template get<NaturalBirthDeathRate<FP>>();
        const FP xi    = params.template get<LossMaternalImmunityRate<FP>>();
        const FP sigma = params.template get<ProgressionRate<FP>>();
        const FP nu    = params.template get<RecoveryRate<FP>>();
        const FP gamma = params.template get<ImmunityWaningRate<FP>>();
        const FP f2    = params.template get<Beta2Factor<FP>>();
        const FP f3    = params.template get<Beta3Factor<FP>>();
        const FP f4    = params.template get<Beta4Factor<FP>>();

        const FP two_pi = FP(2) * std::numbers::pi_v<FP>;
        const FP beta1  = b0 * (FP(1) + b1 * std::cos(two_pi * (t / FP(365.0)) + phi));
        const FP beta2  = f2 * beta1;
        const FP beta3  = f3 * beta1;
        const FP beta4  = f4 * beta1;

        auto idx = [&](InfectionState s) {
            return static_cast<size_t>(s);
        };

        FP I_total = y[idx(InfectionState::I1)] + y[idx(InfectionState::I2)] + y[idx(InfectionState::I3)] +
                     y[idx(InfectionState::I4)];
        FP R_total = y[idx(InfectionState::R1)] + y[idx(InfectionState::R2)] + y[idx(InfectionState::R3)] +
                     y[idx(InfectionState::R4)];
        FP N = pop.sum();

        // dM
        dydt[idx(InfectionState::MaternalImmune)] = mu * R_total - (xi + mu) * y[idx(InfectionState::MaternalImmune)];

        // dS1
        dydt[idx(InfectionState::S1)] = mu * (N - R_total) + xi * y[idx(InfectionState::MaternalImmune)] -
                                        mu * y[idx(InfectionState::S1)] - beta1 * I_total * y[idx(InfectionState::S1)];

        // dE1..E4
        dydt[idx(InfectionState::E1)] =
            beta1 * I_total * y[idx(InfectionState::S1)] - (mu + sigma) * y[idx(InfectionState::E1)];
        dydt[idx(InfectionState::E2)] =
            beta2 * I_total * y[idx(InfectionState::S2)] - (mu + sigma) * y[idx(InfectionState::E2)];
        dydt[idx(InfectionState::E3)] =
            beta3 * I_total * y[idx(InfectionState::S3)] - (mu + sigma) * y[idx(InfectionState::E3)];
        dydt[idx(InfectionState::E4)] =
            beta4 * I_total * y[idx(InfectionState::S4)] - (mu + sigma) * y[idx(InfectionState::E4)];

        // dI1..I4
        dydt[idx(InfectionState::I1)] = sigma * y[idx(InfectionState::E1)] - (nu + mu) * y[idx(InfectionState::I1)];
        dydt[idx(InfectionState::I2)] = sigma * y[idx(InfectionState::E2)] - (nu + mu) * y[idx(InfectionState::I2)];
        dydt[idx(InfectionState::I3)] = sigma * y[idx(InfectionState::E3)] - (nu + mu) * y[idx(InfectionState::I3)];
        dydt[idx(InfectionState::I4)] = sigma * y[idx(InfectionState::E4)] - (nu + mu) * y[idx(InfectionState::I4)];

        // dR1..R4
        dydt[idx(InfectionState::R1)] = nu * y[idx(InfectionState::I1)] - (mu + gamma) * y[idx(InfectionState::R1)];
        dydt[idx(InfectionState::R2)] = nu * y[idx(InfectionState::I2)] - (mu + gamma) * y[idx(InfectionState::R2)];
        dydt[idx(InfectionState::R3)] = nu * y[idx(InfectionState::I3)] - (mu + gamma) * y[idx(InfectionState::R3)];
        dydt[idx(InfectionState::R4)] = nu * y[idx(InfectionState::I4)] - (mu + gamma) * y[idx(InfectionState::R4)];

        // dS2,S3,S4
        dydt[idx(InfectionState::S2)] = gamma * y[idx(InfectionState::R1)] - mu * y[idx(InfectionState::S2)] -
                                        beta2 * I_total * y[idx(InfectionState::S2)];
        dydt[idx(InfectionState::S3)] = gamma * y[idx(InfectionState::R2)] - mu * y[idx(InfectionState::S3)] -
                                        beta3 * I_total * y[idx(InfectionState::S3)];
        dydt[idx(InfectionState::S4)] = gamma * (y[idx(InfectionState::R3)] + y[idx(InfectionState::R4)]) -
                                        mu * y[idx(InfectionState::S4)] - beta4 * I_total * y[idx(InfectionState::S4)];
    }
};

} // namespace omseirs4
} // namespace mio

#endif // ODE_MSEIRS4_MODEL_H
