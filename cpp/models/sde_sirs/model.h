/* 
* Copyright (C) 2020-2024 MEmilio
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

#ifndef SDESIR_MODEL_H
#define SDESIR_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "sde_sirs/infection_state.h"
#include "sde_sirs/parameters.h"
#include "memilio/utils/random_number_generator.h"
#include <iostream>
namespace mio
{
namespace ssirs
{

/********************
    * define the model *
    ********************/

class Model : public CompartmentalModel<InfectionState, Populations<InfectionState>, Parameters>
{
    using Base = CompartmentalModel<InfectionState, mio::Populations<InfectionState>, Parameters>;

public:
    Model()
        : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
    {   
    }

    void get_derivatives_stoch(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt, double dt) const 
    {
        auto& params     = this->parameters;
        double coeffStoI = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                           params.get<TransmissionProbabilityOnContact>() / populations.get_total();
        
        RandomNumberGenerator rng = mio::RandomNumberGenerator();
        double si = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double ir = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double rs = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);                
        
        dydt[(size_t)InfectionState::Susceptible] =
            -coeffStoI * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected]
            + (1.0 / params.get<TimeImmune>()) * y[(size_t)InfectionState::Recovered] 
            - sqrt(coeffStoI * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected]) / sqrt(dt) * si
            + sqrt((1.0 / params.get<TimeImmune>()) * y[(size_t)InfectionState::Recovered]) / sqrt(dt) * rs;
        dydt[(size_t)InfectionState::Infected] =
            coeffStoI * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected]
            - (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected] 
            + sqrt(coeffStoI * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected]) / sqrt(dt) * si
            - sqrt((1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected]) / sqrt(dt) * ir; 
        dydt[(size_t)InfectionState::Recovered] =
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected]
            - (1.0 / params.get<TimeImmune>()) * y[(size_t)InfectionState::Recovered]
            + sqrt((1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected]) / sqrt(dt) * ir
            - sqrt((1.0 / params.get<TimeImmune>()) * y[(size_t)InfectionState::Recovered]) / sqrt(dt) * rs;
    }

private:
    
};

} // namespace ssirs
} // namespace mio

#endif // SDESIRS_MODEL_H
