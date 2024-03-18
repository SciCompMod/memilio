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

#include "memilio/compartments/flow_model.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "sde_sir/infection_state.h"
#include "sde_sir/parameters.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/type_list.h"
#include <iostream>
namespace mio
{
namespace ssir
{

/********************
    * define the model *
    ********************/
using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Infected>,
                       Flow<InfectionState::Infected,     InfectionState::Recovered>>;

class Model : public FlowModel<InfectionState, Populations<InfectionState>, Parameters, Flows>
{
    using Base = FlowModel<InfectionState, mio::Populations<InfectionState>, Parameters, Flows>;

public:
    Model()
        : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
    {   
    }
    void get_flows(Eigen::Ref<const Eigen::VectorXd> , Eigen::Ref<const Eigen::VectorXd> , double , 
                         Eigen::Ref<Eigen::VectorXd> ) const {}
    void get_flows_stochastic(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t, double dt,
                         Eigen::Ref<Eigen::VectorXd> flows) const 
    {
        auto& params     = this->parameters;
        double coeffStoI = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                           params.get<TransmissionProbabilityOnContact>() / populations.get_total();
        

        RandomNumberGenerator rng = mio::RandomNumberGenerator();
        double si = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);
        double ir = mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(rng, 0.0, 1.0);



        flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>()] =
            coeffStoI * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected]
            + sqrt(coeffStoI * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected]) / sqrt(dt) * si;
        flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>()] = 
            std::min(flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>()], y[(size_t)InfectionState::Susceptible] / dt);

        flows[get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>()] =
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected]
            + sqrt((1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected]) / sqrt(dt) * ir;
        flows[get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>()] =
            std::min(flows[get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>()], y[(size_t)InfectionState::Infected] / dt);    
    }


private:
    
};

} // namespace ssir
} // namespace mio

#endif // ODESIR_MODEL_H
