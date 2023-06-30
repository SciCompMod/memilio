/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef SEAIR_MODEL_H
#define SEAIR_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h" // IWYU pragma: keep
#include "ode_seair/infection_state.h"
#include "ode_seair/parameters.h"
#include "ad/ad.hpp"

namespace mio
{
namespace oseair
{

/********************
    * define the model *
    ********************/
template<typename FP=ad::gt1s<double>::type>
class Model : public CompartmentalModel<InfectionState, Populations<InfectionState>,
                                        Parameters, FP>
{
    using Base = CompartmentalModel<InfectionState, mio::Populations<InfectionState>,
                                    Parameters, FP>;

public:
    Model()
        : Base(Populations<FP>({InfectionState::Count}, 0.), ParameterSet<FP>())
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::Matrix<ad::gt1s<double>::type,Eigen::Dynamic,1>> /* pop */,
                         Eigen::Ref<const Eigen::Matrix<ad::gt1s<double>::type,Eigen::Dynamic,1>> y, ad::gt1s<double>::type /* t */,
                         Eigen::Ref<Eigen::Matrix<ad::gt1s<double>::type,Eigen::Dynamic,1>> dydt) const override
    {
        auto& params     = this->parameters;
        /*        double coeffStoE = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                                   params.get<TransmissionProbabilityOnContact>() / populations.get_total();

        dydt[(size_t)InfectionState::Susceptible] =
            -coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Exposed] =
            coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected] -
            (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed];
        dydt[(size_t)InfectionState::Infected] =
            (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed] -
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Recovered] =
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected]; */

        auto& alpha_a = params.template get<AlphaA>();
        auto& alpha_i = params.template get<AlphaI>();
        auto& kappa = params.template get<Kappa>();
        auto& beta = params.template get<Beta>();
        auto& mu = params.template get<Mu>();
        auto& t_latent_inverse = params.template get<TLatentInverse>();
        auto& rho = params.template get<Rho>();
        auto& gamma = params.template get<Gamma>();



        const auto& s = y[(size_t)InfectionState::Susceptible];
        const auto& e = y[(size_t)InfectionState::Exposed];
        const auto& a = y[(size_t)InfectionState::Asymptomatic];
        const auto& i = y[(size_t)InfectionState::Infected];
        const auto& r = y[(size_t)InfectionState::Recovered];



        dydt[(size_t)InfectionState::Susceptible] = -alpha_a * s * a - alpha_i * s * i + gamma * r;
        dydt[(size_t)InfectionState::Exposed] = alpha_a  * s * a + alpha_i * s * i - t_latent_inverse * e;
        dydt[(size_t)InfectionState::Asymptomatic] = t_latent_inverse * e - kappa * a - rho * a;
        dydt[(size_t)InfectionState::Infected] = kappa * a - beta * i - mu * i;
        dydt[(size_t)InfectionState::Recovered] = rho * a + beta * i - gamma * r;
        dydt[(size_t)InfectionState::Perished] = mu * i;
        dydt[(size_t)InfectionState::ObjectiveFunction] = -alpha_i  -alpha_a + 0.1 * kappa;
    }
};

} // namespace oseir
} // namespace mio

#endif // SEAIR_MODEL_H
