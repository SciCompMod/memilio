/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#ifndef SEIR_MODEL_H
#define SEIR_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "seir/infection_state.h"
#include "seir/parameters.h"

namespace mio
{
namespace seir
{

/********************
 * define the model *
 ********************/

class Model : public CompartmentalModel<Populations<InfectionState>, ParametersBase>
{
    using Base = CompartmentalModel<mio::Populations<InfectionState>, ParametersBase>;
    using Po = Base::Populations;
    using Pa = Base::ParameterSet;

public:
    Model()
        : Base(Po({Index<InfectionState>((size_t)InfectionState::Count)}, 0.), Pa())
    {
#if !USE_DERIV_FUNC
        //S to E
        this->add_flow(std::make_tuple(InfectionState::S), std::make_tuple(InfectionState::E),
                       [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           ScalarType cont_freq_eff = params.contact_frequency.get_matrix_at(t)(0, 0);

                           //TODO: We should probably write a static Po::get_total_from function
                           ScalarType divN = 1.0 / (Po::get_from(y, InfectionState::S) + Po::get_from(y, InfectionState::E) +
                                                    Po::get_from(y, InfectionState::I) + Po::get_from(y, InfectionState::R));
                           return cont_freq_eff * Po::get_from(y, InfectionState::S) * Po::get_from(pop, InfectionState::I) *
                                  divN;
                       });

        //E to I
        this->add_flow(std::make_tuple(InfectionState::E), std::make_tuple(InfectionState::I),
                       [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                           return p.times.get_incubation_inv() * Po::get_from(y, InfectionState::E);
                       });

        //I to R
        this->add_flow(std::make_tuple(InfectionState::I), std::make_tuple(InfectionState::R),
                       [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                           return p.times.get_infectious_inv() * Po::get_from(y, InfectionState::I);
                       });
#endif
    }

#if USE_DERIV_FUNC
    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop,
                         Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        auto& params = this->parameters;
        double S2E_coeff = params.get<ContactFrequency>().get_matrix_at(t)(0, 0) * params.get<TransmissionRisk>()
                / populations.get_total();

        dydt[(size_t)InfectionState::S] = -S2E_coeff * y[(size_t)InfectionState::S] * pop[(size_t)InfectionState::I];
        dydt[(size_t)InfectionState::E] = S2E_coeff * y[(size_t)InfectionState::S] * pop[(size_t)InfectionState::I] -
                                    params.get<StageTimeIncubationInv>() * y[(size_t)InfectionState::E];
        dydt[(size_t)InfectionState::I] = params.get<StageTimeIncubationInv>() * y[(size_t)InfectionState::E] -
                                    params.get<StageTimeInfectiousInv>() * y[(size_t)InfectionState::I];
        dydt[(size_t)InfectionState::R] = params.get<StageTimeInfectiousInv>() * y[(size_t)InfectionState::I];
    }

#endif // USE_DERIV_FUNC
};

} // namespace seir
} // namespace mio

#endif // SEIR_MODEL_H
