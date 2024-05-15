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

#ifndef ODESIR_MODEL_H
#define ODESIR_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/parameters.h"

namespace mio
{
namespace osir
{

/********************
    * define the model *
    ********************/

class Model : public CompartmentalModel<InfectionState, Populations<AgeGroup, InfectionState>, Parameters>
{
    using Base = CompartmentalModel<InfectionState, mio::Populations<AgeGroup, InfectionState>, Parameters>;

public:
    Model(int num_agegroups)
        : Base(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        auto params                              = this->parameters;
        AgeGroup n_agegroups                     = params.get_num_groups();
        ContactMatrixGroup const& contact_matrix = params.get<ContactPatterns>();

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            size_t Ii = this->populations.get_flat_index({i, InfectionState::Infected});
            size_t Ri = this->populations.get_flat_index({i, InfectionState::Recovered});

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {

                size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
                size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
                size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

                double Nj = pop[Sj] + pop[Ij] + pop[Rj];

                double coeffStoI = contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                   static_cast<Eigen::Index>((size_t)j)) *
                                   params.get<TransmissionProbabilityOnContact>()[i] / Nj;

                dydt[Si] += -coeffStoI * y[Si] * pop[Ij];
                dydt[Ii] += coeffStoI * y[Si] * pop[Ij];
            }
            dydt[Ii] -= (1.0 / params.get<TimeInfected>()[i]) * y[Ii];
            dydt[Ri] = (1.0 / params.get<TimeInfected>()[i]) * y[Ii];

        }
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Model");
        obj.add_element("Parameters", parameters);
        obj.add_element("Populations", populations);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Model> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Model");
        auto par = obj.expect_element("Parameters", Tag<ParameterSet>{});
        auto pop = obj.expect_element("Populations", Tag<Populations>{});
        return apply(
            io,
            [](auto&& par_, auto&& pop_) {
                return Model{pop_, par_};
            },
            par, pop);
    }
};

} // namespace osir
} // namespace mio

#endif // ODESIR_MODEL_H
