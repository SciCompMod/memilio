/*
* Copyright (C) 2020-2026 MEmilio
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

#include "memilio/compartments/flow_model.h"
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

using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Infected>,
                       Flow<InfectionState::Infected, InfectionState::Recovered>>;

template <typename FP>
class Model
    : public mio::FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>
{
    using Base =
        mio::FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;

    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    Model(int num_agegroups)
        : Base(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    // void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
    //                      Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
    // {
    //     auto params                                  = this->parameters;
    //     AgeGroup n_agegroups                         = params.get_num_groups();
    //     ContactMatrixGroup<FP> const& contact_matrix = params.template get<ContactPatterns<FP>>();

    //     for (auto i = AgeGroup(0); i < n_agegroups; i++) {

    //         size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
    //         size_t Ii = this->populations.get_flat_index({i, InfectionState::Infected});
    //         size_t Ri = this->populations.get_flat_index({i, InfectionState::Recovered});

    //         for (auto j = AgeGroup(0); j < n_agegroups; j++) {

    //             size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
    //             size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
    //             size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

    //             const FP Nj    = pop[Sj] + pop[Ij] + pop[Rj];
    //             const FP divNj = (Nj < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nj);

    //             FP coeffStoI = contact_matrix.get_matrix_at(SimulationTime<FP>(t))(
    //                                static_cast<Eigen::Index>((size_t)i), static_cast<Eigen::Index>((size_t)j)) *
    //                            params.template get<TransmissionProbabilityOnContact<FP>>()[i] * divNj;

    //             dydt[Si] += -coeffStoI * y[Si] * pop[Ij];
    //             dydt[Ii] += coeffStoI * y[Si] * pop[Ij];
    //         }
    //         dydt[Ii] -= (1.0 / params.template get<TimeInfected<FP>>()[i]) * y[Ii];
    //         dydt[Ri] = (1.0 / params.template get<TimeInfected<FP>>()[i]) * y[Ii];
    //     }
    // }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const Index<AgeGroup> age_groups = reduce_index<Index<AgeGroup>>(this->populations.size());
        const auto& params               = this->parameters;

        for (auto i : make_index_range(age_groups)) {
            const size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            const size_t Ii = this->populations.get_flat_index({i, InfectionState::Infected});

            for (auto j : make_index_range(age_groups)) {
                const size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
                const size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
                const size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

                const FP Nj        = pop[Sj] + pop[Ij] + pop[Rj];
                const FP divNj     = (Nj < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nj);
                const FP coeffStoI = params.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(
                                         SimulationTime<FP>(t))(i.get(), j.get()) *
                                     params.template get<TransmissionProbabilityOnContact<FP>>()[i] * divNj;

                flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(i)] +=
                    coeffStoI * y[Si] * pop[Ij];
            }
            flows[Base::template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>(i)] =
                (1.0 / params.template get<TimeInfected<FP>>()[i]) * y[Ii];
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
        obj.add_element("Parameters", this->parameters);
        obj.add_element("Populations", this->populations);
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
