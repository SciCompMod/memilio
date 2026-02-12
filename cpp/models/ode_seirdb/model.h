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
#ifndef SEIRDB_MODEL_H
#define SEIRDB_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/math/interpolation.h"
#include "memilio/utils/time_series.h"
#include "ode_seirdb/infection_state.h"
#include "ode_seirdb/parameters.h"

GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wshadow")
#include <Eigen/Dense>
GCC_CLANG_DIAGNOSTIC(pop)

namespace mio
{
namespace oseirdb
{

/********************
 * define the model *
 ********************/

// clang-format off
using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Exposed>,
                       Flow<InfectionState::Exposed, InfectionState::Infected>,
                       Flow<InfectionState::Infected, InfectionState::Recovered>,
                       Flow<InfectionState::Infected, InfectionState::Dead>,
                       Flow<InfectionState::Dead, InfectionState::Buried>>;
// clang-format on
template <typename FP>
class Model
    : public FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>
{
    using Base = FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>;

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

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const Index<AgeGroup> age_groups = reduce_index<Index<AgeGroup>>(this->populations.size());
        const auto& params               = this->parameters;

        for (auto i : make_index_range(age_groups)) {
            const size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            const size_t Ei = this->populations.get_flat_index({i, InfectionState::Exposed});
            const size_t Ii = this->populations.get_flat_index({i, InfectionState::Infected});
            const size_t Di = this->populations.get_flat_index({i, InfectionState::Dead});

            for (auto j : make_index_range(age_groups)) {
                const size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
                const size_t Ej = this->populations.get_flat_index({j, InfectionState::Exposed});
                const size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
                const size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});
                const size_t Dj = this->populations.get_flat_index({j, InfectionState::Dead});
                const size_t Bj = this->populations.get_flat_index({j, InfectionState::Buried});

                const FP Nj           = pop[Sj] + pop[Ej] + pop[Ij] + pop[Rj] + pop[Dj] + pop[Bj];
                const FP divNj        = (Nj < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nj);
                const FP contact_rate = params.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(
                                            SimulationTime<FP>(t))(i.get(), j.get()) *
                                        divNj;
                const FP coeffStoE     = contact_rate * params.template get<TransmissionProbabilityOnContact<FP>>()[i];
                const FP coeffStoEDead = contact_rate * params.template get<TransmissionProbabilityFromDead<FP>>()[i];

                flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>(i)] +=
                    y[Si] * (coeffStoE * pop[Ij] + coeffStoEDead * pop[Dj]);
            }
            flows[Base::template get_flat_flow_index<InfectionState::Exposed, InfectionState::Infected>(i)] =
                (1.0 / params.template get<TimeExposed<FP>>()[i]) * y[Ei];

            const FP inv_time_infected = 1.0 / params.template get<TimeInfected<FP>>()[i];
            const FP prob_recover      = params.template get<ProbabilityToRecover<FP>>()[i];
            const FP prob_die          = 1.0 - prob_recover;
            flows[Base::template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>(i)] =
                inv_time_infected * prob_recover * y[Ii];
            flows[Base::template get_flat_flow_index<InfectionState::Infected, InfectionState::Dead>(i)] =
                inv_time_infected * prob_die * y[Ii];
            flows[Base::template get_flat_flow_index<InfectionState::Dead, InfectionState::Buried>(i)] =
                (1.0 / params.template get<TimeToBurial<FP>>()[i]) * y[Di];
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

} // namespace oseirdb
} // namespace mio

#endif // SEIRDB_MODEL_H
