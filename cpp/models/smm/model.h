/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: René Schmieding, Julia Bicker, Kilian Volmer
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

#ifndef MIO_SMM_MODEL_H
#define MIO_SMM_MODEL_H

#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "smm/parameters.h"
#include "memilio/compartments/compartmental_model.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/geography/regions.h"

namespace mio
{
namespace smm
{

template <class T, std::size_t>
using age_group = T;

/**
 * @brief Stochastic Metapopulation Model.
 * @tparam regions Number of regions.
 * @tparam Status An infection state enum.
 */
template <typename FP, size_t regions, class Status, size_t... Groups>
class Model : public mio::CompartmentalModel<
                  FP, Status, mio::Populations<FP, mio::regions::Region, Status, age_group<mio::AgeGroup, Groups>...>,
                  ParametersBase<FP, Status, age_group<mio::AgeGroup, Groups>...>>
{
    using Base =
        mio::CompartmentalModel<FP, Status,
                                mio::Populations<FP, mio::regions::Region, Status, age_group<mio::AgeGroup, Groups>...>,
                                ParametersBase<FP, Status, age_group<mio::AgeGroup, Groups>...>>;
    using Index = mio::Index<mio::regions::Region, Status, age_group<mio::AgeGroup, Groups>...>;

public:
    Model()
        : Base(typename Base::Populations(
                   {static_cast<mio::regions::Region>(regions), Status::Count, static_cast<mio::AgeGroup>(Groups)...}),
               typename Base::ParameterSet())
    {
    }

    /**
     * @brief Calculate the current rate of the given adoption.
     * @param[in] rate An adoption rate from this model.
     * @param[in] x The current state of the model.
     * @return Current value of the adoption rate.
     */
    FP evaluate(const AdoptionRate<FP, Status, age_group<mio::AgeGroup, Groups>...>& rate,
                const Eigen::VectorXd& x) const
    {
        const auto& pop       = this->populations;
        const auto index_from = std::apply(
            [&](auto&&... args) {
                return Index{rate.region, rate.from, std::forward<decltype(args)>(args)...};
            },
            rate.group_indices);
        const auto source = pop.get_flat_index(index_from); // Why is here rate.from used? KV
        // determine order and calculate rate
        if (rate.influences.size() == 0) { // first order adoption
            return rate.factor * x[source];
        }
        else { // second order adoption
            // accumulate influences
            FP influences = 0.0;
            for (size_t i = 0; i < rate.influences.size(); i++) {
                FP N = 0; // Welches N brauchen wir hier??

                for (size_t s = 0; s < static_cast<size_t>(Status::Count); ++s) {
                    const auto index = std::apply(
                        [&](auto&&... args) {
                            return Index{rate.influences[i].region.value_or(rate.region), Status(s),
                                         std::forward<decltype(args)>(args)...};
                        },
                        rate.influences[i].group_indices);
                    N += x[pop.get_flat_index(index)];
                }
                const auto index = std::apply(
                    [&](auto&&... args) {
                        return Index{rate.influences[i].region.value_or(rate.region), rate.influences[i].status,
                                     std::forward<decltype(args)>(args)...};
                    },
                    rate.influences[i].group_indices);
                if (N > 0) {
                    influences += rate.influences[i].factor * x[pop.get_flat_index(index)] / N;
                }
            }
            return rate.factor * x[source] * influences;
        }
    }

    /**
     * @brief Calculate the current rate of the given spatial transition.
     * @param[in] rate A transition rate from this model.
     * @param[in] x The current state of the model.
     * @return Current value of the transition rate.
     */
    FP evaluate(const TransitionRate<FP, Status, age_group<mio::AgeGroup, Groups>...>& rate,
                const Eigen::VectorXd& x) const
    {
        auto index = std::apply(
            [&](auto&&... args) {
                return Index{rate.from, rate.status, std::forward<decltype(args)>(args)...};
            },
            rate.group_indices_from);
        const auto source = this->populations.get_flat_index(index);
        return rate.factor * x[source];
    }

    /**
    * Get the RandomNumberGenerator used by this Model for random events.
    * @return The random number generator.
    */
    mio::RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

private:
    mio::RandomNumberGenerator m_rng; ///< Model's random number generator.
};

} //namespace smm

} // namespace mio

#endif
