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

#include "memilio/utils/index.h"
#include "memilio/utils/index_range.h"
#include "smm/parameters.h"
#include "memilio/compartments/compartmental_model.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/geography/regions.h"
#include <utility>

namespace mio
{
namespace smm
{

template <class T, std::size_t>
using age_group = T;

template <class Status, class Region>
using PopulationIndex = decltype(merge_indices(std::declval<Region>(), std::declval<Status>()));

/**
 * @brief Stochastic Metapopulation Model.
 * @tparam regions Number of regions.
 * @tparam Status An infection state enum.
 */
template <typename FP, class Comp, class StatusT = Comp, class RegionT = mio::regions::Region>
class Model : public mio::CompartmentalModel<FP, Comp, mio::Populations<FP, PopulationIndex<StatusT, RegionT>>,
                                             ParametersBase<FP, StatusT, RegionT>>
{
    using Base = mio::CompartmentalModel<FP, Comp, mio::Populations<FP, PopulationIndex<StatusT, RegionT>>,
                                         ParametersBase<FP, StatusT, RegionT>>;
    static_assert(!Base::Populations::Index::has_duplicates, "Do not use the same Index tag multiple times!");

public:
    using Status = StatusT;
    using Region = RegionT;

    Model(Status status_dimensions, Region region_dimensions)
        : Base(typename Base::Populations(merge_indices(region_dimensions, status_dimensions)),
               typename Base::ParameterSet())
    {
    }

    /**
     * @brief Calculate the current rate of the given adoption.
     * @param[in] rate An adoption rate from this model.
     * @param[in] x The current state of the model.
     * @return Current value of the adoption rate.
     */
    FP evaluate(const AdoptionRate<FP, Status, Region>& rate, const Eigen::VectorXd& x) const
    {
        const auto& pop   = this->populations;
        const auto source = pop.get_flat_index({rate.region, rate.from});
        // determine order and calculate rate
        if (rate.influences.size() == 0) { // first order adoption
            return rate.factor * x[source];
        }
        else { // second order adoption
            FP N = 0;
            for (auto status : make_index_range(reduce_index<Status>(this->populations.size()))) {
                N += x[pop.get_flat_index({rate.region, status})];
            }
            // accumulate influences
            FP influences = 0.0;
            for (size_t i = 0; i < rate.influences.size(); i++) {
                influences +=
                    rate.influences[i].factor * x[pop.get_flat_index({rate.region, rate.influences[i].status})];
            }
            return (N > 0) ? (rate.factor * x[source] * influences / N) : 0;
        }
    }

    /**
     * @brief Calculate the current rate of the given spatial transition.
     * @param[in] rate A transition rate from this model.
     * @param[in] x The current state of the model.
     * @return Current value of the transition rate.
     */
    FP evaluate(const TransitionRate<FP, Status, Region>& rate, const Eigen::VectorXd& x) const
    {
        const auto source = this->populations.get_flat_index({rate.from, rate.status});
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
