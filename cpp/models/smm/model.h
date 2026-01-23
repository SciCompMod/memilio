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

/// @brief The Index type used to define the SMM subpopulations.
template <class Status, class Region>
using PopulationIndex = decltype(concatenate_indices(std::declval<Region>(), std::declval<Status>()));

/**
 * @brief Stochastic Metapopulation Model.
 * The stratification of the population of this model is split between "Status" and "Region". This split is mostly
 * arbitrary, with the important distinction, that for second order rates the reference population
 * (i.e., the N in S' = S * I / N) is calculated by accumulating subpopulations only over the Status, i.e. individuals 
 * only interact with other individuals within the same Region.
 * Hence, the assumption of homogeneous mixing of the population only holds across Status groups within one Region. 
 * Across Regions, no direct interaction is possible (only indirectly, by first transitioning into another Region)
 * @tparam Comp An enum representing the infection states. Must also be contained in Status
 * @tparam Status A MultiIndex allowing to further stratify infection state adoptions.
 * @tparam Region A MultiIndex for "spatially" distinct subpopulations, default is @ref mio::regions::Region.
 */
template <typename FP, class Comp, class Status = Comp, class Region = mio::regions::Region>
class Model : public mio::CompartmentalModel<FP, Comp, mio::Populations<FP, PopulationIndex<Status, Region>>,
                                             ParametersBase<FP, Status, Region>>
{
    using Base = mio::CompartmentalModel<FP, Comp, mio::Populations<FP, PopulationIndex<Status, Region>>,
                                         ParametersBase<FP, Status, Region>>;
    static_assert(!Base::Populations::Index::has_duplicates, "Do not use the same Index tag multiple times!");

public:
    Model(Status status_dimensions, Region region_dimensions)
        : Base(typename Base::Populations(concatenate_indices(region_dimensions, status_dimensions)),
               typename Base::ParameterSet())
    {
    }

    /**
     * @brief Calculate the current rate of the given adoption.
     * @param[in] rate An adoption rate from this model.
     * @param[in] x The current state of the model.
     * @return Current value of the adoption rate.
     */
    FP evaluate(const AdoptionRate<FP, Status, Region>& rate, const Eigen::VectorX<FP>& x) const
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
    FP evaluate(const TransitionRate<FP, Status, Region>& rate, const Eigen::VectorX<FP>& x) const
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
