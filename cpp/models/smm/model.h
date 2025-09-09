/* 
* Copyright (C) 2020-2025 German Aerospace Center (DLR-SC)
*
* Authors: René Schmieding, Julia Bicker
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
#include "smm/parameters.h"
#include "memilio/compartments/compartmental_model.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/geography/regions.h"

namespace mio
{
namespace smm
{

/**
 * @brief Stochastic Metapopulation Model.
 * @tparam regions Number of regions.
 * @tparam Status An infection state enum.
 */
template <size_t regions, class Status, class... Groups>
class Model : public mio::CompartmentalModel<ScalarType, Status,
                                             mio::Populations<ScalarType, mio::regions::Region, Status, Groups...>,
                                             ParametersBase<Status, Groups...>>
{
    using Base = mio::CompartmentalModel<ScalarType, Status,
                                         mio::Populations<ScalarType, mio::regions::Region, Status, Groups...>,
                                         ParametersBase<Status, Groups...>>;

public:
    Model()
        : Base(typename Base::Populations({static_cast<mio::regions::Region>(regions), Status::Count}),
               typename Base::ParameterSet())
    {
    }

    /**
     * @brief Calculate the current rate of the given adoption.
     * @param[in] rate An adoption rate from this model.
     * @param[in] x The current state of the model.
     * @return Current value of the adoption rate.
     */
    ScalarType evaluate(const AdoptionRate<Status, Groups...>& rate, const Eigen::VectorXd& x) const
    {
        const auto& pop   = this->populations;
        const auto source = pop.get_flat_index({rate.region, rate.from}); // Why is here rate.from used? KV
        // determine order and calculate rate
        if (rate.influences.size() == 0) { // first order adoption
            return rate.factor * x[source];
        }
        else { // second order adoption
            ScalarType N = 0;
            for (size_t s = 0; s < static_cast<size_t>(Status::Count); ++s) {
                N += x[pop.get_flat_index({rate.region, Status(s)})];
            }
            // accumulate influences
            ScalarType influences = 0.0;
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
    ScalarType evaluate(const TransitionRate<Status, Groups...>& rate, const Eigen::VectorXd& x) const
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
