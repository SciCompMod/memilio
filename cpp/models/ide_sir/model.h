/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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

#ifndef IDESIR_MODEL_H
#define IDESIR_MODEL_H

#include "ide_sir/parameters.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isir
{

class Model
{
    using ParameterSet = Parameters;

public:
    /**
    * @brief Constructor to create an IDE-SIR model.
    *
    * @param[in, out] transitions_init TimeSeries with the initial values of the number of individuals, 
    *   which transit within one timestep dt from one compartment to another.
    *   Possible transitions are specified in #InfectionTransition%s.
    *   Considered time points should have the distance dt. The last time point determines the start time t0 of the 
    *   simulation. 
    * @param[in] N_init Total population.
    */
    Model(TimeSeries<ScalarType>&& populations_init, ScalarType N_init);

    ScalarType get_totalpop() const;

    // ---- Public parameters. ----
    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    // Attention: populations and transitions do not necessarily have the same number of time points due to the
    // initialization part.
    TimeSeries<ScalarType> populations; ///< TimeSeries containing points of time and the corresponding number of
        // people in defined #InfectionState%s for every AgeGroup.

private:
    // ---- Private parameters. ----
    ScalarType m_N; ///< Vector containing the total population size of the considered region for every AgeGroup.
};

} // namespace isir
} // namespace mio

#endif // IDESIR_MODEL_H
