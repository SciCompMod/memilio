/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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

#ifndef LCT_SECIR_MODEL_H
#define LCT_SECIR_MODEL_H

#include "lct_secir/parameters.h"
#include "lct_secir/infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"



namespace mio
{
namespace lsecir
{
class Model: public CompartmentalModel<InfectionState, Populations<InfectionState>, Parameters>
{
    using ParameterSet = Parameters;

public:
    Model();

    /**
    * @brief Checks constraints on model parameters.
    */
    void check_constraints() const;

    void eval_right_hand_side();

    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    InfectionState InfectionStates;
    /* Attention: m_populations and m_transitions do not necessarily have the same number of time points due to the initialization part. */
    TimeSeries<ScalarType>
        m_populations; ///< TimeSeries containing points of time and the corresponding number of people in defined InfectionState%s.

private:
    ScalarType m_N{0}; ///< Total population size of the considered region.
};

} // namespace lsecir
} // namespace mio

#endif // LCTSECIR_MODEL_H
