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
#include <string>

namespace mio
{
namespace lsecir
{
class Model
{
    using ParameterSet = Parameters;

public:
    Model(Eigen::VectorXd init, const InfectionState InfectionState_init = InfectionState(),
          const ParameterSet& Parameterset_init = ParameterSet());

    /**
    * @brief Checks constraints on model parameters.
    */
    void check_constraints() const;

    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorXd> y, ScalarType t,
                              Eigen::Ref<Eigen::VectorXd> dydt) const;

    TimeSeries<ScalarType> calculate_populations(const TimeSeries<ScalarType>& result) const;

    std::string get_heading_CompartmentsBase() const;
    std::string get_heading_Subcompartments() const;

    Eigen::VectorXd get_initial_values()
    {
        return m_initial_values;
    }

    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    InfectionState m_InfectionStates;

private:
    Eigen::VectorXd m_initial_values;
    ScalarType m_N{0}; ///< Total population size of the considered region.
};

} // namespace lsecir
} // namespace mio

#endif // LCTSECIR_MODEL_H
