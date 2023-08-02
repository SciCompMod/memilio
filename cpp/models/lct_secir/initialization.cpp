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
#include "lct_secir/initialization.h"
#include "lct_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/config.h"

namespace mio
{
namespace lsecir
{
Eigen::VectorXd compute_compartments(TimeSeries<ScalarType>&& init, InfectionStateBase Base, InfectionState InfState,
                                     ScalarType gamma, Eigen::Index idx_IncomingFlow, ScalarType dt)
{
    int n = InfState.get_number(Base);
    Eigen::VectorXd initSubcompartments(n);
    ErlangDistribution erlang(n * gamma, 1);

    //initialize relevant parameters
    ScalarType calc_time;
    Eigen::Index calc_time_index;
    ScalarType sum = 0;

    Eigen::Index num_time_points = init.get_num_time_points();
    for (int j = 0; j < n; j++) {
        erlang.set_parameter(j + 1);
        // determine relevant calculation area and corresponding index
        calc_time       = erlang.get_support_max(dt);
        calc_time_index = (Eigen::Index)std::ceil(calc_time / dt) - 1;
        for (Eigen::Index i = num_time_points - calc_time_index; i < num_time_points; i++) {
            ScalarType state_age = (num_time_points - i) * dt;
            sum += erlang.eval(state_age) * init[i][idx_IncomingFlow];
        }
        initSubcompartments[j] = 1 / (n * gamma) * sum;
        sum                    = 0;
    }

    return initSubcompartments;
}
} // namespace lsecir
} // namespace mio