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

#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"

int main()
{
    std::vector<int> SubcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
    SubcompartmentNumbers[1] = 2;
    SubcompartmentNumbers[2] = 3;
    mio::lsecir::InfectionState InfState(SubcompartmentNumbers);

    ScalarType tmax = 10;
    Eigen::VectorXd init(InfState.get_count());
    init[0]  = 750;
    init[1]  = 30;
    init[2]  = 20;
    init[3]  = 20;
    init[4]  = 10;
    init[5]  = 10;
    init[6]  = 50;
    init[7]  = 50;
    init[8]  = 30;
    init[9]  = 20;
    init[10] = 10;

    mio::lsecir::Model model(std::move(init), InfState);
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(0, tmax, 0.5, model);
    mio::lsecir::print_TimeSeries(result, model.get_heading_Subcompartments());
    mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);
    mio::lsecir::print_TimeSeries(populations, model.get_heading_CompartmentsBase());
}