/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler, Lena Ploetzke
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
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>

int main()
{   
    std::vector<int> SubcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
    SubcompartmentNumbers[1]=2;
    SubcompartmentNumbers[2]=3;
    mio::lsecir::InfectionState InfState(SubcompartmentNumbers);

    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax        = 10;
    ScalarType N           = 1000;
    mio::TimeSeries<ScalarType> init(InfState.Count);
    Vec vec_init(InfState.Count);
    vec_init[0]=750;
    vec_init[1]=30;
    vec_init[2]=20;
    vec_init[3]=20;
    vec_init[4]=10;
    vec_init[5]=10;
    vec_init[6]=50;
    vec_init[7]=50;
    vec_init[8]=30;
    vec_init[9]=20;
    vec_init[10]=10;

    init.add_time_point(0, vec_init);

    mio::lsecir::Model model(std::move(init),N,InfState);
    mio::lsecir::simulate(0, tmax, 0.5, model);
}