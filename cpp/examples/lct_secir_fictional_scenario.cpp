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
#include "lct_secir/initialization.h"
#include "ide_secir/infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include <iostream>

int main()
{
    using Vec          = mio::TimeSeries<ScalarType>::Vector;
    int numTransitions = 10;
    ScalarType dt      = 0.5;
    // Set vector that specifies the number of subcompartments
    std::vector<int> subcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
    subcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Exposed]            = 2;
    subcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 3;
    subcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 5;
    mio::lsecir::InfectionState infectionStates(subcompartmentNumbers);

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    mio::TimeSeries<ScalarType> init(numTransitions);

    // add time points for initialization of transitions
    Vec initTransitions(numTransitions);
    initTransitions[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 15.0;
    initTransitions[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
    initTransitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    initTransitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
    initTransitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 1.0;
    initTransitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
    initTransitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 5.0;
    initTransitions[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 1.0;
    initTransitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 1.0;
    initTransitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 1.0;
    initTransitions                                                                              = initTransitions * dt;

    // add initial time point to time series
    init.add_time_point(-40, initTransitions);
    // add further time points until time 0
    while (init.get_last_time() < 0) {
        //initTransitions *=  1.01;
        init.add_time_point(init.get_last_time() + dt, initTransitions);
    }
    mio::lsecir::Initializer initializer(std::move(init), dt, infectionStates);
    auto erg = initializer.compute_initializationvector(100000, 10);
    for (int i = 0; i < infectionStates.get_count(); ++i) {
        std::cout << erg[i] << ", ";
    }
    std::cout << "\n";
}