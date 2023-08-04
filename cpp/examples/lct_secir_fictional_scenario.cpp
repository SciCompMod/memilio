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
#include "lct_secir/parameters.h"
#include "lct_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/time_series.h"
#include <iostream>

int main()
{
    bool save_result   = true;
    using Vec          = mio::TimeSeries<ScalarType>::Vector;
    using ParameterSet = mio::lsecir::Parameters;

    ScalarType dt = 0.1;

    // Set vector that specifies the number of subcompartments
    std::vector<int> num_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]            = 10;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 10;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = 10;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = 10;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 10;
    mio::lsecir::InfectionState infectionStates(num_subcompartments);

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    int num_transitions = (int)mio::lsecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> init(num_transitions);

    // add time points for initialization of transitions
    Vec init_transitions(num_transitions);
    init_transitions[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 30.0;
    init_transitions[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 30.0;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 30.0 * (1 - 0.09);
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 30.0 * 0.09;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] = 30.0 * (1 - 0.09) * 0.2;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered] =
        30.0 * (1 - 0.09) * (1 - 0.2);
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] = 5.0;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]        = 1.0;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]           = 1.0;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]      = 1.0;
    init_transitions                                                                          = init_transitions * dt;

    // add initial time point to time series
    init.add_time_point(-70, init_transitions);
    // add further time points until time 0
    while (init.get_last_time() < 0) {
        //init_transitions *=  1.01;
        init.add_time_point(init.get_last_time() + dt, init_transitions);
    }

    // Define ParameterSet used for Simulation and initialization.
    ParameterSet parameters;

    // Set Parameters of the model
    parameters.get<mio::lsecir::TimeExposed>()            = 2 * 4.2 - 5.2;
    parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2 * (5.2 - 4.2);
    parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
    parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
    // also possible to change values with setter
    parameters.set<mio::lsecir::TimeInfectedCritical>(7.1);

    parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 9.5));

    parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
    parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
    parameters.set<mio::lsecir::DeathsPerCritical>(0.3);

    mio::lsecir::Initializer initializer(std::move(init), infectionStates, std::move(parameters));
    auto init_compartments = initializer.compute_initializationvector(83000000, 10);
    for (int i = 0; i < infectionStates.get_count(); ++i) {
        std::cout << init_compartments[i] << ", ";
    }
    std::cout << "\n";

    init_compartments[infectionStates.get_firstindex(mio::lsecir::InfectionStateBase::Recovered)] = 0;

    // initialize model
    mio::lsecir::Model model(std::move(init_compartments), infectionStates, std::move(parameters));

    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(0, 10, 0.5, model);
    // calculate the distribution in Infectionstates without subcompartments of the result
    mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);
    // print it
    mio::lsecir::print_TimeSeries(populations, model.get_heading_CompartmentsBase());

    if (save_result) {
        auto save_result_status_subcompartments =
            mio::save_result({result}, {0}, 1, "result_lct_subcompartments_fictional.h5");
        auto save_result_status = mio::save_result({populations}, {0}, 1, "result_lct_fictional.h5");
    }
}