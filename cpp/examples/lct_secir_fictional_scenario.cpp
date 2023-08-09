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
    bool print_result  = true;
    using Vec          = mio::TimeSeries<ScalarType>::Vector;
    using ParameterSet = mio::lsecir::Parameters;

    ScalarType dt_flows         = 0.1;
    ScalarType total_population = 83155031.0;

    // Define ParameterSet used for simulation and initialization.
    ParameterSet parameters;
    parameters.get<mio::lsecir::TimeExposed>()                      = 3.335;
    parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()           = 3.31331;
    parameters.get<mio::lsecir::TimeInfectedSymptoms>()             = 6.94547;
    parameters.get<mio::lsecir::TimeInfectedSevere>()               = 7.28196;
    parameters.get<mio::lsecir::TimeInfectedCritical>()             = 13.066;
    parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.0733271;

    mio::ContactMatrixGroup& contact_matrix = parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2.7463));
    contact_matrix[0].add_damping(0.5, mio::SimulationTime(5.));
    /*contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2 * 2.7463));
    contact_matrix[0].add_damping(0.5, mio::SimulationTime(-10.));
    contact_matrix[0].add_damping(0.5, mio::SimulationTime(4.9));
    contact_matrix[0].add_damping(0., mio::SimulationTime(5.));
    contact_matrix[0].finalize();*/

    parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 1;
    parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.3;
    parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.206901;
    parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.0786429;
    parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.173176;
    parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.217177;

    // The initialization vector for the LCT model is calculated by defining transitions.
    // Create TimeSeries with num_transitions elements.
    int num_transitions = (int)mio::lsecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> init(num_transitions);

    // Add time points for initialization of transitions.
    /* For this example, the intention is to create nearly constant values for SusceptiblesToExposed flow 
    at the beginning of the simulation. Therefore we initalize the flows accordingly constant for 
    SusceptiblesToExposed and derive matching values for the other flows.*/
    // 7-Tage-Inzidenz at 15.10.2020 was 34.1, see https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/Okt_2020/2020-10-15-de.pdf?__blob=publicationFile.
    ScalarType SusceptibleToExposed_const = (34.1 / 7) * total_population / 100000;
    Vec init_transitions(num_transitions);
    init_transitions[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]        = SusceptibleToExposed_const;
    init_transitions[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = SusceptibleToExposed_const;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
        SusceptibleToExposed_const * (1 - parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>());
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
        SusceptibleToExposed_const * parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>();
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        parameters.get<mio::lsecir::SeverePerInfectedSymptoms>();
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        (1 - parameters.get<mio::lsecir::SeverePerInfectedSymptoms>());
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        parameters.get<mio::lsecir::CriticalPerSevere>();
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        (1 - parameters.get<mio::lsecir::CriticalPerSevere>());
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        parameters.get<mio::lsecir::DeathsPerCritical>();
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        (1 - parameters.get<mio::lsecir::DeathsPerCritical>());
    init_transitions = init_transitions * dt_flows;

    // add initial time point to time series
    init.add_time_point(-200, init_transitions);
    // add further time points until time 0
    while (init.get_last_time() < 0) {
        //init_transitions *=  1.01;
        init.add_time_point(init.get_last_time() + dt_flows, init_transitions);
    }

    // Set vector that specifies the number of subcompartments.
    int num_subcompartments = 1;
    std::vector<int> vec_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]            = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = num_subcompartments;
    mio::lsecir::InfectionState infectionStates(vec_subcompartments);

    mio::TimeSeries<ScalarType> init_copy(init);

    // Get initialization vector for LCT model with num_subcompartments subcompartments.
    mio::lsecir::Initializer initializer(std::move(init_copy), infectionStates, std::move(parameters));
    auto init_compartments = initializer.compute_initializationvector(total_population, 10);

    // Initialize model and perform simulation.
    mio::lsecir::Model model(std::move(init_compartments), infectionStates, std::move(parameters));
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(0, 20, 0.5, model);
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);

    if (print_result) {
        mio::lsecir::print_TimeSeries(populations, model.get_heading_CompartmentsBase());
    }
    if (save_result) {
        auto save_result_status_subcompartments =
            mio::save_result({result}, {0}, 1, "result_lct_subcompartments_fictional_1.h5");
        auto save_result_status = mio::save_result({populations}, {0}, 1, "result_lct_fictional_1.h5");
    }

    // Perform an other simulation with a different number of subcompartments.
    // Set vector that specifies the number of subcompartments.
    int num_subcompartments2 = 10;
    std::vector<int> vec_subcompartments2((int)mio::lsecir::InfectionStateBase::Count, 1);
    vec_subcompartments2[(int)mio::lsecir::InfectionStateBase::Exposed]            = num_subcompartments2;
    vec_subcompartments2[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = num_subcompartments2;
    vec_subcompartments2[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = num_subcompartments2;
    vec_subcompartments2[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = num_subcompartments2;
    vec_subcompartments2[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = num_subcompartments2;
    mio::lsecir::InfectionState infectionStates2(vec_subcompartments2);

    // Get initialization vector for LCT model with num_subcompartments subcompartments.
    mio::lsecir::Initializer initializer2(std::move(init), infectionStates2, std::move(parameters));
    auto init_compartments2 = initializer2.compute_initializationvector(total_population, 10);

    // Initialize model and perform simulation.
    mio::lsecir::Model model2(std::move(init_compartments2), infectionStates2, std::move(parameters));
    mio::TimeSeries<ScalarType> result2 = mio::lsecir::simulate(0, 20, 0.5, model2);
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations2 = model2.calculate_populations(result2);

    if (print_result) {
        mio::lsecir::print_TimeSeries(populations2, model2.get_heading_CompartmentsBase());
    }
    if (save_result) {
        auto save_result_status_subcompartments2 =
            mio::save_result({result2}, {0}, 1, "result_lct_subcompartments_fictional_10.h5");
        auto save_result_status2 = mio::save_result({populations2}, {0}, 1, "result_lct_fictional_10.h5");
    }
}