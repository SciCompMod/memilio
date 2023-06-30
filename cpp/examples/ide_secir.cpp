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

#include "ide_secir/model.h" // IWYU pragma: keep
#include "ide_secir/infection_state.h"
#include "ide_secir/simulation.h" // IWYU pragma: keep
#include "memilio/config.h"
#include "memilio/math/eigen.h" // IWYU pragma: keep
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h" // IWYU pragma: keep
#include <iostream>

int main()
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax        = 10;
    ScalarType N           = 10000;
    ScalarType Dead_before = 12;
    ScalarType dt          = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    mio::TimeSeries<ScalarType> init(num_transitions);

    // add time points for initialization of transitions
    Vec vec_init(num_transitions);
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 25.0;
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 1.0;
    // add initial time point to time series
    init.add_time_point(-10, vec_init);
    // add further time points until time 0
    while (init.get_last_time() < 0) {
        vec_init *= 1.01;
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model<double> model(std::move(init), N, Dead_before);

    // model.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Susceptible] = 1000;
    // model.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered]   = 0;

    // Set working parameters
    // Set max_support for all Delay Distributions
    std::vector<ScalarType> vec_max_support((int)mio::isecir::InfectionTransition::Count, 2);
    vec_max_support[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]                 = 3;
    vec_max_support[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] = 4;
    std::vector<mio::isecir::DelayDistribution> vec_delaydistrib(num_transitions, mio::isecir::DelayDistribution());
    for (int i = 0; i < (int)mio::isecir::InfectionTransition::Count; i++) {
        vec_delaydistrib[i].set_max_support(vec_max_support[i]);
    }
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
    model.parameters.get<mio::isecir::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrix);

    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(0.5);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(0.5);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(0.5);

    model.check_constraints(dt);

    // Carry out simulation.
    mio::isecir::Simulation<double> sim(model, 0, dt);
    sim.advance(tmax);

    sim.print_transitions();

    sim.print_compartments();
}
