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

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>

int main()
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax  = 10;
    ScalarType N     = 10000;
    ScalarType Dead0 = 12;
    ScalarType dt    = 1;

    int num_transitions = (int)mio::isecir::InfectionTransitions::Count;

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    mio::TimeSeries<ScalarType> init(num_transitions);

    // add time points for initialization
    Vec vec_init(num_transitions);
    vec_init << 25.0, 15.0, 8.0, 4.0, 1.0, 4.0, 1.0, 1.0, 1.0, 1.0;
    init.add_time_point(-10, vec_init);
    while (init.get_last_time() < 0) {
        vec_init *= 1.01;
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), dt, N, Dead0);

    // // Set working parameters.
    model.parameters.set<mio::isecir::TransitionDistributions>(
        std::vector<mio::isecir::DelayDistribution>(num_transitions, mio::isecir::DelayDistribution()));
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransitions::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransitions::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransitions::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);
    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(0.5);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(0.5);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(0.5);

    // Carry out simulation.
    mio::TimeSeries<ScalarType> compartments = model.simulate(tmax);

    model.print_transitions();
    std::cout << "\n" << std::endl;
    model.print_compartments();
}