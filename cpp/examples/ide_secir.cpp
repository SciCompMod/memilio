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
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>

int main()
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    int tmax     = 10;
    size_t N     = 1000;
    size_t Dead0 = 12;
    double dt    = 1;
    // for SECIHURD model we need 6 transitions for simulation
    size_t num_transitions = 6;

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    mio::TimeSeries<ScalarType> init(num_transitions);
    std::cout << "time points: " << init.get_num_time_points()
              << ", elements: " << init.get_num_elements() << "\n";

    // add time points for initialization
    init.add_time_point<Eigen::VectorXd>(-10, Vec::Constant(num_transitions, 40));
    std::cout << "time points: " << init.get_num_time_points()
              << ", elements: " << init.get_num_elements() << "\n";
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt,
                                        Vec::Constant(num_transitions, init.get_last_value()[0]+1.0));
    }

    // print initial transitions
    std::cout << "# time  |  S -> E  |  E - > C  |  C -> I  |  I -> H  |  H -> U  |  U -> D" << std::endl;
    Eigen::Index num_points = init.get_num_time_points();
    for (Eigen::Index i = 0; i < num_points; ++i) {
        std::cout << init.get_time(i) << "      |  " << init[i][0] << "  |  "
                  << init[i][1] << "  |  " << init[i][2] << "  |  " << init[i][3]
                  << "  |  " << init[i][4] << "  |  " << init[i][5]
                  << std::endl; //[Eigen::Index(InfectionState::S)] << std::endl;
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), dt, N, Dead0);

    // Set working parameters.
    model.parameters.set<mio::isecir::TransitionDistributions>(
        std::vector<mio::isecir::DelayDistribution>(9, mio::isecir::DelayDistribution()));
    model.parameters.set<mio::isecir::TransitionProbabilities>(std::vector<double>(6, 0.5));
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>()= mio::UncertainContactMatrix(contact_matrix);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(1.0);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(1.0);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(1.0);

    // // Carry out simulation.
    model.simulate(tmax);

    std::cout << "Ended simulation"
              << "\n";

    model.print_transitions();
    model.print_compartments();
}