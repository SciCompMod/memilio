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
    using Vec = mio::TimeSeries<double>::Vector;

    int tmax  = 10;
    int N     = 100;
    double dt = 1;
    mio::TimeSeries<double> result(1);

    result.add_time_point<Eigen::VectorXd>(-5.0, Vec::Constant(1, N));
    while (result.get_last_time() < 0) {
        result.add_time_point(result.get_last_time() + dt, Vec::Constant(1, (double)result.get_last_value()[0] - 1));
    }

    std::cout << "# time  |  number of susceptibles" << std::endl;
    Eigen::Index num_points = result.get_num_time_points();
    for (Eigen::Index i = 0; i < num_points; ++i) {
        std::cout << result.get_time(i) << "      |  " << result[i]
                  << std::endl; //[Eigen::Index(InfectionState::S)] << std::endl;
    }

    // Initialize model.
    mio::isecir::Model model(std::move(result), dt, N);

    // Set working parameters.
    // model.parameters.set<mio::iseir::LatencyTime>(3.3);
    // model.parameters.set < mio::isecir::mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    // contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
    // //model.parameters.get<mio::iseir::ContactFrequency>() = mio::UncertainContactMatrix(contact_matrix);

    // // Carry out simulation.
    model.simulate(tmax);
}