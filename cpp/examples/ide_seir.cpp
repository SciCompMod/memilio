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

#include "ide_seir/model.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"

int main()
{
    /**
    * Note: the initial values as well as all other parameters are randomly chosen for this example and are not 
    * intended to characterize the real world.
    * This example has the purpose to show how the IDE SEIR model can be applied. 
    */

    using Vec = mio::TimeSeries<double>::Vector;

    int tmax  = 15;
    int N     = 810000;
    double dt = 0.1;
    mio::TimeSeries<double> result(1);

    /**
    * Construction of the initial TimeSeries with point of times and the corresponding number of susceptibles.  
    * The smallest time should be small enough. See the documentation of the IdeSeirModel constructor for 
    * detailed information. Initial data are chosen randomly.
    */
    result.add_time_point<Eigen::VectorXd>(-15.0, Vec::Constant(1, N * 0.95));
    while (result.get_last_time() < 0) {
        result.add_time_point(result.get_last_time() + dt,
                              Vec::Constant(1, (double)result.get_last_value()[0] + result.get_last_time()));
    }

    // Initialize model.
    mio::iseir::Model model(std::move(result), dt, N);

    // Set working parameters.
    model.parameters.set<mio::iseir::LatencyTime>(3.3);
    model.parameters.set<mio::iseir::InfectiousTime>(8.2);
    model.parameters.set<mio::iseir::TransmissionRisk>(0.015);
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
    // Add damping.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(10.));
    model.parameters.get<mio::iseir::ContactFrequency>() = mio::UncertainContactMatrix(contact_matrix);

    // Carry out simulation.
    model.simulate(tmax);
    // Calculate values for compartments EIR.
    model.calculate_EIR();
    //Print results.
    model.print_result(true);
}
