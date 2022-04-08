/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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

#include "seir_ide/seir_ide.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"

#include <vector>
#include <iostream>

int main()
{
    /*  Note that the initial values as well as all other parameters are randomly chosen for this example and are not intend to depict the real world.
    *   This example has only the purpose to show how the IDE model can be applied. 
    */

    using Vec = mio::TimeSeries<double>::Vector;

    int tmax  = 15;
    int N     = 810000;
    double dt = 0.1;
    mio::TimeSeries<double> result(1);

    /*  construct initial TimeSeries with initial times and related quantity of Susceptibles. 
    *   The first time point should satisfy the condition of the IDE model 
    *   (should be earlier than time -(k-1)*dt with k=std::ceil((infectiousTime+latencyTime)/dt))
     */
    result.add_time_point<Eigen::VectorXd>(-16.5, Vec::Constant(1, N * 0.95));
    while (result.get_last_time() < 0) {
        result.add_time_point(result.get_last_time() + dt,
                              Vec::Constant(1, (double)result.get_last_value()[0] + result.get_last_time() / 10.0));
    }

    // initialize model
    mio::IdeModel model(std::move(result), dt, N);

    // set contact matrix as well as dampings; Note: use effective contacts (quantity of Contacts * probability of infection in case of contact)
    // values randomly chosen here as well (here such that initial reproduction number equals 1)
    mio::ContactMatrix& contact_matrix = model.get_contact_matrix();
    contact_matrix =
        mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 1 / 8.2), Eigen::MatrixXd::Constant(1, 1, 0.5 / 8.2));
    contact_matrix.add_damping(0.7, mio::SimulationTime(3.0));

    // carry out simulation
    model.simulate(tmax);
    // calculate values for compartments EIR as well
    model.calculate_EIR();
    model.print_result(true);
}