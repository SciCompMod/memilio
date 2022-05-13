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

#include "ide_seir/ide_seir.h"
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
    *   The TimeSeries should satisfy the conditions of the IDEmodel. Accordingly, the first time point is set to -15.
    */
    result.add_time_point<Eigen::VectorXd>(-15.0, Vec::Constant(1, N * 0.95));
    while (result.get_last_time() < 0) {
        result.add_time_point(result.get_last_time() + dt,
                              Vec::Constant(1, (double)result.get_last_value()[0] + result.get_last_time()));
    }

    // initialize model
    mio::iseir::IdeModel model(std::move(result), dt, N);

    model.m_parameters.set<mio::iseir::LatencyTime>(3.3);
    model.m_parameters.set<mio::iseir::InfectiousTime>(8.2);
    model.m_parameters.set<mio::iseir::TransmissionRisk>(0.015);
    model.m_parameters.get<mio::iseir::ContactFrequency>() = mio::iseir::ContactFrequency::get_default();

    //carry out simulation
    model.simulate(tmax);
    // calculate values for compartments EIR as well
    model.calculate_EIR();
    model.print_result(true);
}