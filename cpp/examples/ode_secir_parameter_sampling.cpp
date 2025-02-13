/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/model.h"

#include <stdio.h>

// verify manually that the correct distribution is chosen
int main()
{
    /*
     * Real valued variable element
     */
    double mean   = 5;
    double stddev = 1.5;
    double min    = 1;
    double max    = 10;
    // check if constructor works correctly
    mio::ParameterDistributionNormal some_parameter(min, max, mean, stddev);
    // "some parameter",
    // std::make_unique<mio::ParameterDistributionNormal>(mio::ParameterDistributionNormal(min, max, mean, stddev))};

    printf("\n N(%.0f,%.0f)-distribution with sampling only in [%.0f,%.0f]", mean, stddev, min, max);
    int counter[10] = {0};
    for (int i = 0; i < 1000; i++) {
        int rounded = (int)(some_parameter.get_sample(mio::thread_local_rng()) - 1);
        if (rounded >= 0 && rounded < 10) {
            counter[rounded]++;
        }
    }
    double acc = 0;
    for (int i = 0; i < 9; i++) {
        acc += (double)counter[i] / 1000.0;
        printf("\n [%d-%d): %.2f %.2f ", i + 1, i + 2, (double)counter[i] / 1000.0, acc);
    }
    printf("\n");

    // check if constructor works correctly
    printf("\n U(%.0f,%.0f)-distribution", min, max);
    mio::ParameterDistributionUniform some_other_parameter(1.0, 10.0);

    double counter_unif[10] = {0};
    for (int i = 0; i < 1000; i++) {
        int rounded = (int)(some_other_parameter.get_sample(mio::thread_local_rng()) - 1);
        if (rounded >= 0 && rounded < 10) {
            counter_unif[rounded]++;
        }
    }
    acc = 0;
    for (int i = 0; i < 9; i++) {
        acc += (double)counter_unif[i] / 1000.0;
        printf("\n [%d-%d): %.2f %.2f ", i + 1, i + 2, (double)counter_unif[i] / 1000.0, acc);
    }

    /*
     * Contact frequency and dampings variable element
     */
    mio::osecir::Model<double> model(3);
    auto& params = model.parameters;

    mio::AgeGroup nb_groups = params.get_num_groups();
    mio::ContactMatrixGroup cm_group{
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.5))};
    params.get<mio::osecir::ContactPatterns<double>>() = cm_group;

    params.get<mio::osecir::ContactPatterns<double>>().get_dampings().push_back(mio::DampingSampling<double>(
        mio::UncertainValue<double>(0.5), mio::DampingLevel(0), mio::DampingType(0), mio::SimulationTime(30.),
        std::vector<size_t>(1, size_t(0)), Eigen::VectorXd::Constant(Eigen::Index(nb_groups.get()), 1.0)));
    params.get<mio::osecir::ContactPatterns<double>>().get_dampings()[0].get_value().set_distribution(
        mio::ParameterDistributionNormal(0.0, 1.0, 0.5, 0.2));

    params.get<mio::osecir::ContactPatterns<double>>().get_dampings().push_back(mio::DampingSampling<double>(
        mio::UncertainValue<double>(0.3), mio::DampingLevel(1), mio::DampingType(0), mio::SimulationTime(10.),
        std::vector<size_t>(1, size_t(0)), Eigen::VectorXd::Constant(Eigen::Index(nb_groups.get()), 1.0)));
    params.get<mio::osecir::ContactPatterns<double>>().get_dampings()[0].get_value().set_distribution(
        mio::ParameterDistributionNormal(0.0, 1.0, 0.4, 0.05));

    draw_sample(model);
    auto& cfmat_sample = params.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat();

    printf("\n\n Number of dampings: %zu\n", cfmat_sample[0].get_dampings().size());

    //Dampings are sorted automatically by time, therefore the second DampingSamping is at the first position
    printf("\n First damping at %.2f with factor %.2f\n", double(cfmat_sample[0].get_dampings()[0].get_time()),
           cfmat_sample[0].get_dampings()[0].get_coeffs()(0, 0));

    // printout the second damping
    printf("\n Damping at day %.2f\n\t", double(cfmat_sample[0].get_dampings()[1].get_time()));
    std::cout << cfmat_sample[0].get_dampings()[1].get_coeffs() << std::endl;
}
