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
#include "load_test_data.h"
#include "ide_seir/model.h"
#include "ide_seir/parameters.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <gtest/gtest.h>

TEST(ModelTestIdeSeirMin, simulateDefault)
{
    int tmax  = 1;
    double dt = 0.1;

    using Vec = mio::TimeSeries<double>::Vector;

    mio::TimeSeries<double> init(1);
    init.add_time_point<Eigen::VectorXd>(-12.0, Vec::Constant(1, 10.));
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, Vec::Constant(1, 10.));
    }

    mio::iseir::Model model(std::move(init), dt, 10);
    model.simulate(tmax);
    auto result = model.calculate_EIR();

    // TODO: It is actually an inconsistency that the simulation of the IDE model goes to tmax + dt,
    // instead of tmax; this has to be corrected in the future.
    EXPECT_NEAR(result.get_last_time(), (double)tmax + 0.1, 1e-10);
}

class ModelTestIdeSeir : public testing::Test
{
protected:
    virtual void SetUp()
    {
        using Vec = mio::TimeSeries<double>::Vector;

        int N     = 810000;
        double dt = 0.1;
        mio::TimeSeries<double> result(1);
        result.add_time_point<Eigen::VectorXd>(-15.0, Vec::Constant(1, N * 0.95));
        while (result.get_last_time() < 0) {
            result.add_time_point(result.get_last_time() + dt,
                                  Vec::Constant(1, (double)result.get_last_value()[0] + result.get_last_time()));
        }

        model = new mio::iseir::Model(std::move(result), dt, N);

        model->parameters.set<mio::iseir::LatencyTime>(3.3);
        model->parameters.set<mio::iseir::InfectiousTime>(8.2);
        model->parameters.set<mio::iseir::TransmissionRisk>(0.015);
        mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
        contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
        contact_matrix[0].add_damping(0.7, mio::SimulationTime(10.));
        model->parameters.get<mio::iseir::ContactFrequency>() = mio::UncertainContactMatrix(contact_matrix);
    }

    virtual void TearDown()
    {
        delete model;
    }

public:
    mio::iseir::Model* model = nullptr;
};

TEST_F(ModelTestIdeSeir, compareWithPreviousRun)
{

    auto compare = load_test_data_csv<double>("ide-seir-compare.csv");
    model->simulate(15);
    auto sim_result = model->calculate_EIR();

    ASSERT_EQ(compare.size(), static_cast<size_t>(sim_result.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(sim_result.get_num_elements()) + 1) << "at row " << i;
        ASSERT_NEAR(sim_result.get_time(i), compare[i][0], 1e-8) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            ASSERT_NEAR(sim_result.get_value(i)[j - 1], compare[i][j], 1e-8) << " at row " << i;
        }
    }
}
