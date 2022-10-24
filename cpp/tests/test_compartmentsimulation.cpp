/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Jan Kleinert, Daniel Abele
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

#include "memilio/compartments/simulation.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

TEST(TestCompartmentSimulation, integrator_uses_model_reference)
{
    struct MockModel {
        Eigen::VectorXd get_initial_values() const
        {
            return Eigen::VectorXd::Zero(1);
        }
        void eval_right_hand_side(const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&,
                                  double, Eigen::Ref<Eigen::VectorXd> dydt) const
        {
            dydt[0] = this->m_dydt;
        }
        double m_dydt = 1.0;
    };

    auto sim = mio::Simulation<MockModel>(MockModel(), 0.0);
    sim.advance(1.0);

    ASSERT_NEAR(sim.get_result().get_last_value()[0], 1.0, 1e-5);

    //modifying the model from the outside should affect the integration result
    sim.get_model().m_dydt = 2.0;
    sim.advance(2.0);

    ASSERT_NEAR(sim.get_result().get_last_value()[0], 3.0, 1e-5);
}
