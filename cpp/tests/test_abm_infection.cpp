/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: David Kerkmann
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
#include "test_abm.h"

TEST(TestInfection, init)
{
    auto params = mio::abm::GlobalInfectionParameters{};

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.4)) // Transition to Infected
        .WillOnce(testing::Return(0.6)) // Transition to Recovered_Infected
        .WillRepeatedly(testing::Return(1.0));

    auto infection = mio::abm::Infection(mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, params,
                                         mio::abm::TimePoint(0), true);

    EXPECT_EQ(infection.get_virus_variant(), mio::abm::VirusVariant::Wildtype);
    EXPECT_EQ(infection.is_detected(), true);

    EXPECT_EQ(infection.get_infection_state(mio::abm::TimePoint(3600)), mio::abm::InfectionState::Recovered_Infected);
    EXPECT_NEAR(infection.get_infectivity(mio::abm::TimePoint(72 * 3600)), 0.2689414213699951, 1e-14);
}
