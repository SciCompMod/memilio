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
#include "ide_secir/infection_state.h"
#include <gtest/gtest.h>

TEST(IdeSecir, infection_transitions)
{
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.size(), mio::isecir::InfectionTransitionsCount);

    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(0),
              std::make_pair(mio::isecir::InfectionState::Susceptible, mio::isecir::InfectionState::Exposed));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(1),
              std::make_pair(mio::isecir::InfectionState::Exposed, mio::isecir::InfectionState::InfectedNoSymptoms));
    EXPECT_EQ(
        mio::isecir::InfectionTransitionsMap.at(2),
        std::make_pair(mio::isecir::InfectionState::InfectedNoSymptoms, mio::isecir::InfectionState::InfectedSymptoms));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(3),
              std::make_pair(mio::isecir::InfectionState::InfectedNoSymptoms, mio::isecir::InfectionState::Recovered));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(4), std::make_pair(mio::isecir::InfectionState::InfectedSymptoms,
                                                                         mio::isecir::InfectionState::InfectedSevere));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(5),
              std::make_pair(mio::isecir::InfectionState::InfectedSymptoms, mio::isecir::InfectionState::Recovered));
    EXPECT_EQ(
        mio::isecir::InfectionTransitionsMap.at(6),
        std::make_pair(mio::isecir::InfectionState::InfectedSevere, mio::isecir::InfectionState::InfectedCritical));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(7),
              std::make_pair(mio::isecir::InfectionState::InfectedSevere, mio::isecir::InfectionState::Recovered));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(8),
              std::make_pair(mio::isecir::InfectionState::InfectedCritical, mio::isecir::InfectionState::Dead));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(9),
              std::make_pair(mio::isecir::InfectionState::InfectedCritical, mio::isecir::InfectionState::Recovered));
}