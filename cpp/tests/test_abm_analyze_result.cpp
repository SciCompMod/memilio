/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Khoa Nguyen
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

TEST(TestAnalyzeResult, ensembleParamsPercentile)
{
    mio::abm::World world1(6);
    mio::abm::World world2(6);

    auto& params1 = world1.parameters;
    params1.get<mio::abm::InfectedToSevere>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}] = 0.1;
    params1.get<mio::abm::SevereToCritical>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}] = 0.2;

    auto& params2 = world2.parameters;
    params2.get<mio::abm::InfectedToSevere>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}] = 0.2;
    params2.get<mio::abm::SevereToCritical>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}] = 0.3;

    auto g1 = std::vector<mio::abm::World>({world1, world2});

    params1.get<mio::abm::InfectedToRecovered>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}] = 0.2;
    params1.get<mio::abm::InfectedToSevere>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]    = 0.3;
    params1.get<mio::abm::SevereToCritical>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]    = 0.4;

    params2.get<mio::abm::InfectedToRecovered>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}] = 0.7;
    params2.get<mio::abm::InfectedToSevere>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]    = 0.4;
    params2.get<mio::abm::SevereToCritical>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]    = 0.5;

    auto g2 = std::vector<mio::abm::World>({world1, world2});

    auto ensemble_params = std::vector<std::vector<mio::abm::World>>({g1, g2});

    auto ensemble_p49_params = mio::abm::ensemble_params_percentile(ensemble_params, 0.49);
    auto ensemble_p51_params = mio::abm::ensemble_params_percentile(ensemble_params, 0.51);

    auto check1 =
        ensemble_p49_params[0]
            .parameters.get<mio::abm::InfectedToSevere>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]
            .value();
    auto check2 =
        ensemble_p49_params[1]
            .parameters.get<mio::abm::InfectedToSevere>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]
            .value();

    EXPECT_EQ(check1, 0.1);
    EXPECT_EQ(check2, 0.2);

    auto check3 =
        ensemble_p51_params[0]
            .parameters.get<mio::abm::InfectedToSevere>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]
            .value();
    auto check4 =
        ensemble_p51_params[1]
            .parameters.get<mio::abm::InfectedToSevere>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]
            .value();

    EXPECT_EQ(check3, 0.3);
    EXPECT_EQ(check4, 0.4);

    auto check5 =
        ensemble_p49_params[0]
            .parameters.get<mio::abm::SevereToCritical>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]
            .value();
    auto check6 =
        ensemble_p49_params[1]
            .parameters.get<mio::abm::SevereToCritical>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]
            .value();

    EXPECT_EQ(check5, 0.2);
    EXPECT_EQ(check6, 0.3);

    auto check7 =
        ensemble_p51_params[0]
            .parameters.get<mio::abm::SevereToCritical>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]
            .value();
    auto check8 =
        ensemble_p51_params[1]
            .parameters.get<mio::abm::SevereToCritical>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Unvaccinated}]
            .value();

    EXPECT_EQ(check7, 0.4);
    EXPECT_EQ(check8, 0.5);
}
