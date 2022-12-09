/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn, Khoa Nguyen
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
#include "abm/location_type.h"
#include "abm/mask.h"
#include "abm/mask_type.h"
#include "abm/person.h"
#include "abm/state.h"
#include "test_abm.h"
#include "gmock/gmock.h"
#include <gtest/gtest.h>

TEST(TestMasks, init)
{
    auto mask = mio::abm::Mask(mio::abm::MaskType::Count);
    ASSERT_EQ(mask.get_time_used().seconds(), 0.);
}

TEST(TestMasks, maskUsage)
{
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto target = mio::abm::Location(mio::abm::LocationType::Work, 0);
    target.set_npi_active(true);
    target.set_required_mask(mio::abm::MaskType::Surgical);

    auto person = mio::abm::Person(home, mio::abm::InfectionState::Count, mio::abm::AgeGroup::Count, {});
    person.get_mask().change_mask(mio::abm::MaskType::Community);

    person.mask_usage(target);

    ASSERT_EQ(person.get_mask().get_type(), mio::abm::MaskType::Surgical);
    ASSERT_TRUE(person.get_wear_mask());

    target.set_npi_active(false);
    person.mask_usage(target);

    ASSERT_FALSE(person.get_wear_mask());

    auto preferences = mio::CustomIndexArray<double, mio::abm::LocationType>({mio::abm::LocationType::Count}, -1.);
    person.set_mask_preferences(preferences);
    person.get_mask().change_mask(mio::abm::MaskType::Community);
    person.mask_usage(target);

    ASSERT_EQ(person.get_mask().get_type(), mio::abm::MaskType::Community);
    ASSERT_FALSE(person.get_wear_mask());

    preferences = mio::CustomIndexArray<double, mio::abm::LocationType>({mio::abm::LocationType::Count}, 1.);
    person.set_mask_preferences(preferences);
    target.set_npi_active(false);
    person.mask_usage(target);

    ASSERT_TRUE(person.get_wear_mask());
}

TEST(TestMasks, MaskProtection)
{
    mio::abm::AgeGroup age =
        mio::abm::AgeGroup(mio::UniformIntDistribution<int>()(0, int(mio::abm::AgeGroup::Count) - 1));
    mio::abm::VaccinationState vaccination_state =
        mio::abm::VaccinationState(mio::UniformIntDistribution<int>()(0, int(mio::abm::VaccinationState::Count) - 1));

    mio::abm::GlobalInfectionParameters params;
    params.set<mio::abm::CarrierToInfected>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::CarrierToInfected>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::CarrierToRecovered>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::CarrierToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::DetectInfection>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::DetectInfection>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::InfectedToSevere>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::InfectedToSevere>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::InfectedToRecovered>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::InfectedToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::SevereToCritical>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::SevereToCritical>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::SevereToRecovered>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::SevereToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::CriticalToDead>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::CriticalToDead>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::CriticalToRecovered>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::CriticalToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::RecoveredToSusceptible>({{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::RecoveredToSusceptible>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::SusceptibleToExposedByCarrier>(
        {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::SusceptibleToExposedByCarrier>()[{age, vaccination_state}] = 1.;
    params.set<mio::abm::SusceptibleToExposedByInfected>(
        {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::SusceptibleToExposedByInfected>()[{age, vaccination_state}] = 1.;

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(2))
        .WillOnce(testing::Return(std::numeric_limits<int>::max()))
        .WillOnce(testing::Return(std::numeric_limits<int>::max()));

    //setup location with some chance of exposure
    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto location = mio::abm::Location(mio::abm::LocationType::School, 0);
    auto person1  = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, params,
                                    vaccination_state);
    home.add_person(person1);
    auto person2 = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, params,
                                    vaccination_state);
    home.add_person(person2);
    auto infected1 = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, params,
                                      vaccination_state);
    home.add_person(infected1);
    infected1.migrate_to(home, location);

    //cache precomputed results
    auto dt = mio::abm::seconds(8640);
    location.begin_step(dt, params);
    // person1 wears a mask
    person1.set_wear_mask(true);
    person1.migrate_to(home, location);
    // person2 does not wear a mask
    person2.set_wear_mask(false);
    person2.migrate_to(home, location);
    location.interact(person1, dt, params);
    location.interact(person2, dt, params);

    // The person should have full protection against an infection
    ASSERT_EQ(person1.get_infection_state(), mio::abm::InfectionState::Susceptible);
    ASSERT_EQ(person2.get_infection_state(), mio::abm::InfectionState::Exposed);
}