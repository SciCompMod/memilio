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
#include "test_abm.h"

TEST(TestLocation, init)
{
    auto location = mio::abm::Location(mio::abm::LocationType::School, 0);
    for (mio::abm::InfectionState i = mio::abm::InfectionState(0); i < mio::abm::InfectionState::Count;
         i                          = mio::abm::InfectionState(size_t(i) + 1)) {
        ASSERT_EQ(location.get_subpopulation(i), 0);
    }
    ASSERT_EQ(print_wrap(location.get_population().get_last_value()),
              print_wrap(mio::TimeSeries<double>::Vector::Zero((size_t)mio::abm::InfectionState::Count)));
}

TEST(TestLocation, initCell)
{
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 2);
    ASSERT_EQ(location.get_cells().size(), 2);
}

TEST(TestLocation, getIndex)
{
    auto location = mio::abm::Location(mio::abm::LocationType::Home, 0);
    ASSERT_EQ((int)location.get_index(), 0);
}

TEST(TestLocation, addRemovePerson)
{
    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0, 0);
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 3);
    auto person1  = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14, {});
    home.add_person(person1);
    person1.migrate_to(home, location, {0, 1});
    auto person2 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age15to34, {});
    home.add_person(person2);
    person2.migrate_to(home, location, {0});
    auto person3 = mio::abm::Person(home, mio::abm::InfectionState::Exposed, mio::abm::AgeGroup::Age35to59, {});
    home.add_person(person3);
    person3.migrate_to(home, location, {0, 1});

    ASSERT_EQ(location.get_subpopulation(mio::abm::InfectionState::Infected), 2);
    ASSERT_EQ(location.get_subpopulation(mio::abm::InfectionState::Exposed), 1);
    ASSERT_EQ(location.get_cells()[0].num_people, 3u);
    ASSERT_EQ(location.get_cells()[1].num_people, 2u);
    ASSERT_EQ(location.get_cells()[2].num_people, 0u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 2u);
    ASSERT_EQ(location.get_cells()[1].num_infected, 1u);
    ASSERT_EQ(location.get_cells()[2].num_infected, 0u);

    location.remove_person(person2);

    ASSERT_EQ(location.get_subpopulation(mio::abm::InfectionState::Infected), 1);
    ASSERT_EQ(location.get_subpopulation(mio::abm::InfectionState::Exposed), 1);
    ASSERT_EQ(location.get_cells()[0].num_people, 2u);
    ASSERT_EQ(location.get_cells()[1].num_people, 2u);
    ASSERT_EQ(location.get_cells()[2].num_people, 0u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 1u);
    ASSERT_EQ(location.get_cells()[1].num_infected, 1u);
    ASSERT_EQ(location.get_cells()[2].num_infected, 0u);
}

TEST(TestLocation, beginStep)
{
    using testing::Return;

    {
        // Test should work identically work with any age.
        mio::abm::AgeGroup age =
            mio::abm::AgeGroup(mio::UniformIntDistribution<int>()(0, int(mio::abm::AgeGroup::Count) - 1));
        mio::abm::VaccinationState vaccination_state = mio::abm::VaccinationState(
            mio::UniformIntDistribution<int>()(0, int(mio::abm::VaccinationState::Count) - 1));

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
        params.set<mio::abm::RecoveredToSusceptible>(
            {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
        params.get<mio::abm::RecoveredToSusceptible>()[{age, vaccination_state}] = 0.5;
        params.set<mio::abm::SusceptibleToExposedByCarrier>(
            {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
        params.get<mio::abm::SusceptibleToExposedByCarrier>()[{age, vaccination_state}] = 0.4;
        params.set<mio::abm::SusceptibleToExposedByInfected>(
            {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
        params.get<mio::abm::SusceptibleToExposedByInfected>()[{age, vaccination_state}] = 0.5;

        //setup location with some chance of exposure
        auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0, 0);
        auto location1 = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 3);
        auto infected1 = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34,
                                          params, vaccination_state);
        home.add_person(infected1);
        infected1.migrate_to(home, location1, {0});
        auto infected2 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age80plus,
                                          params, vaccination_state);
        home.add_person(infected2);
        infected2.migrate_to(home, location1, {0, 1});
        auto infected3 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14,
                                          params, vaccination_state);
        home.add_person(infected3);
        infected3.migrate_to(home, location1, {1});

        //cache precomputed results
        auto dt = mio::abm::seconds(8640);
        location1.begin_step(dt, params);

        EXPECT_NEAR((location1.get_cells()[0].cached_exposure_rate[{age, vaccination_state}]), 0.9, 1e-14);
        EXPECT_NEAR((location1.get_cells()[1].cached_exposure_rate[{age, vaccination_state}]), 1, 1e-14);
        EXPECT_NEAR((location1.get_cells()[2].cached_exposure_rate[{age, vaccination_state}]), 0, 1e-14);
    }

    {
        mio::abm::AgeGroup age =
            mio::abm::AgeGroup(mio::UniformIntDistribution<int>()(0, int(mio::abm::AgeGroup::Count) - 1));
        mio::abm::VaccinationState vaccination_state = mio::abm::VaccinationState(
            mio::UniformIntDistribution<int>()(0, int(mio::abm::VaccinationState::Count) - 1));

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
        params.set<mio::abm::RecoveredToSusceptible>(
            {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
        params.get<mio::abm::RecoveredToSusceptible>()[{age, vaccination_state}] = 0.5;
        params.set<mio::abm::SusceptibleToExposedByCarrier>(
            {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
        params.get<mio::abm::SusceptibleToExposedByCarrier>()[{age, vaccination_state}] = 0.4;
        params.set<mio::abm::SusceptibleToExposedByInfected>(
            {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
        params.get<mio::abm::SusceptibleToExposedByInfected>()[{age, vaccination_state}] = 0.5;

        //setup location with some chance of exposure
        auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
        auto location1 = mio::abm::Location(mio::abm::LocationType::School, 0);
        location1.set_capacity(3, 18);
        auto infected1 = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34,
                                          params, vaccination_state);
        home.add_person(infected1);
        infected1.migrate_to(home, location1);
        auto infected2 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age80plus,
                                          params, vaccination_state);
        home.add_person(infected2);
        infected2.migrate_to(home, location1);
        auto infected3 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14,
                                          params, vaccination_state);
        home.add_person(infected3);
        infected3.migrate_to(home, location1);

        //cache precomputed results
        auto dt = mio::abm::seconds(8640);
        location1.set_capacity_adapted_transmission_risk_flag(true);
        location1.begin_step(dt, params);

        EXPECT_NEAR((location1.get_cached_exposure_rate()[{age, vaccination_state}]), 15.4, 1e-14);
    }
}

TEST(TestLocation, reachCapacity)
{
    using testing::Return;

    auto world     = mio::abm::World();
    auto home_id   = world.add_location(mio::abm::LocationType::Home);
    auto school_id = world.add_location(mio::abm::LocationType::School);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)); // draw random school hour
    // .WillRepeatedly(testing::Return(1.0));

    auto& p1 = world.add_person(home_id, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age5to14);
    auto& p2 = world.add_person(home_id, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14);

    p1.set_assigned_location(school_id);
    p2.set_assigned_location(school_id);
    p1.set_assigned_location(home_id);
    p2.set_assigned_location(home_id);

    auto& home   = world.get_individualized_location(home_id);
    auto& school = world.get_individualized_location(school_id);
    school.set_capacity(1, 66);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

    world.evolve(mio::abm::TimePoint(0) + mio::abm::hours(8), mio::abm::hours(1));

    ASSERT_EQ(p1.get_location_id().type, mio::abm::LocationType::School);
    ASSERT_EQ(p2.get_location_id().type, mio::abm::LocationType::Home); // p2 should not be able to enter the school
    ASSERT_EQ(school.get_population().get_last_value().sum(), 1);
    ASSERT_EQ(home.get_population().get_last_value().sum(), 1);
}

TEST(TestLocation, computeRelativeTransmissionRisk)
{
    using testing::Return;

    mio::abm::AgeGroup age =
        mio::abm::AgeGroup(mio::UniformIntDistribution<int>()(0, int(mio::abm::AgeGroup::Count) - 1));
    mio::abm::VaccinationState vaccination_state =
        mio::abm::VaccinationState(mio::UniformIntDistribution<int>()(0, int(mio::abm::VaccinationState::Count) - 1));

    mio::abm::GlobalInfectionParameters params;

    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    home.set_capacity(4, 264);
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0);
    location.set_capacity(4, 264);

    auto infected1 = mio::abm::Person(home, mio::abm::InfectionState::Carrier, age, params, vaccination_state);
    home.add_person(infected1);
    auto infected2 = mio::abm::Person(home, mio::abm::InfectionState::Carrier, age, params, vaccination_state);
    location.add_person(infected2);

    location.set_capacity_adapted_transmission_risk_flag(true);

    home.compute_relative_transmission_risk();
    location.compute_relative_transmission_risk();

    ASSERT_EQ(location.compute_relative_transmission_risk(), 0.25);
    ASSERT_EQ(home.compute_relative_transmission_risk(), 1.0);
}

TEST(TestLocation, changedState)
{
    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0, 0);
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 1);
    auto p1       = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, {});
    home.add_person(p1);
    p1.migrate_to(home, location, {0});
    auto p2 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age80plus, {});
    home.add_person(p2);
    p2.migrate_to(home, location, {0});
    auto p3 = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age80plus, {});
    home.add_person(p3);
    p3.migrate_to(home, location, {0});

    ASSERT_EQ(location.get_cells()[0].num_carriers, 1u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 1u);
    location.changed_state(p1, mio::abm::InfectionState::Susceptible);
    ASSERT_EQ(location.get_cells()[0].num_carriers, 2u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 1u);
    location.changed_state(p2, mio::abm::InfectionState::Carrier);
    ASSERT_EQ(location.get_cells()[0].num_carriers, 1u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 2u);
    location.changed_state(p3, mio::abm::InfectionState::Infected);
    ASSERT_EQ(location.get_cells()[0].num_carriers, 1u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 1u);
}

TEST(TestLocation, interact)
{
    using testing::Return;

    // Test should work identically work with any age.
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
    params.get<mio::abm::SusceptibleToExposedByCarrier>()[{age, vaccination_state}] = 0.5;
    params.set<mio::abm::SusceptibleToExposedByInfected>(
        {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::SusceptibleToExposedByInfected>()[{age, vaccination_state}] = 0.5;

    //setup location with some chance of exposure
    auto location  = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto infected1 = mio::abm::Person(location, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34,
                                      params, vaccination_state);
    location.add_person(infected1);
    auto infected2 = mio::abm::Person(location, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age80plus,
                                      params, vaccination_state);
    location.add_person(infected2);
    auto infected3 = mio::abm::Person(location, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14,
                                      params, vaccination_state);
    location.add_person(infected3);

    //cache precomputed results
    auto dt = mio::abm::seconds(8640); //0.1 days
    location.begin_step(dt, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    {
        auto susceptible =
            mio::abm::Person(location, mio::abm::InfectionState::Susceptible, age, params, vaccination_state);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(susceptible, dt, params), mio::abm::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.15));
        EXPECT_EQ(location.interact(susceptible, dt, params), mio::abm::InfectionState::Susceptible);
    }

    {
        auto exposed = mio::abm::Person(location, mio::abm::InfectionState::Exposed, age, params, vaccination_state);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(0); //no transitions out of exposed state
        EXPECT_EQ(location.interact(exposed, dt, params), mio::abm::InfectionState::Exposed);
    }

    {
        auto carrier = mio::abm::Person(location, mio::abm::InfectionState::Carrier, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(carrier, dt, params), mio::abm::InfectionState::Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.099));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(carrier, dt, params), mio::abm::InfectionState::Recovered_Carrier);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(carrier, dt, params), mio::abm::InfectionState::Carrier);
    }

    {
        auto infected = mio::abm::Person(location, mio::abm::InfectionState::Infected, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(infected, dt, params), mio::abm::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(infected, dt, params), mio::abm::InfectionState::Infected_Severe);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(infected, dt, params), mio::abm::InfectionState::Infected);
    }

    {
        auto severe =
            mio::abm::Person(location, mio::abm::InfectionState::Infected_Severe, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(severe, dt, params), mio::abm::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(severe, dt, params), mio::abm::InfectionState::Infected_Critical);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(severe, dt, params), mio::abm::InfectionState::Infected_Severe);
    }

    {
        auto critical =
            mio::abm::Person(location, mio::abm::InfectionState::Infected_Critical, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(critical, dt, params), mio::abm::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(critical, dt, params), mio::abm::InfectionState::Dead);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(critical, dt, params), mio::abm::InfectionState::Infected_Critical);
    }

    for (auto&& recovered_state :
         {mio::abm::InfectionState::Recovered_Carrier, mio::abm::InfectionState::Recovered_Infected}) {
        auto recovered = mio::abm::Person(location, recovered_state, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(recovered, dt, params), mio::abm::InfectionState::Susceptible);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(recovered, dt, params), recovered_state);
    }

    //setup location with 2 cells with some chance of exposure
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0, 0);
    auto location2 = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 2);
    auto infected4 = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, params,
                                      vaccination_state);
    home.add_person(infected4);
    infected4.migrate_to(home, location2, {0, 1});
    auto infected5 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age80plus, params,
                                      vaccination_state);
    home.add_person(infected5);
    infected5.migrate_to(home, location2, {0});
    auto infected6 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14, params,
                                      vaccination_state);
    home.add_person(infected6);
    infected6.migrate_to(home, location2, {1});

    //cache precomputed results
    location2.begin_step(dt, params);

    {
        auto susceptible =
            mio::abm::Person(home, mio::abm::InfectionState::Susceptible, age, params, vaccination_state);
        home.add_person(susceptible);
        susceptible.migrate_to(home, location2, {0, 1});
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location2.interact(susceptible, dt, params), mio::abm::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(2).WillOnce(Return(0.2)).WillOnce(Return(0.07));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location2.interact(susceptible, dt, params), mio::abm::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(2).WillOnce(Return(0.15)).WillOnce(Return(0.15));
        EXPECT_EQ(location2.interact(susceptible, dt, params), mio::abm::InfectionState::Susceptible);
    }
}

TEST(TestLocation, setCapacity)
{
    auto location = mio::abm::Location(mio::abm::LocationType::Home, 0);
    location.set_capacity(4, 200);
    ASSERT_EQ(location.get_capacity().persons, 4);
    ASSERT_EQ(location.get_capacity().volume, 200);
}

TEST(TestLocation, addSubpopulationsTimepoint)
{
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 3);
    auto person1  = mio::abm::Person(location, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14, {});
    location.add_person(person1);
    auto person2 = mio::abm::Person(location, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age15to34, {});
    location.add_person(person2);
    auto person3 = mio::abm::Person(location, mio::abm::InfectionState::Exposed, mio::abm::AgeGroup::Age35to59, {});
    location.add_person(person3);

    auto t1 = mio::abm::TimePoint(0) + mio::abm::hours(7);
    location.add_subpopulations_timepoint(t1);
    auto v1 = location.get_population().get_value(1);
    // Check whether the number of persons in infected state at the location is correct
    ASSERT_EQ(v1[size_t(mio::abm::InfectionState::Infected)], 2);

    auto t2 = mio::abm::TimePoint(0) + mio::abm::hours(14);
    person1.set_infection_state(mio::abm::InfectionState::Infected_Critical);
    location.changed_state(person1, mio::abm::InfectionState::Infected);
    location.add_subpopulations_timepoint(t2);
    auto v2 = location.get_population().get_value(2);
    // Check whether the number of persons in infected state at the location is correct
    ASSERT_EQ(v2[size_t(mio::abm::InfectionState::Infected)], 1);

    auto t3 = mio::abm::TimePoint(0) + mio::abm::hours(24);
    person3.set_infection_state(mio::abm::InfectionState::Infected);
    location.changed_state(person3, mio::abm::InfectionState::Exposed);
    location.add_subpopulations_timepoint(t3);
    auto v3 = location.get_population().get_value(3);
    // Check whether the number of persons in infected state at the location is correct
    ASSERT_EQ(v3[size_t(mio::abm::InfectionState::Infected)], 2);

    // Check total number of subpopulation is correct.
    ASSERT_EQ(location.get_population().get_num_time_points(), 4);
    for (auto&& v_iter : location.get_population()) {
        ASSERT_EQ(v_iter.sum(), 3);
    }
    ASSERT_EQ(location.get_population().get_time(1), 7);
    ASSERT_EQ(location.get_population().get_time(2), 14);
    ASSERT_EQ(location.get_population().get_time(3), 24);
}

TEST(TestLocation, setRequiredMask)
{
    auto location = mio::abm::Location(mio::abm::LocationType::Home, 0);
    ASSERT_EQ(location.get_required_mask(), mio::abm::MaskType::Community);

    location.set_required_mask(mio::abm::MaskType::FFP2);
    ASSERT_EQ(location.get_required_mask(), mio::abm::MaskType::FFP2);
}

TEST(TestLocation, setNPIActive)
{
    auto location = mio::abm::Location(mio::abm::LocationType::Home, 0);
    location.set_npi_active(false);
    ASSERT_FALSE(location.get_npi_active());

    location.set_npi_active(true);
    ASSERT_TRUE(location.get_npi_active());
}