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
#include "abm/infection.h"
#include "abm_helpers.h"
#include <memory>

TEST(TestLocation, init)
{
    auto location = mio::abm::Location(mio::abm::LocationType::School, 0);
    for (mio::abm::InfectionState i = mio::abm::InfectionState(0); i < mio::abm::InfectionState::Count;
         i                          = mio::abm::InfectionState(size_t(i) + 1)) {
        ASSERT_EQ(location.get_subpopulation(mio::abm::TimePoint(0), i), 0);
    }
    location.initialize_subpopulations(mio::abm::TimePoint(0));
    ASSERT_EQ(print_wrap(location.get_subpopulations().get_last_value()),
              print_wrap(mio::TimeSeries<double>::Vector::Zero((size_t)mio::abm::InfectionState::Count)));
    EXPECT_EQ(location.get_number_persons(), 0);
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
    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0, 1);
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 3);

    auto person1 = make_test_person(home, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::InfectedSymptoms);
    auto person2 = make_test_person(home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms);
    auto person3 = make_test_person(home, mio::abm::AgeGroup::Age35to59, mio::abm::InfectionState::Exposed);

    home.add_person(person1, {0});
    home.add_person(person2, {0});
    home.add_person(person3, {0});

    person1.migrate_to(location, {0, 1});
    person2.migrate_to(location, {0});
    person3.migrate_to(location, {0, 1});

    auto t = mio::abm::TimePoint(0);
    ASSERT_EQ(home.get_number_persons(), 0u);
    ASSERT_EQ(location.get_subpopulation(t, mio::abm::InfectionState::InfectedSymptoms), 2);
    ASSERT_EQ(location.get_subpopulation(t, mio::abm::InfectionState::Exposed), 1);
    ASSERT_EQ(location.get_cells()[0].m_persons.size(), 3u);
    ASSERT_EQ(location.get_cells()[1].m_persons.size(), 2u);
    ASSERT_EQ(location.get_cells()[2].m_persons.size(), 0u);

    location.remove_person(person2);

    EXPECT_EQ(location.get_number_persons(), 2u);
    ASSERT_EQ(location.get_subpopulation(t, mio::abm::InfectionState::InfectedSymptoms), 1);
    ASSERT_EQ(location.get_subpopulation(t, mio::abm::InfectionState::Exposed), 1);
    ASSERT_EQ(location.get_cells()[0].m_persons.size(), 2u);
    ASSERT_EQ(location.get_cells()[1].m_persons.size(), 2u);
    ASSERT_EQ(location.get_cells()[2].m_persons.size(), 0u);
}

TEST(TestLocation, CacheExposureRate)
{
    using testing::Return;

    {
        mio::abm::AgeGroup age =
            mio::abm::AgeGroup(mio::UniformIntDistribution<int>()(0, int(mio::abm::AgeGroup::Count) - 1));
        mio::abm::VirusVariant variant =
            mio::abm::VirusVariant(mio::UniformIntDistribution<int>()(0, int(mio::abm::VirusVariant::Count) - 1));

        auto t  = mio::abm::TimePoint(0);
        auto dt = mio::abm::seconds(10000);

        mio::abm::GlobalInfectionParameters params;

        // setup a location with some chance of exposure
        auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0, 1);
        auto location  = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 3);
        auto infected1 = mio::abm::Person(home, age);
        infected1.add_new_infection(
            mio::abm::Infection(variant, age, params, t, mio::abm::InfectionState::InfectedNoSymptoms));
        infected1.migrate_to(location, {0});
        auto infected2 = mio::abm::Person(home, age);
        infected2.add_new_infection(
            mio::abm::Infection(variant, age, params, t, mio::abm::InfectionState::InfectedNoSymptoms));
        infected2.migrate_to(location, {0, 1});

        //cache precomputed results
        location.cache_exposure_rates(t, dt);

        EXPECT_NEAR((location.get_cells()[0].m_cached_exposure_rate_contacts[{variant, age}]), 0.015015859523894731,
                    1e-14);
        EXPECT_NEAR((location.get_cells()[0].m_cached_exposure_rate_air[{variant}]), 0.015015859523894731, 1e-14);
        EXPECT_NEAR((location.get_cells()[1].m_cached_exposure_rate_contacts[{variant, age}]), 0.0075079297619473654,
                    1e-14);
        EXPECT_NEAR((location.get_cells()[1].m_cached_exposure_rate_air[{variant}]), 0.0075079297619473654, 1e-14);
        EXPECT_NEAR((location.get_cells()[2].m_cached_exposure_rate_contacts[{variant, age}]), 0, 1e-14);
        EXPECT_NEAR((location.get_cells()[2].m_cached_exposure_rate_air[{variant}]), 0, 1e-14);

        // should also work with capacities
        location.set_capacity_adapted_transmission_risk_flag(true);
        location.set_capacity(2, 22, 0); // Capacity for Cell 1
        location.set_capacity(2, 22, 1); // Capacity for Cell 2
        location.set_capacity(2, 22, 2); // Capacity for Cell 3
        location.cache_exposure_rates(t, dt);

        EXPECT_NEAR((location.get_cells()[0].m_cached_exposure_rate_air[{variant}]), 0.045047578571684191, 1e-14);
        EXPECT_NEAR((location.get_cells()[1].m_cached_exposure_rate_air[{variant}]), 0.022523789285842095, 1e-14);
        EXPECT_NEAR((location.get_cells()[2].m_cached_exposure_rate_air[{variant}]), 0, 1e-14);
    }
}

TEST(TestLocation, reachCapacity)
{
    using testing::Return;

    auto t      = mio::abm::TimePoint{mio::abm::hours(8).seconds()};
    auto dt     = mio::abm::hours(1);
    auto params = mio::abm::GlobalInfectionParameters{};
    //setup so p1 doesn't transition
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, mio::abm::VaccinationState::Unvaccinated}] =
        2 * dt.days();
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, mio::abm::VaccinationState::Unvaccinated}] =
        2 * dt.days();

    auto world     = mio::abm::World(params);
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
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillRepeatedly(testing::Return(1.0));

    auto& p1 =
        add_test_person(world, home_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::InfectedNoSymptoms);
    auto& p2 = add_test_person(world, home_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::Susceptible);

    auto& home   = world.get_individualized_location(home_id);
    auto& school = world.get_individualized_location(school_id);

    p1.set_assigned_location(school_id);
    p2.set_assigned_location(school_id);
    p1.set_assigned_location(home_id);
    p2.set_assigned_location(home_id);

    school.set_capacity(1, 66);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

    world.evolve(t, dt);

    ASSERT_EQ(p1.get_location(), school);
    ASSERT_EQ(p2.get_location(), home); // p2 should not be able to enter the school
    ASSERT_EQ(school.get_number_persons(), 1);
    ASSERT_EQ(home.get_number_persons(), 1);
}

TEST(TestLocation, computeSpacePerPersonRelative)
{
    using testing::Return;

    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0, 3);
    home.set_capacity(4, 264, 0); // Capacity for Cell 1
    home.set_capacity(2, 132, 1); // Capacity for Cell 2
    home.set_capacity(0, 0, 2); // Capacity for Cell 3

    auto cells = home.get_cells();
    ASSERT_EQ(cells[0].compute_space_per_person_relative(), 0.25);
    ASSERT_EQ(cells[1].compute_space_per_person_relative(), 0.5);
    ASSERT_EQ(cells[2].compute_space_per_person_relative(), 1.);
}

TEST(TestLocation, interact)
{
    using testing::Return;

    // Test should work identically work with any age.
    mio::abm::AgeGroup age =
        mio::abm::AgeGroup(mio::UniformIntDistribution<int>()(0, int(mio::abm::AgeGroup::Count) - 1));
    mio::abm::VirusVariant variant =
        mio::abm::VirusVariant(mio::UniformIntDistribution<int>()(0, int(mio::abm::VirusVariant::Count) - 1));

    auto t  = mio::abm::TimePoint(0);
    auto dt = mio::abm::seconds(8640); //0.1 days

    mio::abm::GlobalInfectionParameters params;
    params.set_default<mio::abm::ViralLoadDistributions>();
    params.get<mio::abm::ViralLoadDistributions>()[{variant, age, mio::abm::VaccinationState::Unvaccinated}] = {
        {1., 1.}, {0.0001, 0.0001}, {-0.0001, -0.0001}};
    params.set_default<mio::abm::InfectivityDistributions>();
    params.get<mio::abm::InfectivityDistributions>()[{variant, age}] = {{1., 1.}, {1., 1.}};

    // set incubtion period to two days so that the newly infected person is still exposed
    params.get<mio::abm::IncubationPeriod>()[{variant, age, mio::abm::VaccinationState::Unvaccinated}] = 2.;

    //setup location with some chance of exposure
    auto location  = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto infected1 = make_test_person(location, mio::abm::AgeGroup::Age15to34,
                                      mio::abm::InfectionState::InfectedNoSymptoms, t, params);
    auto infected2 = make_test_person(location, mio::abm::AgeGroup::Age80plus,
                                      mio::abm::InfectionState::InfectedSymptoms, t, params);
    auto infected3 =
        make_test_person(location, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::InfectedSymptoms, t, params);

    location.add_person(infected1, {0});
    location.add_person(infected2, {0});
    location.add_person(infected3, {0});

    //cache precomputed results
    location.cache_exposure_rates(t, dt);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    auto susceptible = make_test_person(location, age, mio::abm::InfectionState::Susceptible, t, params);
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.5));
    location.interact(susceptible, t, dt, params);
    EXPECT_EQ(susceptible.get_infection_state(t + dt), mio::abm::InfectionState::Susceptible);

    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
    location.interact(susceptible, t, dt, params);
    EXPECT_EQ(susceptible.get_infection_state(t + dt), mio::abm::InfectionState::Exposed);
}

TEST(TestLocation, setCapacity)
{
    auto location = mio::abm::Location(mio::abm::LocationType::Home, 0);
    location.set_capacity(4, 200);
    ASSERT_EQ(location.get_capacity().persons, (uint32_t)4);
    ASSERT_EQ(location.get_capacity().volume, (uint32_t)200);
}

TEST(TestLocation, storeSubpopulations)
{
    auto t      = mio::abm::TimePoint(0);
    auto dt     = mio::abm::days(7);
    auto params = mio::abm::GlobalInfectionParameters{};

    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 3);

    //setup: p1 goes from Infected to Recovered, p2 stays in Infected and p3 goes from Exposed to InfectedNoSymptoms to Recovered
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age5to14,
                                                         mio::abm::VaccinationState::Unvaccinated}] = 1.5 * dt.days();

    params.get<mio::abm::InfectedSymptomsToRecovered>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, mio::abm::VaccinationState::Unvaccinated}] =
        5 * dt.days();
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
                                                      mio::abm::VaccinationState::Unvaccinated}] = 5 * dt.days();

    params.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age35to59,
                                              mio::abm::VaccinationState::Unvaccinated}] = 0.4 * dt.days();
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age35to59, mio::abm::VaccinationState::Unvaccinated}] =
        1.8 * dt.days();

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;

    // mock person 1
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillOnce(testing::Return(0.6)) // transition to Recovered
        .WillRepeatedly(testing::Return(1.0));

    auto person1 =
        make_test_person(location, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::InfectedSymptoms, t, params);
    location.add_person(person1, {0});

    // mock person 2 not needed due to high setup of transition times
    auto person2 = make_test_person(location, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms,
                                    t, params);
    location.add_person(person2, {0});

    // mock person 3
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillOnce(testing::Return(0.6)) // transition to Recovered
        .WillRepeatedly(testing::Return(1.0));
    auto person3 =
        make_test_person(location, mio::abm::AgeGroup::Age35to59, mio::abm::InfectionState::Exposed, t, params);
    location.add_person(person3, {0});

    location.initialize_subpopulations(t);
    auto t1 = t + dt;
    location.store_subpopulations(t1);
    auto v1 = location.get_subpopulations().get_value(1);
    // Check whether the number of persons in infected state at the location is correct
    ASSERT_EQ(v1[size_t(mio::abm::InfectionState::InfectedSymptoms)], 2);
    ASSERT_EQ(v1[size_t(mio::abm::InfectionState::InfectedNoSymptoms)], 1);

    auto t2 = t1 + dt;
    location.store_subpopulations(t2);
    auto v2 = location.get_subpopulations().get_value(2);
    // Check whether the number of persons in infected state at the location is correct
    ASSERT_EQ(v2[size_t(mio::abm::InfectionState::InfectedSymptoms)], 1);
    ASSERT_EQ(v2[size_t(mio::abm::InfectionState::Recovered)], 1);
    ASSERT_EQ(v2[size_t(mio::abm::InfectionState::InfectedNoSymptoms)], 1);

    auto t3 = t2 + mio::abm::days(10);
    location.store_subpopulations(t3);
    auto v3 = location.get_subpopulations().get_value(3);
    // Check whether the number of persons in infected state at the location is correct
    ASSERT_EQ(v3[size_t(mio::abm::InfectionState::InfectedSymptoms)], 1);
    ASSERT_EQ(v3[size_t(mio::abm::InfectionState::Recovered)], 2);

    // Check total number of subpopulation is correct.
    ASSERT_EQ(location.get_subpopulations().get_num_time_points(), 4);
    for (auto&& v_iter : location.get_subpopulations()) {
        ASSERT_EQ(v_iter.sum(), 3);
    }
    ASSERT_EQ(location.get_subpopulations().get_time(1), 7);
    ASSERT_EQ(location.get_subpopulations().get_time(2), 14);
    ASSERT_EQ(location.get_subpopulations().get_time(3), 24);
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
