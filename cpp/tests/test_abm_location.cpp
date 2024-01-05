/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn, Khoa Nguyen, Carlotta Gerstein
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
#include "abm/person.h"
#include "abm_helpers.h"
#include "memilio/utils/random_number_generator.h"
#include <memory>

TEST(TestLocation, init)
{
    mio::abm::Location location(mio::abm::LocationType::School, 0, num_age_groups);
    for (mio::abm::InfectionState i = mio::abm::InfectionState(0); i < mio::abm::InfectionState::Count;
         i                          = mio::abm::InfectionState(size_t(i) + 1)) {
        ASSERT_EQ(location.get_subpopulation(mio::abm::TimePoint(0), i), 0);
    }
    EXPECT_EQ(location.get_number_persons(), 0);
}

TEST(TestLocation, copyLocation)
{
    auto location = mio::abm::Location(mio::abm::LocationType::School, 0, num_age_groups);
    auto person   = make_test_person(location, age_group_5_to_14, mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(location.get_number_persons(), 0);
    location.add_person(person);
    EXPECT_EQ(location.get_number_persons(), 1);

    auto copied_location = location.copy_location_without_persons(num_age_groups);
    ASSERT_EQ(copied_location.get_type(), mio::abm::LocationType::School);
    ASSERT_EQ(copied_location.get_index(), location.get_index());
    ASSERT_EQ(copied_location.get_cells().size(), location.get_cells().size());
    EXPECT_EQ(copied_location.get_number_persons(), 0);
}

TEST(TestLocation, initCell)
{
    mio::abm::Location location(mio::abm::LocationType::PublicTransport, 0, 6, 2);
    ASSERT_EQ(location.get_cells().size(), 2);
}

TEST(TestLocation, getIndex)
{
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    ASSERT_EQ((int)location.get_index(), 0);
}

TEST(TestLocation, addRemovePerson)
{
    mio::abm::Location home(mio::abm::LocationType::Home, 0, 6, 1);
    mio::abm::Location location(mio::abm::LocationType::PublicTransport, 0, 6, 3);

    auto person1 = make_test_person(home, age_group_5_to_14, mio::abm::InfectionState::InfectedSymptoms);
    auto person2 = make_test_person(home, age_group_5_to_14, mio::abm::InfectionState::InfectedSymptoms);
    auto person3 = make_test_person(home, age_group_35_to_59, mio::abm::InfectionState::Exposed);

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

    auto rng = mio::RandomNumberGenerator();

    mio::AgeGroup age =
        mio::AgeGroup(mio::UniformIntDistribution<int>::get_instance()(rng, 0, int(num_age_groups - 1)));
    mio::abm::VirusVariant variant = mio::abm::VirusVariant(
        mio::UniformIntDistribution<int>::get_instance()(rng, 0, int(mio::abm::VirusVariant::Count) - 1));

    auto t  = mio::abm::TimePoint(0);
    auto dt = mio::abm::seconds(10000);

    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);

    // setup a location with some chance of exposure
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups, 1);
    mio::abm::Location location(mio::abm::LocationType::PublicTransport, 0, num_age_groups, 3);
    auto infected1     = mio::abm::Person(rng, home, age);
    auto rng_infected1 = mio::abm::Person::RandomNumberGenerator(rng, infected1);
    infected1.add_new_infection(
        mio::abm::Infection(rng_infected1, variant, age, params, t, mio::abm::InfectionState::InfectedNoSymptoms));
    infected1.migrate_to(location, {0});
    auto infected2     = mio::abm::Person(rng, home, age);
    auto rng_infected2 = mio::abm::Person::RandomNumberGenerator(rng, infected2);
    infected2.add_new_infection(
        mio::abm::Infection(rng_infected2, variant, age, params, t, mio::abm::InfectionState::InfectedNoSymptoms));
    infected2.migrate_to(location, {0, 1});

    //cache precomputed results
    location.cache_exposure_rates(t, dt, num_age_groups);

    EXPECT_NEAR((location.get_cells()[0].m_cached_exposure_rate_contacts[{variant, age}]), 0.015015859523894731, 1e-14);
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
    location.cache_exposure_rates(t, dt, num_age_groups);

    EXPECT_NEAR((location.get_cells()[0].m_cached_exposure_rate_air[{variant}]), 0.045047578571684191, 1e-14);
    EXPECT_NEAR((location.get_cells()[1].m_cached_exposure_rate_air[{variant}]), 0.022523789285842095, 1e-14);
    EXPECT_NEAR((location.get_cells()[2].m_cached_exposure_rate_air[{variant}]), 0, 1e-14);
}

TEST(TestLocation, reachCapacity)
{
    using testing::Return;

    auto t     = mio::abm::TimePoint{mio::abm::hours(8).seconds()};
    auto dt    = mio::abm::hours(1);
    auto world = mio::abm::World(num_age_groups);

    //setup so p1 doesn't do transition
    world.parameters
        .get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        2 * dt.days();
    world.parameters
        .get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        2 * dt.days();

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

    auto& p1 = add_test_person(world, home_id, age_group_5_to_14, mio::abm::InfectionState::InfectedNoSymptoms);
    auto& p2 = add_test_person(world, home_id, age_group_5_to_14, mio::abm::InfectionState::Susceptible);

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

    mio::abm::Location home(mio::abm::LocationType::Home, 0, 6, 3);
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

    auto rng = mio::RandomNumberGenerator();

    // Test should work identically work with any age.
    mio::AgeGroup age =
        mio::AgeGroup(mio::UniformIntDistribution<int>::get_instance()(rng, 0, int(num_age_groups - 1)));
    mio::abm::VirusVariant variant = mio::abm::VirusVariant(
        mio::UniformIntDistribution<int>::get_instance()(rng, 0, int(mio::abm::VirusVariant::Count) - 1));

    auto t  = mio::abm::TimePoint(0);
    auto dt = mio::abm::seconds(8640); //0.1 days

    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    params.set_default<mio::abm::ViralLoadDistributions>(num_age_groups);
    params.get<mio::abm::ViralLoadDistributions>()[{variant, age}] = {{1., 1.}, {0.0001, 0.0001}, {-0.0001, -0.0001}};
    params.set_default<mio::abm::InfectivityDistributions>(num_age_groups);
    params.get<mio::abm::InfectivityDistributions>()[{variant, age}] = {{1., 1.}, {1., 1.}};

    // set incubtion period to two days so that the newly infected person is still exposed
    params.get<mio::abm::IncubationPeriod>()[{variant, age}] = 2.;

    //setup location with some chance of exposure
    mio::abm::Location location(mio::abm::LocationType::Work, 0, num_age_groups);
    auto infected1 =
        make_test_person(location, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms, t, params);
    auto infected2 =
        make_test_person(location, age_group_80_plus, mio::abm::InfectionState::InfectedSymptoms, t, params);
    auto infected3 =
        make_test_person(location, age_group_5_to_14, mio::abm::InfectionState::InfectedSymptoms, t, params);

    location.add_person(infected1, {0});
    location.add_person(infected2, {0});
    location.add_person(infected3, {0});

    //cache precomputed results
    location.cache_exposure_rates(t, dt, num_age_groups);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    auto susceptible = make_test_person(location, age, mio::abm::InfectionState::Susceptible, t, params);
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.5));
    auto person_rng = mio::abm::Person::RandomNumberGenerator(rng, susceptible);
    location.interact(person_rng, susceptible, t, dt, params);
    EXPECT_EQ(susceptible.get_infection_state(t + dt), mio::abm::InfectionState::Susceptible);

    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
    location.interact(person_rng, susceptible, t, dt, params);
    EXPECT_EQ(susceptible.get_infection_state(t + dt), mio::abm::InfectionState::Exposed);
}

TEST(TestLocation, setCapacity)
{
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    location.set_capacity(4, 200);
    ASSERT_EQ(location.get_capacity().persons, (uint32_t)4);
    ASSERT_EQ(location.get_capacity().volume, (uint32_t)200);
}

TEST(TestLocation, setRequiredMask)
{
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    ASSERT_EQ(location.get_required_mask(), mio::abm::MaskType::Community);

    location.set_required_mask(mio::abm::MaskType::FFP2);
    ASSERT_EQ(location.get_required_mask(), mio::abm::MaskType::FFP2);
}

TEST(TestLocation, setNPIActive)
{
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    location.set_npi_active(false);
    ASSERT_FALSE(location.get_npi_active());

    location.set_npi_active(true);
    ASSERT_TRUE(location.get_npi_active());
}

TEST(TestLocation, getGeographicalLocation)
{
    auto location                                        = mio::abm::Location(mio::abm::LocationType::Home, 0);
    mio::abm::GeographicalLocation geographical_location = {10.5100470359749, 52.2672785559812};
    location.set_geographical_location(geographical_location);

    ASSERT_EQ(location.get_geographical_location(), geographical_location);
}
