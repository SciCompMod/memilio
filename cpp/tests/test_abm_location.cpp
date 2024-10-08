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

#include "abm/location_id.h"
#include "abm/parameters.h"
#include "abm/person.h"
#include "abm/model.h"
#include "abm_helpers.h"
#include "memilio/utils/random_number_generator.h"

TEST(TestLocation, initCell)
{
    mio::abm::Location location(mio::abm::LocationType::PublicTransport, 0, 6, 2);
    EXPECT_EQ(location.get_cells().size(), 2);
}

TEST(TestLocation, getId)
{
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    ASSERT_EQ(location.get_id(), mio::abm::LocationId(0));
}

TEST(TestLocation, reachCapacity)
{
    using testing::Return;

    auto t     = mio::abm::TimePoint{mio::abm::hours(8).seconds()};
    auto dt    = mio::abm::hours(1);
    auto model = mio::abm::Model(num_age_groups);

    //setup so p1 doesn't do transition
    model.parameters
        .get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        2 * dt.days();
    model.parameters
        .get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        2 * dt.days();
    model.parameters.get<mio::abm::AgeGroupGotoSchool>().set_multiple({age_group_5_to_14}, true);
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    auto home_id   = model.add_location(mio::abm::LocationType::Home);
    auto school_id = model.add_location(mio::abm::LocationType::School);

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

    auto p1 = add_test_person(model, home_id, age_group_5_to_14, mio::abm::InfectionState::InfectedNoSymptoms);
    auto p2 = add_test_person(model, home_id, age_group_5_to_14, mio::abm::InfectionState::Susceptible);

    model.get_person(p1).set_assigned_location(mio::abm::LocationType::School, school_id);
    model.get_person(p2).set_assigned_location(mio::abm::LocationType::School, school_id);
    model.get_person(p1).set_assigned_location(mio::abm::LocationType::Home, home_id);
    model.get_person(p2).set_assigned_location(mio::abm::LocationType::Home, home_id);

    model.get_location(school_id).set_capacity(1, 66);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

    model.evolve(t, dt);

    ASSERT_EQ(model.get_person(p1).get_location(), school_id);
    ASSERT_EQ(model.get_person(p2).get_location(), home_id); // p2 should not be able to enter the school
    ASSERT_EQ(model.get_number_persons(school_id), 1);
    ASSERT_EQ(model.get_number_persons(home_id), 1);
}

TEST(TestLocation, computeSpacePerPersonRelative)
{
    using testing::Return;

    mio::abm::Location home(mio::abm::LocationType::Home, 0, 6, 3);
    home.set_capacity(4, 264, 0); // Capacity for Cell 1
    home.set_capacity(2, 132, 1); // Capacity for Cell 2
    home.set_capacity(0, 0, 2); // Capacity for Cell 3

    auto cells = home.get_cells();
    EXPECT_EQ(cells[0].compute_space_per_person_relative(), 0.25);
    EXPECT_EQ(cells[1].compute_space_per_person_relative(), 0.5);
    EXPECT_EQ(cells[2].compute_space_per_person_relative(), 1.);
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
    std::vector<mio::abm::Person> local_population{infected1, infected2, infected3};

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    auto susceptible = make_test_person(location, age, mio::abm::InfectionState::Susceptible, t, params);
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.5));
    auto person_rng = mio::abm::PersonalRandomNumberGenerator(rng, susceptible);
    interact_testing(person_rng, susceptible, location, local_population, t, dt, params);
    EXPECT_EQ(susceptible.get_infection_state(t + dt), mio::abm::InfectionState::Susceptible);

    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
    interact_testing(person_rng, susceptible, location, local_population, t, dt, params);
    EXPECT_EQ(susceptible.get_infection_state(t + dt), mio::abm::InfectionState::Exposed);
}

TEST(TestLocation, setCapacity)
{
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    location.set_capacity(4, 200);
    EXPECT_EQ(location.get_capacity().persons, (uint32_t)4);
    EXPECT_EQ(location.get_capacity().volume, (uint32_t)200);
}

TEST(TestLocation, setRequiredMask)
{
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    EXPECT_EQ(location.get_required_mask(), mio::abm::MaskType::None);

    location.set_required_mask(mio::abm::MaskType::FFP2);
    EXPECT_EQ(location.get_required_mask(), mio::abm::MaskType::FFP2);
}

TEST(TestLocation, getGeographicalLocation)
{
    auto location                                        = mio::abm::Location(mio::abm::LocationType::Home, 0);
    mio::abm::GeographicalLocation geographical_location = {10.5100470359749, 52.2672785559812};
    location.set_geographical_location(geographical_location);

    EXPECT_EQ(location.get_geographical_location(), geographical_location);
}
