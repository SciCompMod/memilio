/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "abm_helpers.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "random_number_test.h"

using TestLocation = RandomNumberTest;

/**
 * @brief Test that initializing a location with cells correctly creates the given number of cells.
 */
TEST_F(TestLocation, initCell)
{
    // Create a location of type PublicTransport with 2 cells.
    mio::abm::Location location(mio::abm::LocationType::PublicTransport, 0, 6, 2);
    // Verify that the number of cells created is as expected (2).
    EXPECT_EQ(location.get_cells().size(), 2);
}

/**
 * @brief Test that a location correctly returns its ID.
 */
TEST_F(TestLocation, getId)
{
    // Create a location of type Home with an ID of 0.
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    // Verify that the location's ID is correctly set to 0.
    EXPECT_EQ(location.get_id(), mio::abm::LocationId(0));
}

/**
 * @brief Test that the computation of space per person relative to capacity works correctly.
 */
TEST_F(TestLocation, computeSpacePerPersonRelative)
{
    using testing::Return;

    // Create a location of type Home with 3 cells.
    mio::abm::Location home(mio::abm::LocationType::Home, 0, 6, 3);
    home.set_capacity(4, 264, 0); // Capacity for Cell 1
    home.set_capacity(2, 132, 1); // Capacity for Cell 2
    home.set_capacity(0, 0, 2); // Capacity for Cell 3

    auto cells = home.get_cells();

    // Verify the space per person relative for each cell.
    EXPECT_EQ(cells[0].compute_space_per_person_relative(), 0.25);
    EXPECT_EQ(cells[1].compute_space_per_person_relative(), 0.5);
    EXPECT_EQ(cells[2].compute_space_per_person_relative(), 1.);
}

/**
 * @brief Test the interaction between infected and susceptible persons at a location.
 */
TEST_F(TestLocation, interact)
{
    using testing::Return;

    // Test should work identically work with any age.
    mio::AgeGroup age = mio::AgeGroup(this->random_integer(0, int(num_age_groups - 1)));
    mio::abm::VirusVariant variant =
        mio::abm::VirusVariant(this->random_integer(0, int(mio::abm::VirusVariant::Count) - 1));

    auto t  = mio::abm::TimePoint(0);
    auto dt = mio::abm::seconds(8640); //0.1 days

    // Setup model parameters for viral loads and infectivity distributions.
    // Setup model parameters for viral loads and infectivity distributions.
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    params.set_default<mio::abm::ViralLoadDistributions>(num_age_groups);
    params.get<mio::abm::ViralLoadDistributions>()[{variant, age}] = {{1., 1.}, {0.0001, 0.0001}, {-0.0001, -0.0001}};
    params.set_default<mio::abm::InfectivityDistributions>(num_age_groups);
    params.get<mio::abm::InfectivityDistributions>()[{variant, age}] = {{1., 1.}, {1., 1.}};

    // Set incubtion period to two days so that the newly infected person is still exposed
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::LogNormalDistribution<double>>>> mock_logNorm_dist;
    EXPECT_CALL(mock_logNorm_dist.get_mock(), invoke).WillRepeatedly(testing::Return(2));

    // Setup location with some chance of exposure
    mio::abm::Location location(mio::abm::LocationType::Work, 0, num_age_groups);
    auto infected1 = make_test_person(this->get_rng(), location, age_group_15_to_34,
                                      mio::abm::InfectionState::InfectedNoSymptoms, t, params);
    auto infected2 = make_test_person(this->get_rng(), location, age_group_80_plus,
                                      mio::abm::InfectionState::InfectedSymptoms, t, params);
    auto infected3 = make_test_person(this->get_rng(), location, age_group_5_to_14,
                                      mio::abm::InfectionState::InfectedSymptoms, t, params);
    std::vector<mio::abm::Person> local_population{infected1, infected2, infected3};

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    // Create a susceptible person and test interaction.
    auto susceptible =
        make_test_person(this->get_rng(), location, age, mio::abm::InfectionState::Susceptible, t, params);
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.5)); // Probability of no infection
    auto person_rng = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), susceptible);
    interact_testing(person_rng, susceptible, location, local_population, t, dt, params);
    EXPECT_EQ(susceptible.get_infection_state(t + dt), mio::abm::InfectionState::Susceptible);

    // Test with a higher probability of infection leading to an exposed state.
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05)); // Probability of infection
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0)); // Choose infection state
    interact_testing(person_rng, susceptible, location, local_population, t, dt, params);
    EXPECT_EQ(susceptible.get_infection_state(t + dt), mio::abm::InfectionState::Exposed);
}

/**
 * @brief Test setting and getting the capacity of a location.
 */
TEST_F(TestLocation, setCapacity)
{
    // Create a location of type Home.
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    // Set the capacity of the location and verify.
    location.set_capacity(4, 200);
    EXPECT_EQ(location.get_capacity().persons, (uint32_t)4);
    EXPECT_EQ(location.get_capacity().volume, (uint32_t)200);
}

/**
 * @brief Test setting and getting the required mask type at a location.
 */
TEST_F(TestLocation, setRequiredMask)
{
    // Create a location of type Home.
    mio::abm::Location location(mio::abm::LocationType::Home, 0, num_age_groups);
    // Verify that the default required mask is set to None.
    EXPECT_EQ(location.get_required_mask(), mio::abm::MaskType::None);

    // Set a new required mask type and verify the change.
    location.set_required_mask(mio::abm::MaskType::FFP2);
    EXPECT_EQ(location.get_required_mask(), mio::abm::MaskType::FFP2);
}

/**
 * @brief Test setting and getting the geographical location of a location.
 */
TEST_F(TestLocation, getGeographicalLocation)
{
    // Create a location of type Home.
    auto location = mio::abm::Location(mio::abm::LocationType::Home, 0);
    // Set a geographical location for the location.
    mio::abm::GeographicalLocation geographical_location = {10.5100470359749, 52.2672785559812};
    location.set_geographical_location(geographical_location);
    // Verify that the set geographical location matches the expected values.
    EXPECT_EQ(location.get_geographical_location(), geographical_location);
}
