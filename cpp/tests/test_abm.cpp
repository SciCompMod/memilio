/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth
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
#include "abm/abm.h"
#include "abm/migration_rules.h"
#include "abm/lockdown_rules.h"
#include "memilio/math/eigen_util.h"
#include "matchers.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <memory>

TEST(TestLocation, init)
{
    auto location = mio::Location(mio::LocationType::School, 0);
    for (mio::InfectionState i = mio::InfectionState(0); i < mio::InfectionState::Count;
         i                     = mio::InfectionState(size_t(i) + 1)) {
        ASSERT_EQ(location.get_subpopulation(i), 0);
    }
    ASSERT_EQ(print_wrap(location.get_subpopulations()),
              print_wrap(Eigen::VectorXi::Zero(Eigen::Index(mio::InfectionState::Count))));
}

TEST(TestLocation, initCell)
{
    auto location = mio::Location(mio::LocationType::PublicTransport, 0, 2);
    ASSERT_EQ(location.get_cells().size(), 2);
}

TEST(TestLocation, GetIndex)
{
    auto location = mio::Location(mio::LocationType::Home, 0);
    ASSERT_EQ((int)location.get_index(), 0);
}

TEST(TestLocation, addRemovePerson)
{
    auto home     = mio::Location(mio::LocationType::Home, 0, 0);
    auto location = mio::Location(mio::LocationType::PublicTransport, 0, 3);
    auto person1  = mio::Person(home, mio::InfectionState::Infected, mio::AbmAgeGroup::Age5to14, {});
    home.add_person(person1);
    person1.migrate_to(home, location, {0, 1});
    auto person2 = mio::Person(home, mio::InfectionState::Infected, mio::AbmAgeGroup::Age15to34, {});
    home.add_person(person2);
    person2.migrate_to(home, location, {0});
    auto person3 = mio::Person(home, mio::InfectionState::Exposed, mio::AbmAgeGroup::Age35to59, {});
    home.add_person(person3);
    person3.migrate_to(home, location, {0, 1});

    ASSERT_EQ(location.get_subpopulation(mio::InfectionState::Infected), 2);
    ASSERT_EQ(location.get_subpopulation(mio::InfectionState::Exposed), 1);
    ASSERT_EQ(location.get_cells()[0].num_people, 3u);
    ASSERT_EQ(location.get_cells()[1].num_people, 2u);
    ASSERT_EQ(location.get_cells()[2].num_people, 0u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 2u);
    ASSERT_EQ(location.get_cells()[1].num_infected, 1u);
    ASSERT_EQ(location.get_cells()[2].num_infected, 0u);

    location.remove_person(person2);

    ASSERT_EQ(location.get_subpopulation(mio::InfectionState::Infected), 1);
    ASSERT_EQ(location.get_subpopulation(mio::InfectionState::Exposed), 1);
    ASSERT_EQ(location.get_cells()[0].num_people, 2u);
    ASSERT_EQ(location.get_cells()[1].num_people, 2u);
    ASSERT_EQ(location.get_cells()[2].num_people, 0u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 1u);
    ASSERT_EQ(location.get_cells()[1].num_infected, 1u);
    ASSERT_EQ(location.get_cells()[2].num_infected, 0u);
}

TEST(TestLocation, setTestingScheme)
{
    auto location = mio::Location(mio::LocationType::Home, 0);
    location.set_testing_scheme(mio::days(5), 0.9);
    ASSERT_EQ(location.get_testing_scheme().get_interval(), mio::days(5));
    ASSERT_EQ(location.get_testing_scheme().get_probability(), 0.9);
}

/**
 * mock of the generator function of DistributionAdapter<DistT>.
 * can't be used directly as a generator function because it is not copyable.
 * see MockDistributionRef
 */
template <class DistT>
struct MockDistribution {
    using Distribution = DistT;
    // using invoke() instead of operator() because operators cant be mocked in the GMock framework */
    MOCK_METHOD(typename Distribution::ResultType, invoke, (const typename Distribution::ParamType&), ());
};

/**
 * reference wrapper of a MockDistribution object.
 * Mocks are not copyable but the generator function of a distribution must be copyable.
 * This wrapper is copyable and all copies redirect invocations to a shared underlying mock
 * so it can be used as a generator function.
 */
template <class MockDistribution>
struct MockDistributionRef {
    using Distribution = typename MockDistribution::Distribution;
    typename Distribution::ResultType operator()(const typename Distribution::ParamType& p)
    {
        return mock->invoke(p);
    }
    std::shared_ptr<MockDistribution> mock = std::make_shared<MockDistribution>();
};

/**
 * Replaces the generator function in the static instance of DistributionAdapter with a mock.
 * On construction sets the generator and on destruction restores the previous generator.
 */
template <class MockDistribution>
struct ScopedMockDistribution {
    using Distribution = typename MockDistribution::Distribution;
    /**
     * constructor replaces the generator function with a mock
     */
    ScopedMockDistribution()
    {
        old = Distribution::get_instance().get_generator();
        Distribution::get_instance().set_generator(mock_ref);
    }
    ~ScopedMockDistribution()
    {
        Distribution::get_instance().set_generator(old);
    }
    MockDistribution& get_mock()
    {
        return *mock_ref.mock;
    }

    MockDistributionRef<MockDistribution> mock_ref;
    typename Distribution::GeneratorFunction old;
};

TEST(TestPerson, init)
{
    using testing::Return;

    auto location = mio::Location(mio::LocationType::Work, 0);
    //setup rng mock so the time_until_carrier is 1.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformIntDistribution<int>>>>
        mock_uniform_int_dist;
    EXPECT_CALL(mock_uniform_int_dist.get_mock(), invoke).Times(testing::AtLeast(1)).WillRepeatedly(Return(0));

    auto person = mio::Person(location, mio::InfectionState::Exposed, mio::AbmAgeGroup::Age60to79, {});
    location.add_person(person);

    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Exposed);
    ASSERT_EQ(person.get_location_id().index, location.get_index());
    ASSERT_EQ(person.get_location_id().type, location.get_type());
    ASSERT_EQ(person.get_person_id(), mio::INVALID_PERSON_ID);

    auto person2 = mio::Person(location, mio::InfectionState::Exposed, mio::AbmAgeGroup::Age60to79, {},
                               mio::VaccinationState::Unvaccinated, 0);
    ASSERT_EQ(person2.get_person_id(), 0u);

    mio::TimeSpan dt = mio::hours(1);
    person.interact(dt, {}, location, {});
    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Carrier);
}

TEST(TestPerson, migrate)
{
    auto home   = mio::Location(mio::LocationType::Home, 0, 0);
    auto loc1   = mio::Location(mio::LocationType::PublicTransport, 0, 1);
    auto loc2   = mio::Location(mio::LocationType::School, 0);
    auto loc3   = mio::Location(mio::LocationType::PublicTransport, 0, 2);
    auto person = mio::Person(home, mio::InfectionState::Recovered_Carrier, mio::AbmAgeGroup::Age0to4, {});
    home.add_person(person);
    person.migrate_to(home, loc1, {0});

    ASSERT_EQ(person.get_location_id().index, loc1.get_index());
    ASSERT_EQ(person.get_location_id().type, loc1.get_type());
    ASSERT_EQ(loc1.get_subpopulation(mio::InfectionState::Recovered_Carrier), 1);
    ASSERT_EQ(home.get_subpopulation(mio::InfectionState::Recovered_Carrier), 0);
    ASSERT_EQ(loc1.get_cells()[0].num_people, 1u);

    person.migrate_to(loc1, loc2);

    ASSERT_EQ(person.get_location_id().index, loc2.get_index());
    ASSERT_EQ(person.get_location_id().type, loc2.get_type());
    ASSERT_EQ(loc2.get_subpopulation(mio::InfectionState::Recovered_Carrier), 1);
    ASSERT_EQ(loc1.get_subpopulation(mio::InfectionState::Recovered_Carrier), 0);
    ASSERT_EQ(loc1.get_cells()[0].num_people, 0u);

    person.migrate_to(loc2, loc3, {0, 1});

    ASSERT_EQ(loc3.get_cells()[0].num_people, 1u);
    ASSERT_EQ(loc3.get_cells()[1].num_people, 1u);
    ASSERT_EQ(person.get_cells().size(), 2);
    ASSERT_EQ(person.get_cells()[0], 0u);
    ASSERT_EQ(person.get_cells()[1], 1u);
}

TEST(TestPerson, setGetAssignedLocation)
{
    auto location = mio::Location(mio::LocationType::Work, 2);
    auto person   = mio::Person(location, mio::InfectionState::Recovered_Carrier, mio::AbmAgeGroup::Age60to79, {});
    person.set_assigned_location(location);
    ASSERT_EQ((int)person.get_assigned_location_index(mio::LocationType::Work), 2);

    person.set_assigned_location({4, mio::LocationType::Work});
    ASSERT_EQ((int)person.get_assigned_location_index(mio::LocationType::Work), 4);
}

TEST(TestWorld, findLocation)
{
    auto world     = mio::World();
    auto home_id   = world.add_location(mio::LocationType::Home);
    auto school_id = world.add_location(mio::LocationType::School);
    auto work_id   = world.add_location(mio::LocationType::Work);
    auto person    = mio::Person(home_id, mio::InfectionState::Recovered_Carrier, mio::AbmAgeGroup::Age60to79,
                              world.get_global_infection_parameters());
    auto& home     = world.get_individualized_location(home_id);
    auto& school   = world.get_individualized_location(school_id);
    auto& work     = world.get_individualized_location(work_id);
    person.set_assigned_location(home);
    person.set_assigned_location(school);
    person.set_assigned_location({0, mio::LocationType::Work});

    ASSERT_EQ(world.find_location(mio::LocationType::Work, person), &work);
    ASSERT_EQ(world.find_location(mio::LocationType::School, person), &school);
    ASSERT_EQ(world.find_location(mio::LocationType::Home, person), &home);
}

TEST(TestLocation, beginStep)
{
    using testing::Return;

    // Test should work identically work with any age.
    mio::AbmAgeGroup age = mio::AbmAgeGroup(mio::UniformIntDistribution<int>()(0, int(mio::AbmAgeGroup::Count) - 1));
    mio::VaccinationState vaccination_state =
        mio::VaccinationState(mio::UniformIntDistribution<int>()(0, int(mio::VaccinationState::Count) - 1));

    mio::GlobalInfectionParameters params;
    params.set<mio::CarrierToInfected>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::CarrierToInfected>()[{age, vaccination_state}] = 0.5;
    params.set<mio::CarrierToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::CarrierToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::DetectInfection>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::DetectInfection>()[{age, vaccination_state}] = 0.5;
    params.set<mio::InfectedToSevere>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::InfectedToSevere>()[{age, vaccination_state}] = 0.5;
    params.set<mio::InfectedToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::InfectedToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::SevereToCritical>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::SevereToCritical>()[{age, vaccination_state}] = 0.5;
    params.set<mio::SevereToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::SevereToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::CriticalToDead>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::CriticalToDead>()[{age, vaccination_state}] = 0.5;
    params.set<mio::CriticalToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::CriticalToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::RecoveredToSusceptible>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::RecoveredToSusceptible>()[{age, vaccination_state}] = 0.5;
    params.set<mio::SusceptibleToExposedByCarrier>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::SusceptibleToExposedByCarrier>()[{age, vaccination_state}] = 0.4;
    params.set<mio::SusceptibleToExposedByInfected>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::SusceptibleToExposedByInfected>()[{age, vaccination_state}] = 0.5;

    //setup location with some chance of exposure
    auto home      = mio::Location(mio::LocationType::Home, 0, 0);
    auto location1 = mio::Location(mio::LocationType::PublicTransport, 0, 3);
    auto infected1 =
        mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34, params, vaccination_state);
    home.add_person(infected1);
    infected1.migrate_to(home, location1, {0});
    auto infected2 =
        mio::Person(home, mio::InfectionState::Infected, mio::AbmAgeGroup::Age80plus, params, vaccination_state);
    home.add_person(infected2);
    infected2.migrate_to(home, location1, {0, 1});
    auto infected3 =
        mio::Person(home, mio::InfectionState::Infected, mio::AbmAgeGroup::Age5to14, params, vaccination_state);
    home.add_person(infected3);
    infected3.migrate_to(home, location1, {1});

    //cache precomputed results
    auto dt = mio::seconds(8640);
    location1.begin_step(dt, params);

    ASSERT_TRUE(std::abs(location1.get_cells()[0].cached_exposure_rate[{age, vaccination_state}] - 0.9) < 0.001);
    ASSERT_TRUE(std::abs(location1.get_cells()[1].cached_exposure_rate[{age, vaccination_state}] - 1) < 0.001);
    ASSERT_TRUE(std::abs(location1.get_cells()[2].cached_exposure_rate[{age, vaccination_state}]) < 0.001);
}

TEST(TestLocation, changedState)
{
    auto home     = mio::Location(mio::LocationType::Home, 0, 0);
    auto location = mio::Location(mio::LocationType::PublicTransport, 0, 1);
    auto p1       = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34, {});
    home.add_person(p1);
    p1.migrate_to(home, location, {0});
    auto p2 = mio::Person(home, mio::InfectionState::Infected, mio::AbmAgeGroup::Age80plus, {});
    home.add_person(p2);
    p2.migrate_to(home, location, {0});
    auto p3 = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age80plus, {});
    home.add_person(p3);
    p3.migrate_to(home, location, {0});

    ASSERT_EQ(location.get_cells()[0].num_carriers, 1u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 1u);
    location.changed_state(p1, mio::InfectionState::Susceptible);
    ASSERT_EQ(location.get_cells()[0].num_carriers, 2u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 1u);
    location.changed_state(p2, mio::InfectionState::Carrier);
    ASSERT_EQ(location.get_cells()[0].num_carriers, 1u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 2u);
    location.changed_state(p3, mio::InfectionState::Infected);
    ASSERT_EQ(location.get_cells()[0].num_carriers, 1u);
    ASSERT_EQ(location.get_cells()[0].num_infected, 1u);
}

TEST(TestLocation, interact)
{
    using testing::Return;

    // Test should work identically work with any age.
    mio::AbmAgeGroup age = mio::AbmAgeGroup(mio::UniformIntDistribution<int>()(0, int(mio::AbmAgeGroup::Count) - 1));
    mio::VaccinationState vaccination_state =
        mio::VaccinationState(mio::UniformIntDistribution<int>()(0, int(mio::VaccinationState::Count) - 1));

    mio::GlobalInfectionParameters params;
    params.set<mio::CarrierToInfected>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::CarrierToInfected>()[{age, vaccination_state}] = 0.5;
    params.set<mio::CarrierToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::CarrierToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::DetectInfection>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::DetectInfection>()[{age, vaccination_state}] = 0.5;
    params.set<mio::InfectedToSevere>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::InfectedToSevere>()[{age, vaccination_state}] = 0.5;
    params.set<mio::InfectedToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::InfectedToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::SevereToCritical>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::SevereToCritical>()[{age, vaccination_state}] = 0.5;
    params.set<mio::SevereToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::SevereToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::CriticalToDead>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::CriticalToDead>()[{age, vaccination_state}] = 0.5;
    params.set<mio::CriticalToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::CriticalToRecovered>()[{age, vaccination_state}] = 0.5;
    params.set<mio::RecoveredToSusceptible>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::RecoveredToSusceptible>()[{age, vaccination_state}] = 0.5;
    params.set<mio::SusceptibleToExposedByCarrier>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::SusceptibleToExposedByCarrier>()[{age, vaccination_state}] = 0.5;
    params.set<mio::SusceptibleToExposedByInfected>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.});
    params.get<mio::SusceptibleToExposedByInfected>()[{age, vaccination_state}] = 0.5;

    //setup location with some chance of exposure
    auto location = mio::Location(mio::LocationType::Work, 0);
    auto infected1 =
        mio::Person(location, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34, params, vaccination_state);
    location.add_person(infected1);
    auto infected2 =
        mio::Person(location, mio::InfectionState::Infected, mio::AbmAgeGroup::Age80plus, params, vaccination_state);
    location.add_person(infected2);
    auto infected3 =
        mio::Person(location, mio::InfectionState::Infected, mio::AbmAgeGroup::Age5to14, params, vaccination_state);
    location.add_person(infected3);

    //cache precomputed results
    auto dt = mio::seconds(8640); //0.1 days
    location.begin_step(dt, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    {
        auto susceptible = mio::Person(location, mio::InfectionState::Susceptible, age, params, vaccination_state);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(susceptible, dt, params), mio::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.15));
        EXPECT_EQ(location.interact(susceptible, dt, params), mio::InfectionState::Susceptible);
    }

    {
        auto exposed = mio::Person(location, mio::InfectionState::Exposed, age, params, vaccination_state);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(0); //no transitions out of exposed state
        EXPECT_EQ(location.interact(exposed, dt, params), mio::InfectionState::Exposed);
    }

    {
        auto carrier = mio::Person(location, mio::InfectionState::Carrier, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(carrier, dt, params), mio::InfectionState::Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.099));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(carrier, dt, params), mio::InfectionState::Recovered_Carrier);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(carrier, dt, params), mio::InfectionState::Carrier);
    }

    {
        auto infected = mio::Person(location, mio::InfectionState::Infected, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(infected, dt, params), mio::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(infected, dt, params), mio::InfectionState::Infected_Severe);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(infected, dt, params), mio::InfectionState::Infected);
    }

    {
        auto severe = mio::Person(location, mio::InfectionState::Infected_Severe, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(severe, dt, params), mio::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(severe, dt, params), mio::InfectionState::Infected_Critical);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(severe, dt, params), mio::InfectionState::Infected_Severe);
    }

    {
        auto critical = mio::Person(location, mio::InfectionState::Infected_Critical, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(critical, dt, params), mio::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(critical, dt, params), mio::InfectionState::Dead);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(critical, dt, params), mio::InfectionState::Infected_Critical);
    }

    for (auto&& recovered_state : {mio::InfectionState::Recovered_Carrier, mio::InfectionState::Recovered_Infected}) {
        auto recovered = mio::Person(location, recovered_state, age, params, vaccination_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(recovered, dt, params), mio::InfectionState::Susceptible);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(recovered, dt, params), recovered_state);
    }

    //setup location with 2 cells with some chance of exposure
    auto home      = mio::Location(mio::LocationType::Home, 0, 0);
    auto location2 = mio::Location(mio::LocationType::PublicTransport, 0, 2);
    auto infected4 =
        mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34, params, vaccination_state);
    home.add_person(infected4);
    infected4.migrate_to(home, location2, {0, 1});
    auto infected5 =
        mio::Person(home, mio::InfectionState::Infected, mio::AbmAgeGroup::Age80plus, params, vaccination_state);
    home.add_person(infected5);
    infected5.migrate_to(home, location2, {0});
    auto infected6 =
        mio::Person(home, mio::InfectionState::Infected, mio::AbmAgeGroup::Age5to14, params, vaccination_state);
    home.add_person(infected6);
    infected6.migrate_to(home, location2, {1});

    //cache precomputed results
    location2.begin_step(dt, params);

    {
        auto susceptible = mio::Person(home, mio::InfectionState::Susceptible, age, params, vaccination_state);
        home.add_person(susceptible);
        susceptible.migrate_to(home, location2, {0, 1});
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location2.interact(susceptible, dt, params), mio::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(2).WillOnce(Return(0.2)).WillOnce(Return(0.07));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location2.interact(susceptible, dt, params), mio::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(2).WillOnce(Return(0.15)).WillOnce(Return(0.15));
        EXPECT_EQ(location2.interact(susceptible, dt, params), mio::InfectionState::Susceptible);
    }
}

TEST(TestPerson, get_tested)
{
    using testing::Return;

    auto loc         = mio::Location(mio::LocationType::Home, 0);
    auto infected    = mio::Person(loc, mio::InfectionState::Infected, mio::AbmAgeGroup::Age15to34, {});
    auto susceptible = mio::Person(loc, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, {});

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(4)
        .WillOnce(Return(0.4))
        .WillOnce(Return(0.95))
        .WillOnce(Return(0.6))
        .WillOnce(Return(0.999));
    ASSERT_EQ(infected.get_tested({0.9, 0.99}), true);
    ASSERT_EQ(infected.get_tested({0.9, 0.99}), false);
    ASSERT_EQ(susceptible.get_tested({0.9, 0.99}), false);
    ASSERT_EQ(susceptible.get_tested({0.9, 0.99}), true);
    ASSERT_EQ(susceptible.get_time_since_negative_test(), mio::days(0));
}

TEST(TestPerson, getCells)
{
    auto home     = mio::Location(mio::LocationType::Home, 0, 0);
    auto location = mio::Location(mio::LocationType::PublicTransport, 0, 2);
    auto person   = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34, {});
    home.add_person(person);
    person.migrate_to(home, location, {0, 1});
    ASSERT_EQ(person.get_cells().size(), 2);
}

TEST(TestPerson, interact)
{
    using testing::Return;

    auto infection_parameters = mio::GlobalInfectionParameters();
    auto loc                  = mio::Location(mio::LocationType::Home, 0);
    auto person = mio::Person(loc, mio::InfectionState::Infected, mio::AbmAgeGroup::Age15to34, infection_parameters);
    loc.add_person(person);
    auto dt = mio::seconds(8640); //0.1 days
    loc.begin_step(dt, {});

    //setup rng mock so the person has a state transition
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    person.interact(dt, infection_parameters, loc, {});
    EXPECT_EQ(person.get_infection_state(), mio::InfectionState::Recovered_Infected);
    EXPECT_EQ(loc.get_subpopulation(mio::InfectionState::Recovered_Infected), 1);
    EXPECT_EQ(loc.get_subpopulation(mio::InfectionState::Infected), 0);
}

TEST(TestPerson, interact_exposed)
{
    using testing::Return;

    auto infection_parameters = mio::GlobalInfectionParameters();
    infection_parameters.set<mio::IncubationPeriod>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 2.});

    //setup location with some chance of exposure
    auto loc       = mio::Location(mio::LocationType::Work, 0);
    auto infected1 = mio::Person(loc, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34, infection_parameters);
    loc.add_person(infected1);
    auto infected2 = mio::Person(loc, mio::InfectionState::Infected, mio::AbmAgeGroup::Age5to14, infection_parameters);
    loc.add_person(infected2);
    auto infected3 = mio::Person(loc, mio::InfectionState::Infected, mio::AbmAgeGroup::Age60to79, infection_parameters);
    loc.add_person(infected3);
    auto person = mio::Person(loc, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, infection_parameters);
    loc.add_person(person);
    loc.begin_step(mio::hours(1), {});

    //setup rng mock so the person becomes exposed
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.49));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    //person becomes exposed
    person.interact(mio::hours(12), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Exposed);
    EXPECT_EQ(loc.get_subpopulation(mio::InfectionState::Exposed), 1);
    EXPECT_EQ(loc.get_subpopulation(mio::InfectionState::Carrier), 1);
    EXPECT_EQ(loc.get_subpopulation(mio::InfectionState::Infected), 2);

    //person becomes a carrier after the incubation time runs out, not random
    person.interact(mio::hours(12), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Exposed);

    person.interact(mio::hours(12), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Exposed);

    person.interact(mio::hours(24), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Exposed);

    person.interact(mio::hours(1), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Carrier);
    EXPECT_EQ(loc.get_subpopulation(mio::InfectionState::Exposed), 0);
    EXPECT_EQ(loc.get_subpopulation(mio::InfectionState::Carrier), 2);
    EXPECT_EQ(loc.get_subpopulation(mio::InfectionState::Infected), 2);
}

TEST(TestWorld, init)
{
    auto world = mio::World();
    for (uint32_t i = 0; i < (uint32_t)mio::LocationType::Count; i++) {
        ASSERT_THAT(world.get_locations()[i], testing::ElementsAre());
    }
    ASSERT_THAT(world.get_persons(), testing::ElementsAre());
}

TEST(TestWorld, addLocation)
{
    auto world      = mio::World();
    auto school_id1 = world.add_location(mio::LocationType::School);
    auto school_id2 = world.add_location(mio::LocationType::School);
    auto work_id    = world.add_location(mio::LocationType::Work);
    auto home_id    = world.add_location(mio::LocationType::Home);

    ASSERT_EQ((int)school_id1.index, 0);
    ASSERT_EQ((int)school_id2.index, 1);

    auto& school1 = world.get_individualized_location(school_id1);
    auto& school2 = world.get_individualized_location(school_id2);
    auto& work    = world.get_individualized_location(work_id);
    auto& home    = world.get_individualized_location(home_id);

    ASSERT_EQ(world.get_locations().size(), (uint32_t)mio::LocationType::Count);
    ASSERT_EQ(world.get_locations()[(uint32_t)mio::LocationType::School].size(), 2);

    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::LocationType::School][0], &school1);
    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::LocationType::School][1], &school2);
    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::LocationType::Work][0], &work);
    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::LocationType::Home][0], &home);
}

TEST(TestWorld, addPerson)
{
    auto world    = mio::World();
    auto location = world.add_location(mio::LocationType::School);

    auto& p1 = world.add_person(location, mio::InfectionState::Recovered_Carrier);
    auto& p2 = world.add_person(location, mio::InfectionState::Exposed);

    ASSERT_EQ(world.get_persons().size(), 2);
    ASSERT_EQ(&world.get_persons()[0], &p1);
    ASSERT_EQ(&world.get_persons()[1], &p2);
}

TEST(TestWorld, getSubpopulationCombined)
{
    auto world   = mio::World();
    auto school1 = world.add_location(mio::LocationType::School);
    auto school2 = world.add_location(mio::LocationType::School);
    auto school3 = world.add_location(mio::LocationType::School);
    world.add_person(school1, mio::InfectionState::Carrier);
    world.add_person(school1, mio::InfectionState::Susceptible);
    world.add_person(school2, mio::InfectionState::Susceptible);
    world.add_person(school2, mio::InfectionState::Susceptible);
    world.add_person(school3, mio::InfectionState::Carrier);

    ASSERT_EQ(world.get_subpopulation_combined(mio::InfectionState::Susceptible, mio::LocationType::School), 3);
    ASSERT_EQ(world.get_subpopulation_combined(mio::InfectionState::Carrier, mio::LocationType::School), 2);
}

TEST(TestWorld, evolveStateTransition)
{
    using testing::Return;

    auto world     = mio::World();
    auto location1 = world.add_location(mio::LocationType::School);
    auto& p1       = world.add_person(location1, mio::InfectionState::Carrier);
    auto& p2       = world.add_person(location1, mio::InfectionState::Susceptible);
    auto location2 = world.add_location(mio::LocationType::Work);
    auto& p3       = world.add_person(location2, mio::InfectionState::Infected);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);

    //setup mock so only p2 transitions
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::AtLeast(3))
        .WillOnce(Return(0.51))
        .WillOnce(Return(0.04))
        .WillOnce(Return(0.6))
        .WillRepeatedly(Return(1.0));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    world.evolve(mio::TimePoint(0), mio::hours(1));

    EXPECT_EQ(p1.get_infection_state(), mio::InfectionState::Carrier);
    EXPECT_EQ(p2.get_infection_state(), mio::InfectionState::Exposed);
    EXPECT_EQ(p3.get_infection_state(), mio::InfectionState::Infected);
}

TEST(TestMigrationRules, quarantine)
{
    auto home   = mio::Location(mio::LocationType::Home, 0);
    auto work   = mio::Location(mio::LocationType::Work, 0);
    auto p_inf1 = mio::Person(work, {mio::InfectionState::Infected, true}, mio::AbmAgeGroup::Age15to34, {});
    auto t      = mio::TimePoint(12346);
    auto dt     = mio::hours(1);
    p_inf1.set_assigned_location(home);

    ASSERT_EQ(mio::go_to_quarantine(p_inf1, t, dt, {}), mio::LocationType::Home);

    auto p_inf2 = mio::Person(work, mio::InfectionState::Infected, mio::AbmAgeGroup::Age15to34, {});
    p_inf2.set_assigned_location(home);
    ASSERT_EQ(mio::go_to_quarantine(p_inf2, t, dt, {}), mio::LocationType::Work);
}

TEST(TestMigrationRules, hospital)
{
    auto home  = mio::Location(mio::LocationType::Home, 0);
    auto p_inf = mio::Person(home, mio::InfectionState::Infected_Severe, mio::AbmAgeGroup::Age15to34, {});
    auto t     = mio::TimePoint(12346);
    auto dt    = mio::hours(1);

    ASSERT_EQ(mio::go_to_hospital(p_inf, t, dt, {}), mio::LocationType::Hospital);

    auto p_car = mio::Person(home, mio::InfectionState::Infected, mio::AbmAgeGroup::Age15to34, {});
    ASSERT_EQ(mio::go_to_hospital(p_car, t, dt, {}), mio::LocationType::Home);
}

TEST(TestMigrationRules, icu)
{
    auto hospital = mio::Location(mio::LocationType::Hospital, 0);
    auto p_hosp   = mio::Person(hospital, mio::InfectionState::Infected_Critical, mio::AbmAgeGroup::Age15to34, {});
    auto t        = mio::TimePoint(12346);
    auto dt       = mio::hours(1);

    ASSERT_EQ(mio::go_to_icu(p_hosp, t, dt, {}), mio::LocationType::ICU);

    auto work   = mio::Location(mio::LocationType::Work, 0);
    auto p_work = mio::Person(work, mio::InfectionState::Infected, mio::AbmAgeGroup::Age15to34, {});
    ASSERT_EQ(mio::go_to_icu(p_work, t, dt, {}), mio::LocationType::Work);
}

TEST(TestMigrationRules, recover)
{
    auto hospital = mio::Location(mio::LocationType::Hospital, 0);
    auto p_rec    = mio::Person(hospital, mio::InfectionState::Recovered_Infected, mio::AbmAgeGroup::Age60to79, {});
    auto p_inf    = mio::Person(hospital, mio::InfectionState::Infected_Severe, mio::AbmAgeGroup::Age60to79, {});
    auto t        = mio::TimePoint(12346);
    auto dt       = mio::hours(1);

    ASSERT_EQ(mio::return_home_when_recovered(p_rec, t, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::return_home_when_recovered(p_inf, t, dt, {}), mio::LocationType::Hospital);
}

TEST(TestTestingScheme, init)
{
    auto tests = mio::TestingScheme(mio::days(7), 0.8);
    ASSERT_EQ(tests.get_interval(), mio::days(7));
    ASSERT_EQ(tests.get_probability(), 0.8);

    tests.set_interval(mio::days(2));
    ASSERT_EQ(tests.get_interval(), mio::days(2));
}

TEST(TestTestingScheme, runScheme)
{
    auto loc     = mio::Location(mio::LocationType::Home, 0);
    auto person1 = mio::Person(loc, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age5to14, {});
    auto person2 = mio::Person(loc, mio::InfectionState::Recovered_Carrier, mio::AbmAgeGroup::Age5to14, {});
    auto testing = mio::TestingScheme(mio::days(5), 0.9);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.5));
    testing.run_scheme(person1, {});

    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.5));
    testing.run_scheme(person2, {});

    ASSERT_EQ(person2.get_time_since_negative_test(), mio::days(0));
    ASSERT_EQ(person1.is_in_quarantine(), true);
    ASSERT_EQ(person2.is_in_quarantine(), false);
}

TEST(TestWorld, evolveMigration)
{
    using testing::Return;

    auto world       = mio::World();
    auto home_id     = world.add_location(mio::LocationType::Home);
    auto school_id   = world.add_location(mio::LocationType::School);
    auto work_id     = world.add_location(mio::LocationType::Work);
    auto hospital_id = world.add_location(mio::LocationType::Hospital);

    auto& p1 = world.add_person(home_id, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34);
    auto& p2 = world.add_person(home_id, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14);
    auto& p3 = world.add_person(home_id, mio::InfectionState::Infected_Severe, mio::AbmAgeGroup::Age5to14);
    auto& p4 = world.add_person(hospital_id, mio::InfectionState::Recovered_Infected, mio::AbmAgeGroup::Age5to14);
    auto& p5 = world.add_person(home_id, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34);
    p1.set_assigned_location(school_id);
    p2.set_assigned_location(school_id);
    p1.set_assigned_location(work_id);
    p2.set_assigned_location(work_id);
    p1.set_assigned_location(home_id);
    p2.set_assigned_location(home_id);
    p3.set_assigned_location(home_id);
    p4.set_assigned_location(home_id);
    p3.set_assigned_location(hospital_id);
    p5.set_assigned_location(school_id);
    p5.set_assigned_location(work_id);
    p5.set_assigned_location(home_id);

    mio::TripList& data = world.get_trip_list();
    mio::Trip trip1(p1.get_person_id(), mio::TimePoint(0) + mio::hours(9), work_id, home_id);
    mio::Trip trip2(p2.get_person_id(), mio::TimePoint(1) + mio::hours(9), school_id, home_id);
    mio::Trip trip3(p5.get_person_id(), mio::TimePoint(2) + mio::hours(9), school_id, work_id);
    data.add_trip(trip1);
    data.add_trip(trip2);
    data.add_trip(trip3);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

    world.evolve(mio::TimePoint(0) + mio::hours(8), mio::hours(2));

    auto& school   = world.get_individualized_location(school_id);
    auto& work     = world.get_individualized_location(work_id);
    auto& home     = world.get_individualized_location(home_id);
    auto& hospital = world.get_individualized_location(hospital_id);

    EXPECT_EQ(p1.get_location_id().type, mio::LocationType::Work);
    EXPECT_EQ(p2.get_location_id().type, mio::LocationType::School);
    EXPECT_EQ(p3.get_location_id().type, mio::LocationType::Hospital);
    EXPECT_EQ(p4.get_location_id().type, mio::LocationType::Home);
    EXPECT_EQ(p5.get_location_id().type, mio::LocationType::Home);
    EXPECT_EQ(school.get_subpopulations().sum(), 1);
    EXPECT_EQ(work.get_subpopulations().sum(), 1);
    EXPECT_EQ(home.get_subpopulations().sum(), 2);
    EXPECT_EQ(hospital.get_subpopulations().sum(), 1);
}

TEST(TestSimulation, advance_random)
{
    auto world     = mio::World();
    auto location1 = world.add_location(mio::LocationType::School);
    auto location2 = world.add_location(mio::LocationType::School);
    auto& p1       = world.add_person(location1, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age5to14);
    auto& p2       = world.add_person(location1, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14);
    auto& p3       = world.add_person(location2, mio::InfectionState::Infected, mio::AbmAgeGroup::Age5to14);
    auto& p4       = world.add_person(location2, mio::InfectionState::Infected, mio::AbmAgeGroup::Age5to14);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);
    p4.set_assigned_location(location2);

    auto sim = mio::AbmSimulation(mio::TimePoint(0), std::move(world));

    sim.advance(mio::TimePoint(0) + mio::hours(50));
    ASSERT_EQ(sim.get_result().get_num_time_points(), 51);
    ASSERT_THAT(sim.get_result().get_times(), ElementsAreLinspace(0.0, 50.0 / 24.0, 51));
    for (auto&& v : sim.get_result()) {
        ASSERT_EQ(v.sum(), 4);
    }
}

TEST(TestDiscreteDistribution, generate)
{
    using namespace mio;
    auto distribution = mio::DiscreteDistribution<size_t>();

    std::vector<double> weights;
    for (size_t i = 0; i < 50; i++) {
        weights = {};
        ASSERT_EQ(distribution(weights), 0);

        weights = {0.5};
        ASSERT_EQ(distribution(weights), 0);

        weights = {0.5, 1.3, 0.1, 0.4, 0.3};
        auto d  = distribution(weights);
        ASSERT_GE(d, 0);
        ASSERT_LE(d, 4);
    }
}
