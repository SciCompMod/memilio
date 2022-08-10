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
#include "abm/age.h"
#include "abm/location_type.h"
#include "abm/migration_rules.h"
#include "abm/lockdown_rules.h"
#include "memilio/math/eigen_util.h"
#include "matchers.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <memory>

TEST(TestLocation, init)
{
    auto location = mio::abm::Location(mio::abm::LocationType::School, 0);
    for (mio::abm::InfectionState i = mio::abm::InfectionState(0); i < mio::abm::InfectionState::Count;
         i                          = mio::abm::InfectionState(size_t(i) + 1)) {
        ASSERT_EQ(location.get_subpopulation(i), 0);
    }
    ASSERT_EQ(print_wrap(location.get_subpopulations()),
              print_wrap(Eigen::VectorXi::Zero(Eigen::Index(mio::abm::InfectionState::Count))));
}

TEST(TestLocation, initCell)
{
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 2);
    ASSERT_EQ(location.get_cells().size(), 2);
}

TEST(TestLocation, GetIndex)
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

    auto location = mio::abm::Location(mio::abm::LocationType::Work, 0);
    //setup rng mock so the time_until_carrier is 1.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformIntDistribution<int>>>>
        mock_uniform_int_dist;
    EXPECT_CALL(mock_uniform_int_dist.get_mock(), invoke).Times(testing::AtLeast(1)).WillRepeatedly(Return(0));

    auto person = mio::abm::Person(location, mio::abm::InfectionState::Exposed, mio::abm::AgeGroup::Age60to79, {});
    location.add_person(person);

    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);
    ASSERT_EQ(person.get_location_id().index, location.get_index());
    ASSERT_EQ(person.get_location_id().type, location.get_type());
    ASSERT_EQ(person.get_person_id(), mio::abm::INVALID_PERSON_ID);

    auto person2 = mio::abm::Person(location, mio::abm::InfectionState::Exposed, mio::abm::AgeGroup::Age60to79, {},
                                    mio::abm::VaccinationState::Unvaccinated, 0);
    ASSERT_EQ(person2.get_person_id(), 0u);

    mio::abm::TimeSpan dt = mio::abm::hours(1);
    person.interact(dt, {}, location);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Carrier);
}

TEST(TestPerson, migrate)
{
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0, 0);
    auto loc1   = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 1);
    auto loc2   = mio::abm::Location(mio::abm::LocationType::School, 0);
    auto loc3   = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 2);
    auto person = mio::abm::Person(home, mio::abm::InfectionState::Recovered_Carrier, mio::abm::AgeGroup::Age0to4, {});
    home.add_person(person);
    person.migrate_to(home, loc1, {0});

    ASSERT_EQ(person.get_location_id().index, loc1.get_index());
    ASSERT_EQ(person.get_location_id().type, loc1.get_type());
    ASSERT_EQ(loc1.get_subpopulation(mio::abm::InfectionState::Recovered_Carrier), 1);
    ASSERT_EQ(home.get_subpopulation(mio::abm::InfectionState::Recovered_Carrier), 0);
    ASSERT_EQ(loc1.get_cells()[0].num_people, 1u);

    person.migrate_to(loc1, loc2);

    ASSERT_EQ(person.get_location_id().index, loc2.get_index());
    ASSERT_EQ(person.get_location_id().type, loc2.get_type());
    ASSERT_EQ(loc2.get_subpopulation(mio::abm::InfectionState::Recovered_Carrier), 1);
    ASSERT_EQ(loc1.get_subpopulation(mio::abm::InfectionState::Recovered_Carrier), 0);
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
    auto location = mio::abm::Location(mio::abm::LocationType::Work, 2);
    auto person =
        mio::abm::Person(location, mio::abm::InfectionState::Recovered_Carrier, mio::abm::AgeGroup::Age60to79, {});
    person.set_assigned_location(location);
    ASSERT_EQ((int)person.get_assigned_location_index(mio::abm::LocationType::Work), 2);

    person.set_assigned_location({4, mio::abm::LocationType::Work});
    ASSERT_EQ((int)person.get_assigned_location_index(mio::abm::LocationType::Work), 4);
}

TEST(TestWorld, findLocation)
{
    auto world     = mio::abm::World();
    auto home_id   = world.add_location(mio::abm::LocationType::Home);
    auto school_id = world.add_location(mio::abm::LocationType::School);
    auto work_id   = world.add_location(mio::abm::LocationType::Work);
    auto person  = mio::abm::Person(home_id, mio::abm::InfectionState::Recovered_Carrier, mio::abm::AgeGroup::Age60to79,
                                    world.get_global_infection_parameters());
    auto& home   = world.get_individualized_location(home_id);
    auto& school = world.get_individualized_location(school_id);
    auto& work   = world.get_individualized_location(work_id);
    person.set_assigned_location(home);
    person.set_assigned_location(school);
    person.set_assigned_location({0, mio::abm::LocationType::Work});

    ASSERT_EQ(world.find_location(mio::abm::LocationType::Work, person), &work);
    ASSERT_EQ(world.find_location(mio::abm::LocationType::School, person), &school);
    ASSERT_EQ(world.find_location(mio::abm::LocationType::Home, person), &home);
}

TEST(TestLocation, beginStep)
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
    params.get<mio::abm::SusceptibleToExposedByCarrier>()[{age, vaccination_state}] = 0.4;
    params.set<mio::abm::SusceptibleToExposedByInfected>(
        {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 0.});
    params.get<mio::abm::SusceptibleToExposedByInfected>()[{age, vaccination_state}] = 0.5;

    //setup location with some chance of exposure
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0, 0);
    auto location1 = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 3);
    auto infected1 = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, params,
                                      vaccination_state);
    home.add_person(infected1);
    infected1.migrate_to(home, location1, {0});
    auto infected2 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age80plus, params,
                                      vaccination_state);
    home.add_person(infected2);
    infected2.migrate_to(home, location1, {0, 1});
    auto infected3 = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14, params,
                                      vaccination_state);
    home.add_person(infected3);
    infected3.migrate_to(home, location1, {1});

    //cache precomputed results
    auto dt = mio::abm::seconds(8640);
    location1.begin_step(dt, params);

    ASSERT_TRUE(std::abs(location1.get_cells()[0].cached_exposure_rate[{age, vaccination_state}] - 0.9) < 0.001);
    ASSERT_TRUE(std::abs(location1.get_cells()[1].cached_exposure_rate[{age, vaccination_state}] - 1) < 0.001);
    ASSERT_TRUE(std::abs(location1.get_cells()[2].cached_exposure_rate[{age, vaccination_state}]) < 0.001);
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

TEST(TestPerson, quarantine)
{
    using testing::Return;

    auto infection_parameters = mio::abm::GlobalInfectionParameters();
    auto home                 = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto work                 = mio::abm::Location(mio::abm::LocationType::Work, 0);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));

    auto person = mio::abm::Person(home, mio::abm::InfectionProperties(mio::abm::InfectionState::Infected, true),
                                   mio::abm::AgeGroup::Age15to34, infection_parameters);
    home.add_person(person);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(7);
    auto dt        = mio::abm::hours(1);
    ASSERT_EQ(mio::abm::go_to_work(person, t_morning, dt, {}), mio::abm::LocationType::Home);
    //setup rng mock so the person has a state transition to Recovered_Infected
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.04));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
    person.interact(dt, infection_parameters, home);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Recovered_Infected);
    ASSERT_EQ(mio::abm::go_to_work(person, t_morning, dt, {}), mio::abm::LocationType::Work);
}

TEST(TestPerson, get_tested)
{
    using testing::Return;

    auto loc         = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto infected    = mio::abm::Person(loc, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age15to34, {});
    auto susceptible = mio::abm::Person(loc, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, {});

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
    ASSERT_EQ(susceptible.get_time_since_negative_test(), mio::abm::days(0));
}

TEST(TestPerson, getCells)
{
    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0, 0);
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 2);
    auto person   = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, {});
    home.add_person(person);
    person.migrate_to(home, location, {0, 1});
    ASSERT_EQ(person.get_cells().size(), 2);
}

TEST(TestPerson, interact)
{
    using testing::Return;

    auto infection_parameters = mio::abm::GlobalInfectionParameters();
    auto loc                  = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto person =
        mio::abm::Person(loc, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age15to34, infection_parameters);
    loc.add_person(person);
    auto dt = mio::abm::seconds(8640); //0.1 days
    loc.begin_step(dt, {});

    //setup rng mock so the person has a state transition
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    person.interact(dt, infection_parameters, loc);
    EXPECT_EQ(person.get_infection_state(), mio::abm::InfectionState::Recovered_Infected);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Recovered_Infected), 1);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Infected), 0);
}

TEST(TestPerson, interact_exposed)
{
    using testing::Return;

    auto infection_parameters = mio::abm::GlobalInfectionParameters();
    infection_parameters.set<mio::abm::IncubationPeriod>(
        {{mio::abm::AgeGroup::Count, mio::abm::VaccinationState::Count}, 2.});

    //setup location with some chance of exposure
    auto loc = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto infected1 =
        mio::abm::Person(loc, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, infection_parameters);
    loc.add_person(infected1);
    auto infected2 =
        mio::abm::Person(loc, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14, infection_parameters);
    loc.add_person(infected2);
    auto infected3 =
        mio::abm::Person(loc, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age60to79, infection_parameters);
    loc.add_person(infected3);
    auto person = mio::abm::Person(loc, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34,
                                   infection_parameters);
    loc.add_person(person);
    loc.begin_step(mio::abm::hours(1), {});

    //setup rng mock so the person becomes exposed
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.49));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    //person becomes exposed
    person.interact(mio::abm::hours(12), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Exposed), 1);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Carrier), 1);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Infected), 2);

    //person becomes a carrier after the incubation time runs out, not random
    person.interact(mio::abm::hours(12), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);

    person.interact(mio::abm::hours(12), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);

    person.interact(mio::abm::hours(24), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);

    person.interact(mio::abm::hours(1), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Carrier);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Exposed), 0);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Carrier), 2);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Infected), 2);
}

TEST(TestWorld, init)
{
    auto world = mio::abm::World();
    for (uint32_t i = 0; i < (uint32_t)mio::abm::LocationType::Count; i++) {
        ASSERT_THAT(world.get_locations()[i], testing::ElementsAre());
    }
    ASSERT_THAT(world.get_persons(), testing::ElementsAre());
}

TEST(TestWorld, addLocation)
{
    auto world      = mio::abm::World();
    auto school_id1 = world.add_location(mio::abm::LocationType::School);
    auto school_id2 = world.add_location(mio::abm::LocationType::School);
    auto work_id    = world.add_location(mio::abm::LocationType::Work);
    auto home_id    = world.add_location(mio::abm::LocationType::Home);

    ASSERT_EQ((int)school_id1.index, 0);
    ASSERT_EQ((int)school_id2.index, 1);

    auto& school1 = world.get_individualized_location(school_id1);
    auto& school2 = world.get_individualized_location(school_id2);
    auto& work    = world.get_individualized_location(work_id);
    auto& home    = world.get_individualized_location(home_id);

    ASSERT_EQ(world.get_locations().size(), (uint32_t)mio::abm::LocationType::Count);
    ASSERT_EQ(world.get_locations()[(uint32_t)mio::abm::LocationType::School].size(), 2);

    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::abm::LocationType::School][0], &school1);
    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::abm::LocationType::School][1], &school2);
    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::abm::LocationType::Work][0], &work);
    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::abm::LocationType::Home][0], &home);
}

TEST(TestWorld, addPerson)
{
    auto world    = mio::abm::World();
    auto location = world.add_location(mio::abm::LocationType::School);

    auto& p1 = world.add_person(location, mio::abm::InfectionState::Recovered_Carrier);
    auto& p2 = world.add_person(location, mio::abm::InfectionState::Exposed);

    ASSERT_EQ(world.get_persons().size(), 2);
    ASSERT_EQ(&world.get_persons()[0], &p1);
    ASSERT_EQ(&world.get_persons()[1], &p2);
}

TEST(TestWorld, getSubpopulationCombined)
{
    auto world   = mio::abm::World();
    auto school1 = world.add_location(mio::abm::LocationType::School);
    auto school2 = world.add_location(mio::abm::LocationType::School);
    auto school3 = world.add_location(mio::abm::LocationType::School);
    world.add_person(school1, mio::abm::InfectionState::Carrier);
    world.add_person(school1, mio::abm::InfectionState::Susceptible);
    world.add_person(school2, mio::abm::InfectionState::Susceptible);
    world.add_person(school2, mio::abm::InfectionState::Susceptible);
    world.add_person(school3, mio::abm::InfectionState::Carrier);

    ASSERT_EQ(world.get_subpopulation_combined(mio::abm::InfectionState::Susceptible, mio::abm::LocationType::School),
              3);
    ASSERT_EQ(world.get_subpopulation_combined(mio::abm::InfectionState::Carrier, mio::abm::LocationType::School), 2);
}

TEST(TestWorld, evolveStateTransition)
{
    using testing::Return;

    auto world     = mio::abm::World();
    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto& p1       = world.add_person(location1, mio::abm::InfectionState::Carrier);
    auto& p2       = world.add_person(location1, mio::abm::InfectionState::Susceptible);
    auto location2 = world.add_location(mio::abm::LocationType::Work);
    auto& p3       = world.add_person(location2, mio::abm::InfectionState::Infected);
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

    world.evolve(mio::abm::TimePoint(0), mio::abm::hours(1));

    EXPECT_EQ(p1.get_infection_state(), mio::abm::InfectionState::Carrier);
    EXPECT_EQ(p2.get_infection_state(), mio::abm::InfectionState::Exposed);
    EXPECT_EQ(p3.get_infection_state(), mio::abm::InfectionState::Infected);
}

TEST(TestMigrationRules, student_goes_to_school)
{
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));

    auto home    = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_child = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14, {});
    auto p_adult = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, {});

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(7);
    auto t_weekend = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(7);
    auto dt        = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_school(p_child, t_morning, dt, {}), mio::abm::LocationType::School);
    ASSERT_EQ(mio::abm::go_to_school(p_adult, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child, t_weekend, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, students_go_to_school_in_different_times)
{
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        //Mocking the random values will define at what time the student should go to school, i.e:
        // random is in [0,1/3] -> goes to school at 6
        // random is in [1/3,2/3] -> goes to school at 7
        // random is in [2/3,1.] -> goes to school at 8
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.8))
        .WillOnce(testing::Return(0.8))
        .WillOnce(testing::Return(0.8))
        .WillOnce(testing::Return(0.8))
        .WillRepeatedly(testing::Return(1.0));

    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_child_goes_to_school_at_6 =
        mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14, {});
    auto p_child_goes_to_school_at_8 =
        mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14, {});

    auto t_morning_6 = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8 = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt          = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_6, t_morning_6, dt, {}), mio::abm::LocationType::School);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_6, t_morning_8, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_8, t_morning_6, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_8, t_morning_8, dt, {}), mio::abm::LocationType::School);
}

TEST(TestMigrationRules, students_go_to_school_in_different_times_with_smaller_time_steps)
{
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        //Mocking the random values will define at what time the student should go to school, i.e:
        // random is in [0,1/6] -> goes to school at 6
        // random is in [1/6,2/6] -> goes to school at 6:30
        // random is in [2/6,3/6] -> goes to school at 7:00
        // random is in [3/6,4/6] -> goes to school at 7:30
        // random is in [4/6,5/6] -> goes to school at 8:00
        // random is in [5/6,6/6] -> goes to school at 8:30
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillRepeatedly(testing::Return(1.0));

    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_child_goes_to_school_at_6 =
        mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14, {});
    auto p_child_goes_to_school_at_8_30 =
        mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14, {});

    auto t_morning_6    = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8_30 = mio::abm::TimePoint(0) + mio::abm::hours(8) + mio::abm::seconds(1800);
    auto dt             = mio::abm::seconds(1800);

    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_6, t_morning_6, dt, {}), mio::abm::LocationType::School);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_6, t_morning_8_30, dt, {}),
              mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_8_30, t_morning_6, dt, {}),
              mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_8_30, t_morning_8_30, dt, {}),
              mio::abm::LocationType::School);
}

TEST(TestMigrationRules, school_return)
{
    auto school  = mio::abm::Location(mio::abm::LocationType::School, 0);
    auto p_child = mio::abm::Person(school, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14, {});

    auto t  = mio::abm::TimePoint(0) + mio::abm::hours(15);
    auto dt = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_school(p_child, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work)
{
    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillRepeatedly(testing::Return(1.0));

    auto p_retiree = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age60to79, {});
    auto p_adult   = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, {});

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt        = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_work(p_retiree, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t_night, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work_with_non_dividable_timespan)
{
    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillRepeatedly(testing::Return(1.0));

    auto p_retiree = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age60to79, {});
    auto p_adult   = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, {});

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt        = mio::abm::minutes(53);

    ASSERT_EQ(mio::abm::go_to_work(p_retiree, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t_night, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, workers_go_to_work_in_different_times)
{
    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))

        .WillRepeatedly(testing::Return(1.0));

    auto p_adult_goes_to_work_at_6 =
        mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, {});
    auto p_adult_goes_to_work_at_8 =
        mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, {});

    auto t_morning_6 = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8 = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night     = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt          = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_6, t_morning_6, dt, {}), mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_6, t_morning_8, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_6, t_night, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_8, t_morning_6, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_8, t_morning_8, dt, {}), mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_8, t_night, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, work_return)
{
    auto work    = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto p_adult = mio::abm::Person(work, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age35to59, {});
    auto t       = mio::abm::TimePoint(0) + mio::abm::hours(17);
    auto dt      = mio::abm::hours(1);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, quarantine)
{
    auto t  = mio::abm::TimePoint(12346);
    auto dt = mio::abm::hours(1);

    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto work     = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto hospital = mio::abm::Location(mio::abm::LocationType::Hospital, 0);

    auto p_inf1 = mio::abm::Person(work, {mio::abm::InfectionState::Infected, true}, mio::abm::AgeGroup::Age15to34, {});
    ASSERT_EQ(mio::abm::go_to_quarantine(p_inf1, t, dt, {}),
              mio::abm::LocationType::Home); //detected infected person quarantines at home

    auto p_inf2 =
        mio::abm::Person(work, {mio::abm::InfectionState::Infected, false}, mio::abm::AgeGroup::Age15to34, {});
    ASSERT_EQ(mio::abm::go_to_quarantine(p_inf2, t, dt, {}),
              mio::abm::LocationType::Work); //undetected infected person does not quaratine

    auto p_inf3 = mio::abm::Person(hospital, {mio::abm::InfectionState::Infected_Severe, true},
                                   mio::abm::AgeGroup::Age15to34, {});
    ASSERT_EQ(mio::abm::go_to_quarantine(p_inf3, t, dt, {}),
              mio::abm::LocationType::Hospital); //detected infected person does not leave hospital to quarantine
}

TEST(TestMigrationRules, hospital)
{
    auto home  = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_inf = mio::abm::Person(home, mio::abm::InfectionState::Infected_Severe, mio::abm::AgeGroup::Age15to34, {});
    auto t     = mio::abm::TimePoint(12346);
    auto dt    = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_hospital(p_inf, t, dt, {}), mio::abm::LocationType::Hospital);

    auto p_car = mio::abm::Person(home, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age15to34, {});
    ASSERT_EQ(mio::abm::go_to_hospital(p_car, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, go_shopping)
{
    auto hospital = mio::abm::Location(mio::abm::LocationType::Hospital, 0);
    auto p_hosp   = mio::abm::Person(hospital, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age0to4, {});
    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_home   = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age60to79, {});

    auto t_weekday = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto t_sunday  = mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(9);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(1);
    auto dt        = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_shop(p_hosp, t_weekday, dt, {}), mio::abm::LocationType::Hospital);
    ASSERT_EQ(mio::abm::go_to_shop(p_home, t_sunday, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_shop(p_home, t_night, dt, {}), mio::abm::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_shop(p_home, t_weekday, dt, {}), mio::abm::LocationType::BasicsShop);
}

TEST(TestMigrationRules, shop_return)
{
    auto t  = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto dt = mio::abm::hours(1);

    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto shop = mio::abm::Location(mio::abm::LocationType::BasicsShop, 0);
    auto p    = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, {});
    home.add_person(p);
    p.migrate_to(home, shop);
    p.interact(dt, {}, shop); //person only returns home after some time passed

    ASSERT_EQ(mio::abm::go_to_shop(p, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, go_event)
{
    auto work   = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto p_work = mio::abm::Person(work, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age35to59, {});
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_home = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age60to79, {});

    auto t_weekday  = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(20);
    auto t_saturday = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(10);
    auto t_night    = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(1);
    auto dt         = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_event(p_work, t_weekday, dt, {}), mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_event(p_home, t_night, dt, {}), mio::abm::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_event(p_home, t_weekday, dt, {}), mio::abm::LocationType::SocialEvent);

    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_event(p_home, t_saturday, dt, {}), mio::abm::LocationType::SocialEvent);
}

TEST(TestMigrationRules, event_return)
{
    auto t  = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(21);
    auto dt = mio::abm::hours(3);

    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto shop = mio::abm::Location(mio::abm::LocationType::SocialEvent, 0);
    auto p    = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, {});
    home.add_person(p);
    p.migrate_to(home, shop);
    p.interact(dt, {}, shop);

    ASSERT_EQ(mio::abm::go_to_event(p, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestLockdownRules, school_closure)
{
    auto t         = mio::abm::TimePoint(0);
    auto dt        = mio::abm::hours(1);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto school    = mio::abm::Location(mio::abm::LocationType::School, 0);

    //setup rng mock so one person is home schooled and the other goes to school
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.2))
        .WillOnce(testing::Return(0.2))
        .WillOnce(testing::Return(0.2))
        .WillOnce(testing::Return(0.2))
        .WillRepeatedly(testing::Return(1.0));

    auto p1 = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age5to14, {});
    p1.set_assigned_location(home);
    p1.set_assigned_location(school);
    auto p2 = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age5to14, {});
    p2.set_assigned_location(home);
    p2.set_assigned_location(school);
    mio::abm::MigrationParameters params;

    mio::abm::set_school_closure(t, 0.7, params);

    ASSERT_EQ(mio::abm::go_to_school(p1, t_morning, dt, params), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p2, t_morning, dt, params), mio::abm::LocationType::School);
}

TEST(TestLockdownRules, school_opening)
{
    auto t_closing = mio::abm::TimePoint(0);
    auto t_opening = mio::abm::TimePoint(0) + mio::abm::days(1);
    auto dt        = mio::abm::hours(1);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(7);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto school    = mio::abm::Location(mio::abm::LocationType::School, 0);
    //setup rng mock so the person is homeschooled in case of lockdown
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));
    auto p = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age5to14, {});
    p.set_assigned_location(home);
    p.set_assigned_location(school);
    mio::abm::MigrationParameters params;

    mio::abm::set_school_closure(t_closing, 1., params);
    mio::abm::set_school_closure(t_opening, 0., params);

    ASSERT_EQ(mio::abm::go_to_school(p, t_morning, dt, params), mio::abm::LocationType::School);
}

TEST(TestLockdownRules, home_office)
{
    auto t         = mio::abm::TimePoint(0);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt        = mio::abm::hours(1);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto work      = mio::abm::Location(mio::abm::LocationType::Work, 0);
    mio::abm::MigrationParameters params;

    mio::abm::set_home_office(t, 0.4, params);

    //setup rng mock so one person goes to work and the other works at home
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(4))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillRepeatedly(testing::Return(1.0));

    auto person1 = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, {});
    auto person2 = mio::abm::Person(home, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34, {});
    person1.set_assigned_location(home);
    person1.set_assigned_location(work);
    person2.set_assigned_location(home);
    person2.set_assigned_location(work);

    ASSERT_EQ(mio::abm::go_to_work(person1, t_morning, dt, params), mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_work(person2, t_morning, dt, params), mio::abm::LocationType::Home);
}

TEST(TestLockdownRules, no_home_office)
{
    auto t_closing = mio::abm::TimePoint(0);
    auto t_opening = mio::abm::TimePoint(0) + mio::abm::days(1);
    auto dt        = mio::abm::hours(1);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(8);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto work      = mio::abm::Location(mio::abm::LocationType::Work, 0);

    //setup rng mock so the person works in home office
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillRepeatedly(testing::Return(1.0));

    auto p = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, {});
    p.set_assigned_location(home);
    p.set_assigned_location(work);
    mio::abm::MigrationParameters params;

    mio::abm::set_home_office(t_closing, 0.5, params);
    mio::abm::set_home_office(t_opening, 0., params);

    ASSERT_EQ(mio::abm::go_to_work(p, t_morning, dt, params), mio::abm::LocationType::Work);
}

TEST(TestLockdownRules, social_event_closure)
{
    auto t         = mio::abm::TimePoint(0);
    auto dt        = mio::abm::hours(1);
    auto t_evening = mio::abm::TimePoint(0) + mio::abm::hours(19);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto event     = mio::abm::Location(mio::abm::LocationType::SocialEvent, 0);
    auto p         = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age5to14, {});
    p.set_assigned_location(home);
    p.set_assigned_location(event);
    mio::abm::MigrationParameters params;

    mio::abm::close_social_events(t, 1, params);

    ASSERT_EQ(mio::abm::go_to_event(p, t_evening, dt, params), mio::abm::LocationType::Home);
}

TEST(TestLockdownRules, social_events_opening)
{
    auto t_closing = mio::abm::TimePoint(0);
    auto t_opening = mio::abm::TimePoint(0) + mio::abm::days(1);
    auto dt        = mio::abm::hours(1);
    auto t_evening = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(19);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto event     = mio::abm::Location(mio::abm::LocationType::SocialEvent, 0);
    auto p         = mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age5to14, {});
    p.set_assigned_location(event);
    p.set_assigned_location(home);
    mio::abm::MigrationParameters params;

    mio::abm::close_social_events(t_closing, 1, params);
    mio::abm::close_social_events(t_opening, 0, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_event(p, t_evening, dt, params), mio::abm::LocationType::SocialEvent);
}

TEST(TestMigrationRules, icu)
{
    auto hospital = mio::abm::Location(mio::abm::LocationType::Hospital, 0);
    auto p_hosp =
        mio::abm::Person(hospital, mio::abm::InfectionState::Infected_Critical, mio::abm::AgeGroup::Age15to34, {});
    auto t  = mio::abm::TimePoint(12346);
    auto dt = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_icu(p_hosp, t, dt, {}), mio::abm::LocationType::ICU);

    auto work   = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto p_work = mio::abm::Person(work, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age15to34, {});
    ASSERT_EQ(mio::abm::go_to_icu(p_work, t, dt, {}), mio::abm::LocationType::Work);
}

TEST(TestMigrationRules, recover)
{
    auto hospital = mio::abm::Location(mio::abm::LocationType::Hospital, 0);
    auto p_rec =
        mio::abm::Person(hospital, mio::abm::InfectionState::Recovered_Infected, mio::abm::AgeGroup::Age60to79, {});
    auto p_inf =
        mio::abm::Person(hospital, mio::abm::InfectionState::Infected_Severe, mio::abm::AgeGroup::Age60to79, {});
    auto t  = mio::abm::TimePoint(12346);
    auto dt = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::return_home_when_recovered(p_rec, t, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::return_home_when_recovered(p_inf, t, dt, {}), mio::abm::LocationType::Hospital);
}

TEST(TestWorld, evolveMigration)
{
    using testing::Return;

    {
        auto world     = mio::abm::World();
        auto home_id   = world.add_location(mio::abm::LocationType::Home);
        auto school_id = world.add_location(mio::abm::LocationType::School);
        auto work_id   = world.add_location(mio::abm::LocationType::Work);

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>>
            mock_uniform_dist;
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
            .Times(testing::AtLeast(2))
            .WillOnce(testing::Return(0.8))
            .WillOnce(testing::Return(0.8))
            .WillOnce(testing::Return(0.8))
            .WillOnce(testing::Return(0.8))
            .WillOnce(testing::Return(0.8))
            .WillOnce(testing::Return(0.8))
            .WillOnce(testing::Return(0.8))
            .WillOnce(testing::Return(0.8))
            .WillRepeatedly(testing::Return(1.0));

        auto& p1 = world.add_person(home_id, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34);
        auto& p2 = world.add_person(home_id, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14);

        p1.set_assigned_location(school_id);
        p2.set_assigned_location(school_id);
        p1.set_assigned_location(work_id);
        p2.set_assigned_location(work_id);
        p1.set_assigned_location(home_id);
        p2.set_assigned_location(home_id);

        auto& school = world.get_individualized_location(school_id);
        auto& work   = world.get_individualized_location(work_id);

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
            mock_exponential_dist;
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

        world.evolve(mio::abm::TimePoint(0) + mio::abm::hours(8), mio::abm::hours(1));

        EXPECT_EQ(p1.get_location_id().type, mio::abm::LocationType::Work);
        EXPECT_EQ(p2.get_location_id().type, mio::abm::LocationType::School);
        EXPECT_EQ(school.get_subpopulations().sum(), 1);
        EXPECT_EQ(work.get_subpopulations().sum(), 1);
    }

    {
        auto world = mio::abm::World();
        world.use_migration_rules(false);

        auto home_id     = world.add_location(mio::abm::LocationType::Home);
        auto event_id    = world.add_location(mio::abm::LocationType::SocialEvent);
        auto work_id     = world.add_location(mio::abm::LocationType::Work);
        auto hospital_id = world.add_location(mio::abm::LocationType::Hospital);

        auto& p1 = world.add_person(home_id, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34);
        auto& p2 = world.add_person(home_id, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14);
        auto& p3 = world.add_person(home_id, mio::abm::InfectionState::Infected_Severe, mio::abm::AgeGroup::Age5to14);
        auto& p4 =
            world.add_person(hospital_id, mio::abm::InfectionState::Recovered_Infected, mio::abm::AgeGroup::Age5to14);
        auto& p5 = world.add_person(home_id, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age15to34);
        p1.set_assigned_location(event_id);
        p2.set_assigned_location(event_id);
        p1.set_assigned_location(work_id);
        p2.set_assigned_location(work_id);
        p1.set_assigned_location(home_id);
        p2.set_assigned_location(home_id);
        p3.set_assigned_location(home_id);
        p4.set_assigned_location(home_id);
        p3.set_assigned_location(hospital_id);
        p5.set_assigned_location(event_id);
        p5.set_assigned_location(work_id);
        p5.set_assigned_location(home_id);

        mio::abm::TripList& data = world.get_trip_list();
        mio::abm::Trip trip1(p1.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), work_id, home_id);
        mio::abm::Trip trip2(p2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, home_id);
        mio::abm::Trip trip3(p5.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, work_id);
        data.add_trip(trip1);
        data.add_trip(trip2);
        data.add_trip(trip3);

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
            mock_exponential_dist;
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

        world.evolve(mio::abm::TimePoint(0) + mio::abm::hours(8), mio::abm::hours(2));

        auto& event    = world.get_individualized_location(event_id);
        auto& work     = world.get_individualized_location(work_id);
        auto& home     = world.get_individualized_location(home_id);
        auto& hospital = world.get_individualized_location(hospital_id);

        EXPECT_EQ(p1.get_location_id().type, mio::abm::LocationType::Work);
        EXPECT_EQ(p2.get_location_id().type, mio::abm::LocationType::SocialEvent);
        EXPECT_EQ(p3.get_location_id().type, mio::abm::LocationType::Hospital);
        EXPECT_EQ(p4.get_location_id().type, mio::abm::LocationType::Home);
        EXPECT_EQ(p5.get_location_id().type, mio::abm::LocationType::Home);
        EXPECT_EQ(event.get_subpopulations().sum(), 1);
        EXPECT_EQ(work.get_subpopulations().sum(), 1);
        EXPECT_EQ(home.get_subpopulations().sum(), 2);
        EXPECT_EQ(hospital.get_subpopulations().sum(), 1);
    }
}

TEST(TestSimulation, advance_random)
{
    auto world     = mio::abm::World();
    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto location2 = world.add_location(mio::abm::LocationType::School);
    auto& p1       = world.add_person(location1, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age5to14);
    auto& p2       = world.add_person(location1, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14);
    auto& p3       = world.add_person(location2, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14);
    auto& p4       = world.add_person(location2, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);
    p4.set_assigned_location(location2);

    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(world));

    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(50));
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

TEST(TestTestingRule, addremoveandevaluateTestRule)
{
    auto world   = mio::abm::World();
    auto home_id = world.add_location(mio::abm::LocationType::Home);
    auto work_id = world.add_location(mio::abm::LocationType::Work);
    auto person  = mio::abm::Person(home_id, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age15to34,
                                    world.get_global_infection_parameters());
    auto& home   = world.get_individualized_location(home_id);
    auto& work   = world.get_individualized_location(work_id);
    person.set_assigned_location(home);
    person.set_assigned_location(work);

    auto testing_rule = mio::abm::TestingRule({}, {}, {});
    testing_rule.add_infection_state(mio::abm::InfectionState::Infected);
    testing_rule.add_infection_state(mio::abm::InfectionState::Carrier);
    testing_rule.add_location_type(mio::abm::LocationType::Home);
    testing_rule.add_location_type(mio::abm::LocationType::Work);

    ASSERT_EQ(testing_rule.evaluate(person, work), true);
    ASSERT_EQ(testing_rule.evaluate(person, home), true);

    testing_rule.remove_infection_state(mio::abm::InfectionState::Infected);
    ASSERT_EQ(testing_rule.evaluate(person, home), false);

    testing_rule.add_infection_state(mio::abm::InfectionState::Infected);
    testing_rule.remove_location_type(mio::abm::LocationType::Home);
    ASSERT_EQ(testing_rule.evaluate(person, home), false);
}

TEST(TestTestingScheme, init)
{
    std::vector<mio::abm::InfectionState> test_infection_states1 = {mio::abm::InfectionState::Infected,
                                                                    mio::abm::InfectionState::Carrier};
    std::vector<mio::abm::LocationType> test_location_types1     = {mio::abm::LocationType::Home,
                                                                    mio::abm::LocationType::Work};

    auto testing_rule1 = mio::abm::TestingRule({}, test_location_types1, test_infection_states1);
    std::vector<mio::abm::TestingRule> testing_rules = {testing_rule1};

    const auto testing_frequency = mio::abm::days(1);
    const auto start_date        = mio::abm::TimePoint(0);
    const auto end_date          = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability       = 1.0;
    const auto test_type         = mio::abm::PCRTest();

    auto testing_scheme =
        mio::abm::TestingScheme(testing_rules, testing_frequency, start_date, end_date, probability, test_type);

    ASSERT_EQ(testing_scheme.is_active(), false);
    testing_scheme.update_activity_status(mio::abm::TimePoint(10));
    ASSERT_EQ(testing_scheme.is_active(), true);
    testing_scheme.update_activity_status(mio::abm::TimePoint(60 * 60 * 24 * 3 + 200));
    ASSERT_EQ(testing_scheme.is_active(), false);
}

TEST(TestTestingScheme, runScheme)
{
    std::vector<mio::abm::InfectionState> test_infection_states1 = {mio::abm::InfectionState::Infected,
                                                                    mio::abm::InfectionState::Carrier};
    std::vector<mio::abm::LocationType> test_location_types1     = {mio::abm::LocationType::Home,
                                                                    mio::abm::LocationType::Work};

    auto testing_rule1 = mio::abm::TestingRule({}, test_location_types1, test_infection_states1);
    std::vector<mio::abm::TestingRule> testing_rules = {testing_rule1};

    const auto testing_frequency = mio::abm::days(1);
    const auto start_date        = mio::abm::TimePoint(0);
    const auto end_date          = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability       = 1.0;
    const auto test_type         = mio::abm::PCRTest();

    auto testing_scheme =
        mio::abm::TestingScheme(testing_rules, testing_frequency, start_date, end_date, probability, test_type);

    std::vector<mio::abm::InfectionState> test_infection_states2 = {mio::abm::InfectionState::Recovered_Carrier};
    std::vector<mio::abm::LocationType> test_location_types2     = {mio::abm::LocationType::Home};
    auto testing_rule2 = mio::abm::TestingRule({}, test_location_types2, test_infection_states2);
    testing_scheme.add_testing_rule(testing_rule2);

    auto loc_home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto loc_work = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto person1  = mio::abm::Person(loc_home, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age15to34, {});
    auto person2 =
        mio::abm::Person(loc_home, mio::abm::InfectionState::Recovered_Carrier, mio::abm::AgeGroup::Age15to34, {});

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(4))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.5));
    testing_scheme.run_scheme(person1, loc_home);
    testing_scheme.run_scheme(person2, loc_work);

    ASSERT_EQ(person1.is_in_quarantine(), true);
    ASSERT_EQ(person2.is_in_quarantine(), false);
    ASSERT_EQ(person2.get_time_since_negative_test(), mio::abm::days(0));
}

TEST(TestWorldTestingRule, testAddingAndUpdatingAndRunningTestingSchemes)
{

    auto world   = mio::abm::World();
    auto home_id = world.add_location(mio::abm::LocationType::Home);
    auto work_id = world.add_location(mio::abm::LocationType::Work);
    auto person  = mio::abm::Person(home_id, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age15to34,
                                    world.get_global_infection_parameters());
    auto& home   = world.get_individualized_location(home_id);
    auto& work   = world.get_individualized_location(work_id);
    person.set_assigned_location(home);
    person.set_assigned_location(work);

    auto testing_rule = mio::abm::TestingRule({}, {}, {});
    testing_rule.add_infection_state(mio::abm::InfectionState::Infected);
    testing_rule.add_infection_state(mio::abm::InfectionState::Carrier);
    testing_rule.add_location_type(mio::abm::LocationType::Home);
    testing_rule.add_location_type(mio::abm::LocationType::Work);

    const auto testing_frequency = mio::abm::days(1);
    const auto start_date        = mio::abm::TimePoint(20);
    const auto end_date          = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability       = 1.0;
    const auto test_type         = mio::abm::PCRTest();

    auto testing_scheme =
        mio::abm::TestingScheme({testing_rule}, testing_frequency, start_date, end_date, probability, test_type);

    world.add_testing_scheme(testing_scheme);
    auto current_time = mio::abm::TimePoint(0);
    ASSERT_EQ(world.run_testing_schemes(person, work), true);
    current_time = mio::abm::TimePoint(30);
    world.update_testing_scheme_activity_status(current_time);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(4))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.4));
    ASSERT_EQ(world.run_testing_schemes(person, work), false);
}