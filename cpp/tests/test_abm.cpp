/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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

TEST(TestLocation, GetIndex)
{
    auto location = mio::Location(mio::LocationType::Home, 0);
    ASSERT_EQ((int)location.get_index(), 0);
}

TEST(TestLocation, addRemovePerson)
{
    auto location = mio::Location(mio::LocationType::Home, 0);
    auto person1  = mio::Person(location, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14, {});
    location.add_person(person1);
    auto person2 = mio::Person(location, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, {});
    location.add_person(person2);
    auto person3 = mio::Person(location, mio::InfectionState::Exposed, mio::AbmAgeGroup::Age35to59, {});
    location.add_person(person3);

    ASSERT_EQ(location.get_subpopulation(mio::InfectionState::Susceptible), 2);
    ASSERT_EQ(location.get_subpopulation(mio::InfectionState::Exposed), 1);

    location.remove_person(person2);

    ASSERT_EQ(location.get_subpopulation(mio::InfectionState::Susceptible), 1);
    ASSERT_EQ(location.get_subpopulation(mio::InfectionState::Exposed), 1);
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
    EXPECT_CALL(mock_uniform_int_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    auto person = mio::Person(location, mio::InfectionState::Exposed, mio::AbmAgeGroup::Age60to79, {});
    location.add_person(person);

    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Exposed);
    ASSERT_EQ(person.get_location_id().index, location.get_index());
    ASSERT_EQ(person.get_location_id().type, location.get_type());
    ASSERT_EQ(person.get_age(), mio::AbmAgeGroup::Age60to79);

    mio::TimeSpan dt = mio::hours(1);
    person.interact(dt, {}, location, {});
    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Carrier);
}

TEST(TestPerson, migrate)
{
    auto loc1   = mio::Location(mio::LocationType::Work, 0);
    auto loc2   = mio::Location(mio::LocationType::School, 0);
    auto person = mio::Person(loc1, mio::InfectionState::Recovered_Carrier, mio::AbmAgeGroup::Age0to4, {});
    loc1.add_person(person);

    person.migrate_to(loc1, loc2);

    ASSERT_EQ(person.get_location_id().index, loc2.get_index());
    ASSERT_EQ(person.get_location_id().type, loc2.get_type());
    ASSERT_EQ(loc2.get_subpopulation(mio::InfectionState::Recovered_Carrier), 1);
    ASSERT_EQ(loc1.get_subpopulation(mio::InfectionState::Recovered_Carrier), 0);
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
        mio::Person(location, mio::InfectionState::Carrier, vaccination_state, mio::AbmAgeGroup::Age15to34, params);
    location.add_person(infected1);
    auto infected2 =
        mio::Person(location, mio::InfectionState::Infected, vaccination_state, mio::AbmAgeGroup::Age80plus, params);
    location.add_person(infected2);
    auto infected3 =
        mio::Person(location, mio::InfectionState::Infected, vaccination_state, mio::AbmAgeGroup::Age5to14, params);
    location.add_person(infected3);

    //cache precomputed results
    auto dt = mio::seconds(8640); //0.1 days
    location.begin_step(dt, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    {
        auto susceptible = mio::Person(location, mio::InfectionState::Susceptible, vaccination_state, age, params);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(susceptible, dt, params), mio::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.15));
        EXPECT_EQ(location.interact(susceptible, dt, params), mio::InfectionState::Susceptible);
    }

    {
        auto exposed = mio::Person(location, mio::InfectionState::Exposed, vaccination_state, age, params);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(0); //no transitions out of exposed state
        EXPECT_EQ(location.interact(exposed, dt, params), mio::InfectionState::Exposed);
    }

    {
        auto carrier = mio::Person(location, mio::InfectionState::Carrier, vaccination_state, age, params);

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
        auto infected = mio::Person(location, mio::InfectionState::Infected, vaccination_state, age, params);

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
        auto severe = mio::Person(location, mio::InfectionState::Infected_Severe, vaccination_state, age, params);

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
        auto critical = mio::Person(location, mio::InfectionState::Infected_Critical, vaccination_state, age, params);

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
        auto recovered = mio::Person(location, recovered_state, vaccination_state, age, params);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(recovered, dt, params), mio::InfectionState::Susceptible);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(recovered, dt, params), recovered_state);
    }
}

TEST(TestPerson, quarantine)
{
    using testing::Return;

    auto infection_parameters = mio::GlobalInfectionParameters();
    auto home                 = mio::Location(mio::LocationType::Home, 0);
    auto work                 = mio::Location(mio::LocationType::Work, 0);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));

    auto person = mio::Person(home, mio::InfectionProperties(mio::InfectionState::Infected, true),
                              mio::AbmAgeGroup::Age15to34, infection_parameters);
    home.add_person(person);

    auto t_morning = mio::TimePoint(0) + mio::hours(7);
    auto dt        = mio::hours(1);

    ASSERT_EQ(mio::go_to_work(person, t_morning, dt, {}), mio::LocationType::Home);

    //setup rng mock so the person has a state transition to Recovered_Infected
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.04));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    person.interact(dt, infection_parameters, home, {});
    ASSERT_EQ(person.get_infection_state(), mio::InfectionState::Recovered_Infected);
    ASSERT_EQ(mio::go_to_work(person, t_morning, dt, {}), mio::LocationType::Work);
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

    auto home    = mio::Location(mio::LocationType::Home, 0);
    auto p_child = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14, {});
    auto p_adult = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, {});

    auto t_morning = mio::TimePoint(0) + mio::hours(7);
    auto t_weekend = mio::TimePoint(0) + mio::days(5) + mio::hours(7);
    auto dt        = mio::hours(1);

    ASSERT_EQ(mio::go_to_school(p_child, t_morning, dt, {}), mio::LocationType::School);
    ASSERT_EQ(mio::go_to_school(p_adult, t_morning, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_school(p_child, t_weekend, dt, {}), mio::LocationType::Home);
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

    auto home = mio::Location(mio::LocationType::Home, 0);
    auto p_child_goes_to_school_at_6 =
        mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14, {});
    auto p_child_goes_to_school_at_8 =
        mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14, {});

    auto t_morning_6 = mio::TimePoint(0) + mio::hours(6);
    auto t_morning_8 = mio::TimePoint(0) + mio::hours(8);
    auto dt          = mio::hours(1);

    ASSERT_EQ(mio::go_to_school(p_child_goes_to_school_at_6, t_morning_6, dt, {}), mio::LocationType::School);
    ASSERT_EQ(mio::go_to_school(p_child_goes_to_school_at_6, t_morning_8, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_school(p_child_goes_to_school_at_8, t_morning_6, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_school(p_child_goes_to_school_at_8, t_morning_8, dt, {}), mio::LocationType::School);
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

    auto home = mio::Location(mio::LocationType::Home, 0);
    auto p_child_goes_to_school_at_6 =
        mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14, {});
    auto p_child_goes_to_school_at_8_30 =
        mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14, {});

    auto t_morning_6    = mio::TimePoint(0) + mio::hours(6);
    auto t_morning_8_30 = mio::TimePoint(0) + mio::hours(8) + mio::seconds(1800);
    auto dt             = mio::seconds(1800);

    ASSERT_EQ(mio::go_to_school(p_child_goes_to_school_at_6, t_morning_6, dt, {}), mio::LocationType::School);
    ASSERT_EQ(mio::go_to_school(p_child_goes_to_school_at_6, t_morning_8_30, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_school(p_child_goes_to_school_at_8_30, t_morning_6, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_school(p_child_goes_to_school_at_8_30, t_morning_8_30, dt, {}), mio::LocationType::School);
}

TEST(TestMigrationRules, school_return)
{
    auto school  = mio::Location(mio::LocationType::School, 0);
    auto p_child = mio::Person(school, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14, {});

    auto t  = mio::TimePoint(0) + mio::hours(15);
    auto dt = mio::hours(1);

    ASSERT_EQ(mio::go_to_school(p_child, t, dt, {}), mio::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work)
{
    auto home = mio::Location(mio::LocationType::Home, 0);
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

    auto p_retiree = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age60to79, {});
    auto p_adult   = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, {});

    auto t_morning = mio::TimePoint(0) + mio::hours(8);
    auto t_night   = mio::TimePoint(0) + mio::days(1) + mio::hours(4);
    auto dt        = mio::hours(1);

    ASSERT_EQ(mio::go_to_work(p_retiree, t_morning, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_work(p_adult, t_morning, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_work(p_adult, t_night, dt, {}), mio::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work_with_non_dividable_timespan)
{
    auto home = mio::Location(mio::LocationType::Home, 0);
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

    auto p_retiree = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age60to79, {});
    auto p_adult   = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, {});

    auto t_morning = mio::TimePoint(0) + mio::hours(8);
    auto t_night   = mio::TimePoint(0) + mio::days(1) + mio::hours(4);
    auto dt        = mio::minutes(53);

    ASSERT_EQ(mio::go_to_work(p_retiree, t_morning, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_work(p_adult, t_morning, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_work(p_adult, t_night, dt, {}), mio::LocationType::Home);
}

TEST(TestMigrationRules, workers_go_to_work_in_different_times)
{
    auto home = mio::Location(mio::LocationType::Home, 0);
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
        mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, {});
    auto p_adult_goes_to_work_at_8 =
        mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, {});

    auto t_morning_6 = mio::TimePoint(0) + mio::hours(6);
    auto t_morning_8 = mio::TimePoint(0) + mio::hours(8);
    auto t_night     = mio::TimePoint(0) + mio::days(1) + mio::hours(4);
    auto dt          = mio::hours(1);

    ASSERT_EQ(mio::go_to_work(p_adult_goes_to_work_at_6, t_morning_6, dt, {}), mio::LocationType::Work);
    ASSERT_EQ(mio::go_to_work(p_adult_goes_to_work_at_6, t_morning_8, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_work(p_adult_goes_to_work_at_6, t_night, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_work(p_adult_goes_to_work_at_8, t_morning_6, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_work(p_adult_goes_to_work_at_8, t_morning_8, dt, {}), mio::LocationType::Work);
    ASSERT_EQ(mio::go_to_work(p_adult_goes_to_work_at_8, t_night, dt, {}), mio::LocationType::Home);
}

TEST(TestMigrationRules, work_return)
{
    auto work    = mio::Location(mio::LocationType::Work, 0);
    auto p_adult = mio::Person(work, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age35to59, {});

    auto t  = mio::TimePoint(0) + mio::hours(17);
    auto dt = mio::hours(1);

    ASSERT_EQ(mio::go_to_work(p_adult, t, dt, {}), mio::LocationType::Home);
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

TEST(TestMigrationRules, go_shopping)
{
    auto hospital = mio::Location(mio::LocationType::Hospital, 0);
    auto p_hosp   = mio::Person(hospital, mio::InfectionState::Infected, mio::AbmAgeGroup::Age0to4, {});
    auto home     = mio::Location(mio::LocationType::Home, 0);
    auto p_home   = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age60to79, {});

    auto t_weekday = mio::TimePoint(0) + mio::days(4) + mio::hours(9);
    auto t_sunday  = mio::TimePoint(0) + mio::days(6) + mio::hours(9);
    auto t_night   = mio::TimePoint(0) + mio::days(4) + mio::hours(1);
    auto dt        = mio::hours(1);

    ASSERT_EQ(mio::go_to_shop(p_hosp, t_weekday, dt, {}), mio::LocationType::Hospital);
    ASSERT_EQ(mio::go_to_shop(p_home, t_sunday, dt, {}), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_shop(p_home, t_night, dt, {}), mio::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::go_to_shop(p_home, t_weekday, dt, {}), mio::LocationType::BasicsShop);
}

TEST(TestMigrationRules, shop_return)
{
    auto t  = mio::TimePoint(0) + mio::days(4) + mio::hours(9);
    auto dt = mio::hours(1);

    auto home = mio::Location(mio::LocationType::Home, 0);
    auto shop = mio::Location(mio::LocationType::BasicsShop, 0);
    auto p    = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34, {});
    home.add_person(p);
    p.migrate_to(home, shop);
    p.interact(dt, {}, shop, {}); //person only returns home after some time passed

    ASSERT_EQ(mio::go_to_shop(p, t, dt, {}), mio::LocationType::Home);
}

TEST(TestMigrationRules, go_event)
{
    auto work   = mio::Location(mio::LocationType::Work, 0);
    auto p_work = mio::Person(work, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age35to59, {});
    auto home   = mio::Location(mio::LocationType::Home, 0);
    auto p_home = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age60to79, {});

    auto t_weekday  = mio::TimePoint(0) + mio::days(4) + mio::hours(20);
    auto t_saturday = mio::TimePoint(0) + mio::days(5) + mio::hours(10);
    auto t_night    = mio::TimePoint(0) + mio::days(5) + mio::hours(1);
    auto dt         = mio::hours(1);

    ASSERT_EQ(mio::go_to_event(p_work, t_weekday, dt, {}), mio::LocationType::Work);
    ASSERT_EQ(mio::go_to_event(p_home, t_night, dt, {}), mio::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::go_to_event(p_home, t_weekday, dt, {}), mio::LocationType::SocialEvent);

    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::go_to_event(p_home, t_saturday, dt, {}), mio::LocationType::SocialEvent);
}

TEST(TestMigrationRules, event_return)
{
    auto t  = mio::TimePoint(0) + mio::days(4) + mio::hours(21);
    auto dt = mio::hours(3);

    auto home = mio::Location(mio::LocationType::Home, 0);
    auto shop = mio::Location(mio::LocationType::SocialEvent, 0);
    auto p    = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34, {});
    home.add_person(p);
    p.migrate_to(home, shop);
    p.interact(dt, {}, shop, {});

    ASSERT_EQ(mio::go_to_event(p, t, dt, {}), mio::LocationType::Home);
}

TEST(TestLockdownRules, school_closure)
{
    auto t         = mio::TimePoint(0);
    auto dt        = mio::hours(1);
    auto t_morning = mio::TimePoint(0) + mio::hours(6);
    auto home      = mio::Location(mio::LocationType::Home, 0);
    auto school    = mio::Location(mio::LocationType::School, 0);

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

    auto p1 = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age5to14, {});
    p1.set_assigned_location(home);
    p1.set_assigned_location(school);
    auto p2 = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age5to14, {});
    p2.set_assigned_location(home);
    p2.set_assigned_location(school);
    mio::AbmMigrationParameters params;

    mio::set_school_closure(t, 0.7, params);

    ASSERT_EQ(mio::go_to_school(p1, t_morning, dt, params), mio::LocationType::Home);
    ASSERT_EQ(mio::go_to_school(p2, t_morning, dt, params), mio::LocationType::School);
}

TEST(TestLockdownRules, school_opening)
{
    auto t_closing = mio::TimePoint(0);
    auto t_opening = mio::TimePoint(0) + mio::days(1);
    auto dt        = mio::hours(1);
    auto t_morning = mio::TimePoint(0) + mio::days(1) + mio::hours(7);
    auto home      = mio::Location(mio::LocationType::Home, 0);
    auto school    = mio::Location(mio::LocationType::School, 0);
    //setup rng mock so the person is homeschooled in case of lockdown
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));
    auto p = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age5to14, {});
    p.set_assigned_location(home);
    p.set_assigned_location(school);
    mio::AbmMigrationParameters params;

    mio::set_school_closure(t_closing, 1., params);
    mio::set_school_closure(t_opening, 0., params);

    ASSERT_EQ(mio::go_to_school(p, t_morning, dt, params), mio::LocationType::School);
}

TEST(TestLockdownRules, home_office)
{
    auto t         = mio::TimePoint(0);
    auto t_morning = mio::TimePoint(0) + mio::hours(8);
    auto dt        = mio::hours(1);
    auto home      = mio::Location(mio::LocationType::Home, 0);
    auto work      = mio::Location(mio::LocationType::Work, 0);
    mio::AbmMigrationParameters params;

    mio::set_home_office(t, 0.4, params);

    //setup rng mock so one person goes to work and the other works at home
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(4))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillRepeatedly(testing::Return(1.0));

    auto person1 = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, {});
    auto person2 = mio::Person(home, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age15to34, {});
    person1.set_assigned_location(home);
    person1.set_assigned_location(work);
    person2.set_assigned_location(home);
    person2.set_assigned_location(work);

    ASSERT_EQ(mio::go_to_work(person1, t_morning, dt, params), mio::LocationType::Work);
    ASSERT_EQ(mio::go_to_work(person2, t_morning, dt, params), mio::LocationType::Home);
}

TEST(TestLockdownRules, no_home_office)
{
    auto t_closing = mio::TimePoint(0);
    auto t_opening = mio::TimePoint(0) + mio::days(1);
    auto dt        = mio::hours(1);
    auto t_morning = mio::TimePoint(0) + mio::days(1) + mio::hours(8);
    auto home      = mio::Location(mio::LocationType::Home, 0);
    auto work      = mio::Location(mio::LocationType::Work, 0);

    //setup rng mock so the person works in home office
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillRepeatedly(testing::Return(1.0));

    auto p = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34, {});
    p.set_assigned_location(home);
    p.set_assigned_location(work);
    mio::AbmMigrationParameters params;

    mio::set_home_office(t_closing, 0.5, params);
    mio::set_home_office(t_opening, 0., params);

    ASSERT_EQ(mio::go_to_work(p, t_morning, dt, params), mio::LocationType::Work);
}

TEST(TestLockdownRules, social_event_closure)
{
    auto t         = mio::TimePoint(0);
    auto dt        = mio::hours(1);
    auto t_evening = mio::TimePoint(0) + mio::hours(19);
    auto home      = mio::Location(mio::LocationType::Home, 0);
    auto event     = mio::Location(mio::LocationType::SocialEvent, 0);
    auto p         = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age5to14, {});
    p.set_assigned_location(home);
    p.set_assigned_location(event);
    mio::AbmMigrationParameters params;

    mio::close_social_events(t, 1, params);

    ASSERT_EQ(mio::go_to_event(p, t_evening, dt, params), mio::LocationType::Home);
}

TEST(TestLockdownRules, social_events_opening)
{
    auto t_closing = mio::TimePoint(0);
    auto t_opening = mio::TimePoint(0) + mio::days(1);
    auto dt        = mio::hours(1);
    auto t_evening = mio::TimePoint(0) + mio::days(1) + mio::hours(19);
    auto home      = mio::Location(mio::LocationType::Home, 0);
    auto event     = mio::Location(mio::LocationType::SocialEvent, 0);
    auto p         = mio::Person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age5to14, {});
    p.set_assigned_location(event);
    p.set_assigned_location(home);
    mio::AbmMigrationParameters params;

    mio::close_social_events(t_closing, 1, params);
    mio::close_social_events(t_opening, 0, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::go_to_event(p, t_evening, dt, params), mio::LocationType::SocialEvent);
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

    auto world     = mio::World();
    auto home_id   = world.add_location(mio::LocationType::Home);
    auto school_id = world.add_location(mio::LocationType::School);
    auto work_id   = world.add_location(mio::LocationType::Work);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
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

    auto& p1 = world.add_person(home_id, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age15to34);
    auto& p2 = world.add_person(home_id, mio::InfectionState::Susceptible, mio::AbmAgeGroup::Age5to14);
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

    world.evolve(mio::TimePoint(0) + mio::hours(8), mio::hours(1));

    EXPECT_EQ(p1.get_location_id().type, mio::LocationType::Work);
    EXPECT_EQ(p2.get_location_id().type, mio::LocationType::School);
    EXPECT_EQ(school.get_subpopulations().sum(), 1);
    EXPECT_EQ(work.get_subpopulations().sum(), 1);
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

TEST(TestWorldBuilder, averageContacts)
{
    Eigen::VectorXd x_sol(5);
    x_sol << 30, 2, 10, 22, 100;
    Eigen::MatrixXd contact_matrix(2, 2);
    Eigen::MatrixXd average_contacs(2, 2);
    contact_matrix << 24, 7, 11.6666666, 19.3333333;
    mio::compute_average_contact_matrix(average_contacs, x_sol, 2);
    ASSERT_TRUE(std::abs(average_contacs(0, 0) - contact_matrix(0, 0)) < 0.00001);
    ASSERT_TRUE(std::abs(average_contacs(0, 1) - contact_matrix(0, 1)) < 0.00001);
    ASSERT_TRUE(std::abs(average_contacs(1, 0) - contact_matrix(1.0)) < 0.00001);
    ASSERT_TRUE(std::abs(average_contacs(1, 1) - contact_matrix(1, 1)) < 0.00001);
}

TEST(TestWorldBuilder, findOptimalLocations)
{
    Eigen::VectorXd people(2);
    people << 40, 24;
    Eigen::MatrixXd contact_matrix(2, 2);
    contact_matrix << 24, 7, 11.6666666, 19.3333333;
    int num_locs              = 2;
    Eigen::VectorXd size_locs = Eigen::VectorXd::Constant(2, 22);
    Eigen::VectorXd x_sol     = mio::find_optimal_locations(people, num_locs, contact_matrix, size_locs);
    Eigen::MatrixXd average_contacs(2, 2);
    mio::compute_average_contact_matrix(average_contacs, x_sol, 2);
    ASSERT_TRUE(std::abs(average_contacs(0, 0) - contact_matrix(0, 0)) < 0.0001);
    ASSERT_TRUE(std::abs(average_contacs(0, 1) - contact_matrix(0, 1)) < 0.0001);
    ASSERT_TRUE(std::abs(average_contacs(1, 0) - contact_matrix(1.0)) < 0.0001);
    ASSERT_TRUE(std::abs(average_contacs(1, 1) - contact_matrix(1, 1)) < 0.0001);
}

TEST(TestWorldBuilder, create_locations)
{
    //create environment with only 2 age groups and small contact matrix
    Eigen::MatrixXd M(2, 2);
    M(0, 0)           = 24;
    M(0, 1)           = 7;
    M(1, 0)           = 11.666666666;
    M(1, 1)           = 19.333333333;
    uint32_t num_locs = 2;

    auto world = mio::World();
    auto home  = world.add_location(mio::LocationType::Home);
    for (int i = 0; i < 40; i++) {
        auto& p1 = world.add_person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age0to4);
        p1.set_assigned_location(home);
    }
    for (int i = 0; i < 24; i++) {
        auto& p1 = world.add_person(home, mio::InfectionState::Carrier, mio::AbmAgeGroup::Age5to14);
        p1.set_assigned_location(home);
    }
    Eigen::VectorXd size_locs(2);
    size_locs << 22, 22;

    //assign social event
    mio::create_locations(num_locs, mio::LocationType::SocialEvent, world, M, size_locs);

    //count people of each age group at the locations
    Eigen::VectorXi counter = Eigen::VectorXi::Zero(4);
    for (auto& p : world.get_persons()) {
        int index = p.get_assigned_location_index(mio::LocationType::SocialEvent);
        if (index == 0 && (size_t)p.get_age() == 0) {
            counter(0)++;
        }
        else if (index == 0 && (size_t)p.get_age() == 1) {
            counter(1)++;
        }
        else if (index == 1 && (size_t)p.get_age() == 0) {
            counter(2)++;
        }
        else if (index == 1 && (size_t)p.get_age() == 1) {
            counter(3)++;
        }
    }
    if (counter(0) > counter(2)) {
        ASSERT_EQ(counter(0), 30);
        ASSERT_EQ(counter(1), 2);
        ASSERT_EQ(counter(2), 10);
        ASSERT_EQ(counter(3), 22);
    }
    else {
        ASSERT_EQ(counter(2), 30);
        ASSERT_EQ(counter(3), 2);
        ASSERT_EQ(counter(0), 10);
        ASSERT_EQ(counter(1), 22);
    }
}
