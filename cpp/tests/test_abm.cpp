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
#include "epidemiology/abm/abm.h"
#include "epidemiology/abm/migration_rules.h"
#include "epidemiology/abm/lockdown_rules.h"
#include "epidemiology/utils/eigen_util.h"
#include "matchers.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <memory>

TEST(TestLocation, init)
{
    auto location = epi::Location(epi::LocationType::School, 0);
    for (epi::InfectionState i = epi::InfectionState(0); i < epi::InfectionState::Count;
         i                     = epi::InfectionState(size_t(i) + 1)) {
        ASSERT_EQ(location.get_subpopulation(i), 0);
    }
    ASSERT_EQ(print_wrap(location.get_subpopulations()),
              print_wrap(Eigen::VectorXi::Zero(Eigen::Index(epi::InfectionState::Count))));
}

TEST(TestLocation, GetIndex)
{
    auto location = epi::Location(epi::LocationType::Home, 0);
    ASSERT_EQ(location.get_index(), 0);
}

TEST(TestLocation, addRemovePerson)
{
    auto location = epi::Location(epi::LocationType::Home, 0);
    auto person1  = epi::Person(location, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age5to14, {}, {});
    location.add_person(person1);
    auto person2 = epi::Person(location, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34, {}, {});
    location.add_person(person2);
    auto person3 = epi::Person(location, epi::InfectionState::Exposed, epi::AbmAgeGroup::Age35to59, {}, {});
    location.add_person(person3);

    ASSERT_EQ(location.get_subpopulation(epi::InfectionState::Susceptible), 2);
    ASSERT_EQ(location.get_subpopulation(epi::InfectionState::Exposed), 1);

    location.remove_person(person2);

    ASSERT_EQ(location.get_subpopulation(epi::InfectionState::Susceptible), 1);
    ASSERT_EQ(location.get_subpopulation(epi::InfectionState::Exposed), 1);
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

    auto location = epi::Location(epi::LocationType::Work, 0);
    //setup rng mock so the time_until_carrier is 1.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformIntDistribution<int>>>>
        mock_uniform_int_dist;
    EXPECT_CALL(mock_uniform_int_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
    
    auto person = epi::Person(location, epi::InfectionState::Exposed, epi::AbmAgeGroup::Age60to79, {}, {});
    location.add_person(person);
    
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);
    ASSERT_EQ(person.get_location_id().index, location.get_index());
    ASSERT_EQ(person.get_location_id().type, location.get_type());
    ASSERT_EQ(person.get_age(), epi::AbmAgeGroup::Age60to79);

    epi::TimeSpan dt = epi::hours(1);
    person.interact(dt, {}, location, {});
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Carrier);
}

TEST(TestPerson, migrate)
{
    auto loc1 = epi::Location(epi::LocationType::Work, 0);
    auto loc2 = epi::Location(epi::LocationType::School, 0);
    auto person    = epi::Person(loc1, epi::InfectionState::Recovered_Carrier, epi::AbmAgeGroup::Age0to4, {}, {});
    loc1.add_person(person);

    person.migrate_to(loc1, loc2);

    ASSERT_EQ(person.get_location_id().index, loc2.get_index());
    ASSERT_EQ(person.get_location_id().type, loc2.get_type());
    ASSERT_EQ(loc2.get_subpopulation(epi::InfectionState::Recovered_Carrier), 1);
    ASSERT_EQ(loc1.get_subpopulation(epi::InfectionState::Recovered_Carrier), 0);
}

TEST(TestPerson, setGetAssignedLocation)
{
    auto location = epi::Location(epi::LocationType::Work, 2);
    auto person = epi::Person(location, epi::InfectionState::Recovered_Carrier, epi::AbmAgeGroup::Age60to79, {}, {});
    person.set_assigned_location(location);
    ASSERT_EQ(person.get_assigned_location_index(epi::LocationType::Work), 2);

    person.set_assigned_location({4, epi::LocationType::Work});
    ASSERT_EQ(person.get_assigned_location_index(epi::LocationType::Work), 4);
}

TEST(TestWorld, findLocation)
{
    auto world = epi::World();
    auto home_id = world.add_location(epi::LocationType::Home);
    auto school_id = world.add_location(epi::LocationType::School);
    auto work_id = world.add_location(epi::LocationType::Work);
    auto person = epi::Person(home_id, epi::InfectionState::Recovered_Carrier, epi::AbmAgeGroup::Age60to79, world.get_global_infection_parameters(), world.get_global_testing_parameters());
    auto& home =world.get_individualized_location(home_id);
    auto& school =world.get_individualized_location(school_id);
    auto& work =world.get_individualized_location(work_id);
    person.set_assigned_location(home);
    person.set_assigned_location(school);
    person.set_assigned_location({0, epi::LocationType::Work});

    ASSERT_EQ(world.find_location(epi::LocationType::Work, person), &work);
    ASSERT_EQ(world.find_location(epi::LocationType::School, person), &school);
    ASSERT_EQ(world.find_location(epi::LocationType::Home, person), &home);
}

TEST(TestLocation, interact)
{
    using testing::Return;

    //test should work identically work with any age
    epi::AbmAgeGroup age = epi::AbmAgeGroup(epi::UniformIntDistribution<int>()(0, int(epi::AbmAgeGroup::Count) - 1));
    epi::GlobalInfectionParameters params;
    params.set<epi::CarrierToInfected>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::CarrierToInfected>()[{age}] = 0.5;
    params.set<epi::CarrierToRecovered>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::CarrierToRecovered>()[{age}] = 0.5;
    params.set<epi::DetectInfection>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::DetectInfection>()[{age}] = 0.5;
    params.set<epi::InfectedToSevere>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::InfectedToSevere>()[{age}] = 0.5;
    params.set<epi::InfectedToRecovered>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::InfectedToRecovered>()[{age}] = 0.5;
    params.set<epi::SevereToCritical>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::SevereToCritical>()[{age}] = 0.5;
    params.set<epi::SevereToRecovered>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::SevereToRecovered>()[{age}] = 0.5;
    params.set<epi::CriticalToDead>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::CriticalToDead>()[{age}] = 0.5;
    params.set<epi::CriticalToRecovered>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::CriticalToRecovered>()[{age}] = 0.5;
    params.set<epi::RecoveredToSusceptible>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::RecoveredToSusceptible>()[{age}] = 0.5;
    params.set<epi::SusceptibleToExposedByCarrier>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::SusceptibleToExposedByCarrier>()[{age}] = 0.5;
    params.set<epi::SusceptibleToExposedByInfected>({{epi::AbmAgeGroup::Count}, 0.});
    params.get<epi::SusceptibleToExposedByInfected>()[{age}] = 0.5;

    //setup location with some chance of exposure
    auto location  = epi::Location(epi::LocationType::Work, 0);
    auto infected1 = epi::Person(location, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age15to34, params, {});
    location.add_person(infected1);
    auto infected2 = epi::Person(location, epi::InfectionState::Infected, epi::AbmAgeGroup::Age80plus, params, {});
    location.add_person(infected2);
    auto infected3 = epi::Person(location, epi::InfectionState::Infected, epi::AbmAgeGroup::Age5to14, params, {});
    location.add_person(infected3);
    

    //cache precomputed results
    auto dt = epi::seconds(8640); //0.1 days
    location.begin_step(dt, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    {
        auto susceptible = epi::Person(location, epi::InfectionState::Susceptible, age, params, {});
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(susceptible, dt, params), epi::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.15));
        EXPECT_EQ(location.interact(susceptible, dt, params), epi::InfectionState::Susceptible);
    }

    {
        auto exposed = epi::Person(location, epi::InfectionState::Exposed, age, params, {});
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(0); //no transitions out of exposed state
        EXPECT_EQ(location.interact(exposed, dt, params), epi::InfectionState::Exposed);
    }

    {
        auto carrier = epi::Person(location, epi::InfectionState::Carrier, age, params, {});

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(carrier, dt, params), epi::InfectionState::Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.099));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(carrier, dt, params), epi::InfectionState::Recovered_Carrier);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(carrier, dt, params), epi::InfectionState::Carrier);
    }

    {
        auto infected = epi::Person(location, epi::InfectionState::Infected, age, params, {});

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(infected, dt, params), epi::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(infected, dt, params), epi::InfectionState::Infected_Severe);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(infected, dt, params), epi::InfectionState::Infected);
    }

    {
        auto severe = epi::Person(location, epi::InfectionState::Infected_Severe, age, params, {});

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(severe, dt, params), epi::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(severe, dt, params), epi::InfectionState::Infected_Critical);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(severe, dt, params), epi::InfectionState::Infected_Severe);
    }

    {
        auto critical = epi::Person(location, epi::InfectionState::Infected_Critical, age, params, {});

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(critical, dt, params), epi::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(critical, dt, params), epi::InfectionState::Dead);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(critical, dt, params), epi::InfectionState::Infected_Critical);
    }

    for (auto&& recovered_state : {epi::InfectionState::Recovered_Carrier, epi::InfectionState::Recovered_Infected}) {
        auto recovered = epi::Person(location, recovered_state, age, params, {});

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(recovered, dt, params), epi::InfectionState::Susceptible);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(recovered, dt, params), recovered_state);
    }
}

TEST(TestPerson, quarantine)
{
    using testing::Return;

    auto infection_parameters = epi::GlobalInfectionParameters();
    auto home = epi::Location(epi::LocationType::Home, 0);
    auto work = epi::Location(epi::LocationType::Work, 0);

    //setup rng mock so the person will be tested positive
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformDistribution<double>>>>
        mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(4).WillOnce(Return(0.04))
                                                            .WillOnce(Return(0.04))
                                                            .WillOnce(Return(0.04))
                                                            .WillOnce(Return(0.04));

    auto person   = epi::Person(home, epi::InfectionState::Infected, epi::AbmAgeGroup::Age15to34, infection_parameters, {});
    home.add_person(person);

    auto t_morning  = epi::TimePoint(0) + epi::hours(8);
    auto dt = epi::hours(1);

    ASSERT_EQ(epi::go_to_work(person, t_morning, dt, {}), epi::LocationType::Home);

    //setup rng mock so the person has a state transition to Recovered_Infected
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.04));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    person.interact(dt, infection_parameters, home, {});
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Recovered_Infected);
    ASSERT_EQ(epi::go_to_work(person, t_morning, dt, {}), epi::LocationType::Work);
}

TEST(TestPerson, interact)
{
    using testing::Return;

    auto infection_parameters = epi::GlobalInfectionParameters();
    auto loc = epi::Location(epi::LocationType::Home, 0);
    auto person   = epi::Person(loc, epi::InfectionState::Infected, epi::AbmAgeGroup::Age15to34, infection_parameters, {});
    loc.add_person(person);
    auto dt = epi::seconds(8640); //0.1 days
    loc.begin_step(dt, {});

    //setup rng mock so the person has a state transition
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    person.interact(dt, infection_parameters, loc, {});
    EXPECT_EQ(person.get_infection_state(), epi::InfectionState::Recovered_Infected);
    EXPECT_EQ(loc.get_subpopulation(epi::InfectionState::Recovered_Infected), 1);
    EXPECT_EQ(loc.get_subpopulation(epi::InfectionState::Infected), 0);
}

TEST(TestPerson, interact_exposed)
{
    using testing::Return;

    auto infection_parameters = epi::GlobalInfectionParameters();
    infection_parameters.set<epi::IncubationPeriod>({{epi::AbmAgeGroup::Count}, 2.});

    //setup location with some chance of exposure
    auto loc  = epi::Location(epi::LocationType::Work, 0);
    auto infected1 = epi::Person(loc, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age15to34, infection_parameters, {});
    loc.add_person(infected1);
    auto infected2 = epi::Person(loc, epi::InfectionState::Infected, epi::AbmAgeGroup::Age5to14, infection_parameters, {});
    loc.add_person(infected2);
    auto infected3 = epi::Person(loc, epi::InfectionState::Infected, epi::AbmAgeGroup::Age60to79, infection_parameters, {});
    loc.add_person(infected3);
    auto person = epi::Person(loc, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34, infection_parameters, {});
    loc.add_person(person);
    loc.begin_step(epi::hours(1), {});

    //setup rng mock so the person becomes exposed
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.49));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    //person becomes exposed
    person.interact(epi::hours(12), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);
    EXPECT_EQ(loc.get_subpopulation(epi::InfectionState::Exposed), 1);
    EXPECT_EQ(loc.get_subpopulation(epi::InfectionState::Carrier), 1);
    EXPECT_EQ(loc.get_subpopulation(epi::InfectionState::Infected), 2);

    //person becomes a carrier after the incubation time runs out, not random
    person.interact(epi::hours(12), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);

    person.interact(epi::hours(12), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);

    person.interact(epi::hours(24), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);

    person.interact(epi::hours(1), infection_parameters, loc, {});
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Carrier);
    EXPECT_EQ(loc.get_subpopulation(epi::InfectionState::Exposed), 0);
    EXPECT_EQ(loc.get_subpopulation(epi::InfectionState::Carrier), 2);
    EXPECT_EQ(loc.get_subpopulation(epi::InfectionState::Infected), 2);
}

TEST(TestWorld, init)
{
    auto world = epi::World();
    for (uint32_t i=0; i<(uint32_t)epi::LocationType::Count; i++){
        ASSERT_THAT(world.get_locations()[i], testing::ElementsAre());
    }
    ASSERT_THAT(world.get_persons(), testing::ElementsAre());
}

TEST(TestWorld, addLocation)
{
    auto world   = epi::World();
    auto school_id1 = world.add_location(epi::LocationType::School);
    auto school_id2 = world.add_location(epi::LocationType::School);
    auto work_id   = world.add_location(epi::LocationType::Work);
    auto home_id   = world.add_location(epi::LocationType::Home);

    ASSERT_EQ(school_id1.index, 0);
    ASSERT_EQ(school_id2.index, 1);

    auto& school1 = world.get_individualized_location(school_id1);
    auto& school2 = world.get_individualized_location(school_id2);
    auto& work = world.get_individualized_location(work_id);
    auto& home = world.get_individualized_location(home_id);

    ASSERT_EQ(world.get_locations().size(), (uint32_t)epi::LocationType::Count);
    ASSERT_EQ(world.get_locations()[(uint32_t)epi::LocationType::School].size(), 2);
    
    ASSERT_EQ(&world.get_locations()[(uint32_t)epi::LocationType::School][0], &school1);
    ASSERT_EQ(&world.get_locations()[(uint32_t)epi::LocationType::School][1], &school2);
    ASSERT_EQ(&world.get_locations()[(uint32_t)epi::LocationType::Work][0], &work);
    ASSERT_EQ(&world.get_locations()[(uint32_t)epi::LocationType::Home][0], &home);
}

TEST(TestWorld, addPerson)
{
    auto world = epi::World();
    auto location = world.add_location(epi::LocationType::School);

    auto& p1 = world.add_person(location, epi::InfectionState::Recovered_Carrier);
    auto& p2 = world.add_person(location, epi::InfectionState::Exposed);

    ASSERT_EQ(world.get_persons().size(), 2);
    ASSERT_EQ(&world.get_persons()[0], &p1);
    ASSERT_EQ(&world.get_persons()[1], &p2);
}

TEST(TestWorld, getSubpopulationCombined)
{
    auto world = epi::World();
    auto school1 = world.add_location(epi::LocationType::School);
    auto school2 = world.add_location(epi::LocationType::School);
    auto school3 = world.add_location(epi::LocationType::School);
    world.add_person(school1, epi::InfectionState::Carrier);
    world.add_person(school1, epi::InfectionState::Susceptible);
    world.add_person(school2, epi::InfectionState::Susceptible);
    world.add_person(school2, epi::InfectionState::Susceptible);
    world.add_person(school3, epi::InfectionState::Carrier);

    ASSERT_EQ(world.get_subpopulation_combined(epi::InfectionState::Susceptible, epi::LocationType::School), 3);
    ASSERT_EQ(world.get_subpopulation_combined(epi::InfectionState::Carrier, epi::LocationType::School), 2);
}

TEST(TestWorld, evolveStateTransition)
{
    using testing::Return;

    auto world      = epi::World();
    auto location1 = world.add_location(epi::LocationType::School);
    auto &p1 = world.add_person(location1, epi::InfectionState::Carrier);
    auto &p2 = world.add_person(location1, epi::InfectionState::Susceptible);
    auto location2 = world.add_location(epi::LocationType::Work);
    auto &p3 = world.add_person(location2, epi::InfectionState::Infected);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);


    //setup mock so only p2 transitions
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::AtLeast(3))
        .WillOnce(Return(0.51))
        .WillOnce(Return(0.04))
        .WillOnce(Return(0.6))
        .WillRepeatedly(Return(1.0));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    world.evolve(epi::TimePoint(0), epi::hours(1));

    EXPECT_EQ(p1.get_infection_state(), epi::InfectionState::Carrier);
    EXPECT_EQ(p2.get_infection_state(), epi::InfectionState::Exposed);
    EXPECT_EQ(p3.get_infection_state(), epi::InfectionState::Infected);
}

TEST(TestMigrationRules, school)
{
    auto home    = epi::Location(epi::LocationType::Home, 0);
    auto p_child = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age5to14, {}, {});
    auto p_adult = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34, {}, {});

    auto t_morning  = epi::TimePoint(0) + epi::hours(8);
    auto t_weekend = epi::TimePoint(0) + epi::days(5) + epi::hours(8);
    auto dt = epi::hours(1);

    ASSERT_EQ(epi::go_to_school(p_child, t_morning, dt, {}), epi::LocationType::School);
    ASSERT_EQ(epi::go_to_school(p_adult, t_morning, dt, {}), epi::LocationType::Home);
    ASSERT_EQ(epi::go_to_school(p_child, t_weekend, dt, {}), epi::LocationType::Home);
}

TEST(TestMigrationRules, school_return)
{
    auto school  = epi::Location(epi::LocationType::School, 0);
    auto p_child = epi::Person(school, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age5to14, {}, {});

    auto t  = epi::TimePoint(0) + epi::hours(15);
    auto dt = epi::hours(1);

    ASSERT_EQ(epi::go_to_school(p_child, t, dt, {}), epi::LocationType::Home);
}

TEST(TestMigrationRules, work)
{
    auto home      = epi::Location(epi::LocationType::Home, 0);
    auto p_retiree = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age60to79, {}, {});
    auto p_adult   = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34, {}, {});

    auto t_morning  = epi::TimePoint(0) + epi::hours(8);
    auto t_night = epi::TimePoint(0) + epi::days(1) + epi::hours(4);
    auto dt = epi::hours(1);

    ASSERT_EQ(epi::go_to_work(p_retiree, t_morning, dt, {}), epi::LocationType::Home);
    ASSERT_EQ(epi::go_to_work(p_adult, t_morning, dt, {}), epi::LocationType::Work);
    ASSERT_EQ(epi::go_to_work(p_adult, t_night, dt, {}), epi::LocationType::Home);
}

TEST(TestMigrationRules, work_return)
{
    auto work    = epi::Location(epi::LocationType::Work, 0);
    auto p_adult = epi::Person(work, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age35to59, {}, {});

    auto t  = epi::TimePoint(0) + epi::hours(17);
    auto dt = epi::hours(1);

    ASSERT_EQ(epi::go_to_work(p_adult, t, dt, {}), epi::LocationType::Home);
}

TEST(TestMigrationRules, hospital)
{
    auto home  = epi::Location(epi::LocationType::Home, 0);
    auto p_inf = epi::Person(home, epi::InfectionState::Infected_Severe, epi::AbmAgeGroup::Age15to34, {}, {});
    auto t     = epi::TimePoint(12346);
    auto dt    = epi::hours(1);

    ASSERT_EQ(epi::go_to_hospital(p_inf, t, dt, {}), epi::LocationType::Hospital);

    auto p_car = epi::Person(home, epi::InfectionState::Infected, epi::AbmAgeGroup::Age15to34, {}, {});
    ASSERT_EQ(epi::go_to_hospital(p_car, t, dt, {}), epi::LocationType::Home);
}

TEST(TestMigrationRules, icu)
{
    auto hospital = epi::Location(epi::LocationType::Hospital, 0);
    auto p_hosp   = epi::Person(hospital, epi::InfectionState::Infected_Critical, epi::AbmAgeGroup::Age15to34, {}, {});
    auto t        = epi::TimePoint(12346);
    auto dt       = epi::hours(1);

    ASSERT_EQ(epi::go_to_icu(p_hosp, t, dt, {}), epi::LocationType::ICU);

    auto work   = epi::Location(epi::LocationType::Work, 0);
    auto p_work = epi::Person(work, epi::InfectionState::Infected, epi::AbmAgeGroup::Age15to34, {}, {});
    ASSERT_EQ(epi::go_to_icu(p_work, t, dt, {}), epi::LocationType::Work);
}

TEST(TestMigrationRules, recover)
{
    auto hospital = epi::Location(epi::LocationType::Hospital, 0);
    auto p_rec    = epi::Person(hospital, epi::InfectionState::Recovered_Infected, epi::AbmAgeGroup::Age60to79, {}, {});
    auto p_inf    = epi::Person(hospital, epi::InfectionState::Infected_Severe, epi::AbmAgeGroup::Age60to79, {}, {});
    auto t        = epi::TimePoint(12346);
    auto dt       = epi::hours(1);

    ASSERT_EQ(epi::return_home_when_recovered(p_rec, t, dt, {}), epi::LocationType::Home);
    ASSERT_EQ(epi::return_home_when_recovered(p_inf, t, dt, {}), epi::LocationType::Hospital);
}

TEST(TestMigrationRules, go_shopping)
{
    auto hospital = epi::Location(epi::LocationType::Hospital, 0);
    auto p_hosp   = epi::Person(hospital, epi::InfectionState::Infected, epi::AbmAgeGroup::Age0to4, {}, {});
    auto home     = epi::Location(epi::LocationType::Home, 0);
    auto p_home   = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age60to79, {}, {});

    auto t_weekday = epi::TimePoint(0) + epi::days(4) + epi::hours(9);
    auto t_sunday  = epi::TimePoint(0) + epi::days(6) + epi::hours(9);
    auto t_night   = epi::TimePoint(0) + epi::days(4) + epi::hours(1);
    auto dt        = epi::hours(1);

    ASSERT_EQ(epi::go_to_shop(p_hosp, t_weekday, dt, {}), epi::LocationType::Hospital);
    ASSERT_EQ(epi::go_to_shop(p_home, t_sunday, dt, {}), epi::LocationType::Home);
    ASSERT_EQ(epi::go_to_shop(p_home, t_night, dt, {}), epi::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(epi::go_to_shop(p_home, t_weekday, dt, {}), epi::LocationType::BasicsShop);
}

TEST(TestMigrationRules, shop_return)
{
    auto t  = epi::TimePoint(0) + epi::days(4) + epi::hours(9);
    auto dt = epi::hours(1);

    auto home = epi::Location(epi::LocationType::Home, 0);
    auto shop = epi::Location(epi::LocationType::BasicsShop, 0);
    auto p    = epi::Person(home, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age15to34, {}, {});
    home.add_person(p);
    p.migrate_to(home, shop);
    p.interact(dt, {}, shop, {}); //person only returns home after some time passed

    ASSERT_EQ(epi::go_to_shop(p, t, dt, {}), epi::LocationType::Home);
}

TEST(TestMigrationRules, go_event)
{
    auto work   = epi::Location(epi::LocationType::Work, 0);
    auto p_work = epi::Person(work, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age35to59, {}, {});
    auto home   = epi::Location(epi::LocationType::Home, 0);
    auto p_home = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age60to79, {}, {});

    auto t_weekday  = epi::TimePoint(0) + epi::days(4) + epi::hours(20);
    auto t_saturday = epi::TimePoint(0) + epi::days(5) + epi::hours(10);
    auto t_night    = epi::TimePoint(0) + epi::days(5) + epi::hours(1);
    auto dt         = epi::hours(1);

    ASSERT_EQ(epi::go_to_event(p_work, t_weekday, dt, {}), epi::LocationType::Work);
    ASSERT_EQ(epi::go_to_event(p_home, t_night, dt, {}), epi::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(epi::go_to_event(p_home, t_weekday, dt, {}), epi::LocationType::SocialEvent);

    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(epi::go_to_event(p_home, t_saturday, dt, {}), epi::LocationType::SocialEvent);
}

TEST(TestMigrationRules, event_return)
{
    auto t  = epi::TimePoint(0) + epi::days(4) + epi::hours(21);
    auto dt = epi::hours(3);

    auto home = epi::Location(epi::LocationType::Home, 0);
    auto shop = epi::Location(epi::LocationType::SocialEvent, 0);
    auto p    = epi::Person(home, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age15to34, {}, {});
    home.add_person(p);
    p.migrate_to(home, shop);
    p.interact(dt, {}, shop, {});

    ASSERT_EQ(epi::go_to_event(p, t, dt, {}), epi::LocationType::Home);
}

TEST(TestLockdownRules, school_closure)
{
    auto t  = epi::TimePoint(0);
    auto dt = epi::hours(1);
    auto t_morning  = epi::TimePoint(0) + epi::hours(8);
    auto home = epi::Location(epi::LocationType::Home, 0);
    auto school = epi::Location(epi::LocationType::School, 0);

    //setup rng mock so one person is home schooled and the other goes to school
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformDistribution<double>>>>
        mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(4))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.2))
        .WillOnce(testing::Return(0.2))
        .WillRepeatedly(testing::Return(1.0)); 

    auto p1    = epi::Person(home, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age5to14, {}, {});
    p1.set_assigned_location(home);
    p1.set_assigned_location(school);
    auto p2    = epi::Person(home, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age5to14, {}, {});
    p2.set_assigned_location(home);
    p2.set_assigned_location(school);
    epi::AbmMigrationParameters params;

    epi::set_school_closure(t, 0.7, params);

    ASSERT_EQ(epi::go_to_school(p1, t_morning, dt, params), epi::LocationType::Home);
    ASSERT_EQ(epi::go_to_school(p2, t_morning, dt, params), epi::LocationType::School);
}

TEST(TestLockdownRules, school_opening)
{
    auto t_closing  = epi::TimePoint(0);
    auto t_opening  = epi::TimePoint(0) + epi::days(1);
    auto dt = epi::hours(1);
    auto t_morning  = epi::TimePoint(0) + epi::days(1) + epi::hours(8);
    auto home = epi::Location(epi::LocationType::Home, 0);
    auto school = epi::Location(epi::LocationType::School, 0);
    //setup rng mock so the person is homeschooled in case of lockdown
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformDistribution<double>>>>
        mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.5))
        .WillRepeatedly(testing::Return(1.0)); 
    auto p    = epi::Person(home, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age5to14, {}, {});
    p.set_assigned_location(home);
    p.set_assigned_location(school);
    epi::AbmMigrationParameters params;

    epi::set_school_closure(t_closing, 1., params);
    epi::set_school_closure(t_opening, 0., params);

    ASSERT_EQ(epi::go_to_school(p, t_morning, dt, params), epi::LocationType::School);
}

TEST(TestLockdownRules, home_office)
{
    auto t = epi::TimePoint(0);
    auto t_morning  = epi::TimePoint(0) + epi::hours(8);
    auto dt = epi::hours(1);
    auto home      = epi::Location(epi::LocationType::Home, 0);
    auto work      = epi::Location(epi::LocationType::Work, 0);
    epi::AbmMigrationParameters params;

    epi::set_home_office(t, 0.4, params);
    
    //setup rng mock so one person goes to work and the other works at home
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformDistribution<double>>>>
        mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(4))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillRepeatedly(testing::Return(1.0));        

    auto person1 = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34, {}, {});
    auto person2 = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34, {}, {});
    person1.set_assigned_location(home);
    person1.set_assigned_location(work);
    person2.set_assigned_location(home);
    person2.set_assigned_location(work);

    ASSERT_EQ(epi::go_to_work(person1, t_morning, dt, params), epi::LocationType::Work);
    ASSERT_EQ(epi::go_to_work(person2, t_morning, dt, params), epi::LocationType::Home);
    
}

TEST(TestLockdownRules, no_home_office)
{
    auto t_closing  = epi::TimePoint(0);
    auto t_opening  = epi::TimePoint(0) + epi::days(1);
    auto dt = epi::hours(1);
    auto t_morning  = epi::TimePoint(0) + epi::days(1) + epi::hours(8);
    auto home = epi::Location(epi::LocationType::Home, 0);
    auto work = epi::Location(epi::LocationType::Work, 0);


    //setup rng mock so the person works in home office
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformDistribution<double>>>>
        mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillRepeatedly(testing::Return(1.0)); 

    auto p    = epi::Person(home, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age15to34, {}, {});
    p.set_assigned_location(home);
    p.set_assigned_location(work);
    epi::AbmMigrationParameters params;

    epi::set_home_office(t_closing, 0.5, params);
    epi::set_home_office(t_opening, 0., params);

    ASSERT_EQ(epi::go_to_work(p, t_morning, dt, params), epi::LocationType::Work);
}

TEST(TestLockdownRules, social_event_closure)
{
    auto t  = epi::TimePoint(0);
    auto dt = epi::hours(1);
    auto t_evening  = epi::TimePoint(0) + epi::hours(19);
    auto home = epi::Location(epi::LocationType::Home, 0);
    auto event = epi::Location(epi::LocationType::SocialEvent, 0);
    auto p    = epi::Person(home, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age5to14, {}, {});
    p.set_assigned_location(home);
    p.set_assigned_location(event);
    epi::AbmMigrationParameters params;

    epi::close_social_events(t, 1, params);

    ASSERT_EQ(epi::go_to_event(p, t_evening, dt, params), epi::LocationType::Home);
}

TEST(TestLockdownRules, social_events_opening)
{
    auto t_closing  = epi::TimePoint(0);
    auto t_opening  = epi::TimePoint(0) + epi::days(1);
    auto dt = epi::hours(1);
    auto t_evening  = epi::TimePoint(0) + epi::days(1) + epi::hours(19);
    auto home = epi::Location(epi::LocationType::Home, 0);
    auto event = epi::Location(epi::LocationType::SocialEvent, 0);
    auto p    = epi::Person(home, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age5to14, {}, {});
    p.set_assigned_location(event);
    p.set_assigned_location(home);
    epi::AbmMigrationParameters params;

    epi::close_social_events(t_closing, 1, params);
    epi::close_social_events(t_opening, 0, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(epi::go_to_event(p, t_evening, dt, params), epi::LocationType::SocialEvent);
}

TEST(TestWorld, evolveMigration)
{
    using testing::Return;

    auto world      = epi::World();
    auto home_id = world.add_location(epi::LocationType::Home);
    auto school_id = world.add_location(epi::LocationType::School);
    auto work_id = world.add_location(epi::LocationType::Work);
    auto &p1 = world.add_person(home_id, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age15to34);
    auto &p2 = world.add_person(home_id, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age5to14);
    p1.set_assigned_location(school_id);
    p2.set_assigned_location(school_id);
    p1.set_assigned_location(work_id);
    p2.set_assigned_location(work_id);
    p1.set_assigned_location(home_id);
    p2.set_assigned_location(home_id);

    auto& school = world.get_individualized_location(school_id);
    auto& work = world.get_individualized_location(work_id);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

    world.evolve(epi::TimePoint(0) + epi::hours(8), epi::hours(1));

    EXPECT_EQ(p1.get_location_id().type, epi::LocationType::Work);
    EXPECT_EQ(p2.get_location_id().type, epi::LocationType::School);
    EXPECT_EQ(school.get_subpopulations().sum(), 1);
    EXPECT_EQ(work.get_subpopulations().sum(), 1);
}

TEST(TestSimulation, advance_random)
{
    auto world      = epi::World();
    auto location1 = world.add_location(epi::LocationType::School);
    auto location2 = world.add_location(epi::LocationType::School);
    auto &p1 = world.add_person(location1, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age5to14);
    auto &p2 = world.add_person(location1, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age5to14);
    auto &p3 = world.add_person(location2, epi::InfectionState::Infected, epi::AbmAgeGroup::Age5to14);
    auto &p4 = world.add_person(location2, epi::InfectionState::Infected, epi::AbmAgeGroup::Age5to14);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);
    p4.set_assigned_location(location2);

    auto sim = epi::AbmSimulation(epi::TimePoint(0), std::move(world));

    sim.advance(epi::TimePoint(0) + epi::hours(50));
    ASSERT_EQ(sim.get_result().get_num_time_points(), 51);
    ASSERT_THAT(sim.get_result().get_times(), ElementsAreLinspace(0.0, 50.0 / 24.0, 51));
    for (auto&& v : sim.get_result()) {
        ASSERT_EQ(v.sum(), 4);
    }
}

TEST(TestDiscreteDistribution, generate)
{
    using namespace epi;
    auto distribution = epi::DiscreteDistribution<size_t>();

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
