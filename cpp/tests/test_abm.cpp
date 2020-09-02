#include "epidemiology/abm/abm.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <memory>

TEST(TestLocation, init)
{
    auto location = epi::Location(epi::LocationType::School);
    for (epi::InfectionState i = epi::InfectionState(0); i < epi::InfectionState::Count;
         i                     = epi::InfectionState(size_t(i) + 1)) {
        ASSERT_EQ(location.get_subpopulation(i), 0);
    }
    ASSERT_EQ(epi::print_wrap(location.get_subpopulations()),
              epi::print_wrap(Eigen::VectorXi::Zero(Eigen::Index(epi::InfectionState::Count))));
}

TEST(TestLocation, addRemovePerson)
{
    auto location = epi::Location(epi::LocationType::Home);
    auto person1  = epi::Person(location, epi::InfectionState::Susceptible);
    location.add_person(person1);
    auto person2 = epi::Person(location, epi::InfectionState::Susceptible);
    location.add_person(person2);
    auto person3 = epi::Person(location, epi::InfectionState::Exposed);
    location.add_person(person3);

    ASSERT_EQ(location.get_subpopulation(epi::InfectionState::Susceptible), 2);
    ASSERT_EQ(location.get_subpopulation(epi::InfectionState::Exposed), 1);

    location.remove_person(person2);

    ASSERT_EQ(location.get_subpopulation(epi::InfectionState::Susceptible), 1);
    ASSERT_EQ(location.get_subpopulation(epi::InfectionState::Exposed), 1);
}

TEST(TestPerson, init)
{
    auto location = epi::Location(epi::LocationType::Work);
    auto person   = epi::Person(location, epi::InfectionState::Recovered_Carrier);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Recovered_Carrier);
    ASSERT_EQ(&person.get_location(), &location);
}

TEST(TestPerson, migrate)
{
    auto location1 = epi::Location(epi::LocationType::Work);
    auto location2 = epi::Location(epi::LocationType::School);
    auto person    = epi::Person(location1, epi::InfectionState::Recovered_Carrier);
    location1.add_person(person);

    person.migrate_to(location2);

    ASSERT_EQ(&person.get_location(), &location2);
    ASSERT_EQ(location2.get_subpopulation(epi::InfectionState::Recovered_Carrier), 1);
    ASSERT_EQ(location1.get_subpopulation(epi::InfectionState::Recovered_Carrier), 0);
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

TEST(TestLocation, interact)
{
    using testing::Return;

    //setup location with some chance of exposure
    auto location  = epi::Location(epi::LocationType::Work);
    auto infected1 = epi::Person(location, epi::InfectionState::Carrier);
    location.add_person(infected1);
    auto infected2 = epi::Person(location, epi::InfectionState::Infected_Detected);
    location.add_person(infected2);
    auto infected3 = epi::Person(location, epi::InfectionState::Infected_Undetected);
    location.add_person(infected3);
    location.begin_step(0.1, {});

    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformDistribution<double>>>> mock_uniform_dist;

    {
        auto susceptible = epi::Person(location, epi::InfectionState::Susceptible);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.5));
        EXPECT_EQ(location.interact(susceptible, 0.1, {}), epi::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.15));
        EXPECT_EQ(location.interact(susceptible, 0.1, {}), epi::InfectionState::Susceptible);
    }

    {
        auto exposed = epi::Person(location, epi::InfectionState::Exposed);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(0); //no transitions out of exposed state
        EXPECT_EQ(location.interact(exposed, 0.1, {}), epi::InfectionState::Exposed);
    }

    {
        auto carrier = epi::Person(location, epi::InfectionState::Carrier);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1));
        EXPECT_EQ(location.interact(carrier, 0.1, {}), epi::InfectionState::Infected_Detected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.4));
        EXPECT_EQ(location.interact(carrier, 0.1, {}), epi::InfectionState::Infected_Undetected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.099));
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.7));
        EXPECT_EQ(location.interact(carrier, 0.1, {}), epi::InfectionState::Recovered_Carrier);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(carrier, 0.1, {}), epi::InfectionState::Carrier);
    }

    for (auto&& infected_state : {epi::InfectionState::Infected_Detected, epi::InfectionState::Infected_Undetected}) {
        auto infected = epi::Person(location, infected_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.3));
        EXPECT_EQ(location.interact(infected, 0.1, {}), epi::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.6));
        EXPECT_EQ(location.interact(infected, 0.1, {}), epi::InfectionState::Dead);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(infected, 0.1, {}), infected_state);
    }

    for (auto&& recovered_state : {epi::InfectionState::Recovered_Carrier, epi::InfectionState::Recovered_Infected}) {
        auto recovered = epi::Person(location, recovered_state);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.5));
        EXPECT_EQ(location.interact(recovered, 0.1, {}), epi::InfectionState::Susceptible);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(recovered, 0.1, {}), recovered_state);
    }
}

TEST(TestPerson, interact)
{
    using testing::Return;

    auto location = epi::Location(epi::LocationType::Home);
    auto person   = epi::Person(location, epi::InfectionState::Infected_Detected);
    location.add_person(person);
    location.begin_step(0.1, {});

    //setup rng mock so the person has a state transition
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.3));

    auto infection_parameters = epi::GlobalInfectionParameters();
    person.interact(0.1, infection_parameters);
    EXPECT_EQ(person.get_infection_state(), epi::InfectionState::Recovered_Infected);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Recovered_Infected), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Detected), 0);
}

TEST(TestPerson, interact_exposed)
{
    using testing::Return;

    //setup location with some chance of exposure
    auto location  = epi::Location(epi::LocationType::Work);
    auto infected1 = epi::Person(location, epi::InfectionState::Carrier);
    location.add_person(infected1);
    auto infected2 = epi::Person(location, epi::InfectionState::Infected_Detected);
    location.add_person(infected2);
    auto infected3 = epi::Person(location, epi::InfectionState::Infected_Undetected);
    location.add_person(infected3);
    auto person = epi::Person(location, epi::InfectionState::Susceptible);
    location.add_person(person);
    location.begin_step(0.1, {});

    auto infection_parameters              = epi::GlobalInfectionParameters();
    infection_parameters.incubation_period = 2;

    //setup rng mock so the person becomes exposed
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.49));
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.5));

    //person becomes exposed
    person.interact(0.5, infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Exposed), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Carrier), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Detected), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Undetected), 1);

    //person becomes a carrier after the incubation time runs out, not random
    person.interact(0.5, infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);

    person.interact(0.5, infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);

    person.interact(1.0, infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);

    person.interact(0.1, infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Carrier);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Exposed), 0);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Carrier), 2);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Detected), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Undetected), 1);
}

TEST(TestWorld, init)
{
    auto world = epi::World();
    ASSERT_THAT(world.get_locations(), testing::ElementsAre());
    ASSERT_THAT(world.get_persons(), testing::ElementsAre());
}

TEST(TestWorld, addLocation)
{
    auto world   = epi::World();
    auto& school = world.add_location(epi::LocationType::School);
    auto& work   = world.add_location(epi::LocationType::Work);
    auto& home   = world.add_location(epi::LocationType::Home);

    ASSERT_EQ(world.get_locations().size(), 3);
    ASSERT_EQ(&world.get_locations()[0], &school);
    ASSERT_EQ(&world.get_locations()[1], &work);
    ASSERT_EQ(&world.get_locations()[2], &home);
}

TEST(TestWorld, addPerson)
{
    auto world     = epi::World();
    auto& location = world.add_location(epi::LocationType::School);

    auto& p1 = world.add_person(location, epi::InfectionState::Recovered_Carrier);
    auto& p2 = world.add_person(location, epi::InfectionState::Exposed);

    ASSERT_EQ(world.get_persons().size(), 2);
    ASSERT_EQ(&world.get_persons()[0], &p1);
    ASSERT_EQ(&world.get_persons()[1], &p2);
}

TEST(TestWorld, evolve)
{
    using testing::Return;

    auto world      = epi::World();
    auto& location1 = world.add_location(epi::LocationType::School);
    auto& p1        = world.add_person(location1, epi::InfectionState::Carrier);
    auto& p2        = world.add_person(location1, epi::InfectionState::Susceptible);
    auto& location2 = world.add_location(epi::LocationType::School);
    auto& p3        = world.add_person(location2, epi::InfectionState::Infected_Detected);

    //setup mock so only p2 transitions
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::AtLeast(3))
        .WillOnce(Return(0.51))
        .WillOnce(Return(0.4))
        .WillOnce(Return(0.6))
        .WillRepeatedly(Return(1e10)); //no random migration, not yet implemented properly
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.5));

    world.evolve(0.5);

    EXPECT_EQ(world.get_persons()[0].get_infection_state(), epi::InfectionState::Carrier);
    EXPECT_EQ(world.get_persons()[1].get_infection_state(), epi::InfectionState::Exposed);
    EXPECT_EQ(world.get_persons()[2].get_infection_state(), epi::InfectionState::Infected_Detected);
}

template <class T>
auto ElementsAreLinspace(T b, T e, size_t num_points)
{
    assert(num_points >= 2);

    std::vector<T> values(num_points);
    for (size_t i = 0; i < num_points; i++) {
        values[i] = i * (e - b) / (num_points - 1);
    }
    return testing::ElementsAreArray(values);
}

TEST(TestSimulation, advance_random)
{
    auto world      = epi::World();
    auto& location1 = world.add_location(epi::LocationType::School);
    auto& p1        = world.add_person(location1, epi::InfectionState::Carrier);
    auto& p2        = world.add_person(location1, epi::InfectionState::Susceptible);
    auto& location2 = world.add_location(epi::LocationType::School);
    auto& p3        = world.add_person(location2, epi::InfectionState::Infected_Detected);
    auto& p4        = world.add_person(location2, epi::InfectionState::Infected_Undetected);

    auto sim = epi::AbmSimulation(0, std::move(world));

    sim.advance(50);
    ASSERT_EQ(sim.get_result().get_num_time_points(), 51);
    ASSERT_THAT(sim.get_result().get_times(), ElementsAreLinspace(0.0, 50.0, 51));
    for (auto&& v : sim.get_result()) {
        ASSERT_EQ(v.sum(), 4);
    }
}
