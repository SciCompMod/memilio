#include "epidemiology/abm/abm.h"
#include "epidemiology/abm/migration_rules.h"
#include "epidemiology/utils/eigen_util.h"
#include "matchers.h"
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
    ASSERT_EQ(print_wrap(location.get_subpopulations()),
              print_wrap(Eigen::VectorXi::Zero(Eigen::Index(epi::InfectionState::Count))));
}

TEST(TestLocation, addRemovePerson)
{
    auto location = epi::Location(epi::LocationType::Home);
    auto person1  = epi::Person(location, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age5to14);
    location.add_person(person1);
    auto person2 = epi::Person(location, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34);
    location.add_person(person2);
    auto person3 = epi::Person(location, epi::InfectionState::Exposed, epi::AbmAgeGroup::Age35to59);
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
    auto person   = epi::Person(location, epi::InfectionState::Recovered_Carrier, epi::AbmAgeGroup::Age60to79);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Recovered_Carrier);
    ASSERT_EQ(&person.get_location(), &location);
    ASSERT_EQ(person.get_age(), epi::AbmAgeGroup::Age60to79);
}

TEST(TestPerson, migrate)
{
    auto location1 = epi::Location(epi::LocationType::Work);
    auto location2 = epi::Location(epi::LocationType::School);
    auto person    = epi::Person(location1, epi::InfectionState::Recovered_Carrier, epi::AbmAgeGroup::Age0to4);
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
    auto infected1 = epi::Person(location, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age15to34);
    location.add_person(infected1);
    auto infected2 = epi::Person(location, epi::InfectionState::Infected_Detected, epi::AbmAgeGroup::Age80plus);
    location.add_person(infected2);
    auto infected3 = epi::Person(location, epi::InfectionState::Infected_Undetected, epi::AbmAgeGroup::Age5to14);
    location.add_person(infected3);

    //test should work identically work with any age
    epi::AbmAgeGroup age = epi::AbmAgeGroup(epi::UniformIntDistribution<int>()(0, int(epi::AbmAgeGroup::Count) - 1));
    epi::GlobalInfectionParameters params;
    params.carrier_to_infected = 0;
    params.carrier_to_infected[age] = 0.5;
    params.carrier_to_recovered = 0;
    params.carrier_to_recovered[age] = 0.5;
    params.detect_infection = 0;
    params.detect_infection[age] = 0.5;
    params.infected_to_dead = 0;
    params.infected_to_dead[age] = 0.5;
    params.infected_to_recovered = 0;
    params.infected_to_recovered[age] = 0.5;
    params.recovered_to_susceptible = 0;
    params.recovered_to_susceptible[age] = 0.5;
    params.susceptible_to_exposed_by_carrier = 0;
    params.susceptible_to_exposed_by_carrier[age] = 0.5;
    params.susceptible_to_exposed_by_infected = 0;
    params.susceptible_to_exposed_by_infected[age] = 0.5;

    //cache precomputed results
    auto dt = epi::seconds(8640); //0.1 days
    location.begin_step(dt, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    {
        auto susceptible = epi::Person(location, epi::InfectionState::Susceptible, age);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(susceptible, dt, params), epi::InfectionState::Exposed);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.15));
        EXPECT_EQ(location.interact(susceptible, dt, params), epi::InfectionState::Susceptible);
    }

    {
        auto exposed = epi::Person(location, epi::InfectionState::Exposed, age);
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(0); //no transitions out of exposed state
        EXPECT_EQ(location.interact(exposed, dt, params), epi::InfectionState::Exposed);
    }

    {
        auto carrier = epi::Person(location, epi::InfectionState::Carrier, age);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.05));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(carrier, dt, params), epi::InfectionState::Infected_Detected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(carrier, dt, params), epi::InfectionState::Infected_Undetected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.099));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(2));
        EXPECT_EQ(location.interact(carrier, dt, params), epi::InfectionState::Recovered_Carrier);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(carrier, dt, {}), epi::InfectionState::Carrier);
    }

    for (auto&& infected_state : {epi::InfectionState::Infected_Detected, epi::InfectionState::Infected_Undetected}) {
        auto infected = epi::Person(location, infected_state, age);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(infected, dt, params), epi::InfectionState::Recovered_Infected);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(1));
        EXPECT_EQ(location.interact(infected, dt, params), epi::InfectionState::Dead);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.1001));
        EXPECT_EQ(location.interact(infected, dt, params), infected_state);
    }

    for (auto&& recovered_state : {epi::InfectionState::Recovered_Carrier, epi::InfectionState::Recovered_Infected}) {
        auto recovered = epi::Person(location, recovered_state, age);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
        EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
        EXPECT_EQ(location.interact(recovered, dt, params), epi::InfectionState::Susceptible);

        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.11));
        EXPECT_EQ(location.interact(recovered, dt, params), recovered_state);
    }
}

TEST(TestPerson, interact)
{
    using testing::Return;

    auto location = epi::Location(epi::LocationType::Home);
    auto person   = epi::Person(location, epi::InfectionState::Infected_Detected, epi::AbmAgeGroup::Age15to34);
    location.add_person(person);
    auto dt = epi::seconds(8640); //0.1 days
    location.begin_step(dt, {});

    //setup rng mock so the person has a state transition
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    auto infection_parameters = epi::GlobalInfectionParameters();
    person.interact(dt, infection_parameters);
    EXPECT_EQ(person.get_infection_state(), epi::InfectionState::Recovered_Infected);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Recovered_Infected), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Detected), 0);
}

TEST(TestPerson, interact_exposed)
{
    using testing::Return;

    //setup location with some chance of exposure
    auto location  = epi::Location(epi::LocationType::Work);
    auto infected1 = epi::Person(location, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age15to34);
    location.add_person(infected1);
    auto infected2 = epi::Person(location, epi::InfectionState::Infected_Detected, epi::AbmAgeGroup::Age5to14);
    location.add_person(infected2);
    auto infected3 = epi::Person(location, epi::InfectionState::Infected_Undetected, epi::AbmAgeGroup::Age60to79);
    location.add_person(infected3);
    auto person = epi::Person(location, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34);
    location.add_person(person);
    location.begin_step(epi::hours(1), {});

    auto infection_parameters              = epi::GlobalInfectionParameters();
    infection_parameters.incubation_period = 2;

    //setup rng mock so the person becomes exposed
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.49));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    //person becomes exposed
    person.interact(epi::hours(12), infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Exposed), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Carrier), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Detected), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Undetected), 1);

    //person becomes a carrier after the incubation time runs out, not random
    person.interact(epi::hours(12), infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);

    person.interact(epi::hours(12), infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);

    person.interact(epi::hours(24), infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);

    person.interact(epi::hours(1), infection_parameters);
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

TEST(TestWorld, evolveStateTransition)
{
    using testing::Return;

    auto world      = epi::World();
    auto& location1 = world.add_location(epi::LocationType::School);
    world.add_person(location1, epi::InfectionState::Carrier);
    world.add_person(location1, epi::InfectionState::Susceptible);
    auto& location2 = world.add_location(epi::LocationType::Work);
    world.add_person(location2, epi::InfectionState::Infected_Detected);

    //setup mock so only p2 transitions
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(3)
        .WillOnce(Return(0.51))
        .WillOnce(Return(0.04))
        .WillOnce(Return(0.6));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    world.evolve(epi::TimePoint(0), epi::hours(1));

    EXPECT_EQ(world.get_persons()[0].get_infection_state(), epi::InfectionState::Carrier);
    EXPECT_EQ(world.get_persons()[1].get_infection_state(), epi::InfectionState::Exposed);
    EXPECT_EQ(world.get_persons()[2].get_infection_state(), epi::InfectionState::Infected_Detected);
}

TEST(TestMigrationRules, school)
{
    auto home = epi::Location(epi::LocationType::Home);
    auto p_child = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age5to14);
    auto p_adult = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34);

    auto t = epi::TimePoint(0) + epi::hours(8);
    auto dt = epi::hours(1);

    ASSERT_EQ(epi::go_to_school(p_child, t, dt, {}), epi::LocationType::School);
    ASSERT_EQ(epi::go_to_school(p_adult, t, dt, {}), epi::LocationType::Home);    
}

TEST(TestMigrationRules, school_return)
{
    auto school = epi::Location(epi::LocationType::School);
    auto p_child = epi::Person(school, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age5to14);

    auto t = epi::TimePoint(0) + epi::hours(15);
    auto dt = epi::hours(1);

    ASSERT_EQ(epi::go_to_school(p_child, t, dt, {}), epi::LocationType::Home);
}

TEST(TestMigrationRules, work)
{
    auto home = epi::Location(epi::LocationType::Home);
    auto p_retiree = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age60to79);
    auto p_adult = epi::Person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age15to34);

    auto t = epi::TimePoint(0) + epi::hours(8);
    auto dt = epi::hours(1);

    ASSERT_EQ(epi::go_to_work(p_retiree, t, dt, {}), epi::LocationType::Home);
    ASSERT_EQ(epi::go_to_work(p_adult, t, dt, {}), epi::LocationType::Work);
}

TEST(TestMigrationRules, work_return)
{
    auto work = epi::Location(epi::LocationType::Work);
    auto p_adult = epi::Person(work, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age35to59);

    auto t = epi::TimePoint(0) + epi::hours(17);
    auto dt = epi::hours(1);

    ASSERT_EQ(epi::go_to_work(p_adult, t, dt, {}), epi::LocationType::Home);
}

TEST(TestWorld, evolveMigration)
{
    using testing::Return;

    auto world      = epi::World();
    auto& home = world.add_location(epi::LocationType::Home);
    world.add_person(home, epi::InfectionState::Carrier, epi::AbmAgeGroup::Age15to34);
    world.add_person(home, epi::InfectionState::Susceptible, epi::AbmAgeGroup::Age5to14);
    auto& school = world.add_location(epi::LocationType::School);
    auto& work = world.add_location(epi::LocationType::Work);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<epi::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

    world.evolve(epi::TimePoint(0) + epi::hours(8), epi::hours(1));

    EXPECT_EQ(world.get_persons()[0].get_location().get_type(), epi::LocationType::Work);
    EXPECT_EQ(world.get_persons()[1].get_location().get_type(), epi::LocationType::School);
    EXPECT_EQ(school.get_subpopulations().sum(), 1);
    EXPECT_EQ(work.get_subpopulations().sum(), 1);
}

TEST(TestSimulation, advance_random)
{
    auto world      = epi::World();
    auto& location1 = world.add_location(epi::LocationType::School);
    world.add_person(location1, epi::InfectionState::Carrier);
    world.add_person(location1, epi::InfectionState::Susceptible);
    auto& location2 = world.add_location(epi::LocationType::School);
    world.add_person(location2, epi::InfectionState::Infected_Detected);
    world.add_person(location2, epi::InfectionState::Infected_Undetected);

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

TEST(TestDependentParameter, init)
{
    enum class E
    {
        E1 = 0, E2, Count,
    };
    epi::DependentParameter<E> p1(3.0);
    ASSERT_EQ(p1[E::E1], 3.0);
    ASSERT_EQ(p1[E::E2], 3.0);

    auto v = (Eigen::ArrayXd(2) << 1.0, 2.0).finished();
    epi::DependentParameter<E> p2(v);
    ASSERT_EQ(p2[E::E1], 1.0);
    ASSERT_EQ(p2[E::E2], 2.0);
    ASSERT_THAT(print_wrap(p2.array().matrix()), print_wrap(v.matrix()));
}

TEST(TestDependentParameter, assign)
{
    enum class E
    {
        E1 = 0, E2, Count,
    };
    epi::DependentParameter<E> p;

    p = 5.0;
    ASSERT_EQ(p[E::E1], 5.0);
    ASSERT_EQ(p[E::E2], 5.0);

    auto v = (Eigen::ArrayXd(2) << 1.0, 2.0).finished();
    p = v;
    ASSERT_EQ(p[E::E1], 1.0);
    ASSERT_EQ(p[E::E2], 2.0);
}