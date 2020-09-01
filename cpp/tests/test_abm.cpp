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

//mocking a rng and using it in e.g. uniform_int_distribution is actually not allowed.
//the standard requires the generator to be actually random, which the mock is not.
//the resulting behaviour is undefined, standard compliant implementations may e.g. loop forever
//possible solution: mock the distribution instead of the generator
struct MockRng {
    MOCK_METHOD(epi::RandomNumberGenerator::result_type, invoke, (), ());
};
template <class MockRng = MockRng>
struct MockRngRef {
    epi::RandomNumberGenerator::result_type operator()()
    {
        return mock->invoke();
    }
    std::shared_ptr<MockRng> mock = std::make_shared<MockRng>();
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

    MockRngRef<testing::StrictMock<MockRng>> mock_rng;
    epi::thread_local_rng().generator = mock_rng;

    {
        auto susceptible = epi::Person(location, epi::InfectionState::Susceptible);
        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.99 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.5 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(susceptible, 0.1, {}), epi::InfectionState::Exposed);
        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.1 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.5 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(susceptible, 0.1, {}), epi::InfectionState::Susceptible);
    }

    {
        auto exposed = epi::Person(location, epi::InfectionState::Exposed);
        EXPECT_CALL(*mock_rng.mock, invoke).Times(0); //no transitions out of exposed state
        EXPECT_EQ(location.interact(exposed, 0.1, {}), epi::InfectionState::Exposed);
    }

    {
        auto carrier = epi::Person(location, epi::InfectionState::Carrier);

        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.99 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.2 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(carrier, 0.1, {}), epi::InfectionState::Infected_Detected);
        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.99 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.4 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(carrier, 0.1, {}), epi::InfectionState::Infected_Undetected);
        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.99 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.6 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(carrier, 0.1, {}), epi::InfectionState::Recovered_Carrier);
        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.5 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.5 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(carrier, 0.1, {}), epi::InfectionState::Carrier);
    }

    for (auto&& infected_state : {epi::InfectionState::Infected_Detected, epi::InfectionState::Infected_Undetected}) {
        auto infected = epi::Person(location, infected_state);

        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.99 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.4 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(infected, 0.1, {}), epi::InfectionState::Recovered_Infected);
        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.99 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.6 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(infected, 0.1, {}), epi::InfectionState::Dead);
        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.5 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.5 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(infected, 0.1, {}), infected_state);
    }

    for (auto&& recovered_state : {epi::InfectionState::Recovered_Carrier, epi::InfectionState::Recovered_Infected}) {
        auto recovered = epi::Person(location, recovered_state);

        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.99 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.5 * epi::RandomNumberGenerator::max()));
        EXPECT_EQ(location.interact(recovered, 0.1, {}), epi::InfectionState::Susceptible);
        EXPECT_CALL(*mock_rng.mock, invoke)
            .Times(2)
            .WillOnce(Return(0.5 * epi::RandomNumberGenerator::max()))
            .WillOnce(Return(0.5 * epi::RandomNumberGenerator::max()));
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
    MockRngRef<testing::StrictMock<MockRng>> mock_rng;
    epi::thread_local_rng().generator = mock_rng;
    EXPECT_CALL(*mock_rng.mock, invoke)
        .Times(2)
        .WillOnce(Return(0.99 * epi::RandomNumberGenerator::max()))
        .WillOnce(Return(0.4 * epi::RandomNumberGenerator::max()));

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
    MockRngRef<testing::StrictMock<MockRng>> mock_rng;
    epi::thread_local_rng().generator = mock_rng;
    EXPECT_CALL(*mock_rng.mock, invoke)
        .Times(2)
        .WillOnce(Return(0.99 * epi::RandomNumberGenerator::max()))
        .WillOnce(Return(0.4 * epi::RandomNumberGenerator::max()));

    //person becomes a carrier after the incubation time runs out
    person.interact(0.1, infection_parameters);
    ASSERT_EQ(person.get_infection_state(), epi::InfectionState::Exposed);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Exposed), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Carrier), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Detected), 1);
    EXPECT_EQ(location.get_subpopulation(epi::InfectionState::Infected_Undetected), 1);

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

    MockRngRef<testing::StrictMock<MockRng>> mock_rng;
    epi::thread_local_rng().generator = mock_rng;
    EXPECT_CALL(*mock_rng.mock, invoke)
        .Times(testing::AtLeast(6))
        .WillOnce(Return(0.1 * epi::RandomNumberGenerator::max())) //don't transition
        .WillOnce(Return(0.4 * epi::RandomNumberGenerator::max()))
        .WillOnce(Return(0.9 * epi::RandomNumberGenerator::max())) //transition to exposed
        .WillOnce(Return(0.4 * epi::RandomNumberGenerator::max()))
        .WillOnce(Return(0.1 * epi::RandomNumberGenerator::max())) //don't transition
        .WillOnce(Return(0.4 * epi::RandomNumberGenerator::max()))
        .WillRepeatedly(Return(
            epi::RandomNumberGenerator::max() -
            1)); //no migrations (not properly implemented yet); can't return max() because of a bug(?) in the STL

    world.evolve(0.5);

    EXPECT_EQ(world.get_persons()[0].get_infection_state(), epi::InfectionState::Carrier);
    EXPECT_EQ(world.get_persons()[1].get_infection_state(), epi::InfectionState::Exposed);
    EXPECT_EQ(world.get_persons()[2].get_infection_state(), epi::InfectionState::Infected_Detected);
}

TEST(TestSimulation, advance)
{
    auto world      = epi::World();
    auto& location1 = world.add_location(epi::LocationType::School);
    auto& p1        = world.add_person(location1, epi::InfectionState::Carrier);
    auto& p2        = world.add_person(location1, epi::InfectionState::Susceptible);
    auto& location2 = world.add_location(epi::LocationType::School);
    auto& p3        = world.add_person(location2, epi::InfectionState::Infected_Detected);
    auto& p4        = world.add_person(location2, epi::InfectionState::Infected_Undetected);

    auto sim = epi::AbmSimulation(0, std::move(world));

    MockRngRef<testing::StrictMock<MockRng>> mock_rng;
    epi::thread_local_rng().generator = mock_rng;
    EXPECT_CALL(*mock_rng.mock, invoke)
        .Times(testing::AtLeast(40))
        .WillRepeatedly(testing::Return(epi::RandomNumberGenerator::max() - 1));

    sim.advance(5);
    ASSERT_EQ(sim.get_result().get_num_time_points(), 6);
    ASSERT_THAT(sim.get_result().get_times(), testing::ElementsAre(0., 1., 2., 3., 4., 5.));
    for (auto&& v : sim.get_result()) {
        ASSERT_EQ(v.sum(), 4);
    }
}
