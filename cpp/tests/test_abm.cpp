#include "test_abm.h"

void add_test_infection(mio::abm::Person& p, mio::abm::InfectionState infection_state, mio::abm::TimePoint t,
                        mio::abm::GlobalInfectionParameters params)
{
    // compute start time of Infection to match InfectionState at desired TimePoint
    // always takes a path with the Infected status back to Exposed, if not starting with RecoveredCarrier
    // important: To create the correct path, mocks are used. A person will always recover.
    mio::abm::TimePoint inf_start(t);
    mio::abm::InfectionState start_state(infection_state);
    //setup rng mock so the person has a state transition to infection_state
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<ScalarType>>>>
        mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(testing::AnyNumber()).WillRepeatedly(testing::Return(1.0));

    if (start_state == mio::abm::InfectionState::Recovered_Infected) {
        inf_start =
            inf_start -
            mio::abm::TimeSpan(params.get<mio::abm::InfectedToRecovered>()[{
                static_cast<mio::abm::VirusVariant>(0), p.get_age(), mio::abm::VaccinationState::Unvaccinated}]);
        start_state = mio::abm::InfectionState::Infected;
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.6)).RetiresOnSaturation();
    }
    if (start_state == mio::abm::InfectionState::Recovered_Carrier) {
        inf_start =
            inf_start -
            mio::abm::TimeSpan(params.get<mio::abm::CarrierToRecovered>()[{
                static_cast<mio::abm::VirusVariant>(0), p.get_age(), mio::abm::VaccinationState::Unvaccinated}]);
        start_state = mio::abm::InfectionState::Carrier;
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.6)).RetiresOnSaturation();
    }
    if (start_state == mio::abm::InfectionState::Infected_Critical) {
        inf_start =
            inf_start -
            mio::abm::TimeSpan(params.get<mio::abm::SevereToCritical>()[{
                static_cast<mio::abm::VirusVariant>(0), p.get_age(), mio::abm::VaccinationState::Unvaccinated}]);
        start_state = mio::abm::InfectionState::Infected_Severe;
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.4)).RetiresOnSaturation();
    }
    if (start_state == mio::abm::InfectionState::Infected_Severe) {
        inf_start =
            inf_start -
            mio::abm::TimeSpan(params.get<mio::abm::InfectedToSevere>()[{
                static_cast<mio::abm::VirusVariant>(0), p.get_age(), mio::abm::VaccinationState::Unvaccinated}]);
        start_state = mio::abm::InfectionState::Infected;
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.4)).RetiresOnSaturation();
    }
    if (start_state == mio::abm::InfectionState::Infected) {
        inf_start =
            inf_start -
            mio::abm::TimeSpan(params.get<mio::abm::CarrierToInfected>()[{
                static_cast<mio::abm::VirusVariant>(0), p.get_age(), mio::abm::VaccinationState::Unvaccinated}]);
        start_state = mio::abm::InfectionState::Carrier;
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.4)).RetiresOnSaturation();
    }
    if (start_state == mio::abm::InfectionState::Carrier) {
        inf_start =
            inf_start -
            mio::abm::TimeSpan(params.get<mio::abm::IncubationPeriod>()[{
                static_cast<mio::abm::VirusVariant>(0), p.get_age(), mio::abm::VaccinationState::Unvaccinated}]);
    }

    p.add_new_infection(mio::abm::Infection(static_cast<mio::abm::VirusVariant>(0), p.get_age(), params, inf_start));
    assert(p.get_infection_state(t) == infection_state &&
           "Error in adding infection. Desired infection state does not match state of infection.");
}

mio::abm::Person make_test_person(mio::abm::Location& location, mio::abm::AgeGroup age_group,
                                  mio::abm::InfectionState infection_state, mio::abm::TimePoint t,
                                  mio::abm::GlobalInfectionParameters params)
{
    mio::abm::Person p = mio::abm::Person(location, age_group);
    if (infection_state != mio::abm::InfectionState::Susceptible) {
        add_test_infection(p, infection_state, t, params);
    }
    return p;
}

mio::abm::Person& add_test_person(mio::abm::World& world, mio::abm::LocationId loc_id, mio::abm::AgeGroup age,
                                  mio::abm::InfectionState infection_state, mio::abm::TimePoint t,
                                  mio::abm::GlobalInfectionParameters params)
{
    mio::abm::Person& p = world.add_person(loc_id, age);
    if (infection_state != mio::abm::InfectionState::Susceptible) {
        add_test_infection(p, infection_state, t, params);
    }
    return p;
}