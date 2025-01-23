/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: David Kerkmann, Khoa Nguyen
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

#include "abm/location_type.h"
#include "abm/person.h"
#include "abm_helpers.h"
#include "random_number_test.h"

using TestInfection = RandomNumberTest;

/**
 * @brief Test infection initialization, viral load setup, and transitions.
 */
TEST_F(TestInfection, init)
{
    auto params             = mio::abm::Parameters(num_age_groups);
    auto virus_variant_test = mio::abm::VirusVariant::Wildtype;
    auto age_group_test     = age_group_15_to_34;
    mio::abm::Location loc(mio::abm::LocationType::Hospital, 0);

    //set up a personal RNG for infections
    //uses uniformdistribution but result doesn't matter, so init before the mock
    auto counter = mio::Counter<uint32_t>(0);
    auto prng    = mio::abm::PersonalRandomNumberGenerator(this->get_rng().get_key(), mio::abm::PersonId(0), counter);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(7))
        .WillOnce(testing::Return(0.4)) // Transition to Infected
        .WillOnce(testing::Return(0.6)) // Transition to Recovered
        .WillOnce(testing::Return(params.get<mio::abm::ViralLoadDistributions>()[{virus_variant_test, age_group_test}]
                                      .viral_load_peak.params.a())) // Viral load draws
        .WillOnce(testing::Return(params.get<mio::abm::ViralLoadDistributions>()[{virus_variant_test, age_group_test}]
                                      .viral_load_incline.params.a()))
        .WillOnce(testing::Return(params.get<mio::abm::ViralLoadDistributions>()[{virus_variant_test, age_group_test}]
                                      .viral_load_decline.params.a()))
        .WillOnce(testing::Return(params.get<mio::abm::InfectivityDistributions>()[{virus_variant_test, age_group_test}]
                                      .infectivity_alpha.params.a())) // Infectivity draws
        .WillOnce(testing::Return(params.get<mio::abm::InfectivityDistributions>()[{virus_variant_test, age_group_test}]
                                      .infectivity_beta.params.a()))
        .WillOnce(testing::Return(0.1)) // Transition to Infected
        .WillOnce(testing::Return(0.1)) // Transition to Recovered
        .WillOnce(testing::Return(params.get<mio::abm::ViralLoadDistributions>()[{virus_variant_test, age_group_test}]
                                      .viral_load_peak.params.a())) // Viral load draws
        .WillOnce(testing::Return(params.get<mio::abm::ViralLoadDistributions>()[{virus_variant_test, age_group_test}]
                                      .viral_load_incline.params.a()))
        .WillOnce(testing::Return(params.get<mio::abm::ViralLoadDistributions>()[{virus_variant_test, age_group_test}]
                                      .viral_load_decline.params.a()))
        .WillOnce(testing::Return(params.get<mio::abm::InfectivityDistributions>()[{virus_variant_test, age_group_test}]
                                      .infectivity_alpha.params.a())) // Infectivity draws
        .WillOnce(testing::Return(params.get<mio::abm::InfectivityDistributions>()[{virus_variant_test, age_group_test}]
                                      .infectivity_beta.params.a()))
        .WillRepeatedly(testing::Return(1.0));

    auto infection = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_15_to_34, params,
                                         mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed,
                                         {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, true);

    // Test virus variant and detection status
    EXPECT_EQ(infection.get_virus_variant(), mio::abm::VirusVariant::Wildtype);
    EXPECT_EQ(infection.is_detected(), true);
    // Test state transitions based on time
    EXPECT_EQ(infection.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1) - mio::abm::seconds(1)),
              mio::abm::InfectionState::Exposed);
    EXPECT_EQ(infection.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1)),
              mio::abm::InfectionState::InfectedNoSymptoms);
    // Test infectivity at a specific time point
    EXPECT_NEAR(infection.get_infectivity(mio::abm::TimePoint(0) + mio::abm::days(3)), 0.2689414213699951, 1e-14);

    // Test infection with previous exposure and recovery state transition.
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_test,
                                                      virus_variant_test}] =
        mio::TimeSeriesFunctor<ScalarType>{mio::TimeSeriesFunctorType::LinearInterpolation, {{0, 0.91}, {30, 0.81}}};
    params.get<mio::abm::HighViralLoadProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_test,
                                                           virus_variant_test}] =
        mio::TimeSeriesFunctor<ScalarType>{mio::TimeSeriesFunctorType::LinearInterpolation, {{0, 0.91}, {30, 0.81}}};
    auto infection_w_previous_exp =
        mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_test, params, mio::abm::TimePoint(0),
                            mio::abm::InfectionState::InfectedSymptoms,
                            {mio::abm::ProtectionType::GenericVaccine, mio::abm::TimePoint(0)}, true);
    // Test infection state transition
    EXPECT_EQ(
        infection_w_previous_exp.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1) - mio::abm::seconds(1)),
        mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(infection_w_previous_exp.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1)),
              mio::abm::InfectionState::Recovered);
    // Test infectivity at a specific time point
    EXPECT_NEAR(infection_w_previous_exp.get_infectivity(mio::abm::TimePoint(0) + mio::abm::days(3)),
                0.45760205922564895, 1e-14);
}

/**
 * @brief Test getInfectionState function at different time points.
 */
TEST_F(TestInfection, getInfectionState)
{
    auto counter = mio::Counter<uint32_t>(0);
    auto prng    = mio::abm::PersonalRandomNumberGenerator(this->get_rng().get_key(), mio::abm::PersonId(0), counter);
    auto params     = mio::abm::Parameters(num_age_groups);
    auto t          = mio::abm::TimePoint(0);

    // Initialize infection in Exposed state
    auto infection = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_15_to_34, params, t,
                                         mio::abm::InfectionState::Exposed,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, true);

    // Test infection state at different time points
    EXPECT_EQ(infection.get_infection_state(t), mio::abm::InfectionState::Exposed);
    EXPECT_EQ(infection.get_infection_state(t - mio::abm::TimeSpan(1)), mio::abm::InfectionState::Susceptible);
}

/**
 * @brief Test infection state forward transitions.
 */
TEST_F(TestInfection, drawInfectionCourseForward)
{
    auto counter = mio::Counter<uint32_t>(0);
    auto prng    = mio::abm::PersonalRandomNumberGenerator(this->get_rng().get_key(), mio::abm::PersonId(0), counter);
    auto params     = mio::abm::Parameters(num_age_groups);
    auto t          = mio::abm::TimePoint(0);

    // Mock recovery transition
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 1;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(6)) // First five draws for viral load
        .WillRepeatedly(testing::Return(0.8)); // Sixth draw: Recovered draw in drawInfectionCourseForward
    auto infection1 = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_15_to_34, params, t,
                                         mio::abm::InfectionState::InfectedCritical,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, true);
    // Test state transitions from Critical to Recovered
    EXPECT_EQ(infection1.get_infection_state(t), mio::abm::InfectionState::InfectedCritical);
    EXPECT_EQ(infection1.get_infection_state(t + mio::abm::days(1)), mio::abm::InfectionState::Recovered);

    // Mock death transition
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 1;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(6)) // First five draws for viral load
        .WillRepeatedly(testing::Return(0.2)); // Sixth draw: Dead draw in drawInfectionCourseForward
    auto infection2 = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_15_to_34, params, t,
                                          mio::abm::InfectionState::InfectedSevere,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, true);
    EXPECT_EQ(infection2.get_infection_state(t), mio::abm::InfectionState::InfectedSevere);
    EXPECT_EQ(infection2.get_infection_state(t + mio::abm::days(1)), mio::abm::InfectionState::Dead);
}

/**
 * @brief Test infection state backward transitions.
 */
TEST_F(TestInfection, drawInfectionCourseBackward)
{
    auto counter = mio::Counter<uint32_t>(0);
    auto prng    = mio::abm::PersonalRandomNumberGenerator(this->get_rng().get_key(), mio::abm::PersonId(0), counter);
    auto t       = mio::abm::TimePoint(1);
    auto dt      = mio::abm::days(1);
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);

    // Time to go from all infected states to recover is 1 day (dt).
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]             = 1;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]           = 1;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]   = 1;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 1;

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(14))
        .WillOnce(testing::Return(0.1)) // Transition to InfectedNoSymptoms
        .WillOnce(testing::Return(0.1))
        .WillOnce(testing::Return(0.1))
        .WillOnce(testing::Return(0.1))
        .WillOnce(testing::Return(0.1))
        .WillOnce(testing::Return(0.1))
        .WillOnce(testing::Return(0.3)) // Transition to InfectedSymptoms
        .WillOnce(testing::Return(0.3))
        .WillOnce(testing::Return(0.3))
        .WillOnce(testing::Return(0.3))
        .WillOnce(testing::Return(0.3))
        .WillOnce(testing::Return(0.3))
        .WillOnce(testing::Return(0.6)) // Transition to InfectedSevere
        .WillRepeatedly(testing::Return(0.9)); // Transition to InfectedCritical

    auto infection1 = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_60_to_79, params,
                                          mio::abm::TimePoint(t + dt), mio::abm::InfectionState::Recovered,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, false);
    auto infection2 = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_60_to_79, params,
                                          mio::abm::TimePoint(t + dt), mio::abm::InfectionState::Recovered,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, false);
    auto infection3 = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_60_to_79, params,
                                          mio::abm::TimePoint(t + dt), mio::abm::InfectionState::Recovered,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, false);
    auto infection4 = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_60_to_79, params,
                                          mio::abm::TimePoint(t + dt), mio::abm::InfectionState::Recovered,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, false);

    // Validate infection state progression backward.
    EXPECT_EQ(infection1.get_infection_state(t), mio::abm::InfectionState::InfectedNoSymptoms);
    EXPECT_EQ(infection2.get_infection_state(t), mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(infection3.get_infection_state(t), mio::abm::InfectionState::InfectedSevere);
    EXPECT_EQ(infection4.get_infection_state(t), mio::abm::InfectionState::InfectedCritical);
}

/**
 * @brief Test personal protective factors for infection and severity.
 */
TEST_F(TestInfection, getPersonalProtectiveFactor)
{
    const ScalarType eps = 1e-4;

    auto location = mio::abm::Location(mio::abm::LocationType::School, 0, num_age_groups);
    auto person   = mio::abm::Person(this->get_rng(), location.get_type(), location.get_id(), age_group_15_to_34);
    person.add_new_vaccination(mio::abm::ProtectionType::GenericVaccine, mio::abm::TimePoint(0));
    auto latest_protection = person.get_latest_protection();

    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Test default parameter functions
    auto defaut_infection_protection       = params.get<mio::abm::InfectionProtectionFactor>()[{
        mio::abm::ProtectionType::GenericVaccine, mio::AgeGroup(0), mio::abm::VirusVariant::Wildtype}](0);
    auto defaut_severity_protection        = params.get<mio::abm::SeverityProtectionFactor>()[{
        mio::abm::ProtectionType::GenericVaccine, mio::AgeGroup(0), mio::abm::VirusVariant::Wildtype}](0);
    auto defaut_high_viral_load_protection = params.get<mio::abm::HighViralLoadProtectionFactor>()[{
        mio::abm::ProtectionType::GenericVaccine, mio::AgeGroup(0), mio::abm::VirusVariant::Wildtype}](0);
    EXPECT_NEAR(defaut_infection_protection, 0, eps);
    EXPECT_NEAR(defaut_severity_protection, 0, eps);
    EXPECT_NEAR(defaut_high_viral_load_protection, 0, eps);
    // Test linear interpolation with one node
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, person.get_age(),
                                                       mio::abm::VirusVariant::Wildtype}] =
        mio::TimeSeriesFunctor<ScalarType>{mio::TimeSeriesFunctorType::LinearInterpolation, {{2, 0.91}}};
    auto t = mio::abm::TimePoint(6 * 24 * 60 * 60);

    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.91, eps);
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, person.get_age(),
                                                       mio::abm::VirusVariant::Wildtype}] =
        mio::TimeSeriesFunctor<ScalarType>{mio::TimeSeriesFunctorType::LinearInterpolation, {{2, 0.91}, {30, 0.81}}};
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, person.get_age(),
                                                      mio::abm::VirusVariant::Wildtype}] =
        mio::TimeSeriesFunctor<ScalarType>{mio::TimeSeriesFunctorType::LinearInterpolation, {{2, 0.91}, {30, 0.81}}};
    params.get<mio::abm::HighViralLoadProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, person.get_age(),
                                                           mio::abm::VirusVariant::Wildtype}] =

        mio::TimeSeriesFunctor<ScalarType>{mio::TimeSeriesFunctorType::LinearInterpolation, {{2, 0.91}, {30, 0.81}}};

    // Test Parameter InfectionProtectionFactor and get_protection_factor()
    t                                = mio::abm::TimePoint(0) + mio::abm::days(2);
    auto infection_protection_factor = params.get<mio::abm::InfectionProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.time.days());
    EXPECT_NEAR(infection_protection_factor, 0.91, eps);
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.91, eps);

    t                           = mio::abm::TimePoint(0) + mio::abm::days(15);
    infection_protection_factor = params.get<mio::abm::InfectionProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.time.days());
    EXPECT_NEAR(infection_protection_factor, 0.8635, eps);
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.8635, eps);

    t                           = mio::abm::TimePoint(0) + mio::abm::days(40);
    infection_protection_factor = params.get<mio::abm::InfectionProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.time.days());
    EXPECT_NEAR(infection_protection_factor, 0.81, eps);
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.81, eps);

    // Test Parameter SeverityProtectionFactor
    t                               = mio::abm::TimePoint(0) + mio::abm::days(2);
    auto severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.time.days());
    EXPECT_NEAR(severity_protection_factor, 0.91, eps);

    t                          = mio::abm::TimePoint(0) + mio::abm::days(15);
    severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.time.days());
    EXPECT_NEAR(severity_protection_factor, 0.8635, eps);

    t                          = mio::abm::TimePoint(0) + mio::abm::days(40);
    severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.time.days());
    EXPECT_NEAR(severity_protection_factor, 0.81, eps);

    // Test Parameter HighViralLoadProtectionFactor
    t                                 = mio::abm::TimePoint(0) + mio::abm::days(2);
    auto high_viral_protection_factor = params.get<mio::abm::HighViralLoadProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](t.days());
    EXPECT_NEAR(high_viral_protection_factor, 0.91, eps);
    t                            = mio::abm::TimePoint(0) + mio::abm::days(15);
    high_viral_protection_factor = params.get<mio::abm::HighViralLoadProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](t.days());
    EXPECT_NEAR(high_viral_protection_factor, 0.8635, eps);
    t                            = mio::abm::TimePoint(0) + mio::abm::days(40);
    high_viral_protection_factor = params.get<mio::abm::HighViralLoadProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](t.days());
    EXPECT_NEAR(high_viral_protection_factor, 0.81, eps);
}
