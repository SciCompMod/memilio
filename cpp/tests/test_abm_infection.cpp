/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: David Kerkmann, Khoa Nguyen
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
#include "abm/parameters.h"
#include "abm/person.h"
#include "abm/person_id.h"
#include "abm_helpers.h"
#include "memilio/utils/abstract_parameter_distribution.h"
#include "memilio/utils/parameter_distributions.h"
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

    params.get<mio::abm::ViralShedFactor>()[{virus_variant_test, age_group_test}] =
        mio::ParameterDistributionUniform(0.1, 0.2);

    //set up a personal RNG for infections
    //uses uniformdistribution but result doesn't matter, so init before the mock
    auto counter = mio::Counter<uint32_t>(0);
    auto prng    = mio::abm::PersonalRandomNumberGenerator(this->get_rng().get_key(), 0, counter);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::LogNormalDistribution<double>>>>
        mock_logNormal_dist;

    //Distribution for state transitions
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(5))
        // 1st infection
        .WillOnce(testing::Return(0.4)) // Transition to Infected
        .WillOnce(testing::Return(0.6)) // Transition to Recovered
        .WillOnce(testing::Return(params.get<mio::abm::ViralShedFactor>()[{virus_variant_test, age_group_test}]
                                      .params()[0])) // Viral Shed Factor
        // 2nd infection
        .WillOnce(testing::Return(1.0)) // Transition to Recovered
        .WillOnce(testing::Return(params.get<mio::abm::ViralShedFactor>()[{virus_variant_test, age_group_test}]
                                      .params()[0])) // Viral Shed Factor
        .WillRepeatedly(testing::Return(1.0));

    //Distribution for stay times
    EXPECT_CALL(mock_logNormal_dist.get_mock(), invoke)
        // 1st infection
        .WillOnce(testing::Return(1.)) // IncubationTime
        .WillOnce(testing::Return(1.)) // TimeInfectedNoSymptomsToSymptoms
        .WillOnce(testing::Return(1.)) // TimeInfectedSymptomsToRecovered
        // 2nd infection
        .WillOnce(testing::Return(1.0)) // TimeInfectedNoSymptomsToSymptoms
        .WillOnce(testing::Return(1.0)) // IncubationTime
        .WillOnce(testing::Return(1.0)) // TimeInfectedSymptomsToRecovered
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
    // Test viral shed at a specific time point
    EXPECT_NEAR(infection.get_viral_shed(mio::abm::TimePoint(0) + mio::abm::days(3)), 0.02689414213699951, 1e-14);

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
    // Test viral shed at a specific time point
    EXPECT_NEAR(infection_w_previous_exp.get_viral_shed(mio::abm::TimePoint(0) + mio::abm::days(3)),
                9.1105119440064545e-05, 1e-14);
}

/**
 * @brief Test getInfectionState function at different time points.
 */
TEST_F(TestInfection, getInfectionState)
{
    auto counter = mio::Counter<uint32_t>(0);
    auto prng    = mio::abm::PersonalRandomNumberGenerator(this->get_rng().get_key(), 0, counter);
    auto params  = mio::abm::Parameters(num_age_groups);
    auto t       = mio::abm::TimePoint(0);

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
    auto prng    = mio::abm::PersonalRandomNumberGenerator(this->get_rng().get_key(), 0, counter);
    auto params  = mio::abm::Parameters(num_age_groups);
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.;
    auto t = mio::abm::TimePoint(0);

    // Mock state transition times
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::LogNormalDistribution<double>>>>
        mock_logNormal_dist;
    EXPECT_CALL(mock_logNormal_dist.get_mock(), invoke)
        .Times(testing::Exactly(13)) // 5 draws for infection 1 and 4 for infection 2 and 4 draws for infection3
        .WillRepeatedly(testing::Return(1)); // All times will be 1 day
    // Mock viral load draws and infection paths
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(
            7)) // 6 viral load draws per infection and 1 infection path draw per infection (for infection1 and 2) and 2 infection paths draws for infection3
        .WillRepeatedly(testing::Return(
            0.55)); // is necessary for infection 1 to recover, infection 2 to die and infection 3 to turn severe

    auto infection1 = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_15_to_34, params, t,
                                          mio::abm::InfectionState::InfectedCritical,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, true);
    // Test state transitions from Critical to Recovered
    EXPECT_EQ(infection1.get_infection_state(t), mio::abm::InfectionState::InfectedCritical);
    EXPECT_EQ(infection1.get_infection_state(t + mio::abm::days(1)), mio::abm::InfectionState::Recovered);

    auto infection2 = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_15_to_34, params, t,
                                          mio::abm::InfectionState::InfectedSevere,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, true);
    EXPECT_EQ(infection2.get_infection_state(t), mio::abm::InfectionState::InfectedSevere);
    EXPECT_EQ(infection2.get_infection_state(t + mio::abm::days(1)), mio::abm::InfectionState::Dead);

    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.6;
    auto infection3 = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, age_group_15_to_34, params, t,
                                          mio::abm::InfectionState::InfectedSymptoms,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, true);
    EXPECT_EQ(infection3.get_infection_state(t), mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(infection3.get_infection_state(t + mio::abm::days(1)), mio::abm::InfectionState::InfectedSevere);
}

/**
 * @brief Test infection state backward transitions.
 */
TEST_F(TestInfection, drawInfectionCourseBackward)
{
    auto counter                = mio::Counter<uint32_t>(0);
    auto prng                   = mio::abm::PersonalRandomNumberGenerator(this->get_rng().get_key(), 0, counter);
    auto t                      = mio::abm::TimePoint(1);
    auto dt                     = mio::abm::days(1);
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);

    auto virus_variant_test = mio::abm::VirusVariant::Wildtype;
    auto age_group_test     = age_group_60_to_79;

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::LogNormalDistribution<double>>>>
        mock_logNormal_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(7))
        // 1st infection
        .WillOnce(testing::Return(0.6)) // Transition to InfectedNoSymptoms
        .WillOnce(testing::Return(params.get<mio::abm::ViralShedFactor>()[{virus_variant_test, age_group_test}]
                                      .params()[0])) // Viral Shed Factor
        // 2nd infection
        .WillOnce(testing::Return(0.4)) // Transition to InfectedSymptoms
        .WillOnce(testing::Return(params.get<mio::abm::ViralShedFactor>()[{virus_variant_test, age_group_test}]
                                      .params()[0])) // Viral Shed Factor
        // 3rd infection
        .WillOnce(testing::Return(0.2)) // Transition to InfectedSevere
        .WillOnce(testing::Return(params.get<mio::abm::ViralShedFactor>()[{virus_variant_test, age_group_test}]
                                      .params()[0])) // Viral Shed Factor
        // 4th infection
        .WillOnce(testing::Return(0.05)) // Transition to InfectedCritical
        .WillRepeatedly(testing::Return(1.0));

    EXPECT_CALL(mock_logNormal_dist.get_mock(), invoke)
        .Times(testing::AtLeast(10))
        // 1st infection
        .WillOnce(testing::Return(1.0)) // TimeInfectedNoSymptomsToRecovered
        .WillOnce(testing::Return(1.0)) // TimeExposedToNoSymptoms
        // 2nd infection
        .WillOnce(testing::Return(1.0)) // TimeInfectedSymptomsToRecovered
        .WillOnce(testing::Return(1.0)) // TimeInfectedNoSymptomsToSymptoms
        .WillOnce(testing::Return(1.0)) // TimeExposedToNoSymptoms
        // 3rd infection
        .WillOnce(testing::Return(1.0)) // TimeInfectedSevereToRecovered
        .WillOnce(testing::Return(1.0)) // TimeInfectedSymptomsToSevere
        .WillOnce(testing::Return(1.0)) // TimeInfectedNoSymptomsToSymptoms
        .WillOnce(testing::Return(1.0)) // TimeExposedToNoSymptoms
        // 4th infection
        .WillOnce(testing::Return(1.0)) //TimeInfectedCriticalToRecovered
        .WillRepeatedly(testing::Return(1.0));

    auto infection1 = mio::abm::Infection(prng, virus_variant_test, age_group_test, params, mio::abm::TimePoint(t + dt),
                                          mio::abm::InfectionState::Recovered,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, false);
    auto infection2 = mio::abm::Infection(prng, virus_variant_test, age_group_test, params, mio::abm::TimePoint(t + dt),
                                          mio::abm::InfectionState::Recovered,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, false);
    auto infection3 = mio::abm::Infection(prng, virus_variant_test, age_group_test, params, mio::abm::TimePoint(t + dt),
                                          mio::abm::InfectionState::Recovered,
                                          {mio::abm::ProtectionType::NoProtection, mio::abm::TimePoint(0)}, false);
    auto infection4 = mio::abm::Infection(prng, virus_variant_test, age_group_test, params, mio::abm::TimePoint(t + dt),
                                          mio::abm::InfectionState::Recovered,
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
    auto person   = mio::abm::Person(this->get_rng(), location.get_type(), location.get_id(), location.get_model_id(),
                                     age_group_15_to_34);
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
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](t.days() -
                                                                                       latest_protection.time.days());
    EXPECT_NEAR(infection_protection_factor, 0.91, eps);
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.91, eps);

    t                           = mio::abm::TimePoint(0) + mio::abm::days(15);
    infection_protection_factor = params.get<mio::abm::InfectionProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](t.days() -
                                                                                       latest_protection.time.days());
    EXPECT_NEAR(infection_protection_factor, 0.8635, eps);
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.8635, eps);

    t                           = mio::abm::TimePoint(0) + mio::abm::days(40);
    infection_protection_factor = params.get<mio::abm::InfectionProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](t.days() -
                                                                                       latest_protection.time.days());
    EXPECT_NEAR(infection_protection_factor, 0.81, eps);
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.81, eps);

    // Test Parameter SeverityProtectionFactor
    t                               = mio::abm::TimePoint(0) + mio::abm::days(2);
    auto severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](t.days() -
                                                                                       latest_protection.time.days());
    EXPECT_NEAR(severity_protection_factor, 0.91, eps);

    t                          = mio::abm::TimePoint(0) + mio::abm::days(15);
    severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](t.days() -
                                                                                       latest_protection.time.days());
    EXPECT_NEAR(severity_protection_factor, 0.8635, eps);

    t                          = mio::abm::TimePoint(0) + mio::abm::days(40);
    severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.type, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](t.days() -
                                                                                       latest_protection.time.days());
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
