/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: David Kerkmann
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
#include "memilio/math/interpolation.h"
#include "memilio/utils/random_number_generator.h"
#include "abm_helpers.h"

TEST(TestInfection, init)
{
    auto params             = mio::abm::Parameters(num_age_groups);
    auto virus_variant_test = mio::abm::VirusVariant::Wildtype;
    auto age_group_test     = age_group_15_to_34;

    //set up a personal RNG for infections
    //uses uniformdistribution but result doesn't matter, so init before the mock
    mio::abm::Location loc(mio::abm::LocationType::Hospital, 0);
    auto counter = mio::Counter<uint32_t>(0);
    auto rng     = mio::abm::PersonalRandomNumberGenerator(mio::Key<uint64_t>{0}, mio::abm::PersonId(0), counter);

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

    auto infection = mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, age_group_15_to_34, params,
                                         mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed, {}, true);

    EXPECT_EQ(infection.get_virus_variant(), mio::abm::VirusVariant::Wildtype);
    EXPECT_EQ(infection.is_detected(), true);

    EXPECT_EQ(infection.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1) - mio::abm::seconds(1)),
              mio::abm::InfectionState::Exposed);
    EXPECT_EQ(infection.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1)),
              mio::abm::InfectionState::InfectedNoSymptoms);
    EXPECT_NEAR(infection.get_infectivity(mio::abm::TimePoint(0) + mio::abm::days(3)), 0.2689414213699951, 1e-14);

    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ExposureType::GenericVaccine, age_group_test,
                                                      virus_variant_test}] = mio::TimeDependentParameterFunctor{
        mio::TimeDependentParameterFunctor::Type::LinearInterpolation, {{0, 0.91}, {30, 0.81}}};
    params.get<mio::abm::HighViralLoadProtectionFactor>() = mio::TimeDependentParameterFunctor{
        mio::TimeDependentParameterFunctor::Type::LinearInterpolation, {{0, 0.91}, {30, 0.81}}};
    auto infection_w_previous_exp =
        mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, age_group_test, params, mio::abm::TimePoint(0),
                            mio::abm::InfectionState::InfectedSymptoms,
                            {mio::abm::ExposureType::GenericVaccine, mio::abm::TimePoint(0)}, true);
    EXPECT_EQ(
        infection_w_previous_exp.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1) - mio::abm::seconds(1)),
        mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(infection_w_previous_exp.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1)),
              mio::abm::InfectionState::Recovered);
    EXPECT_NEAR(infection_w_previous_exp.get_infectivity(mio::abm::TimePoint(0) + mio::abm::days(3)),
                0.45760205922564895, 1e-14);
}

TEST(TestInfection, getInfectionState)
{
    auto counter   = mio::Counter<uint32_t>(0);
    auto rng       = mio::abm::PersonalRandomNumberGenerator(mio::Key<uint64_t>{0}, mio::abm::PersonId(0), counter);
    auto params    = mio::abm::Parameters(num_age_groups);
    auto t         = mio::abm::TimePoint(0);
    auto infection = mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, age_group_15_to_34, params, t,
                                         mio::abm::InfectionState::Exposed, {}, true);
    EXPECT_EQ(infection.get_infection_state(t), mio::abm::InfectionState::Exposed);
    EXPECT_EQ(infection.get_infection_state(t - mio::abm::TimeSpan(1)), mio::abm::InfectionState::Susceptible);
}

TEST(TestInfection, drawInfectionCourseBackward)
{
    auto counter = mio::Counter<uint32_t>(0);
    auto rng     = mio::abm::PersonalRandomNumberGenerator(mio::Key<uint64_t>{0}, mio::abm::PersonId(0), counter);

    auto t                      = mio::abm::TimePoint(1);
    auto dt                     = mio::abm::days(1);
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

    auto infection1 = mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, age_group_60_to_79, params,
                                          mio::abm::TimePoint(t + dt), mio::abm::InfectionState::Recovered,
                                          {mio::abm::ExposureType::NoProtection, mio::abm::TimePoint(0)}, false);
    auto infection2 = mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, age_group_60_to_79, params,
                                          mio::abm::TimePoint(t + dt), mio::abm::InfectionState::Recovered,
                                          {mio::abm::ExposureType::NoProtection, mio::abm::TimePoint(0)}, false);
    auto infection3 = mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, age_group_60_to_79, params,
                                          mio::abm::TimePoint(t + dt), mio::abm::InfectionState::Recovered,
                                          {mio::abm::ExposureType::NoProtection, mio::abm::TimePoint(0)}, false);
    auto infection4 = mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, age_group_60_to_79, params,
                                          mio::abm::TimePoint(t + dt), mio::abm::InfectionState::Recovered,
                                          {mio::abm::ExposureType::NoProtection, mio::abm::TimePoint(0)}, false);

    EXPECT_EQ(infection1.get_infection_state(t), mio::abm::InfectionState::InfectedNoSymptoms);
    EXPECT_EQ(infection2.get_infection_state(t), mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(infection3.get_infection_state(t), mio::abm::InfectionState::InfectedSevere);
    EXPECT_EQ(infection4.get_infection_state(t), mio::abm::InfectionState::InfectedCritical);
}

TEST(TestInfection, getPersonalProtectiveFactor)
{
    const ScalarType eps = 1e-4;
    auto rng             = mio::RandomNumberGenerator();

    auto location = mio::abm::Location(mio::abm::LocationType::School, 0, num_age_groups);
    auto person   = mio::abm::Person(rng, location.get_type(), location.get_id(), age_group_15_to_34);
    person.add_new_vaccination(mio::abm::ExposureType::GenericVaccine, mio::abm::TimePoint(0));
    auto latest_protection = person.get_latest_protection();

    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Test default parameter functions
    auto defaut_infection_protection = params.get<mio::abm::InfectionProtectionFactor>()[{
        mio::abm::ExposureType::GenericVaccine, mio::AgeGroup(0), mio::abm::VirusVariant::Wildtype}](0);
    auto defaut_severity_protection  = params.get<mio::abm::SeverityProtectionFactor>()[{
        mio::abm::ExposureType::GenericVaccine, mio::AgeGroup(0), mio::abm::VirusVariant::Wildtype}](0);
    EXPECT_NEAR(defaut_infection_protection, 0, eps);
    EXPECT_NEAR(defaut_severity_protection, 0, eps);

    // Test linear interpolation with one node
    // mio::set_log_level(mio::LogLevel::critical); //this throws an error either way
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ExposureType::GenericVaccine, person.get_age(),
                                                       mio::abm::VirusVariant::Wildtype}] =
        mio::TimeDependentParameterFunctor{mio::TimeDependentParameterFunctor::Type::LinearInterpolation, {{2, 0.91}}};
    auto t = mio::abm::TimePoint(6 * 24 * 60 * 60);
    // TODO: Discuss: Assumption of interpolation in TDPF is that the function is constant with value at front/back entry outside of [front, back] time range. This works with one node as well and prints no errors
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.91, eps);
    // mio::set_log_level(mio::LogLevel::warn); //this throws an error either way
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ExposureType::GenericVaccine, person.get_age(),
                                                       mio::abm::VirusVariant::Wildtype}] =
        mio::TimeDependentParameterFunctor{mio::TimeDependentParameterFunctor::Type::LinearInterpolation,
                                           {{2, 0.91}, {30, 0.81}}};
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ExposureType::GenericVaccine, person.get_age(),
                                                      mio::abm::VirusVariant::Wildtype}] =
        mio::TimeDependentParameterFunctor{mio::TimeDependentParameterFunctor::Type::LinearInterpolation,
                                           {{2, 0.91}, {30, 0.81}}};
    params.get<mio::abm::HighViralLoadProtectionFactor>() = mio::TimeDependentParameterFunctor{
        mio::TimeDependentParameterFunctor::Type::LinearInterpolation, {{2, 0.91}, {30, 0.81}}};

    // Test Parameter InfectionProtectionFactor and get_protection_factor()
    t                                = mio::abm::TimePoint(0) + mio::abm::days(2);
    auto infection_protection_factor = params.get<mio::abm::InfectionProtectionFactor>()[{
        latest_protection.first, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.second.days());
    EXPECT_NEAR(infection_protection_factor, 0.91, eps);
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.91, eps);

    t                           = mio::abm::TimePoint(0) + mio::abm::days(15);
    infection_protection_factor = params.get<mio::abm::InfectionProtectionFactor>()[{
        latest_protection.first, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.second.days());
    EXPECT_NEAR(infection_protection_factor, 0.8635, eps);
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.8635, eps);

    t                           = mio::abm::TimePoint(0) + mio::abm::days(40);
    infection_protection_factor = params.get<mio::abm::InfectionProtectionFactor>()[{
        latest_protection.first, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.second.days());
    EXPECT_NEAR(infection_protection_factor, 0.81,
                eps); // TODO: why was this 0? should there be an instant falloff after last data point?
    EXPECT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.81,
                eps); // TODO: why was this 0? should there be an instant falloff after last data point?

    // Test Parameter SeverityProtectionFactor
    t                               = mio::abm::TimePoint(0) + mio::abm::days(2);
    auto severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.first, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.second.days());
    EXPECT_NEAR(severity_protection_factor, 0.91, eps);

    t                          = mio::abm::TimePoint(0) + mio::abm::days(15);
    severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.first, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.second.days());
    EXPECT_NEAR(severity_protection_factor, 0.8635, eps);

    t                          = mio::abm::TimePoint(0) + mio::abm::days(40);
    severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.first, age_group_15_to_34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.second.days());
    EXPECT_NEAR(severity_protection_factor, 0.81,
                eps); // TODO: why was this 0? should there be an instant falloff after last data point?

    // Test Parameter HighViralLoadProtectionFactor
    t = mio::abm::TimePoint(0) + mio::abm::days(2);
    EXPECT_NEAR(params.get<mio::abm::HighViralLoadProtectionFactor>()(t.days()), 0.91, eps);
    t = mio::abm::TimePoint(0) + mio::abm::days(15);
    EXPECT_NEAR(params.get<mio::abm::HighViralLoadProtectionFactor>()(t.days()), 0.8635, eps);
    t = mio::abm::TimePoint(0) + mio::abm::days(40);
    EXPECT_NEAR(params.get<mio::abm::HighViralLoadProtectionFactor>()(t.days()), 0.81,
                eps); // TODO: why was this 0? should there be an instant falloff after last data point?
}
