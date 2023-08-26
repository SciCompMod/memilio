/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "abm_helpers.h"

TEST(TestInfection, init)
{
    auto params             = mio::abm::GlobalInfectionParameters{};
    auto virus_variant_test = mio::abm::VirusVariant::Wildtype;
    auto age_group_test     = mio::abm::AgeGroup::Age15to34;

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
        .WillRepeatedly(testing::Return(1.0));

    auto infection = mio::abm::Infection(mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, params,
                                         mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed, {}, true);

    EXPECT_EQ(infection.get_virus_variant(), mio::abm::VirusVariant::Wildtype);
    EXPECT_EQ(infection.is_detected(), true);

    EXPECT_EQ(infection.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1) - mio::abm::seconds(1)),
              mio::abm::InfectionState::Exposed);
    EXPECT_EQ(infection.get_infection_state(mio::abm::TimePoint(0) + mio::abm::days(1)),
              mio::abm::InfectionState::InfectedNoSymptoms);
    EXPECT_NEAR(infection.get_infectivity(mio::abm::TimePoint(0) + mio::abm::days(3)), 0.2689414213699951, 1e-14);
}

TEST(TestInfection, getInfectionState)
{
    auto t = mio::abm::TimePoint(0);
    auto infection =
        mio::abm::Infection(mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
                            mio::abm::GlobalInfectionParameters{}, t, mio::abm::InfectionState::Exposed, {}, true);
    EXPECT_EQ(infection.get_infection_state(t), mio::abm::InfectionState::Exposed);
    EXPECT_EQ(infection.get_infection_state(t - mio::abm::TimeSpan(1)), mio::abm::InfectionState::Susceptible);
}

TEST(TestInfection, getPersonalProtectiveFactor)
{
    auto location = mio::abm::Location(mio::abm::LocationType::School, 0);
    auto person   = mio::abm::Person(location, mio::abm::AgeGroup::Age15to34);

    mio::abm::GlobalInfectionParameters params = mio::abm::GlobalInfectionParameters();

    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ExposureType::GenericVaccine, person.get_age(),
                                                       mio::abm::VirusVariant::Wildtype}] =
        [](ScalarType days) -> ScalarType {
        return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>({{2, 0.91}, {30, 0.81}}, days);
    };
    params.get<mio::abm::HighViralLoadProtectionFactor>() = [](ScalarType days) -> ScalarType {
        return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>({{2, 0.91}, {30, 0.81}}, days);
    };
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ExposureType::GenericVaccine, person.get_age(),
                                                      mio::abm::VirusVariant::Wildtype}] =
        [](ScalarType days) -> ScalarType {
        return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>({{2, 0.91}, {30, 0.81}}, days);
    };
    person.add_new_vaccination(mio::abm::ExposureType::GenericVaccine, mio::abm::TimePoint(0));
    auto latest_protection = person.get_latest_protection();

    auto t = mio::abm::TimePoint(2 * 24 * 60 * 60);
    ASSERT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0, 0.0001);
    t = mio::abm::TimePoint(15 * 24 * 60 * 60);
    ASSERT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0.8635, 0.0001);
    t = mio::abm::TimePoint(40 * 24 * 60 * 60);
    ASSERT_NEAR(person.get_protection_factor(t, mio::abm::VirusVariant::Wildtype, params), 0, 0.0001);

    t                               = mio::abm::TimePoint(2 * 24 * 60 * 60);
    auto severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.first, mio::abm::AgeGroup::Age15to34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.second.days());
    ASSERT_NEAR(severity_protection_factor, 0, 0.0001);
    t                          = mio::abm::TimePoint(15 * 24 * 60 * 60);
    severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.first, mio::abm::AgeGroup::Age15to34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.second.days());
    ASSERT_NEAR(severity_protection_factor, 0.8635, 0.0001);
    t                          = mio::abm::TimePoint(40 * 24 * 60 * 60);
    severity_protection_factor = params.get<mio::abm::SeverityProtectionFactor>()[{
        latest_protection.first, mio::abm::AgeGroup::Age15to34, mio::abm::VirusVariant::Wildtype}](
        t.days() - latest_protection.second.days());
    ASSERT_NEAR(severity_protection_factor, 0, 0.0001);

    t = mio::abm::TimePoint(2 * 24 * 60 * 60);
    ASSERT_NEAR(params.get<mio::abm::HighViralLoadProtectionFactor>()(t.days()), 0, 0.0001);
    t = mio::abm::TimePoint(15 * 24 * 60 * 60);
    ASSERT_NEAR(params.get<mio::abm::HighViralLoadProtectionFactor>()(t.days()), 0.8635, 0.0001);
    t = mio::abm::TimePoint(40 * 24 * 60 * 60);
    ASSERT_NEAR(params.get<mio::abm::HighViralLoadProtectionFactor>()(t.days()), 0, 0.0001);
}
