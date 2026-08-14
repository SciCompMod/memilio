/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Nils Wassmuth
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
#include "abm/pcr_surveillance.h"
#include "abm/model.h"
#include "abm_helpers.h"
#include "random_number_test.h"

using TestPcrSurveillance = RandomNumberTest;

namespace
{
// Assay parameters tuned so infected Persons test positive and uninfected Persons test
// true-negative deterministically without needing to mock the underlying distributions.
mio::abm::PcrAssayParameters deterministic_assay()
{
    mio::abm::PcrAssayParameters params;
    params.sample_adequacy = 1.0;
    params.cutoff          = 1e9;
    params.specificity     = 1.0;
    return params;
}

// A surveillance-only design that samples the whole (tiny test) population every day.
// Diagnostic care-seeking is disabled (p_careseek = 0).
mio::abm::TestingBudget survey_design(mio::TimeSpan reporting_delay, int test_hour)
{
    mio::abm::TestingBudget design;
    design.budget_fraction  = 1.0;
    design.test_period_days = 1;
    design.test_hour        = test_hour;
    design.reporting_delay  = reporting_delay;
    return design;
}
} // namespace

TEST_F(TestPcrSurveillance, inactiveByDefault)
{
    mio::abm::Model model(num_age_groups);
    auto home = model.add_location(mio::abm::LocationType::Home);
    add_test_person(model, home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);

    mio::abm::PcrSurveillance surveillance;
    surveillance.run(mio::abm::TimePoint(0) + mio::abm::hours(8), mio::abm::hours(1), model);

    EXPECT_TRUE(surveillance.get_resolved_records().empty());
}

TEST_F(TestPcrSurveillance, surveillanceTestsNonIsolatedPopulation)
{
    mio::abm::Model model(num_age_groups);
    auto home     = model.add_location(mio::abm::LocationType::Home);
    const auto t0 = mio::abm::TimePoint(0) + mio::abm::hours(8);

    // Two infected and one uninfected Person, none isolating.
    add_test_person(model, home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms, t0);
    add_test_person(model, home, age_group_35_to_59, mio::abm::InfectionState::InfectedSymptoms, t0);
    add_test_person(model, home, age_group_60_to_79, mio::abm::InfectionState::Susceptible, t0);

    mio::abm::PcrSurveillance surveillance;
    surveillance.configure(survey_design(mio::abm::hours(0), t0.hour_of_day()), deterministic_assay(), 0.0);

    surveillance.run(t0, mio::abm::hours(1), model); // samples all 3; nothing resolved yet
    EXPECT_TRUE(surveillance.get_resolved_records().empty());

    surveillance.run(t0 + mio::abm::hours(1), mio::abm::hours(1), model); // resolves the pending tests
    const auto& records = surveillance.get_resolved_records();
    ASSERT_EQ(records.size(), 3u);

    int positives = 0;
    for (const auto& r : records) {
        EXPECT_EQ(r.source, mio::abm::TestSource::Survey);
        positives += r.positive ? 1 : 0;
    }
    EXPECT_EQ(positives, 2) << "the two infected Persons must test positive, the uninfected one negative";
}

TEST_F(TestPcrSurveillance, positiveResultIsDeferredByReportingDelay)
{
    mio::abm::Model model(num_age_groups);
    auto home = model.add_location(mio::abm::LocationType::Home);

    const auto t0 = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto person_id =
        add_test_person(model, home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms, t0);
    auto& person = model.get_person(person_id);
    person.set_compliance(mio::abm::InterventionType::Isolation, 1.0);

    mio::abm::PcrSurveillance surveillance;
    surveillance.configure(survey_design(mio::abm::hours(24), t0.hour_of_day()), deterministic_assay(), 0.0);

    surveillance.run(t0, mio::abm::hours(1), model);
    EXPECT_FALSE(person.is_in_quarantine(t0, model.parameters)) << "isolation must not start at sample time";

    surveillance.run(t0 + mio::abm::hours(23), mio::abm::hours(1), model);
    EXPECT_FALSE(person.is_in_quarantine(t0 + mio::abm::hours(23), model.parameters))
        << "isolation must not start before the reporting delay has elapsed";

    surveillance.run(t0 + mio::abm::hours(24), mio::abm::hours(1), model);
    EXPECT_TRUE(person.is_in_quarantine(t0 + mio::abm::hours(24), model.parameters))
        << "isolation must start once the reporting delay has elapsed";
}

TEST_F(TestPcrSurveillance, isolatedPersonsAreExcludedFromSurvey)
{
    mio::abm::Model model(num_age_groups);
    auto home     = model.add_location(mio::abm::LocationType::Home);
    const auto t0 = mio::abm::TimePoint(0) + mio::abm::hours(8);

    auto isolating_id = add_test_person(model, home, age_group_15_to_34, mio::abm::InfectionState::Susceptible, t0);
    auto& isolating   = model.get_person(isolating_id);
    isolating.set_compliance(mio::abm::InterventionType::Isolation, 1.0);
    auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), isolating);
    isolating.start_isolation(rng, t0);
    ASSERT_TRUE(isolating.is_in_quarantine(t0, model.parameters));

    auto eligible_id = add_test_person(model, home, age_group_35_to_59, mio::abm::InfectionState::Susceptible, t0);

    mio::abm::PcrSurveillance surveillance;
    surveillance.configure(survey_design(mio::abm::hours(0), t0.hour_of_day()), deterministic_assay(), 0.0);

    surveillance.run(t0, mio::abm::hours(1), model); // samples only the non-isolated Person
    surveillance.run(t0 + mio::abm::hours(1), mio::abm::hours(1), model); // resolves it

    const auto& records = surveillance.get_resolved_records();
    ASSERT_EQ(records.size(), 1u) << "an isolated Person must not be part of the surveillance frame";
    EXPECT_EQ(records.front().person_id, eligible_id);
}

TEST_F(TestPcrSurveillance, surveillanceRespectsCadence)
{
    mio::abm::Model model(num_age_groups);
    auto home     = model.add_location(mio::abm::LocationType::Home);
    const auto t0 = mio::abm::TimePoint(0) + mio::abm::hours(8);

    add_test_person(model, home, age_group_15_to_34, mio::abm::InfectionState::Susceptible, t0);

    mio::abm::TestingBudget design = survey_design(mio::abm::hours(1), t0.hour_of_day());
    design.test_period_days        = 2; // sample only on days 0, 2, 4, ...

    mio::abm::PcrSurveillance surveillance;
    surveillance.configure(design, deterministic_assay(), 0.0);

    // Day 0: sampled, then resolved.
    surveillance.run(t0, mio::abm::hours(1), model);
    surveillance.run(t0 + mio::abm::hours(1), mio::abm::hours(1), model);
    EXPECT_EQ(surveillance.get_resolved_records().size(), 1u) << "day 0 is a surveillance day";

    // Day 1: off-cadence, nothing to resolve.
    const auto day1 = t0 + mio::abm::days(1);
    surveillance.run(day1, mio::abm::hours(1), model);
    surveillance.run(day1 + mio::abm::hours(1), mio::abm::hours(1), model);
    EXPECT_TRUE(surveillance.get_resolved_records().empty()) << "day 1 is off-cadence, no sample taken";

    // Day 2: sampled again.
    const auto day2 = t0 + mio::abm::days(2);
    surveillance.run(day2, mio::abm::hours(1), model);
    surveillance.run(day2 + mio::abm::hours(1), mio::abm::hours(1), model);
    EXPECT_EQ(surveillance.get_resolved_records().size(), 1u) << "day 2 is a surveillance day again";
}

TEST_F(TestPcrSurveillance, diagnosticCareSeekingTestsSymptomaticOnly)
{
    mio::abm::Model model(num_age_groups);
    auto home     = model.add_location(mio::abm::LocationType::Home);
    const auto t0 = mio::abm::TimePoint(0) + mio::abm::hours(8);

    auto symptomatic_id = add_test_person(model, home, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms, t0);
    auto& symptomatic   = model.get_person(symptomatic_id);
    symptomatic.set_compliance(mio::abm::InterventionType::Isolation, 1.0);
    add_test_person(model, home, age_group_35_to_59, mio::abm::InfectionState::InfectedNoSymptoms, t0);

    mio::abm::TestingBudget design;
    design.budget_fraction  = 0.0;
    design.test_hour        = t0.hour_of_day();
    design.reporting_delay  = mio::abm::hours(1);

    mio::abm::PcrSurveillance surveillance;
    surveillance.configure(design, deterministic_assay(), /*p_careseek=*/1.0);

    surveillance.run(t0, mio::abm::hours(1), model); // symptomatic Person seeks a test
    surveillance.run(t0 + mio::abm::hours(1), mio::abm::hours(1), model);

    const auto& records = surveillance.get_resolved_records();
    ASSERT_EQ(records.size(), 1u) << "only the symptomatic Person seeks a diagnostic test";
    EXPECT_EQ(records.front().person_id, symptomatic_id);
    EXPECT_EQ(records.front().source, mio::abm::TestSource::Diagnostic);
    EXPECT_TRUE(records.front().positive);
    EXPECT_TRUE(symptomatic.is_in_quarantine(t0 + mio::abm::hours(1), model.parameters))
        << "a resolved positive diagnostic test must start isolation";
}
