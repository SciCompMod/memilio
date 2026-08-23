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
#ifndef MIO_ABM_PCR_SURVEILLANCE_H
#define MIO_ABM_PCR_SURVEILLANCE_H

#include "abm/location_type.h"
#include "abm/person_id.h"
#include "abm/time.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace mio
{
namespace abm
{

class Model;

/**
 * @brief Parameters of the per-agent PCR test-outcome model.
 * Infected: Ci = a - b * Vi(t) + eps, eps ~ N(0, sigma^2), where Vi is the log10 viral
 * load. The test is positive iff the sample is adequate (Bernoulli sample_adequacy) and Ci < cutoff.
 * Not infected: a contamination sample occurs with probability (1 - specificity),
 * Ci ~ N(contam_mean, contam_sd), and is reported positive under the same rule Ci < cutoff,
 * so an above-cutoff contamination sample reads negative, just like any above-cutoff sample.
 * contam_mean sits 2 SD below cutoff, so nearly all contamination samples clear it and
 * (1 - specificity) is a a close approximation of the reported false-positive rate
 */
struct PcrAssayParameters {
    ScalarType a = 40.0; ///< CT value at zero viral load.
    ScalarType b = 3.32192809489; ///< CT decrease per unit of log10 viral load (1/log10(2)).
    ScalarType sigma = 1.0; ///< Standard deviation of the CT measurement noise.
    ScalarType cutoff = 35.0; ///< Positivity cutoff c.
    ScalarType specificity = 0.99; ///< Probability of a true negative for an uninfected Person.
    ScalarType contam_mean = 33.0; ///< Mean CT of a false-positive sample (F_contam); 2 SD below cutoff.
    ScalarType contam_sd = 1.0; ///< Standard deviation of F_contam.
    ScalarType sample_adequacy = 0.95; ///< Probability s(V) that a sample from an infected Person is adequate.
};

/**
 * @brief Which testing process produced a result.
 */
enum class TestSource : std::uint32_t {
    Survey     = 0, ///< Active random surveillance
    Diagnostic = 1, ///< Demand-driven symptomatic care-seeking, suppresses the epidemic.
};

/**
 * @brief The active-surveillance testing "design" d: how much to sample and how often.
 */
struct TestingBudget {
    /// Fraction of the population sampled at each surveillance event.
    /// 0 disables surveillance.
    ScalarType event_budget_fraction = 0.0;
    int test_period_days = 1; ///< Surveillance samples every test_period_days days (cadence). 1 = daily.
    int test_hour = 8; ///< Hour of day at which sampling occurs.
    TimeSpan reporting_delay = days(1); ///< Delay between taking a sample and acting on/reporting its result.
    /// Minimum time between two PCR tests (Survey or Diagnostic) of the same Person. Since this is
    /// always >= reporting_delay by construction of the default, it also rules out testing someone
    /// while a previous result is still pending; the two conditions are nonetheless checked
    /// independently in PcrSurveillance so the guarantee holds regardless of how these are configured.
    TimeSpan min_retest_gap = days(2);
};

/**
 * @brief One resolved PCR test.
 */
struct PcrTestRecord {
    PersonId person_id;
    AgeGroup age;
    LocationType location_type; ///< Location of the Person at the time the sample was taken.
    TimePoint test_time; ///< TimePoint the sample was taken.
    TimePoint report_time; ///< TimePoint the result is acted upon/reported.
    ScalarType ct_value; ///< Reported CT value.
    bool positive;
    TestSource source; ///< Which process produced this test.
};

/**
 * @brief Two independent PCR testing processes with delayed reporting.
 *
 * - Active surveillance: every TestingBudget::test_period_days days at TestingBudget::test_hour,
 *   a random `event_budget_fraction * population` sample of the non-isolated population is
 *   PCR-tested
 * - Diagnostic care-seeking: each day at test_hour, every symptomatic, non-isolated Person seeks
 *   a test with per-day hazard `p_careseek`. This suppresses the epidemic through isolation.
 *
 * A positive result starts a fixed-duration isolation (subject to the Person's Isolation
 * compliance) once TestingBudget::reporting_delay has elapsed.
 */
class PcrSurveillance
{
public:
    PcrSurveillance() = default;

    /// @brief (Re-)configure the surveillance design, assay parameters, and diagnostic care-seeking hazard.
    void configure(const TestingBudget& budget, const PcrAssayParameters& assay, ScalarType p_careseek)
    {
        m_budget     = budget;
        m_assay      = assay;
        m_p_careseek = p_careseek;
    }

    const TestingBudget& get_budget() const
    {
        return m_budget;
    }

    const PcrAssayParameters& get_assay_parameters() const
    {
        return m_assay;
    }

    /**
     * @brief Advance testing processes by one Model step.
     */
    void run(TimePoint t, TimeSpan dt, Model& model);

    /**
     * @brief Get the PcrTestRecord%s resolved by the most recent call to run().
     */
    const std::vector<PcrTestRecord>& get_resolved_records() const
    {
        return m_resolved_records;
    }

private:
    /// @brief Time window of the most recent PCR test (either source) taken from one Person.
    struct LastTest {
        TimePoint test_time; ///< TimePoint the sample was taken.
        TimePoint report_time; ///< TimePoint the result is/was acted upon.
    };

    void resolve_pending(TimePoint t, Model& model);
    void run_surveillance(TimePoint t, Model& model); ///< periodic random sample of the non-isolated population
    void run_diagnostic(TimePoint t, Model& model); ///< daily care-seeking hazard over symptomatic non-isolated Persons
    void test_person(PersonId id, TimePoint t, Model& model, TestSource source);
    /// @brief Whether a Person may be tested (by either source) at time t: not while a previous
    /// result of theirs is still pending, and not within m_budget.min_retest_gap of their last
    /// sample -- so nobody is ever double-tested the same day or tested on back-to-back days.
    bool is_eligible_for_test(PersonId id, TimePoint t) const;

    TestingBudget m_budget{};
    PcrAssayParameters m_assay{};
    ScalarType m_p_careseek = 0.0; ///< Per-day probability a symptomatic non-isolated Person seeks a diagnostic test.
    std::vector<PcrTestRecord> m_pending; ///< Tests taken but not yet reported.
    std::vector<PcrTestRecord> m_resolved_records; ///< Results resolved by the most recent run().
    std::unordered_map<std::uint64_t, LastTest> m_last_test; ///< Most recent test window per Person, keyed by PersonId.
};

} // namespace abm
} // namespace mio

#endif // MIO_ABM_PCR_SURVEILLANCE_H
