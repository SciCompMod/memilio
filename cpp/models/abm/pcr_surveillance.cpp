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
#include "abm/test_type.h"
#include "memilio/utils/random_number_generator.h"

#include <algorithm>
#include <cmath>

namespace mio
{
namespace abm
{

void PcrSurveillance::run(TimePoint t, TimeSpan /*dt*/, Model& model)
{
    const bool surveillance_on = m_budget.event_budget_fraction > 0.0;
    const bool diagnostic_on   = m_p_careseek > 0.0;
    if (!surveillance_on && !diagnostic_on) {
        return;
    }

    resolve_pending(t, model);

    if (t.hour_of_day() != m_budget.test_hour) {
        return;
    }
    if (diagnostic_on) {
        run_diagnostic(t, model);
    }
    if (surveillance_on) {
        const int period = std::max(1, m_budget.test_period_days);
        const int day    = static_cast<int>(std::floor(t.days()));
        if (day % period == 0) {
            run_surveillance(t, model);
        }
    }
}

void PcrSurveillance::resolve_pending(TimePoint t, Model& model)
{
    m_resolved_records.clear();
    auto it = m_pending.begin();
    while (it != m_pending.end()) {
        if (it->report_time <= t) {
            Person& person = model.get_person(it->person_id);
            auto rng       = PersonalRandomNumberGenerator(model.get_rng(), person);
            if (it->positive) {
                // Fixed-duration isolation, subject to the Person's Isolation compliance
                // (checked inside start_isolation).
                person.start_isolation(rng, t);
            }
            person.add_test_result(t, TestType::PCR, it->positive, it->ct_value);
            m_resolved_records.push_back(*it);
            it = m_pending.erase(it);
        }
        else {
            ++it;
        }
    }
}

void PcrSurveillance::run_surveillance(TimePoint t, Model& model)
{
    // Random sample of the non-isolated, alive, currently test-eligible population
    std::vector<PersonId> frame;
    for (auto& person : model.get_persons()) {
        if (person.is_in_quarantine(t, model.parameters) ||
            person.get_infection_state(t) == InfectionState::Dead ||
            !is_eligible_for_test(person.get_id(), t)) {
            continue;
        }
        frame.push_back(person.get_id());
    }

    const auto population_size = model.get_persons().size();
    const auto n_sample =
        static_cast<size_t>(std::llround(static_cast<double>(population_size) * m_budget.event_budget_fraction));

    std::shuffle(frame.begin(), frame.end(), mio::thread_local_rng());
    const size_t n_take = std::min(n_sample, frame.size());
    for (size_t i = 0; i < n_take; ++i) {
        test_person(frame[i], t, model, TestSource::Survey);
    }
}

void PcrSurveillance::run_diagnostic(TimePoint t, Model& model)
{
    // Each symptomatic, non-isolated Person seeks a test today with hazard
    // m_p_careseek. False Positives are not modeled.
    for (auto& person : model.get_persons()) {
        if (person.is_in_quarantine(t, model.parameters) ||
            person.get_infection_state(t) != InfectionState::InfectedSymptoms ||
            !is_eligible_for_test(person.get_id(), t)) {
            continue;
        }
        auto rng = PersonalRandomNumberGenerator(model.get_rng(), person);
        if (UniformDistribution<ScalarType>::get_instance()(rng) < m_p_careseek) {
            test_person(person.get_id(), t, model, TestSource::Diagnostic);
        }
    }
}

void PcrSurveillance::test_person(PersonId id, TimePoint t, Model& model, TestSource source)
{
    Person& person      = model.get_person(id);
    auto rng            = PersonalRandomNumberGenerator(model.get_rng(), person);
    auto [positive, ct] = person.get_tested_pcr(rng, t, m_assay);
    const TimePoint report_time = t + m_budget.reporting_delay;
    m_pending.push_back(
        PcrTestRecord{id, person.get_age(), person.get_location_type(), t, report_time, ct, positive, source});
    m_last_test[id.get()] = LastTest{t, report_time};
}

bool PcrSurveillance::is_eligible_for_test(PersonId id, TimePoint t) const
{
    auto it = m_last_test.find(id.get());
    if (it == m_last_test.end()) {
        return true;
    }
    // Not while a previous test has a pending result, and not within
    // min_retest_gap of the last sample. As a consequence, rules out same-day 
    // double-testing across streams and back-to-back testing on consecutive days.
    return t >= it->second.report_time && (t - it->second.test_time) >= m_budget.min_retest_gap;
}

} // namespace abm
} // namespace mio
