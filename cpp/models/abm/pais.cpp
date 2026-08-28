/*
* Copyright (C) 2020-2026 MEmilio
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

#include "abm/pais.h"
#include "abm/person.h"
#include "abm/random_events.h"

namespace mio
{
namespace abm
{

void PAIS::update_severity(const Parameters& params, PersonalRandomNumberGenerator& rng, TimePoint t, TimeSpan dt)
{
    if (severity.empty() || t <= severity.back().first) {
        return; // only update if PAIS has already started and the last update was before t
    }
    const auto current = severity.back().second;
    // Collect the transition rates out of the current state. The diagonal is skipped, since staying in
    // the current state is not an event of its own: it is the complement of all transitions happening.
    std::pair<PAISState, ScalarType> transition_rates[static_cast<uint32_t>(PAISState::Count) - 1];
    size_t num_transitions = 0;
    for (auto&& v : enum_members<PAISState>()) {
        if (v == current) {
            continue;
        }
        transition_rates[num_transitions++] = {
            v, params.get<PAISTransitionMatrix>()(static_cast<Eigen::Index>(current),
                                                  static_cast<Eigen::Index>(v))};
    }
    auto severity_new = random_transition(rng, current, dt, transition_rates);
    if (severity_new != severity.back().second) {
        this->severity.push_back({t, severity_new}); // only update if there is a change in severity
    }
}

void PAIS::init_or_refresh(const Parameters& params, Person& p, PersonalRandomNumberGenerator& rng,
                           const Infection& inf, TimePoint t)
{
    // get highest InfectionState of the new infection
    auto highest_state = inf.get_highest_infection_state();
    if (highest_state.second != InfectionState::Dead) {
        // if the Person already had an active PAIS and gets a reinfection, refresh the PAIS status
        if (get_severity(t) != PAISState::Count) {
            add_new_severity(t, highest_state);
        }
        else {
            // base probability of developing PAIS based on age, sex, virus variant and number of vaccinations
            ScalarType pais_prob = params.get<PAISProbability>()[{inf.get_virus_variant(), p.get_age(), p.get_sex(),
                                                                  get_vaccination_class(p.get_vaccinations().size())}];
            // increase probability of developing PAIS if the Person had a severe acute infection or worse
            if (highest_state.second == InfectionState::InfectedSevere ||
                highest_state.second == InfectionState::InfectedCritical) {
                pais_prob *= params.get<PAISProbabilitySeverityFactor>()[{
                    inf.get_virus_variant(), get_vaccination_class(p.get_vaccinations().size())}];
            }
            // reduce probability of developing PAIS if the Person has not had PAIS after an earlier infection
            if (!p.get_infections().empty() && get_severity(t) == PAISState::Count) {
                pais_prob *= params.get<PAISProtectionAtSecondInfection>()[{
                    inf.get_virus_variant(), get_vaccination_class(p.get_vaccinations().size())}];
            }

            // Draw only if the outcome is not already determined. PAISProbability is zero by default, so
            // this keeps the random number stream of a model that does not use PAIS unchanged.
            if (pais_prob > 0.0) {
                auto& uniform_dist = UniformDistribution<ScalarType>::get_instance();
                if (uniform_dist(rng) < pais_prob) {
                    TimePoint time_recovered = inf.get_infection_state_start_date(InfectionState::Recovered);
                    add_new_severity(time_recovered, highest_state);
                }
            }
        }
    }
}

void PAIS::add_new_severity(TimePoint t, std::pair<TimePoint, InfectionState> highest_state)
{
    PAISState severity_new;
    if (highest_state.second == InfectionState::InfectedSevere ||
        highest_state.second == InfectionState::InfectedCritical) {
        severity_new = PAISState::Severe;
    }
    else {
        severity_new = PAISState::Medium;
    }
    if (severity.empty() || (t > severity.back().first && severity_new != severity.back().second)) {
        this->severity.push_back({t, severity_new});
    }
}

PAISState PAIS::get_severity(TimePoint t) const
{
    if (severity.empty()) {
        return PAISState::Count;
    }
    if (t < severity[0].first) {
        return PAISState::Count;
    }

    auto it = std::upper_bound(severity.begin(), severity.end(), t,
                               [](const TimePoint& s, const std::pair<TimePoint, PAISState>& state) {
                                   return state.first > s;
                               });
    return std::prev(it)->second;
}

} // namespace abm
} // namespace mio
