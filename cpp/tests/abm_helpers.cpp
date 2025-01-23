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
#include "abm_helpers.h"
#include "abm/person.h"
#include "abm/person_id.h"
#include "memilio/utils/random_number_generator.h"

mio::abm::Person make_test_person(mio::RandomNumberGenerator& rng, mio::abm::Location& location, mio::AgeGroup age,
                                  mio::abm::InfectionState infection_state, mio::abm::TimePoint t,
                                  mio::abm::Parameters params, mio::abm::PersonId id)
{
    assert(age.get() < params.get_num_groups());
    mio::abm::Person p(rng, location.get_type(), location.get_id(), location.get_model_id(), age, id);
    if (infection_state != mio::abm::InfectionState::Susceptible) {
        auto rng_p = mio::abm::PersonalRandomNumberGenerator(rng, p);
        p.add_new_infection(
            mio::abm::Infection(rng_p, static_cast<mio::abm::VirusVariant>(0), age, params, t, infection_state));
    }
    return p;
}

mio::abm::PersonId add_test_person(mio::abm::Model& model, mio::abm::LocationId loc_id, mio::AgeGroup age,
                                   mio::abm::InfectionState infection_state, mio::abm::TimePoint t)
{
    return model.add_person(make_test_person(model.get_rng(), model.get_location(loc_id), age, infection_state, t,
                                             model.parameters, static_cast<uint64_t>(model.get_persons().size())));
}

void interact_testing(mio::abm::PersonalRandomNumberGenerator& personal_rng, mio::abm::Person& person,
                      const mio::abm::Location& location, const std::vector<mio::abm::Person>& local_population,
                      const mio::abm::TimePoint t, const mio::abm::TimeSpan dt,
                      const mio::abm::Parameters& global_parameters)
{
    // allocate and initialize air exposures with 0
    mio::abm::AirExposureRates local_air_exposure;
    local_air_exposure.resize({mio::abm::CellIndex(location.get_cells().size()), mio::abm::VirusVariant::Count});
    std::for_each(local_air_exposure.begin(), local_air_exposure.end(), [](auto& r) {
        r = 0.0;
    });
    // allocate and initialize contact exposures with 0
    mio::abm::ContactExposureRates local_contact_exposure;
    local_contact_exposure.resize({mio::abm::CellIndex(location.get_cells().size()), mio::abm::VirusVariant::Count,
                                   mio::AgeGroup(global_parameters.get_num_groups())});
    std::for_each(local_contact_exposure.begin(), local_contact_exposure.end(), [](auto& r) {
        r = 0.0;
    });
    // caclculate current exposures
    for (const mio::abm::Person& p : local_population) {
        add_exposure_contribution(local_air_exposure, local_contact_exposure, p, location, t, dt);
    }
    // run interaction
    mio::abm::interact(personal_rng, person, location, local_air_exposure, local_contact_exposure, t, dt,
                       global_parameters);
}
