#include "functions.h"

#include "abm/random_events.h"
#include "abm/infection.h"

namespace mio
{

namespace abm
{

void interact(Person& person, Location& location, TimePoint t, TimeSpan dt, const Parameters& global_parameters,
              Person::RandomNumberGenerator& personal_rng)
{
    if (person.get_infection_state(t) == InfectionState::Susceptible) {
        auto& local_parameters = location.get_infection_parameters();
        // TODO: we need to define what a cell is used for, as the loop may lead to incorrect results for multiple cells
        auto age_receiver          = person.get_age();
        ScalarType mask_protection = person.get_mask_protective_factor(global_parameters);
        assert(person.get_cells().size() && "Person is in multiple cells. Interact logic is incorrect at the moment.");
        for (auto cell_index :
             person.get_cells()) { // TODO: the logic here is incorrect in case a person is in multiple cells
            std::pair<VirusVariant, ScalarType> local_indiv_trans_prob[static_cast<uint32_t>(VirusVariant::Count)];
            for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
                VirusVariant virus = static_cast<VirusVariant>(v);
                ScalarType local_indiv_trans_prob_v =
                    (std::min(local_parameters.get<MaximumContacts>(),
                              location.transmission_contacts_per_day(cell_index, virus, age_receiver,
                                                                     global_parameters.get_num_groups())) +
                     location.transmission_air_per_day(cell_index, virus, global_parameters)) *
                    (1 - mask_protection) * dt.days() * (1 - person.get_protection_factor(t, virus, global_parameters));

                local_indiv_trans_prob[v] = std::make_pair(virus, local_indiv_trans_prob_v);
            }
            VirusVariant virus =
                random_transition(personal_rng, VirusVariant::Count, dt,
                                  local_indiv_trans_prob); // use VirusVariant::Count for no virus submission
            if (virus != VirusVariant::Count) {
                person.add_new_infection(Infection(personal_rng, virus, age_receiver, global_parameters, t + dt / 2,
                                                   mio::abm::InfectionState::Exposed, person.get_latest_protection(),
                                                   false)); // Starting time in first approximation
            }
        }
    }
    person.add_time_at_location(dt);
}

bool migrate(Person& person, Location& destination, const std::vector<uint32_t>& cells, TransportMode mode)
{
    if (person.get_location() != destination.get_id()) {
        person.set_location(destination);
        person.get_cells() = cells;
        person.set_last_transport_mode(mode);

        return true;
    }
    else {
        return false;
    }
}

} // namespace abm

} // namespace mio