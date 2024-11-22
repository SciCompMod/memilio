#include "abm/mobility_rules.h"
#include "abm/person_id.h"
#include "memilio/utils/mioomp.h"
#include "models/abm/model.h"
#include <iostream>
#include <string>

namespace mio
{
namespace abm
{

struct ModelDefinition {
    ModelDefinition()
        : m_local_population_cache()
        , m_air_exposure_rates_cache()
        , m_contact_exposure_rates_cache()
        , m_is_local_population_cache_valid(false)
        , m_are_exposure_caches_valid(false)
        , m_exposure_caches_need_rebuild(true)
    {
    }

    void evolve(Model& m, TimePoint t, TimeSpan dt)
    {
        begin_step(m, t, dt);
        log_info("ABM Model interaction.");
        interaction(m, t, dt);
        log_info("ABM Model mobility.");
        perform_mobility(m, t, dt);
    }

    void interaction(Model& m, TimePoint t, TimeSpan dt)
    {
        const uint32_t num_persons = static_cast<uint32_t>(m.get_persons().size());
        PRAGMA_OMP(parallel for)
        for (uint32_t i = 0; i < num_persons; ++i) {
            const PersonId person_id(i);
            if (!m_are_exposure_caches_valid) {
                // checking caches is only needed for external calls
                // during simulation (i.e. in evolve()), the caches are computed in begin_step
                compute_exposure_caches(m, t, dt);
                m_are_exposure_caches_valid = true;
            }
            auto personal_rng = PersonalRandomNumberGenerator(m.get_rng(), m.get_person(person_id));
            mio::abm::interact(personal_rng, m.get_person(person_id), m.get_location(person_id),
                               m_air_exposure_rates_cache[m.get_location(person_id).get_id().get()],
                               m_contact_exposure_rates_cache[m.get_location(person_id).get_id().get()], t, dt,
                               m.parameters);
        }
    }

    void perform_mobility(Model& m, TimePoint t, TimeSpan dt)
    {
        const uint32_t num_persons = static_cast<uint32_t>(m.get_persons().size());
        PRAGMA_OMP(parallel for)
        for (uint32_t i = 0; i < num_persons; ++i) {
            const PersonId person_id(i);
            Person& person    = m.get_person(person_id);
            auto personal_rng = PersonalRandomNumberGenerator(m.get_rng(), person);

            auto try_mobility_rule = [&](auto rule) -> bool {
                // run mobility rule and check if change of location can actually happen
                auto target_type                  = rule(personal_rng, person, t, dt, m.parameters);
                const Location& target_location   = m.get_location(m.find_location(target_type, person_id));
                const LocationId current_location = person.get_location();

                // the Person cannot move if they do not wear mask as required at targeted location
                if (target_location.is_mask_required() && !person.is_compliant(personal_rng, InterventionType::Mask)) {
                    return false;
                }
                // the Person cannot move if the capacity of targeted Location is reached
                if (target_location.get_id() == current_location ||
                    get_number_persons(m, target_location.get_id()) >= target_location.get_capacity().persons) {
                    return false;
                }
                // the Person cannot move if the performed TestingStrategy is positive
                if (!m.get_testing_strategy().run_strategy(personal_rng, person, target_location, t)) {
                    return false;
                }
                // update worn mask to target location's requirements
                if (target_location.is_mask_required()) {
                    // if the current MaskProtection level is lower than required, the Person changes mask
                    if (m.parameters.get<MaskProtection>()[person.get_mask().get_type()] <
                        m.parameters.get<MaskProtection>()[target_location.get_required_mask()]) {
                        person.set_mask(target_location.get_required_mask(), t);
                    }
                }
                else {
                    person.set_mask(MaskType::None, t);
                }
                // all requirements are met, move to target location
                m.change_location(person_id, target_location.get_id());
                return true;
            };

            // run mobility rules one after the other if the corresponding location type exists
            // shortcutting of bool operators ensures the rules stop after the first rule is applied
            if (m.use_mobility_rules()) {
                (m.has_locations({LocationType::Cemetery}) && try_mobility_rule(&get_buried)) ||
                    (m.has_locations({LocationType::Home}) && try_mobility_rule(&return_home_when_recovered)) ||
                    (m.has_locations({LocationType::Hospital}) && try_mobility_rule(&go_to_hospital)) ||
                    (m.has_locations({LocationType::ICU}) && try_mobility_rule(&go_to_icu)) ||
                    (m.has_locations({LocationType::School, LocationType::Home}) && try_mobility_rule(&go_to_school)) ||
                    (m.has_locations({LocationType::Work, LocationType::Home}) && try_mobility_rule(&go_to_work)) ||
                    (m.has_locations({LocationType::BasicsShop, LocationType::Home}) &&
                     try_mobility_rule(&go_to_shop)) ||
                    (m.has_locations({LocationType::SocialEvent, LocationType::Home}) &&
                     try_mobility_rule(&go_to_event)) ||
                    (m.has_locations({LocationType::Home}) && try_mobility_rule(&go_to_quarantine));
            }
            else {
                // no daily routine mobility, just infection related
                (m.has_locations({LocationType::Cemetery}) && try_mobility_rule(&get_buried)) ||
                    (m.has_locations({LocationType::Home}) && try_mobility_rule(&return_home_when_recovered)) ||
                    (m.has_locations({LocationType::Hospital}) && try_mobility_rule(&go_to_hospital)) ||
                    (m.has_locations({LocationType::ICU}) && try_mobility_rule(&go_to_icu)) ||
                    (m.has_locations({LocationType::Home}) && try_mobility_rule(&go_to_quarantine));
            }
        }

        // check if a person makes a trip
        bool weekend     = t.is_weekend();
        size_t num_trips = m.get_trip_list().num_trips(weekend);

        for (; m.get_trip_list().get_current_index() < num_trips &&
               m.get_trip_list().get_next_trip_time(weekend).seconds() < (t + dt).time_since_midnight().seconds();
             m.get_trip_list().increase_index()) {
            auto& trip        = m.get_trip_list().get_next_trip(weekend);
            auto& person      = m.get_person(trip.person_id);
            auto personal_rng = PersonalRandomNumberGenerator(m.get_rng(), person);
            // skip the trip if the person is in quarantine or is dead
            if (person.is_in_quarantine(t, m.parameters) || person.get_infection_state(t) == InfectionState::Dead) {
                continue;
            }
            auto& target_location = m.get_location(trip.destination);
            // skip the trip if the Person wears mask as required at targeted location
            if (target_location.is_mask_required() && !person.is_compliant(personal_rng, InterventionType::Mask)) {
                continue;
            }
            // skip the trip if the performed TestingStrategy is positive
            if (!m.get_testing_strategy().run_strategy(personal_rng, person, target_location, t)) {
                continue;
            }
            // all requirements are met, move to target location
            m.change_location(person.get_id(), target_location.get_id(), trip.trip_mode);
            // update worn mask to target location's requirements
            if (target_location.is_mask_required()) {
                // if the current MaskProtection level is lower than required, the Person changes mask
                if (m.parameters.get<MaskProtection>()[person.get_mask().get_type()] <
                    m.parameters.get<MaskProtection>()[target_location.get_required_mask()]) {
                    person.set_mask(target_location.get_required_mask(), t);
                }
            }
            else {
                person.set_mask(MaskType::None, t);
            }
        }
        if (((t).days() < std::floor((t + dt).days()))) {
            m.get_trip_list().reset_index();
        }
    }

    void build_compute_local_population_cache(const Model& m) const
    {
        PRAGMA_OMP(single)
        {
            const size_t num_locations = m.get_locations().size();
            const size_t num_persons   = m.get_persons().size();
            m_local_population_cache.resize(num_locations);
            PRAGMA_OMP(taskloop)
            for (size_t i = 0; i < num_locations; i++) {
                m_local_population_cache[i] = 0;
            } // implicit taskloop barrier
            PRAGMA_OMP(taskloop)
            for (size_t i = 0; i < num_persons; i++) {
                ++m_local_population_cache[m.get_persons()[i].get_location().get()];
            } // implicit taskloop barrier
        } // implicit single barrier
    }

    void build_exposure_caches(const Model& m)
    {
        PRAGMA_OMP(single)
        {
            const size_t num_locations = m.get_locations().size();
            m_air_exposure_rates_cache.resize(num_locations);
            m_contact_exposure_rates_cache.resize(num_locations);
            PRAGMA_OMP(taskloop)
            for (size_t i = 0; i < num_locations; i++) {
                m_air_exposure_rates_cache[i].resize(
                    {CellIndex(m.get_locations()[i].get_cells().size()), VirusVariant::Count});
                m_contact_exposure_rates_cache[i].resize({CellIndex(m.get_locations()[i].get_cells().size()),
                                                          VirusVariant::Count,
                                                          AgeGroup(m.parameters.get_num_groups())});
            } // implicit taskloop barrier
            m_are_exposure_caches_valid    = false;
            m_exposure_caches_need_rebuild = false;
        } // implicit single barrier
    }

    void compute_exposure_caches(const Model& m, TimePoint t, TimeSpan dt)
    {
        PRAGMA_OMP(single)
        {
            // if cache shape was changed (e.g. by add_location), rebuild it
            if (m_exposure_caches_need_rebuild) {
                build_exposure_caches(m);
            }
            // use these const values to help omp recognize that the for loops are bounded
            const auto num_locations = m.get_locations().size();
            const auto num_persons   = m.get_persons().size();

            // 1) reset all cached values
            // Note: we cannot easily reuse values, as they are time dependant (get_infection_state)
            PRAGMA_OMP(taskloop)
            for (size_t i = 0; i < num_locations; ++i) {
                const auto index         = i;
                auto& local_air_exposure = m_air_exposure_rates_cache[index];
                std::for_each(local_air_exposure.begin(), local_air_exposure.end(), [](auto& r) {
                    r = 0.0;
                });
                auto& local_contact_exposure = m_contact_exposure_rates_cache[index];
                std::for_each(local_contact_exposure.begin(), local_contact_exposure.end(), [](auto& r) {
                    r = 0.0;
                });
            } // implicit taskloop barrier
            // here is an implicit (and needed) barrier from parallel for

            // 2) add all contributions from each person
            PRAGMA_OMP(taskloop)
            for (size_t i = 0; i < num_persons; ++i) {
                const Person& person = m.get_persons()[i];
                const auto location  = person.get_location().get();
                mio::abm::add_exposure_contribution(m_air_exposure_rates_cache[location],
                                                    m_contact_exposure_rates_cache[location], person,
                                                    m.get_location(person.get_id()), t, dt);
            } // implicit taskloop barrier
        } // implicit single barrier
    }

    void begin_step(Model& m, TimePoint t, TimeSpan dt)
    {
        m.get_testing_strategy().update_activity_status(t);

        if (!m_is_local_population_cache_valid) {
            build_compute_local_population_cache(m);
            m_is_local_population_cache_valid = true;
        }
        compute_exposure_caches(m, t, dt);
        m_are_exposure_caches_valid = true;
    }

    size_t get_number_persons(const Model& m, LocationId location) const
    {
        if (!m_is_local_population_cache_valid) {
            build_compute_local_population_cache(m);
        }
        return m_local_population_cache[location.get()];
    }

    mutable Eigen::Matrix<std::atomic_int_fast32_t, Eigen::Dynamic, 1>
        m_local_population_cache; ///< Current number of Persons in a given location.
    Eigen::Matrix<AirExposureRates, Eigen::Dynamic, 1>
        m_air_exposure_rates_cache; ///< Cache for local exposure through droplets in #transmissions/day.
    Eigen::Matrix<ContactExposureRates, Eigen::Dynamic, 1>
        m_contact_exposure_rates_cache; ///< Cache for local exposure through contacts in #transmissions/day.
    bool m_is_local_population_cache_valid = false;
    bool m_are_exposure_caches_valid       = false;
    bool m_exposure_caches_need_rebuild    = true;
};

} // namespace abm
} // namespace mio
