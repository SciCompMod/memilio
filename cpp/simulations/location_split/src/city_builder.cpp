#include "city_builder.h"
#include "constants.h"
#include "parameter_setter.h"

#include "abm/location_type.h"
#include "memilio/utils/random_number_generator.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

namespace
{
/// @brief Index of an AgeGroup in the six age groups of this simulation.
inline size_t age_index(mio::AgeGroup age)
{
    return (size_t)age.get();
}
} // namespace

mio::abm::Model CityBuilder::build_model(const CityConfig& config, const mio::RandomNumberGenerator& rng)
{
    auto model      = mio::abm::Model(num_age_groups);
    model.get_rng() = rng;

    const auto infra = config.infrastructure();

    const auto hospitals = create_locations(model, mio::abm::LocationType::Hospital, infra.num_hospitals);
    const auto icus      = create_locations(model, mio::abm::LocationType::ICU, infra.num_icus);
    const auto households =
        create_locations(model, mio::abm::LocationType::Home,
                         std::accumulate(infra.num_households_hh_size.begin(), infra.num_households_hh_size.end(), 0));
    const auto workplaces   = create_locations(model, mio::abm::LocationType::Work, infra.num_workplaces);
    const auto prim_schools = create_locations(model, mio::abm::LocationType::School, infra.num_elementary_schools);
    const auto sec_schools  = create_locations(model, mio::abm::LocationType::School, infra.num_secondary_schools);
    const auto shops        = create_locations(model, mio::abm::LocationType::BasicsShop, infra.num_stores);
    const auto events       = create_locations(model, mio::abm::LocationType::SocialEvent, infra.num_events);

    create_and_assign_people(model, config, infra, households, workplaces, prim_schools, sec_schools, shops, events,
                             hospitals, icus);

    // The contact rates have to be set after the Location%s exist.
    set_local_parameters(model, config.region);
    set_parameters(model.parameters);

    return model;
}

std::vector<mio::abm::LocationId> CityBuilder::create_locations(mio::abm::Model& model, mio::abm::LocationType type,
                                                                int count)
{
    std::vector<mio::abm::LocationId> ids;
    ids.reserve(static_cast<size_t>(std::max(0, count)));
    for (int i = 0; i < count; ++i) {
        ids.push_back(model.add_location(type));
    }
    return ids;
}

void CityBuilder::create_and_assign_people(mio::abm::Model& model, const CityConfig& config,
                                           const CityParameters::CityInfrastructure& infra,
                                           const std::vector<mio::abm::LocationId>& households,
                                           const std::vector<mio::abm::LocationId>& workplaces,
                                           const std::vector<mio::abm::LocationId>& prim_schools,
                                           const std::vector<mio::abm::LocationId>& sec_schools,
                                           const std::vector<mio::abm::LocationId>& shops,
                                           const std::vector<mio::abm::LocationId>& events,
                                           const std::vector<mio::abm::LocationId>& hospitals,
                                           const std::vector<mio::abm::LocationId>& icus)
{
    // Persons are created household by household and then assigned to the remaining location types
    // following these rules:
    //  1. Every household that contains a person of age group 0-4 or 5-14 also contains at least one
    //     person of an older age group.
    //  2. The elementary schools are filled with persons of age group 5-14 first.
    //  3. The secondary schools are filled with the remaining persons of age group 5-14 and then with
    //     persons of age group 15-34.
    //  4. The workplaces are filled with the remaining persons of age groups 15-34 and 35-59, and with
    //     persons of age group 60+ if any workplaces are left.
    //  5. Everybody is assigned to the single hospital and the single ICU.
    //  6. Everybody is assigned to a social event, everybody outside of age group 0-4 also to a shop.

    auto& model_rng      = model.get_rng();
    const auto params    = config.region_parameters();
    const auto age_count = CityParameters::age_vector(config.total_population, params);

    std::vector<std::vector<mio::AgeGroup>> temp_households;
    temp_households.reserve(households.size());
    std::vector<int> remaining_ages = age_count;

    const int max_household_size = static_cast<int>(infra.num_households_hh_size.size());
    for (int hh_size = 1; hh_size <= max_household_size; ++hh_size) {
        const int num_households = infra.num_households_hh_size[hh_size - 1];

        for (int hh = 0; hh < num_households; ++hh) {
            std::vector<mio::AgeGroup> household;
            bool needs_parent = false;

            for (int person = 0; person < hh_size; ++person) {
                std::vector<int> available_ages;
                for (int age = 0; age < (int)remaining_ages.size(); ++age) {
                    if (remaining_ages[age] > 0) {
                        available_ages.push_back(age);
                    }
                }
                if (available_ages.empty()) {
                    break;
                }

                const auto draw = mio::UniformIntDistribution<size_t>::get_instance()(model_rng, size_t(0),
                                                                                      available_ages.size() - 1);
                const int selected_age = available_ages[draw];

                household.push_back(static_cast<mio::AgeGroup>(selected_age));
                remaining_ages[selected_age]--;

                if (selected_age == 0 || selected_age == 1) {
                    needs_parent = true;
                }
            }

            // A household of children only gets one of its children replaced by an adult.
            if (needs_parent) {
                const bool has_adult = std::any_of(household.begin(), household.end(), [](mio::AgeGroup age) {
                    return age != age_group_0_to_4 && age != age_group_5_to_14;
                });

                if (!has_adult) {
                    for (size_t i = 0; i < household.size(); ++i) {
                        const auto current_age = household[i];
                        if (current_age == age_group_0_to_4 || current_age == age_group_5_to_14) {
                            remaining_ages[age_index(current_age)]++;

                            std::vector<int> adult_ages;
                            for (int age = 2; age < (int)remaining_ages.size(); ++age) {
                                if (remaining_ages[age] > 0) {
                                    adult_ages.push_back(age);
                                }
                            }
                            if (!adult_ages.empty()) {
                                const auto adult_draw = mio::UniformIntDistribution<size_t>::get_instance()(
                                    model_rng, size_t(0), adult_ages.size() - 1);
                                const int selected_adult = adult_ages[adult_draw];
                                household[i]             = static_cast<mio::AgeGroup>(selected_adult);
                                remaining_ages[selected_adult]--;
                            }
                            else {
                                // No adult left in the pool, put the child back.
                                remaining_ages[age_index(current_age)]--;
                            }
                            break;
                        }
                    }
                }
            }

            temp_households.push_back(household);
        }
    }

    // Shuffle so that the household size is not correlated with the order the persons are created in.
    std::shuffle(temp_households.begin(), temp_households.end(), model_rng);

    std::vector<int> persons_per_household(households.size(), 0);
    for (size_t household_index = 0; household_index < temp_households.size(); ++household_index) {
        const auto household_id = households[household_index];
        for (const auto& age : temp_households[household_index]) {
            const auto person_id = model.add_person(household_id, age);
            model.assign_location(person_id, household_id);
            model.assign_location(person_id, hospitals[0]);
            model.assign_location(person_id, icus[0]);
            persons_per_household[household_index]++;
        }
    }

    // Sanity check on the realized household size distribution.
    std::vector<int> household_sizes(infra.num_households_hh_size.size(), 0);
    for (const auto count : persons_per_household) {
        if (count >= 1 && count <= max_household_size) {
            household_sizes[count - 1]++;
        }
    }
    for (size_t i = 0; i < household_sizes.size(); ++i) {
        if (household_sizes[i] != infra.num_households_hh_size[i]) {
            std::cerr << "Warning: expected " << infra.num_households_hh_size[i] << " households of size " << (i + 1)
                      << ", but built " << household_sizes[i] << ".\n";
        }
    }

    // Schools are filled with the 5-14 year olds first. If there are not enough of them, the
    // remaining places are taken by 15-34 year olds.
    const int number_of_5_to_14_year_olds = age_count[1];
    int prim_overhang                     = 0;
    int sec_overhang                      = 0;
    if (number_of_5_to_14_year_olds >= infra.num_persons_elementary_schools) {
        if (number_of_5_to_14_year_olds <
            infra.num_persons_secondary_schools + infra.num_persons_elementary_schools) {
            sec_overhang = infra.num_persons_secondary_schools + infra.num_persons_elementary_schools -
                           number_of_5_to_14_year_olds;
        }
    }
    else {
        prim_overhang = infra.num_persons_elementary_schools - number_of_5_to_14_year_olds;
        sec_overhang  = infra.num_persons_secondary_schools;
    }

    int amount_of_60plus_that_work =
        std::max(0, infra.num_worker - sec_overhang - age_count[2] - age_count[3]);

    int count_elementary_schools = 0;
    int count_secondary_schools  = 0;
    int count_shops              = 0;
    int count_events             = 0;
    int count_workers            = 0;

    for (auto& person : model.get_persons()) {
        const auto person_id = person.get_id();
        const auto age       = person.get_age();

        if (age != age_group_0_to_4) {
            model.assign_location(person_id, shops[count_shops % shops.size()]);
            count_shops++;
        }
        model.assign_location(person_id, events[count_events % events.size()]);
        count_events++;

        if (age == age_group_5_to_14) {
            if (count_elementary_schools < infra.num_persons_elementary_schools) {
                model.assign_location(person_id, prim_schools[count_elementary_schools % prim_schools.size()]);
                count_elementary_schools++;
                continue;
            }
            if (count_secondary_schools < infra.num_persons_secondary_schools) {
                model.assign_location(person_id, sec_schools[count_secondary_schools % sec_schools.size()]);
                count_secondary_schools++;
                continue;
            }
        }
        else if ((prim_overhang > 0 || sec_overhang > 0) && age == age_group_15_to_34) {
            if (prim_overhang > 0) {
                model.assign_location(person_id, prim_schools[count_elementary_schools % prim_schools.size()]);
                count_elementary_schools++;
                prim_overhang--;
                continue;
            }
            model.assign_location(person_id, sec_schools[count_secondary_schools % sec_schools.size()]);
            count_secondary_schools++;
            sec_overhang--;
            continue;
        }
        else if (age == age_group_35_to_59 || age == age_group_15_to_34) {
            if (count_workers < infra.num_worker) {
                model.assign_location(person_id, workplaces[count_workers % workplaces.size()]);
                count_workers++;
                continue;
            }
        }
        else if (age == age_group_60_to_79 || age == age_group_80_plus) {
            if (amount_of_60plus_that_work > 0) {
                model.assign_location(person_id, workplaces[count_workers % workplaces.size()]);
                count_workers++;
                amount_of_60plus_that_work--;
                continue;
            }
        }
    }

    // Person%s that have no Location of some type are detached from that type explicitly.
    //
    // Model::perform_mobility only checks whether the *model id* assigned for the target LocationType
    // matches the Model, and that id defaults to 0, which is also the id of this Model. A Person
    // without, say, a Work would therefore pass that check and the Model would dereference an invalid
    // LocationId. Setting the model id of the unused types to an id no Model has makes the mobility
    // rules skip them.
    constexpr int detached_model_id = -1;
    for (auto& person : model.get_persons()) {
        for (uint32_t type = 0; type < (uint32_t)mio::abm::LocationType::Count; ++type) {
            const auto location_type = mio::abm::LocationType(type);
            if (person.get_assigned_location(location_type) == mio::abm::LocationId::invalid_id()) {
                person.set_assigned_location(location_type, mio::abm::LocationId::invalid_id(),
                                             detached_model_id);
            }
        }
    }

    // Sanity checks on the assignment.
    int assigned_workers  = 0;
    int school_attendees  = 0;
    const auto invalid_id = mio::abm::LocationId::invalid_id();

    for (const auto& person : model.get_persons()) {
        const bool has_school = person.get_assigned_location(mio::abm::LocationType::School) != invalid_id;
        const bool has_work   = person.get_assigned_location(mio::abm::LocationType::Work) != invalid_id;

        school_attendees += has_school ? 1 : 0;
        assigned_workers += has_work ? 1 : 0;

        if (person.get_assigned_location(mio::abm::LocationType::Home) == invalid_id) {
            std::cerr << "Error: person " << person.get_id().get() << " has no assigned home.\n";
        }
        if (has_school && has_work) {
            std::cerr << "Error: person " << person.get_id().get() << " has both a work and a school assigned.\n";
        }
        if (person.get_assigned_location(mio::abm::LocationType::SocialEvent) == invalid_id) {
            std::cerr << "Error: person " << person.get_id().get() << " has no assigned social event.\n";
        }
        if (person.get_assigned_location(mio::abm::LocationType::BasicsShop) == invalid_id &&
            person.get_age() != age_group_0_to_4) {
            std::cerr << "Error: person " << person.get_id().get() << " has no assigned shop.\n";
        }
        if (person.get_assigned_location(mio::abm::LocationType::Hospital) == invalid_id) {
            std::cerr << "Error: person " << person.get_id().get() << " has no assigned hospital.\n";
        }
        if (person.get_assigned_location(mio::abm::LocationType::ICU) == invalid_id) {
            std::cerr << "Error: person " << person.get_id().get() << " has no assigned ICU.\n";
        }
    }

    if (assigned_workers != infra.num_worker) {
        std::cerr << "Warning: expected " << infra.num_worker << " workers, but assigned " << assigned_workers
                  << ".\n";
    }
    const int expected_students = infra.num_persons_elementary_schools + infra.num_persons_secondary_schools;
    if (school_attendees != expected_students) {
        std::cerr << "Warning: expected " << expected_students << " students, but assigned " << school_attendees
                  << ".\n";
    }
}

void CityBuilder::print_city_summary(const CityConfig& config)
{
    const auto params = config.region_parameters();
    const auto infra  = config.infrastructure();
    const auto ages   = CityParameters::age_vector(config.total_population, params);

    const int total_households =
        std::accumulate(infra.num_households_hh_size.begin(), infra.num_households_hh_size.end(), 0);
    const int total_locations = total_households + infra.num_workplaces + infra.num_elementary_schools +
                                infra.num_secondary_schools + infra.num_hospitals + infra.num_icus +
                                infra.num_stores + infra.num_events;

    static const char* const age_labels[] = {"0-4  ", "5-14 ", "15-34", "35-59", "60-79", "80+  "};

    std::cout << "\n" << std::string(50, '=') << "\n";
    std::cout << "  CITY SUMMARY (" << CityParameters::region_to_string(config.region) << ")\n";
    std::cout << std::string(50, '=') << "\n\n";

    std::cout << "POPULATION\n";
    std::cout << "  Total population: " << config.total_population << "\n";
    for (size_t i = 0; i < params.age_distribution.size(); ++i) {
        std::cout << "  " << age_labels[i] << " years: " << std::fixed << std::setprecision(1)
                  << (params.age_distribution[i] * 100) << "% (" << ages[i] << " people)\n";
    }

    std::cout << "\nHOUSING\n";
    std::cout << "  Households: " << total_households << "\n";
    std::cout << "  Average household size: " << std::setprecision(2)
              << static_cast<double>(config.total_population) / total_households << "\n";
    for (size_t i = 0; i < infra.num_households_hh_size.size(); ++i) {
        std::cout << "  " << (i + 1) << " person: " << infra.num_households_hh_size[i] << "\n";
    }

    std::cout << "\nEMPLOYMENT\n";
    std::cout << "  Workplaces: " << infra.num_workplaces << "\n";
    std::cout << "  Workers: " << infra.num_worker << "\n";
    std::cout << "  Average employees per workplace: " << std::setprecision(1)
              << static_cast<double>(infra.num_worker) / infra.num_workplaces << "\n";

    std::cout << "\nEDUCATION\n";
    std::cout << "  Elementary schools: " << infra.num_elementary_schools << " ("
              << infra.num_persons_elementary_schools << " students)\n";
    std::cout << "  Secondary schools: " << infra.num_secondary_schools << " ("
              << infra.num_persons_secondary_schools << " students)\n";

    std::cout << "\nHEALTHCARE, RETAIL AND SOCIAL\n";
    std::cout << "  Hospitals: " << infra.num_hospitals << ", ICUs: " << infra.num_icus << "\n";
    std::cout << "  Stores: " << infra.num_stores << "\n";
    std::cout << "  Social event venues: " << infra.num_events << "\n";

    std::cout << "\nINFRASTRUCTURE DENSITY\n";
    std::cout << "  Total locations: " << total_locations << "\n";
    std::cout << "  Locations per 1000 people: " << std::setprecision(1)
              << static_cast<double>(total_locations) * 1000.0 / config.total_population << "\n";

    std::cout << "\n" << std::string(50, '=') << "\n\n";
}

void CityBuilder::save_city_to_file(const CityConfig& config, const std::string& filename)
{
    const auto params = config.region_parameters();
    const auto infra  = config.infrastructure();
    const auto ages   = CityParameters::age_vector(config.total_population, params);

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Error: could not open " << filename << " for writing.\n";
        return;
    }

    static const char* const age_keys[] = {"0_to_4", "5_to_14", "15_to_34", "35_to_59", "60_to_79", "80_plus"};

    ofs << "key,value\n";
    ofs << "region," << CityParameters::region_to_string(config.region) << "\n";

    // Population.
    ofs << "total_population," << config.total_population << "\n";
    for (size_t i = 0; i < params.age_distribution.size(); ++i) {
        ofs << "age_group_" << age_keys[i] << "," << ages[i] << "\n";
    }
    for (size_t i = 0; i < params.age_distribution.size(); ++i) {
        ofs << "age_percentage_" << age_keys[i] << "," << std::fixed << std::setprecision(2)
            << (params.age_distribution[i] * 100) << "\n";
    }

    // Households.
    const int total_households =
        std::accumulate(infra.num_households_hh_size.begin(), infra.num_households_hh_size.end(), 0);
    ofs << "total_households," << total_households << "\n";
    ofs << "average_household_size," << std::setprecision(2)
        << static_cast<double>(config.total_population) / total_households << "\n";
    for (size_t i = 0; i < infra.num_households_hh_size.size(); ++i) {
        ofs << "households_size_" << (i + 1) << "," << infra.num_households_hh_size[i] << "\n";
        ofs << "households_size_" << (i + 1) << "_percentage," << std::setprecision(2)
            << (static_cast<double>(infra.num_households_hh_size[i]) / total_households * 100) << "\n";
    }

    // Employment.
    ofs << "num_workplaces," << infra.num_workplaces << "\n";
    ofs << "num_workers," << infra.num_worker << "\n";
    ofs << "employment_rate," << std::setprecision(3) << params.ratios.employment_rate << "\n";
    ofs << "average_employees_per_workplace," << std::setprecision(1)
        << static_cast<double>(infra.num_worker) / infra.num_workplaces << "\n";

    // Education.
    ofs << "num_elementary_schools," << infra.num_elementary_schools << "\n";
    ofs << "num_secondary_schools," << infra.num_secondary_schools << "\n";
    ofs << "elementary_school_students," << infra.num_persons_elementary_schools << "\n";
    ofs << "secondary_school_students," << infra.num_persons_secondary_schools << "\n";
    ofs << "school_age_population_5_14," << ages[1] << "\n";

    // Healthcare, retail, social.
    ofs << "num_hospitals," << infra.num_hospitals << "\n";
    ofs << "num_icus," << infra.num_icus << "\n";
    ofs << "num_stores," << infra.num_stores << "\n";
    ofs << "shoppers," << (config.total_population - ages[0]) << "\n";
    ofs << "num_events," << infra.num_events << "\n";

    // Density.
    const int total_locations = total_households + infra.num_workplaces + infra.num_elementary_schools +
                                infra.num_secondary_schools + infra.num_hospitals + infra.num_icus +
                                infra.num_stores + infra.num_events;
    ofs << "total_locations," << total_locations << "\n";
    ofs << "locations_per_1000_people," << std::setprecision(1)
        << static_cast<double>(total_locations) * 1000.0 / config.total_population << "\n";

    std::cout << "City statistics saved to " << filename << "\n";
}
