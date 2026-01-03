#include "../include/city_builder.h"
#include "../include/constants.h"
#include "../include/parameter_setter.h"
#include "../include/city_parameters.h"
#include <random>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <iostream>
#include <iomanip>

mio::abm::World CityBuilder::build_world(const CityConfig& config, const mio::RandomNumberGenerator& rng)
{
    auto world      = mio::abm::World(num_age_groups);
    world.get_rng() = rng; // Set the random number generator for the world

    // Get infrastructure configuration based on German demographic data
    auto infra = config.infrastructure();

    // Healthcare facilities
    std::vector<mio::abm::LocationId> hospitals;
    std::vector<mio::abm::LocationId> icus;

    for (int i = 0; i < infra.num_hospitals; ++i) {
        hospitals.push_back(world.add_location(mio::abm::LocationType::Hospital));
    }
    for (int i = 0; i < infra.num_icus; ++i) {
        icus.push_back(world.add_location(mio::abm::LocationType::ICU));
    }

    // Create all location types based on German infrastructure ratios
    auto households = create_households(
        world, std::accumulate(infra.num_households_hh_size.begin(), infra.num_households_hh_size.end(), 0));
    auto workplaces   = create_workplaces(world, infra.num_workplaces);
    auto prim_schools = create_schools(world, infra.num_elementary_schools);
    auto sec_schools  = create_schools(world, infra.num_secondary_schools);
    auto shops        = create_shops(world, infra.num_stores);
    auto events       = create_events(world, infra.num_events);

    // Assign people to locations using German demographic distribution
    create_and_assign_people_to_locations(world, households, workplaces, prim_schools, sec_schools, shops, events,
                                          hospitals, icus, config.total_population, infra.num_households_hh_size,
                                          infra);

    set_local_parameters_ger(world);
    set_parameters(world.parameters);

    return world;
}

std::vector<mio::abm::LocationId> CityBuilder::create_households(mio::abm::World& world, int num_households)
{
    std::vector<mio::abm::LocationId> household_ids;
    household_ids.reserve(num_households);

    for (int i = 0; i < num_households; ++i) {
        auto household_id = world.add_location(mio::abm::LocationType::Home);
        household_ids.push_back(household_id);
    }

    return household_ids;
}

std::vector<mio::abm::LocationId> CityBuilder::create_workplaces(mio::abm::World& world, int num_workplaces)
{
    std::vector<mio::abm::LocationId> workplace_ids;
    workplace_ids.reserve(num_workplaces);

    for (int i = 0; i < num_workplaces; ++i) {
        auto workplace_id = world.add_location(mio::abm::LocationType::Work);
        workplace_ids.push_back(workplace_id);
    }

    return workplace_ids;
}

std::vector<mio::abm::LocationId> CityBuilder::create_schools(mio::abm::World& world, int num_schools)
{
    std::vector<mio::abm::LocationId> school_ids;
    school_ids.reserve(num_schools);

    for (int i = 0; i < num_schools; ++i) {
        auto school_id = world.add_location(mio::abm::LocationType::School);
        school_ids.push_back(school_id);
    }

    return school_ids;
}

std::vector<mio::abm::LocationId> CityBuilder::create_events(mio::abm::World& world, int num_events)
{
    std::vector<mio::abm::LocationId> event_ids;
    event_ids.reserve(num_events);

    // Create large events (concert halls, stadiums)
    for (int i = 0; i < num_events; ++i) {
        auto event_id = world.add_location(mio::abm::LocationType::SocialEvent);
        event_ids.push_back(event_id);
    }

    return event_ids;
}

std::vector<mio::abm::LocationId> CityBuilder::create_shops(mio::abm::World& world, int num_shops)
{
    std::vector<mio::abm::LocationId> shop_ids;
    shop_ids.reserve(num_shops);

    // Create grocery stores
    for (int i = 0; i < num_shops; ++i) {
        auto shop_id = world.add_location(mio::abm::LocationType::BasicsShop);
        shop_ids.push_back(shop_id);
    }

    return shop_ids;
}

std::vector<int> CityBuilder::create_age_vector(int total_population)
{
    std::vector<int> age_vector;
    age_vector.resize(num_age_groups);

    // Use German age distribution from 2023 census data
    for (size_t i = 0; i < CityParameters::AGE_DISTRIBUTION.size(); ++i) {
        int age_group_population = static_cast<int>(total_population * CityParameters::AGE_DISTRIBUTION[i]);
        age_vector.at(i)         = age_group_population;
    }

    // see if total population matches the expected distribution
    int total_age_population = std::accumulate(age_vector.begin(), age_vector.end(), 0);
    if (total_age_population != total_population) {
        std::cerr << "Warning: Total population does not match expected distribution. "
                  << "Expected: " << total_population << ", Actual: " << total_age_population << "\n";
    }

    return age_vector;
}

int CityBuilder::ageGroupTInt6(mio::AgeGroup age_group)
{
    if (age_group == age_group_0_to_4) {
        return 0;
    }
    else if (age_group == age_group_5_to_14) {
        return 1;
    }
    else if (age_group == age_group_15_to_34) {
        return 2;
    }
    else if (age_group == age_group_35_to_59) {
        return 3;
    }
    else if (age_group == age_group_60_to_79) {
        return 4;
    }
    else if (age_group == age_group_80_plus) {
        return 5;
    }
    else {
        return -1;
    }
}

void CityBuilder::create_and_assign_people_to_locations(
    mio::abm::World& world, const std::vector<mio::abm::LocationId>& household_locations,
    const std::vector<mio::abm::LocationId>& workplaces, const std::vector<mio::abm::LocationId>& prim_schools,
    const std::vector<mio::abm::LocationId>& sec_schools, const std::vector<mio::abm::LocationId>& shops,
    const std::vector<mio::abm::LocationId>& events, const std::vector<mio::abm::LocationId>& hospitals,
    const std::vector<mio::abm::LocationId>& icus, int total_population, std::vector<int>& hh_per_size,
    const CityParameters::CityInfrastructure infra)
{

    // We need to create and assign persons.
    // We just want to create them (with an age) and assign them (to a home, perhaps a school, perhaps work, the only icu, the only hospital, one event, one shop.

    // With these rules:

    //    We have 5 age groups:
    //    0: 0-4 years
    //    1: 5-14 years
    //    2: 15-34 years
    //    3: 35-64 years
    //    4: 65+ years
    // 1.
    // Every household with persons in age group 0 and 1 needs to be paired with at least one other person in an age group 2 or 3. (At least but randomly.)
    // 2.
    // We assign maximum amount of num_elementary_schools to persons in age group 1
    // 3.
    // We assign maximum amount of num_sec_schools to persons first in age group 1 and then to age group 2
    // 4.
    // We assign the amount of working people randomly between the remaining persons in age group 2 and group 3, if there are any left, we assign them to group 4
    // 5.
    // Everybody gets assigned to the same ICU and Hospital.
    // 6.
    // Everybody gets assigned to a random social event and everybody whos not in age group 0 to a random Basic shop.

    auto world_rng  = world.get_rng();
    auto age_vector = create_age_vector(total_population);

    std::vector<std::vector<mio::AgeGroup>> temp_households;
    std::vector<int> remaining_ages = age_vector;

    for (int hh_size = 1; hh_size <= 5; ++hh_size) {
        int num_households = hh_per_size[hh_size - 1];

        for (int hh = 0; hh < num_households; ++hh) {
            std::vector<mio::AgeGroup> household;

            // Check if we need to enforce the child-parent rule
            bool needs_parent = false;

            // First, randomly assign ages to fill the household
            for (int person = 0; person < hh_size; ++person) {
                // Find available age groups
                std::vector<int> available_ages;
                for (int age = 0; age < (int)remaining_ages.size(); ++age) {
                    if (remaining_ages[age] > 0) {
                        available_ages.push_back(age);
                    }
                }

                if (available_ages.empty())
                    break;

                // Randomly select an age group
                auto age         = mio::UniformIntDistribution<size_t>::get_instance()(world_rng, size_t(0),
                                                                               available_ages.size() - 1);
                int selected_age = available_ages[age];

                household.push_back(static_cast<mio::AgeGroup>(selected_age));
                remaining_ages[selected_age]--;

                // Check if this person is a child (age 0 or 1)
                if (selected_age == 0 || selected_age == 1) {
                    needs_parent = true;
                }
            }

            // If household has children but no adults, replace one person with an adult
            if (needs_parent) {
                bool has_adult = false;
                for (const auto& age : household) {
                    if (age != age_group_0_to_4 && age != age_group_5_to_14) {
                        has_adult = true;
                        break;
                    }
                }

                if (!has_adult) {
                    // Find a non-child to replace
                    for (int i = 0; i < (int)household.size(); ++i) {
                        auto current_age = household[i];
                        if ((current_age == age_group_0_to_4) || (current_age == age_group_5_to_14)) {
                            // Return this age to the pool
                            remaining_ages[ageGroupTInt6(current_age)]++;

                            // Find an available adult age (2-5)
                            std::vector<int> adult_ages;
                            for (int age = 2; age <= 5; ++age) {
                                if (remaining_ages[age] > 0) {
                                    adult_ages.push_back(age);
                                }
                            }

                            if (!adult_ages.empty()) {
                                auto adult_dist = mio::UniformIntDistribution<size_t>::get_instance()(
                                    world_rng, size_t(0), adult_ages.size() - 1);
                                int selected_adult = adult_ages[adult_dist];
                                household[i]       = static_cast<mio::AgeGroup>(selected_adult);
                                remaining_ages[selected_adult]--;
                            }
                            break;
                        }
                    }
                }
            }

            temp_households.push_back(household);
        }
    }
    // We shuffle the households to ensure randomness in assignment
    std::shuffle(temp_households.begin(), temp_households.end(), world_rng);
    // We now have a vector of households, each containing a vector of AgeGroups.
    // Each household is a vector of AgeGroups, where each AgeGroup corresponds to a person in that household.
    // For example, a household might look like this:
    // [
    //     [AgeGroup::AgeGroup_0_to_4, AgeGroup::AgeGroup_5_to_14],
    //     [AgeGroup::AgeGroup_15_to_34, AgeGroup::AgeGroup_35_to_59],
    //     [AgeGroup::AgeGroup_60_to_79, AgeGroup::AgeGroup_80_plus],
    //     ...
    // ]
    // Each inner vector represents a household, and each AgeGroup represents a person in that household
    // We will now create persons for each household and assign them to the household locations.

    int household_index = 0;
    for (const auto& household : temp_households) {
        // Create a new household location
        auto household_id = household_locations[household_index];
        // Create persons for each age in the household
        for (const auto& age : household) {
            // Create a new person with the given age
            auto& person = world.add_person(household_id, age);
            // Assign the person to the household location
            person.set_assigned_location(household_id);
            person.set_assigned_location(hospitals[0]);
            person.set_assigned_location(icus[0]);
        }
        // Increment the household index
        household_index++;
    }

    // Quick check if we created the right amount household sizes:
    std::vector<int> household_sizes(5, 0);
    for (const auto& loc : world.get_locations()) {
        if (loc.get_type() != mio::abm::LocationType::Home) {
            continue; // Only count households
        }
        household_sizes[loc.get_persons().size() - 1]++;
    }

    for (size_t i = 0; i < household_sizes.size(); ++i) {
        if (household_sizes[i] != hh_per_size[i]) {
            std::cerr << "Error: Expected " << hh_per_size[i] << " households of size " << (i + 1) << ", but found "
                      << household_sizes[i] << ".\n";
        }
    }

    // Finally we assign schools, workplaces, shops, events to the persons based on their age groups and according to the rules:
    int count_elementary_schools = 0;
    int count_secondary_schools  = 0;
    int count_shops              = 0;
    int count_events             = 0;
    int count_workers            = 0;

    int number_of_4_to_14_years_old = age_vector[1];
    int prim_overhang               = 0;
    int sec_overhang                = 0;
    if (number_of_4_to_14_years_old >= infra.num_persons_elementary_schools) {
        if (number_of_4_to_14_years_old >= infra.num_persons_secondary_schools + infra.num_persons_elementary_schools) {
            prim_overhang = 0;
            sec_overhang  = 0;
        }
        else {
            prim_overhang = 0;
            sec_overhang  = infra.num_persons_secondary_schools + infra.num_persons_elementary_schools -
                           number_of_4_to_14_years_old;
        }
    }
    else {
        prim_overhang = infra.num_persons_elementary_schools - number_of_4_to_14_years_old;
        sec_overhang  = infra.num_persons_secondary_schools;
    }

    int amount_of_60plus_that_work = infra.num_worker - sec_overhang - age_vector[2] - age_vector[3];
    if (amount_of_60plus_that_work < 0) {
        amount_of_60plus_that_work = 0; // No one over 60 works
    }

    for (auto& person : world.get_persons()) {

        if (person.get_age() != age_group_0_to_4) {
            // Assign to a random shop if not in age group 0-4
            person.set_assigned_location(shops[count_shops % infra.num_stores]);
            count_shops++;
        }
        // Assign to a random event
        person.set_assigned_location(events[count_events % infra.num_events]);
        count_events++;

        if (person.get_age() == age_group_5_to_14) {
            // Assign to a school if in age group 0-14
            if (count_elementary_schools < infra.num_persons_elementary_schools) {
                person.set_assigned_location(prim_schools[count_elementary_schools % infra.num_elementary_schools]);
                count_elementary_schools++;
                continue;
            }
            else if (count_secondary_schools < infra.num_persons_secondary_schools) {
                person.set_assigned_location(sec_schools[count_secondary_schools % infra.num_secondary_schools]);
                count_secondary_schools++;
                continue;
            }
        }
        else if ((prim_overhang > 0 || sec_overhang > 0) && person.get_age() == age_group_15_to_34) {
            if (prim_overhang > 0) {
                // Assign to a primary school if there are overhangs
                person.set_assigned_location(prim_schools[count_elementary_schools % infra.num_elementary_schools]);
                count_elementary_schools++;
                prim_overhang--;
                continue;
            }
            else if (sec_overhang > 0) {
                // Assign to a secondary school if there are overhangs
                person.set_assigned_location(sec_schools[count_secondary_schools % infra.num_secondary_schools]);
                count_secondary_schools++;
                sec_overhang--;
                continue;
            }
        }

        else if (person.get_age() == age_group_35_to_59 || person.get_age() == age_group_15_to_34) {
            // Assign to a workplace if in age group 35-59
            if (count_workers < infra.num_worker) {
                person.set_assigned_location(workplaces[count_workers % infra.num_workplaces]);
                count_workers++;
                continue;
            }
        }
        else if (person.get_age() == age_group_60_to_79 || person.get_age() == age_group_80_plus) {
            // Assign to a workplace if in age group 60-79 or 80+
            if (amount_of_60plus_that_work > 0) {
                person.set_assigned_location(workplaces[count_workers % infra.num_workplaces]);
                count_workers++;
                amount_of_60plus_that_work--;
                continue;
            }
        }
    }

    int workerss         = 0;
    int school_attendees = 0;
    // Check if we did everything correctly
    for (auto& person : world.get_persons()) {
        // We count amount of workers and school attendees

        if (person.get_assigned_locations()[1] != std::numeric_limits<uint32_t>::max()) {
            school_attendees++;
        }
        if (person.get_assigned_locations()[2] != std::numeric_limits<uint32_t>::max()) {
            workerss++;
        }

        if (person.get_assigned_locations()[0] == std::numeric_limits<uint32_t>::max()) {
            //Home
            std::cerr << "Error: Person " << person.get_person_id() << " has no assigned home location.\n";
        }
        if (person.get_assigned_locations()[1] != std::numeric_limits<uint32_t>::max() &&
            person.get_assigned_locations()[2] != std::numeric_limits<uint32_t>::max()) {
            std::cerr << "Error: Person " << person.get_person_id() << " has an wokr and school place assigned.\n";
        }

        if (person.get_assigned_locations()[3] == std::numeric_limits<uint32_t>::max()) {
            //Event
            std::cerr << "Error: Person " << person.get_person_id() << " has no assigned event location.\n";
        }
        if (person.get_assigned_locations()[4] == std::numeric_limits<uint32_t>::max() &&
            person.get_age() != age_group_0_to_4) {
            //Shop
            std::cerr << "Error: Person " << person.get_person_id() << " has no assigned shop location.\n";
        }

        if (person.get_assigned_locations()[5] == std::numeric_limits<uint32_t>::max()) {
            //Hospital
            std::cerr << "Error: Person " << person.get_person_id() << " has no assigned hospital location.\n";
        }
        if (person.get_assigned_locations()[6] == std::numeric_limits<uint32_t>::max()) {
            //ICU
            std::cerr << "Error: Person " << person.get_person_id() << " has no assigned ICU location.\n";
        }
    }

    // Check if we assigned the right amount of workers and school attendees
    if (workerss != infra.num_worker) {
        std::cerr << "Error: Expected " << infra.num_worker << " workers, but found " << workerss << ".\n";
    }
    if (school_attendees != infra.num_persons_elementary_schools + infra.num_persons_secondary_schools) {
        std::cerr << "Error: Expected " << (infra.num_persons_elementary_schools + infra.num_persons_secondary_schools)
                  << " school attendees, but found " << school_attendees << ".\n";
    }
}

void CityBuilder::print_city_summary(const CityConfig& config)
{
    auto infra = config.infrastructure();

    std::cout << "\n" << std::string(50, '=') << "\n";
    std::cout << "     GERMAN CITY SIMULATION SUMMARY\n";
    std::cout << std::string(50, '=') << "\n\n";

    // Population Overview
    std::cout << "📊 POPULATION OVERVIEW\n";
    std::cout << "Total Population: " << config.total_population << "\n";
    std::cout << "Expected Age Distribution (based on German 2023 census):\n";
    std::cout << "  • 0-4 years:   " << std::fixed << std::setprecision(1)
              << (CityParameters::AGE_DISTRIBUTION[0] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[0]) << " people)\n";
    std::cout << "  • 5-14 years:  " << (CityParameters::AGE_DISTRIBUTION[1] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[1]) << " people)\n";
    std::cout << "  • 15-34 years: " << (CityParameters::AGE_DISTRIBUTION[2] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[2]) << " people)\n";
    std::cout << "  • 35-59 years: " << (CityParameters::AGE_DISTRIBUTION[3] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[3]) << " people)\n";
    std::cout << "  • 60-79 years: " << (CityParameters::AGE_DISTRIBUTION[4] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[4]) << " people)\n";
    std::cout << "  • 80+ years:   " << (CityParameters::AGE_DISTRIBUTION[5] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[5]) << " people)\n\n";

    // Housing
    std::cout << "🏠 HOUSING\n";
    int total_households = std::accumulate(infra.num_households_hh_size.begin(), infra.num_households_hh_size.end(), 0);
    std::cout << "Households: " << total_households << "\n";
    std::cout << "Average household size: " << std::setprecision(2)
              << static_cast<double>(config.total_population) / total_households << " people\n";
    std::cout << "Households by size:\n";
    for (size_t i = 0; i < infra.num_households_hh_size.size(); ++i) {
        std::cout << "  • " << (i + 1) << " person: " << infra.num_households_hh_size[i] << "\n";
    }

    // Employment
    std::cout << "💼 EMPLOYMENT\n";
    std::cout << "Workplaces: " << infra.num_workplaces << "\n";
    std::cout << "Average employees per workplace: " << std::setprecision(1)
              << static_cast<double>(config.total_population) * CityParameters::InfrastructureRatios::EMPLOYMENT_RATE /
                     infra.num_workplaces
              << "\n";
    std::cout << "Expected employment rate: " << std::setprecision(1)
              << (CityParameters::InfrastructureRatios::EMPLOYMENT_RATE * 100) << "%\n\n";

    // Education
    std::cout << "🎓 EDUCATION\n";
    std::cout << "Elementary schools: " << infra.num_elementary_schools << "\n";
    std::cout << "Secondary schools: " << infra.num_secondary_schools << "\n";
    std::cout << "Total schools: " << (infra.num_elementary_schools + infra.num_secondary_schools) << "\n";
    std::cout << "Students per elementary school: ~" << std::setprecision(0) << infra.num_persons_elementary_schools
              << "\n";
    std::cout << "Students per secondary school: ~" << std::setprecision(0) << infra.num_persons_secondary_schools
              << "\n\n";

    // Healthcare
    std::cout << "🏥 HEALTHCARE\n";
    std::cout << "Hospitals: " << infra.num_hospitals << "\n";
    std::cout << "ICU units: " << infra.num_icus << "\n";
    std::cout << "People per hospital: " << std::setprecision(0)
              << static_cast<double>(config.total_population) / infra.num_hospitals << "\n";
    std::cout << "People per ICU: " << std::setprecision(0)
              << static_cast<double>(config.total_population) / infra.num_icus << "\n\n";

    // Retail & Services
    std::cout << "🛒 RETAIL & SERVICES\n";
    std::cout << "Stores: " << infra.num_stores << "\n";
    std::cout << "People per store: " << std::setprecision(0)
              << static_cast<double>(config.total_population) / infra.num_stores << "\n\n";

    // Social & Entertainment
    std::cout << "🎉 SOCIAL & ENTERTAINMENT\n";
    std::cout << "Event venues: " << infra.num_events << "\n";
    std::cout << "People per event venue: " << std::setprecision(0)
              << static_cast<double>(config.total_population) / infra.num_events << "\n\n";

    // Infrastructure Summary
    std::cout << "📈 INFRASTRUCTURE DENSITY\n";
    std::cout << "Total locations: "
              << (total_households + infra.num_workplaces + infra.num_elementary_schools + infra.num_secondary_schools +
                  infra.num_hospitals + infra.num_icus + infra.num_stores + infra.num_events)
              << "\n";
    std::cout << "Locations per 1000 people: " << std::setprecision(1)
              << static_cast<double>(total_households + infra.num_workplaces + infra.num_elementary_schools +
                                     infra.num_secondary_schools + infra.num_hospitals + infra.num_icus +
                                     infra.num_stores + infra.num_events) *
                     1000.0 / config.total_population
              << "\n\n";

    std::cout << "📍 DATA SOURCES\n";
    std::cout << "• German Federal Statistical Office (Destatis) 2023\n";
    std::cout << "• German Trade Association (HDE) 2023\n";

    std::cout << "\n" << std::string(50, '=') << "\n\n";
}

void CityBuilder::save_city_to_file(const CityConfig& config, const std::string& filename)
{
    auto infra = config.infrastructure();
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Error: Could not open file " << filename << " for writing\n";
        return;
    }

    // Write header
    ofs << "City Simulation Statistics\n";
    ofs << "=========================\n";

    // Population Statistics
    ofs << "POPULATION STATISTICS\n";
    ofs << "--------------------\n";
    ofs << "total_population," << config.total_population << "\n";

    // Age group distribution
    ofs << "age_group_0_to_4," << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[0])
        << "\n";
    ofs << "age_group_5_to_14," << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[1])
        << "\n";
    ofs << "age_group_15_to_34," << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[2])
        << "\n";
    ofs << "age_group_35_to_59," << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[3])
        << "\n";
    ofs << "age_group_60_to_79," << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[4])
        << "\n";
    ofs << "age_group_80_plus," << static_cast<int>(config.total_population * CityParameters::AGE_DISTRIBUTION[5])
        << "\n";

    // Age percentages
    ofs << "age_percentage_0_to_4," << std::fixed << std::setprecision(2) << (CityParameters::AGE_DISTRIBUTION[0] * 100)
        << "\n";
    ofs << "age_percentage_5_to_14," << (CityParameters::AGE_DISTRIBUTION[1] * 100) << "\n";
    ofs << "age_percentage_15_to_34," << (CityParameters::AGE_DISTRIBUTION[2] * 100) << "\n";
    ofs << "age_percentage_35_to_59," << (CityParameters::AGE_DISTRIBUTION[3] * 100) << "\n";
    ofs << "age_percentage_60_to_79," << (CityParameters::AGE_DISTRIBUTION[4] * 100) << "\n";
    ofs << "age_percentage_80_plus," << (CityParameters::AGE_DISTRIBUTION[5] * 100) << "\n\n";

    // Household Statistics
    ofs << "HOUSEHOLD STATISTICS\n";
    ofs << "-------------------\n";
    int total_households = std::accumulate(infra.num_households_hh_size.begin(), infra.num_households_hh_size.end(), 0);
    ofs << "total_households," << total_households << "\n";
    ofs << "average_household_size," << std::setprecision(2)
        << static_cast<double>(config.total_population) / total_households << "\n";

    // Household size distribution
    for (size_t i = 0; i < infra.num_households_hh_size.size(); ++i) {
        ofs << "households_size_" << (i + 1) << "," << infra.num_households_hh_size[i] << "\n";
        ofs << "households_size_" << (i + 1) << "_percentage," << std::setprecision(2)
            << (static_cast<double>(infra.num_households_hh_size[i]) / total_households * 100) << "\n";
    }
    ofs << "\n";

    // Employment Statistics
    ofs << "EMPLOYMENT STATISTICS\n";
    ofs << "--------------------\n";
    ofs << "num_workplaces," << infra.num_workplaces << "\n";
    ofs << "num_workers," << infra.num_worker << "\n";
    ofs << "employment_rate," << std::setprecision(2) << (CityParameters::InfrastructureRatios::EMPLOYMENT_RATE * 100)
        << "\n";
    ofs << "average_employees_per_workplace," << std::setprecision(1)
        << static_cast<double>(infra.num_worker) / infra.num_workplaces << "\n";

    // Calculate workers by age group
    auto age_vector     = create_age_vector(config.total_population);
    int workers_15_34   = std::min(age_vector[2], infra.num_worker);
    int workers_35_59   = std::min(age_vector[3], std::max(0, infra.num_worker - workers_15_34));
    int workers_60_plus = std::max(0, infra.num_worker - workers_15_34 - workers_35_59);

    ofs << "workers_age_15_34," << workers_15_34 << "\n";
    ofs << "workers_age_35_59," << workers_35_59 << "\n";
    ofs << "workers_age_60_plus," << workers_60_plus << "\n\n";

    // Education Statistics
    ofs << "EDUCATION STATISTICS\n";
    ofs << "-------------------\n";
    ofs << "num_elementary_schools," << infra.num_elementary_schools << "\n";
    ofs << "num_secondary_schools," << infra.num_secondary_schools << "\n";
    ofs << "total_schools," << (infra.num_elementary_schools + infra.num_secondary_schools) << "\n";
    ofs << "elementary_school_students," << infra.num_persons_elementary_schools << "\n";
    ofs << "secondary_school_students," << infra.num_persons_secondary_schools << "\n";
    ofs << "total_students," << (infra.num_persons_elementary_schools + infra.num_persons_secondary_schools) << "\n";
    ofs << "students_per_elementary_school," << std::setprecision(1)
        << static_cast<double>(infra.num_persons_elementary_schools) / infra.num_elementary_schools << "\n";
    ofs << "students_per_secondary_school," << std::setprecision(1)
        << static_cast<double>(infra.num_persons_secondary_schools) / infra.num_secondary_schools << "\n";
    ofs << "school_age_population_5_14," << age_vector[1] << "\n\n";

    // Healthcare Statistics
    ofs << "HEALTHCARE STATISTICS\n";
    ofs << "--------------------\n";
    ofs << "num_hospitals," << infra.num_hospitals << "\n";
    ofs << "num_icus," << infra.num_icus << "\n";
    ofs << "people_per_hospital," << std::setprecision(0)
        << static_cast<double>(config.total_population) / infra.num_hospitals << "\n";
    ofs << "people_per_icu," << std::setprecision(0) << static_cast<double>(config.total_population) / infra.num_icus
        << "\n\n";

    // Retail & Services Statistics
    ofs << "RETAIL & SERVICES STATISTICS\n";
    ofs << "---------------------------\n";
    ofs << "num_stores," << infra.num_stores << "\n";
    ofs << "people_per_store," << std::setprecision(0)
        << static_cast<double>(config.total_population) / infra.num_stores << "\n";
    ofs << "shoppers," << (config.total_population - age_vector[0]) << "\n"; // All except 0-4 age group
    ofs << "shoppers_per_store," << std::setprecision(1)
        << static_cast<double>(config.total_population - age_vector[0]) / infra.num_stores << "\n\n";

    // Social & Entertainment Statistics
    ofs << "SOCIAL & ENTERTAINMENT STATISTICS\n";
    ofs << "--------------------------------\n";
    ofs << "num_events," << infra.num_events << "\n";
    ofs << "people_per_event_venue," << std::setprecision(0)
        << static_cast<double>(config.total_population) / infra.num_events << "\n\n";

    // Infrastructure Density
    ofs << "INFRASTRUCTURE DENSITY\n";
    ofs << "---------------------\n";
    int total_locations = total_households + infra.num_workplaces + infra.num_elementary_schools +
                          infra.num_secondary_schools + infra.num_hospitals + infra.num_icus + infra.num_stores +
                          infra.num_events;
    ofs << "total_locations," << total_locations << "\n";
    ofs << "locations_per_1000_people," << std::setprecision(1)
        << static_cast<double>(total_locations) * 1000.0 / config.total_population << "\n";
    ofs << "households_per_1000_people," << std::setprecision(1)
        << static_cast<double>(total_households) * 1000.0 / config.total_population << "\n";
    ofs << "workplaces_per_1000_people," << std::setprecision(1)
        << static_cast<double>(infra.num_workplaces) * 1000.0 / config.total_population << "\n";
    ofs << "schools_per_1000_people," << std::setprecision(1)
        << static_cast<double>(infra.num_elementary_schools + infra.num_secondary_schools) * 1000.0 /
               config.total_population
        << "\n\n";

    // Contact Network Potential
    ofs << "CONTACT NETWORK POTENTIAL\n";
    ofs << "------------------------\n";
    ofs << "max_workplace_contacts," << infra.num_worker << "\n";
    ofs << "max_school_contacts," << (infra.num_persons_elementary_schools + infra.num_persons_secondary_schools)
        << "\n";
    ofs << "max_shopping_contacts," << (config.total_population - age_vector[0]) << "\n";
    ofs << "max_social_event_contacts," << config.total_population << "\n";
    ofs << "non_working_non_school_population,"
        << (config.total_population - infra.num_worker - infra.num_persons_elementary_schools -
            infra.num_persons_secondary_schools)
        << "\n\n";

    // German Demographic Ratios Used
    ofs << "GERMAN DEMOGRAPHIC RATIOS USED\n";
    ofs << "------------------------------\n";
    ofs << "employment_rate_used," << std::setprecision(3) << CityParameters::InfrastructureRatios::EMPLOYMENT_RATE
        << "\n";

    // Infrastructure ratios (you may need to expose these from CityParameters)
    ofs << "households_per_person," << std::setprecision(4)
        << static_cast<double>(total_households) / config.total_population << "\n";
    ofs << "workplaces_per_person," << std::setprecision(4)
        << static_cast<double>(infra.num_workplaces) / config.total_population << "\n";
    ofs << "stores_per_person," << std::setprecision(4)
        << static_cast<double>(infra.num_stores) / config.total_population << "\n";
    ofs << "events_per_person," << std::setprecision(4)
        << static_cast<double>(infra.num_events) / config.total_population << "\n\n";

    // Summary Statistics for Quick Reference
    ofs << "SUMMARY STATISTICS\n";
    ofs << "-----------------\n";
    ofs << "population_density_category,";
    if (config.total_population <= 500)
        ofs << "small_town\n";
    else if (config.total_population <= 2000)
        ofs << "medium_town\n";
    else
        ofs << "large_town\n";

    ofs << "average_persons_per_household," << std::setprecision(2)
        << static_cast<double>(config.total_population) / total_households << "\n";
    ofs << "worker_to_population_ratio," << std::setprecision(3)
        << static_cast<double>(infra.num_worker) / config.total_population << "\n";
    ofs << "student_to_population_ratio," << std::setprecision(3)
        << static_cast<double>(infra.num_persons_elementary_schools + infra.num_persons_secondary_schools) /
               config.total_population
        << "\n";

    ofs.close();
    std::cout << "City statistics saved to " << filename << "\n";
}