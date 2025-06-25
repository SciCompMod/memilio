#include "../include/world_creator.h"
#include "../include/constants.h"
#include "../include/parameter_setter.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "abm/abm.h"

std::map<uint32_t, bool> read_infection_data(const std::string& filename)
{
    std::map<uint32_t, bool> infected_status;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return infected_status;
    }

    std::string line;
    std::getline(file, line); // Skip header

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        uint32_t id;
        std::string table;
        int infected;

        if (iss >> id >> table >> infected) {
            infected_status[id] = (infected == 1);
        }
    }

    file.close();
    return infected_status;
}

mio::abm::World create_world_from_file(const std::string& infection_data_file, int number_of_persons)
{
    // Create world with 6 age groups
    auto world = mio::abm::World(num_age_groups);

    // Set infection parameters
    set_parameters(world.parameters);

    // Read infection data
    auto infected_status = read_infection_data(infection_data_file);

    // Create a hospital and ICU for severe cases
    auto hospital = world.add_location(mio::abm::LocationType::Hospital);
    auto icu      = world.add_location(mio::abm::LocationType::ICU);

    // Create common locations
    auto event  = world.add_location(mio::abm::LocationType::SocialEvent);
    auto shop   = world.add_location(mio::abm::LocationType::BasicsShop);
    auto school = world.add_location(mio::abm::LocationType::School);
    auto work   = world.add_location(mio::abm::LocationType::Work);

    // Create homes (one per table) for the initial persons from the file
    std::map<std::string, mio::abm::LocationId> table_homes;
    for (const auto& [id, infected] : infected_status) {
        // Create home location for each person
        auto home                       = world.add_location(mio::abm::LocationType::Home);
        table_homes[std::to_string(id)] = home;
    }

    // Create people from the infection data and assign them to locations
    std::vector<double> weights_age      = {0.05, 0.25, 0.2, 0.4, 0.07, 0.03};
    std::vector<double> weight_home_size = {0.2, 0.3, 0.2, 0.2, 0.1};
    int initial_person_count             = 0;

    for (const auto& [id, infected] : infected_status) {
        mio::AgeGroup age =
            mio::AgeGroup(mio::DiscreteDistribution<size_t>::get_instance()(world.get_rng(), weights_age));

        // Create person and assign them to their home
        auto& person = world.add_person(table_homes[std::to_string(id)], age);
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        person.set_assigned_location(event);
        person.set_assigned_location(world.get_individualized_location(table_homes[std::to_string(id)]));
        person.set_assigned_location(shop);

        // Assign school or work based on age
        if (age == age_group_5_to_14) {
            person.set_assigned_location(school);
        }
        if (age == age_group_15_to_34 || age == age_group_35_to_59) {
            person.set_assigned_location(work);
        }

        // If infected according to data, set them as exposed
        if (infected) {
            auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world.parameters, mio::abm::TimePoint(0),
                                                         mio::abm::InfectionState::Exposed));
        }

        initial_person_count++;
    }

    // If a higher number of persons is requested, create additional people
    if (number_of_persons > initial_person_count) {
        // Create additional homes with 1-5 persons per home
        std::vector<mio::abm::LocationId> additional_homes;
        std::vector<int> persons_per_home;

        // First create homes and determine how many people in each
        int remaining_persons = number_of_persons - initial_person_count;
        while (remaining_persons > 0) {
            // Randomly decide how many people in this home (1-5)
            int persons_in_home =
                std::min(remaining_persons, 1 + static_cast<int>(mio::DiscreteDistribution<size_t>::get_instance()(
                                                    world.get_rng(), weight_home_size)));

            // Create a new home
            auto home = world.add_location(mio::abm::LocationType::Home);
            additional_homes.push_back(home);
            persons_per_home.push_back(persons_in_home);

            remaining_persons -= persons_in_home;
        }

        // Now create the people and assign them to homes
        for (size_t i = 0; i < additional_homes.size(); ++i) {
            for (int j = 0; j < persons_per_home[i]; ++j) {
                // Randomly determine age based on weights
                mio::AgeGroup age =
                    mio::AgeGroup(mio::DiscreteDistribution<size_t>::get_instance()(world.get_rng(), weights_age));

                // Create person and assign them to their home
                auto& person = world.add_person(additional_homes[i], age);

                // Assign to common locations
                person.set_assigned_location(world.get_individualized_location(additional_homes[i]));
                person.set_assigned_location(hospital);
                person.set_assigned_location(icu);
                person.set_assigned_location(event);
                person.set_assigned_location(shop);

                // Assign school or work based on age
                if (age == age_group_5_to_14) {
                    person.set_assigned_location(school);
                }
                if (age == age_group_15_to_34 || age == age_group_35_to_59) {
                    person.set_assigned_location(work);
                }
            }
        }
    }

    set_local_parameters(world);
    std::cout << "Created world with " << world.get_persons().size() << " people" << std::endl;
    std::cout << "Number of homes: "
              << std::count_if(world.get_locations().begin(), world.get_locations().end(),
                               [](const mio::abm::Location& loc) {
                                   return loc.get_type() == mio::abm::LocationType::Home;
                               })
              << std::endl;

    return world;
}