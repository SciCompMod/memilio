#include "../include/event_simulator.h"
#include "../include/constants.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <parameter_setter.h>

mio::IOResult<std::string> EventSimulator::get_panvadere_file_for_event_type(EventType type)
{
    auto it = EVENT_TYPE_TO_FILE.find(type);
    if (it != EVENT_TYPE_TO_FILE.end()) {
        return mio::success(it->second);
    }

    // No fallback return error handling
    std::cerr << "Error: No Panvadere file defined for event type " << event_type_to_string(type) << std::endl;
    std::cerr << "Please define a Panvadere file for this event type in EVENT_TYPE_TO_FILE." << std::endl;

    return mio::IOStatus(std::make_error_code(std::errc::invalid_argument), "No Panvadere file defined for event type");
}

std::string EventSimulationConfig::get_panvadere_file() const
{
    auto result = EventSimulator::get_panvadere_file_for_event_type(type);
    if (result.has_error()) {
        // Fallback to default if there's an error
        std::cerr << "Error getting Panvadere file for event type " << EventSimulator::event_type_to_string(type)
                  << ": " << result.error().message() << std::endl;
        std::cerr << "Using default infections file: ./data/default/infections.txt" << std::endl;
        // Return default file path
        return "./data/default/infections.txt";
    }
    return result.value();
}

mio::IOResult<double> EventSimulator::calculate_infection_parameter_k(const EventSimulationConfig& config,
                                                                      mio::abm::World& city,
                                                                      std::map<uint32_t, uint32_t>& event_map)
{

    // Check how many infectious persons there are in the panvadere file
    BOOST_OUTCOME_TRY(auto panvadere_file, get_panvadere_file_for_event_type(config.type));

    // Read Panvadere infection data
    BOOST_OUTCOME_TRY(auto panvadere_data, read_panv_file_restaurant(panvadere_file));

    if (panvadere_data.empty()) {
        return mio::failure(mio::StatusCode::InvalidFileFormat, "No valid infection data found in Panvadere file");
    }
    // Count infected persons in the Panvadere data
    size_t infected_count = 0;
    for (const auto& entry : panvadere_data) {
        if (std::get<2>(entry)) { // Check if infected status is true
            infected_count++;
        }
    }

    // Now we need to simulate the event world and search for a suitable K parameter with which the simulation has the same number of infected persons as in the Panvadere data

    // Rudimentary grid search for K parameter
    double k_min = 1.0; // Minimum K value
    double k_max = 20.0; // Maximum K value

    double k_step                  = 1.0; // Step size for K value
    double best_k                  = k_min;
    size_t best_infected_count     = 0;
    const int num_runs             = 200; // Number of runs for averaging
    double best_avg_infected_count = -1.0;

    for (double k = k_min; k <= k_max; k += k_step) {
        size_t total_infected_for_k = 0;
        for (int i = 0; i < num_runs; ++i) {
            BOOST_OUTCOME_TRY(auto calculation_world, create_event_world(config, city, event_map));
            set_parameters(calculation_world.parameters);

            double ratio = 1.0;
            if (config.type == EventType::Restaurant_Table_Equals_Half_Household ||
                config.type == EventType::Restaurant_Table_Equals_Household) {
                // Set restaurant-specific parameters
                ratio = 3.0 / config.event_duration_hours; // 3 hours divided by event duration
            }
            else if (config.type == EventType::WorkMeeting_Many_Meetings ||
                     config.type == EventType::WorkMeeting_Few_Meetings) {
                ratio = 8.0 / config.event_duration_hours; // 3 hours divided by event duration
            }
            set_local_parameters_event(calculation_world, ratio);

            calculation_world.use_migration_rules(false); // Disable migration rules for this simulation
            auto t0 = mio::abm::TimePoint(0); // Start time per simulation
            auto tmax =
                mio::abm::TimePoint(0) + mio::abm::hours(config.event_duration_hours); // End time per simulation

            // Set infection parameter K
            calculation_world.parameters.get<mio::abm::InfectionRateFromViralShed>() = k;

            auto sim = mio::abm::Simulation(t0, std::move(calculation_world));
            sim.advance(tmax);

            // Count infected persons in the simulation
            size_t current_infected_count = 0;
            for (const auto& person : sim.get_world().get_persons()) {
                if (person.is_infected(tmax)) {
                    current_infected_count++;
                }
            }
            total_infected_for_k += current_infected_count;
        }

        double avg_infected_count = static_cast<double>(total_infected_for_k) / num_runs;

        // Check if this K value gives a better match to the Panvadere data
        if (best_avg_infected_count < 0 ||
            std::abs(avg_infected_count - static_cast<double>(infected_count)) <
                std::abs(best_avg_infected_count - static_cast<double>(infected_count))) {
            best_k                  = k;
            best_avg_infected_count = avg_infected_count;
            best_infected_count     = static_cast<size_t>(std::round(avg_infected_count));
            mio::unused(best_infected_count); // Avoid unused variable warning
        }
    }

    return mio::success(best_k);
}

mio::IOResult<std::vector<std::tuple<uint32_t, std::string, bool>>>
EventSimulator::read_panv_file_restaurant(const std::string& filename)
{
    std::vector<std::tuple<uint32_t, std::string, bool>> infected_status;

    // Check if filename is empty
    if (filename.empty()) {
        return mio::failure(mio::StatusCode::InvalidValue, "Filename cannot be empty");
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return mio::failure(mio::StatusCode::InvalidFileFormat, "Failed to open file: " + filename);
    }

    std::string line;
    size_t line_number = 0;

    // Skip header line (pedestrianId	table	infected)
    if (std::getline(file, line)) {
        line_number++;
        // Verify it's actually a header by checking for expected column names
        if (line.find("pedestrianId") == std::string::npos || line.find("table") == std::string::npos ||
            line.find("infected") == std::string::npos) {
            std::cerr << "Warning: Expected header line not found. Processing first line as data." << std::endl;
            // Reset file to beginning if header verification fails
            file.clear();
            file.seekg(0);
            line_number = 0;
        }
    }

    while (std::getline(file, line)) {
        line_number++;

        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos) {
            continue;
        }

        // Parse tab-delimited data
        std::istringstream iss(line);
        std::string pedestrian_id_str, table_str, infected_str;

        // Extract fields using tab delimiter
        if (std::getline(iss, pedestrian_id_str, '\t') && std::getline(iss, table_str, '\t') &&
            std::getline(iss, infected_str)) {

            try {
                // Parse pedestrian ID
                uint32_t pedestrian_id = static_cast<uint32_t>(std::stoul(pedestrian_id_str));

                // Parse table (remove any trailing whitespace)
                table_str.erase(table_str.find_last_not_of(" \t\r\n") + 1);

                // Parse infection status
                infected_str.erase(infected_str.find_last_not_of(" \t\r\n") + 1);
                int infected_int = std::stoi(infected_str);

                if (infected_int < 0 || infected_int > 1) {
                    std::cerr << "Warning: Invalid infection status " << infected_int << " for pedestrian "
                              << pedestrian_id << " at line " << line_number << ", treating as not infected"
                              << std::endl;
                    infected_status.emplace_back(pedestrian_id, table_str, false);
                }
                else {
                    bool is_infected = (infected_int == 1);
                    infected_status.emplace_back(pedestrian_id, table_str, is_infected);
                }
            }
            catch (const std::exception& e) {
                std::cerr << "Warning: Failed to parse line " << line_number << ": '" << line << "' - " << e.what()
                          << std::endl;
                continue;
            }
        }
        else {
            std::cerr << "Warning: Failed to parse line " << line_number << " (expected 3 tab-separated fields): '"
                      << line << "'" << std::endl;
        }
    }

    file.close();

    if (infected_status.empty()) {
        std::cerr << "Warning: No valid infection data found in file: " << filename << std::endl;
        return mio::failure(mio::StatusCode::InvalidFileFormat, "No valid data found in file");
    }

    // std::cout << "Successfully read " << infected_status.size() << " persons from " << filename << std::endl;

    return mio::success(infected_status);
}

mio::IOResult<std::vector<uint32_t>>
EventSimulator::initialize_from_event_simulation(mio::RandomNumberGenerator rng, EventType event_type,
                                                 std::map<uint32_t, uint32_t>& event_map)
{
    // Map to households based on event type
    std::vector<uint32_t> household_infections;

    switch (event_type) {
    case EventType::Restaurant_Table_Equals_Half_Household:
    case EventType::Restaurant_Table_Equals_Random:
    case EventType::Restaurant_Table_Equals_Household: {
        // BOOST_OUTCOME_TRY(auto panvadere_file, get_panvadere_file_for_event_type(event_type));

        // // Read Panvadere infection data
        // BOOST_OUTCOME_TRY(auto panvadere_data, read_panv_file_restaurant(panvadere_file));

        // size_t amount_of_infected = 0;
        // for (const auto& entry : panvadere_data) {
        //     bool is_infected = std::get<2>(entry);
        //     if (is_infected) {
        //         amount_of_infected++;
        //     }
        // }

        // auto event_map_copy = event_map;
        // //shuffle the event_map_copy to randomize the household mapping
        // std::vector<std::pair<uint32_t, uint32_t>> event_map_vector(event_map_copy.begin(), event_map_copy.end());
        // std::shuffle(event_map_vector.begin(), event_map_vector.end(), rng);

        // // take amount_of_infected households from the shuffled event_map_vector
        // for (size_t i = 0; i < amount_of_infected && i < event_map_vector.size(); ++i) {
        //     household_infections.push_back(event_map_vector[i].second);
        // }

        std::vector<uint32_t> return_vector;
        mio::unused(rng); // Avoid unused variable warning
        return_vector.reserve(13); // Reserve space for 13 households -> this is one of the worst case scenarios
        return_vector = {0, 4, 9, 13, 18, 20, 25, 1, 5, 10, 14, 21, 26, 27};
        for (size_t i = 0; i < return_vector.size(); ++i) {
            household_infections.push_back(event_map[return_vector[i]]); // Map to household IDs
        }

        break;
    }
    case EventType::WorkMeeting_Many_Meetings:
        break;
    case EventType::WorkMeeting_Few_Meetings:
        break;
    default:
        // Map events to households
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Event type not supported for Panvadere data mapping: " + event_type_to_string(event_type));
        break;
    }

    return mio::success(household_infections);
}

mio::IOResult<std::vector<uint32_t>> EventSimulator::initialize_from_panvadere(EventType event_type,
                                                                               std::map<uint32_t, uint32_t>& event_map)
{

    // Map to households based on event type
    std::vector<uint32_t> household_infections;

    // Get the appropriate file for this event type

    switch (event_type) {
    case EventType::Restaurant_Table_Equals_Household:
    case EventType::Restaurant_Table_Equals_Random:
    case EventType::Restaurant_Table_Equals_Half_Household: {
        BOOST_OUTCOME_TRY(auto panvadere_file, get_panvadere_file_for_event_type(event_type));

        // Read Panvadere infection data
        BOOST_OUTCOME_TRY(auto panvadere_data, read_panv_file_restaurant(panvadere_file));

        for (const auto& entry : panvadere_data) {
            uint32_t person_id = std::get<0>(entry);
            std::string table  = std::get<1>(entry);
            bool is_infected   = std::get<2>(entry);

            // Map person ID to household ID
            if (is_infected) {
                // Check if the table is already mapped to a household
                auto it = event_map.find(person_id);
                if (it != event_map.end()) {
                    // If already mapped, add to household infections
                    household_infections.push_back(it->second);
                }
                else {
                    return mio::failure(mio::StatusCode::InvalidValue, "Person ID " + std::to_string(person_id) +
                                                                           " not found in event map for table " +
                                                                           table);
                }
            }
        }
        break;
    }
    case EventType::WorkMeeting_Many_Meetings:
        break;
    case EventType::WorkMeeting_Few_Meetings:
        break;
    default:
        // Map events to households
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Event type not supported for Panvadere data mapping: " + event_type_to_string(event_type));
        break;
    }

    return mio::success(household_infections);
}

std::string EventSimulator::event_type_to_string(EventType type)
{
    switch (type) {
    case EventType::Restaurant_Table_Equals_Household:
        return "RestaurantEqualHH";
    case EventType::Restaurant_Table_Equals_Half_Household:
        return "RestaurantNotEqualHH";
    case EventType::WorkMeeting_Many_Meetings:
        return "WorkMeetingMany";
    case EventType::WorkMeeting_Few_Meetings:
        return "WorkMeetingFew";
    default:
        return "Unknown";
    }
}

std::string EventSimulator::simulation_type_to_string(SimType type)
{
    switch (type) {
    case SimType::Panvadere:
        return "Panvadere";
    case SimType::Memilio:
        return "Memilio";
    case SimType::Both:
        return "Both";
    default:
        return "Unknown";
    }
}

mio::IOResult<std::map<uint32_t, uint32_t>> EventSimulator::map_restaurant_tables_to_households(mio::abm::World& city)
{
    // This funciton maps a simulation id to a household id in the world simulation.

    // First we read in the infection data from the Panvadere simulation
    BOOST_OUTCOME_TRY(auto panvadere_file,
                      get_panvadere_file_for_event_type(EventType::Restaurant_Table_Equals_Household));
    BOOST_OUTCOME_TRY(auto panvadere_data, read_panv_file_restaurant(panvadere_file));

    // For this we want to map each table to a full household. For this we look at a table and assign the first household with that
    // amount of people to the table. If the table has more than 5 people we devide it into two households: one household for the first 5 people and one for the rest.

    //First we calculate how many person are at each table and save the ids for each table
    std::vector<std::pair<std::string, uint32_t>> tables;
    std::vector<std::tuple<std::string, std::vector<uint32_t>>> table_person_ids; // Store person IDs for each table
    for (const auto& entry : panvadere_data) {
        const std::string& table = std::get<1>(entry);
        // Find the table in the vector or create a new entry
        auto it = std::find_if(tables.begin(), tables.end(), [&table](const std::pair<std::string, uint32_t>& entry) {
            return entry.first == table;
        });
        if (it != tables.end()) {
            it->second++; // Increment count for existing table
            // Also add the person ID to the table_person_ids vector
            auto person_id = std::get<0>(entry);
            auto person_it = std::find_if(table_person_ids.begin(), table_person_ids.end(),
                                          [&table](const std::tuple<std::string, std::vector<uint32_t>>& entry) {
                                              return std::get<0>(entry) == table;
                                          });
            if (person_it != table_person_ids.end()) {
                std::get<1>(*person_it).push_back(person_id); // Add person ID
            }
            else {
                // If the table is not found, create a new entry with the person ID, this should not happen
                std::cerr << "Warning: Table " << table << " not found in table_person_ids, creating new entry."
                          << std::endl;
                table_person_ids.emplace_back(table, std::vector<uint32_t>{person_id}); // Add new entry with person ID
            }
        }
        else {
            tables.emplace_back(table, 1); // Add new table with count 1
            // Also add the person ID to the table_person_ids vector
            table_person_ids.emplace_back(table, std::vector<uint32_t>{std::get<0>(entry)});
        }
    }

    // Secondly we now assign each table a vector of household sizes
    std::map<std::string, std::vector<uint32_t>> table_household_sizes;
    for (const auto& [table, count] : tables) {
        if (count <= 5) {
            // If the table has 5 or less people, we assign it to one household
            table_household_sizes[table].push_back(count);
        }
        else if (count > 5 && count <= 10) {
            // If the table has more than 5 people, we divide it into two households
            uint32_t household_size1 = 5; // First household takes 5 people
            uint32_t household_size2 = count - 5; // Second household takes the rest
            table_household_sizes[table].push_back(household_size1);
            table_household_sizes[table].push_back(household_size2);
        }
        else {
            // If the table has more than 10 people, we divide it into multiple households of 5 people each
            uint32_t num_households = count / 5;
            for (uint32_t i = 0; i < num_households; ++i) {
                table_household_sizes[table].push_back(5);
            }
            if (count % 5 != 0) {
                // Add remaining people as a smaller household
                table_household_sizes[table].push_back(count % 5);
            }
        }
    }

    // Now we search for the households in the city world and assign them to the tables
    std::map<uint32_t, uint32_t> household_map; // Maps panv id to household id in the world simulation
    std::vector<uint32_t> household_ids_reserved; // Store household IDs which have been assigned to tables
    for (const auto& table_entry : table_household_sizes) {
        const std::string& table                     = table_entry.first;
        const std::vector<uint32_t>& household_sizes = table_entry.second;
        int table_global_counter                     = 0;
        for (const auto& household_size : household_sizes) {
            // Find a household in the city world with the same size which hasnt been assigned yet
            auto household_it =
                std::find_if(city.get_locations().begin(), city.get_locations().end(),
                             [&household_size, &household_ids_reserved](const auto& location) {
                                 return location.get_persons().size() == household_size &&
                                        std::find(household_ids_reserved.begin(), household_ids_reserved.end(),
                                                  location.get_index()) == household_ids_reserved.end();
                             });
            // add to reserved household ids
            household_ids_reserved.push_back(household_it->get_index());
            if (household_it != city.get_locations().end()) {
                // If we found a household, map the table to the household id
                for (size_t i = 0; i < household_size; i++) {
                    // Get the person ID from the table_person_ids vector
                    auto person_it =
                        std::find_if(table_person_ids.begin(), table_person_ids.end(), [&table](const auto& entry) {
                            return std::get<0>(entry) == table;
                        });
                    if (person_it != table_person_ids.end() && i < std::get<1>(*person_it).size()) {
                        uint32_t person_id = std::get<1>(*person_it)[table_global_counter++];
                        household_map[person_id] =
                            household_it->get_persons().at(i).get()->get_person_id(); // Map person ID to household ID
                    }
                    else {
                        std::cerr << "Warning: Not enough persons for table " << table << " with size "
                                  << household_size << ". Skipping this person." << std::endl;
                    }
                }
            }
            else {
                // If we didn't find a household, we can either skip this table or create a new household
                std::cerr << "Warning: No household found for table " << table << " with size " << household_size
                          << ". Skipping this table." << std::endl;
            }
        }
    }

    // Summarize the household map
    // std::cout << "Household map created with " << household_map.size() << " entries." << std::endl;
    // for (const auto& [person_id_pre, person_id_sim] : household_map) {
    //     std::cout << "Person ID Pre: " << person_id_pre << " Person ID Sim: " << person_id_sim << " HouseholdID: "
    //               << city.get_person(person_id_sim).get_assigned_location_index(mio::abm::LocationType::Home)
    //               << std::endl;
    // }

    return mio::success(household_map);
}

mio::IOResult<std::map<uint32_t, uint32_t>>
EventSimulator::map_two_person_per_household_restaurant_tables_to_households(mio::abm::World& city)
{
    // This funciton maps a simulation id to a household id in the world simulation.

    // First we read in the infection data from the Panvadere simulation
    BOOST_OUTCOME_TRY(auto panvadere_file,
                      get_panvadere_file_for_event_type(EventType::Restaurant_Table_Equals_Household));
    BOOST_OUTCOME_TRY(auto panvadere_data, read_panv_file_restaurant(panvadere_file));

    // For this we want to map each two seats to a new household. For this we look at a table and assign the first household with that
    // amount of people to the table. If the table has more than 5 people we devide it into two households: one household for the first 5 people and one for the rest.

    //First we calculate how many person are at each table and save the ids for each table
    std::vector<std::pair<std::string, uint32_t>> tables;
    std::vector<std::tuple<std::string, std::vector<uint32_t>>> table_person_ids; // Store person IDs for each table
    for (const auto& entry : panvadere_data) {
        const std::string& table = std::get<1>(entry);
        // Find the table in the vector or create a new entry
        auto it = std::find_if(tables.begin(), tables.end(), [&table](const std::pair<std::string, uint32_t>& entry) {
            return entry.first == table;
        });
        if (it != tables.end()) {
            it->second++; // Increment count for existing table
            // Also add the person ID to the table_person_ids vector
            auto person_id = std::get<0>(entry);
            auto person_it = std::find_if(table_person_ids.begin(), table_person_ids.end(),
                                          [&table](const std::tuple<std::string, std::vector<uint32_t>>& entry) {
                                              return std::get<0>(entry) == table;
                                          });
            if (person_it != table_person_ids.end()) {
                std::get<1>(*person_it).push_back(person_id); // Add person ID
            }
            else {
                // If the table is not found, create a new entry with the person ID, this should not happen
                std::cerr << "Warning: Table " << table << " not found in table_person_ids, creating new entry."
                          << std::endl;
                table_person_ids.emplace_back(table, std::vector<uint32_t>{person_id}); // Add new entry with person ID
            }
        }
        else {
            tables.emplace_back(table, 1); // Add new table with count 1
            // Also add the person ID to the table_person_ids vector
            table_person_ids.emplace_back(table, std::vector<uint32_t>{std::get<0>(entry)});
        }
    }

    // Secondly we now assign each table a vector of household sizes, we want to assign two seats at each table to a household
    std::map<std::string, std::vector<uint32_t>> table_household_sizes;
    for (const auto& [table, count] : tables) {
        if (count <= 2) {
            // If the table has 2 or less people, we assign it to one household
            table_household_sizes[table].push_back(count);
        }
        else if (count > 2 && count <= 4) {
            // If the table has more than 2 people, we divide it into two households
            uint32_t household_size1 = 2; // First household takes 2 people
            uint32_t household_size2 = count - 2; // Second household takes the rest
            table_household_sizes[table].push_back(household_size1);
            table_household_sizes[table].push_back(household_size2);
        }
        else {
            // If the table has more than 4 people, we divide it into multiple households of 2 people each
            uint32_t num_households = count / 2;
            for (uint32_t i = 0; i < num_households; ++i) {
                table_household_sizes[table].push_back(2);
            }
            if (count % 2 != 0) {
                // Add remaining person as a smaller household
                table_household_sizes[table].push_back(count % 2);
            }
        }
    }

    // Now we search for the households in the city world and assign them to the tables
    std::map<uint32_t, uint32_t> household_map; // Maps panv id to household id in the world simulation
    std::vector<uint32_t> household_ids_reserved; // Store household IDs which have been assigned to tables
    for (const auto& table_entry : table_household_sizes) {
        const std::string& table                     = table_entry.first;
        const std::vector<uint32_t>& household_sizes = table_entry.second;
        int table_global_counter                     = 0;
        for (const auto& household_size : household_sizes) {
            // Find a household in the city world with the same size which hasnt been assigned yet
            auto household_it =
                std::find_if(city.get_locations().begin(), city.get_locations().end(),
                             [&household_size, &household_ids_reserved](const auto& location) {
                                 return location.get_persons().size() > household_size &&
                                        std::find(household_ids_reserved.begin(), household_ids_reserved.end(),
                                                  location.get_index()) == household_ids_reserved.end();
                             });
            // add to reserved household ids
            household_ids_reserved.push_back(household_it->get_index());
            if (household_it != city.get_locations().end()) {
                // If we found a household, map the table to the household id
                for (size_t i = 0; i < household_size; i++) {
                    // Get the person ID from the table_person_ids vector
                    auto person_it =
                        std::find_if(table_person_ids.begin(), table_person_ids.end(), [&table](const auto& entry) {
                            return std::get<0>(entry) == table;
                        });
                    if (person_it != table_person_ids.end() && i < std::get<1>(*person_it).size()) {
                        uint32_t person_id = std::get<1>(*person_it)[table_global_counter++];
                        household_map[person_id] =
                            household_it->get_persons().at(i).get()->get_person_id(); // Map person ID to household ID
                    }
                    else {
                        std::cerr << "Warning: Not enough persons for table " << table << " with size "
                                  << household_size << ". Skipping this person." << std::endl;
                    }
                }
            }
            else {
                // If we didn't find a household, we can either skip this table or create a new household
                std::cerr << "Warning: No household found for table " << table << " with size " << household_size
                          << ". Skipping this table." << std::endl;
            }
        }
    }

    // // Summarize the household map
    // std::cout << "Household map created with " << household_map.size() << " entries." << std::endl;
    // for (const auto& [person_id_pre, person_id_sim] : household_map) {
    //     std::cout << "Person ID Pre: " << person_id_pre << " Person ID Sim: " << person_id_sim << " HouseholdID: "
    //               << city.get_person(person_id_sim).get_assigned_location_index(mio::abm::LocationType::Home)
    //               << std::endl;
    // }

    return mio::success(household_map);
}

mio::IOResult<std::map<uint32_t, uint32_t>>
EventSimulator::map_random_restaurant_tables_to_households(mio::abm::World& city)
{
    // This funciton maps a simulation id to a household id in the world simulation.

    // First we read in the infection data from the Panvadere simulation
    BOOST_OUTCOME_TRY(auto panvadere_file,
                      get_panvadere_file_for_event_type(EventType::Restaurant_Table_Equals_Household));
    BOOST_OUTCOME_TRY(auto panvadere_data, read_panv_file_restaurant(panvadere_file));

    // For this we want to map each seat to a random household. For this we look at a table and assign the first household with that
    // amount of people to the table. If the table has more than 5 people we devide it into two households: one household for the first 5 people and one for the rest.

    // First we just need each person ID
    std::vector<uint32_t> person_ids;
    for (const auto& entry : panvadere_data) {
        uint32_t person_id = std::get<0>(entry);
        // Add person ID to the vector
        person_ids.push_back(person_id);
    }

    // Get all households in the city world
    std::vector<std::tuple<mio::abm::LocationId, uint32_t>> households;
    for (const auto& location : city.get_locations()) {
        if (location.get_type() == mio::abm::LocationType::Home) {
            households.push_back(std::make_tuple(mio::abm::LocationId({location.get_index(), location.get_type()}),
                                                 location.get_persons().size())); // Store household index
        }
    }

    //Delete all households which have less than 2 persons, as we want to assign at least two persons to a household
    households.erase(std::remove_if(households.begin(), households.end(),
                                    [](const std::tuple<mio::abm::LocationId, uint32_t>& household) {
                                        return std::get<1>(household) < 4; // Remove households with less than 2 persons
                                    }),
                     households.end());

    //Shuffle the households to get a random order
    std::shuffle(households.begin(), households.end(), city.get_rng());

    // Now we assign each person to a random household
    std::map<uint32_t, uint32_t> household_map;
    int counter_hh = 0; // Counter for households, to assign each person to a household
    for (const auto& person_id : person_ids) {
        auto household_persons   = city.get_individualized_location(std::get<0>(households[counter_hh])).get_persons();
        household_map[person_id] = household_persons[0].get()->get_person_id(); // Map person ID to household ID
        counter_hh++; // Increment household counter
    }

    return mio::success(household_map);
}

mio::IOResult<std::map<uint32_t, uint32_t>>
EventSimulator::map_work_meeting_to_households(const std::map<uint32_t, bool>& panvadere_data)
{
    std::map<uint32_t, uint32_t> household_infections;

    // TODO: Implement mapping logic
    // For work meeting: one person per household in each room
    // Distribute according to household proportion in population

    for (const auto& [person_id, is_infected] : panvadere_data) {
        if (is_infected) {
            // Each person represents a different household
            uint32_t household_id              = person_id;
            household_infections[household_id] = true;
        }
    }

    return mio::success(household_infections);
}

mio::IOResult<mio::abm::World> EventSimulator::create_event_world(const EventSimulationConfig& config,
                                                                  mio::abm::World& city,
                                                                  std::map<uint32_t, uint32_t>& event_map)
{
    // We take the persons out of the city world and create a new world for the event simulation.
    // There we just need one locaiton for the event, which is the event location.
    // All the persons are directly assigned to this location and dont move.
    // Then we need one contagious person, which is the first person in the event_map.

    auto world = mio::abm::World(num_age_groups);
    set_parameters(world.parameters);

    mio::abm::LocationId loc_id; // Location ID for the event location
    if (config.type == EventType::Restaurant_Table_Equals_Household ||
        config.type == EventType::Restaurant_Table_Equals_Half_Household) {
        loc_id = world.add_location(mio::abm::LocationType::SocialEvent); // Add one event location
    }
    else if (config.type == EventType::WorkMeeting_Many_Meetings ||
             config.type == EventType::WorkMeeting_Few_Meetings) {
        loc_id = world.add_location(mio::abm::LocationType::Work); // Add one work meeting location
    }
    // We at least need one house, hospital and icu
    auto home     = world.add_location(mio::abm::LocationType::Home); // Add one home location
    auto hospital = world.add_location(mio::abm::LocationType::Hospital); // Add one hospital
    auto icu      = world.add_location(mio::abm::LocationType::ICU); // Add one ICU location

    // Add persons to the event world
    for (const auto& [person_id_before, person_id_sim] : event_map) {
        // Get the person from the city world
        auto& person = city.get_person(person_id_sim);

        // Create a new person in the event world
        auto& added_person = world.add_person(loc_id, person.get_age());

        // assign the person to the event location
        added_person.set_assigned_location(home);
        added_person.set_assigned_location(hospital);
        added_person.set_assigned_location(icu);
        added_person.set_assigned_location(loc_id); // Assign to the event location
    }

    // Set the first person as infected
    if (!event_map.empty()) {
        auto& first_person = world.get_persons()[0]; // Get the first person in the event world
        auto prng          = mio::abm::Person::RandomNumberGenerator(world.get_rng(), first_person);
        first_person.add_new_infection(mio::abm::Infection(prng, mio::abm::VirusVariant::Alpha, first_person.get_age(),
                                                           world.parameters, mio::abm::TimePoint(0),
                                                           mio::abm::InfectionState::InfectedNoSymptoms));
    }

    return mio::success(std::move(world));
}

mio::IOResult<std::map<uint32_t, uint32_t>> EventSimulator::map_events_to_persons(const mio::abm::World& city,
                                                                                  EventType event_type)
{
    switch (event_type) {
    case EventType::Restaurant_Table_Equals_Household: {
        BOOST_OUTCOME_TRY(auto household_map,
                          EventSimulator::map_restaurant_tables_to_households(const_cast<mio::abm::World&>(city)));
        return mio::success(household_map);
    }
    case EventType::Restaurant_Table_Equals_Half_Household: {
        BOOST_OUTCOME_TRY(auto household_map,
                          EventSimulator::map_two_person_per_household_restaurant_tables_to_households(
                              const_cast<mio::abm::World&>(city)));
        return mio::success(household_map);
    }
    case EventType::Restaurant_Table_Equals_Random: {
        BOOST_OUTCOME_TRY(auto household_map, EventSimulator::map_random_restaurant_tables_to_households(
                                                  const_cast<mio::abm::World&>(city)));
        return mio::success(household_map);
    }
    case EventType::WorkMeeting_Many_Meetings:
        return mio::failure(mio::StatusCode::InvalidValue, "WorkMeeting_Many_Meetings is not implemented yet");
    case EventType::WorkMeeting_Few_Meetings:
        return mio::failure(mio::StatusCode::InvalidValue, "WorkMeeting_Few_Meetings is not implemented yet");
    default:
        std::cerr << "Error: Unknown event type " << EventSimulator::event_type_to_string(event_type) << std::endl;
        // Return an error result
        return mio::failure(mio::StatusCode::InvalidValue, "Unknown event type");
    }
}