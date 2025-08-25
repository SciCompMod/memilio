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
    std::vector<std::tuple<uint32_t, uint32_t, bool>> panvadere_data_work;
    std::vector<std::tuple<uint32_t, std::string, bool>> panvadere_data_restaurant;
    int event_duration_hours;
    if (config.type == EventType::Restaurant_Table_Equals_Household ||
        config.type == EventType::Restaurant_Table_Equals_Half_Household ||
        config.type == EventType::Restaurant_Table_Equals_Random) {
        // For restaurant events, read the restaurant-specific Panvadere file
        BOOST_OUTCOME_TRY(panvadere_data_restaurant, read_panv_file_restaurant(panvadere_file));
        event_duration_hours = 2; // Default event duration for restaurant events
    }
    else if (config.type == EventType::WorkMeeting_Many_Meetings ||
             config.type == EventType::WorkMeeting_Baseline_Meetings) {
        // For work meeting events, read the work-specific Panvadere file
        BOOST_OUTCOME_TRY(panvadere_data_work, read_panv_file_work(panvadere_file));
        event_duration_hours = 8; // Default event duration for work meetings
    }
    else {
        std::cerr << "Error: Unsupported event type for Panvadere data reading: "
                  << EventSimulator::event_type_to_string(config.type) << std::endl;
        return mio::failure(mio::StatusCode::InvalidValue, "Unsupported event type for Panvadere data reading: " +
                                                               EventSimulator::event_type_to_string(config.type));
    }

    // Count infected persons in the Panvadere data
    size_t infected_count = 0;
    if (config.type == EventType::Restaurant_Table_Equals_Household ||
        config.type == EventType::Restaurant_Table_Equals_Half_Household ||
        config.type == EventType::Restaurant_Table_Equals_Random) {
        // For restaurant events, count infected persons in the restaurant data
        for (const auto& entry : panvadere_data_restaurant) {
            if (std::get<2>(entry)) { // Check if infected status is true
                infected_count++;
            }
        }
    }
    else if (config.type == EventType::WorkMeeting_Many_Meetings ||
             config.type == EventType::WorkMeeting_Baseline_Meetings) {
        for (const auto& entry : panvadere_data_work) {
            if (std::get<2>(entry)) { // Check if infected status is true
                infected_count++;
            }
        }
    }

    // Now we need to simulate the event world and search for a suitable K parameter with which the simulation has the same number of infected persons as in the Panvadere data

    // Rudimentary grid search for K parameter
    double k_min = 1.0; // Minimum K value
    double k_max = 50.0; // Maximum K value

    double k_step                  = 0.2; // Step size for K value
    double best_k                  = k_min;
    size_t best_infected_count     = 0;
    const int num_runs             = 100; // Number of runs for averaging
    double best_avg_infected_count = -1.0;

    for (double k = k_min; k <= k_max; k += k_step) {
        size_t total_infected_for_k = 0;
        for (int i = 0; i < num_runs; ++i) {

            BOOST_OUTCOME_TRY(auto calculation_world, create_event_world(config, city, event_map));
            auto run_rng_counter =
                mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(i), mio::Counter<uint32_t>(0));
            calculation_world.get_rng().set_counter(run_rng_counter);

            set_parameters(calculation_world.parameters);

            for (auto& person : calculation_world.get_persons()) {
                // Reset the infection state for each person
                person.get_rng_counter() = mio::Counter<uint32_t>(i * 10000 + i); // Ensure unique counter for each run
            }

            double ratio = 1.0;
            if (config.type == EventType::Restaurant_Table_Equals_Half_Household ||
                config.type == EventType::Restaurant_Table_Equals_Household ||
                config.type == EventType::Restaurant_Table_Equals_Random) {
                // Set restaurant-specific parameters
                ratio = 4.0 * 4.0 / event_duration_hours; // 4 hours divided by event duration
            }
            else if (config.type == EventType::WorkMeeting_Many_Meetings ||
                     config.type == EventType::WorkMeeting_Baseline_Meetings) {
                ratio = 7.0 / event_duration_hours; // 7 hours divided by event duration
            }
            set_local_parameters_event(calculation_world, ratio);

            calculation_world.use_migration_rules(false); // Disable migration rules for this simulation
            auto t0 = mio::abm::TimePoint(0); // Start time per simulation

            auto tmax = mio::abm::TimePoint(0) + mio::abm::hours(event_duration_hours); // End time per simulation

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
            std::cout << "Found better K: " << best_k << " with average infected count: " << best_avg_infected_count
                      << " (target: " << infected_count << ")" << std::endl;
        }
    }

    return mio::success(best_k);
}

mio::IOResult<std::vector<std::tuple<uint32_t, uint32_t, bool>>>
EventSimulator::read_panv_file_work(const std::string& filename)
{
    std::vector<std::tuple<uint32_t, uint32_t, bool>> infected_status;

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
        if (line.find("pedestrianId") == std::string::npos || line.find("workRoomId") == std::string::npos ||
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
        std::string pedestrian_id_str, workroom_id_int, infected_str;

        // Extract fields using tab delimiter
        if (std::getline(iss, pedestrian_id_str, '\t') && std::getline(iss, workroom_id_int, '\t') &&
            std::getline(iss, infected_str)) {

            try {
                // Parse pedestrian ID
                uint32_t pedestrian_id = static_cast<uint32_t>(std::stoul(pedestrian_id_str));

                // Parse workroom (remove any trailing whitespace)
                workroom_id_int.erase(workroom_id_int.find_last_not_of(" \t\r\n") + 1);
                uint32_t workroom_id = static_cast<uint32_t>(std::stoi(workroom_id_int));

                // Parse infection status
                infected_str.erase(infected_str.find_last_not_of(" \t\r\n") + 1);
                int infected_int = std::stoi(infected_str);

                if (infected_int < 0 || infected_int > 1) {
                    std::cerr << "Warning: Invalid infection status " << infected_int << " for pedestrian "
                              << pedestrian_id << " at line " << line_number << ", treating as not infected"
                              << std::endl;
                    infected_status.emplace_back(pedestrian_id, workroom_id, false);
                }
                else {
                    bool is_infected = (infected_int == 1);
                    infected_status.emplace_back(pedestrian_id, workroom_id, is_infected);
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
EventSimulator::initialize_from_event_simulation(EventType event_type, std::map<uint32_t, uint32_t>& event_map)
{
    // Map to households based on event type
    std::vector<uint32_t> household_infections;

    switch (event_type) {
    case EventType::Restaurant_Table_Equals_Half_Household: {
        std::vector<uint32_t> return_vector;
        return_vector.reserve(9); // Reserve space for 9 households
        return_vector = {1, 10, 19, 28, 37, 47, 55, 64, 89};
        for (size_t i = 0; i < return_vector.size(); ++i) {
            household_infections.push_back(event_map[return_vector[i]]); // Map to household IDs
        }

        break;
    }
    case EventType::Restaurant_Table_Equals_Random: {
        std::vector<uint32_t> return_vector;
        return_vector.reserve(9); // Reserve space for 9 households
        return_vector = {1, 10, 19, 28, 37, 47, 55, 64, 89};
        for (size_t i = 0; i < return_vector.size(); ++i) {
            household_infections.push_back(event_map[return_vector[i]]); // Map to household IDs
        }

        break;
    }

    case EventType::Restaurant_Table_Equals_Household: {
        std::vector<uint32_t> return_vector;
        return_vector.reserve(9); // Reserve space for 9 households
        return_vector = {1, 10, 19, 28, 37, 47, 55, 64, 89};
        for (size_t i = 0; i < return_vector.size(); ++i) {
            household_infections.push_back(event_map[return_vector[i]]); // Map to household IDs
        }

        break;
    }
    case EventType::WorkMeeting_Many_Meetings: {
        std::vector<uint32_t> return_vector;
        return_vector.reserve(6); // Reserve space for 6 households
        return_vector = {1, 5, 9, 14, 18, 22};
        for (size_t i = 0; i < return_vector.size(); ++i) {
            household_infections.push_back(event_map[return_vector[i]]); // Map to household IDs
        }
        break;
    }
    case EventType::WorkMeeting_Baseline_Meetings: {
        std::vector<uint32_t> return_vector;
        return_vector.reserve(3); // Reserve space for 3 households
        return_vector = {1, 12, 23};
        for (size_t i = 0; i < return_vector.size(); ++i) {
            household_infections.push_back(event_map[return_vector[i]]); // Map to household IDs
        }
        break;
    }
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
    BOOST_OUTCOME_TRY(auto panvadere_file, get_panvadere_file_for_event_type(event_type));

    switch (event_type) {
    case EventType::Restaurant_Table_Equals_Household:
    case EventType::Restaurant_Table_Equals_Random:
    case EventType::Restaurant_Table_Equals_Half_Household: {
        // Read Panvadere infection data
        BOOST_OUTCOME_TRY(auto panvadere_data, read_panv_file_restaurant(panvadere_file));
        for (const auto& entry : panvadere_data) {
            uint32_t person_id = std::get<0>(entry);
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
                    return mio::failure(mio::StatusCode::InvalidValue, "Person ID " + std::to_string(person_id));
                }
            }
        }
        break;
    }
    case EventType::WorkMeeting_Many_Meetings:
    case EventType::WorkMeeting_Baseline_Meetings: {
        // Read Panvadere infection data
        BOOST_OUTCOME_TRY(auto panvadere_data, read_panv_file_work(panvadere_file));
        for (const auto& entry : panvadere_data) {
            uint32_t person_id = std::get<0>(entry);
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
                                                                           " not found in event map for table ");
                }
            }
        }
        break;
    }
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
    case EventType::WorkMeeting_Baseline_Meetings:
        return "WorkMeetingBaseline";
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
                                 return location.get_persons().size() >= household_size &&
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

    std::cout << "Household map created with " << household_map.size() << " entries." << std::endl;
    for (const auto& [person_id_pre, person_id_sim] : household_map) {
        std::cout << "Person ID Pre: " << person_id_pre << " Person ID Sim: " << person_id_sim << " HouseholdID: "
                  << city.get_person(person_id_sim).get_assigned_location_index(mio::abm::LocationType::Home)
                  << "House has "
                  << city
                         .get_individualized_location(mio::abm::LocationId{
                             city.get_person(person_id_sim).get_assigned_location_index(mio::abm::LocationType::Home),
                             mio::abm::LocationType::Home})
                         .get_persons()
                         .size()
                  << " persons." << std::endl;
    }
    return mio::success(household_map);
}

mio::IOResult<std::map<uint32_t, uint32_t>>
EventSimulator::map_three_person_per_household_restaurant_tables_to_households(mio::abm::World& city)
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
        if (count <= 3) {
            // If the table has 3 or less people, we assign it to one household
            table_household_sizes[table].push_back(count);
        }
        else if (count > 3 && count <= 4) {
            // If the table has more than 3 people, we divide it into two households
            uint32_t household_size1 = 2; // First household takes 2 people
            uint32_t household_size2 = count - 2; // Second household takes the rest
            table_household_sizes[table].push_back(household_size1);
            table_household_sizes[table].push_back(household_size2);
        }
        else {
            // If the table has more than 6 people, we divide it into multiple households of 2 people each
            uint32_t num_households = count / 3;
            for (uint32_t i = 0; i < num_households; ++i) {
                table_household_sizes[table].push_back(3);
            }
            if (count % 3 != 0) {
                // Add remaining person as a smaller household
                table_household_sizes[table].push_back(count % 3);
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
                                 return location.get_persons().size() >= household_size &&
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
    std::cout << "Household map created with " << household_map.size() << " entries." << std::endl;
    for (const auto& [person_id_pre, person_id_sim] : household_map) {
        std::cout << "Person ID Pre: " << person_id_pre << " Person ID Sim: " << person_id_sim << " HouseholdID: "
                  << city.get_person(person_id_sim).get_assigned_location_index(mio::abm::LocationType::Home)
                  << std::endl;
    }

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
                                        return std::get<1>(household) < 2; // Remove households with less than 2 persons
                                    }),
                     households.end());

    // Now we assign each person to a random household
    std::map<uint32_t, uint32_t> household_map;
    int counter_hh = 0; // Counter for households, to assign each person to a household
    for (const auto& person_id : person_ids) {
        auto household_persons   = city.get_individualized_location(std::get<0>(households[counter_hh])).get_persons();
        household_map[person_id] = household_persons[0].get()->get_person_id(); // Map person ID to household ID
        counter_hh++; // Increment household counter
    }

    // std::cout << "Household map created with " << household_map.size() << " entries." << std::endl;
    // for (const auto& [person_id_pre, person_id_sim] : household_map) {
    //     std::cout << "Person ID Pre: " << person_id_pre << " Person ID Sim: " << person_id_sim << " HouseholdID: "
    //               << city.get_person(person_id_sim).get_assigned_location_index(mio::abm::LocationType::Home)
    //               << "House has "
    //               << city
    //                      .get_individualized_location(mio::abm::LocationId{
    //                          city.get_person(person_id_sim).get_assigned_location_index(mio::abm::LocationType::Home),
    //                          mio::abm::LocationType::Home})
    //                      .get_persons()
    //                      .size()
    //               << " persons." << std::endl;
    // }

    return mio::success(household_map);
}

mio::IOResult<std::map<uint32_t, uint32_t>> EventSimulator::map_work_meeting_to_households(mio::abm::World& city,
                                                                                           EventType event_type)
{
    std::map<uint32_t, uint32_t> household_infections;
    // For this we divide the persons in the simualtion event to the same workplaces. E.g. every 4 persons work at the same workplace.
    // For now we assign 1 workteam 2 rooms, this means that 1,2,9,10,11 = 1 team, 3,4,12,13,14 = 2 team, 5,6,15,16,17 = 3 team, 7,8,18,19,20 = 4 team, 21,22,23,24,25,26 = 5 team

    std::string panvadere_file;
    switch (event_type) // Use a switch statement to handle different event types
    {
    case EventType::WorkMeeting_Many_Meetings: {
        BOOST_OUTCOME_TRY(auto panvadere_file_read,
                          get_panvadere_file_for_event_type(EventType::WorkMeeting_Many_Meetings));
        panvadere_file = panvadere_file_read;
        break;
    }
    case EventType::WorkMeeting_Baseline_Meetings: {
        BOOST_OUTCOME_TRY(auto panvadere_file_read,
                          get_panvadere_file_for_event_type(EventType::WorkMeeting_Baseline_Meetings));
        panvadere_file = panvadere_file_read;
        break;
    }
    default:
        break;
    }

    std::vector<uint32_t> work_teams = {0, 1, 2, 3, 4}; // Work teams for the event simulation
    std::vector<std::vector<uint32_t>> work_persons(work_teams.size());
    work_persons[0] = {1, 2, 3, 4}; // Work team 1 Room 1 and 2
    work_persons[1] = {5, 6, 7, 8}; // Work team 2 Room 3 and 4
    work_persons[2] = {9, 10, 11, 12, 13, 14}; // Work team 3 Room 5 and 6
    work_persons[3] = {15, 16, 17, 18, 19, 20}; // Work team 4 Room 7 and 8
    work_persons[4] = {21, 22, 23, 24, 25, 26}; // Work team 5 Room 9 and 10
    std::vector<uint32_t> workplace_id_per_work_team(work_teams.size()); // Workplace IDs for each work team

    int current_work_team = 0; // Current work team index
    for (auto locations : city.get_locations()) {
        if (current_work_team >= (int)work_teams.size()) {
            break; // Reset to the first work team if we exceed the number of work teams
        }
        if (locations.get_type() == mio::abm::LocationType::Work) {
            workplace_id_per_work_team[current_work_team++] =
                locations.get_index(); // Assign the workplace ID to the work team
        }
    }

    std::vector<std::vector<uint32_t>> work_persons_global(work_teams.size());

    // Now we go through the populatoin, and assign the first persons which have a workplace index to the work team
    for (auto person : city.get_persons()) {
        auto workplace_index =
            person.get_assigned_location_index(mio::abm::LocationType::Work); // Get the workplace index of the person
        if (std::find(workplace_id_per_work_team.begin(), workplace_id_per_work_team.end(), workplace_index) !=
            workplace_id_per_work_team.end()) {
            auto index_of_work_team = std::distance(
                workplace_id_per_work_team.begin(),
                std::find(workplace_id_per_work_team.begin(), workplace_id_per_work_team.end(), workplace_index));
            // Add the person to the work team if there are not already 5 persons assigned to the work team
            if (work_persons_global[index_of_work_team].size() < work_persons[index_of_work_team].size()) {
                work_persons_global[index_of_work_team].push_back(person.get_person_id());
            }
        }
    }
    for (size_t i = 0; i < work_persons_global.size(); ++i) {
        // Now we map the work persons to the household infections
        for (size_t j = 0; j < work_persons_global[i].size(); ++j) {
            household_infections[work_persons[i][j]] =
                work_persons_global[i][j]; // Map the person ID to the work team ID
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
    world.get_rng().seed(city.get_rng().get_seeds()); // Copy the random number generator seeds from the city world
    set_parameters(world.parameters);

    mio::abm::LocationId loc_id; // Location ID for the event location
    if (config.type == EventType::Restaurant_Table_Equals_Household ||
        config.type == EventType::Restaurant_Table_Equals_Half_Household ||
        config.type == EventType::Restaurant_Table_Equals_Random) {
        loc_id = world.add_location(mio::abm::LocationType::SocialEvent); // Add one event location
    }
    else if (config.type == EventType::WorkMeeting_Many_Meetings ||
             config.type == EventType::WorkMeeting_Baseline_Meetings) {
        loc_id = world.add_location(mio::abm::LocationType::Work); // Add one work meeting location
    }
    // We at least need one house, hospital and icu
    auto home     = world.add_location(mio::abm::LocationType::Home); // Add one home location
    auto hospital = world.add_location(mio::abm::LocationType::Hospital); // Add one hospital
    auto icu      = world.add_location(mio::abm::LocationType::ICU); // Add one ICU location

    // Add persons to the event world
    // bool first_person_age = true; // Flag to check if we are adding the first person
    for (const auto& [person_id_before, person_id_sim] : event_map) {
        mio::unused(person_id_before, person_id_sim); // Avoid unused variable warning
        // Create a new person in the event world
        auto& added_person = world.add_person(loc_id, age_group_35_to_59);

        // assign the person to the event location
        added_person.set_assigned_location(home);
        added_person.set_assigned_location(hospital);
        added_person.set_assigned_location(icu);
        added_person.set_assigned_location(loc_id); // Assign to the event location
    }

    // Set the first person as infected
    if (!event_map.empty()) {
        auto& first_person = world.get_persons()[0];
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
                          EventSimulator::map_three_person_per_household_restaurant_tables_to_households(
                              const_cast<mio::abm::World&>(city)));
        return mio::success(household_map);
    }
    case EventType::Restaurant_Table_Equals_Random: {
        BOOST_OUTCOME_TRY(auto household_map, EventSimulator::map_random_restaurant_tables_to_households(
                                                  const_cast<mio::abm::World&>(city)));
        return mio::success(household_map);
    }
    case EventType::WorkMeeting_Many_Meetings:
    case EventType::WorkMeeting_Baseline_Meetings: {
        BOOST_OUTCOME_TRY(auto household_map, EventSimulator::map_work_meeting_to_households(
                                                  const_cast<mio::abm::World&>(city), event_type));
        return mio::success(household_map);
    }
    default:
        std::cerr << "Error: Unknown event type " << EventSimulator::event_type_to_string(event_type) << std::endl;
        // Return an error result
        return mio::failure(mio::StatusCode::InvalidValue, "Unknown event type");
    }
}