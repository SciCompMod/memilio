#include "../include/event_simulator.h"
#include "../include/constants.h"
#include <fstream>
#include <iostream>
#include <sstream>

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

mio::IOResult<double> EventSimulator::calculate_infection_parameter_k(const EventSimulationConfig& config)
{
    // Validate configuration first

    // Calculate K parameter based on event type and duration
    double base_k_value = 1.0;

    switch (config.type) {
    case EventType::Restaurant_Table_Equals_Household:
        base_k_value = 20.0; // Higher transmission in close dining
        break;
    case EventType::Restaurant_Table_Equals_Half_Household:
        base_k_value = 1.1; // Slightly lower than full household tables
        break;
    case EventType::WorkMeeting_Many_Meetings:
        base_k_value = 0.8; // Lower transmission in meeting room
        break;
    case EventType::WorkMeeting_Few_Meetings:
        base_k_value = 0.9; // Slightly higher with fewer, longer meetings
        break;
    }

    // Adjust K value based on event duration (longer events = higher transmission)
    double duration_factor  = 1.0 + (config.event_duration_hours - 2) * 0.1; // 10% increase per extra hour
    double adjusted_k_value = base_k_value * duration_factor * config.infection_parameter_k;

    return mio::success(adjusted_k_value);
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

    std::cout << "Successfully read " << infected_status.size() << " persons from " << filename << std::endl;

    return mio::success(infected_status);
}

mio::IOResult<std::map<uint32_t, bool>> EventSimulator::initialize_from_panvadere(EventType event_type)
{
    // Get the appropriate file for this event type
    BOOST_OUTCOME_TRY(auto panvadere_file, get_panvadere_file_for_event_type(event_type));

    // Read Panvadere infection data
    BOOST_OUTCOME_TRY(auto panvadere_data, read_panv_file_restaurant(panvadere_file));

    // Map to households based on event type
    std::map<uint32_t, bool> household_infections;

    switch (event_type) {
    case EventType::Restaurant_Table_Equals_Household:
    case EventType::Restaurant_Table_Equals_Half_Household:
    case EventType::WorkMeeting_Many_Meetings:
    case EventType::WorkMeeting_Few_Meetings:
    default:
        // Map events to households
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Event type not supported for Panvadere data mapping: " + event_type_to_string(event_type));
        break;
    }

    return mio::success(household_infections);
}

mio::IOResult<std::map<uint32_t, bool>>
EventSimulator::initialize_from_event_simulation(const EventSimulationConfig& config, const mio::abm::World& city)
{
    // Create event-specific world
    BOOST_OUTCOME_TRY(auto event_world, create_event_world(config));

    // Simulate event transmission
    auto infected_people = simulate_event_transmission(config, city);

    return mio::success(infected_people);
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

    std::map<uint32_t, uint32_t> household_map;

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
    std::map<uint32_t, std::vector<std::tuple<uint32_t, uint32_t>>> household_map; // Maps a table to
    for (const auto& [table, household_sizes] : table_household_sizes) {
        for (const auto& household_size : household_sizes) {
            // Find a household in the city world with the same size
            auto household_it = std::find_if(city.get_locations().begin(), city.get_locations().end(),
                                             [&household_size](const auto& location) {
                                                 return location.get_persons().size() == household_size;
                                             });
            if (household_it != city.get_locations().end()) {
                // If we found a household, map the table to the household id
                uint32_t household_id       = household_it->get_index(); // Get the household id from the location index
                household_map[household_id] = table; // Map household id to table
            }
            else {
                // If we didn't find a household, we can either skip this table or create a new household
                std::cerr << "Warning: No household found for table " << table << " with size " << household_size
                          << ". Skipping this table." << std::endl;
            }
        }
    }

    mio::unused(city); // Prevent unused warning for the city world
    return mio::success(household_map);
}

std::map<uint32_t, bool> EventSimulator::map_work_meeting_to_households(const std::map<uint32_t, bool>& panvadere_data)
{
    std::map<uint32_t, bool> household_infections;

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

    return household_infections;
}

std::map<uint32_t, bool> EventSimulator::map_choir_to_households(const std::map<uint32_t, bool>& panvadere_data)
{
    // This function is kept for potential future use but not currently used
    // with the current EventType enum
    std::map<uint32_t, bool> household_infections;

    for (const auto& [person_id, is_infected] : panvadere_data) {
        if (is_infected) {
            uint32_t household_id              = person_id / 2; // Families might attend together
            household_infections[household_id] = true;
        }
    }

    return household_infections;
}

mio::IOResult<mio::abm::World> EventSimulator::create_event_world(const EventSimulationConfig& config)
{
    // TODO: Create a specialized world for event simulation
    // This would be a smaller world focused on the specific event type
    mio::unused(config); // Prevent unused warning
    auto world = mio::abm::World(num_age_groups);

    // Add event location
    // auto event_location = world.add_location(mio::abm::LocationType::SocialEvent);

    // Add people to event
    // for (int i = 0; i < config.num_persons; ++i) {
    //     auto person_id = size_t(i);
    //     auto age_group = mio::AgeGroup(i % num_age_groups);
    //     auto person =
    //         mio::abm::Person(world.get_individualized_location(mio::abm::LocationType::Home, person_id), age_group);
    //     world.add_person(person_id, std::move(person));
    // }

    return mio::success(std::move(world));
}

std::map<uint32_t, bool> EventSimulator::simulate_event_transmission(const EventSimulationConfig& config,
                                                                     const mio::abm::World& event_world)
{
    std::map<uint32_t, bool> infected_people;

    mio::unused(event_world, config); // Prevent unused warning

    // TODO: Run event simulation to determine who gets infected
    // 1. Place one initially infected person
    // 2. Run simulation for event duration
    // 3. Return list of infected people

    // For now, simple random infection
    // for (const auto& [person_id, person] : event_world.get_persons()) {
    //     // Simple placeholder: 10% infection rate
    //     if (person_id.get() % 10 == 0) {
    //         infected_people[person_id.get()] = true;
    //     }
    // }

    return infected_people;
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
    case EventType::Restaurant_Table_Equals_Half_Household:
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Restaurant_Table_Equals_Half_Household is not implemented yet");
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