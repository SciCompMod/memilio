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
    if (!config.is_valid()) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Invalid configuration: infection_parameter_k must be > 0 and <= 10, "
                            "event_duration_hours must be > 0 and <= 24");
    }

    // Calculate K parameter based on event type and duration
    double base_k_value = 1.0;

    switch (config.type) {
    case EventType::Restaurant_Table_Equals_Household:
        base_k_value = 1.2; // Higher transmission in close dining
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

mio::IOResult<std::map<uint32_t, bool>> EventSimulator::read_infection_data(const std::string& filename)
{
    std::map<uint32_t, bool> infected_status;

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

    // Skip header if present
    if (std::getline(file, line)) {
        line_number++;
    }

    while (std::getline(file, line)) {
        line_number++;
        if (line.empty())
            continue; // Skip empty lines

        std::istringstream iss(line);
        uint32_t id;
        std::string table;
        int infected;

        if (iss >> id >> table >> infected) {
            if (infected < 0 || infected > 1) {
                std::cerr << "Warning: Invalid infection status " << infected << " at line " << line_number
                          << ", treating as not infected" << std::endl;
                infected_status[id] = false;
            }
            else {
                infected_status[id] = (infected == 1);
            }
        }
        else {
            std::cerr << "Warning: Failed to parse line " << line_number << ": '" << line << "'" << std::endl;
        }
    }

    file.close();

    if (infected_status.empty()) {
        std::cerr << "Warning: No valid infection data found in file: " << filename << std::endl;
    }

    return mio::success(infected_status);
}

mio::IOResult<std::map<uint32_t, bool>> EventSimulator::initialize_from_panvadere(EventType event_type)
{
    // Get the appropriate file for this event type
    BOOST_OUTCOME_TRY(auto panvadere_file, get_panvadere_file_for_event_type(event_type));

    // Read Panvadere infection data
    BOOST_OUTCOME_TRY(auto panvadere_data, read_infection_data(panvadere_file));

    // Map to households based on event type
    std::map<uint32_t, bool> household_infections;

    switch (event_type) {
    case EventType::Restaurant_Table_Equals_Household:
    case EventType::Restaurant_Table_Equals_Half_Household:
        household_infections = map_restaurant_tables_to_households(panvadere_data);
        break;
    case EventType::WorkMeeting_Many_Meetings:
    case EventType::WorkMeeting_Few_Meetings:
        household_infections = map_work_meeting_to_households(panvadere_data);
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
        return "Restaurant (Table = Household)";
    case EventType::Restaurant_Table_Equals_Half_Household:
        return "Restaurant (Table = Half Household)";
    case EventType::WorkMeeting_Many_Meetings:
        return "Work Meeting (Many Meetings)";
    case EventType::WorkMeeting_Few_Meetings:
        return "Work Meeting (Few Meetings)";
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

std::map<uint32_t, bool>
EventSimulator::map_restaurant_tables_to_households(const std::map<uint32_t, bool>& panvadere_data)
{
    std::map<uint32_t, bool> household_infections;

    // TODO: Implement mapping logic
    // For restaurant: each table represents a household
    // Interchange households - half of each household at each table

    for (const auto& [person_id, is_infected] : panvadere_data) {
        if (is_infected) {
            // Map person to household (simple: person_id / 4 = household_id)
            uint32_t household_id              = person_id / 4;
            household_infections[household_id] = true;
        }
    }

    return household_infections;
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