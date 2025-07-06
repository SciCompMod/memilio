#pragma once
#include "memilio/io/result_io.h"
#include <string>
#include <map>
#include <unordered_map>
#include "abm/world.h"
#include "abm/location_type.h"
#include "config_validation.h"

enum class EventType
{
    Restaurant_Table_Equals_Household,
    Restaurant_Table_Equals_Half_Household,
    WorkMeeting_Many_Meetings,
    WorkMeeting_Few_Meetings
};

enum class SimType
{
    Panvadere,
    Memilio,
    Both
};

// Mapping from EventType to corresponding Panvadere data file
const std::unordered_map<EventType, std::string> EVENT_TYPE_TO_FILE = {
    {EventType::Restaurant_Table_Equals_Household,
     std::string(Config::DEFAULT_BASE_DIR) + "/data/restaurant_scenario/infections.txt"},
    {EventType::Restaurant_Table_Equals_Half_Household,
     std::string(Config::DEFAULT_BASE_DIR) + "/data/restaurant_scenario/infections.txt"},
    {EventType::WorkMeeting_Many_Meetings,
     std::string(Config::DEFAULT_BASE_DIR) + "/data/work_meeting_many/infections.txt"},
    {EventType::WorkMeeting_Few_Meetings,
     std::string(Config::DEFAULT_BASE_DIR) + "/data/work_meeting_few/infections.txt"}};

struct EventSimulationConfig {
    EventType type;
    double infection_parameter_k = 1.0;
    int event_duration_hours     = 2;

    // Get the appropriate panvadere file for this event type
    std::string get_panvadere_file() const;
};

/**
 * @brief EventSimulator class for handling event-based infection simulations
 * 
 * This class provides functionality to:
 * - Initialize infections from Panvadere simulation data
 * - Calculate infection parameters for different event types
 * - Map event participants to household infections
 */
class EventSimulator
{
public:
    // === Configuration and Validation ===

    /**
     * @brief Calculate infection parameter K for a specific event configuration
     * @param config Event simulation configuration
     * @return IOResult containing the calculated K parameter
     */
    static mio::IOResult<double> calculate_infection_parameter_k(const EventSimulationConfig& config);

    // === Initialization Methods ===

    /**
     * @brief Initialize infections from Panvadere simulation data
     * @param event_type Type of event to process
     * @return IOResult containing household infection mapping
     */
    static mio::IOResult<std::map<uint32_t, bool>> initialize_from_panvadere(EventType event_type);

    /**
     * @brief Initialize infections from event simulation
     * @param config Event simulation configuration
     * @param city The city world for context
     * @return IOResult containing infected person mapping
     */
    static mio::IOResult<std::map<uint32_t, bool>> initialize_from_event_simulation(const EventSimulationConfig& config,
                                                                                    const mio::abm::World& city);

    // === Utility Functions ===

    /**
     * @brief Map events to persons (placeholder for future implementation)
     */
    static mio::IOResult<std::map<uint32_t, uint32_t>> map_events_to_persons(const mio::abm::World& city,
                                                                             EventType event_type);

    /**
     * @brief Convert EventType to human-readable string
     */
    static std::string event_type_to_string(EventType type);

    /**
     * @brief Convert SimType to human-readable string
     */
    static std::string simulation_type_to_string(SimType type);

    /**
     * @brief Get the panvadere file path for a specific event type
     */
    static mio::IOResult<std::string> get_panvadere_file_for_event_type(EventType type);

    /**
     * @brief Read infection data from a Panvadere file
     */
    static mio::IOResult<std::vector<std::tuple<uint32_t, std::string, bool>>>
    read_panv_file_restaurant(const std::string& filename);

    /**
     * @brief Read infection data from a Panvadere file
     */
    static mio::IOResult<std::map<uint32_t, bool>> read_panv_file_work(const std::string& filename);

    /**
     * @brief Map restaurant tables to households for infection tracking
     */
    static mio::IOResult<std::map<uint32_t, uint32_t>> map_restaurant_tables_to_households(mio::abm::World& city);

private:
    static std::map<uint32_t, bool> map_work_meeting_to_households(const std::map<uint32_t, bool>& panvadere_data);
    static std::map<uint32_t, bool> map_choir_to_households(const std::map<uint32_t, bool>& panvadere_data);

    static mio::IOResult<mio::abm::World> create_event_world(const EventSimulationConfig& config);
    static std::map<uint32_t, bool> simulate_event_transmission(const EventSimulationConfig& config,
                                                                const mio::abm::World& city);
};