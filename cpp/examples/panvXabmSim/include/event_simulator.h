#pragma once
#include "abm/abm.h"
#include "memilio/io/result_io.h"
#include "simulation_runner.h" // For SimType enum
#include <string>
#include <map>
#include <unordered_map>

enum class EventType
{
    Restaurant_Table_Equals_Household,
    Restaurant_Table_Equals_Half_Household,
    WorkMeeting_Many_Meetings,
    WorkMeeting_Few_Meetings
};

// Mapping from EventType to corresponding Panvadere data file
const std::unordered_map<EventType, std::string> EVENT_TYPE_TO_FILE = {
    {EventType::Restaurant_Table_Equals_Household, "./data/restaurant_scenario/infections.txt"},
    {EventType::Restaurant_Table_Equals_Half_Household, "./data/restaurant_scenario/infections.txt"},
    {EventType::WorkMeeting_Many_Meetings, ""},
    {EventType::WorkMeeting_Few_Meetings, ""}};

struct EventSimulationConfig {
    EventType type;
    double infection_parameter_k = 1.0;
    int event_duration_hours     = 2;

    // Get the appropriate panvadere file for this event type
    std::string get_panvadere_file() const;
};

class EventSimulator
{
public:
    static mio::IOResult<double> calculate_infection_parameter_k(const EventSimulationConfig& config);

    static mio::IOResult<std::map<uint32_t, bool>> initialize_from_panvadere(EventType event_type);

    static mio::IOResult<std::map<uint32_t, bool>> initialize_from_event_simulation(const EventSimulationConfig& config,
                                                                                    const mio::abm::World& city);

    static mio::IOResult<std::map<uint32_t, bool>> map_events_to_persons(const mio::abm::World& city,
                                                                         EventType event_type);

    static std::string event_type_to_string(EventType type);

    static std::string simulation_type_to_string(SimType type);

    // Get the panvadere file path for a specific event type
    static mio::IOResult<std::string> get_panvadere_file_for_event_type(EventType type);

private:
    static std::map<uint32_t, bool> map_restaurant_tables_to_households(const std::map<uint32_t, bool>& panvadere_data);
    static std::map<uint32_t, bool> map_work_meeting_to_households(const std::map<uint32_t, bool>& panvadere_data);
    static std::map<uint32_t, bool> map_choir_to_households(const std::map<uint32_t, bool>& panvadere_data);

    static mio::IOResult<mio::abm::World> create_event_world(const EventSimulationConfig& config);
    static std::map<uint32_t, bool> simulate_event_transmission(const mio::abm::World& event_world,
                                                                const EventSimulationConfig& config,
                                                                const mio::abm::World& city);
};