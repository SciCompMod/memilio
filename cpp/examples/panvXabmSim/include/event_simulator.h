#pragma once
#include "abm/abm.h"
#include "memilio/io/result_io.h"
#include <string>
#include <map>

enum class EventType
{
    Restaurant_Table_Equals_Household,
    Restaurant_Table_Equals_Half_Household,
    WorkMeeting_Many_Meetings,
    WorkMeeting_Few_Meetings
};

struct EventSimulationConfig {
    EventType type;
    std::string panvadere_file;
    double infection_parameter_k = 1.0;
    bool use_panvadere_init      = true;
    int event_duration_hours     = 2;
    int max_attendees            = 50;
};

class EventSimulator
{
public:
    // Step 2: Berechne K-Parameter für Event
    static mio::IOResult<double> calculate_infection_parameter_k(const EventSimulationConfig& config);

    // Step 3.1: Initialisierung von Panvadere
    static mio::IOResult<std::map<uint32_t, bool>> initialize_from_panvadere(const std::string& panvadere_file,
                                                                             EventType event_type);

    // Step 3.2: Initialisierung durch eigene Event-Simulation
    static mio::IOResult<std::map<uint32_t, bool>> initialize_from_event_simulation(const EventSimulationConfig& config,
                                                                                    const mio::abm::World& city);

    static std::string event_type_to_string(EventType type);

private:
    static std::map<uint32_t, bool> map_restaurant_tables_to_households(const std::map<uint32_t, bool>& panvadere_data);
    static std::map<uint32_t, bool> map_work_meeting_to_households(const std::map<uint32_t, bool>& panvadere_data);
    static std::map<uint32_t, bool> map_choir_to_households(const std::map<uint32_t, bool>& panvadere_data);

    static mio::IOResult<mio::abm::World> create_event_world(const EventSimulationConfig& config);
    static std::map<uint32_t, bool> simulate_event_transmission(const mio::abm::World& event_world,
                                                                const EventSimulationConfig& config);
};