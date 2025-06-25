#pragma once
#include "abm/abm.h"
#include "memilio/io/result_io.h"
#include <vector>
#include <map>

struct CityConfig {
    int total_population = 1000;
    int num_households = 250;    // ~4 people per household
    int num_workplaces = 50;     // ~20 people per workplace
    int num_schools = 3;         // Elementary, Middle, High
    int num_large_events = 2;    // Concert hall, stadium
    int num_small_events = 5;    // Restaurants, bars
};

class CityBuilder {
public:
    static mio::IOResult<mio::abm::World> create_city_from_census(const CityConfig& config);
    
private:
    static std::vector<mio::abm::LocationId> create_households(mio::abm::World& world, int num_households);
    static std::vector<mio::abm::LocationId> create_workplaces(mio::abm::World& world, int num_workplaces);
    static std::vector<mio::abm::LocationId> create_schools(mio::abm::World& world, int num_schools);
    static std::vector<mio::abm::LocationId> create_events(mio::abm::World& world, int large, int small);
    static mio::IOResult<void> assign_people_to_locations(mio::abm::World& world, 
                                                         const std::vector<mio::abm::LocationId>& households,
                                                         const std::vector<mio::abm::LocationId>& workplaces,
                                                         const std::vector<mio::abm::LocationId>& schools,
                                                         int total_population);
    static mio::AgeGroup assign_age_group_from_demographics(int person_index);
};