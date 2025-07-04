#pragma once
#include "memilio/io/result_io.h"
#include <vector>
#include <map>
#include <random>
#include "abm/world.h"
#include "abm/location.h"
#include "abm/person.h"
#include "city_parameters.h"

// Configuration for representative German city simulation
struct CityConfig {
    int total_population = 1000;

    // Infrastructure will be calculated based on German demographic data
    CityParameters::CityInfrastructure infrastructure() const
    {
        return CityParameters::CityInfrastructure::calculate(total_population);
    }
};

class CityBuilder
{
public:
    static mio::IOResult<mio::abm::World> build_world(const CityConfig& config);
    static void print_city_summary(const CityConfig& config);

private:
    static std::vector<mio::abm::LocationId> create_households(mio::abm::World& world, int num_households);
    static std::vector<mio::abm::LocationId> create_workplaces(mio::abm::World& world, int num_workplaces);
    static std::vector<mio::abm::LocationId> create_schools(mio::abm::World& world, int num_schools);
    static std::vector<mio::abm::LocationId> create_shops(mio::abm::World& world, int grocery_stores, int pharmacies,
                                                          int general_stores);
    static std::vector<mio::abm::LocationId> create_events(mio::abm::World& world, int large, int small);
    static mio::IOResult<void> assign_people_to_locations(
        mio::abm::World& world, const std::vector<mio::abm::LocationId>& households,
        const std::vector<mio::abm::LocationId>& workplaces, const std::vector<mio::abm::LocationId>& schools,
        const std::vector<mio::abm::LocationId>& shops, const std::vector<mio::abm::LocationId>& events,
        const std::vector<mio::abm::LocationId>& hospitals, const std::vector<mio::abm::LocationId>& icus,
        int total_population);
    static mio::AgeGroup assign_age_group_from_demographics(mio::RandomNumberGenerator& gen);
    static mio::abm::LocationId assign_household(const std::vector<mio::abm::LocationId>& households,
                                                 mio::RandomNumberGenerator& gen);
    static bool should_attend_school(const mio::AgeGroup& age_group, mio::RandomNumberGenerator& gen);
    static bool should_be_employed(const mio::AgeGroup& age_group, mio::RandomNumberGenerator& gen);
};