#pragma once
#include "memilio/io/result_io.h"
#include <vector>
#include <map>
#include <random>
#include "abm/world.h"
#include "abm/location.h"
#include "abm/person.h"
#include "city_parameters.h"
#include "config_validation.h"

// Configuration for representative German city simulation
struct CityConfig {
    int total_population = Config::DEFAULT_POPULATION;

    // Infrastructure will be calculated based on German demographic data
    CityParameters::CityInfrastructure infrastructure() const
    {
        return CityParameters::CityInfrastructure::calculate(total_population);
    }
};

class CityBuilder
{
public:
    static mio::IOResult<mio::abm::World> build_world(const CityConfig& config, const mio::RandomNumberGenerator& rng);
    static void print_city_summary(const CityConfig& config);

private:
    static std::vector<mio::abm::LocationId> create_households(mio::abm::World& world, int num_households);
    static std::vector<mio::abm::LocationId> create_workplaces(mio::abm::World& world, int num_workplaces);
    static std::vector<mio::abm::LocationId> create_schools(mio::abm::World& world, int num_schools);
    static std::vector<mio::abm::LocationId> create_shops(mio::abm::World& world, int num_shops);
    static std::vector<mio::abm::LocationId> create_events(mio::abm::World& world, int num_events);
    static std::vector<int> create_age_vector(int total_population);
    static mio::IOResult<void> create_and_assign_people_to_locations(
        mio::abm::World& world, const std::vector<mio::abm::LocationId>& households,
        const std::vector<mio::abm::LocationId>& workplaces, const std::vector<mio::abm::LocationId>& prim_schools,
        const std::vector<mio::abm::LocationId>& sec_schools, const std::vector<mio::abm::LocationId>& shops,
        const std::vector<mio::abm::LocationId>& events, const std::vector<mio::abm::LocationId>& hospitals,
        const std::vector<mio::abm::LocationId>& icus, int total_population, std::vector<int>& hh_per_size,
        const CityParameters::CityInfrastructure config);
    static int ageGroupTInt6(mio::AgeGroup age_group);
};