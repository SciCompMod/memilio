#ifndef LOCATION_SPLIT_CITY_BUILDER_H
#define LOCATION_SPLIT_CITY_BUILDER_H

#include "city_parameters.h"
#include "defaults.h"

#include "abm/location_id.h"
#include "abm/model.h"
#include "memilio/utils/random_number_generator.h"

#include <string>
#include <vector>

/**
 * @brief Configuration of the synthetic city that the simulation is run on.
 */
struct CityConfig {
    int total_population          = Config::DEFAULT_POPULATION;
    CityParameters::Region region = CityParameters::Region::Germany;

    /// @brief Demographic parameters of the configured region.
    CityParameters::RegionParameters region_parameters() const
    {
        return CityParameters::get_region_parameters(region);
    }

    /// @brief Infrastructure derived from the population size and the region.
    CityParameters::CityInfrastructure infrastructure() const
    {
        return CityParameters::CityInfrastructure::calculate(total_population, region_parameters());
    }
};

/**
 * @brief Builds a synthetic city as an abm::Model from demographic and infrastructure parameters.
 */
class CityBuilder
{
public:
    /**
     * @brief Build the Model of a city.
     * @param[in] config Configuration of the city.
     * @param[in] rng Random number generator used for the Model and for drawing the households.
     */
    static mio::abm::Model build_model(const CityConfig& config, const mio::RandomNumberGenerator& rng);

    /// @brief Print a human readable summary of the city to stdout.
    static void print_city_summary(const CityConfig& config);

    /// @brief Write the derived city statistics as a csv-like key/value file.
    static void save_city_to_file(const CityConfig& config, const std::string& filename);

private:
    static std::vector<mio::abm::LocationId> create_locations(mio::abm::Model& model, mio::abm::LocationType type,
                                                              int count);
    static void create_and_assign_people(mio::abm::Model& model, const CityConfig& config,
                                         const CityParameters::CityInfrastructure& infra,
                                         const std::vector<mio::abm::LocationId>& households,
                                         const std::vector<mio::abm::LocationId>& workplaces,
                                         const std::vector<mio::abm::LocationId>& prim_schools,
                                         const std::vector<mio::abm::LocationId>& sec_schools,
                                         const std::vector<mio::abm::LocationId>& shops,
                                         const std::vector<mio::abm::LocationId>& events,
                                         const std::vector<mio::abm::LocationId>& hospitals,
                                         const std::vector<mio::abm::LocationId>& icus);
};

#endif // LOCATION_SPLIT_CITY_BUILDER_H
