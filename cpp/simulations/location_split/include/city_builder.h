#ifndef LOCATION_SPLIT_CITY_BUILDER_H
#define LOCATION_SPLIT_CITY_BUILDER_H

#include "city_parameters.h"
#include "defaults.h"
#include "parameter_setter.h"

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
    /// Pupils per school room. 0 keeps every school as a single Location. Otherwise every school is
    /// split into ceil(pupils / school_room_size) rooms of near equal size, each a School Location of
    /// its own, and every pupil is assigned to one room at random for the whole simulation.
    int school_room_size = 0;
    /// Relative transmission weight per location type, see ContactWeights.
    ContactWeights contact_weights{};

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
    /**
     * @brief Split every School into rooms of at most room_size pupils.
     * The pupils of a school are shuffled with a stream that is independent of the Model's generator,
     * so the city and the initial infections drawn after the build are the same as without the split.
     * The first room reuses the school's Location, the others are new School Location%s.
     */
    static void split_schools_into_rooms(mio::abm::Model& model, int room_size,
                                         const mio::RandomNumberGenerator& rng);
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
