#ifndef LOCATION_SPLIT_CITY_PARAMETERS_H
#define LOCATION_SPLIT_CITY_PARAMETERS_H

#include <string>
#include <tuple>
#include <utility>
#include <vector>

/**
 * @file city_parameters.h
 * @brief Demographic and infrastructure parameters used to build a representative city.
 *
 * The parameters of a region are bundled in RegionParameters so that the same CityBuilder can
 * build a German, French or US city. On the panvXabmSim branch these lived in three headers that
 * all reused the namespace CityParameters and therefore could not be used together.
 *
 * Sources per region are documented at the corresponding factory function in city_parameters.cpp.
 */
namespace CityParameters
{

/// @brief Regions for which demographic parameters and contact matrices are available.
enum class Region
{
    Germany,
    France,
    UnitedStates
};

/// @brief Convert a Region to the string accepted by the --region command line option.
std::string region_to_string(Region region);

/// @brief Parse the string accepted by the --region command line option. Returns false if unknown.
bool region_from_string(const std::string& name, Region& region);

/**
 * @brief Infrastructure ratios of a region.
 *
 * All rates refer to the total population unless stated otherwise.
 */
struct InfrastructureRatios {
    double employment_rate; ///< Share of the working age population that is employed.
    double people_per_workplace; ///< Average number of employees per workplace.
    double school_rate; ///< Share of the total population that attends school.
    double max_students_per_elementary_school;
    double max_students_per_secondary_school;
    double ratio_elementary_to_secondary_school; ///< Secondary students per elementary student.
    double retail_stores_per_1000_people;
    double people_per_event; ///< Average number of people per social event location.
};

/**
 * @brief All demographic parameters of one region.
 */
struct RegionParameters {
    Region region;
    /// Shares of the age groups 0-4, 5-14, 15-34, 35-59, 60-79, 80+. Six entries.
    std::vector<double> age_distribution;
    /// Shares of the household sizes 1, 2, 3, 4, 5+. Five entries.
    std::vector<double> household_size_distribution;
    InfrastructureRatios ratios;
};

/// @brief Parameters of a representative German city.
RegionParameters germany();
/// @brief Parameters of a representative French city.
RegionParameters france();
/// @brief Parameters of a representative US city.
RegionParameters united_states();
/// @brief Parameters of the given region.
RegionParameters get_region_parameters(Region region);

/**
 * @brief Infrastructure of a city, derived from its population size and its RegionParameters.
 */
struct CityInfrastructure {
    std::vector<int> num_households_hh_size; ///< Number of households per household size 1..5.
    int num_workplaces;
    int num_elementary_schools;
    int num_secondary_schools;
    int num_hospitals;
    int num_icus;
    int num_stores;
    int num_events;
    int num_persons_elementary_schools;
    int num_persons_secondary_schools;
    int num_worker;

    /**
     * @brief Derive the infrastructure of a city.
     * @param[in] population Total population of the city.
     * @param[in] params Demographic parameters of the region.
     */
    static CityInfrastructure calculate(int population, const RegionParameters& params);
};

/// @brief Number of persons per age group for the given population size.
std::vector<int> age_vector(int population, const RegionParameters& params);

/// @brief Number of households per household size 1..5 that add up to the given population.
std::vector<int> household_sizes(int population, const RegionParameters& params);

} // namespace CityParameters

#endif // LOCATION_SPLIT_CITY_PARAMETERS_H
