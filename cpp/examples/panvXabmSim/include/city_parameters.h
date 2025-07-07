#pragma once
#include <vector>
#include <map>

/**
 * @file city_parameters.h
 * @brief Parameters for building a representative German city
 * 
 * This file contains demographic and infrastructure parameters based on German statistics
 * to create realistic city simulations.
 * 
 * Sources:
 * - German Federal Statistical Office (Destatis) 2023: Population by age groups
 * - German Federal Statistical Office: Household statistics 2023
 * - German Federal Statistical Office: Labor force statistics 2023
 * - German Federal Statistical Office: Education statistics 2023
 * - German Trade Association (HDE): Retail structure statistics 2023
 * - German Hotel and Restaurant Association (DEHOGA): Hospitality statistics 2023
 * - German Federal Statistical Office: Healthcare infrastructure 2023
 */

namespace CityParameters
{

/**
 * @brief German age distribution based on 2023 census data
 * Age groups: 0-4, 5-14, 15-34, 35-59, 60-79, 80+
 * Source: Destatis 2023 population statistics
 */
const std::vector<double> GERMAN_AGE_DISTRIBUTION = {
    0.047, // 0-4 years: 4.7%
    0.096, // 5-14 years: 9.6%
    0.245, // 15-34 years: 24.5%
    0.338, // 35-59 years: 33.8%
    0.214, // 60-79 years: 21.4%
    0.060 // 80+ years: 6.0%
};

/**
 * @brief German household size distribution
 * Source: Destatis Mikrozensus 2023
 */
const std::vector<double> HOUSEHOLD_SIZE_DISTRIBUTION = {
    0.432, // 1-person households: 43.2%
    0.336, // 2-person households: 33.6%
    0.117, // 3-person households: 11.7%
    0.095, // 4-person households: 9.5%
    0.020 // 5+ person households: 2.0%
};

/**
 * @brief Average household size in Germany
 * Source: Destatis 2023
 */
const double AVERAGE_HOUSEHOLD_SIZE = 1.95;

/**
 * @brief Infrastructure ratios per 1000 people
 * These ratios are based on German national averages
 */
struct InfrastructureRatios {
    // Employment and workplaces
    // Source: Destatis labor force statistics 2023
    static constexpr double EMPLOYMENT_RATE      = 0.638; // 63.8% of population aged 15-64
    static constexpr double PEOPLE_PER_WORKPLACE = 25.0; // Average employees per workplace

    // Education
    // Source: Destatis education statistics 2023
    static constexpr double STUDENTS_PER_ELEMENTARY_SCHOOL = 180.0; // Primary schools
    static constexpr double STUDENTS_PER_SECONDARY_SCHOOL  = 450.0; // Secondary schools
    static constexpr double ELEMENTARY_TO_SECONDARY_RATIO  = 3.5; // 3.5 elementary per secondary

    // Retail and services
    // Source: HDE retail statistics 2023
    static constexpr double PEOPLE_PER_GROCERY_STORE = 2000.0; // Basic necessities
    static constexpr double PEOPLE_PER_PHARMACY      = 3500.0; // Pharmacies
    static constexpr double PEOPLE_PER_GENERAL_STORE = 1500.0; // General retail

    // Social and event locations
    // Source: DEHOGA hospitality statistics 2023
    static constexpr double PEOPLE_PER_RESTAURANT  = 400.0; // Restaurants/cafes
    static constexpr double PEOPLE_PER_BAR         = 800.0; // Bars/pubs
    static constexpr double PEOPLE_PER_LARGE_EVENT = 50000.0; // Concert halls, stadiums
    static constexpr double PEOPLE_PER_SMALL_EVENT = 2000.0; // Community centers, clubs
};

/**
 * @brief School attendance rates by age group
 * Source: Destatis education statistics 2023
 */
const std::map<int, double> SCHOOL_ATTENDANCE_RATES = {
    {0, 0.0}, // 0-4 years: 0% (some kindergarten, but not modeled)
    {1, 1.0}, // 5-14 years: 100% school attendance
    {2, 0.0}, // 15-34 years: 0% in education/training
    {3, 0.0}, // 35-59 years: 0% in continuing education
    {4, 0.0}, // 60-79 years: 0%
    {5, 0.0} // 80+ years: 0%
};

/**
 * @brief Employment rates by age group
 * Source: Destatis labor force statistics 2023
 */
const std::map<int, double> EMPLOYMENT_RATES = {
    {0, 0.0}, // 0-4 years: 0%
    {1, 0.0}, // 5-14 years: 0%
    {2, 1.0}, // 15-34 years: 78%
    {3, 1.0}, // 35-59 years: 85%
    {4, 0.0}, // 60-79 years: 32% (part-time, retirement transition)
    {5, 0.0} // 80+ years: 0%
};

/**
 * @brief Calculate city infrastructure based on population size
 * @param population Total population of the city
 * @return Infrastructure configuration for the city
 */
struct CityInfrastructure {
    int num_households;
    int num_workplaces;
    int num_elementary_schools;
    int num_secondary_schools;
    int num_hospitals;
    int num_icus;
    int num_grocery_stores;
    int num_pharmacies;
    int num_general_stores;
    int num_restaurants;
    int num_bars;
    int num_large_events;
    int num_small_events;

    static CityInfrastructure calculate(int population)
    {
        CityInfrastructure infra;

        // Households based on average household size
        infra.num_households = static_cast<int>(population / AVERAGE_HOUSEHOLD_SIZE);

        // Workplaces
        int working_population = static_cast<int>(population * InfrastructureRatios::EMPLOYMENT_RATE);
        infra.num_workplaces =
            std::max(1, static_cast<int>(working_population / InfrastructureRatios::PEOPLE_PER_WORKPLACE));

        // Schools
        int school_age_population =
            static_cast<int>(population * (GERMAN_AGE_DISTRIBUTION[1] + GERMAN_AGE_DISTRIBUTION[2] * 0.45));
        infra.num_elementary_schools = std::max(
            1, static_cast<int>(school_age_population * 0.7 / InfrastructureRatios::STUDENTS_PER_ELEMENTARY_SCHOOL));
        infra.num_secondary_schools = std::max(
            1, static_cast<int>(infra.num_elementary_schools / InfrastructureRatios::ELEMENTARY_TO_SECONDARY_RATIO));

        // Healthcare
        // infra.num_hospitals = std::max(1, static_cast<int>(population / InfrastructureRatios::PEOPLE_PER_HOSPITAL_BED));
        // infra.num_icus      = std::max(1, static_cast<int>(population / InfrastructureRatios::PEOPLE_PER_ICU_BED));
        infra.num_hospitals = 1;
        infra.num_icus      = 1; // Simplified for this model, could be adjusted based on real ratios
        // Retail
        infra.num_grocery_stores =
            std::max(1, static_cast<int>(population / InfrastructureRatios::PEOPLE_PER_GROCERY_STORE));
        infra.num_pharmacies = std::max(1, static_cast<int>(population / InfrastructureRatios::PEOPLE_PER_PHARMACY));
        infra.num_general_stores =
            std::max(1, static_cast<int>(population / InfrastructureRatios::PEOPLE_PER_GENERAL_STORE));

        // Social venues
        infra.num_restaurants = std::max(1, static_cast<int>(population / InfrastructureRatios::PEOPLE_PER_RESTAURANT));
        infra.num_bars        = std::max(1, static_cast<int>(population / InfrastructureRatios::PEOPLE_PER_BAR));
        infra.num_large_events =
            std::max(1, static_cast<int>(population / InfrastructureRatios::PEOPLE_PER_LARGE_EVENT));
        infra.num_small_events =
            std::max(1, static_cast<int>(population / InfrastructureRatios::PEOPLE_PER_SMALL_EVENT));

        return infra;
    }
};

} // namespace CityParameters
