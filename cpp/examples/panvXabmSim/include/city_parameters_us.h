#pragma once
#include <vector>
#include <map>

/**
 * @file city_parameters_us.h
 * @brief Parameters for building a representative US city
 * 
 * This file contains demographic and infrastructure parameters based on US statistics
 * to create realistic city simulations.
 * 
 * Sources:
 * - U.S. Census Bureau: Annual Estimates of the Resident Population 2024
 *   https://www.census.gov/data/tables/time-series/demo/popest/2020s-national-detail.html
 * - U.S. Census Bureau: American Community Survey (ACS) 2023
 *   https://www.census.gov/programs-surveys/acs
 * - Bureau of Labor Statistics: Employment Situation 2024
 *   https://www.bls.gov/news.release/empsit.toc.htm
 * - National Center for Education Statistics (NCES): Public Elementary and Secondary School Statistics 2023
 *   https://nces.ed.gov/programs/digest/
 * - U.S. Census Bureau: Annual Retail Trade Survey 2022
 *   https://www.census.gov/retail/index.html
 */

namespace CityParameters
{

/**
 * @brief US age distribution based on 2024 data
 * Age groups: 0-4, 5-14, 15-34, 35-59, 60-79, 80+
 */
const std::vector<double> US_AGE_DISTRIBUTION = {
    0.061, // 0-4 years: 6.1%
    0.128, // 5-14 years: 12.8%
    0.267, // 15-34 years: 26.7%
    0.317, // 35-59 years: 31.7%
    0.181, // 60-79 years: 18.1%
    0.046 // 80+ years: 4.6%
};

/**
 * @brief US household size distribution
 * Source: U.S. Census Bureau - American Community Survey (ACS) 2023
 */
const std::vector<double> HOUSEHOLD_SIZE_DISTRIBUTION = {
    0.283, // 1-person households: 28.3%
    0.349, // 2-person households: 34.9%
    0.150, // 3-person households: 15.0%
    0.127, // 4-person households: 12.7%
    0.091 // 5+ person households: 9.1%
};

/**
 * @brief Infrastructure ratios per 1000 people
 * These ratios are based on US national averages
 */
struct InfrastructureRatios {
    // Employment and workplaces
    // Source: Bureau of Labor Statistics - Employment Situation 2024
    static constexpr double EMPLOYMENT_RATE      = 0.738; // 73.8% for ages 16-64
    static constexpr double PEOPLE_PER_WORKPLACE = 12.5; // Average employees per workplace (larger than Europe)

    // Education
    // Source: National Center for Education Statistics (NCES) 2023
    static constexpr double SCHOOL_RATE                          = 0.155; // 15.5% of population (ages 5-17)
    static constexpr double MAX_STUDENTS_PER_ELEMENTARY_SCHOOL   = 470; // Primary schools (K-5)
    static constexpr double MAX_STUDENTS_PER_SECONDARY_SCHOOL    = 850; // Secondary schools (middle + high schools)
    static constexpr double RATIO_ELEMENTARY_TO_SECONDARY_SCHOOL = 1.65; // Ratio of elementary to secondary students

    // Retail and services
    // Source: U.S. Census Bureau - Annual Retail Trade Survey 2022
    static constexpr double amount_of_retail_stores_per_1000_people = 3.5;

    // Social and event locations
    static constexpr double PEOPLE_PER_EVENT = 15.0;
};

/**
 * @brief Calculate city infrastructure based on population size
 * @param population Total population of the city
 * @return Infrastructure configuration for the city
 */
struct CityInfrastructure {
    std::vector<int> num_households_hh_size;
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

    std::vector<int> calc_household_sizes(int population) const
    {
        const int n_bins = HOUSEHOLD_SIZE_DISTRIBUTION.size();

        // 1. Estimate household count based on average household size
        double avg_household_size = 0.0;
        for (int i = 0; i < n_bins; ++i)
            avg_household_size += (i + 1) * HOUSEHOLD_SIZE_DISTRIBUTION[i];

        int estimated_total_households = static_cast<int>(std::round(population / avg_household_size));

        // 2. Calculate raw household numbers per bin
        std::vector<double> raw_households(n_bins);
        for (int i = 0; i < n_bins; ++i)
            raw_households[i] = estimated_total_households * HOUSEHOLD_SIZE_DISTRIBUTION[i];

        // 3. Round and fix totals
        std::vector<int> households_by_size(n_bins);
        for (int i = 0; i < n_bins; ++i)
            households_by_size[i] = static_cast<int>(std::round(raw_households[i]));

        // 4. Ensure population total is correct
        int total_people = 0;
        for (int i = 0; i < n_bins; ++i)
            total_people += households_by_size[i] * (i + 1);

        int diff = population - total_people;
        // Fix small rounding errors by adding/subtracting people in 5-person households
        while (diff != 0) {
            int idx = (diff > 0) ? 4 : 0; // 5-person or 1-person
            households_by_size[idx] += (diff > 0) ? 1 : -1;
            diff += (diff > 0) ? -(idx + 1) : (idx + 1);
        }

        return households_by_size;
    }

    std::pair<int, int> calc_num_workplaces_and_worker(int population) const
    {
        // Calculate number of workplaces based on employment rate and average employees per workplace
        std::vector<int> age_vector(US_AGE_DISTRIBUTION.size());
        for (size_t i = 0; i < CityParameters::US_AGE_DISTRIBUTION.size(); ++i) {
            int age_group_population = static_cast<int>(population * CityParameters::US_AGE_DISTRIBUTION[i]);
            age_vector.at(i)         = age_group_population;
        }

        int n_potential_worker = age_vector[2] + age_vector[3] + ((4 / 20) * age_vector[4]);
        auto n_w = static_cast<int>(std::round(n_potential_worker * InfrastructureRatios::EMPLOYMENT_RATE));
        int n_wp = static_cast<int>(std::round(n_w / InfrastructureRatios::PEOPLE_PER_WORKPLACE));
        return std::make_pair(n_wp, n_w);
    }

    std::tuple<int, int, int, int> calc_num_elem_and_sec_schools(int population) const
    {
        // Calculate number of elementary schools based on school rate and max students per school
        int total_students          = static_cast<int>(population * InfrastructureRatios::SCHOOL_RATE);
        int num_elementary_students = static_cast<int>(
            std::round(total_students / (1.0 + InfrastructureRatios::RATIO_ELEMENTARY_TO_SECONDARY_SCHOOL)));
        int num_secondary_students = total_students - num_elementary_students;
        return std::make_tuple(
            static_cast<int>(static_cast<int>(
                std::ceil(num_elementary_students / InfrastructureRatios::MAX_STUDENTS_PER_ELEMENTARY_SCHOOL))),
            static_cast<int>(static_cast<int>(
                std::ceil(num_secondary_students / InfrastructureRatios::MAX_STUDENTS_PER_SECONDARY_SCHOOL))),
            num_elementary_students, num_secondary_students);
    }

    int calc_stores(int population) const
    {
        // Calculate number of retail stores based on average stores per 1000 people
        int n_stores = static_cast<int>(
            std::ceil(population * InfrastructureRatios::amount_of_retail_stores_per_1000_people / 1000.0));
        return n_stores;
    }

    int calc_num_events(int population) const
    {
        // Calculate number of events based on people per event and chance to attend
        int n_events = static_cast<int>(std::ceil(population / InfrastructureRatios::PEOPLE_PER_EVENT));
        return n_events;
    }

    static CityInfrastructure calculate(int population)
    {
        CityInfrastructure infra;

        // Calculate household sizes
        infra.num_households_hh_size = infra.calc_household_sizes(population);

        // Calculate workplaces
        auto [num_workplaces, num_worker] = infra.calc_num_workplaces_and_worker(population);
        infra.num_workplaces              = num_workplaces;
        infra.num_worker                  = num_worker;

        // Calculate schools
        auto [num_elem_schools, num_sec_schools, num_elem_students, num_sec_students] =
            infra.calc_num_elem_and_sec_schools(population);
        infra.num_persons_elementary_schools = num_elem_students;
        infra.num_persons_secondary_schools  = num_sec_students;
        infra.num_elementary_schools         = num_elem_schools;
        infra.num_secondary_schools          = num_sec_schools;

        // Calculate hospitals and ICUs
        infra.num_hospitals = 1; // 1 hospital per 100,000 people
        infra.num_icus      = 1;

        // Calculate stores
        infra.num_stores = infra.calc_stores(population);

        // Calculate events
        infra.num_events = infra.calc_num_events(population);

        return infra;
    }
};

} // namespace CityParameters
