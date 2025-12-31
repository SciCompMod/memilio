#pragma once
#include <vector>
#include <map>

/**
 * @file city_parameters_france.h
 * @brief Parameters for building a representative French city
 * 
 * This file contains demographic and infrastructure parameters based on French statistics
 * to create realistic city simulations.
 * 
 * Sources:
 * - INSEE (Institut national de la statistique et des études économiques) 2024: Population par âge
 *   https://www.insee.fr/en/statistiques/6040016
 * - INSEE: Composition des ménages 2020
 *   https://www.insee.fr/fr/statistiques/2381486
 * - Ministère de l'Éducation nationale: Repères et références statistiques 2023
 * - FCD (Fédération du Commerce et de la Distribution) 2023
 * - Employment rate https://www.insee.fr/en/outil-interactif/5543645/tableau/50_MTS/51_EPA
 */

namespace CityParameters
{

/**
 * @brief French age distribution based on 2024 data
 * Age groups: 0-4, 5-14, 15-34, 35-59, 60-79, 80+
 */
const std::vector<double> AGE_DISTRIBUTION = {
    0.052, // 0-4 years: 5.2%
    0.122, // 5-14 years: 12.2%
    0.235, // 15-34 years: 23.4%
    0.319, // 35-59 years: 31.9%
    0.211, // 60-79 years: 21.0%
    0.061 // 80+ years: 6.1%
};

/**
 * @brief French household size distribution
 * Source: INSEE - Composition des ménages 2020
 */
const std::vector<double> HOUSEHOLD_SIZE_DISTRIBUTION = {
    0.384, // 1-person households: 38.4%
    0.324, // 2-person households: 32.4%
    0.130, // 3-person households: 13.0%
    0.108, // 4-person households: 10.8%
    0.054 // 5+ person households: 5.4%
};

/**
 * @brief Infrastructure ratios per 1000 people
 * These ratios are based on French national averages
 */
struct InfrastructureRatios {
    // Employment and workplaces
    // Source: INSEE - Taux d'emploi 2023
    static constexpr double EMPLOYMENT_RATE =
        0.744; // 74.4% for ages 15-64 https://www.insee.fr/en/outil-interactif/5543645/tableau/50_MTS/51_EPA
    static constexpr double PEOPLE_PER_WORKPLACE = 10.0; // Average employees per workplace

    // Education
    // https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.education.gouv.fr/media/159309/download&ved=2ahUKEwirla_lneaRAxXZR_EDHT4kFdIQFnoECBoQAQ&usg=AOvVaw3qgtpLr0ugXgQG_5NtfTWf
    static constexpr double SCHOOL_RATE                        = 0.151; //15.1% of population
    static constexpr double MAX_STUDENTS_PER_ELEMENTARY_SCHOOL = 200; // Primary schools (écoles élémentaires)
    static constexpr double MAX_STUDENTS_PER_SECONDARY_SCHOOL = 300; // Secondary schools (collèges: ~500, lycées: ~540)
    static constexpr double RATIO_ELEMENTARY_TO_SECONDARY_SCHOOL = 1.65; // Ratio of elementary to secondary students

    // Retail and services
    // Source: FCD (Fédération du Commerce et de la Distribution) 2023
    static constexpr double amount_of_retail_stores_per_1000_people = 4.2;

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
        std::vector<int> age_vector(AGE_DISTRIBUTION.size());
        for (size_t i = 0; i < CityParameters::AGE_DISTRIBUTION.size(); ++i) {
            int age_group_population = static_cast<int>(population * CityParameters::AGE_DISTRIBUTION[i]);
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
