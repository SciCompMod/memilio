/**
 * @file german_city_example.cpp
 * @brief Example demonstrating how to create a representative German city simulation
 * 
 * This example shows how to use the enhanced city builder with German demographic data
 * to create realistic city simulations at different scales.
 */

#include "../include/city_builder.h"
#include "../include/city_parameters.h"
#include <iostream>
#include <iomanip>

void print_city_stats(const CityConfig& config)
{
    auto infra = config.infrastructure();

    std::cout << "\n=== German City Simulation Statistics ===\n";
    std::cout << "Population: " << config.total_population << "\n\n";

    std::cout << "Households: " << infra.num_households << " (avg " << std::fixed << std::setprecision(1)
              << static_cast<double>(config.total_population) / infra.num_households << " people/household)\n";

    std::cout << "Workplaces: " << infra.num_workplaces << "\n";
    std::cout << "Schools: " << infra.num_elementary_schools << " elementary + " << infra.num_secondary_schools
              << " secondary\n";

    std::cout << "Healthcare: " << infra.num_hospitals << " hospitals, " << infra.num_icus << " ICUs\n";

    std::cout << "Retail: " << infra.num_grocery_stores << " grocery stores, " << infra.num_pharmacies
              << " pharmacies, " << infra.num_general_stores << " general stores\n";

    std::cout << "Social venues: " << infra.num_restaurants << " restaurants, " << infra.num_bars << " bars\n";

    std::cout << "Events: " << infra.num_large_events << " large events, " << infra.num_small_events
              << " small events\n";

    std::cout << "========================================\n\n";
}

int main()
{
    // Example 1: Small town (10,000 people)
    std::cout << "Creating small German town...\n";
    CityConfig small_town_config;
    small_town_config.total_population = 10000;
    print_city_stats(small_town_config);

    auto small_town_result = CityBuilder::build_world(small_town_config);
    if (small_town_result) {
        std::cout << "Small town created successfully!\n";
    }
    else {
        std::cout << "Failed to create small town: " << small_town_result.error().formatted_message() << "\n";
    }

    // Example 2: Medium city (100,000 people)
    std::cout << "\nCreating medium German city...\n";
    CityConfig medium_city_config;
    medium_city_config.total_population = 100000;
    print_city_stats(medium_city_config);

    auto medium_city_result = CityBuilder::build_world(medium_city_config);
    if (medium_city_result) {
        std::cout << "Medium city created successfully!\n";
    }
    else {
        std::cout << "Failed to create medium city: " << medium_city_result.error().formatted_message() << "\n";
    }

    // Example 3: Large city (1,000,000 people)
    std::cout << "\nCreating large German city...\n";
    CityConfig large_city_config;
    large_city_config.total_population = 1000000;
    print_city_stats(large_city_config);

    auto large_city_result = CityBuilder::build_world(large_city_config);
    if (large_city_result) {
        std::cout << "Large city created successfully!\n";
    }
    else {
        std::cout << "Failed to create large city: " << large_city_result.error().formatted_message() << "\n";
    }

    // Example 4: Scale comparable to Germany (80,000,000 people)
    std::cout << "\nCreating Germany-scale simulation...\n";
    CityConfig germany_config;
    germany_config.total_population = 80000000;
    print_city_stats(germany_config);

    std::cout << "Note: For 80 million people, you would typically want to use\n";
    std::cout << "a distributed simulation approach rather than a single city model.\n";

    return 0;
}
