#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>
#include <string>
#include <map>
#include "abm/abm.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/uncertain_value.h"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "abm/vaccine.h"
#include "abm/common_abm_loggers.h"
#include "memilio/utils/miompi.h"
#include "memilio/io/binary_serializer.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "abm/household.h"
#include "memilio/utils/random_number_generator.h"

// Number of age groups and their definitions
size_t num_age_groups        = 6;
const auto age_group_0_to_4   = mio::AgeGroup(0);
const auto age_group_5_to_14  = mio::AgeGroup(1);
const auto age_group_15_to_34 = mio::AgeGroup(2);
const auto age_group_35_to_59 = mio::AgeGroup(3);
const auto age_group_60_to_79 = mio::AgeGroup(4);
const auto age_group_80_plus  = mio::AgeGroup(5);

/**
 * Helper function to convert log-normal parameters from mean and std to mu and sigma
 */
std::pair<double, double> get_my_and_sigma(std::pair<double, double> mean_and_std)
{
    auto mean    = mean_and_std.first;
    auto stddev  = mean_and_std.second;
    double my    = log(mean * mean / sqrt(mean * mean + stddev * stddev));
    double sigma = sqrt(log(1 + stddev * stddev / (mean * mean)));
    return {my, sigma};
}

/**
 * Sets the parameters for the disease model
 */
void set_parameters(mio::abm::Parameters& params)
{
    auto incubation_period_my_sigma          = get_my_and_sigma({4.5, 1.5});
    params.get<mio::abm::IncubationPeriod>() = {incubation_period_my_sigma.first, incubation_period_my_sigma.second};

    auto InfectedNoSymptoms_to_symptoms_my_sigma             = get_my_and_sigma({1.1, 0.9});
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>() = {InfectedNoSymptoms_to_symptoms_my_sigma.first,
                                                                InfectedNoSymptoms_to_symptoms_my_sigma.second};

    auto TimeInfectedNoSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() = {TimeInfectedNoSymptomsToRecovered_my_sigma.first,
                                                                 TimeInfectedNoSymptomsToRecovered_my_sigma.second};

    auto TimeInfectedSymptomsToSevere_my_sigma           = get_my_and_sigma({6.6, 4.9});
    params.get<mio::abm::TimeInfectedSymptomsToSevere>() = {TimeInfectedSymptomsToSevere_my_sigma.first,
                                                            TimeInfectedSymptomsToSevere_my_sigma.second};

    auto TimeInfectedSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedSymptomsToRecovered>() = {TimeInfectedSymptomsToRecovered_my_sigma.first,
                                                               TimeInfectedSymptomsToRecovered_my_sigma.second};

    auto TimeInfectedSevereToCritical_my_sigma           = get_my_and_sigma({1.5, 2.0});
    params.get<mio::abm::TimeInfectedSevereToCritical>() = {TimeInfectedSevereToCritical_my_sigma.first,
                                                            TimeInfectedSevereToCritical_my_sigma.second};

    auto TimeInfectedSevereToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedSevereToRecovered>() = {TimeInfectedSevereToRecovered_my_sigma.first,
                                                             TimeInfectedSevereToRecovered_my_sigma.second};

    auto TimeInfectedCriticalToDead_my_sigma           = get_my_and_sigma({10.7, 4.8});
    params.get<mio::abm::TimeInfectedCriticalToDead>() = {TimeInfectedCriticalToDead_my_sigma.first,
                                                          TimeInfectedCriticalToDead_my_sigma.second};

    auto TimeInfectedCriticalToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedCriticalToRecovered>() = {TimeInfectedCriticalToRecovered_my_sigma.first,
                                                               TimeInfectedCriticalToRecovered_my_sigma.second};

    // Set percentage parameters
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]  = 0.75;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = 0.75;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.8;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.8;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.8;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.8;

    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.01;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]  = 0.01;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.02;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.07;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.3;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.6;

    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]  = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.17;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.6;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.8;

    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.055;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]  = 0.055;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.14;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.28;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.55;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.75;

    // Set infection parameters
    params.get<mio::abm::InfectionRateFromViralShed>()[{mio::abm::VirusVariant::Wildtype}] = 5.6;
    
    // Set contact parameter (people of same age group meet more often)
    // params.get<mio::abm::AgeGroupGotoSocialEvent>() = true;
}


/**
 * Read infection data from file and return a map of infected person IDs
 */
std::map<uint32_t, bool> read_infection_data(const std::string& filename)
{
    std::map<uint32_t, bool> infected_status;
    std::map<std::string, std::vector<uint32_t>> tables;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return infected_status;
    }
    
    std::string line;
    // Skip header line
    std::getline(file, line);
    
    // Read each line
    while (std::getline(file, line)) {
        // Parse tab-separated values
        std::istringstream iss(line);
        uint32_t id;
        std::string table;
        int infected;
        
        if (iss >> id >> table >> infected) {
            infected_status[id] = (infected == 1);
            tables[table].push_back(id);
        }
    }
    
    file.close();
    return infected_status;
}

/**
 * Creates a restaurant simulation world from infection data file
 */
mio::abm::World create_restaurant_world(const std::string& infection_data_file)
{
    // Create world with 6 age groups
    auto world = mio::abm::World(num_age_groups);
    
    // Set infection parameters
    set_parameters(world.parameters);
    
    // Read infection data
    auto infected_status = read_infection_data(infection_data_file);
    
    // Create tables/households map
    std::map<std::string, std::vector<uint32_t>> tables;
    
    // Create one restaurant (social event location)
    auto restaurant = world.add_location(mio::abm::LocationType::SocialEvent);
    
    // Create a hospital and ICU for severe cases
    auto hospital = world.add_location(mio::abm::LocationType::Hospital);
    auto icu = world.add_location(mio::abm::LocationType::ICU);
    
    // Create homes (one per table)
    std::map<std::string, mio::abm::LocationId> table_homes;
    for (const auto& [id, infected] : infected_status) {
        // Create home location for each person (or group them by table)
        auto home = world.add_location(mio::abm::LocationType::Home);
        table_homes[std::to_string(id)] = home;
    }
    
    // Create people and assign them to locations
    std::vector<double> weights = {0.05, 0.05, 0.4, 0.4, 0.07, 0.03};
    for (const auto& [id, infected] : infected_status) {
        mio::AgeGroup age = mio::AgeGroup(
            mio::DiscreteDistribution<size_t>::get_instance()(
                world.get_rng(), weights
            )
        );
        
        // Create person and assign them to their home
        auto& person = world.add_person(table_homes[std::to_string(id)], age);
        
        // Assign the person to the restaurant, hospital and ICU
        person.set_assigned_location(restaurant);
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        
        // If infected according to data, set them as exposed
        if (infected) {
            auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
            person.add_new_infection(mio::abm::Infection(
                rng, 
                mio::abm::VirusVariant::Wildtype, 
                person.get_age(),
                world.parameters, 
                mio::abm::TimePoint(0), 
                mio::abm::InfectionState::InfectedNoSymptoms
            ));
        }
        
        // Keep track of people at each table for contact matrix
        tables[std::to_string(id)].push_back(person.get_person_id());
    }
    
    // Print an overview of the simulation
    std::cout << "Created restaurant simulation with " << world.get_persons().size() << " people." << std::endl;
    std::cout << "Infection data read from: " << infection_data_file << std::endl;
    std::cout << "Number of tables: " << tables.size() << std::endl;
    std::cout << "Infected people: " << infected_status.size() << std::endl;
    std::cout << "Tables and their members and their infection state:" << std::endl;
    for (const auto& [table, members] : tables) {
        std::cout << "  Table " << table << ": ";
        for (const auto& member : members) {
            std::cout << member << " ";
            auto person = world.get_person(member);
            if (person.is_infected(mio::abm::TimePoint(0))) {
                std::cout << "(" << static_cast<uint32_t>(person.get_infection_state(mio::abm::TimePoint(0))) << ") ";
            } else {
                std::cout << "(not found) ";
            }        
        }
        std::cout << std::endl;
    }
    std::cout << "Simulation created successfully." << std::endl;
    
    return world;
}

/**
 * Template function for reading a world from a text file (to be implemented later)
 */
template<typename DataFormat>
mio::abm::World read_world_from_file(const std::string& filename)
{
    mio::unused(filename);
    // This is a template function that can be implemented later
    // For now, it creates a default world
    auto world = mio::abm::World(num_age_groups);
    set_parameters(world.parameters);
    return world;
}

int main(int argc, char** argv)
{
    // Default infection data file path
    std::string infection_data_file = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabm/PanVadere/restaurant/simulation_runs/lu-2020_no_airflow/infections.txt";
    
    // Allow overriding the file path via command line
    if (argc > 1) {
        infection_data_file = argv[1];
    }
    
    std::cout << "Creating restaurant simulation from: " << infection_data_file << std::endl;
    
    // Create the simulation world
    auto world = create_restaurant_world(infection_data_file);
    
    // Set up simulation timeframe
    auto t0 = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(10);
    
    // Initialize simulation
    auto sim = mio::abm::Simulation(t0, std::move(world));
    
    // Create a history object to store time series of infection states
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)
    };
    
    // Run simulation and collect data
    sim.advance(tmax, historyTimeSeries);
    
    // Write results to file
    std::ofstream outfile("panvXabm_results.txt");
    std::get<0>(historyTimeSeries.get_log())
        .print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
    
    std::cout << "Results written to panvXabm_results.txt" << std::endl;
    
    return 0;
}