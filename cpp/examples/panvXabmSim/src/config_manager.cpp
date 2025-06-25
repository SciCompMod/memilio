#include "../include/config_manager.h"
#include <iostream>

SimulationRunner::SimulationConfig ConfigManager::get_default_config()
{
    SimulationRunner::SimulationConfig config;
    config.infection_data_file =
        "/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/cpp/examples/panvXabmSim/data/restaurant/simulation_runs/"
        "lu-2020_airflow_inlet_right_outlet_left/infections.txt";
    config.input_dir       = "/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/cpp/examples/panvXabmSim/results";
    config.n_persons       = 1000;
    config.simulation_days = 20;
    return config;
}

SimulationRunner::SimulationConfig ConfigManager::load_config_from_args(int argc, char* argv[])
{
    auto config = get_default_config();

    // Simple argument parsing - extend as needed
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--persons" && i + 1 < argc) {
            config.n_persons = std::stoi(argv[++i]);
        }
        else if (arg == "--days" && i + 1 < argc) {
            config.simulation_days = std::stoi(argv[++i]);
        }
        else if (arg == "--input-file" && i + 1 < argc) {
            config.infection_data_file = argv[++i];
        }
        else if (arg == "--output-dir" && i + 1 < argc) {
            config.input_dir = argv[++i];
        }
        else if (arg == "--help") {
            print_usage();
            exit(0);
        }
    }

    return config;
}

void ConfigManager::print_usage()
{
    std::cout << "Usage: panvXabm [options]\n"
              << "Options:\n"
              << "  --persons N          Number of persons to add (default: 1000)\n"
              << "  --days N             Simulation days (default: 20)\n"
              << "  --input-file FILE    Infection data file path\n"
              << "  --output-dir DIR     Output directory path\n"
              << "  --help               Show this help message\n";
}