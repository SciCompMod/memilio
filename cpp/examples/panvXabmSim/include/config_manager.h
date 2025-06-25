#pragma once
#include "simulation_runner.h"
#include <string>

class ConfigManager {
public:
    static SimulationRunner::SimulationConfig get_default_config();
    static SimulationRunner::SimulationConfig load_config_from_args(int argc, char* argv[]);
    static void print_usage();
};