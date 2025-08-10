#pragma once

namespace Config
{
// Validation limits
constexpr int MIN_RUNS = 1;
constexpr int MAX_RUNS = 1000;

constexpr int MIN_DAYS = 1;
constexpr int MAX_DAYS = 30;

constexpr int MIN_POPULATION = 1000;
constexpr int MAX_POPULATION = 10000;

constexpr double MIN_INFECTION_K = 0.01;
constexpr double MAX_INFECTION_K = 100.0;

// Default values
constexpr int DEFAULT_RUNS               = 5;
constexpr int DEFAULT_DAYS               = 7;
constexpr int DEFAULT_POPULATION         = 1000;
constexpr double DEFAULT_INFECTION_K     = 20.0;
constexpr int DEFAULT_EVENT_HOURS        = 2;
constexpr const char* DEFAULT_BASE_DIR   = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim";
constexpr const char* DEFAULT_OUTPUT_DIR = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results";
} // namespace Config
