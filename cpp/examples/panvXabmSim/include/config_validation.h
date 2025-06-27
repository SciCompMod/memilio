#pragma once

namespace Config
{
// Validation limits
constexpr int MIN_RUNS = 1;
constexpr int MAX_RUNS = 1000;

constexpr int MIN_DAYS = 1;
constexpr int MAX_DAYS = 365;

constexpr int MIN_POPULATION = 100;
constexpr int MAX_POPULATION = 1000000;

constexpr double MIN_INFECTION_K = 0.1;
constexpr double MAX_INFECTION_K = 10.0;

constexpr int MIN_EVENT_HOURS = 1;
constexpr int MAX_EVENT_HOURS = 24;

// Default values
constexpr int DEFAULT_RUNS               = 1;
constexpr int DEFAULT_DAYS               = 1;
constexpr int DEFAULT_POPULATION         = 1000;
constexpr double DEFAULT_INFECTION_K     = 1.0;
constexpr int DEFAULT_EVENT_HOURS        = 2;
constexpr const char* DEFAULT_OUTPUT_DIR = "/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/cpp/examples/panvXabmSim";
} // namespace Config
