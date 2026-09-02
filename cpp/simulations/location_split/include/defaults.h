#ifndef LOCATION_SPLIT_DEFAULTS_H
#define LOCATION_SPLIT_DEFAULTS_H

namespace Config
{
// Validation limits.
constexpr int MIN_RUNS = 1;
constexpr int MAX_RUNS = 1000;

constexpr int MIN_DAYS = 1;
constexpr int MAX_DAYS = 365;

constexpr int MIN_POPULATION = 100;
constexpr int MAX_POPULATION = 10000000;

constexpr double MIN_INFECTION_K = 0.01;
constexpr double MAX_INFECTION_K = 100.0;

// Default values.
constexpr int DEFAULT_RUNS               = 5;
constexpr int DEFAULT_DAYS               = 7;
constexpr int DEFAULT_POPULATION         = 1000;
constexpr int DEFAULT_INITIAL_INFECTIONS = 10;
constexpr double DEFAULT_INFECTION_K     = 22.6;
constexpr const char* DEFAULT_OUTPUT_DIR = "results";
} // namespace Config

#endif // LOCATION_SPLIT_DEFAULTS_H
