#pragma once

#include <string>
#include <utility>
#include <vector>

#include "models/ode_seair/model.h"

// ============================================================================
// --- Configuration Structs -------------------------------------------------
// ============================================================================

struct ProblemSettings {
    int numControlIntervals = 20;
    int pcResolution = 5;
    double t0 = 0.0;
    double tmax = 100.0;
    int numIntervals = numControlIntervals * pcResolution;

    int integratorResolution = 5;
    double dt = (tmax - t0) / (numIntervals * integratorResolution);

    double N = 327'167'434; // total US population

    std::vector<std::tuple<std::string, std::pair<double, double>, double>> controlBounds = {
        {"SocialDistancing", {0.05, 0.5}, 0.1},
        {"Quarantined",      {0.01, 0.3}, 0.05},
        {"TestingRate",      {0.15, 0.3}, 0.2}
    };
    
    // Vector of path constraint bounds: { name, {lower, upper} }
    std::vector<std::pair<std::string, std::pair<double, double>>> pathConstraints = {
        {"Infections", {0.0, 1'000'000}}
    };

    int numControls = static_cast<int>(controlBounds.size());
    int numPathConstraints = static_cast<int>(pathConstraints.size());
};