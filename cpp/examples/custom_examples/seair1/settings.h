#pragma once

#include <string>
#include <utility>
#include <vector>
#include <cstddef>

constexpr double INF = 1e19;

#include "models/ode_seair/model.h"
#include "memilio/compartments/simulation.h"

// ============================================================================
// --- Configuration Structs --------------------------------------------------
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

    // Vector of control bounds for various parameters: { name, {lower, upper}, initial value }
    std::vector<std::tuple<std::string, std::pair<double, double>, double>> controlBounds = {
        {"SocialDistancing", {0.05, 0.5}, 0.1},
        {"Quarantined",      {0.01, 0.3}, 0.05},
        {"TestingRate",      {0.15, 0.3}, 0.2}
    };
    
    // Vector of path constraint bounds: { name, {lower, upper} }
    std::vector<std::pair<std::string, std::pair<double, double>>> pathConstraints = {
        {"Infections", {0.0, 1'000'000}},
        {"Dead", {0.0, INF}}
    };

    int numControls = static_cast<int>(controlBounds.size());
    int numPathConstraints = static_cast<int>(pathConstraints.size());
};


// ============================================================================
// --- Configure Objective and Constraints ------------------------------------
// ============================================================================

template <typename FP>
FP objective_function(const ProblemSettings& settings, const FP* parameters, size_t n) {
    FP objective = 0.0;
    for (size_t controlIndex = 0; controlIndex < settings.numControlIntervals; controlIndex++){
        FP socialDistancing = parameters[controlIndex * settings.numControls + 0];
        FP quarantined      = parameters[controlIndex * settings.numControls + 1];
        FP testingRate      = parameters[controlIndex * settings.numControls + 2];

        objective += settings.pcResolution * (-socialDistancing - quarantined + 0.1 * testingRate);
    }
    return objective;
}

template <typename FP>
void constraint_functions(const ProblemSettings& settings, const FP* parameters, size_t n, FP* constraints, size_t m) {

    std::vector<FP> grid(settings.numIntervals + 1);
    FP grid_spacing = (settings.tmax - settings.t0) / settings.numIntervals;
    for (size_t i = 0; i < grid.size(); ++i) {
        grid[i] = settings.t0 + i * grid_spacing;
    }

    mio::oseair::Model<FP> model;
    set_initial_values(model, settings);

    for (size_t controlIndex = 0; controlIndex < settings.numControlIntervals; controlIndex++){

        FP socialDistancing = parameters[controlIndex * settings.numControls + 0];
        FP quarantined      = parameters[controlIndex * settings.numControls + 1];
        FP testingRate      = parameters[controlIndex * settings.numControls + 2];
        model.parameters.template get<mio::oseair::SocialDistancing<FP>>() = socialDistancing;
        model.parameters.template get<mio::oseair::Quarantined<FP>>()      = quarantined;
        model.parameters.template get<mio::oseair::TestingRate<FP>>()      = testingRate;

        for (size_t i = 0; i < settings.pcResolution ; i++)
        {
            size_t gridindex = controlIndex * settings.pcResolution + i;

            auto result = mio::simulate<FP, mio::oseair::Model<FP>>(
                grid[gridindex], grid[gridindex + 1], settings.dt, model
            );

            for (int j = 0; j < (int)mio::oseair::InfectionState::Count; ++j) {
                model.populations[mio::oseair::InfectionState(j)] = result.get_last_value()[j];
            }

            for (size_t constraintIndex = 0; constraintIndex < settings.numPathConstraints; constraintIndex++) {
                const auto& constraint = settings.pathConstraints[constraintIndex];
                constraints[gridindex * settings.numPathConstraints + 0] = 
                    result.get_last_value()[(int)mio::oseair::InfectionState::Infected];
                constraints[gridindex * settings.numPathConstraints + 1] = 
                    result.get_last_value()[(int)mio::oseair::InfectionState::Dead];
            }
        }
    }
}

template <typename FP>
void set_initial_values(mio::oseair::Model<FP>& model, const ProblemSettings& settings)
{
    using IS = mio::oseair::InfectionState;

    const std::array<std::pair<IS, FP>, 6> init = {{
        {IS::Susceptible,  0.9977558755803503 * settings.N},
        {IS::Exposed,      0.0003451395725394549 * settings.N},
        {IS::Asymptomatic, 0.00037846880968213874 * settings.N},
        {IS::Infected,     337072.0},
        {IS::Recovered,    17448.0},
        {IS::Dead,         9619.0}
    }};

    for (const auto& [state, value] : init) {
        model.populations[{mio::Index<IS>(state)}] = value;
    }
}


template<typename FP>
void save_solution(
    const ProblemSettings& settings,
    size_t n, const FP* x, const FP* z_L, const FP* z_U,
    size_t m, const FP* g, const FP* lambda, FP obj_value
) {
    // Load the Data in a .csv file

    std::ofstream outputFile("result.csv");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }
    
    outputFile << "Time,Infected,Recovered,Dead,Susceptible,Exposed,Asymptomatic\n";
    mio::oseair::Model<double> model;
    set_initial_values(model, settings);

    std::vector<double> grid(settings.numIntervals + 1);
    double grid_spacing = (settings.tmax - settings.t0) / settings.numIntervals;
    for (size_t i = 0; i < grid.size(); ++i) {
        grid[i] = settings.t0 + i * grid_spacing;
    }

    for (size_t controlIndex = 0; controlIndex < settings.numControlIntervals; controlIndex++){
        double socialDistancing = x[controlIndex * settings.numControls + 0];
        double quarantined      = x[controlIndex * settings.numControls + 1];
        double testingRate      = x[controlIndex * settings.numControls + 2];
        model.parameters.template get<mio::oseair::SocialDistancing<double>>() = socialDistancing;
        model.parameters.template get<mio::oseair::Quarantined<double>>()      = quarantined;
        model.parameters.template get<mio::oseair::TestingRate<double>>()      = testingRate;

        for (size_t i = 0; i < settings.pcResolution; i++)
        {
            size_t gridindex = controlIndex * settings.pcResolution + i;
            auto result = mio::simulate<double, mio::oseair::Model<double>>(
                grid[gridindex], grid[gridindex + 1], settings.dt, model
            );

            for (int j = 0; j < (int)mio::oseair::InfectionState::Count; ++j) {
                model.populations[mio::oseair::InfectionState(j)] = result.get_last_value()[j];
            }
            outputFile << grid[gridindex] << "," << result.get_last_value()[(int)mio::oseair::InfectionState::Infected] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Recovered] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Dead] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Susceptible] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Exposed] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Asymptomatic] << "\n";
        }
    }
    outputFile.close();
    std::cout << "Solution saved to output.csv\n";


    // Next save the optimal control values
    std::ofstream controlFile("controls.csv");
    if (!controlFile.is_open()) {
        std::cerr << "Error opening control output file!" << std::endl;
        return;
    }
    controlFile << "ControlInterval,SocialDistancing,Quarantined,TestingRate\n";
    controlFile << std::scientific << std::setprecision(6);

    for (size_t controlIndex = 0; controlIndex < settings.numControlIntervals; controlIndex++){
        controlFile << controlIndex * settings.pcResolution << "," 
                    << x[controlIndex * settings.numControls + 0] << ","
                    << x[controlIndex * settings.numControls + 1] << ","
                    << x[controlIndex * settings.numControls + 2] << "\n";
    }
    // Save the last control interval 
    size_t lastControlIndex = settings.numControlIntervals - 1;
    controlFile << (lastControlIndex+1) * settings.pcResolution << ","
                << x[lastControlIndex * settings.numControls + 0] << ","
                << x[lastControlIndex * settings.numControls + 1] << ","
                << x[lastControlIndex * settings.numControls + 2] << "\n";

    controlFile.close();
    std::cout << "Control values saved to control.csv\n";
}
