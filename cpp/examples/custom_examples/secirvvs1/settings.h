#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <iomanip>



#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>
#include <cstddef>
#include <stdexcept>

#include "models/ode_secirvvs/model.h"
#include "memilio/compartments/simulation.h"

// ------------------------------ //
// Handle InfectionStates Parsing //
// ------------------------------ //

using InfectionState = mio::osecirvvs::InfectionState;

inline std::string infection_state_to_string(InfectionState state) {
    switch (state) {
        case InfectionState::SusceptibleNaive: return "SusceptibleNaive";
        case InfectionState::SusceptiblePartialImmunity: return "SusceptiblePartialImmunity";
        case InfectionState::ExposedNaive: return "ExposedNaive";
        case InfectionState::ExposedPartialImmunity: return "ExposedPartialImmunity";
        case InfectionState::ExposedImprovedImmunity: return "ExposedImprovedImmunity";
        case InfectionState::InfectedNoSymptomsNaive: return "InfectedNoSymptomsNaive";
        case InfectionState::InfectedNoSymptomsPartialImmunity: return "InfectedNoSymptomsPartialImmunity";
        case InfectionState::InfectedNoSymptomsImprovedImmunity: return "InfectedNoSymptomsImprovedImmunity";
        case InfectionState::InfectedNoSymptomsNaiveConfirmed: return "InfectedNoSymptomsNaiveConfirmed";
        case InfectionState::InfectedNoSymptomsPartialImmunityConfirmed: return "InfectedNoSymptomsPartialImmunityConfirmed";
        case InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed: return "InfectedNoSymptomsImprovedImmunityConfirmed";
        case InfectionState::InfectedSymptomsNaive: return "InfectedSymptomsNaive";
        case InfectionState::InfectedSymptomsPartialImmunity: return "InfectedSymptomsPartialImmunity";
        case InfectionState::InfectedSymptomsImprovedImmunity: return "InfectedSymptomsImprovedImmunity";
        case InfectionState::InfectedSymptomsNaiveConfirmed: return "InfectedSymptomsNaiveConfirmed";
        case InfectionState::InfectedSymptomsPartialImmunityConfirmed: return "InfectedSymptomsPartialImmunityConfirmed";
        case InfectionState::InfectedSymptomsImprovedImmunityConfirmed: return "InfectedSymptomsImprovedImmunityConfirmed";
        case InfectionState::InfectedSevereNaive: return "InfectedSevereNaive";
        case InfectionState::InfectedSeverePartialImmunity: return "InfectedSeverePartialImmunity";
        case InfectionState::InfectedSevereImprovedImmunity: return "InfectedSevereImprovedImmunity";
        case InfectionState::InfectedCriticalNaive: return "InfectedCriticalNaive";
        case InfectionState::InfectedCriticalPartialImmunity: return "InfectedCriticalPartialImmunity";
        case InfectionState::InfectedCriticalImprovedImmunity: return "InfectedCriticalImprovedImmunity";
        case InfectionState::SusceptibleImprovedImmunity: return "SusceptibleImprovedImmunity";
        case InfectionState::DeadNaive: return "DeadNaive";
        case InfectionState::DeadPartialImmunity: return "DeadPartialImmunity";
        case InfectionState::DeadImprovedImmunity: return "DeadImprovedImmunity";
        case InfectionState::Count: return "Count";
        default: return "Unknown";
    }
}

inline std::vector<InfectionState> query_infection_states(const std::string& state_query) {
    static const std::vector<std::pair<std::string, InfectionState>> all_states = {
        {"SusceptibleNaive", InfectionState::SusceptibleNaive},
        {"SusceptiblePartialImmunity", InfectionState::SusceptiblePartialImmunity},
        {"ExposedNaive", InfectionState::ExposedNaive},
        {"ExposedPartialImmunity", InfectionState::ExposedPartialImmunity},
        {"ExposedImprovedImmunity", InfectionState::ExposedImprovedImmunity},
        {"InfectedNoSymptomsNaive", InfectionState::InfectedNoSymptomsNaive},
        {"InfectedNoSymptomsPartialImmunity", InfectionState::InfectedNoSymptomsPartialImmunity},
        {"InfectedNoSymptomsImprovedImmunity", InfectionState::InfectedNoSymptomsImprovedImmunity},
        {"InfectedNoSymptomsNaiveConfirmed", InfectionState::InfectedNoSymptomsNaiveConfirmed},
        {"InfectedNoSymptomsPartialImmunityConfirmed", InfectionState::InfectedNoSymptomsPartialImmunityConfirmed},
        {"InfectedNoSymptomsImprovedImmunityConfirmed", InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed},
        {"InfectedSymptomsNaive", InfectionState::InfectedSymptomsNaive},
        {"InfectedSymptomsPartialImmunity", InfectionState::InfectedSymptomsPartialImmunity},
        {"InfectedSymptomsImprovedImmunity", InfectionState::InfectedSymptomsImprovedImmunity},
        {"InfectedSymptomsNaiveConfirmed", InfectionState::InfectedSymptomsNaiveConfirmed},
        {"InfectedSymptomsPartialImmunityConfirmed", InfectionState::InfectedSymptomsPartialImmunityConfirmed},
        {"InfectedSymptomsImprovedImmunityConfirmed", InfectionState::InfectedSymptomsImprovedImmunityConfirmed},
        {"InfectedSevereNaive", InfectionState::InfectedSevereNaive},
        {"InfectedSeverePartialImmunity", InfectionState::InfectedSeverePartialImmunity},
        {"InfectedSevereImprovedImmunity", InfectionState::InfectedSevereImprovedImmunity},
        {"InfectedCriticalNaive", InfectionState::InfectedCriticalNaive},
        {"InfectedCriticalPartialImmunity", InfectionState::InfectedCriticalPartialImmunity},
        {"InfectedCriticalImprovedImmunity", InfectionState::InfectedCriticalImprovedImmunity},
        {"SusceptibleImprovedImmunity", InfectionState::SusceptibleImprovedImmunity},
        {"DeadNaive", InfectionState::DeadNaive},
        {"DeadPartialImmunity", InfectionState::DeadPartialImmunity},
        {"DeadImprovedImmunity", InfectionState::DeadImprovedImmunity}
    };

    // Split query into separate subqueries
    std::vector<std::string> and_group_tokens;
    std::vector<std::string> or_group_tokens;

    std::istringstream iss(state_query);
    std::string word;
    while (iss >> word) {
        if (word.find('+') != std::string::npos) {
            std::istringstream substream(word);
            std::string sub;
            while (std::getline(substream, sub, '+')) {
                and_group_tokens.push_back(sub);
            }
        } else {
            or_group_tokens.push_back(word);
        }
    }

    std::vector<InfectionState> result;

    if (!and_group_tokens.empty()) {
        for (const auto& [name, state] : all_states) {
            bool match_all = std::all_of(and_group_tokens.begin(), and_group_tokens.end(), [&](const std::string& token) {
                return name.find(token) != std::string::npos;
            });
            if (match_all) {
                result.push_back(state);
            }
        }
    }

    if (!or_group_tokens.empty()) {
        for (const auto& [name, state] : all_states) {
            bool match_any = std::any_of(or_group_tokens.begin(), or_group_tokens.end(), [&](const std::string& token) {
                return name.find(token) != std::string::npos;
            });
            if (match_any && std::find(result.begin(), result.end(), state) == result.end()) {
                result.push_back(state);
            }
        }
    }

    if (result.empty()) {
        throw std::invalid_argument("No matching InfectionState for query: " + state_query);
    }

    return result;
}

// ------------------------------------------ //
// Store Information in ProblemSettings Class //
// ------------------------------------------ //

class ProblemSettings {
public:

    ProblemSettings(
        int num_control_intervals,
        int control_interval_resolution,
        double t_start,
        double t_end,
        int integrator_resolution,
        std::vector<std::tuple<std::string, std::pair<double, double>, double>> control_bounds,
        std::vector<std::pair<std::string, std::pair<double, double>>> path_constraints,
        std::vector<std::pair<std::string, std::pair<double, double>>> terminal_constraints,
        int num_age_groups,
        int num_graph_nodes
    )
        : numControlIntervals_(num_control_intervals)
        , controlIntervalResolution_(control_interval_resolution)
        , t0_(t_start)
        , tmax_(t_end)
        , integratorResolution_(integrator_resolution)
        , controlBounds_(std::move(control_bounds))
        , pathConstraints_(std::move(path_constraints))
        , terminalConstraints_(std::move(terminal_constraints))
        , numAgeGroups_(num_age_groups)
        , numGraphNodes_(num_graph_nodes)
    {
        numIntervals_ = numControlIntervals_ * controlIntervalResolution_;
        dt_ = (tmax_ - t0_) / (numIntervals_ * integratorResolution_);

        numControls_ = static_cast<int>(controlBounds_.size());
        numPathConstraints_ = static_cast<int>(pathConstraints_.size());
        numTerminalConstraints_ = static_cast<int>(terminalConstraints_.size());
    }

    // Getters (can be extended with setters if needed)
    int numControlIntervals() const { return numControlIntervals_; }
    int controlIntervalResolution() const { return controlIntervalResolution_; }
    double t0() const { return t0_; }
    double tmax() const { return tmax_; }
    int numIntervals() const { return numIntervals_; }
    int integratorResolution() const { return integratorResolution_; }
    double dt() const { return dt_; }

    int numControls() const { return numControls_; }
    int numPathConstraints() const { return numPathConstraints_; }
    int numTerminalConstraints() const { return numTerminalConstraints_; }

    const auto& controlBounds() const { return controlBounds_; }
    const auto& pathConstraints() const { return pathConstraints_; }
    const auto& terminalConstraints() const { return terminalConstraints_; }

    int numAgeGroups() const { return numAgeGroups_; }
    int numGraphNodes() const { return numGraphNodes_; }

    // Print function for debugging
    void print() const {
        std::cout << "Problem Settings:\n"
                  << "  Control Intervals: " << numControlIntervals_ << "\n"
                  << "  Control Interval Resolution: " << controlIntervalResolution_ << "\n"
                  << "  t0: " << t0_ << "\n"
                  << "  tmax: " << tmax_ << "\n"
                  << "  Num Intervals: " << numIntervals_ << "\n"
                  << "  Integrator Resolution: " << integratorResolution_ << "\n"
                  << "  dt: " << dt_ << "\n"
                  << "  Number of Controls: " << numControls_ << "\n"
                  << "  Number of Path Constraints: " << numPathConstraints_ << "\n"
                  << "  Number of Terminal Constraints: " << numTerminalConstraints_ << "\n"
                  << "  Number of Agegroups: " << numAgeGroups_ << "\n"
                  << "  Number of Graph Nodes: " << numGraphNodes_ << "\n";
    }

private:
    int numControlIntervals_;
    int controlIntervalResolution_;
    double t0_;
    double tmax_;
    int numIntervals_;
    int integratorResolution_;
    double dt_;

    std::vector<std::tuple<std::string, std::pair<double, double>, double>> controlBounds_;
    std::vector<std::pair<std::string, std::pair<double, double>>> pathConstraints_;
    std::vector<std::pair<std::string, std::pair<double, double>>> terminalConstraints_;

    int numControls_;
    int numPathConstraints_;
    int numTerminalConstraints_;

    int numAgeGroups_;
    int numGraphNodes_;
};

// ============================================================================
// --- Configure Objective and Constraints ------------------------------------
// ============================================================================

template <typename FP>
std::vector<FP> make_time_grid(FP t0, FP tmax, size_t num_intervals) {
    std::vector<FP> grid(num_intervals + 1);
    FP grid_spacing = (tmax - t0) / num_intervals;
    for (size_t i = 0; i <= num_intervals; i++) {
        grid[i] = t0 + i * grid_spacing;
    }
    return grid;
}


template <typename FP>
FP objective_function(const ProblemSettings& settings, const FP* parameters, size_t n) {
    FP objective = 0.0;

    std::vector<FP> timeSteps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.numIntervals());

    mio::osecirvvs::Model<FP> model(settings.numAgeGroups());
    set_initial_values(model, settings);

    return objective;
}


template <typename FP>
void update_path_constraints(
    const mio::osecirvvs::Model<FP>& model, const ProblemSettings& settings, std::vector<FP>& max_path_values)
{
    assert(max_path_values.size() == settings.numPathConstraints());

    for (size_t i = 0; i < settings.numPathConstraints(); ++i) {
        const auto& constraint = settings.pathConstraints()[i];
        auto states = query_infection_states(constraint.first);

        FP value = 0.0;
        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (const auto& state : states) {
                value += model.populations[{age_group, state}];
            }
        }
        max_path_values[i] = std::max(max_path_values[i], value);
    }
}

template <typename FP>
void update_terminal_constraints(
    const mio::osecirvvs::Model<FP>& model, const ProblemSettings& settings, std::vector<FP>& terminal_values)
{
    assert(terminal_values.size() == settings.numTerminalConstraints());

    for (size_t i = 0; i < settings.numTerminalConstraints(); ++i) {
        const auto& constraint = settings.terminalConstraints()[i];
        auto states = query_infection_states(constraint.first);

        FP value = 0.0;
        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (const auto& state : states) {
                value += model.populations[{age_group, state}];
            }
        }
        terminal_values[i] = value;
    }
}

template <typename FP>
void constraint_functions(const ProblemSettings& settings, const FP* parameters, size_t n, FP* constraints, size_t m) {
    // ----------- // 
    // Setup Model //
    // ----------- // 
    mio::osecirvvs::Model<FP> model(settings.numAgeGroups());
    set_initial_values(model, settings);

    assert(settings.numPathConstraints() + settings.numTerminalConstraints() == m);

    std::vector<FP> max_path_values(settings.numPathConstraints(), FP(0.0));
    std::vector<FP> terminal_values(settings.numTerminalConstraints(), FP(0.0));

    std::vector<FP> timeSteps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.numIntervals());
    size_t num_infection_states = static_cast<size_t>(InfectionState::Count);

    // --------------------------------------------- //
    // Initial Constraint Gathering at t = t0 Before //
    // --------------------------------------------- //
    update_path_constraints(model, settings, max_path_values);

    for (size_t controlIndex = 0; controlIndex < settings.numControlIntervals(); controlIndex++) {

        // ------------------------ //
        // Gather Control Parameter //
        // ------------------------ //
        // FP socialDistancing = parameters[controlIndex * settings.numControls() + static_cast<int>(string_to_control("SocialDistancing"))];
        // FP quarantined = parameters[controlIndex * settings.numControls() + static_cast<int>(string_to_control("Quarantined"))];
        // FP testingRate = parameters[controlIndex * settings.numControls() + static_cast<int>(string_to_control("TestingRate"))];

        // model.parameters.template get<mio::oseair::SocialDistancing<FP>>() = socialDistancing;
        // model.parameters.template get<mio::oseair::Quarantined<FP>>()      = quarantined;
        // model.parameters.template get<mio::oseair::TestingRate<FP>>()      = testingRate;

        for (size_t substep = 0; substep < settings.controlIntervalResolution(); substep++) {

            // -------------- //
            // Run Simulation //
            // -------------- //
            size_t timeStepIndex = controlIndex * settings.controlIntervalResolution() + substep;

            mio::TimeSeries<FP> result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(
                timeSteps[timeStepIndex], timeSteps[timeStepIndex + 1], settings.dt(), model
            );
            const auto& final_state = result.get_last_value();

            for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
                for (size_t state_index = 0; state_index < num_infection_states; state_index++) {
                    size_t idx = age_group.get() * num_infection_states + state_index;
                    model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
                }
            }

            // ----------------------------- // 
            // Gather Constraint Information //
            // ----------------------------- // 
            update_path_constraints(model, settings, max_path_values);
        }
    }

    update_terminal_constraints(model, settings, terminal_values);

    // -------------------------- //
    // Parse Constraints to IPOPT //
    // -------------------------- //
    for (size_t constraintIndex = 0; constraintIndex < settings.numPathConstraints(); constraintIndex++) {
        constraints[constraintIndex] = max_path_values[constraintIndex];
    }

    size_t terminal_constraint_offset = settings.numPathConstraints();
    for (size_t constraintIndex = 0; constraintIndex < settings.numTerminalConstraints(); constraintIndex++) {
        constraints[terminal_constraint_offset + constraintIndex] = terminal_values[constraintIndex];
    }

}



template <typename FP>
void set_initial_values(mio::osecirvvs::Model<FP>& model, const ProblemSettings& settings)
{
    const std::array<std::pair<InfectionState, FP>, 27> initial_population_distribution = {{
        {InfectionState::SusceptibleNaive, 0.0}, // Will be automatically adjusted to match total population size
        {InfectionState::SusceptiblePartialImmunity, 7.0},
        {InfectionState::ExposedNaive, 10.0},
        {InfectionState::ExposedPartialImmunity, 12.0},
        {InfectionState::ExposedImprovedImmunity, 11.0},
        {InfectionState::InfectedNoSymptomsNaive, 13.0},
        {InfectionState::InfectedNoSymptomsPartialImmunity, 14.0},
        {InfectionState::InfectedNoSymptomsImprovedImmunity, 15.0},
        {InfectionState::InfectedNoSymptomsNaiveConfirmed, 13.0},
        {InfectionState::InfectedNoSymptomsPartialImmunityConfirmed, 14.0},
        {InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed, 15.0},
        {InfectionState::InfectedSymptomsNaive, 5.0},
        {InfectionState::InfectedSymptomsPartialImmunity, 6.0},
        {InfectionState::InfectedSymptomsImprovedImmunity, 7.0},
        {InfectionState::InfectedSymptomsNaiveConfirmed, 5.0},
        {InfectionState::InfectedSymptomsPartialImmunityConfirmed, 6.0},
        {InfectionState::InfectedSymptomsImprovedImmunityConfirmed, 7.0},
        {InfectionState::InfectedSevereNaive, 8.0},
        {InfectionState::InfectedSeverePartialImmunity, 2.0},
        {InfectionState::InfectedSevereImprovedImmunity, 1.0},
        {InfectionState::InfectedCriticalNaive, 3.0},
        {InfectionState::InfectedCriticalPartialImmunity, 4.0},
        {InfectionState::InfectedCriticalImprovedImmunity, 5.0},
        {InfectionState::SusceptibleImprovedImmunity, 6.0},
        {InfectionState::DeadNaive, 0.0},
        {InfectionState::DeadPartialImmunity, 0.0},
        {InfectionState::DeadImprovedImmunity, 0.0},
    }};

    // Initialize populations from the initial distribution
    for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
        for (const auto& [infection_state, population] : initial_population_distribution) {
            model.populations[{age_group, infection_state}] = population;
        }
    }

    // Fill the remaining portion as SusceptibleNaive to reach the total of 1000
    const FP total_population_size = 1'000.0;
    for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
        model.populations.template set_difference_from_group_total<mio::AgeGroup>(
            {age_group, InfectionState::SusceptibleNaive}, total_population_size
        );
    }

    // Define constants and simulation dimensions
    const size_t num_time_steps = static_cast<size_t>(settings.tmax());
    const mio::AgeGroup age_0(0);

    // Create references to model parameter components
    auto& params = model.parameters;
    auto& partial_vacc = params.template get<mio::osecirvvs::DailyPartialVaccinations<FP>>();
    auto& full_vacc    = params.template get<mio::osecirvvs::DailyFullVaccinations<FP>>();
    auto& contact_matrix = params.template get<mio::osecirvvs::ContactPatterns<FP>>().get_cont_freq_mat();

    // General settings
    params.template get<mio::osecirvvs::ICUCapacity<FP>>() = 100.0;
    params.template get<mio::osecirvvs::TestAndTraceCapacity<FP>>() = 0.0143;
    params.template get<mio::osecirvvs::DynamicNPIsImplementationDelay<FP>>() = 7;
    params.template get<mio::osecirvvs::Seasonality<FP>>() = 0.2;

    // Resize vaccination schedules
    partial_vacc.resize(mio::SimulationDay(num_time_steps + 1));
    full_vacc.resize(mio::SimulationDay(num_time_steps + 1));

    // Assign daily vaccinations
    const size_t daily_vaccinations = 10;
    for (size_t i = 0; i <= num_time_steps; ++i) {
        FP vacc = static_cast<FP>(i * daily_vaccinations);
        partial_vacc[{age_0, mio::SimulationDay(i)}] = vacc;
        full_vacc[{age_0, mio::SimulationDay(i)}] = vacc;
    }

    // Contact matrix setup
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    contact_matrix[0].add_damping(0.5, mio::SimulationTime(5.0));

    // Durations (in days)
    params.template get<mio::osecirvvs::TimeExposed<FP>>()[age_0] = 3.33;
    params.template get<mio::osecirvvs::TimeInfectedNoSymptoms<FP>>()[age_0] = 1.87;
    params.template get<mio::osecirvvs::TimeInfectedSymptoms<FP>>()[age_0] = 7.0;
    params.template get<mio::osecirvvs::TimeInfectedSevere<FP>>()[age_0] = 6.0;
    params.template get<mio::osecirvvs::TimeInfectedCritical<FP>>()[age_0] = 7.0;

    // Transmission probabilities
    params.template get<mio::osecirvvs::TransmissionProbabilityOnContact<FP>>()[age_0] = 0.15;
    params.template get<mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>>()[age_0] = 0.5;
    params.template get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>>()[age_0] = 0.0;
    params.template get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>>()[age_0] = 0.4;

    // Progression probabilities
    params.template get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>>()[age_0] = 0.2;
    params.template get<mio::osecirvvs::SeverePerInfectedSymptoms<FP>>()[age_0] = 0.1;
    params.template get<mio::osecirvvs::CriticalPerSevere<FP>>()[age_0] = 0.1;
    params.template get<mio::osecirvvs::DeathsPerCritical<FP>>()[age_0] = 0.1;

    // Immunity effects
    params.template get<mio::osecirvvs::ReducExposedPartialImmunity<FP>>()[age_0] = 0.8;
    params.template get<mio::osecirvvs::ReducExposedImprovedImmunity<FP>>()[age_0] = 0.331;
    params.template get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>>()[age_0] = 0.65;
    params.template get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>>()[age_0] = 0.243;
    params.template get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[age_0] = 0.1;
    params.template get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[age_0] = 0.091;
    params.template get<mio::osecirvvs::ReducTimeInfectedMild<FP>>()[age_0] = 0.9;

    // Enforce constraints after assignment
    model.apply_constraints();
}


template<typename FP>
void save_solution(
    const ProblemSettings& settings,
    size_t n, const FP* x, const FP* z_L, const FP* z_U,
    size_t m, const FP* g, const FP* lambda, FP obj_value
) {
    // ----------- // 
    // Setup Model //
    // ----------- // 
    mio::osecirvvs::Model<FP> model(settings.numAgeGroups());
    set_initial_values(model, settings);

    assert(settings.numPathConstraints() + settings.numTerminalConstraints() == m);

    std::vector<FP> timeSteps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.numIntervals());
    size_t num_infection_states = static_cast<size_t>(InfectionState::Count);

    // Open CSV file and write header
    std::ofstream csv_file("population_time_series.csv");
    csv_file << "Time";
    for (size_t i = 0; i < num_infection_states; ++i) {
        csv_file << "," << infection_state_to_string(static_cast<InfectionState>(i));
    }
    csv_file << "\n";

    for (size_t controlIndex = 0; controlIndex < settings.numControlIntervals(); controlIndex++) {
        for (size_t substep = 0; substep < settings.controlIntervalResolution(); substep++) {
            size_t timeStepIndex = controlIndex * settings.controlIntervalResolution() + substep;

            // Simulate next step
            mio::TimeSeries<FP> result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(
                timeSteps[timeStepIndex], timeSteps[timeStepIndex + 1], settings.dt(), model
            );
            const auto& final_state = result.get_last_value();

            // Update model state
            for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
                for (size_t state_index = 0; state_index < num_infection_states; state_index++) {
                    size_t idx = age_group.get() * num_infection_states + state_index;
                    model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
                }
            }

            // Sum over all age groups per infection state
            std::vector<FP> total_by_state(num_infection_states, FP(0.0));
            for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
                for (size_t state_index = 0; state_index < num_infection_states; state_index++) {
                    InfectionState state = static_cast<InfectionState>(state_index);
                    total_by_state[state_index] += model.populations[{age_group, state}];
                }
            }

            // Write to CSV
            csv_file << std::fixed << std::setprecision(6) << timeSteps[timeStepIndex];
            for (const auto& val : total_by_state) {
                csv_file << "," << val;
            }
            csv_file << "\n";
        }
    }

    // Close CSV file
    csv_file.close();
}