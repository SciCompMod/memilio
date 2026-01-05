// clang-format off

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
#include <filesystem>
#include <string>
#include "memilio/io/io.h"
#include "memilio/io/mobility_io.h"

#include "models/ode_secirvvs/model.h"
#include "memilio/compartments/simulation.h"

// ------------------------------------------ //
// Store Information in ProblemSettings Class //
// ------------------------------------------ //

enum class PathConstraintMode {Individual, GlobalMax};

class ProblemSettings {
public:
    ProblemSettings(
        int num_control_intervals,
        int pc_resolution,
        double t0,
        double tmax,
        int integrator_resolution,
        PathConstraintMode pc_mode,
        std::vector<std::tuple<std::string, std::pair<double, double>, double>> control_bounds,
        std::vector<std::pair<std::string, std::pair<double, double>>> path_constraints,
        std::vector<std::pair<std::string, std::pair<double, double>>> terminal_constraints
    )
        : num_control_intervals_(num_control_intervals)
        , pc_resolution_(pc_resolution)
        , t0_(t0)
        , tmax_(tmax)
        , integrator_resolution_(integrator_resolution)
        , pc_mode_(pc_mode)
        , control_bounds_(std::move(control_bounds))
        , path_constraints_(std::move(path_constraints))
        , terminal_constraints_(std::move(terminal_constraints))
    {
        num_intervals_ = num_control_intervals_ * pc_resolution_;
        dt_ = (tmax_ - t0_) / (num_intervals_ * integrator_resolution_);

        num_controls_ = static_cast<int>(control_bounds_.size());
        num_path_constraints_ = static_cast<int>(path_constraints_.size());
        num_terminal_constraints_ = static_cast<int>(terminal_constraints_.size());
    }

    int numControlIntervals() const { return num_control_intervals_; }
    int pcResolution() const { return pc_resolution_; }
    double t0() const { return t0_; }
    double tmax() const { return tmax_; }
    int numIntervals() const { return num_intervals_; }
    int integratorResolution() const { return integrator_resolution_; }
    double dt() const { return dt_; }
    PathConstraintMode pathConstraintMode() const { return pc_mode_; }

    int numControls() const { return num_controls_; }
    int numPathConstraints() const { return num_path_constraints_; }
    int numTerminalConstraints() const { return num_terminal_constraints_; }
    const auto& controlBounds() const { return control_bounds_; }
    const auto& pathConstraints() const { return path_constraints_; }
    const auto& terminalConstraints() const { return terminal_constraints_; }

private:
    int num_control_intervals_;
    int pc_resolution_;
    double t0_;
    double tmax_;
    int num_intervals_;
    int integrator_resolution_;
    double dt_;
    double infinity_;
    PathConstraintMode pc_mode_;

    std::vector<std::tuple<std::string, std::pair<double, double>, double>> control_bounds_;
    std::vector<std::pair<std::string, std::pair<double, double>>> path_constraints_;
    std::vector<std::pair<std::string, std::pair<double, double>>> terminal_constraints_;

    int num_controls_;
    int num_path_constraints_;
    int num_terminal_constraints_;
};

class SimulationSettings {
public:
    SimulationSettings(
        int num_age_groups,
        std::filesystem::path data_directory,
        std::filesystem::path initialization_file
    )
        : num_age_groups_(num_age_groups)
        , data_directory_(data_directory)
        , initialization_file_(initialization_file)
    {}

    int numAgeGroups() const { return num_age_groups_; }
    std::filesystem::path dataDirectory() const { return data_directory_; }
    std::filesystem::path initializationFile() const { return initialization_file_; }

private:
    int num_age_groups_;
    std::filesystem::path data_directory_;
    std::filesystem::path initialization_file_;
};


// ------------------------------ //
// Handle InfectionStates Parsing //
// ------------------------------ //

using InfectionState = mio::osecirvvs::InfectionState;

inline std::string infection_state_to_string(InfectionState state) {
    switch (state) {
        case InfectionState::SusceptibleNaive: return "SusceptibleNaive";
        case InfectionState::SusceptiblePartialImmunity: return "SusceptiblePartialImmunity";
        case InfectionState::SusceptibleImprovedImmunity: return "SusceptibleImprovedImmunity";
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
        case InfectionState::DeadNaive: return "DeadNaive";
        case InfectionState::DeadPartialImmunity: return "DeadPartialImmunity";
        case InfectionState::DeadImprovedImmunity: return "DeadImprovedImmunity";
        case InfectionState::Count: return "Count";
        default: return "InfectionState Unknown";
    }
}

inline std::vector<InfectionState> query_infection_states(const std::string& state_query) {
    static const std::vector<std::pair<std::string, InfectionState>> all_states = {
        {"SusceptibleNaive", InfectionState::SusceptibleNaive},
        {"SusceptiblePartialImmunity", InfectionState::SusceptiblePartialImmunity},
        {"SusceptibleImprovedImmunity", InfectionState::SusceptibleImprovedImmunity},
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


// ------------------------- //
// Handle Control Parameters //
// ------------------------- //

enum class ControlParameter {
    SchoolClosure,
    HomeOffice,
    PhysicalDistancingSchool,
    PhysicalDistancingWork,
    PhysicalDistancingOther,
    Count
};

inline ControlParameter string_to_control(const std::string& control_name) {
    static const std::vector<std::pair<std::string, ControlParameter>> control_pairs = {
        {"SchoolClosure", ControlParameter::SchoolClosure},
        {"HomeOffice", ControlParameter::HomeOffice},
        {"PhysicalDistancingSchool", ControlParameter::PhysicalDistancingSchool},
        {"PhysicalDistancingWork", ControlParameter::PhysicalDistancingWork},
        {"PhysicalDistancingOther", ControlParameter::PhysicalDistancingOther},
    };
    auto it = std::find_if(
        control_pairs.begin(), control_pairs.end(),
        [&](const auto& pair) {
            return pair.first == control_name;
        }
    );
    if (it != control_pairs.end()) {
        return it->second;
    }
    throw std::runtime_error("Invalid control name: " + control_name);
}


// ============================================================
// --- Configure the Model ------------------------------------
// ============================================================

template <typename FP>
void set_initial_values(
    mio::osecirvvs::Model<FP>& model, const ProblemSettings& problem_settings, const SimulationSettings& simulation_settings
) {
    auto& params = model.parameters;
    mio::AgeGroup num_age_groups = params.get_num_groups();

    params.template get<mio::osecirvvs::ICUCapacity<FP>>() = 100;
    params.template get<mio::osecirvvs::TestAndTraceCapacity<FP>>() = 0.0143;
    params.template get<mio::osecirvvs::Seasonality<FP>>() = 0.2;
    params.template get<mio::osecirvvs::DynamicNPIsImplementationDelay<FP>>() = 7;

    size_t tmax_days = static_cast<size_t>(problem_settings.tmax());
    params.template get<mio::osecirvvs::DailyPartialVaccinations<FP>>().resize(mio::SimulationDay(tmax_days + 1));
    params.template get<mio::osecirvvs::DailyFullVaccinations<FP>>().resize(mio::SimulationDay(tmax_days + 1));

    size_t daily_vaccinations = 10;
    for (mio::AgeGroup age_group = 0; age_group < num_age_groups; age_group++) {
        for (size_t current_day = 0; current_day <= tmax_days; current_day++) {
            mio::Index<mio::AgeGroup, mio::SimulationDay> index(age_group, mio::SimulationDay(current_day));
            FP num_vaccinations = static_cast<FP>(current_day * daily_vaccinations);
            params.template get<mio::osecirvvs::DailyPartialVaccinations<FP>>()[index] = num_vaccinations;
            params.template get<mio::osecirvvs::DailyFullVaccinations<FP>>()[index] = num_vaccinations;
        }
    }

    auto set_all_groups = [&](auto Tag, FP value) {
        auto& age_group_params = params.template get<decltype(Tag)>();  
        for (mio::AgeGroup age_group = 0; age_group < num_age_groups; age_group++) {
            age_group_params[age_group] = value;
        }
    };

    // — times (all groups same value)
    set_all_groups(mio::osecirvvs::TimeExposed<FP>{}, 3.33);
    set_all_groups(mio::osecirvvs::TimeInfectedNoSymptoms<FP>{}, 1.87);
    set_all_groups(mio::osecirvvs::TimeInfectedSymptoms<FP>{}, 7);
    set_all_groups(mio::osecirvvs::TimeInfectedSevere<FP>{}, 6);
    set_all_groups(mio::osecirvvs::TimeInfectedCritical<FP>{}, 7);
    // — probabilities (all groups same value)
    set_all_groups(mio::osecirvvs::TransmissionProbabilityOnContact<FP>{}, 0.15);
    set_all_groups(mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>{}, 0.5);
    set_all_groups(mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>{}, 0.0);
    set_all_groups(mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>{}, 0.4);
    set_all_groups(mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>{}, 0.2);
    set_all_groups(mio::osecirvvs::SeverePerInfectedSymptoms<FP>{}, 0.1);
    set_all_groups(mio::osecirvvs::CriticalPerSevere<FP>{}, 0.1);
    set_all_groups(mio::osecirvvs::DeathsPerCritical<FP>{}, 0.1);
    // — immunity reductions (all groups same value)
    set_all_groups(mio::osecirvvs::ReducExposedPartialImmunity<FP>{}, 0.8);
    set_all_groups(mio::osecirvvs::ReducExposedImprovedImmunity<FP>{}, 0.331);
    set_all_groups(mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>{}, 0.65);
    set_all_groups(mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>{}, 0.243);
    set_all_groups(mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>{}, 0.1);
    set_all_groups(mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>{}, 0.091);
    set_all_groups(mio::osecirvvs::ReducTimeInfectedMild<FP>{}, 0.9);

    std::vector<std::string> contact_locations = {"home", "school_pf_eig", "work", "other"};

    mio::ContactMatrixGroup<FP> contact_matrices(contact_locations.size(), num_age_groups.get());
    for (size_t i = 0; i < contact_locations.size(); ++i) {
        std::string location = contact_locations[i];
        std::filesystem::path baseline_file_path = simulation_settings.dataDirectory() / "contacts" / ("baseline_" + location + ".txt");

        mio::IOResult<Eigen::MatrixX<ScalarType>> baseline_result = mio::read_mobility_plain(baseline_file_path.string());
        if (!baseline_result) {
            std::cerr << "Error reading file: " << baseline_result.error().message() << std::endl;
        }

        contact_matrices[i].get_baseline() = baseline_result.value().cast<FP>();
        contact_matrices[i].get_minimum() = 
            Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic>::Zero(
                num_age_groups.get(), num_age_groups.get()
            );
    }

    params.template get<mio::osecirvvs::ContactPatterns<FP>>() = 
        mio::UncertainContactMatrix<FP>(std::move(contact_matrices));

    // — initial populations from file
    Eigen::MatrixXd initial_population;
    {
        std::ifstream ifs(simulation_settings.initializationFile());
        if (!ifs) {
            throw std::runtime_error("Cannot open initialization file: " + simulation_settings.initializationFile().string());
        }
        initial_population.resize(num_age_groups.get(), static_cast<size_t>(InfectionState::Count));
        for (size_t row = 0; row < num_age_groups.get(); row++) {
            for (size_t column = 0; column < static_cast<size_t>(InfectionState::Count); column++) {
                ifs >> initial_population(row, column);
                if (!ifs) {
                    throw std::runtime_error("Error reading init file at row " + std::to_string(row) + ", col " + std::to_string(column));
                }
            }
        }
    }

    for (size_t column = 0; column < initial_population.cols(); column++) {
        for (size_t row = 0; row < initial_population.rows(); row++){
            mio::Index<mio::AgeGroup, InfectionState> index{
                mio::AgeGroup(row), InfectionState(column)
            };
            model.populations[index] = static_cast<FP>(initial_population(row, column));
        }
    }
    
    model.apply_constraints();
}


// ============================================================================
// --- Configure Objective and Constraints ------------------------------------
// ============================================================================

template <typename FP>
std::vector<FP> make_time_grid(FP t0, FP tmax, size_t num_intervals) {
    // ------------------------------------------------------------------- //
    // Create a uniform distribution of points in the specified timeframe. //
    // ------------------------------------------------------------------ //
    std::vector<FP> grid(num_intervals + 1);
    FP grid_spacing = (tmax - t0) / num_intervals;
    for (size_t i = 0; i <= num_intervals; i++) {
        grid[i] = t0 + i * grid_spacing;
    }
    return grid;
}

template <typename FP>
void update_path_values(
    const mio::osecirvvs::Model<FP>& model, const ProblemSettings& settings, int interval, std::vector<FP>& path_values
) {
    // ------------------------------------------------------------------ //
    // Consider the model to be at the 'interval' simulation step.        //
    // Query and fill the constraints of the model into 'path_values'     //
    // ------------------------------------------------------------------ //
    int num_intervals = settings.numIntervals();
    int num_path_constraints = settings.numPathConstraints();
    assert( static_cast<int>(path_values.size()) == num_intervals * num_path_constraints);

    for (int constraintIndex = 0; constraintIndex < num_path_constraints; constraintIndex++) {
        const auto& constraint = settings.pathConstraints()[constraintIndex];
        auto states = query_infection_states(std::get<0>(constraint));

        FP value = 0.0;
        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (const auto& state : states) {
                value += model.populations[{age_group, state}].value();
            }
        }
        path_values[interval * num_path_constraints + constraintIndex] = value;
    }
}

template <typename FP>
void update_terminal_values(
    const mio::osecirvvs::Model<FP>& model, const ProblemSettings& settings, std::vector<FP>& terminal_values)
{
    // ------------------------------------------------------------------ //
    // Consider the model to be at the final simulation step.             //
    // Query and fill the constraints of the model into 'terminal_values' //
    // ------------------------------------------------------------------ //
    int num_terminal_constraints = settings.numTerminalConstraints();
    assert(static_cast<int>(terminal_values.size()) == num_terminal_constraints);

    for (size_t i = 0; i < num_terminal_constraints; ++i) {
        const auto& constraint = settings.terminalConstraints()[i];
        auto states = query_infection_states(constraint.first);

        FP value = 0.0;
        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (const auto& state : states) {
                value += model.populations[{age_group, state}].value();
            }
        }
        terminal_values[i] = value;
    }
}

template <typename FP>
void fill_constraints(
    const ProblemSettings& settings, FP* constraints, std::vector<FP>& path_values, std::vector<FP>& terminal_values)
{
    // --------------------------------------------------------------------- //
    // Parse the constraints provided in 'path_values' and 'terminal_values' //
    // into the IPOPT format 'FP* constraints'                               //
    // --------------------------------------------------------------------- //
    int num_intervals = settings.numIntervals();
    int num_path_constraints = settings.numPathConstraints();
    int num_terminal_constraints = settings.numTerminalConstraints();
    // ----------------------- //
    // Update Path Constraints //
    // ----------------------- //
    int effective_pc_nodes = 0;
    switch (settings.pathConstraintMode()) {
        case PathConstraintMode::Individual:
            effective_pc_nodes = num_intervals;
            for (size_t i = 0; i < path_values.size(); ++i) {
                constraints[i] = path_values[i];
            }
            break;
        case PathConstraintMode::GlobalMax:
            effective_pc_nodes = 1;
            for (size_t constraintIndex = 0; constraintIndex < num_path_constraints; constraintIndex++) {
                FP max_value = path_values[constraintIndex];
                for (size_t interval = 1; interval < num_intervals; interval++) {
                    size_t idx = interval * num_path_constraints + constraintIndex;
                    if (path_values[idx] > max_value) {
                        max_value = path_values[idx];
                    }
                }
                constraints[constraintIndex] = max_value;
            }
            break;
        default:
            throw std::runtime_error("Unsupported path constraint mode");
    }
    // --------------------------- //
    // Update Terminal Constraints //
    // --------------------------- //
    int terminal_offset = effective_pc_nodes * num_path_constraints;
    for (size_t constraintIndex = 0; constraintIndex < num_terminal_constraints; constraintIndex++) {
        constraints[terminal_offset + constraintIndex] = terminal_values[constraintIndex];
    }
}


enum class Intervention{Home,SchoolClosure,HomeOffice,GatheringBanFacilitiesClosure,PhysicalDistanceAndMasks,SeniorAwareness,Count};
enum class InterventionLevel{Main,PhysicalDistanceAndMasks,SeniorAwareness,Holidays,Count};
enum class ContactLocation{Home,School,Work,Other,Count};

template <typename FP>
FP objective_function(
    mio::osecirvvs::Model<FP> model, const ProblemSettings& settings, const FP* parameters, size_t n) 
{
    // ------------------------------------------------------------------ //
    // Evaluate the objective function of the model.                      //
    // Step 1. Define dampings based on 'const FP* parameters'.           //
    // Step 2. Evaluate the objective function based on                   //
    //         the parameters and the infection states in the simulation. // 
    // ------------------------------------------------------------------ //
    int num_control_intervals = settings.numControlIntervals();
    int pc_resolution = settings.pcResolution();
    int num_intervals = settings.numIntervals();
    int num_controls = settings.numControls();
    int num_path_constraints = settings.numPathConstraints();
    int num_terminal_constraints = settings.numTerminalConstraints();
    double dt = settings.dt();

    FP objective = 0.0;

    std::vector<FP> timeSteps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.numIntervals());
    size_t num_infection_states = static_cast<size_t>(InfectionState::Count);

    // --------------- //
    // Define Dampings //
    // --------------- //
    auto damping_helper = [=](
        mio::SimulationTime<FP> time, FP min, FP max, 
        mio::DampingLevel damping_level, mio::DampingType damping_type, 
        const std::vector<size_t> location, Eigen::VectorX<FP> group_weights
    ) {
        auto damping_value = mio::UncertainValue<FP>(0.5 * (max + min));
        damping_value.set_distribution(mio::ParameterDistributionUniform(ad::value(min), ad::value(max)));
        return mio::DampingSampling<FP>(damping_value, damping_level, damping_type, time, location, group_weights);
    };

    auto group_weights_all = Eigen::VectorX<FP>::Constant(model.parameters.get_num_groups().get(), 1.0);

    auto set_school_closure = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::Main)),
                                            mio::DampingType(size_t(Intervention::SchoolClosure)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto set_home_office = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::Main)),
                                            mio::DampingType(size_t(Intervention::HomeOffice)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto set_physical_distancing_school = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(size_t(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto set_physical_distancing_work = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(size_t(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto set_physical_distancing_other = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(size_t(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Other)}, group_weights_all);
    };

    // --------------------------------------- //
    // Create dampings from control parameters //
    // --------------------------------------- //
    mio::UncertainContactMatrix<FP>& contacts = model.parameters.template get<mio::osecirvvs::ContactPatterns<FP>>();
    std::vector<mio::DampingSampling<FP>>& contact_dampings = contacts.get_dampings();

    contact_dampings.clear();
    for (size_t controlIndex = 0; controlIndex < num_control_intervals; controlIndex++) {

        auto param_at = [&](const std::string& name) {
            return parameters[static_cast<size_t>(string_to_control(name)) + controlIndex * num_controls];
        };
        FP school_closure             = param_at("SchoolClosure");
        FP home_office                = param_at("HomeOffice");
        FP physical_distancing_school = param_at("PhysicalDistancingSchool");
        FP physical_distancing_work   = param_at("PhysicalDistancingWork");
        FP physical_distancing_other  = param_at("PhysicalDistancingOther");

        size_t timeStepIndex = controlIndex * pc_resolution;
        mio::SimulationTime<FP> time(timeSteps[timeStepIndex]);

        contact_dampings.push_back(set_school_closure(time, 1.0 * school_closure));
        contact_dampings.push_back(set_home_office(time, 0.25 * home_office));
        contact_dampings.push_back(set_physical_distancing_school(time, 0.25 * physical_distancing_school));
        contact_dampings.push_back(set_physical_distancing_work(time, 0.25 * physical_distancing_work));
        contact_dampings.push_back(set_physical_distancing_other(time, 0.35 * physical_distancing_other));
    }
    contacts.make_matrix();

    // ---------------- //
    // Start simulation //
    // ---------------- //
    for (size_t interval = 0; interval < num_intervals; interval++) {
        FP abs_tol = 1e-12;
        FP rel_tol = 1e-8;
        FP dt_min = std::numeric_limits<ScalarType>::min();
        FP dt_max = dt;
        auto integrator = std::make_unique<
            mio::ControlledStepperWrapper<FP, boost::numeric::odeint::runge_kutta_fehlberg78>
        >(abs_tol, rel_tol, dt_min, dt_max);

        const int controlIndex = interval / pc_resolution;
        auto param_at = [&](const std::string& name) {
            return parameters[static_cast<size_t>(string_to_control(name)) + controlIndex * num_controls];
        };
        FP school_closure             = param_at("SchoolClosure");
        FP home_office                = param_at("HomeOffice");
        FP physical_distancing_school = param_at("PhysicalDistancingSchool");
        FP physical_distancing_work   = param_at("PhysicalDistancingWork");
        FP physical_distancing_other  = param_at("PhysicalDistancingOther");

        mio::TimeSeries<FP> result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(
            timeSteps[interval], timeSteps[interval+1], dt, model, std::move(integrator)
        );
        const auto& final_state = result.get_last_value();

        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (size_t state_index = 0; state_index < num_infection_states; state_index++) {
                size_t idx = age_group.get() * num_infection_states + state_index;
                model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
            }
        }

        // ------------------ //
        // Objective Function //
        // ------------------ //
        objective += (school_closure + home_office + 
                     physical_distancing_school + physical_distancing_work + physical_distancing_other) / num_intervals;
    }

    return objective;
}


template <typename FP>
void constraint_functions(
    mio::osecirvvs::Model<FP> model, const ProblemSettings& settings, 
    const FP* parameters, size_t n, FP* constraints, size_t m
) {
    // ----------------------------------------------------------------- //
    // Evaluate the constraints on the model.                            //
    // Step 1. Define dampings based on 'const FP* parameters'.          //
    // Step 2. Gather information in 'path_values' in 'terminal_values'. //
    // Step 3. Fill information into 'FP* constraints'.                  //
    // ----------------------------------------------------------------- //
    int num_control_intervals = settings.numControlIntervals();
    int pc_resolution = settings.pcResolution();
    int num_intervals = settings.numIntervals();
    int num_controls = settings.numControls();
    int num_path_constraints = settings.numPathConstraints();
    int num_terminal_constraints = settings.numTerminalConstraints();
    double dt = settings.dt();

    std::vector<FP> path_values(num_intervals * num_path_constraints, FP(0.0));
    std::vector<FP> terminal_values(num_terminal_constraints, FP(0.0));

    std::vector<FP> timeSteps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.numIntervals());
    size_t num_infection_states = static_cast<size_t>(InfectionState::Count);

    // --------------- //
    // Define Dampings //
    // --------------- //
    auto damping_helper = [=](
        mio::SimulationTime<FP> time, FP min, FP max, 
        mio::DampingLevel damping_level, mio::DampingType damping_type, 
        const std::vector<size_t> location, Eigen::VectorX<FP> group_weights
    ) {
        auto damping_value = mio::UncertainValue<FP>(0.5 * (max + min));
        damping_value.set_distribution(mio::ParameterDistributionUniform(ad::value(min), ad::value(max)));
        return mio::DampingSampling<FP>(damping_value, damping_level, damping_type, time, location, group_weights);
    };

    auto group_weights_all = Eigen::VectorX<FP>::Constant(model.parameters.get_num_groups().get(), 1.0);

    auto set_school_closure = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::Main)),
                                            mio::DampingType(size_t(Intervention::SchoolClosure)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto set_home_office = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::Main)),
                                            mio::DampingType(size_t(Intervention::HomeOffice)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto set_physical_distancing_school = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(size_t(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto set_physical_distancing_work = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(size_t(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto set_physical_distancing_other = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(size_t(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Other)}, group_weights_all);
    };

    // --------------------------------------- //
    // Create dampings from control parameters //
    // --------------------------------------- //
    mio::UncertainContactMatrix<FP>& contacts = model.parameters.template get<mio::osecirvvs::ContactPatterns<FP>>();
    std::vector<mio::DampingSampling<FP>>& contact_dampings = contacts.get_dampings();

    for (size_t controlIndex = 0; controlIndex < num_control_intervals; controlIndex++) {

        auto param_at = [&](const std::string& name) {
            return parameters[static_cast<size_t>(string_to_control(name)) + controlIndex * num_controls];
        };
        FP school_closure             = param_at("SchoolClosure");
        FP home_office                = param_at("HomeOffice");
        FP physical_distancing_school = param_at("PhysicalDistancingSchool");
        FP physical_distancing_work   = param_at("PhysicalDistancingWork");
        FP physical_distancing_other  = param_at("PhysicalDistancingOther");

        size_t timeStepIndex = controlIndex * pc_resolution;
        mio::SimulationTime<FP> time(timeSteps[timeStepIndex]);

        contact_dampings.push_back(set_school_closure(time, 1.0 * school_closure));
        contact_dampings.push_back(set_home_office(time, 0.25 * home_office));
        contact_dampings.push_back(set_physical_distancing_school(time, 0.25 * physical_distancing_school));
        contact_dampings.push_back(set_physical_distancing_work(time, 0.25 * physical_distancing_work));
        contact_dampings.push_back(set_physical_distancing_other(time, 0.35 * physical_distancing_other));
    }
    contacts.make_matrix();

    // ----------------- //
    // Define integrator //
    // ----------------- //
    FP abs_tol = 1e-12;
    FP rel_tol = 1e-8;
    FP dt_min = std::numeric_limits<ScalarType>::min();
    FP dt_max = dt;
    auto integrator = std::make_unique<
        mio::ControlledStepperWrapper<FP, boost::numeric::odeint::runge_kutta_fehlberg78>
    >(abs_tol, rel_tol, dt_min, dt_max);

    // ---------------- //
    // Start simulation //
    // ---------------- //
    update_path_values(model, settings, 0, path_values);
    for (size_t interval = 0; interval < num_intervals; interval++) {

        mio::TimeSeries<FP> result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(
            timeSteps[interval], timeSteps[interval+1], dt, model, std::move(integrator)
        );
        const auto& final_state = result.get_last_value();

        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (size_t state_index = 0; state_index < num_infection_states; state_index++) {
                size_t idx = age_group.get() * num_infection_states + state_index;
                model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
            }
        }
        update_path_values(model, settings, interval, path_values);
    }
    update_terminal_values(model, settings, terminal_values);

    // ---------------------------------------- //
    // Parse constraints in the format of IPOPT //
    // ---------------------------------------- //
    fill_constraints(settings, constraints, path_values, terminal_values);
}



template<typename FP>
void save_solution(
    mio::osecirvvs::Model<FP> model, const ProblemSettings& settings,
    size_t n, const FP* x, const FP* z_L, const FP* z_U,
    size_t m, const FP* g, const FP* lambda, FP obj_value)
{
    switch (settings.pathConstraintMode()) {
        case PathConstraintMode::Individual:
            std::cout << "\nIndividual Path Constraints:\n";
            for (size_t i = 0; i < settings.numPathConstraints(); ++i) {
                std::cout<< settings.pathConstraints()[i].first << ": " << std::endl;
                for(size_t j = 0; j < settings.numIntervals(); ++j) {
                    std::cout <<  g[i + j * settings.numPathConstraints()];
                    if(j < settings.numIntervals() - 1) {
                        std::cout << ", ";
                    } else {
                        std::cout << std::endl;
                    }
                }
            }
            std::cout<< "\nTerminal Constraints:\n";
            for (size_t i = 0; i < settings.numTerminalConstraints(); ++i) {
                std::cout << settings.terminalConstraints()[i].first << ": " << g[i + settings.numPathConstraints() * settings.numIntervals()] << std::endl;
            }
            break;
        case PathConstraintMode::GlobalMax:
            std::cout<< "\nGlobal Max Path Constraints:\n";
            for (size_t i = 0; i < settings.numPathConstraints(); ++i) {
                std::cout << settings.pathConstraints()[i].first << ": " << g[i] << std::endl;
            }
            std::cout<< "\nTerminal Constraints:\n";
            for (size_t i = 0; i < settings.numTerminalConstraints(); ++i) {
                std::cout << settings.terminalConstraints()[i].first << ": " << g[i + settings.numPathConstraints()] << std::endl;
            }
            break;
        default:
            throw std::runtime_error("Unsupported path constraint mode");
    }

    int num_control_intervals = settings.numControlIntervals();
    int pc_resolution = settings.pcResolution();
    int num_intervals = settings.numIntervals();
    int num_controls = settings.numControls();
    int num_path_constraints = settings.numPathConstraints();
    int num_terminal_constraints = settings.numTerminalConstraints();
    double dt = settings.dt();

    {
        std::ofstream control_file("control_parameters.csv");
        control_file << "Time";
        for (const auto& control : settings.controlBounds()) {
            control_file << "," << std::get<0>(control);
        }
        control_file << "\n";

        double t0 = settings.t0();
        double tmax = settings.tmax();
        int num_control_intervals = settings.numControlIntervals();
        double control_interval_duration = (tmax - t0) / num_control_intervals;

        for (int i = 0; i < num_control_intervals; ++i) {
            double time = t0 + i * control_interval_duration;
            control_file << std::fixed << std::setprecision(6) << time;
            
            for (int j = 0; j < settings.numControls(); ++j) {
                size_t idx = i * settings.numControls() + j;
                control_file << "," << x[idx];
            }
            control_file << "\n";
        }

        // duplicate the last entry at the end of the file
        double time = tmax;
        control_file << std::fixed << std::setprecision(6) << time;
        for (int j = 0; j < settings.numControls(); ++j) {
            size_t idx = (num_control_intervals - 1) * settings.numControls() + j;
            control_file << "," << x[idx];
        }
        control_file << "\n";

        control_file.close();
    }


    std::vector<FP> timeSteps   = make_time_grid<FP>(settings.t0(), settings.tmax(), num_intervals);
    size_t num_infection_states = static_cast<size_t>(InfectionState::Count);


    auto damping_helper = [=](
        mio::SimulationTime<FP> time, FP min, FP max, 
        mio::DampingLevel damping_level, mio::DampingType damping_type, 
        const std::vector<size_t> location, Eigen::VectorX<FP> group_weights
    ) {
        auto damping_value = mio::UncertainValue<FP>(0.5 * (max + min));
        damping_value.set_distribution(mio::ParameterDistributionUniform(ad::value(min), ad::value(max)));
        return mio::DampingSampling<FP>(damping_value, damping_level, damping_type, time, location, group_weights);
    };

    auto group_weights_all = Eigen::VectorX<FP>::Constant(model.parameters.get_num_groups().get(), 1.0);

    auto set_school_closure = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::Main)),
                                            mio::DampingType(size_t(Intervention::SchoolClosure)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto set_home_office = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::Main)),
                                            mio::DampingType(size_t(Intervention::HomeOffice)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto set_physical_distancing_school = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(size_t(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto set_physical_distancing_work = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(size_t(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto set_physical_distancing_other = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(size_t(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(size_t(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Other)}, group_weights_all);
    };

    mio::UncertainContactMatrix<FP>& contacts = model.parameters.template get<mio::osecirvvs::ContactPatterns<FP>>();
    std::vector<mio::DampingSampling<FP>>& contact_dampings = contacts.get_dampings();

    for (size_t controlIndex = 0; controlIndex < num_control_intervals; controlIndex++) {

        auto param_at = [&](const std::string& name) {
            return x[static_cast<size_t>(string_to_control(name)) + controlIndex * num_controls];
        };
        FP school_closure             = param_at("SchoolClosure");
        FP home_office                = param_at("HomeOffice");
        FP physical_distancing_school = param_at("PhysicalDistancingSchool");
        FP physical_distancing_work   = param_at("PhysicalDistancingWork");
        FP physical_distancing_other  = param_at("PhysicalDistancingOther");

        size_t timeStepIndex = controlIndex * pc_resolution;
        mio::SimulationTime<FP> time(timeSteps[timeStepIndex]);

        contact_dampings.push_back(set_school_closure(time, 1.0 * school_closure));
        contact_dampings.push_back(set_home_office(time, 0.25 * home_office));
        contact_dampings.push_back(set_physical_distancing_school(time, 0.25 * physical_distancing_school));
        contact_dampings.push_back(set_physical_distancing_work(time, 0.25 * physical_distancing_work));
        contact_dampings.push_back(set_physical_distancing_other(time, 0.35 * physical_distancing_other));
    }
    contacts.make_matrix();

    // Open CSV file and write header
    std::ofstream csv_file("population_time_series.csv");
    csv_file << "Time";
    for (size_t i = 0; i < num_infection_states; ++i) {
        csv_file << "," << infection_state_to_string(static_cast<InfectionState>(i));
    }
    csv_file << "\n";

    for (size_t controlIndex = 0; controlIndex < num_control_intervals; controlIndex++) {
        for (size_t substep = 0; substep < pc_resolution; substep++) {
            size_t timeStepIndex = controlIndex * pc_resolution + substep;

            // Simulate next step
            mio::TimeSeries<FP> result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(
                timeSteps[timeStepIndex], timeSteps[timeStepIndex + 1], settings.dt(), model);
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
