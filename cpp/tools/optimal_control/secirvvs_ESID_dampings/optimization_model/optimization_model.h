#pragma once

#include "boost/filesystem.hpp"

#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "models/ode_secirvvs/model.h"

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <any>
#include <typeindex>
#include <mutex>

class OptimizationModel
{
public:
    OptimizationModel(const boost::filesystem::path& data_directory, size_t simulation_days, size_t num_age_groups);

    boost::filesystem::path data_directory() const;
    size_t simulation_days() const;
    size_t num_age_groups() const;

    // Retrieves (or generates) a graph model of the simulation.
    // Uses template-based caching since reading in models can be time-consuming.
    template <typename FP>
    mio::Graph<mio::SimulationNode<FP, mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>> get_graph_model() const;

private:
    boost::filesystem::path m_data_directory;
    size_t m_simulation_days;
    size_t m_num_age_groups;

    // Creates an artificial (synthetic) graph model with preset parameters and population data.
    template <typename FP>
    mio::Graph<mio::SimulationNode<FP, mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>>
    create_artificial_graph_model() const;

    struct Cache {
        std::unordered_map<std::type_index, std::any> data;
        std::mutex mutex;
    };
    std::shared_ptr<Cache> m_cache;
};

// Retrieves a graph model from the cache or builds a new one if necessary.
template <typename FP>
auto OptimizationModel::get_graph_model() const
    -> mio::Graph<mio::SimulationNode<FP, mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>>
{
    const auto key = std::type_index(typeid(FP));
    {
        // Check if cached
        std::scoped_lock lock(m_cache->mutex);
        auto it = m_cache->data.find(key);
        if (it != m_cache->data.end()) {
            return std::any_cast<
                mio::Graph<mio::SimulationNode<FP, mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>>>(it->second);
        }
    }
    auto graph = create_artificial_graph_model<FP>();
    {
        // Cache the result
        std::scoped_lock lock(m_cache->mutex);
        m_cache->data.emplace(key, graph);
    }
    return graph;
}

// ------------------------------------------------------------ //
// Creates a synthetic graph model for testing the optimization //
template <typename FP>
mio::Graph<mio::SimulationNode<FP, mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>>
OptimizationModel::create_artificial_graph_model() const
{
    // --------------------------- //
    // Define the model parameters //
    constexpr size_t num_age_groups = 6;
    mio::osecirvvs::Model<FP> model(num_age_groups);
    auto& params = model.parameters;

    params.template get<mio::osecirvvs::ICUCapacity<FP>>()                    = 100;
    params.template get<mio::osecirvvs::TestAndTraceCapacity<FP>>()           = 0.0143;
    params.template get<mio::osecirvvs::Seasonality<FP>>()                    = 0.2;
    params.template get<mio::osecirvvs::DynamicNPIsImplementationDelay<FP>>() = 7;

    params.template get<mio::osecirvvs::DailyPartialVaccinations<FP>>().resize(
        mio::SimulationDay(m_simulation_days + 1));
    params.template get<mio::osecirvvs::DailyFullVaccinations<FP>>().resize(mio::SimulationDay(m_simulation_days + 1));

    size_t daily_vaccinations = 10;
    for (size_t age_group = 0; age_group < num_age_groups; age_group++) {
        for (size_t current_day = 0; current_day <= m_simulation_days; current_day++) {
            mio::Index<mio::AgeGroup, mio::SimulationDay> index{mio::AgeGroup(age_group),
                                                                mio::SimulationDay(current_day)};
            FP num_vaccinations = static_cast<FP>(current_day * daily_vaccinations);
            params.template get<mio::osecirvvs::DailyPartialVaccinations<FP>>()[index] = num_vaccinations;
            params.template get<mio::osecirvvs::DailyFullVaccinations<FP>>()[index]    = num_vaccinations;
        }
    }

    auto set_all_groups = [&](auto Tag, FP value) {
        auto& age_group_params = params.template get<decltype(Tag)>();
        for (size_t age_group = 0; age_group < num_age_groups; age_group++) {
            age_group_params[mio::AgeGroup(age_group)] = value;
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

    // --------------------------- //
    // Define the contact matrices //
    int num_contact_locations = 4;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline_home;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline_school_pf_eig;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline_work;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline_other;

    // clang-format off
    baseline_home << 
        0.4413, 0.4504, 1.2383, 0.8033, 0.0494, 0.0017, 
        0.0485, 0.7616, 0.6532, 1.1614, 0.0256, 0.0013,
        0.1800, 0.1795, 0.8806, 0.6413, 0.0429, 0.0032, 
        0.0495, 0.2639, 0.5189, 0.8277, 0.0679, 0.0014, 
        0.0087, 0.0394, 0.1417, 0.3834, 0.7064, 0.0447, 
        0.0292, 0.0648, 0.1248, 0.4179, 0.3497, 0.1544;

    baseline_school_pf_eig << 
        2.9964, 0.2501, 0.9132, 0.2509, 0.0000, 0.0000, 
        0.2210, 1.9155, 0.2574, 0.2863, 0.0070, 0.0000, 
        0.0324, 0.3598, 1.2613, 0.1854, 0.0041, 0.0000, 
        0.1043, 0.4794, 1.1886, 0.1836, 0.0052, 0.0022, 
        0.0000, 0.1150, 0.0000, 0.0000, 0.0168, 0.0000, 
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;

    baseline_other << 
        0.5170, 0.3997, 0.7957, 0.9958, 0.3239, 0.0428, 
        0.0632, 0.9121, 0.3254, 0.4731, 0.2355, 0.0148,
        0.0336, 0.1604, 1.7529, 0.8622, 0.1440, 0.0077, 
        0.0204, 0.1444, 0.5738, 1.2127, 0.3433, 0.0178,
        0.0371, 0.0393, 0.4171, 0.9666, 0.7495, 0.0257, 
        0.0791, 0.0800, 0.3480, 0.5588, 0.2769, 0.0180;

    baseline_work << 
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
        0.0000, 0.0127, 1.7570, 1.6050, 0.0133, 0.0000, 
        0.0000, 0.0020, 1.0311, 2.3166, 0.0098, 0.0000, 
        0.0000, 0.0002, 0.0194, 0.0325, 0.0003, 0.0000, 
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;
    // clang-format on

    Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum_home;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum_school_pf_eig;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum_work;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum_other;

    minimum_home.setZero();
    minimum_school_pf_eig.setZero();
    minimum_work.setZero();
    minimum_other.setZero();

    mio::ContactMatrixGroup<FP> contact_matrices(num_contact_locations, num_age_groups);

    contact_matrices[0].get_baseline() = baseline_home;
    contact_matrices[1].get_baseline() = baseline_school_pf_eig;
    contact_matrices[2].get_baseline() = baseline_work;
    contact_matrices[3].get_baseline() = baseline_other;

    contact_matrices[0].get_minimum() = minimum_home;
    contact_matrices[1].get_minimum() = minimum_school_pf_eig;
    contact_matrices[2].get_minimum() = minimum_work;
    contact_matrices[3].get_minimum() = minimum_other;

    params.template get<mio::osecirvvs::ContactPatterns<FP>>() =
        mio::UncertainContactMatrix<FP>(std::move(contact_matrices));

    // ----------------------------- //
    // Define the initial population //
    // clang-format off
    size_t num_infection_states = static_cast<size_t>(mio::osecirvvs::InfectionState::Count);
    Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic> population(num_age_groups, num_infection_states);
    population <<
        /* ------------- */
        /* AgeGroup: 0-4 */
        2753030.497430,
        0.0, // Susceptible
        11.708820, 0.0, 2.212815, // Exposed
        7.578866, 0.0, 4.084619, 0.0, 0.0, 0.0, // Infected No Symptoms
        22.853384, 0.0, 2.146616, 0.0, 0.0, 0.0, // Infected Symptoms
        0.150378, 0.0, 0.001865, // Infected Severe
        0.150441, 0.0, 0.005732, // Infected Critical
        1051258.609035, // Recovered
        64.000000, 0.0, 0.0, // Dead
        /* ------------- */
        /* AgeGroup: 5-14 */
        3326384.779416, 0.0, // Susceptible
        5.015423, 0.0, 3.287046, // Exposed
        2.945401, 0.0, 5.907620, 0.0, 0.0, 0.0, // Infected No Symptoms
        5.168372, 0.0, 1.831628, 0.0, 0.0, 0.0, // Infected Symptoms
        0.023957, 0.0, 0.001460, // Infected Severe
        0.138450, 0.0, 0.017723, // Infected Critical
        4506996.883505, // Recovered
        45.000000, 0.0, 0.0, // Dead
        /* ------------- */
        /* AgeGroup: 15-34 */
        7474316.734280, 0.0, // Susceptible
        15.185091, 0.0, 9.963930, // Exposed
        10.119921, 0.0, 19.347768, 0.0, 0.0, 0.0, // Infected No Symptoms
        27.662806, 0.0, 9.908623, 0.0, 0.0, 0.0, // Infected Symptoms
        0.839570, 0.0, 0.044404, // Infected Severe
        0.349355, 0.0, 0.046282, // Infected Critical
        11247005.297970, // Recovered
        510.0, 0.0, 0.0, // Dead
        /* ------------- */
        /* AgeGroup: 35-59 */
        13079641.709818, 0.0, // Susceptible
        28.087697, 0.0, 14.523549, // Exposed
        23.028794, 0.0, 36.771572, 0.0, 0.0, 0.0, // Infected No Symptoms
        74.950401, 0.0, 21.335313, 0.0, 0.0, 0.0, // Infected Symptoms
        7.139552, 0.0, 0.287452, // Infected Severe
        2.329300, 0.0, 0.231930, // Infected Critical
        15095393.104624, // Recovered
        8936.000000, 0.0, 0.0, // Dead
        /* ------------- */
        /* AgeGroup: 60-79 */
        12370660.206757, 0.0, // Susceptible
        66.165784, 0.0, 13.006396, // Exposed
        51.098574, 0.0, 30.752489, 0.0, 0.0, 0.0, // Infected No Symptoms
        162.557513, 0.0, 16.728201, 0.0, 0.0, 0.0, // Infected Symptoms
        38.893189, 0.0, 0.564814, // Infected Severe
        13.248394, 0.0, 0.494789, // Infected Critical
        5239420.683100, // Recovered
        56248.571429, 0.0, 0.0, // Dead
        /* ------------- */
        /* AgeGroup: 80+ */
        5614932.848062, 0.0, // Susceptible
        61.140294, 0.0, 9.149288, // Exposed
        46.298860, 0.0, 21.109925, 0.0, 0.0, 0.0, // Infected No Symptoms
        137.630764, 0.0, 11.369236, 0.0, 0.0, 0.0, // Infected Symptoms
        59.433799, 0.0, 0.673454, // Infected Severe
        24.245028, 0.0, 0.742577, // Infected Critical
        1695860.958713, // Recovered
        121825.428571, 0.0, 0.0; // Dead
    // clang-format on

    for (size_t column = 0; column < static_cast<size_t>(population.cols()); column++) {
        for (size_t row = 0; row < static_cast<size_t>(population.rows()); row++) {
            mio::Index<mio::AgeGroup, mio::osecirvvs::InfectionState> index{mio::AgeGroup(row),
                                                                            mio::osecirvvs::InfectionState(column)};
            model.populations[index] = population(row, column);
        }
    }

    model.apply_constraints();

    // ---------------------- //
    // Define the graph model //
    mio::Graph<mio::SimulationNode<FP, mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>> graph_simulation;

    graph_simulation.add_node(0, model);

    return graph_simulation;
}
