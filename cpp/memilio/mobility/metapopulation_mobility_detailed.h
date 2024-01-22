/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Henrik Zunker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef METAPOPULATION_MOBILITY_DETAILED_H
#define METAPOPULATION_MOBILITY_DETAILED_H

#include "memilio/compartments/parameter_studies.h"
#include "memilio/epidemiology/simulation_day.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "memilio/math/eigen_util.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/math/euler.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/dynamic_npis.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/date.h"
#include "memilio/mobility/graph.h"
#include "memilio/io/mobility_io.h"

#include "boost/filesystem.hpp"

#include <cassert>
#include <string>

namespace mio
{

template <class NodePropertyT>
struct NodeDetailed : Node<NodePropertyT> {
    template <class... Args>
    NodeDetailed(int node_id, Args&&... args)
        : Node<NodePropertyT>(node_id, std::forward<Args>(args)...)
        , stay_duration(0.5)
        , mobility(std::forward<Args>(args)...)
    {
    }

    template <typename Model>
    NodeDetailed(int node_id, double duration, Model property_arg, Model mobility_arg, double m_t0,
                 double m_dt_integration)
        : Node<NodePropertyT>(node_id, property_arg, m_t0, m_dt_integration)
        , stay_duration(duration)
        , mobility(mobility_arg, m_t0, m_dt_integration)
    {
    }

    template <class... Args>
    NodeDetailed(int node_id, double duration, Args&&... args)
        : Node<NodePropertyT>(node_id, std::forward<Args>(args)...)
        , stay_duration(duration)
        , mobility(std::forward<Args>(args)...)
    {
    }

    NodeDetailed(int node_id, double duration, NodePropertyT property_arg, NodePropertyT mobility_pt_arg)
        : Node<NodePropertyT>(node_id, property_arg)
        , stay_duration(duration)
        , mobility(mobility_pt_arg)
    {
    }

    double stay_duration;
    NodePropertyT mobility;
};

template <class EdgePropertyT>
struct EdgeDetailed : Edge<EdgePropertyT> {
    template <class... Args>
    EdgeDetailed(size_t start, size_t end, Args&&... args)
        : Edge<EdgePropertyT>(start, end, std::forward<Args>(args)...)
        , traveltime(0.)
        , path{static_cast<int>(start), static_cast<int>(end)}
    {
    }

    template <class... Args>
    EdgeDetailed(size_t start, size_t end, double t_travel, Args&&... args)
        : Edge<EdgePropertyT>(start, end, std::forward<Args>(args)...)
        , traveltime(t_travel)
        , path{static_cast<int>(start), static_cast<int>(end)}
    {
    }

    template <class... Args>
    EdgeDetailed(size_t start, size_t end, double t_travel, std::vector<int> path_mobility, Args&&... args)
        : Edge<EdgePropertyT>(start, end, std::forward<Args>(args)...)
        , traveltime(t_travel)
        , path(path_mobility)
    {
    }
    double traveltime;
    std::vector<int> path;
};

template <class NodePropertyT, class EdgePropertyT>
class GraphDetailed : public Graph<NodePropertyT, EdgePropertyT>
{
public:
    using Graph<NodePropertyT, EdgePropertyT>::Graph;

    /**
     * @brief range of nodes
     */
    auto nodes()
    {
        return make_range(begin(m_nodes), end(m_nodes));
    }

    /**
     * @brief range of nodes
     */
    auto nodes() const
    {
        return make_range(begin(m_nodes), end(m_nodes));
    };

    /**
     * @brief range of edges
     */
    auto edges()
    {
        return make_range(begin(m_edges), end(m_edges));
    }

    /**
     * @brief range of edges
     */
    auto edges() const
    {
        return make_range(begin(m_edges), end(m_edges));
    }

    /**
     * @brief range of edges going out from a specific node
     */
    auto out_edges(size_t node_idx)
    {
        return out_edges(begin(m_edges), end(m_edges), node_idx);
    }

    /**
     * @brief range of edges going out from a specific node
     */
    auto out_edges(size_t node_idx) const
    {
        return out_edges(begin(m_edges), end(m_edges), node_idx);
    }

    template <class ModelType>
    NodeDetailed<NodePropertyT>& add_node(int id, double duration_stay, ModelType& model1, ModelType& model2)
    {
        m_nodes.emplace_back(id, duration_stay, model1, model2);
        return m_nodes.back();
    }

    template <class ModelType>
    NodeDetailed<NodePropertyT>& add_node(int id, double duration_stay, ModelType& model1, ModelType& model2,
                                          double m_t0, double m_dt_integration)
    {
        m_nodes.emplace_back(id, duration_stay, model1, model2, m_t0, m_dt_integration);
        return m_nodes.back();
    }

    EdgeDetailed<EdgePropertyT>& add_edge(size_t start_node_idx, size_t end_node_idx, double traveltime,
                                          mio::MigrationParameters& args)
    {
        assert(m_nodes.size() > start_node_idx && m_nodes.size() > end_node_idx);
        return *insert_sorted_replace(m_edges,
                                      EdgeDetailed<EdgePropertyT>(start_node_idx, end_node_idx, traveltime, args),
                                      [](auto&& e1, auto&& e2) {
                                          return e1.start_node_idx == e2.start_node_idx
                                                     ? e1.end_node_idx < e2.end_node_idx
                                                     : e1.start_node_idx < e2.start_node_idx;
                                      });
    }

    template <class... Args>
    EdgeDetailed<EdgePropertyT>& add_edge(size_t start_node_idx, size_t end_node_idx, double traveltime,
                                          std::vector<int> path, Args&&... args)
    {
        assert(m_nodes.size() > start_node_idx && m_nodes.size() > end_node_idx);
        return *insert_sorted_replace(
            m_edges,
            EdgeDetailed<EdgePropertyT>(start_node_idx, end_node_idx, traveltime, path, std::forward<Args>(args)...),
            [](auto&& e1, auto&& e2) {
                return e1.start_node_idx == e2.start_node_idx ? e1.end_node_idx < e2.end_node_idx
                                                              : e1.start_node_idx < e2.start_node_idx;
            });
    }

private:
    std::vector<NodeDetailed<NodePropertyT>> m_nodes;
    std::vector<EdgeDetailed<EdgePropertyT>> m_edges;
};

/**
 * @brief Sets the graph nodes for counties or districts.
 * Reads the node ids which could refer to districts or counties and the epidemiological
 * data from json files and creates one node for each id. Every node contains a model.
 * @param[in] params Model Parameters that are used for every node.
 * @param[in] start_date Start date for which the data should be read.
 * @param[in] end_data End date for which the data should be read.
 * @param[in] data_dir Directory that contains the data files.
 * @param[in] population_data_path Path to json file containing the population data.
 * @param[in] stay_times_data_path Path to txt file containing the stay times for the considered local entities.
 * @param[in] is_node_for_county Specifies whether the node ids should be county ids (true) or district ids (false).
 * @param[in, out] params_graph Graph whose nodes are set by the function.
 * @param[in] read_func Function that reads input data for german counties and sets Model compartments.
 * @param[in] node_func Function that returns the county ids.
 * @param[in] scaling_factor_inf Factor of confirmed cases to account for undetected cases in each county.
 * @param[in] scaling_factor_icu Factor of ICU cases to account for underreporting.
 * @param[in] tnt_capacity_factor Factor for test and trace capacity.
 * @param[in] num_days Number of days to be simulated; required to load data for vaccinations during the simulation.
 * @param[in] export_time_series If true, reads data for each day of simulation and writes it in the same directory as the input files.
 */
template <class TestAndTrace, class ContactPattern, class Model, class MigrationParams, class Parameters,
          class ReadFunction, class NodeIdFunction>
IOResult<void> set_nodes_detailed(const Parameters& params, Date start_date, Date end_date, const fs::path& data_dir,
                                  const std::string& population_data_path, const std::string& stay_times_data_path,
                                  bool is_node_for_county, GraphDetailed<Model, MigrationParams>& params_graph,
                                  ReadFunction&& read_func, NodeIdFunction&& node_func,
                                  const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                  double tnt_capacity_factor, int num_days = 0, bool export_time_series = false)
{

    BOOST_OUTCOME_TRY(duration_stay, mio::read_duration_stay(stay_times_data_path));
    BOOST_OUTCOME_TRY(node_ids, node_func(population_data_path, is_node_for_county));

    std::vector<Model> nodes(node_ids.size(), Model(int(size_t(params.get_num_groups()))));

    for (auto& node : nodes) {
        node.parameters = params;
    }
    BOOST_OUTCOME_TRY(read_func(nodes, start_date, node_ids, scaling_factor_inf, scaling_factor_icu, data_dir.string(),
                                num_days, export_time_series));

    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {

        auto tnt_capacity = nodes[node_idx].populations.get_total() * tnt_capacity_factor;

        //local parameters
        auto& tnt_value = nodes[node_idx].parameters.template get<TestAndTrace>();
        tnt_value       = mio::UncertainValue(0.5 * (1.2 * tnt_capacity + 0.8 * tnt_capacity));
        tnt_value.set_distribution(mio::ParameterDistributionUniform(0.8 * tnt_capacity, 1.2 * tnt_capacity));

        //holiday periods
        auto id              = int(mio::regions::CountyId(node_ids[node_idx]));
        auto holiday_periods = mio::regions::get_holidays(mio::regions::get_state_id(id), start_date, end_date);
        auto& contacts       = nodes[node_idx].parameters.template get<ContactPattern>();
        contacts.get_school_holidays() =
            std::vector<std::pair<mio::SimulationTime, mio::SimulationTime>>(holiday_periods.size());
        std::transform(
            holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(), [=](auto& period) {
                return std::make_pair(mio::SimulationTime(mio::get_offset_in_days(period.first, start_date)),
                                      mio::SimulationTime(mio::get_offset_in_days(period.second, start_date)));
            });

        //uncertainty in populations
        for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = mio::Index<typename Model::Compartments>(0); j < Model::Compartments::Count; ++j) {
                auto& compartment_value = nodes[node_idx].populations[{i, j}];
                compartment_value =
                    mio::UncertainValue(0.5 * (1.1 * double(compartment_value) + 0.9 * double(compartment_value)));
                compartment_value.set_distribution(mio::ParameterDistributionUniform(0.9 * double(compartment_value),
                                                                                     1.1 * double(compartment_value)));
            }
        }

        // Add mobility node
        auto mobility = nodes[node_idx];
        mobility.populations.set_total(0);

        params_graph.add_node(node_ids[node_idx], duration_stay((Eigen::Index)node_idx), nodes[node_idx], mobility);
    }
    return success();
}

/**
 * @brief Sets the graph edges.
 * Reads the commuting matrices, travel times and paths from data and creates one edge for each pair of nodes.
 * @param[in] travel_times_path Path to txt file containing the travel times between counties.
 * @param[in] mobility_data_path Path to txt file containing the commuting matrices.
 * @param[in] travelpath_path Path to txt file containing the paths between counties.
 * @param[in, out] params_graph Graph whose nodes are set by the function.
 * @param[in] migrating_compartments Compartments that commute.
 * @param[in] contact_locations_size Number of contact locations.
 * @param[in] commuting_weights Vector with a commuting weight for every AgeGroup.
 */
template <class ContactLocation, class Model, class MigrationParams, class MigrationCoefficientGroup,
          class InfectionState>
IOResult<void>
set_edges_detailed(const std::string& travel_times_path, const std::string mobility_data_path,
                   const std::string& travelpath_path, GraphDetailed<Model, MigrationParams>& params_graph,
                   std::initializer_list<InfectionState>& migrating_compartments, size_t contact_locations_size,
                   std::vector<ScalarType> commuting_weights = std::vector<ScalarType>{},
                   ScalarType theshold_edges                 = 4e-5)
{
    BOOST_OUTCOME_TRY(mobility_data_commuter, mio::read_mobility_plain(mobility_data_path));
    BOOST_OUTCOME_TRY(travel_times, mio::read_mobility_plain(travel_times_path));
    BOOST_OUTCOME_TRY(path_mobility, mio::read_path_mobility(travelpath_path));

    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            auto& populations = params_graph.nodes()[county_idx_i].property.populations;

            // mobility coefficients have the same number of components as the contact matrices.
            // so that the same NPIs/dampings can be used for both (e.g. more home office => fewer commuters)
            auto mobility_coeffs = MigrationCoefficientGroup(contact_locations_size, populations.numel());
            auto num_age_groups  = (size_t)params_graph.nodes()[county_idx_i].property.parameters.get_num_groups();
            commuting_weights =
                (commuting_weights.size() == 0 ? std::vector<ScalarType>(num_age_groups, 1.0) : commuting_weights);

            //commuters
            auto working_population = 0.0;
            auto min_commuter_age   = mio::AgeGroup(2);
            auto max_commuter_age   = mio::AgeGroup(4); //this group is partially retired, only partially commutes
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                working_population += populations.get_group_total(age) * commuting_weights[size_t(age)];
            }
            auto commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) / working_population;
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_index = populations.get_flat_index({age, compartment});
                    mobility_coeffs[size_t(ContactLocation::Work)].get_baseline()[coeff_index] =
                        commuter_coeff_ij * commuting_weights[size_t(age)];
                }
            }

            auto path = path_mobility[county_idx_i][county_idx_j];
            if (static_cast<size_t>(path[0]) != county_idx_i ||
                static_cast<size_t>(path[path.size() - 1]) != county_idx_j)
                std::cout << "Wrong Path for edge " << county_idx_i << " " << county_idx_j << "\n";

            //only add edges with mobility above thresholds for performance
            if (commuter_coeff_ij > theshold_edges) {
                params_graph.add_edge(county_idx_i, county_idx_j, travel_times(county_idx_i, county_idx_j),
                                      path_mobility[county_idx_i][county_idx_j], std::move(mobility_coeffs));
            }
        }
    }

    return success();
}

// class MigrationEdgeDetailed : public MigrationEdge
// {
// public:
//     template <class Sim>
//     void apply_migration(double t, double dt, SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to, int mode);
// };

// /**
//  * edge functor for migration simulation.
//  * @see MigrationEdge::apply_migration
//  */
// template <class Sim, class MigrationEdgeDetailed>
// void apply_migration(double t, double dt, MigrationEdgeDetailed& migrationEdgeDetailed, SimulationNode<Sim>& node_from,
//                      SimulationNode<Sim>& node_to, int mode)
// {
//     migrationEdgeDetailed.apply_migration(t, dt, node_from, node_to, mode);
// }

template <class S>
class ParameterStudyDetailed : public ParameterStudy<S>
{
public:
    /**
    * The type of simulation of a single node of the graph.
    */
    using Simulation = S;
    /**
    * The Graph type that stores the parametes of the simulation.
    * This is the input of ParameterStudies.
    */
    using ParametersGraphDetailed = mio::GraphDetailed<typename Simulation::Model, mio::MigrationParameters>;
    /**
    * The Graph type that stores simulations and their results of each run.
    * This is the output of ParameterStudies for each run.
    */
    using SimulationGraphDetailed = mio::GraphDetailed<mio::SimulationNode<Simulation>, mio::MigrationEdge>;

    /**
     * create study for graph of compartment models.
     * @param graph graph of parameters
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param graph_sim_dt time step of graph simulation
     * @param num_runs number of runs
     */
    ParameterStudyDetailed(const ParametersGraphDetailed& graph, double t0, double tmax, double graph_sim_dt,
                           size_t num_runs)
        : ParameterStudy<S>(graph, t0, tmax, graph_sim_dt, num_runs)
        , m_graph(graph)
        , m_num_runs(num_runs)
        , m_t0{t0}
        , m_tmax{tmax}
        , m_dt_graph_sim(graph_sim_dt)
    {
    }

    /*
     * @brief Carry out all simulations in the parameter study.
     * Save memory and enable more runs by immediately processing and/or discarding the result.
     * The result processing function is called when a run is finished. It receives the result of the run 
     * (a SimulationGraphDetailed object) and an ordered index. The values returned by the result processing function 
     * are gathered and returned as a list.
     * This function is parallelized if memilio is configured with MEMILIO_ENABLE_MPI.
     * The MPI processes each contribute a share of the runs. The sample function and result processing function 
     * are called in the same process that performs the run. The results returned by the result processing function are 
     * gathered at the root process and returned as a list by the root in the same order as if the programm 
     * were running sequentially. Processes other than the root return an empty list.
     * @param sample_graph Function that receives the ParametersGraph and returns a sampled copy.
     * @param result_processing_function Processing function for simulation results, e.g., output function.
     * @returns At the root process, a list of values per run that have been returned from the result processing function.
     *          At all other processes, an empty list.
     * @tparam SampleGraphFunction Callable type, accepts instance of ParametersGraph.
     * @tparam HandleSimulationResultFunction Callable type, accepts instance of SimulationGraphDetailed and an index of type size_t.
     */
    template <class SampleGraphFunction, class HandleSimulationResultFunction>
    std::vector<std::invoke_result_t<HandleSimulationResultFunction, SimulationGraphDetailed, size_t>>
    run(SampleGraphFunction sample_graph, HandleSimulationResultFunction result_processing_function)
    {
        int num_procs, rank;
#ifdef MEMILIO_ENABLE_MPI
        MPI_Comm_size(mpi::get_world(), &num_procs);
        MPI_Comm_rank(mpi::get_world(), &rank);
#else
        num_procs = 1;
        rank      = 0;
#endif

        //The ParameterDistributions used for sampling parameters use thread_local_rng()
        //So we set our own RNG to be used.
        //Assume that sampling uses the thread_local_rng() and isn't multithreaded
        this->m_rng.synchronize();
        thread_local_rng() = this->m_rng;

        auto run_distribution = this->distribute_runs(m_num_runs, num_procs);
        auto start_run_idx =
            std::accumulate(run_distribution.begin(), run_distribution.begin() + size_t(rank), size_t(0));
        auto end_run_idx = start_run_idx + run_distribution[size_t(rank)];

        std::vector<std::invoke_result_t<HandleSimulationResultFunction, SimulationGraphDetailed, size_t>>
            ensemble_result;
        ensemble_result.reserve(m_num_runs);

        for (size_t run_idx = start_run_idx; run_idx < end_run_idx; run_idx++) {
            log(LogLevel::info, "ParameterStudies: run {}", run_idx);

            //prepare rng for this run by setting the counter to the right offset
            //Add the old counter so that this call of run() produces different results
            //from the previous call
            auto run_rng_counter =
                this->m_rng.get_counter() +
                rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(run_idx), Counter<uint32_t>(0));
            thread_local_rng().set_counter(run_rng_counter);

            //sample
            auto sim = create_sampled_sim(sample_graph);
            log(LogLevel::info, "ParameterStudies: Generated {} random numbers.",
                (thread_local_rng().get_counter() - run_rng_counter).get());

            //perform run
            sim.advance(m_tmax);

            //handle result and store
            ensemble_result.emplace_back(result_processing_function(std::move(sim).get_graph(), run_idx));
        }

        //Set the counter of our RNG so that future calls of run() produce different parameters.
        this->m_rng.set_counter(this->m_rng.get_counter() +
                                rng_totalsequence_counter<uint64_t>(m_num_runs, Counter<uint32_t>(0)));

#ifdef MEMILIO_ENABLE_MPI
        //gather results
        if (rank == 0) {
            for (int src_rank = 1; src_rank < num_procs; ++src_rank) {
                int bytes_size;
                MPI_Recv(&bytes_size, 1, MPI_INT, src_rank, 0, mpi::get_world(), MPI_STATUS_IGNORE);
                ByteStream bytes(bytes_size);
                MPI_Recv(bytes.data(), bytes.data_size(), MPI_BYTE, src_rank, 0, mpi::get_world(), MPI_STATUS_IGNORE);

                auto src_ensemble_results = deserialize_binary(bytes, Tag<decltype(ensemble_result)>{});
                if (!src_ensemble_results) {
                    log_error("Error receiving ensemble results from rank {}.", src_rank);
                }
                std::copy(src_ensemble_results.value().begin(), src_ensemble_results.value().end(),
                          std::back_inserter(ensemble_result));
            }
        }
        else {
            auto bytes      = serialize_binary(ensemble_result);
            auto bytes_size = int(bytes.data_size());
            MPI_Send(&bytes_size, 1, MPI_INT, 0, 0, mpi::get_world());
            MPI_Send(bytes.data(), bytes.data_size(), MPI_BYTE, 0, 0, mpi::get_world());
            ensemble_result.clear(); //only return root process
        }
#endif

        return ensemble_result;
    }

private:
    //sample parameters and create simulation
    template <class SampleGraphFunction>
    mio::GraphSimulationDetailed<SimulationGraphDetailed> create_sampled_sim(SampleGraphFunction sample_graph)
    {
        SimulationGraphDetailed sim_graph;

        auto sampled_graph = sample_graph(m_graph);
        for (auto&& node : sampled_graph.nodes()) {
            sim_graph.add_node(node.id, node.stay_duration, node.property, node.mobility, m_t0, this->m_dt_integration);
        }
        for (auto&& edge : sampled_graph.edges()) {
            sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.traveltime, edge.property);
        }
        return make_migration_sim(m_t0, m_dt_graph_sim, std::move(sim_graph));
    }

private:
    // Stores Graph with the names and ranges of all parameters
    ParametersGraphDetailed m_graph;
    size_t m_num_runs;
    double m_t0;
    double m_tmax;
    double m_dt_graph_sim;
};

/**
 * @brief number of migrated people when they return according to the model.
 * E.g. during the time in the other node, some people who left as susceptible will return exposed.
 * Implemented for general compartmentmodel simulations, overload for your custom model if necessary
 * so that it can be found with argument-dependent lookup, i.e. in the same namespace as the model.
 * @param migrated number of people that migrated as input, number of people that return as output
 * @param sim Simulation that is used for the migration
 * @param integrator Integrator that is used for the estimation. Has to be a one-step scheme.
 * @param total total population in the node that the people migrated to.
 * @param t time of migration
 * @param dt timestep
 */
template <class Sim, class = std::enable_if_t<is_compartment_model_simulation<Sim>::value>>
void update_status_migrated(Eigen::Ref<TimeSeries<double>::Vector> migrated, const Sim& sim, IntegratorCore& integrator,
                            Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt)
{
    auto y0 = migrated.eval();
    auto y1 = migrated;
    integrator.step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            sim.get_model().get_derivatives(total, y, t_, dydt);
        },
        y0, t, dt, y1);
}

template <typename FP>
using Vector = Eigen::Matrix<FP, Eigen::Dynamic, 1>;

/*
    * move migrated people from one node to another.
    * @param migrated number of people that migrated 
    * @param results_from current state of node that people migrated from
    * @param results_to current state of node that people migrated to
*/
template <typename FP>
void move_migrated(Eigen::Ref<Vector<FP>> migrated, Eigen::Ref<Vector<FP>> results_from,
                   Eigen::Ref<Vector<FP>> results_to)
{
    // The performance of the low order integrator often fails for small compartments i.e. compartments are sometimes negative.
    // We handel this by checking for values below 0 in m_migrated and add them to the highest value
    Eigen::Index max_index_migrated;
    for (size_t i = 0; i < migrated.size(); ++i) {
        migrated.maxCoeff(&max_index_migrated);
        if (migrated(i) < 0) {
            migrated(max_index_migrated) -= migrated(i);
            migrated(i) = 0;
        }
        if (results_from(i) < migrated(i)) {
            migrated(max_index_migrated) -= (migrated(i) - results_from(i));
            migrated(i) = results_from(i);
        }
        results_from(i) -= migrated(i);
        results_to(i) += migrated(i);
    }
}

// template <class Sim>
// void MigrationEdge::apply_migration(double t, double dt, SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to,
//                                     int mode) //
// {

//     if (mode == 0) {
//         //normal daily migration
//         m_migrated.add_time_point(
//             t, (node_from.get_last_state().array() * m_parameters.get_coefficients().get_matrix_at(t).array() *
//                 get_migration_factors(node_from, t, node_from.get_last_state()).array())
//                    .matrix());
//         m_return_times.add_time_point(t + dt);
//         move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
//                       node_to.get_result().get_last_value());
//     }
//     // change county of migrated
//     else if (mode == 1) {
//         // update status of migrated before moving to next county
//         IntegratorCore& integrator_node = node_from.get_simulation().get_integrator();
//         update_status_migrated(m_migrated.get_last_value(), node_from.get_simulation(), integrator_node,
//                                node_from.get_result().get_last_value(), t, dt);
//         move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
//                       node_to.get_result().get_last_value());
//     }
//     // option for last time point to remove time points
//     else if (mode == 2) {
//         Eigen::Index idx                = m_return_times.get_num_time_points() - 1;
//         IntegratorCore& integrator_node = node_from.get_simulation().get_integrator();
//         update_status_migrated(m_migrated[idx], node_from.get_simulation(), integrator_node,
//                                node_from.get_result().get_last_value(), t, dt);

//         move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
//                       node_to.get_result().get_last_value());

//         for (Eigen::Index i = m_return_times.get_num_time_points() - 1; i >= 0; --i) {
//             if (m_return_times.get_time(i) <= t) {
//                 m_migrated.remove_time_point(i);
//                 m_return_times.remove_time_point(i);
//             }
//         }
//     }
//     // just update status of migrated
//     else if (mode == 3) {
//         Eigen::Index idx                = m_return_times.get_num_time_points() - 1;
//         IntegratorCore& integrator_node = node_from.get_simulation().get_integrator();
//         update_status_migrated(m_migrated[idx], node_from.get_simulation(), integrator_node,
//                                node_from.get_result().get_last_value(), t, dt);

//         move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
//                       node_from.get_result().get_last_value());
//     }
//     else {
//         std::cout << "Invalid input mode. Should be 0 or 1."
//                   << "\n";
//     }
// }

/**
 * create a migration simulation.
 * After every second time step, for each edge a portion of the population corresponding to the coefficients of the edge
 * moves from one node to the other. In the next timestep, the migrated population return to their "home" node. 
 * Returns are adjusted based on the development in the target node. 
 * @param t0 start time of the simulation
 * @param dt time step between migrations
 * @param graph set up for migration simulation
 * @{
 */
template <class Sim>
GraphSimulationDetailed<GraphDetailed<SimulationNode<Sim>, MigrationEdge>>
make_migration_sim(double t0, double dt, const GraphDetailed<SimulationNode<Sim>, MigrationEdge>& graph)
{
    return make_graph_sim_detailed(t0, dt, graph, &evolve_model<Sim>, &apply_migration<Sim, MigrationEdge>);
}

template <class Sim>
GraphSimulationDetailed<GraphDetailed<SimulationNode<Sim>, MigrationEdge>>
make_migration_sim(double t0, double dt, GraphDetailed<SimulationNode<Sim>, MigrationEdge>&& graph)
{
    return make_graph_sim_detailed(t0, dt, std::move(graph), &evolve_model<Sim>, &apply_migration<Sim>);
}

} // namespace mio

#endif //METAPOPULATION_MOBILITY_DETAILED_H
