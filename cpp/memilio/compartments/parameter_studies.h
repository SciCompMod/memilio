/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#ifndef PARAMETER_STUDIES_H
#define PARAMETER_STUDIES_H

#include "memilio/io/binary_serializer.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/time_series.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"

#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>

#include <type_traits>

// Check if a type has a member called 'stay_duration'
template <typename T, typename = void>
struct has_stay_duration : std::false_type {
};

template <typename T>
struct has_stay_duration<T, std::void_t<decltype(std::declval<T>().stay_duration)>> : std::true_type {
};

// Check if a type has a member called 'travel_time'
template <typename T, typename = void>
struct has_travel_time : std::false_type {
};

template <typename T>
struct has_travel_time<T, std::void_t<decltype(std::declval<T>().travel_time)>> : std::true_type {
};

// Check if a type has a member called 'path'
template <typename T, typename = void>
struct has_path : std::false_type {
};

template <typename T>
struct has_path<T, std::void_t<decltype(std::declval<T>().path)>> : std::true_type {
};

namespace mio
{

/**
 * Class that performs multiple simulation runs with randomly sampled parameters.
 * Can simulate mobility graphs with one simulation in each node or single simulations.
 * @tparam S type of simulation that runs in one node of the graph.
 * @tparam ParametersGraph stores the parameters of the simulation. This is the input of ParameterStudies.
 * @tparam SimulationGraph stores simulations and their results of each run. This is the output of ParameterStudies for each run.
 */
template <class S, class ParametersGraph = Graph<typename S::Model, MobilityParameters<double>>,
          class SimulationGraph = Graph<SimulationNode<S>, MobilityEdge<double>>>
class ParameterStudy
{
public:
    using Simulation = S;

    /**
     * create study for graph of compartment models.
     * @param graph graph of parameters
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param graph_sim_dt time step of graph simulation
     * @param num_runs number of runs
     */
    ParameterStudy(const ParametersGraph& graph, double t0, double tmax, double graph_sim_dt, size_t num_runs)
        : m_graph(graph)
        , m_num_runs(num_runs)
        , m_t0{t0}
        , m_tmax{tmax}
        , m_dt_graph_sim(graph_sim_dt)
    {
    }

    /**
     * create study for graph of compartment models.
     * Creates distributions for all parameters of the models in the graph.
     * @param graph graph of parameters
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param dev_rel relative deviation of the created distributions from the initial value.
     * @param graph_sim_dt time step of graph simulation
     * @param num_runs number of runs
     */
    ParameterStudy(const ParametersGraph& graph, double t0, double tmax, double dev_rel, double graph_sim_dt,
                   size_t num_runs)
        : ParameterStudy(graph, t0, tmax, graph_sim_dt, num_runs)
    {
        for (auto& params_node : m_graph.nodes()) {
            set_params_distributions_normal(params_node, t0, tmax, dev_rel);
        }
    }

    /**
     * @brief Create study for single compartment model.
     * @param model compartment model with initial values
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param num_runs number of runs in ensemble run
     */
    ParameterStudy(typename Simulation::Model const& model, double t0, double tmax, size_t num_runs)
        : ParameterStudy({}, t0, tmax, tmax - t0, num_runs)
    {
        m_graph.add_node(0, model);
    }

    /*
     * @brief Carry out all simulations in the parameter study.
     * Save memory and enable more runs by immediately processing and/or discarding the result.
     * The result processing function is called when a run is finished. It receives the result of the run 
     * (a SimulationGraph object) and an ordered index. The values returned by the result processing function 
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
     * @tparam HandleSimulationResultFunction Callable type, accepts instance of SimulationGraph and an index of type size_t.
     */
    template <class SampleGraphFunction, class HandleSimulationResultFunction>
    std::vector<std::invoke_result_t<HandleSimulationResultFunction, SimulationGraph, size_t>>
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
        m_rng.synchronize();
        thread_local_rng() = m_rng;

        auto run_distribution = distribute_runs(m_num_runs, num_procs);
        auto start_run_idx =
            std::accumulate(run_distribution.begin(), run_distribution.begin() + size_t(rank), size_t(0));
        auto end_run_idx = start_run_idx + run_distribution[size_t(rank)];

        std::vector<std::invoke_result_t<HandleSimulationResultFunction, SimulationGraph, size_t>> ensemble_result;
        ensemble_result.reserve(m_num_runs);

        for (size_t run_idx = start_run_idx; run_idx < end_run_idx; run_idx++) {
            log(LogLevel::info, "ParameterStudies: run {}", run_idx);

            //prepare rng for this run by setting the counter to the right offset
            //Add the old counter so that this call of run() produces different results
            //from the previous call
            auto run_rng_counter = m_rng.get_counter() + rng_totalsequence_counter<uint64_t>(
                                                             static_cast<uint32_t>(run_idx), Counter<uint32_t>(0));
            thread_local_rng().set_counter(run_rng_counter);

            //sample
            auto sim = create_sampled_simulation(sample_graph);
            log(LogLevel::info, "ParameterStudies: Generated {} random numbers.",
                (thread_local_rng().get_counter() - run_rng_counter).get());

            //perform run
            sim.advance(m_tmax);

            //handle result and store
            ensemble_result.emplace_back(result_processing_function(std::move(sim).get_graph(), run_idx));
        }

        //Set the counter of our RNG so that future calls of run() produce different parameters.
        m_rng.set_counter(m_rng.get_counter() + rng_totalsequence_counter<uint64_t>(m_num_runs, Counter<uint32_t>(0)));

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

    /*
     * @brief Carry out all simulations in the parameter study.
     * Convenience function for a few number of runs, but can use more memory because it stores all runs until the end.
     * Unlike the other overload, this function is not MPI-parallel.
     * @return vector of SimulationGraph for each run.
     */
    template <class SampleGraphFunction>
    std::vector<SimulationGraph> run(SampleGraphFunction sample_graph)
    {
        std::vector<SimulationGraph> ensemble_result;
        ensemble_result.reserve(m_num_runs);

        //The ParameterDistributions used for sampling parameters use thread_local_rng()
        //So we set our own RNG to be used.
        //Assume that sampling uses the thread_local_rng() and isn't multithreaded
        thread_local_rng() = m_rng;

        for (size_t i = 0; i < m_num_runs; i++) {
            log(LogLevel::info, "ParameterStudies: run {}", i);

            //prepare rng for this run by setting the counter to the right offset
            //Add the old counter so that this call of run() produces different results
            //from the previous call
            auto run_rng_counter = m_rng.get_counter() +
                                   rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(i), Counter<uint32_t>(0));
            thread_local_rng().set_counter(run_rng_counter);

            auto sim = create_sampled_simulation(sample_graph);
            log(LogLevel::info, "ParameterStudies: Generated {} random numbers.",
                (thread_local_rng().get_counter() - run_rng_counter).get());

            sim.advance(m_tmax);

            ensemble_result.emplace_back(std::move(sim).get_graph());
        }

        //Set the counter of our RNG so that future calls of run() produce different parameters.
        m_rng.set_counter(m_rng.get_counter() + rng_totalsequence_counter<uint64_t>(m_num_runs, Counter<uint32_t>(0)));

        return ensemble_result;
    }

    /*
     * @brief sets the number of Monte Carlo runs
     * @param[in] num_runs number of runs
     */

    void set_num_runs(size_t num_runs)
    {
        m_num_runs = num_runs;
    }

    /*
     * @brief returns the number of Monte Carlo runs
     */
    int get_num_runs() const
    {
        return static_cast<int>(m_num_runs);
    }

    /*
     * @brief sets end point in simulation
     * @param[in] tmax end point in simulation
     */
    void set_tmax(double tmax)
    {
        m_tmax = tmax;
    }

    /*
     * @brief returns end point in simulation
     */
    double get_tmax() const
    {
        return m_tmax;
    }

    void set_t0(double t0)
    {
        m_t0 = t0;
    }

    /*
     * @brief returns start point in simulation
     */
    double get_t0() const
    {
        return m_t0;
    }

    /**
     * Get the input model that the parameter study is run for.
     * Use for single node simulations, use get_model_graph for graph simulations.
     * @{
     */
    const typename Simulation::Model& get_model() const
    {
        return m_graph.nodes()[0].property;
    }
    typename Simulation::Model& get_model()
    {
        return m_graph.nodes()[0].property;
    }
    /** @} */

    /**
     * Get the input graph that the parameter study is run for.
     * Use for graph simulations, use get_model for single node simulations.
     * @{
     */
    const ParametersGraph& get_model_graph() const
    {
        return m_graph;
    }
    ParametersGraph& get_model_graph()
    {
        return m_graph;
    }
    /** @} */

    RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

private:
    template <class SampleGraphFunction>
    auto create_sampled_simulation(SampleGraphFunction sample_graph)
    {
        SimulationGraph sim_graph;

        auto sampled_graph = sample_graph(m_graph);
        for (auto&& node : sampled_graph.nodes()) {
            using PropertyType = typename std::decay<decltype(node.property)>::type;
            add_node_with_properties(sim_graph, node,
                                     std::integral_constant<bool, has_stay_duration<PropertyType>::value>{});
        }
        for (auto&& edge : sampled_graph.edges()) {
            using PropertyType = typename std::decay<decltype(edge.property)>::type;
            add_edge_with_properties(sim_graph, edge, std::integral_constant < bool,
                                     has_travel_time<PropertyType>::value&& has_path<PropertyType>::value > {});
        }

        return make_mobility_sim(m_t0, m_dt_graph_sim, std::move(sim_graph));
    }

    template <typename GraphType, typename NodeType>
    void add_node_with_properties(GraphType& graph, const NodeType& node, std::false_type)
    {
        graph.add_node(node.id, node.property, m_t0, m_dt_integration);
    }

    template <typename GraphType, typename NodeType>
    void add_node_with_properties(GraphType& graph, const NodeType& node, std::true_type)
    {
        graph.add_node(node.id, node.property.base_sim, node.property.mobility_sim, node.property.stay_duration);
    }

    template <typename GraphType, typename EdgeType>
    void add_edge_with_properties(GraphType& graph, const EdgeType& edge, std::false_type)
    {
        graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
    }

    template <typename GraphType, typename EdgeType>
    void add_edge_with_properties(GraphType& graph, const EdgeType& edge, std::true_type)
    {
        graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property.get_parameters(),
                       edge.property.travel_time, edge.property.path);
    }

    std::vector<size_t> distribute_runs(size_t num_runs, int num_procs)
    {
        //evenly distribute runs
        //lower processes do one more run if runs are not evenly distributable
        auto num_runs_local = num_runs / num_procs; //integer division!
        auto remainder      = num_runs % num_procs;

        std::vector<size_t> run_distribution(num_procs);
        std::fill(run_distribution.begin(), run_distribution.begin() + remainder, num_runs_local + 1);
        std::fill(run_distribution.begin() + remainder, run_distribution.end(), num_runs_local);

        return run_distribution;
    }

private:
    // Stores Graph with the names and ranges of all parameters
    ParametersGraph m_graph;

    size_t m_num_runs;

    // Start time (should be the same for all simulations)
    double m_t0;
    // End time (should be the same for all simulations)
    double m_tmax;
    // time step of the graph
    double m_dt_graph_sim;
    // adaptive time step of the integrator (will be corrected if too large/small)
    double m_dt_integration = 0.1;
    //
    RandomNumberGenerator m_rng;
};

} // namespace mio

#endif // PARAMETER_STUDIES_H
