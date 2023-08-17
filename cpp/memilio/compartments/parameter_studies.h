/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include <iterator>
#include <numeric>

namespace mio
{

/**
 * Class that performs multiple simulation runs with randomly sampled parameters.
 * Can simulate migration graphs with one simulation in each node or single simulations.
 * @tparam S type of simulation that runs in one node of the graph.
 */
template <class S>
class ParameterStudy
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
    using ParametersGraph = mio::Graph<typename Simulation::Model, mio::MigrationParameters>;
    /**
    * The Graph type that stores simulations and their results of each run.
    * This is the output of ParameterStudies for each run.
    */
    using SimulationGraph = mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>;

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

        auto run_distribution = distribute_runs(m_num_runs, num_procs);
        auto start_run_idx = std::accumulate(run_distribution.begin(), run_distribution.begin() + size_t(rank), size_t(0));
        auto end_run_idx = start_run_idx + run_distribution[size_t(rank)];

        std::vector<std::invoke_result_t<HandleSimulationResultFunction, SimulationGraph, size_t>> ensemble_result;
        ensemble_result.reserve(m_num_runs);

        for (size_t run_idx = start_run_idx; run_idx < end_run_idx; run_idx++) {
            //ensure reproducible results if the number of runs or MPI process changes.
            //assumptions:
            //- the sample_graph functor uses the RNG provided by our library
            //- the RNG has been freshly seeded/initialized before this call
            //- the seeds are identical on all MPI processes
            //- the block size of the RNG is sufficiently big to cover one run
            //  (when in doubt, use a larger block size; fast-forwarding the RNG is cheap and the period length 
            //   of the mt19937 RNG is huge)
            mio::thread_local_rng().forward_to_block(run_idx);

            auto sim = create_sampled_simulation(sample_graph);
            sim.advance(m_tmax);

            ensemble_result.emplace_back(result_processing_function(std::move(sim).get_graph(), run_idx));
        }

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

        for (size_t i = 0; i < m_num_runs; i++) {
            auto sim = create_sampled_simulation(sample_graph);
            sim.advance(m_tmax);

            ensemble_result.emplace_back(std::move(sim).get_graph());
        }

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

private:
    //sample parameters and create simulation
    template <class SampleGraphFunction>
    mio::GraphSimulation<SimulationGraph> create_sampled_simulation(SampleGraphFunction sample_graph)
    {
        SimulationGraph sim_graph;

        auto sampled_graph = sample_graph(m_graph);
        for (auto&& node : sampled_graph.nodes()) {
            sim_graph.add_node(node.id, node.property, m_t0, m_dt_integration);
        }
        for (auto&& edge : sampled_graph.edges()) {
            sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
        }

        return make_migration_sim(m_t0, m_dt_graph_sim, std::move(sim_graph));
    }

    std::vector<size_t> distribute_runs(size_t num_runs, int num_procs)
    {
        //evenly distribute runs
        //lower processes do one more run if runs are not evenly distributable
        auto num_runs_local = num_runs / num_procs; //integer division!
        auto remainder = num_runs % num_procs;

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
};

} // namespace mio

#endif // PARAMETER_STUDIES_H
