/*
* Copyright (C) 2020-2025 MEmilio
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
#ifndef MIO_COMPARTMENTS_PARAMETER_STUDIES_H
#define MIO_COMPARTMENTS_PARAMETER_STUDIES_H

#include "memilio/io/binary_serializer.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/time_series.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <type_traits>
#include <utility>

namespace mio
{

/**
 * Class that performs multiple simulation runs with randomly sampled parameters.
 * Can simulate mobility graphs with one simulation in each node or single simulations.
 * @tparam S type of simulation that runs in one node of the graph.
 */
template <class SimulationType, class ParameterType, typename TimeType, typename StepType = TimeType>
class ParameterStudy2
{
public:
    using Simulation = SimulationType;
    using Parameters = ParameterType;
    using Time       = TimeType;
    using Step       = StepType;

    // TODO: replacement for "set_params_distributions_normal". Maybe a special ctor for UncertainParameterSet?

    /**
     * create study for graph of compartment models.
     * @param graph graph of parameters
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param graph_sim_dt time step of graph simulation
     * @param num_runs number of runs
     */
    ParameterStudy2(const Parameters& global_parameters, Time t0, Time tmax, Step dt, size_t num_runs)
        : m_parameters(global_parameters)
        , m_num_runs(num_runs)
        , m_t0{t0}
        , m_tmax{tmax}
        , m_dt(dt)
    {
    }

    // /**
    //  * @brief Create study for single compartment model.
    //  * @param model compartment model with initial values
    //  * @param t0 start time of simulations
    //  * @param tmax end time of simulations
    //  * @param num_runs number of runs in ensemble run
    //  */
    // template <class = std::enable_if_t<std::is_same_v<typename Simulation::Model, Parameters>>>
    // ParameterStudy2(typename Simulation::Model const& model, Time t0, Time tmax, Step dt, size_t num_runs)
    //     : ParameterStudy2<Simulation, typename Simulation::Model, Time, Step>(model, t0, tmax, dt, num_runs)
    // {
    //     // TODO how is this supposed to work wrt. model? is this just a special case where ParameterType=Model?
    // }

    /**
     * @brief  
     * @param sample_simulation A function that accepts ParameterType and returns an instance of SimulationType.
     * @param process_simulation_result A function that accepts S 
     */
    template <class CreateSimulationFunction, class ProcessSimulationResultFunction, class... Args>
    std::vector<std::decay_t<std::invoke_result_t<ProcessSimulationResultFunction, Simulation, size_t>>>
    run(CreateSimulationFunction&& create_simulation, ProcessSimulationResultFunction&& process_simulation_result)
    {
        static_assert(std::is_invocable_r_v<Simulation, CreateSimulationFunction, Parameters, Time, Step, size_t>,
                      "Incorrect Type for create_simulation.");
        static_assert(std::is_invocable_v<ProcessSimulationResultFunction, Simulation, size_t>,
                      "Incorrect Type for process_simulation_result.");
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

        std::vector<std::decay_t<std::invoke_result_t<ProcessSimulationResultFunction, Simulation, size_t>>>
            ensemble_result;
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
            Simulation sim =
                create_simulation(std::as_const(m_parameters), std::as_const(m_t0), std::as_const(m_dt), run_idx);
            log(LogLevel::info, "ParameterStudies: Generated {} random numbers.",
                (thread_local_rng().get_counter() - run_rng_counter).get());

            //perform run
            sim.advance(m_tmax);

            //handle result and store
            ensemble_result.emplace_back(process_simulation_result(std::move(sim), run_idx));
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

    /**
     * @brief Carry out all simulations in the parameter study.
     * Convenience function for a few number of runs, but can use more memory because it stores all runs until the end.
     * Unlike the other overload, this function is not MPI-parallel.
     * @return vector of SimulationGraph for each run.
     */
    template <class CreateSimulationFunction>
    std::vector<Simulation> run(CreateSimulationFunction&& create_simulation)
    {
        return run(std::forward<CreateSimulationFunction>(create_simulation), &result_forwarding_function);
    }

    /**
     * @brief returns the number of Monte Carlo runs
     */
    size_t get_num_runs() const
    {
        return m_num_runs;
    }

    /**
     * @brief returns end point in simulation
     */
    Time get_tmax() const
    {
        return m_tmax;
    }

    /**
     * @brief returns start point in simulation
     */
    Time get_t0() const
    {
        return m_t0;
    }
    /**
     * Get the input graph that the parameter study is run for.
     * Use for graph simulations, use get_model for single node simulations.
     * @{
     */
    const Parameters& get_parameters() const
    {
        return m_parameters;
    }
    /** @} */

    RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

private:
    std::vector<size_t> distribute_runs(size_t num_runs, int num_procs)
    {
        assert(num_procs > 0);
        //evenly distribute runs
        //lower processes do one more run if runs are not evenly distributable
        auto num_runs_local = num_runs / num_procs; //integer division!
        auto remainder      = num_runs % num_procs;

        std::vector<size_t> run_distribution(num_procs);
        std::fill(run_distribution.begin(), run_distribution.begin() + remainder, num_runs_local + 1);
        std::fill(run_distribution.begin() + remainder, run_distribution.end(), num_runs_local);

        return run_distribution;
    }

    inline static Simulation&& result_forwarding_function(Simulation&& sim, size_t)
    {
        return std::move(sim);
    }

private:
    // Stores Graph with the names and ranges of all parameters
    ParameterType m_parameters;

    size_t m_num_runs;

    // Start time (should be the same for all simulations)
    Time m_t0;
    // End time (should be the same for all simulations)
    Time m_tmax;
    // adaptive time step of the integrator (will be corrected if too large/small)
    Step m_dt;
    //
    RandomNumberGenerator m_rng;
};

template <class Simulation, class Parameters, typename Time, typename Step = Time>
ParameterStudy2<Simulation, Parameters, Time, Step> make_parameter_study(const Parameters& global_parameters, Time t0,
                                                                         Time tmax, Step dt, size_t num_runs)
{
    return {global_parameters, t0, tmax, dt, num_runs};
}

template <class FP, class Sim>
auto make_parameter_study_graph_ode(const Graph<typename Sim::Model, MobilityParameters<FP>>& global_parameters, FP t0,
                                    FP tmax, FP dt, size_t num_runs)
{
    using SimGraph    = Graph<mio::SimulationNode<ScalarType, Sim>, mio::MobilityEdge<ScalarType>>;
    using SimGraphSim = mio::GraphSimulation<ScalarType, SimGraph, ScalarType, ScalarType>;
    using Params      = Graph<typename Sim::Model, MobilityParameters<FP>>;

    return ParameterStudy2<SimGraphSim, Params, ScalarType>{global_parameters, t0, tmax, dt, num_runs};
}

//sample parameters and create simulation
template <typename FP, class GraphSim>
GraphSim make_sampled_graph_simulation(
    const Graph<typename GraphSim::Graph::NodeProperty::Simulation::Model, MobilityParameters<FP>>& sampled_graph,
    FP t0, FP dt_node_sim, FP dt_graph_sim)
{
    typename GraphSim::Graph sim_graph;

    for (auto&& node : sampled_graph.nodes()) {
        sim_graph.add_node(node.id, node.property, t0, dt_node_sim);
    }
    for (auto&& edge : sampled_graph.edges()) {
        sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
    }

    return make_mobility_sim<FP, typename GraphSim::Graph::NodeProperty::Simulation>(t0, dt_graph_sim,
                                                                                     std::move(sim_graph));
}

} // namespace mio

#endif // MIO_COMPARTMENTS_PARAMETER_STUDIES_H
