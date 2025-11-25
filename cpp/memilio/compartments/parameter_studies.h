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
#include "memilio/io/io.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"

#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

namespace mio
{

/**
 * @brief Class used to perform multiple simulation runs with randomly sampled parameters.
 * Note that the type of simulation is not determined until calling one of the run functions.
 * @tparam ParameterType The parameters used to create simulations.
 * @tparam TimeType The time type used by the simulation, e.g. double or TimePoint.
 * @tparam StepType The time step type used by the simulation, e.g. double or TimeStep. May be the same as TimeType.
 */
template <class ParameterType, typename TimeType, typename StepType = TimeType>
class ParameterStudy
{
public:
    using Parameters = ParameterType;
    using Time       = TimeType;
    using Step       = StepType;

private:
    /// @brief The return type of `create_simulation`. Ensures that the function is invocable.
    template <class CreateSimulationFunction>
        requires std::is_invocable_v<CreateSimulationFunction, Parameters, Time, Step, size_t>
    using SimulationT = std::decay_t<std::invoke_result_t<CreateSimulationFunction, Parameters, Time, Step, size_t>>;
    /// @brief The return type of `process_simulation_result`. Ensures that the function is invocable.
    template <class CreateSimulationFunction, class ProcessSimulationResultFunction>
        requires std::is_invocable_v<ProcessSimulationResultFunction, SimulationT<CreateSimulationFunction>, size_t>
    using ProcessedResultT = std::decay_t<
        std::invoke_result_t<ProcessSimulationResultFunction, SimulationT<CreateSimulationFunction>, size_t>>;
    /// @brief Type returned by run functions. Is void if ProcessedResultT is, otherwise a vector of ProcessedResultT.
    template <class CreateSimulationFunction, class ProcessSimulationResultFunction>
    using EnsembleResultT =
        std::conditional_t<std::is_void_v<ProcessedResultT<CreateSimulationFunction, ProcessSimulationResultFunction>>,
                           void,
                           std::vector<ProcessedResultT<CreateSimulationFunction, ProcessSimulationResultFunction>>>;

public:
    // TODO: replacement for "set_params_distributions_normal"? Maybe a special ctor for UncertainParameterSet?

    /**
     * @brief Create a parameter study with some parameters.
     * The simulation type is determined when calling any "run" member function.
     * @param parameters The parameters used to create simulations.
     * @param t0 Start time of simulations.
     * @param tmax End time of simulations.
     * @param dt Initial time step of simulations.
     * @param num_runs Number of simulations that will be created and run.
     */
    ParameterStudy(const Parameters& parameters, Time t0, Time tmax, Step dt, size_t num_runs)
        : m_parameters(parameters)
        , m_num_runs(num_runs)
        , m_t0{t0}
        , m_tmax{tmax}
        , m_dt(dt)
    {
    }

    /**
     * @brief Run all simulations in serial.
     * @param[in] create_simulation A callable sampling the study's parameters and returning a simulation.
     * @param[in] process_simulation_result (Optional) A callable that takes the simulation and processes its result.
     * @return A vector containing (processed) simulation results for each run, or void if processing returns nothing.
     *
     * Important side effect: Calling this function overwrites seed and counter of thread_local_rng().
     * Use this RNG when sampling parameters in create_simulation.
     *
     * The function signature for create_simulation is
     * `SimulationT(const Parameters& study_parameters, Time t0, Step dt, size_t run_idx)`,
     * where SimulationT is some kind of simulation.
     * The function signature for process_simulation_result is
     * `ProcessedResultT(SimulationT&&, size_t run_index)`,
     * where ProcessedResultT is a (de)serializable result, or void. A void function can be useful if the results
     * should be fully handled during the study, for example, when memory is limited and the results have to be written
     * to a disk.
     * @{
     */
    template <class CreateSimulationFunction, class ProcessSimulationResultFunction>
    EnsembleResultT<CreateSimulationFunction, ProcessSimulationResultFunction>
    run_serial(CreateSimulationFunction&& create_simulation,
               ProcessSimulationResultFunction&& process_simulation_result)
    {
        return run_impl(0, m_num_runs, std::forward<CreateSimulationFunction>(create_simulation),
                        std::forward<ProcessSimulationResultFunction>(process_simulation_result));
    }

    template <class CreateSimulationFunction>
    std::vector<SimulationT<CreateSimulationFunction>> run_serial(CreateSimulationFunction&& create_simulation)
    {
        return run_serial(
            std::forward<CreateSimulationFunction>(create_simulation),
            [](SimulationT<CreateSimulationFunction>&& sim, size_t) -> SimulationT<CreateSimulationFunction>&& {
                return std::move(sim);
            });
    }
    /** @} */

    /**
     * @brief Run all simulations distributed over multiple MPI ranks.
     * @param[in] create_simulation A callable sampling the study's parameters and returning a simulation.
     * @param[in] process_simulation_result A callable that takes the simulation and processes its result.
     * @return A vector that contains processed simulation results for each run, or void if processing returns nothing.
     *
     *
     * Important: Do not forget to use mio::mpi::init and finalize when using this function!
     *
     * Important side effect: Calling this function overwrites seed and counter of thread_local_rng().
     * Use this RNG when sampling parameters in create_simulation.
     *
     * The function signature for create_simulation is
     * `SimulationT(const Parameters& study_parameters, Time t0, Step dt, size_t run_idx)`,
     * where SimulationT is some kind of simulation.
     * The function signature for process_simulation_result is
     * `ProcessedResultT(SimulationT&&, size_t run_index)`,
     * where ProcessedResultT is a (de)serializable result, or void. A void function can be useful if the results
     * should be fully handled during the study, for example, when memory is limited and the results have to be written
     * to a disk.
     * 
     * If MPI is enabled and the results are non-void, all results are gathered on the root rank 0. Other ranks will
     * return an empty vector.
     */
    template <class CreateSimulationFunction, class ProcessSimulationResultFunction>
    EnsembleResultT<CreateSimulationFunction, ProcessSimulationResultFunction>
    run(CreateSimulationFunction&& create_simulation, ProcessSimulationResultFunction&& process_simulation_result)
    {
        using ResultT = EnsembleResultT<CreateSimulationFunction, ProcessSimulationResultFunction>;
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

        std::vector<size_t> run_distribution = distribute_runs(m_num_runs, num_procs);
        size_t start_run_idx =
            std::accumulate(run_distribution.begin(), run_distribution.begin() + size_t(rank), size_t(0));
        size_t end_run_idx = start_run_idx + run_distribution[size_t(rank)];

        if constexpr (std::is_void_v<ResultT>) {
            // if the processor returns nothing, there is nothing to synchronize
            run_impl(start_run_idx, end_run_idx, std::forward<CreateSimulationFunction>(create_simulation),
                     std::forward<ProcessSimulationResultFunction>(process_simulation_result));
            return;
        }
        else {
            auto ensemble_result =
                run_impl(start_run_idx, end_run_idx, std::forward<CreateSimulationFunction>(create_simulation),
                         std::forward<ProcessSimulationResultFunction>(process_simulation_result));

#ifdef MEMILIO_ENABLE_MPI
            //gather results
            if (rank == 0) {
                for (int src_rank = 1; src_rank < num_procs; ++src_rank) {
                    int bytes_size;
                    MPI_Recv(&bytes_size, 1, MPI_INT, src_rank, 0, mpi::get_world(), MPI_STATUS_IGNORE);
                    ByteStream bytes(bytes_size);
                    MPI_Recv(bytes.data(), bytes.data_size(), MPI_BYTE, src_rank, 0, mpi::get_world(),
                             MPI_STATUS_IGNORE);

                    IOResult<ResultT> src_ensemble_results = deserialize_binary(bytes, Tag<ResultT>{});
                    if (!src_ensemble_results) {
                        log_error("Error receiving ensemble results from rank {}.", src_rank);
                    }
                    std::copy(src_ensemble_results.value().begin(), src_ensemble_results.value().end(),
                              std::back_inserter(ensemble_result));
                }
            }
            else {
                ByteStream bytes = serialize_binary(ensemble_result);
                int bytes_size   = int(bytes.data_size());
                MPI_Send(&bytes_size, 1, MPI_INT, 0, 0, mpi::get_world());
                MPI_Send(bytes.data(), bytes.data_size(), MPI_BYTE, 0, 0, mpi::get_world());
                ensemble_result.clear(); //only return root process
            }
#endif

            return ensemble_result;
        }
    }

    /// @brief Return the number of total runs that the study will make.
    size_t get_num_runs() const
    {
        return m_num_runs;
    }

    /// @brief Return the final time point for simulations.
    Time get_tmax() const
    {
        return m_tmax;
    }

    /// @brief Return the initial time point for simulations.
    Time get_t0() const
    {
        return m_t0;
    }

    /// @brief Return the initial step sized used by simulations.
    Time get_dt() const
    {
        return m_dt;
    }

    /**
     * @brief Get the input parameters that each simulation in the study is created from.
     * @{
     */
    const Parameters& get_parameters() const
    {
        return m_parameters;
    }
    Parameters& get_parameters()
    {
        return m_parameters;
    }
    /** @} */

    /// @brief Access the study's random number generator.
    RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

private:
    /// @brief Return the ensemble result vector, or "void" in form of a char.
    template <class T>
    inline auto make_ensemble_result()
    {
        if constexpr (std::is_void_v<T>) {
            return char(0); // placeholder, as we cannot instanciate void
        }
        else {
            T result;
            result.reserve(m_num_runs);
            return result;
        }
    }

    /**
     * @brief Main loop creating and running simulations.
     * @param[in] start_run_idx, end_run_idx Range of indices. Performs one run for each index.
     * @param[in] create_simulation A callable sampling the study's parameters and return a simulation.
     * @param[in] process_simulation_result A callable that takes the simulation and processes its result.
     * @return A vector that contains processed simulation results for each run.
     *
     * Important side effect: Calling this function overwrites seed and counter of thread_local_rng().
     * Use this RNG when sampling parameters in create_simulation.
     */
    template <class CreateSimulationFunction, class ProcessSimulationResultFunction>
    EnsembleResultT<CreateSimulationFunction, ProcessSimulationResultFunction>
    run_impl(size_t start_run_idx, size_t end_run_idx, CreateSimulationFunction&& create_simulation,
             ProcessSimulationResultFunction&& process_simulation_result)
    {
        using ResultT = EnsembleResultT<CreateSimulationFunction, ProcessSimulationResultFunction>;
        assert(start_run_idx <= end_run_idx);
        // this bool (and all code that it is used in) enables using void functions
        constexpr bool should_gather_results = !std::is_void_v<ResultT>;

        // Note that this overwrites seed and counter of thread_local_rng, but it does not replace it.
        thread_local_rng() = m_rng;

        // the result is not used or returned if the processing function returns nothing (i.e. void)
        [[maybe_unused]] auto ensemble_result = make_ensemble_result<ResultT>();

        for (size_t run_idx = start_run_idx; run_idx < end_run_idx; run_idx++) {
            log(LogLevel::info, "ParameterStudies: run {}", run_idx);

            //prepare rng for this run by setting the counter to the right offset
            //Add the old counter so that this call of run() produces different results
            //from the previous call
            Counter<uint64_t> run_rng_counter =
                m_rng.get_counter() +
                rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(run_idx), Counter<uint32_t>(0));
            thread_local_rng().set_counter(run_rng_counter);

            //sample
            SimulationT<CreateSimulationFunction> sim =
                create_simulation(std::as_const(m_parameters), std::as_const(m_t0), std::as_const(m_dt), run_idx);

            // the create_counter is only used for debug logging, hence it is unused in Release builds
            [[maybe_unused]] const uint64_t create_counter = (thread_local_rng().get_counter() - run_rng_counter).get();
            log_debug("ParameterStudy: Generated {} random numbers creating simulation #{}.", create_counter, run_idx);

            //perform run
            sim.advance(m_tmax);

            log_debug("ParameterStudy: Generated {} random numbers running simulation #{}.",
                      run_rng_counter.get() - create_counter, run_idx);

            //handle result
            if constexpr (should_gather_results) {
                ensemble_result.emplace_back(process_simulation_result(std::move(sim), run_idx));
            }
            else {
                process_simulation_result(std::move(sim), run_idx);
            }
        }

        //Set the counter of our RNG so that future calls of run() produce different parameters.
        m_rng.set_counter(m_rng.get_counter() + rng_totalsequence_counter<uint64_t>(m_num_runs, Counter<uint32_t>(0)));

        if constexpr (should_gather_results) {
            return ensemble_result;
        }
        // else: there is nothing to return, because process_simulation_result returns void
    }

    /**
     * @brief Distribute a number of runs over a number of processes.
     * Processes with low ranks get additional runs, if the number is not evenly divisible.
     * @param num_runs The total number of runs.
     * @param num_procs The total number of processes, i.e. the size of MPI_Comm.
     * @return A vector of size num_procs with the number of runs each process should make.
     */
    static std::vector<size_t> distribute_runs(size_t num_runs, int num_procs)
    {
        assert(num_procs > 0);
        //evenly distribute runs
        //lower processes do one more run if runs are not evenly distributable
        size_t num_runs_local = num_runs / num_procs; //integer division!
        size_t remainder      = num_runs % num_procs;

        std::vector<size_t> run_distribution(num_procs);
        std::fill(run_distribution.begin(), run_distribution.begin() + remainder, num_runs_local + 1);
        std::fill(run_distribution.begin() + remainder, run_distribution.end(), num_runs_local);

        return run_distribution;
    }

    ParameterType m_parameters; ///< Stores parameters used to create a simulation for each run.
    size_t m_num_runs; ///< Total number of runs (i.e. simulations) to do when calling "run".
    Time m_t0, m_tmax; ///< Start and end time for the simulations.
    Step m_dt; ///< Initial step size of the simulation. Some integrators may adapt their step size during simulation.
    RandomNumberGenerator m_rng; ///< The random number generator used by the study.
};

/**
 * @brief Create a GraphSimulation from a parameter graph.
 * @param[in] sampled_graph A graph of models as nodes and mobility parameters as edges, with pre-sampled values.
 * @param[in] t0 Start time of the graph simulation.
 * @param[in] dt_node_sim (Initial) time step used by each node in the GraphSimulation.
 * @param[in] dt_graph_sim Time step used by the GraphSimulation itself.
 */
template <typename FP, class Sim>
auto make_sampled_graph_simulation(const Graph<typename Sim::Model, MobilityParameters<FP>>& sampled_graph, FP t0,
                                   FP dt_node_sim, FP dt_graph_sim)
{
    using SimGraph = Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>;

    SimGraph sim_graph;

    for (auto&& node : sampled_graph.nodes()) {
        sim_graph.add_node(node.id, node.property, t0, dt_node_sim);
    }
    for (auto&& edge : sampled_graph.edges()) {
        sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
    }

    return make_mobility_sim<FP, Sim>(t0, dt_graph_sim, std::move(sim_graph));
}

} // namespace mio

#endif // MIO_COMPARTMENTS_PARAMETER_STUDIES_H
