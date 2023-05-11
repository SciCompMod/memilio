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

#include "memilio/utils/time_series.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"

#include <cmath>

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
    using Simulation = S;

    /**
     * create study for graph of compartment models.
     * @param graph graph of parameters
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param graph_sim_dt time step of graph simulation
     * @param num_runs number of runs
     */
    ParameterStudy(const mio::Graph<typename Simulation::Model, mio::MigrationParameters>& graph, double t0,
                   double tmax, double graph_sim_dt, size_t num_runs)
        : m_graph(graph)
        , m_num_runs(num_runs)
        , m_t0{t0}
        , m_tmax{tmax}
        , m_dt_graph_sim(graph_sim_dt)
    {
    }

    ParameterStudy(const mio::Graph<typename Simulation::Model, mio::MigrationParameters>& graph, double t0,
                   double tmax, double dev_rel, double graph_sim_dt, size_t num_runs)
        : m_graph(graph)
        , m_num_runs(num_runs)
        , m_t0{t0}
        , m_tmax{tmax}
        , m_dt_graph_sim(graph_sim_dt)
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
        : m_num_runs{num_runs}
        , m_t0{t0}
        , m_tmax{tmax}
        , m_dt_graph_sim(tmax - t0)
    {
        m_graph.add_node(0, model);
    }

    /*
     * @brief Carry out all simulations in the parameter study.
     * Save memory and enable more runs by immediately processing and/or discarding the result.
     * @param result_processing_function Processing function for simulation results, e.g., output function.
     *                                   Receives the result after each run is completed.
     */
    template <class SampleGraphFunction, class HandleSimulationResultFunction>
    void run(SampleGraphFunction sample_graph, HandleSimulationResultFunction result_processing_function)
    {
        // Iterate over all parameters in the parameter space
        for (size_t i = 0; i < m_num_runs; i++) {
            auto sim = create_sampled_simulation(sample_graph);
            sim.advance(m_tmax);

            result_processing_function(std::move(sim).get_graph());
        }
    }

    /*
     * @brief Carry out all simulations in the parameter study.
     * Convenience function for a few number of runs, but uses a lot of memory.
     * @return vector of results of each run.
     */
    template <class SampleGraphFunction>
    std::vector<mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>> run(SampleGraphFunction sample_graph)
    {
        std::vector<mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>> ensemble_result;
        ensemble_result.reserve(m_num_runs);

        run(sample_graph, [&ensemble_result](auto&& r) {
            ensemble_result.emplace_back(std::move(r));
        });

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
    const Graph<typename Simulation::Model, MigrationParameters>& get_model_graph() const
    {
        return m_graph;
    }
    Graph<typename Simulation::Model, MigrationParameters>& get_model_graph()
    {
        return m_graph;
    }
    /** @} */

private:
    //sample parameters and create simulation
    template <class SampleGraphFunction>
    mio::GraphSimulation<mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>>
    create_sampled_simulation(SampleGraphFunction sample_graph)
    {
        mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge> sim_graph;

        auto sampled_graph = sample_graph(m_graph);
        for (auto&& node : sampled_graph.nodes()) {
            sim_graph.add_node(node.id, node.property, m_t0, m_dt_integration);
        }
        for (auto&& edge : sampled_graph.edges()) {
            sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
        }

        return make_migration_sim(m_t0, m_dt_graph_sim, std::move(sim_graph));
    }

private:
    // Stores Graph with the names and ranges of all parameters
    mio::Graph<typename Simulation::Model, mio::MigrationParameters> m_graph;

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
