#ifndef PARAMETER_STUDIES_H
#define PARAMETER_STUDIES_H

#include "epidemiology/secir/secir.h"
#include "epidemiology/secir/parameter_space.h"
#include "epidemiology/utils/time_series.h"
#include <epidemiology/migration/migration.h>

namespace epi
{

using HandleSimulationResultFunction = std::function<void(const SecirParams&, const TimeSeries<double>&, int node)>;

auto DummyHandleResultFunction = [](const SecirParams& params, const TimeSeries<double>&, int node) {};

/**
 * Class that performs multiple simulation runs with randomly sampled parameters.
 */
class ParameterStudy
{
public:
    /**
     * create study for graph of compartment models.
     * @param graph graph of parameters
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param graph_sim_dt time step of graph simulation
     * @param num_runs number of runs
     */
    ParameterStudy(const epi::Graph<epi::SecirParams, epi::MigrationEdge>& graph,
                   double t0, double tmax, double graph_sim_dt, size_t num_runs);

    /**
     * @brief Create study for single compartment model.
     * @param params SecirParams object 
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param num_runs number of runs in ensemble run
     */
    ParameterStudy(SecirParams const& params, double t0, double tmax, size_t num_runs);

    /**
     * @brief create study for single compartment model with normal distributions.
     * Sets all parameters to normal distribution with specified relative deviation.
     * @param params SecirParams object 
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param dev_rel relative deviation of parameters distributions
     * @param num_runs number of runs in ensemble run
     */
    ParameterStudy(SecirParams const& params, double t0, double tmax, double dev_rel, size_t num_runs);

    /*
     * @brief Carry out all simulations in the parameter study.
     * @param[in] result_processing_function Processing function for simulation results, e.g., output function
     */
    std::vector<epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge>>
    run(HandleSimulationResultFunction result_processing_function = DummyHandleResultFunction);

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

    const SecirParams& get_secir_params() const
    {
        return m_graph.nodes()[0];
    }

    SecirParams& get_secir_params()
    {
        return m_graph.nodes()[0];
    }

private:
    // Stores Graph with the names and ranges of all parameters
    epi::Graph<epi::SecirParams, epi::MigrationEdge> m_graph;

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

inline ParameterStudy::ParameterStudy(const epi::Graph<epi::SecirParams, epi::MigrationEdge>& graph, double t0, double tmax, double graph_sim_dt,
                                      size_t num_runs)
    : m_graph(graph)
    , m_num_runs(num_runs)
    , m_t0{t0}
    , m_tmax{tmax}
    , m_dt_graph_sim(graph_sim_dt)
{
}

inline ParameterStudy::ParameterStudy(SecirParams const& params, double t0, double tmax, size_t num_runs)
    : m_num_runs{num_runs}
    , m_t0{t0}
    , m_tmax{tmax}
    , m_dt_graph_sim(tmax - t0)
{
    m_graph.add_node(params);
}

inline ParameterStudy::ParameterStudy(SecirParams const& params, double t0, double tmax, double dev_rel,
                                      size_t num_runs)
    : ParameterStudy(params, t0, tmax, num_runs)
{
    set_params_distributions_normal(m_graph.nodes()[0], t0, tmax, dev_rel);
}

inline std::vector<epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge>>
ParameterStudy::run(HandleSimulationResultFunction simulation_result_function)
{
    std::vector<epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge>> ensemble_result;

    // Iterate over all parameters in the parameter space
    for (size_t i = 0; i < m_num_runs; i++) {
        epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge> sim_graph;

        for (auto& params_node : m_graph.nodes()) {
            draw_sample(params_node);
            sim_graph.add_node(params_node, m_t0, m_dt_integration);
        }

        for (auto& edge : m_graph.edges()) {
            sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property.coefficients);
        }

        // Call the simulation function
        auto sim = make_migration_sim(m_t0, m_dt_graph_sim, sim_graph);
        sim.advance(m_tmax);

        auto result = sim.get_graph();

        int node_id = 0;
        for (auto& node : result.nodes()) {
            simulation_result_function(node.model.get_params(), node.model.get_result(), node_id);
            node_id++;
        }

        ensemble_result.push_back(result);
    }

    return ensemble_result;
}

} // namespace epi

#endif // PARAMETER_STUDIES_H
