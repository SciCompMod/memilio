#ifndef PARAMETER_STUDIES_H
#define PARAMETER_STUDIES_H

#include "epidemiology/secir/secir.h"
#include "epidemiology/secir/parameter_space.h"
#include "epidemiology/utils/time_series.h"
#include <epidemiology/migration/migration.h>

namespace epi
{
using HandleSimulationResultFunction = std::function<void(const SecirParams&, const TimeSeries<double>&, int node)>;

// The function type for the kind of simulation that we want to run
using secir_simulation_function_t =
    std::function<GraphSimulation<epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge>>(
        double t0, double dt, epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge> sim_graph)>;

auto dummy_func = [](auto& params, const TimeSeries<double>&, int node) {};

// TODO: document class

class ParameterStudy
{
public:
    ParameterStudy(secir_simulation_function_t const& simu_func,
                   epi::Graph<epi::ModelNode<epi::SecirParams>, epi::MigrationEdge> graph, size_t num_runs, double t0,
                   double tmax);

    /* 
     * @brief Constructor from file name
     * @param[in] simu_func simulation function for Secir simulation
     * @param[in] params SecirParams object 
     * @param[in] t0 start time of simulations
     * @param[in] tmax end time of simulations
     * @param[in] num_runs number of runs in ensemble run
     */
    ParameterStudy(secir_simulation_function_t const& simu_func, SecirParams&& params, double t0, double tmax,
                   size_t num_runs);

    /* 
     * @brief Constructor from contact frequency matrix and parameter vector
     * @param[in] simu_func simulation function for Secir simulation
     * @param[in] params SecirParams object 
     * @param[in] t0 start time of simulations
     * @param[in] tmax end time of simulations
     * @param[in] num_runs number of runs in ensemble run
     */
    ParameterStudy(secir_simulation_function_t const& simu_func, SecirParams const& params, double t0, double tmax,
                   size_t num_runs);

    ParameterStudy(secir_simulation_function_t const& simu_func, SecirParams const& params, double t0, double tmax,
                   double dev_rel, size_t num_runs);

    /*
     * @brief Carry out all simulations in the parameter study.
     * @param[in] result_processing_function Processing function for simulation results, e.g., output function
     */
    std::vector<epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge>>
    run(HandleSimulationResultFunction result_processing_function = dummy_func);

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
        return m_graph.nodes()[0].model;
    }

    SecirParams& get_secir_params()
    {
        return m_graph.nodes()[0].model;
    }

private:
    // The path of the file storing the parameter ranges
    std::string parameter_file;

    // Stores Graph with the names and ranges of all parameters
    epi::Graph<epi::ModelNode<epi::SecirParams>, epi::MigrationEdge> m_graph;

    // The function that carries out our simulation
    secir_simulation_function_t simulation_function;

    size_t m_num_runs = 100;

    // Start time (should be the same for all simulations)
    double m_t0 = 0.0;
    // End time (should be the same for all simulations)
    double m_tmax = 400.0;
    // adaptive time step (will be corrected if too large/small)
    double m_dt = 0.1;
};

inline ParameterStudy::ParameterStudy(secir_simulation_function_t const& simu_func,
                                      epi::Graph<epi::ModelNode<epi::SecirParams>, epi::MigrationEdge> graph,
                                      size_t num_runs, double t0, double tmax)
    : simulation_function(simu_func)
    , m_graph(std::move(graph))
    , m_num_runs(num_runs)
    , m_t0{t0}
    , m_tmax{tmax}
{
}

inline ParameterStudy::ParameterStudy(const secir_simulation_function_t& simu_func, SecirParams&& params, double t0,
                                      double tmax, size_t num_runs)
    : simulation_function(simu_func)
    , m_num_runs(num_runs)
    , m_t0{t0}
    , m_tmax{tmax}
{
    m_graph.add_node(params);
}

inline ParameterStudy::ParameterStudy(secir_simulation_function_t const& simu_func, SecirParams const& params,
                                      double t0, double tmax, size_t num_runs)
    : simulation_function{simu_func}
    , m_num_runs{num_runs}
    , m_t0{t0}
    , m_tmax{tmax}
{
    m_graph.add_node(params);
}

inline ParameterStudy::ParameterStudy(secir_simulation_function_t const& simu_func, SecirParams const& params,
                                      double t0, double tmax, double dev_rel, size_t num_runs)
    : simulation_function{simu_func}
    , m_num_runs{num_runs}
    , m_t0{t0}
    , m_tmax{tmax}
{
    m_graph.add_node(params);
    set_params_distributions_normal(m_graph.nodes()[0].model, t0, tmax, dev_rel);
}

inline std::vector<epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge>>
ParameterStudy::run(HandleSimulationResultFunction simulation_result_function)
{
    std::vector<epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge>> ensemble_result;

    // Iterate over all parameters in the parameter space
    for (size_t i = 0; i < m_num_runs; i++) {
        epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge> sim_graph;

        for (auto& node : m_graph.nodes()) {
            SecirParams params_sample = node.model;
            draw_sample(params_sample);
            sim_graph.add_node(params_sample, m_t0);
        }

        for (auto& edge : m_graph.edges()) {
            sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property.coefficients);
        }

        // Call the simulation function
        auto sim = simulation_function(m_t0, m_dt, sim_graph);
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
