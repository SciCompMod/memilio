#ifndef PARAMETER_STUDIES_H
#define PARAMETER_STUDIES_H

#include "epidemiology/secir/secir.h"
#include "epidemiology/secir/parameter_space.h"
#include "epidemiology/utils/time_series.h"
#include <epidemiology/migration/migration.h>
#include <epidemiology/model/simulation.h>

#include <cmath>

namespace epi
{

/**
 * Class that performs multiple simulation runs with randomly sampled parameters.
 */
template <class Model>
class ParameterStudy
{

public:

    using HandleSimulationResultFunction = std::function<void(epi::Graph<epi::ModelNode<epi::Simulation<Model>>, epi::MigrationEdge>)>;


    /**
     * create study for graph of compartment models.
     * @param graph graph of parameters
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param graph_sim_dt time step of graph simulation
     * @param num_runs number of runs
     */
    ParameterStudy(const epi::Graph<Model, epi::MigrationParameters>& graph, double t0, double tmax, double graph_sim_dt,
                   size_t num_runs)
        : m_graph(graph)
        , m_num_runs(num_runs)
        , m_t0{t0}
        , m_tmax{tmax}
        , m_dt_graph_sim(graph_sim_dt)
    {
    }

    ParameterStudy(const epi::Graph<Model, epi::MigrationParameters>& graph, double t0, double tmax,
                   double dev_rel, double graph_sim_dt, size_t num_runs)
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
     * @param params SecirParams object 
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param num_runs number of runs in ensemble run
     */
    ParameterStudy(Model const& model, double t0, double tmax, size_t num_runs)
        : m_num_runs{num_runs}
        , m_t0{t0}
        , m_tmax{tmax}
        , m_dt_graph_sim(tmax - t0)
    {
        m_graph.add_node(0, model);
    }

    /**
     * @brief create study for single compartment model with normal distributions.
     * Sets all parameters to normal distribution with specified relative deviation.
     * @param params SecirParams object 
     * @param t0 start time of simulations
     * @param tmax end time of simulations
     * @param dev_rel relative deviation of parameters distributions
     * @param num_runs number of runs in ensemble run
     */
    ParameterStudy(Model const& model, double t0, double tmax, double dev_rel, size_t num_runs)
        : ParameterStudy(model, t0, tmax, num_runs)
    {
        set_params_distributions_normal(m_graph.nodes()[0].property, t0, tmax, dev_rel);
    }

    /*
     * @brief Carry out all simulations in the parameter study.
     * Save memory and enable more runs by immediately processing and/or discarding the result.
     * @param result_processing_function Processing function for simulation results, e.g., output function.
     *                                   Receives the result after each run is completed.
     */
    void run(HandleSimulationResultFunction result_processing_function)
    {
        // Iterate over all parameters in the parameter space
        for (size_t i = 0; i < m_num_runs; i++) {
            auto sim = create_sampled_simulation();
            sim.advance(m_tmax);

            auto result = sim.get_graph();

            result_processing_function(std::move(result));
        }
    }

    /*
     * @brief Carry out all simulations in the parameter study.
     * Convinience function for a few number of runs, but uses a lot of memory.
     * @return vector of results of each run.
     */
    std::vector<epi::Graph<epi::ModelNode<epi::Simulation<Model>>, epi::MigrationEdge>> run()
    {
        std::vector<epi::Graph<epi::ModelNode<epi::Simulation<Model>>, epi::MigrationEdge>> ensemble_result;
        ensemble_result.reserve(m_num_runs);

        run([&ensemble_result](auto r) {
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

    const Model& get_model() const
    {
        return m_graph.nodes()[0].property;
    }

    Model& get_model()
    {
        return m_graph.nodes()[0].property;
    }

    const Graph<Model, MigrationParameters>& get_secir_model_graph() const
    {
        return m_graph;
    }

    Graph<Model, MigrationParameters>& get_secir_model_graph()
    {
        return m_graph;
    }

private:
    //sample parameters and create simulation
    epi::GraphSimulation<epi::Graph<epi::ModelNode<epi::Simulation<Model>>, epi::MigrationEdge>> create_sampled_simulation()
    {
        epi::Graph<epi::ModelNode<epi::Simulation<Model>>, epi::MigrationEdge> sim_graph;

        //sample global parameters
        auto& shared_params_model = m_graph.nodes()[0].property;
        draw_sample_infection(shared_params_model);
        auto& shared_contacts = shared_params_model.parameters.template get<epi::ContactPatterns>();
        shared_contacts.draw_sample();

        for (auto& params_node : m_graph.nodes()) {
            auto& node_model = params_node.property;

            //sample local parameters
            draw_sample_demographics(params_node.property);

            //copy global parameters
            node_model.parameters = shared_params_model.parameters;

            node_model.apply_constraints();
            
            sim_graph.add_node(params_node.id, node_model, m_t0, m_dt_integration);
        }

        for (auto& edge : m_graph.edges()) {
            auto edge_params = edge.property;
            apply_dampings(edge_params.get_coefficients(), shared_contacts.get_dampings(), [&edge_params](auto& v) {
                return make_migration_damping_vector(edge_params.get_coefficients().get_shape(), v);
            });
            sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge_params);
        }

        return make_migration_sim(m_t0, m_dt_graph_sim, sim_graph);
    }

private:
    // Stores Graph with the names and ranges of all parameters
    epi::Graph<Model, epi::MigrationParameters> m_graph;

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

} // namespace epi

#endif // PARAMETER_STUDIES_H
