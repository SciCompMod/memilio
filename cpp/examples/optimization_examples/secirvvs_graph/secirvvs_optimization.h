#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <optional>
#include <cstddef>

#include "models/ode_secirvvs/model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/graph.h"

#include "integrator.h"
#include "ad_type.h"

#include "controls.h"
#include "constraint_parsing.h"

class SecirvvsOptimization
{
public:
    SecirvvsOptimization(
        const mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::osecirvvs::Model<double>>>,
                         mio::MobilityEdge<double>>& graph,
        double t0, double tmax);

    void addPathConstraint(const std::string& name, const std::pair<double, double>& bounds);
    void addPathConstraint(const std::string& name, const std::pair<double, double>& bounds, size_t node_index);

    void addTerminalConstraint(const std::string& name, const std::pair<double, double>& bounds);
    void addTerminalConstraint(const std::string& name, const std::pair<double, double>& bounds, size_t node_index);

    const mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::osecirvvs::Model<double>>>,
                     mio::MobilityEdge<double>>&
    graph() const;
    double t0() const;
    double tmax() const;

    size_t numControlIntervals() const;
    void setNumControlIntervals(size_t num_control_intervals);
    size_t pcResolution() const;
    void setPcResolution(size_t pc_resolution);

    bool randomStart() const;
    void setRandomStart(bool random_start);

    IntegratorType integratorType() const;
    void setIntegratorType(IntegratorType integrator_type);

    size_t integratorResolution() const;
    void setIntegratorResolution(size_t integrator_resolution);

    ADType adEvalF() const;
    void setAdEvalF(ADType ad_eval_f);

    ADType adEvalJac() const;
    void setAdEvalJac(ADType ad_eval_jac);

    // --- IPOPT Options ---
    int ipopt_verbose() const;
    void setIpoptVerbose(int ipopt_verbose);

    int ipopt_printFrequency() const;
    void setIpoptPrintFrequency(int ipopt_print_frequency);

    bool ipopt_printTimings() const;
    void setIpoptPrintTimings(bool ipopt_print_timings);

    int ipopt_maxIter() const;
    void setIpoptMaxIter(int ipopt_max_iter);

    double ipopt_tolerance() const;
    void setIpoptTolerance(double ipopt_tolerance);

    void configure();

    size_t numGraphNodes() const;
    size_t numGraphEdges() const;

    size_t numIntervals() const;
    double dt() const;

    size_t numControls() const;
    size_t numPathConstraints() const;
    size_t numTerminalConstraints() const;

    std::vector<ConstraintPolicy<double>>& pathConstraints();
    std::vector<ConstraintPolicy<double>>& terminalConstraints();

    void validateConstraints();

    // // --- Solve ---
    // void solve()
    // {

    //     validateConstraints();
    //     setupOptimizationProblem();
    //     // IPOPT solving logic would go here
    //     std::cout << "Starting optimization with " << m_num_control_intervals << " control intervals...\n";
    //     // Actual IPOPT integration would be implemented here
    //     std::cout << "Optimization completed.\n";

    //     // TODO: Implement transcription to NLP, setup IPOPT, and solve
    //     // 1. Discretize controls
    //     // 2. Apply path and terminal constraints
    //     // 3. Configure integrator and IPOPT options
    //     // 4. Call IPOPT
    //     // 5. Extract solution
    // }

private:
    // Problem definition
    const mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::osecirvvs::Model<double>>>,
                     mio::MobilityEdge<double>>& m_graph;
    double m_t0;
    double m_tmax;

    // Vector of control policies for various parameters: { name, {lower, upper}, effectiveness, cost_penalty }
    std::vector<ControlPolicy<double>> m_control_policy;

    // Vector of constraint bounds: { name, {lower, upper}, node_index (optional) }
    std::vector<ConstraintPolicy<double>> m_path_constraints;
    std::vector<ConstraintPolicy<double>> m_terminal_constraints;

    // Discretization parameters
    size_t m_num_control_intervals;
    size_t m_pc_resolution;
    bool m_random_start;

    IntegratorType m_integrator_type;
    size_t m_integrator_resolution;

    ADType m_ad_eval_f;
    ADType m_ad_eval_jac;

    // IPOPT options
    int m_ipopt_verbose;
    int m_ipopt_print_frequency;
    bool m_ipopt_print_timings;
    int m_ipopt_max_iter;
    double m_ipopt_tolerance;

    size_t m_num_graph_nodes;
    size_t m_num_graph_edges;

    size_t m_num_intervals;
    double m_dt;

    size_t m_num_controls;
    size_t m_num_path_constraints;
    size_t m_num_terminal_constraints;
    size_t m_num_graph_simulation_controls;
    size_t m_num_graph_simulation_path_constraints;
    size_t m_num_graph_simulation_terminal_constraints;

    void set_restrictive_dampings(
        mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::osecirvvs::Model<double>>>,
                   mio::MobilityEdge<double>>& graph_simulation);

    void update_path_constraint(auto graph_sim_mobility, std::vector<double>& terminal_constraint_values);
    void update_terminal_constraint(auto graph_sim_mobility, std::vector<double>& terminal_constraint_values);
};
