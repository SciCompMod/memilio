
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"

#include <omp.h>

bool age_groups = true;

void set_contact_matrices(mio::oseir::Parameters<double>& params)
{
    if (age_groups) {
        Eigen::MatrixXd contact_matrix_eigen(6, 6);
        contact_matrix_eigen << 3.9547, 1.1002, 2.9472, 2.05, 0.3733, 0.0445, 0.3327, 3.5892, 1.236, 1.9208, 0.2681,
            0.0161, 0.246, 0.7124, 5.6518, 3.2939, 0.2043, 0.0109, 0.1742, 0.8897, 3.3124, 4.5406, 0.4262, 0.0214,
            0.0458, 0.1939, 0.5782, 1.3825, 1.473, 0.0704, 0.1083, 0.1448, 0.4728, 0.9767, 0.6266, 0.1724;
        mio::ContactMatrixGroup& contact_matrix =
            params.template get<mio::oseir::ContactPatterns<>>().get_cont_freq_mat();
        contact_matrix[0].get_baseline() = contact_matrix_eigen;
    }
    else {
        mio::ContactMatrixGroup& contact_matrix = params.get<mio::oseir::ContactPatterns<>>().get_cont_freq_mat();
        contact_matrix[0].get_baseline().setConstant(7.95);
    }
}

/**
 * Set epidemiological parameters of Sars-CoV-2 for a immune-naive
 * population and wild type variant.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
void set_covid_parameters(mio::oseir::Parameters<double>& params)
{
    params.template set<mio::oseir::TimeExposed<>>(3.335);

    if (age_groups) {
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(0)] = 8.0096875;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(1)] = 8.0096875;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(2)] = 8.2182;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(3)] = 8.1158;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(4)] = 8.033;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(5)] = 7.985;

        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(0)] = 0.03;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(1)] = 0.06;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(2)] = 0.06;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(3)] = 0.06;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(4)] = 0.09;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(5)] = 0.175;
    }
    else {
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(0)] = 0.07333;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(0)]                     = 8.097612257;
    }
}

void set_population_data(mio::oseir::Parameters<double>& params,
                         mio::Graph<mio::SimulationNode<mio::Simulation<ScalarType, mio::oseir::Model<>>>,
                                    mio::MobilityEdge<>>& params_graph,
                         size_t number_regions)
{
    std::vector<mio::oseir::Model<double>> nodes(number_regions,
                                                 mio::oseir::Model(int(size_t(params.get_num_groups()))));

    mio::Populations<double, mio::AgeGroup, mio::oseir::InfectionState> population(
        {params.get_num_groups(), mio::oseir::InfectionState::Count});

    for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
        population[{i, mio::oseir::InfectionState::Susceptible}] = 10000;
    }
    for (auto& node : nodes) {
        node.parameters  = params;
        node.populations = population;
    }
    // for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
    nodes[0].populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] += 100;
    nodes[0].populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] -= 100;
    // }

    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {
        params_graph.add_node(node_idx, nodes[node_idx]);
    }
}

void set_parameters_and_population(mio::Graph<mio::SimulationNode<mio::Simulation<ScalarType, mio::oseir::Model<>>>,
                                              mio::MobilityEdge<>>& params_graph,
                                   size_t number_regions)
{
    int num_age_groups = 1;
    if (age_groups) {
        num_age_groups = 6;
    }

    mio::oseir::Parameters<double> params(num_age_groups);

    set_covid_parameters(params);

    // set contact matrix
    set_contact_matrices(params);

    set_population_data(params, params_graph, number_regions);

    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            double commuter_coeff_ij = 1. / (2 * number_regions);
            if (county_idx_i == county_idx_j) {
                commuter_coeff_ij = 0;
            }
            params_graph.add_edge(
                county_idx_i, county_idx_j,
                Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count * size_t(params.get_num_groups()),
                                          commuter_coeff_ij));
        }
    }
}

void simulate(ScalarType tol, ScalarType tmax)
{
    ScalarType t0         = 0.;
    ScalarType dt         = 0.5;
    size_t number_regions = 100;

    mio::Graph<mio::SimulationNode<mio::Simulation<ScalarType, mio::oseir::Model<>>>, mio::MobilityEdge<>> params_graph;

    set_parameters_and_population(params_graph, number_regions);

    using DefaultIntegratorCore =
        mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>;

    for (auto& node : params_graph.nodes()) {
        node.property.get_simulation().set_integrator(std::make_shared<DefaultIntegratorCore>(tol));
    }

    auto sim = mio::make_mobility_sim(t0, dt, std::move(params_graph));
    sim.advance(tmax);

    auto result_graph = std::move(sim).get_graph();

    std::cout << " \"Regions\": " << number_regions << "," << std::endl;
    std::cout << " \"Steps Hotspot\": " << result_graph.nodes()[0].property.get_result().get_num_time_points() - 1
              << "," << std::endl;
    std::cout << " \"Steps other Regions\": "
              << result_graph.nodes()[99].property.get_result().get_num_time_points() - 1;
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::off);
    const ScalarType tmax = 10;
    ScalarType tol        = 1e-2;

    if (argc > 1) {
        tol = std::stod(argv[1]);
    }

    std::cout << "{ \"Absolute tolerance\": " << tol << ", " << std::endl;
    simulate(tol, tmax);
    std::cout << "}," << std::endl;

    return 0;
}
