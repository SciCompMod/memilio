
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"
#include "memilio/io/result_io.h"
#include "memilio/io/epi_data.h"

/**
 * indices of contact matrix corresponding to locations where contacts occur.
 */
enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

/**
 * Set contact matrices.
 * Reads contact matrices from files in the data directory.
 * @param data_dir data directory.
 * @param params Object that the contact matrices will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::oseir::Parameters<double>& params)
{
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(
            auto&& baseline,
            mio::read_mobility_plain(
                (data_dir / "Germany" / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));
        BOOST_OUTCOME_TRY(
            auto&& minimum,
            mio::read_mobility_plain(
                (data_dir / "Germany" / "contacts" / ("minimum_" + contact_location.second + ".txt")).string()));
        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = minimum;
    }
    params.get<mio::oseir::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    printf("Setting contact matrices successful.\n");
    return mio::success();
}

/**
 * Set epidemiological parameters of Sars-CoV-2 for a immune-naive
 * population and wild type variant.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::oseir::Parameters<double>& params)
{
    params.template set<mio::oseir::TimeExposed<>>(3.335);
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

    printf("Setting epidemiological parameters successful.\n");
    return mio::success();
}

mio::IOResult<std::vector<mio::oseir::Model<double>>>
set_population_data(const fs::path& data_dir, mio::oseir::Parameters<double>& params, std::vector<int> node_ids)
{
    size_t number_regions = node_ids.size();

    std::vector<mio::oseir::Model<double>> nodes(number_regions,
                                                 mio::oseir::Model(int(size_t(params.get_num_groups()))));

    for (auto& node : nodes) {
        node.parameters = params;
    }

    BOOST_OUTCOME_TRY(const auto&& population_data,
                      mio::read_population_data(
                          (data_dir / "Germany" / "pydata" / "county_current_population_nrw.json").string(), true));

    std::vector<std::vector<double>> vnum_population(node_ids.size(),
                                                     std::vector<double>((size_t)params.get_num_groups(), 0.0));

    for (auto&& entry : population_data) {
        auto it = std::find_if(node_ids.begin(), node_ids.end(), [&entry](auto r) {
            return r == 0 ||
                   (entry.county_id && mio::regions::StateId(r) == mio::regions::get_state_id(int(*entry.county_id))) ||
                   (entry.county_id && mio::regions::CountyId(r) == *entry.county_id) ||
                   (entry.district_id && mio::regions::DistrictId(r) == *entry.district_id);
        });
        if (it != node_ids.end()) {
            auto region_idx      = size_t(it - node_ids.begin());
            auto& num_population = vnum_population[region_idx];
            for (size_t age = 0; age < num_population.size(); age++) {
                num_population[age] += entry.population[mio::AgeGroup(age)];
            }
        }
    }

    for (size_t region = 0; region < node_ids.size(); region++) {
        auto num_groups = nodes[region].parameters.get_num_groups();
        for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
            nodes[region].populations.template set_difference_from_group_total<mio::AgeGroup>(
                {i, mio::oseir::InfectionState::Susceptible}, vnum_population[region][size_t(i)]);
        }
    }
    nodes[27].populations[{mio::AgeGroup(4), mio::oseir::InfectionState::Susceptible}] -= 100;
    nodes[27].populations[{mio::AgeGroup(4), mio::oseir::InfectionState::Exposed}] += 100;

    return mio::success(nodes);
}

mio::IOResult<void> run(const fs::path& data_dir, double t0, double tmax, double dt)
{
    mio::set_log_level(mio::LogLevel::off);
    // global parameters
    const int num_age_groups = 6;

    mio::oseir::Parameters<double> params(num_age_groups);

    BOOST_OUTCOME_TRY(set_covid_parameters(params));

    // set contact matrix
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));

    // graph of counties with populations and local parameters
    // and mobility between counties
    mio::Graph<mio::SimulationNode<mio::Simulation<ScalarType, mio::oseir::Model<>>>, mio::MobilityEdge<>> params_graph;

    BOOST_OUTCOME_TRY(
        auto&& node_ids,
        mio::get_node_ids((data_dir / "Germany" / "pydata" / "county_current_population_nrw.json").string(), true,
                          true));

    BOOST_OUTCOME_TRY(auto&& nodes, set_population_data(data_dir, params, node_ids));
    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {
        params_graph.add_node(node_ids[node_idx], nodes[node_idx]);
    }
    printf("Setting population from data successful.\n");

    BOOST_OUTCOME_TRY(
        auto&& mobility_data_commuter,
        mio::read_mobility_plain((data_dir / "Germany" / "mobility" / "commuter_mobility_2022_nrw.txt").string()));
    if (mobility_data_commuter.rows() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_commuter.cols() != Eigen::Index(params_graph.nodes().size())) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Mobility matrices do not have the correct size. You may need to run "
                            "transformMobilitydata.py from pycode memilio epidata package.");
    }

    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            auto& populations = params_graph.nodes()[county_idx_i].property.get_simulation().get_model().populations;

            auto commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) / populations.get_total();
            params_graph.add_edge(
                county_idx_i, county_idx_j,
                Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count * size_t(params.get_num_groups()),
                                          commuter_coeff_ij));
        }
    }

    auto sim = mio::make_mobility_sim(t0, dt, std::move(params_graph));

    printf("Start Simulation\n");
    sim.advance(tmax);

    auto result_graph = std::move(sim).get_graph();
    auto result       = mio::interpolate_simulation_result(result_graph);

    std::vector<int> county_ids(result_graph.nodes().size());
    std::transform(result_graph.nodes().begin(), result_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    auto save_result_status = save_result(result, county_ids, num_age_groups, "graph_result_nrw.h5");

    return mio::success();
}

int main()
{
    const auto t0   = 0.;
    const auto tmax = 100.;
    const auto dt   = 0.5; //time step of mobility, daily mobility every second step

    const std::string& data_dir = "";

    auto result = run(data_dir, t0, tmax, dt);

    return 0;
}
