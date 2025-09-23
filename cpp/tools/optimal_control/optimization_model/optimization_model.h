#pragma once

#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <tuple>
#include <cstddef>
#include <algorithm>

#include "boost/outcome/try.hpp"
#include "boost/outcome/result.hpp"
#include "boost/optional.hpp"
#include "boost/filesystem.hpp"

#include "memilio/utils/date.h"
#include "memilio/io/mobility_io.h"
#include "models/ode_secirvvs/model.h"
#include "models/ode_secirvvs/parameters_io.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"

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
                                                                         {ContactLocation::Other, "other"}
                                                                        };

template <typename FP, class Tag>
void set_parameters(mio::osecirvvs::Parameters<FP>& params, const std::vector<FP>& parameter_values)
{
    std::copy(parameter_values.begin(), parameter_values.end(), params.template get<Tag>());
}

template <typename FP, class Tag>
void set_parameters(mio::osecirvvs::Parameters<FP>& params, double parameter_values)
{
    std::fill(params.template get<Tag>().begin(), params.template get<Tag>().end(), FP(parameter_values));
}
template <class ContactLocation>
class OptimizationModel
{
public:
    template <typename FP>    
    using ContactPatterns = mio::osecirvvs::ContactPatterns<FP>
    template <typename FP>
    using EMBModel = mio::osecirvvs::Model<FP>
    template <typename FP>
    using Graph = mio::Graph<mio::SimulationNode<mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>>


    OptimizationModel(boost::filesystem::path data_directory, double t0, double tmax, bool use_county)
        : m_data_directory(std::move(data_directory))
        , m_t0(t0)
        , m_tmax(tmax)
        , m_use_county(use_county)
    {
    }

    const boost::filesystem::path& data_directory() const
    {
        return m_data_directory;
    }

    double t0() const
    {
        return m_t0;
    }
    double tmax() const
    {
        return m_tmax;
    }
    int num_age_groups() const
    {
        return m_num_age_groups;
    }

    template <typename FP>
    Graph<FP> create_model()
    {
        Graph<FP> graph;
        auto out = set_initial_values<FP>(graph);

        return graph;
    }

    void add_contact_location(Eigen::MatrixXd baseline, Eigen::MatrixXd minimum, std::string location_name)
    {
        m_contact_matrices.emplace_back(std::make_tuple(location_name, baseline, minimum));
    }

    void add_contact_location(Eigen::MatrixXd baseline, std::string location_name)
    {
        auto&& minimum = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic>::Zero(baseline.rows(), baseline.cols());
        m_contact_matrices.emplace_back(std::make_tuple(location_name, baseline, minimum));
    }


private:

    template <typename FP>
    mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters<FP>& params)
    {
        BOOST_OUTCOME_TRY(auto&& parameter_list, mio::read_json((m_data_directory / "parameter_list.json").string()));

        std::map<std::string, std::string> id_to_name{};
        for(auto entry:  parameter_list)
        {
            id_to_name[entry["id"].asString()] = entry["name"].asString();
        }

        // only uses group 0 for all parameters
        BOOST_OUTCOME_TRY(auto&& scenario_data_run, mio::read_json((m_data_directory / "scenario_data_run.json").string()));
        std::map<std::string, double> parameter_values{};
        for(auto parameter:  scenario_data_run["modelParameters"])
        {
            std::string parameter_name = id_to_name[parameter["parameterId"].asString()];
            parameter_values[parameter_name] = 0.5 * (parameter["values"][0]["valueMin"].asDouble() + 
                                                        parameter["values"][0]["valueMax"].asDouble());
        }

        //times
        set_parameters<FP, mio::osecirvvs::TimeExposed<FP>>(params, parameter_values["TimeExposed"]);
        set_parameters<FP, mio::osecirvvs::TimeInfectedNoSymptoms<FP>>(params, parameter_values["TimeInfectedNoSymptoms"]);
        set_parameters<FP, mio::osecirvvs::TimeInfectedSymptoms<FP>>(params, parameter_values["TimeInfectedSymptoms"]);
        set_parameters<FP, mio::osecirvvs::TimeInfectedSevere<FP>>(params, parameter_values["TimeInfectedSevere"]);
        set_parameters<FP, mio::osecirvvs::TimeInfectedCritical<FP>>(params, parameter_values["TimeInfectedCritical"]);

        //probabilities
        set_parameters<FP, mio::osecirvvs::TransmissionProbabilityOnContact<FP>>(params, parameter_values["TransmissionProbabilityOnContact"]);
        set_parameters<FP, mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>>(params, parameter_values["RelativeTransmissionNoSymptoms"]);
        // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
        // depends on incidence and test and trace capacity
        set_parameters<FP, mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>>(params, parameter_values["RiskOfInfectionFromSymptomatic"]);
        set_parameters<FP, mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>>(params, parameter_values["MaxRiskOfInfectionFromSymptomatic"]);
        set_parameters<FP, mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>>(params, parameter_values["RecoveredPerInfectedNoSymptoms"]);
        set_parameters<FP, mio::osecirvvs::SeverePerInfectedSymptoms<FP>>(params, parameter_values["SeverePerInfectedSymptoms"]);
        set_parameters<FP, mio::osecirvvs::CriticalPerSevere<FP>>(params, parameter_values["CriticalPerSevere"]);
        set_parameters<FP, mio::osecirvvs::DeathsPerCritical<FP>>(params, parameter_values["DeathsPerCritical"]);

        set_parameters<FP, mio::osecirvvs::ReducExposedPartialImmunity<FP>>(params, parameter_values["ReducedExposedPartialImmunity"]);
        set_parameters<FP, mio::osecirvvs::ReducExposedImprovedImmunity<FP>>(params, parameter_values["ReducedExposedImprovedImmunity"]);
        set_parameters<FP, mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>>(params, parameter_values["ReducedInfectedSymptomsPartialImmunity"]);
        set_parameters<FP, mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>>(params, parameter_values["ReducedInfectedSymptomsImprovedImmunity"]);
        set_parameters<FP, mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>(params, parameter_values["ReducedInfectedSevereCriticalDeadPartialImmunity"]);
        set_parameters<FP, mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>(params, parameter_values["ReducedInfectedSevereCriticalDeadImprovedImmunity"]);
        set_parameters<FP, mio::osecirvvs::ReducTimeInfectedMild<FP>>(params, parameter_values["ReducedTimeInfectedMild"]);

        BOOST_OUTCOME_TRY(mio::Date start_date, mio::parse_date(scenario_data_run["startDate"].asString()));
        params.template get<mio::osecirvvs::StartDay>() = get_day_in_year(start_date);
        params.template get<mio::osecirvvs::Seasonality<FP>>() = parameter_values["Seasonality"];

        return mio::success();
    }

    template <typename FP>
    mio::IOResult<void> set_contact_matrices(mio::osecirvvs::Parameters<FP>& params)
    {
        auto contact_matrices = mio::ContactMatrixGroup<FP>(contact_locations.size(), size_t(params.get_num_groups()));

        size_t idx = 0;
        for (auto&& contact_location : m_contact_matrices) {
            contact_matrices[idx].get_baseline() = std::get<1>(contact_location).cast<FP>();
            contact_matrices[idx].get_minimum()  = std::get<2>(contact_location).cast<FP>();
            idx++;
        }
        params.template get<ContactPatterns<FP>>() = mio::UncertainContactMatrix<FP>(contact_matrices);

        return mio::success();
    }

    template <typename FP>
    mio::IOResult<void> set_population_data_germany(Graph<FP>& params_graph, mio::osecirvvs::Parameters<FP>& params)
    {
        auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 2.5);
        auto scaling_factor_icu      = 1.0;
        auto tnt_capacity_factor     = 7.5 / 100000.;


        std::vector<mio::osecirvvs::Model<FP>> nodes(1, mio::osecirvvs::Model<FP>(int(size_t(params.get_num_groups()))));
        for (auto& node : nodes) {
            node.parameters = params;
        }

        // std::vector<int> node_ids = {0};
        mio::Date start_date(2022, 12, 1);
        mio::Date end_date(2020, 12, 12);
        std::string pydata_dir = mio::path_join(m_data_directory.string(), "Germany", "pydata");
        BOOST_OUTCOME_TRY(mio::osecirvvs::read_input_data_germany<mio::osecirvvs::Model<FP>>(nodes, start_date, scaling_factor_infected, scaling_factor_icu, pydata_dir));

        params_graph.add_node(0, nodes[0]);

        mio::unused(tnt_capacity_factor);

        return mio::success();
    }

    template <typename FP>
    mio::IOResult<void> set_population_data_county(Graph<FP>& params_graph, mio::osecirvvs::Parameters<FP>& params)
    {
        mio::Graph<EMBModel<FP>, mio::MobilityParameters<FP>> param_graph;

        auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 2.5);
        auto scaling_factor_icu      = 1.0;
        auto tnt_capacity_factor     = 7.5 / 100000.;
        auto mobile_compartments     = {mio::osecirvvs::InfectionState::SusceptibleNaive, mio::osecirvvs::InfectionState::ExposedNaive,
                                        mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive,
                                        mio::osecirvvs::InfectionState::InfectedSymptomsNaive};

        // graph of counties with populations and local parameters
        // and mobility between counties
        const auto& read_function_nodes = mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model<FP>>;
        // const auto& read_function_edges = mio::read_mobility_plain;
        const auto& node_id_function    = mio::get_node_ids;

        const auto& set_node_function =
            mio::set_nodes<mio::osecirvvs::TestAndTraceCapacity<FP>, ContactPatterns<FP>,
                        mio::osecirvvs::Model<FP>, mio::MobilityParameters<FP>, mio::osecirvvs::Parameters<FP>,
                        decltype(read_function_nodes), decltype(node_id_function)>;
        const auto& set_edge_function =
            mio::set_edges<ContactLocation, mio::osecirvvs::Model<FP>, mio::MobilityParameters<FP>,
                        mio::MobilityCoefficientGroup, mio::osecirvvs::InfectionState, decltype(read_function_edges)>;

        mio::Date start_date(2022, 12, 1);
        mio::Date end_date(2020, 12, 12);
        std::string pydata_dir = mio::path_join(m_data_directory.string(), "Germany", "pydata");
        BOOST_OUTCOME_TRY(
            set_node_function(params, start_date, end_date, pydata_dir,
                            mio::path_join(pydata_dir, "county_current_population.json"),
                            true, params_graph, read_function_nodes, node_id_function, scaling_factor_infected,
                            scaling_factor_icu, tnt_capacity_factor, 0, false, true));
        BOOST_OUTCOME_TRY(set_edge_function(data_dir, params_graph, mobile_compartments, contact_locations.size(),
                                            read_function_edges, std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.}));
        
        for (auto& node : param_graph.nodes()) {
            graph.add_node(node.id, node.property);
        }
        for (auto& edge : param_graph.edges()) {
            graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
        }

        return mio::success();
    }

    template <typename FP>
    mio::IOResult<void> set_initial_values(Graph<FP>& params_graph)
    {
        mio::osecirvvs::Parameters<FP> params(m_num_age_groups);

        BOOST_OUTCOME_TRY(set_covid_parameters<FP>(params));
        BOOST_OUTCOME_TRY(set_contact_matrices<FP>(params));
        //BOOST_OUTCOME_TRY(auto&& node_ids, node_func(population_data_path, is_node_for_county, rki_age_groups));
        // BOOST_OUTCOME_TRY(set_synthetic_population_data(model));
        if (m_use_county){
            BOOST_OUTCOME_TRY(set_population_data_germany<FP>(params_graph, params));
        }
        else {
            BOOST_OUTCOME_TRY(set_population_data_county<FP>(params_graph, params));
        }
        // model.apply_constraints();

        return mio::success();
    }

    boost::filesystem::path m_data_directory;
    double m_t0;
    double m_tmax;
    int m_num_age_groups = 6;
    bool m_use_county = true;
    std::vector<std::tuple<std::string, Eigen:MatrixXd, Eigen:MatrixXd>> m_contact_matrices{};
};

template <class ContactLocation>
mio::IOResult<void> add_contact_matrices_from_data(OptimizationModel<ContactLocation>& model_creator, std::map<ContactLocation, std::string> contact_locations)
{
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                        mio::read_mobility_plain(
                            (model_creator.data_directory() / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        model_creator.add_contact_location(contact_locations.second, baseline);
    }
}
