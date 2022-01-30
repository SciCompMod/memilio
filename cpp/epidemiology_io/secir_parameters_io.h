/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kuehn, Lena Ploetzke
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
#ifndef SAVE_PARAMETERS_H
#define SAVE_PARAMETERS_H

#include "epidemiology/utils/eigen_util.h"
#include "epidemiology/utils/graph.h"
#include "epidemiology/migration/migration.h"
#include "epidemiology/secir/parameter_studies.h"
#include "epidemiology/secir/parameter_studies_vaccinated.h"
#include "epidemiology_io/io.h"
#include "epidemiology_io/secir_result_io.h"
#include "epidemiology_io/json_serializer.h"
#include "epidemiology/utils/date.h"
#include "epidemiology/secir/analyze_result.h"

namespace epi
{

/**
 * @brief creates json files for each node in a simulation graph.
 * Creates two files per node: one contains the models and its parameters, one contains the outgoing edges.
 * @param graph Graph which should be written
 * @param directory directory where files should be stored
 * @param ioflags flags that set the behavior of serialization; see epi::IOFlags
 */
template <class Model>
IOResult<void> write_graph(const Graph<Model, MigrationParameters>& graph, const std::string& directory,
                           int ioflags = IOF_None)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");

    std::string abs_path;
    BOOST_OUTCOME_TRY(created, create_directory(directory, abs_path));

    if (created) {
        log_info("Results are stored in {:s}/results.", epi::get_current_dir_name());
    }
    else {
        log_info("Results are stored in {:s}/results. Files from previous "
                 "graph will be "
                 "overwritten",
                 epi::get_current_dir_name());
    }

    //write two files per node
    //one file that contains outgoing edges from the node
    //one file for the model (parameters and population)
    for (auto inode = size_t(0); inode < graph.nodes().size(); ++inode) {
        //node
        auto& node = graph.nodes()[inode];
        BOOST_OUTCOME_TRY(js_node_model, serialize_json(node.property, ioflags));
        Json::Value js_node(Json::objectValue);
        js_node["NodeId"]  = node.id;
        js_node["Model"]   = js_node_model;
        auto node_filename = path_join(abs_path, "GraphNode" + std::to_string(inode) + ".json");
        BOOST_OUTCOME_TRY(write_json(node_filename, js_node));

        //list of edges
        auto out_edges = graph.out_edges(inode);
        if (out_edges.size()) {
            Json::Value js_edges(Json::arrayValue);
            for (auto& e : graph.out_edges(inode)) {
                BOOST_OUTCOME_TRY(js_edge_params, serialize_json(e.property, ioflags));
                Json::Value js_edge{Json::objectValue};
                js_edge["StartNodeIndex"] = Json::UInt64(e.start_node_idx);
                js_edge["EndNodeIndex"]   = Json::UInt64(e.end_node_idx);
                js_edge["Parameters"]     = js_edge_params;
                js_edges.append(std::move(js_edge));
            }
            auto edge_filename = path_join(abs_path, "GraphEdges_node" + std::to_string(inode) + ".json");
            BOOST_OUTCOME_TRY(write_json(edge_filename, js_edges));
        }
    }

    return success();
}

/**
 * @brief reads graph json files and returns a simulation graph.
 * See write_graph for information of expected files.
 * @tparam the type of the simulation model.
 * @param directory directory from where graph should be read.
 */
template <class Model>
IOResult<Graph<Model, MigrationParameters>> read_graph(const std::string& directory, int ioflags = IOF_None)
{
    std::string abs_path;
    if (!file_exists(directory, abs_path)) {
        log_error("Directory {} does not exist.", directory);
        return failure(StatusCode::FileNotFound, directory);
    }

    auto graph = Graph<Model, MigrationParameters>{};

    //read nodes, as many as files are available
    for (auto inode = 0;; ++inode) {
        auto node_filename = path_join(abs_path, "GraphNode" + std::to_string(inode) + ".json");
        if (!file_exists(node_filename, node_filename)) {
            break;
        }
        BOOST_OUTCOME_TRY(js_node, read_json(node_filename));
        if (!js_node["NodeId"].isInt()) {
            log_error("NodeId field must be an integer.");
            return failure(StatusCode::InvalidType, node_filename + ", NodeId must be an integer.");
        }
        auto node_id = js_node["NodeId"].asInt();
        BOOST_OUTCOME_TRY(model, deserialize_json(js_node["Model"], Tag<Model>{}, ioflags));
        graph.add_node(node_id, model);
    }

    //read edges; nodes must already be available for that)
    for (auto inode = size_t(0); inode < graph.nodes().size(); ++inode) {
        //list of edges
        auto edge_filename = path_join(abs_path, "GraphEdges_node" + std::to_string(inode) + ".json");
        BOOST_OUTCOME_TRY(js_edges, read_json(edge_filename));

        for (auto& e : js_edges) {
            auto start_node_idx  = inode;
            auto js_end_node_idx = e["EndNodeIndex"];
            if (!js_end_node_idx.isUInt64()) {
                log_error("EndNodeIndex must be an integer.");
                return failure(StatusCode::InvalidType, edge_filename + ", EndNodeIndex must be an integer.");
            }
            auto end_node_idx = js_end_node_idx.asUInt64();
            if (end_node_idx >= graph.nodes().size()) {
                log_error("EndNodeIndex not in range of number of graph nodes.");
                return failure(StatusCode::OutOfRange,
                               edge_filename + ", EndNodeIndex not in range of number of graph nodes.");
            }
            BOOST_OUTCOME_TRY(parameters, deserialize_json(e["Parameters"], Tag<MigrationParameters>{}, ioflags));
            graph.add_edge(start_node_idx, end_node_idx, parameters);
        }
    }

    return success(graph);
}

/**
 * @brief creates xml file with a single run parameter study with std 0 (used to save parameters of individual runs)
 * @param filename Name of file
 * @param params Secir parameters used during run
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
template <class Simulation>
IOResult<void> write_single_run_result(const int run,
                                       const epi::Graph<epi::SimulationNode<Simulation>, epi::MigrationEdge>& graph)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");

    std::string abs_path;
    BOOST_OUTCOME_TRY(created, create_directory("results", abs_path));

    if (created) {
        log_info("Results are stored in {:s}/results.", epi::get_current_dir_name());
    }
    else if (run == 0) {
        log_info(
            "Directory '{:s}' already exists. Results are stored in {:s}/results. Files from previous runs will be "
            "overwritten",
            epi::get_current_dir_name());
    }

    //write sampled parameters for this run
    //omit edges to save space as they are not sampled
    for (auto inode = size_t(0); inode < graph.nodes().size(); ++inode) {
        auto& node = graph.nodes()[inode];
        BOOST_OUTCOME_TRY(js_node_model, serialize_json(node.property.get_result(), epi::IOF_OmitDistributions));
        Json::Value js_node(Json::objectValue);
        js_node["NodeId"] = node.id;
        js_node["Model"]  = js_node_model;
        auto node_filename =
            path_join(abs_path, "Parameters_run" + std::to_string(run) + "_node" + std::to_string(inode) + ".json");
        BOOST_OUTCOME_TRY(write_json(node_filename, js_node));
    }

    //write results for this run
    std::vector<TimeSeries<double>> all_results;
    std::vector<int> ids;

    ids.reserve(graph.nodes().size());
    all_results.reserve(graph.nodes().size());
    std::transform(graph.nodes().begin(), graph.nodes().end(), std::back_inserter(all_results), [](auto& node) {
        return node.property.get_result();
    });
    std::transform(graph.nodes().begin(), graph.nodes().end(), std::back_inserter(ids), [](auto& node) {
        return node.id;
    });
    BOOST_OUTCOME_TRY(
        save_result(all_results, ids, path_join(abs_path, ("Results_run" + std::to_string(run) + ".h5"))));

    return success();
}

namespace details
{
    /**
     * @brief interpolates age_ranges to param_ranges and saves ratios in interpolation
     * @param age_ranges original age ranges of the data
     * @param interpolation vector of ratios that are aplied to the data of age_ranges
     * @param carry_over boolean vector which indicates whether there is an overflow from one age group to the next while interpolating data
     */
    void interpolate_ages(const std::vector<double>& age_ranges, std::vector<std::vector<double>>& interpolation,
                          std::vector<bool>& carry_over);

    /**
     * @brief reads populations data from RKI
     * @param path Path to RKI file
     * @param id_name Name of region key column
     * @param date Date for which the arrays are initialized
     * @param region vector of keys of the region of interest
     * @param year Specifies year at which the data is read
     * @param month Specifies month at which the data is read
     * @param day Specifies day at which the data is read
     * @param num_* output vector for number of people in the corresponding compartement
     * @param t_* vector average time it takes to get from one compartement to another for each age group
     * @param mu_* vector probabilities to get from one compartement to another for each age group
     */
    IOResult<void>
    read_rki_data(std::string const& path, const std::string& id_name, std::vector<int> const& region, Date date,
                  std::vector<std::vector<double>>& num_exp, std::vector<std::vector<double>>& num_car,
                  std::vector<std::vector<double>>& num_inf, std::vector<std::vector<double>>& num_hosp,
                  std::vector<std::vector<double>>& num_icu, std::vector<std::vector<double>>& num_death,
                  std::vector<std::vector<double>>& num_rec, const std::vector<std::vector<int>>& t_car_to_rec,
                  const std::vector<std::vector<int>>& t_car_to_inf, const std::vector<std::vector<int>>& t_exp_to_car,
                  const std::vector<std::vector<int>>& t_inf_to_rec, const std::vector<std::vector<int>>& t_inf_to_hosp,
                  const std::vector<std::vector<int>>& t_hosp_to_rec,
                  const std::vector<std::vector<int>>& t_hosp_to_icu,
                  const std::vector<std::vector<int>>& t_icu_to_dead, const std::vector<std::vector<int>>& t_icu_to_rec,
                  const std::vector<std::vector<double>>& mu_C_R, const std::vector<std::vector<double>>& mu_I_H,
                  const std::vector<std::vector<double>>& mu_H_U, const std::vector<std::vector<double>>& mu_U_D,
                  const std::vector<double>& scaling_factor_inf);

    IOResult<void> read_rki_recovered_data(std::string const& path, const std::string& id_name,
                                           std::vector<int> const& vregion, Date date,
                                           std::vector<std::vector<double>>& vnum_rec);

    /**
     * @brief sets populations data from RKI into a SecirModel
     * @param model vector of objects in which the data is set
     * @param path Path to RKI file
     * @param id_name Name of region key column
     * @param region vector of keys of the region of interest
     * @param year Specifies year at which the data is read
     * @param month Specifies month at which the data is read
     * @param day Specifies day at which the data is read
     * @param scaling_factor_inf factors by which to scale the confirmed cases of
     * rki data
     */
    template <class Model, class ModelType = InfectionState>
    IOResult<void> set_rki_data(std::vector<Model>& model, const std::string& path, const std::string& id_name,
                                std::vector<int> const& region, Date date,
                                const std::vector<double>& scaling_factor_inf)
    {

        std::vector<double> age_ranges = {5., 10., 20., 25., 20., 20.};
        assert(scaling_factor_inf.size() == age_ranges.size());

        std::vector<std::vector<int>> t_car_to_rec{model.size()}; // R9
        std::vector<std::vector<int>> t_car_to_inf{model.size()}; // R3
        std::vector<std::vector<int>> t_exp_to_car{model.size()}; // R2
        std::vector<std::vector<int>> t_inf_to_rec{model.size()}; // R4
        std::vector<std::vector<int>> t_inf_to_hosp{model.size()}; // R6
        std::vector<std::vector<int>> t_hosp_to_rec{model.size()}; // R5
        std::vector<std::vector<int>> t_hosp_to_icu{model.size()}; // R7
        std::vector<std::vector<int>> t_icu_to_dead{model.size()}; // R10
        std::vector<std::vector<int>> t_icu_to_rec{model.size()};    

        std::vector<std::vector<double>> mu_C_R{model.size()};
        std::vector<std::vector<double>> mu_I_H{model.size()};
        std::vector<std::vector<double>> mu_H_U{model.size()};
        std::vector<std::vector<double>> mu_U_D{model.size()};        

        for (size_t county = 0; county < model.size(); county++) {
            for (size_t group = 0; group < age_ranges.size(); group++) {

                t_car_to_inf[county].push_back(static_cast<int>(
                    2 * (model[county].parameters.template get<epi::IncubationTime>()[(epi::AgeGroup)group] -
                         model[county].parameters.template get<epi::SerialInterval>()[(epi::AgeGroup)group])));
                t_car_to_rec[county].push_back(static_cast<int>(
                    t_car_to_inf[county][group] +
                    0.5 * model[county].parameters.template get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group]));
                t_exp_to_car[county].push_back(static_cast<int>(
                    2 * model[county].parameters.template get<epi::SerialInterval>()[(epi::AgeGroup)group] -
                    model[county].parameters.template get<epi::IncubationTime>()[(epi::AgeGroup)group]));
                t_inf_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group]));
                t_inf_to_hosp[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)group]));
                t_hosp_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)group]));
                t_hosp_to_icu[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)group]));
                t_icu_to_dead[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::ICUToDeathTime>()[(epi::AgeGroup)group]));
                t_icu_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::ICUToHomeTime>()[(epi::AgeGroup)group]));                    

                mu_C_R[county].push_back(
                    model[county].parameters.template get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)group]);
                mu_I_H[county].push_back(
                    model[county].parameters.template get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)group]);
                mu_H_U[county].push_back(
                    model[county].parameters.template get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)group]);
                mu_U_D[county].push_back(
                    model[county].parameters.template get<epi::DeathsPerICU>()[(epi::AgeGroup)group]);                    
            }
        }
        std::vector<std::vector<double>> num_inf(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_death(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_exp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_car(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_hosp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_icu(model.size(), std::vector<double>(age_ranges.size(), 0.0));

        BOOST_OUTCOME_TRY(read_rki_data(path, id_name, region, date, num_exp, num_car, num_inf, num_hosp, num_icu,
                                        num_death, num_rec, t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec,
                                        t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead, t_icu_to_rec, mu_C_R, mu_I_H,
                                        mu_H_U, mu_U_D, scaling_factor_inf));

        for (size_t county = 0; county < model.size(); county++) {
            if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
                size_t num_groups = (size_t)model[county].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[county].populations[{AgeGroup(i), ModelType::Exposed}]      = num_exp[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::Carrier}]      = num_car[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::Infected}]     = num_inf[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::Hospitalized}] = num_hosp[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::Dead}]         = num_death[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::Recovered}] =
                        num_rec[county][i] - model[county].populations[{AgeGroup(i), ModelType::ICU}];
                    model[county].populations[{AgeGroup(i), ModelType::InfTotal}] =
                        model[county].populations[{AgeGroup(i), ModelType::Recovered}] +
                        model[county].populations[{AgeGroup(i), ModelType::Dead}] +
                        model[county].populations[{AgeGroup(i), ModelType::Hospitalized}] +
                        model[county].populations[{AgeGroup(i), ModelType::ICU}];
                }
            }
            else {
                log_warning("No infections reported on date " + std::to_string(date.year) + "-" +
                            std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                            std::to_string(region[county]) + ". Population data has not been set.");
            }
        }
        return success();
    }

    /**
     * @brief sets populations data from RKI into a SecirModel
     * @param model vector of objects in which the data is set
     * @param path Path to RKI file
     * @param id_name Name of region key column
     * @param region vector of keys of the region of interest
     * @param year Specifies year at which the data is read
     * @param month Specifies month at which the data is read
     * @param day Specifies day at which the data is read
     * @param scaling_factor_inf factors by which to scale the confirmed cases of
     * rki data
     */
    template <class Model, class ModelType = InfectionStateV>
    IOResult<void> set_rki_data_vaccmodel(std::vector<Model>& model, const std::string& path, const std::string& id_name,
                                        std::vector<int> const& region, Date date,
                                        const std::vector<double>& scaling_factor_inf)
    {

        std::vector<double> age_ranges = {5., 10., 20., 25., 20., 20.};
        assert(scaling_factor_inf.size() == age_ranges.size());

        std::vector<std::vector<int>> t_car_to_rec{model.size()}; // R9
        std::vector<std::vector<int>> t_car_to_inf{model.size()}; // R3
        std::vector<std::vector<int>> t_exp_to_car{model.size()}; // R2
        std::vector<std::vector<int>> t_inf_to_rec{model.size()}; // R4
        std::vector<std::vector<int>> t_inf_to_hosp{model.size()}; // R6
        std::vector<std::vector<int>> t_hosp_to_rec{model.size()}; // R5
        std::vector<std::vector<int>> t_hosp_to_icu{model.size()}; // R7
        std::vector<std::vector<int>> t_icu_to_dead{model.size()}; // R10
        std::vector<std::vector<int>> t_icu_to_rec{model.size()};

        std::vector<std::vector<double>> mu_C_R{model.size()};
        std::vector<std::vector<double>> mu_I_H{model.size()};
        std::vector<std::vector<double>> mu_H_U{model.size()};
        std::vector<std::vector<double>> mu_U_D{model.size()};        

        std::vector<std::vector<double>> num_inf(model.size());
        std::vector<std::vector<double>> num_death(model.size());
        std::vector<std::vector<double>> num_rec(model.size());
        std::vector<std::vector<double>> num_exp(model.size());
        std::vector<std::vector<double>> num_car(model.size());
        std::vector<std::vector<double>> num_hosp(model.size());
        std::vector<std::vector<double>> num_icu(model.size());

        /*----------- UNVACCINATED -----------*/
        for (size_t county = 0; county < model.size(); county++) {
            for (size_t group = 0; group < age_ranges.size(); group++) {

                num_inf[county].push_back(0);
                num_death[county].push_back(0);
                num_rec[county].push_back(0);
                num_exp[county].push_back(0);
                num_car[county].push_back(0);
                num_hosp[county].push_back(0);
                num_icu[county].push_back(0);

                t_car_to_inf[county].push_back(static_cast<int>(
                    2 * (model[county].parameters.template get<epi::IncubationTime>()[(epi::AgeGroup)group] -
                         model[county].parameters.template get<epi::SerialInterval>()[(epi::AgeGroup)group])));
                t_car_to_rec[county].push_back(static_cast<int>(
                    t_car_to_inf[county][group] +
                    0.5 * model[county].parameters.template get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group]));
                t_exp_to_car[county].push_back(static_cast<int>(
                    2 * model[county].parameters.template get<epi::SerialInterval>()[(epi::AgeGroup)group] -
                    model[county].parameters.template get<epi::IncubationTime>()[(epi::AgeGroup)group]));
                t_inf_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group]));
                t_inf_to_hosp[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)group]));
                t_hosp_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)group]));
                t_hosp_to_icu[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)group]));
                t_icu_to_dead[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::ICUToDeathTime>()[(epi::AgeGroup)group]));
                t_icu_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::ICUToHomeTime>()[(epi::AgeGroup)group]));                    

                mu_C_R[county].push_back(
                    model[county].parameters.template get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)group]);
                mu_I_H[county].push_back(
                    model[county].parameters.template get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)group]);
                mu_H_U[county].push_back(
                    model[county].parameters.template get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)group]);
                mu_U_D[county].push_back(
                    model[county].parameters.template get<epi::DeathsPerICU>()[(epi::AgeGroup)group]);                         
            }
        }

        BOOST_OUTCOME_TRY(read_rki_data(path, id_name, region, date, num_exp, num_car, num_inf, num_hosp, num_icu,
                                        num_death, num_rec, t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec,
                                        t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead, t_icu_to_rec, mu_C_R, mu_I_H,
                                        mu_H_U, mu_U_D, scaling_factor_inf));

        for (size_t county = 0; county < model.size(); county++) {
            if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
                size_t num_groups = (size_t)model[county].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[county].populations[{AgeGroup(i), ModelType::Exposed}]      = num_exp[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::Carrier}]      = num_car[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::Infected}]     = num_inf[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::Hospitalized}] = num_hosp[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::Recovered}]    = num_rec[county][i];
                }
            }
            else {
                log_warning("No infections reported on date " + std::to_string(date.year) + "-" +
                            std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                            std::to_string(region[county]) + ". Population data has not been set.");
            }
        }

        /*----------- PARTIALLY VACCINATED -----------*/
        for (size_t county = 0; county < model.size(); county++) {
            t_car_to_rec[county].clear();
            t_car_to_inf[county].clear();
            t_exp_to_car[county].clear();
            t_inf_to_rec[county].clear();
            t_inf_to_hosp[county].clear();
            t_hosp_to_rec[county].clear();
            t_hosp_to_icu[county].clear();
            t_icu_to_dead[county].clear();
            t_icu_to_rec[county].clear();

            mu_C_R[county].clear();
            mu_I_H[county].clear();
            mu_H_U[county].clear();
            mu_U_D[county].clear();            

            num_inf[county].clear();
            num_death[county].clear();
            num_rec[county].clear();
            num_exp[county].clear();
            num_car[county].clear();
            num_hosp[county].clear();
            num_icu[county].clear();
            for (size_t group = 0; group < age_ranges.size(); group++) {

                num_inf[county].push_back(0);
                num_death[county].push_back(0);
                num_rec[county].push_back(0);
                num_exp[county].push_back(0);
                num_car[county].push_back(0);
                num_hosp[county].push_back(0);
                num_icu[county].push_back(0);

                double reduc_t = model[0].parameters.template get<epi::ReducTime>()[(epi::AgeGroup)group];
                t_car_to_inf[county].push_back(static_cast<int>(
                    2 * (model[county].parameters.template get<epi::IncubationTime>()[(epi::AgeGroup)group] -
                         model[county].parameters.template get<epi::SerialInterval>()[(epi::AgeGroup)group])));
                t_car_to_rec[county].push_back(static_cast<int>(
                    t_car_to_inf[county][group] +
                    0.5 * model[county].parameters.template get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group] *
                        reduc_t));
                t_exp_to_car[county].push_back(static_cast<int>(
                    2 * model[county].parameters.template get<epi::SerialInterval>()[(epi::AgeGroup)group] -
                    model[county].parameters.template get<epi::IncubationTime>()[(epi::AgeGroup)group]));
                t_inf_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group] * reduc_t));
                t_inf_to_hosp[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)group]));
                t_hosp_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)group]));
                t_hosp_to_icu[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)group]));
                t_icu_to_dead[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::ICUToDeathTime>()[(epi::AgeGroup)group]));
                t_icu_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::ICUToHomeTime>()[(epi::AgeGroup)group]));

                mu_C_R[county].push_back(
                    (1 - model[county].parameters.template get<epi::ReducExpInf>()[(epi::AgeGroup)group] *
                             (1 - model[county].parameters.template get<epi::AsymptoticCasesPerInfectious>()[(
                                      epi::AgeGroup)group])));
                mu_I_H[county].push_back(
                    model[county].parameters.template get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)group] *
                    model[county].parameters.template get<epi::ReducInfHosp>()[(epi::AgeGroup)group]);
                mu_H_U[county].push_back(
                    model[county].parameters.template get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)group] *
                    model[county].parameters.template get<epi::ReducInfHosp>()[(epi::AgeGroup)group]);
                mu_U_D[county].push_back(
                    model[county].parameters.template get<epi::DeathsPerICU>()[(epi::AgeGroup)group]);                       
            }
        }

        BOOST_OUTCOME_TRY(read_rki_data(path, id_name, region, date, num_exp, num_car, num_inf, num_hosp, num_icu,
                                        num_death, num_rec, t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec,
                                        t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead, t_icu_to_rec, mu_C_R, mu_I_H,
                                        mu_H_U, mu_U_D, scaling_factor_inf));

        for (size_t county = 0; county < model.size(); county++) {
            if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
                size_t num_groups = (size_t)model[county].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[county].populations[{AgeGroup(i), ModelType::ExposedV1}]      = num_exp[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::CarrierV1}]      = num_car[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::InfectedV1}]     = num_inf[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::HospitalizedV1}] = num_hosp[county][i];
                }
            }
            else {
                log_warning("No infections reported on date " + std::to_string(date.year) + "-" +
                            std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                            std::to_string(region[county]) + ". Population data has not been set.");
            }
        }

        /*----------- FULLY VACCINATED -----------*/
        for (size_t county = 0; county < model.size(); county++) {
            t_car_to_rec[county].clear();
            t_car_to_inf[county].clear();
            t_exp_to_car[county].clear();
            t_inf_to_rec[county].clear();
            t_inf_to_hosp[county].clear();
            t_hosp_to_rec[county].clear();
            t_hosp_to_icu[county].clear();
            t_icu_to_dead[county].clear();
            t_icu_to_rec[county].clear();

            mu_C_R[county].clear();
            mu_I_H[county].clear();
            mu_H_U[county].clear();
            mu_U_D[county].clear();

            num_inf[county].clear();
            num_death[county].clear();
            num_rec[county].clear();
            num_exp[county].clear();
            num_car[county].clear();
            num_hosp[county].clear();
            num_icu[county].clear();
            for (size_t group = 0; group < age_ranges.size(); group++) {

                num_inf[county].push_back(0);
                num_death[county].push_back(0);
                num_rec[county].push_back(0);
                num_exp[county].push_back(0);
                num_car[county].push_back(0);
                num_hosp[county].push_back(0);
                num_icu[county].push_back(0);

                double reduc_t = model[0].parameters.template get<epi::ReducTime>()[(epi::AgeGroup)group];
                t_car_to_inf[county].push_back(static_cast<int>(
                    2 * (model[county].parameters.template get<epi::IncubationTime>()[(epi::AgeGroup)group] -
                         model[county].parameters.template get<epi::SerialInterval>()[(epi::AgeGroup)group])));
                t_car_to_rec[county].push_back(static_cast<int>(
                    t_car_to_inf[county][group] +
                    0.5 * model[county].parameters.template get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group] *
                        reduc_t));
                t_exp_to_car[county].push_back(static_cast<int>(
                    2 * model[county].parameters.template get<epi::SerialInterval>()[(epi::AgeGroup)group] -
                    model[county].parameters.template get<epi::IncubationTime>()[(epi::AgeGroup)group]));
                t_inf_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group] * reduc_t));
                t_inf_to_hosp[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)group]));
                t_hosp_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)group]));
                t_hosp_to_icu[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)group]));
                t_icu_to_dead[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::ICUToDeathTime>()[(epi::AgeGroup)group]));
                t_icu_to_rec[county].push_back(static_cast<int>(
                    model[county].parameters.template get<epi::ICUToHomeTime>()[(epi::AgeGroup)group]));                    

                mu_C_R[county].push_back(
                    (1 - model[county].parameters.template get<epi::ReducImmuneExpInf>()[(epi::AgeGroup)group] *
                             (1 - model[county].parameters.template get<epi::AsymptoticCasesPerInfectious>()[(
                                      epi::AgeGroup)group])));
                mu_I_H[county].push_back(
                    model[county].parameters.template get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)group] *
                    model[county].parameters.template get<epi::ReducImmuneInfHosp>()[(epi::AgeGroup)group]);
                mu_H_U[county].push_back(
                    model[county].parameters.template get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)group] *
                    model[county].parameters.template get<epi::ReducImmuneInfHosp>()[(epi::AgeGroup)group]);
                mu_U_D[county].push_back(
                    model[county].parameters.template get<epi::DeathsPerICU>()[(epi::AgeGroup)group]);                       
            }
        }

        BOOST_OUTCOME_TRY(read_rki_data(path, id_name, region, date, num_exp, num_car, num_inf, num_hosp, num_icu,
                                        num_death, num_rec, t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec,
                                        t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead, t_icu_to_rec, mu_C_R, mu_I_H,
                                        mu_H_U, mu_U_D, scaling_factor_inf));

        for (size_t county = 0; county < model.size(); county++) {
            if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
                size_t num_groups = (size_t)model[county].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[county].populations[{AgeGroup(i), ModelType::ExposedV2}]      = num_exp[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::CarrierV2}]      = num_car[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::InfectedV2}]     = num_inf[county][i];
                    model[county].populations[{AgeGroup(i), ModelType::HospitalizedV2}] = num_hosp[county][i];
                }
            }
            else {
                log_warning("No infections reported on date " + std::to_string(date.year) + "-" +
                            std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                            std::to_string(region[county]) + ". Population data has not been set.");
            }
        } // namespace details
        return success();
    } // namespace epi

    /**
     * @brief reads number of ICU patients from DIVI register into SecirParams
     * @param path Path to DIVI file
     * @param id_name Name of region key column
     * @param vregion Keys of the region of interest
     * @param year Specifies year at which the data is read
     * @param month Specifies month at which the data is read
     * @param day Specifies day at which the data is read
     * @param vnum_icu number of ICU patients
     */
    IOResult<void> read_divi_data(const std::string& path, const std::string& id_name, const std::vector<int>& vregion,
                                  Date date, std::vector<double>& vnum_icu);

    /**
     * @brief sets populations data from DIVI register into Model
     * @param model vector of objects in which the data is set
     * @param path Path to DIVI file
     * @param id_name Name of region key column
     * @param vregion vector of keys of the regions of interest
     * @param year Specifies year at which the data is read
     * @param month Specifies month at which the data is read
     * @param day Specifies day at which the data is read
     * @param scaling_factor_icu factor by which to scale the icu cases of divi data
     */
    template <class Model, class ModelType = InfectionState>
    IOResult<void> set_divi_data(std::vector<Model>& model, const std::string& path, const std::string& id_name,
                                 const std::vector<int>& vregion, Date date, double scaling_factor_icu)
    {
        std::vector<double> sum_mu_I_U(vregion.size(), 0);
        std::vector<std::vector<double>> mu_I_U{model.size()};
        for (size_t region = 0; region < vregion.size(); region++) {
            auto num_groups = model[region].parameters.get_num_groups();
            for (auto i = epi::AgeGroup(0); i < num_groups; i++) {
                sum_mu_I_U[region] += model[region].parameters.template get<epi::ICUCasesPerHospitalized>()[i] *
                                      model[region].parameters.template get<epi::HospitalizedCasesPerInfectious>()[i];
                mu_I_U[region].push_back(
                    model[region].parameters.template get<epi::ICUCasesPerHospitalized>()[i] *
                    model[region].parameters.template get<epi::HospitalizedCasesPerInfectious>()[i]);
            }
        }
        std::vector<double> num_icu(model.size(), 0.0);
        BOOST_OUTCOME_TRY(read_divi_data(path, id_name, vregion, date, num_icu));

        for (size_t region = 0; region < vregion.size(); region++) {
            auto num_groups = model[region].parameters.get_num_groups();
            for (auto i = epi::AgeGroup(0); i < num_groups; i++) {
                model[region].populations[{i, ModelType::ICU}] =
                    scaling_factor_icu * num_icu[region] * mu_I_U[region][(size_t)i] / sum_mu_I_U[region];
            }
        }

        return success();
    }

    /**
     * @brief reads population data from census data
     * @param path Path to RKI file
     * @param id_name Name of region key column
     * @param vregion vector of keys of the regions of interest
     */
    IOResult<std::vector<std::vector<double>>> read_population_data(const std::string& path, const std::string& id_name,
                                                                    const std::vector<int>& vregion);

    /**
     * @brief sets population data from census data
     * @param model vector of objects in which the data is set
     * @param path Path to RKI file
     * @param id_name Name of region key column
     * @param vregion vector of keys of the regions of interest
     */

    template <class Model, class ModelType = InfectionState>
    IOResult<void> set_population_data(std::vector<Model>& model, const std::string& path, const std::string& id_name,
                                       const std::vector<int>& vregion)
    {
        BOOST_OUTCOME_TRY(num_population, read_population_data(path, id_name, vregion));

        for (size_t region = 0; region < vregion.size(); region++) {
            if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
                auto num_groups = model[region].parameters.get_num_groups();
                for (auto i = AgeGroup(0); i < num_groups; i++) {
                    model[region].populations.template set_difference_from_group_total<epi::AgeGroup>(
                        {i, ModelType::Susceptible}, num_population[region][size_t(i)]);
                }
            }
            else {
                log_warning("No population data available for region " + std::to_string(region) +
                            ". Population data has not been set.");
            }
        }

        return success();
    }

    template <class Model>
    IOResult<void> set_population_data_vaccmodel(std::vector<Model>& model, const std::string& path,
                                               const std::string& path_rki, const std::string& id_name,
                                               const std::vector<int>& vregion, Date date)
    {

        BOOST_OUTCOME_TRY(num_population, read_population_data(path, id_name, vregion));

        std::vector<double> age_ranges = {5., 10., 20., 25., 20., 20.};
        std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(age_ranges.size(), 0.0));

        BOOST_OUTCOME_TRY(read_rki_recovered_data(path_rki, id_name, vregion, date, num_rec));

        for (size_t region = 0; region < vregion.size(); region++) {
            if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
                auto num_groups = model[region].parameters.get_num_groups();
                for (auto i = AgeGroup(0); i < num_groups; i++) {

                    double S_v = model[region].parameters.template get<DailyFullVaccination>()[i][0] +
                                 num_rec[region][size_t(i)];
                    double S_pv = model[region].parameters.template get<DailyFirstVaccination>()[i][0] -
                                  model[region].parameters.template get<DailyFullVaccination>()[i][0]; // use std::min with 0
                    double S;
                    if (num_population[region][size_t(i)] - S_pv - S_v < 0.0) {
                        std::cout << "county: " << region << ", AgeGroup: " << size_t(i) << std::endl;
                        std::cout << "num pop: " << num_population[region][size_t(i)] << std::endl;
                        std::cout << "S_pv: " << S_pv << std::endl;
                        std::cout << "S_v (V,rec): " << S_v << " ("
                                  << model[region].parameters.template get<DailyFullVaccination>()[i][0] << ", "
                                  << num_rec[region][size_t(i)] << ")" << std::endl;
                        std::cout << std::endl;
                        S = 0.0;
                        S_v += num_population[region][size_t(i)] - S_pv - S_v;
                    }
                    else {
                        S = num_population[region][size_t(i)] - S_pv - S_v;
                    }

                    double ratio_E  = 1 / (S + S_pv * model[region].parameters.template get<ReducVaccExp>()[i] +
                                          S_v * model[region].parameters.template get<ReducImmuneExp>()[i]);
                    double ratio_C  = 1 / (S + S_pv + S_v);
                    double ratio_I  = 1 / (S + S_pv * model[region].parameters.template get<ReducExpInf>()[i] +
                                          S_v * model[region].parameters.template get<ReducImmuneExpInf>()[i]);
                    double ratio_HU = 1 / (S + S_pv * model[region].parameters.template get<ReducInfHosp>()[i] +
                                           S_v * model[region].parameters.template get<ReducImmuneInfHosp>()[i]);

                    model[region].populations[{i, InfectionStateV::Exposed}] =
                        S * ratio_E * model[region].populations[{i, InfectionStateV::Exposed}];
                    model[region].populations[{i, InfectionStateV::ExposedV1}] =
                        S_pv * model[region].parameters.template get<ReducVaccExp>()[i] * ratio_E *
                        model[region].populations[{i, InfectionStateV::ExposedV1}];
                    model[region].populations[{i, InfectionStateV::ExposedV2}] =
                        S_v * model[region].parameters.template get<ReducImmuneExp>()[i] * ratio_E *
                        model[region].populations[{i, InfectionStateV::ExposedV2}];

                    model[region].populations[{i, InfectionStateV::Carrier}] =
                        S * ratio_C * model[region].populations[{i, InfectionStateV::Carrier}];
                    model[region].populations[{i, InfectionStateV::CarrierV1}] =
                        S_pv * ratio_C * model[region].populations[{i, InfectionStateV::CarrierV1}];
                    model[region].populations[{i, InfectionStateV::CarrierV2}] =
                        S_v * ratio_C * model[region].populations[{i, InfectionStateV::CarrierV2}];

                    model[region].populations[{i, InfectionStateV::Infected}] =
                        S * ratio_I * model[region].populations[{i, InfectionStateV::Infected}];
                    model[region].populations[{i, InfectionStateV::InfectedV1}] =
                        S_pv * model[region].parameters.template get<ReducExpInf>()[i] * ratio_I *
                        model[region].populations[{i, InfectionStateV::InfectedV1}];
                    model[region].populations[{i, InfectionStateV::InfectedV2}] =
                        S_v * model[region].parameters.template get<ReducImmuneExpInf>()[i] * ratio_I *
                        model[region].populations[{i, InfectionStateV::InfectedV2}];

                    model[region].populations[{i, InfectionStateV::Hospitalized}] =
                        S * ratio_HU * model[region].populations[{i, InfectionStateV::Hospitalized}];
                    model[region].populations[{i, InfectionStateV::HospitalizedV1}] =
                        S_pv * model[region].parameters.template get<ReducInfHosp>()[i] * ratio_HU *
                        model[region].populations[{i, InfectionStateV::HospitalizedV1}];
                    model[region].populations[{i, InfectionStateV::HospitalizedV2}] =
                        S_v * model[region].parameters.template get<ReducImmuneInfHosp>()[i] * ratio_HU *
                        model[region].populations[{i, InfectionStateV::HospitalizedV2}];

                    model[region].populations[{i, InfectionStateV::ICUV1}] =
                        S_pv * model[region].parameters.template get<ReducInfHosp>()[i] * ratio_HU *
                        model[region].populations[{i, InfectionStateV::ICU}];
                    model[region].populations[{i, InfectionStateV::ICUV2}] =
                        S_v * model[region].parameters.template get<ReducImmuneInfHosp>()[i] * ratio_HU *
                        model[region].populations[{i, InfectionStateV::ICU}];
                    model[region].populations[{i, InfectionStateV::ICU}] =
                        S * ratio_HU * model[region].populations[{i, InfectionStateV::ICU}];

                    model[region].populations[{i, InfectionStateV::Recovered}] =
                        model[region].parameters.template get<DailyFullVaccination>()[i][0] +
                        model[region].populations[{i, InfectionStateV::Recovered}] -
                        (model[region].populations[{i, InfectionStateV::Infected}] +
                         model[region].populations[{i, InfectionStateV::InfectedV1}] +
                         model[region].populations[{i, InfectionStateV::InfectedV2}] +
                         model[region].populations[{i, InfectionStateV::Hospitalized}] +
                         model[region].populations[{i, InfectionStateV::HospitalizedV1}] +
                         model[region].populations[{i, InfectionStateV::HospitalizedV2}] +
                         model[region].populations[{i, InfectionStateV::ICU}] +
                         model[region].populations[{i, InfectionStateV::ICUV1}] +
                         model[region].populations[{i, InfectionStateV::ICUV1}]) -
                        model[region].populations[{i, InfectionStateV::Dead}];

                    model[region].populations[{i, InfectionStateV::Recovered}] =
                        std::min(S + S_pv + S_v,
                                 std::max(0.0, double(model[region].populations[{i, InfectionStateV::Recovered}])));

                    model[region].populations[{i, InfectionStateV::SusceptibleV1}] =
                        S_pv - model[region].populations[{i, InfectionStateV::ExposedV1}] -
                        model[region].populations[{i, InfectionStateV::InfectedV1}] -
                        model[region].populations[{i, InfectionStateV::HospitalizedV1}] -
                        model[region].populations[{i, InfectionStateV::ICUV1}];

                    model[region].populations.template set_difference_from_group_total<epi::AgeGroup>(
                        {i, InfectionStateV::Susceptible}, num_population[region][size_t(i)]);
                }

                for (auto i = epi::AgeGroup(0); i < AgeGroup(6); i++) {
                    for (auto j = epi::Index<epi::InfectionStateV>(0); j < epi::InfectionStateV::Count; ++j) {
                        if (model[region].populations[{i, j}] < 0) {
                            std::cout << "populations(" << i << ", " << j << ") is negative:"
                                      << model[region].populations[{i, j}] / num_population[region][size_t(i)]
                                      << std::endl;
                        }
                    }
                }
            }
            else {
                log_warning("No population data available for region " + std::to_string(region) +
                            ". Population data has not been set.");
            }
        }

        return success();
    }

    void set_vaccine_data(std::vector<SecirModelV>& model, const std::string& path, const std::string& id_name,
                          const std::vector<int>& vregion);

    void get_vaccine_growth(std::vector<SecirModelV>& model, const std::string& path, int window);

    void set_new_vaccine_data(std::vector<SecirModelV>& model, const std::string& path, const std::string& id_name,
                              const std::vector<int>& vregion);

    void get_new_vaccine_growth(std::vector<SecirModelV>& model, const std::string& path, Date date,
                                const std::string& id_name, const std::vector<int>& vregion, int num_days);
} //namespace details

/**
* @brief Sets the RKI data as in the model initialization step for a 
    arbitrary number of days and exports the corresponding time series
    to compare the simulation output against it.
* @param model vector of objects in which the data is set
* @param data_dir Path to RKI files
* @param results_dir Path to result files
* @param id_name Name of region key column
* @param region vector of keys of the region of interest
* @param year Specifies year at which the data is read
* @param month Specifies month at which the data is read
* @param day Specifies day at which the data is read
* @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
*/
template <class Model, class ModelType = InfectionState>
IOResult<void> set_rki_data_export_timeseries(std::vector<Model>& model, const std::string& data_dir,
                                       const std::string& results_dir, std::vector<int> const& region, Date date,
                                       const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                       int num_days)
{

    std::string id_name            = "ID_County";
    std::vector<double> age_ranges = {5., 10., 20., 25., 20., 20.};
    assert(scaling_factor_inf.size() == age_ranges.size());

    std::vector<std::vector<int>> t_car_to_rec{model.size()}; // R9
    std::vector<std::vector<int>> t_car_to_inf{model.size()}; // R3
    std::vector<std::vector<int>> t_exp_to_car{model.size()}; // R2
    std::vector<std::vector<int>> t_inf_to_rec{model.size()}; // R4
    std::vector<std::vector<int>> t_inf_to_hosp{model.size()}; // R6
    std::vector<std::vector<int>> t_hosp_to_rec{model.size()}; // R5
    std::vector<std::vector<int>> t_hosp_to_icu{model.size()}; // R7
    std::vector<std::vector<int>> t_icu_to_dead{model.size()}; // R10
    std::vector<std::vector<int>> t_icu_to_rec{model.size()}; 

    std::vector<std::vector<double>> mu_C_R{model.size()};
    std::vector<std::vector<double>> mu_I_H{model.size()};
    std::vector<std::vector<double>> mu_H_U{model.size()};
    std::vector<std::vector<double>> mu_U_D{model.size()};    

    std::vector<double> sum_mu_I_U(region.size(), 0);
    std::vector<std::vector<double>> mu_I_U{model.size()};

    for (size_t county = 0; county < model.size(); county++) {
        for (size_t group = 0; group < age_ranges.size(); group++) {

            t_car_to_inf[county].push_back(
                static_cast<int>(2 * (model[county].parameters.template get<IncubationTime>()[AgeGroup(group)] -
                                      model[county].parameters.template get<SerialInterval>()[AgeGroup(group)])));
            t_car_to_rec[county].push_back(
                static_cast<int>(t_car_to_inf[county][group] +
                                 0.5 * model[county].parameters.template get<InfectiousTimeMild>()[AgeGroup(group)]));
            t_exp_to_car[county].push_back(
                static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[AgeGroup(group)] -
                                 model[county].parameters.template get<IncubationTime>()[AgeGroup(group)]));
            t_inf_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<InfectiousTimeMild>()[AgeGroup(group)]));
            t_inf_to_hosp[county].push_back(
                static_cast<int>(model[county].parameters.template get<HomeToHospitalizedTime>()[AgeGroup(group)]));
            t_hosp_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<HospitalizedToHomeTime>()[AgeGroup(group)]));
            t_hosp_to_icu[county].push_back(
                static_cast<int>(model[county].parameters.template get<HospitalizedToICUTime>()[AgeGroup(group)]));
            t_icu_to_dead[county].push_back(
                static_cast<int>(model[county].parameters.template get<ICUToDeathTime>()[AgeGroup(group)]));
            t_icu_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<epi::ICUToHomeTime>()[AgeGroup(group)]));

            mu_C_R[county].push_back(
                model[county].parameters.template get<AsymptoticCasesPerInfectious>()[AgeGroup(group)]);
            mu_I_H[county].push_back(
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)]);
            mu_H_U[county].push_back(model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)]);

            sum_mu_I_U[county] +=
                model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)];
            mu_I_U[county].push_back(
                model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)]);
            mu_U_D[county].push_back(
                model[county].parameters.template get<epi::DeathsPerICU>()[AgeGroup(group)]);                  
        }
    }

    std::vector<TimeSeries<double>> rki_data(
        region.size(), TimeSeries<double>::zero(num_days, (size_t)ModelType::Count * age_ranges.size()));

    for (size_t j = 0; j < static_cast<size_t>(num_days); j++) {
        std::vector<std::vector<double>> num_inf(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_death(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_exp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_car(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_hosp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> dummy_icu(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<double> num_icu(model.size(), 0.0);

        BOOST_OUTCOME_TRY(details::read_rki_data(
            path_join(data_dir, "all_county_age_ma7_rki.json"), id_name, region, date, num_exp, num_car, num_inf,
            num_hosp, dummy_icu, num_death, num_rec, t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec,
            t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead, t_icu_to_rec, mu_C_R, mu_I_H, mu_H_U, mu_U_D, scaling_factor_inf));
        BOOST_OUTCOME_TRY(
            details::read_divi_data(path_join(data_dir, "county_divi_ma7.json"), id_name, region, date, num_icu));
        BOOST_OUTCOME_TRY(num_population, details::read_population_data(
                                              path_join(data_dir, "county_current_population.json"), id_name, region));

        for (size_t i = 0; i < region.size(); i++) {
            for (size_t age = 0; age < age_ranges.size(); age++) {
                rki_data[i][j]((size_t)ModelType::Exposed + (size_t)ModelType::Count * age)      = num_exp[i][age];
                rki_data[i][j]((size_t)ModelType::Carrier + (size_t)ModelType::Count * age)      = num_car[i][age];
                rki_data[i][j]((size_t)ModelType::Infected + (size_t)ModelType::Count * age)     = num_inf[i][age];
                rki_data[i][j]((size_t)ModelType::Hospitalized + (size_t)ModelType::Count * age) = num_hosp[i][age];
                rki_data[i][j]((size_t)ModelType::ICU + (size_t)ModelType::Count * age) =
                    scaling_factor_icu * num_icu[i] * mu_I_U[i][age] / sum_mu_I_U[i];
                rki_data[i][j]((size_t)ModelType::Recovered + (size_t)ModelType::Count * age) = num_rec[i][age];
                rki_data[i][j]((size_t)ModelType::Dead + (size_t)ModelType::Count * age)      = num_death[i][age];
                rki_data[i][j]((size_t)ModelType::Susceptible + (size_t)ModelType::Count * age) =
                    num_population[i][age] - num_exp[i][age] - num_car[i][age] - num_inf[i][age] - num_hosp[i][age] -
                    num_rec[i][age] - num_death[i][age] -
                    rki_data[i][j]((size_t)ModelType::ICU + (size_t)ModelType::Count * age);
            }
        }
        std::cout << "extrapolated rki data for date: " << date.day << "." << date.month << "." << date.year
                  << std::endl;
        date = offset_date_by_days(date, 1);
    }
    BOOST_OUTCOME_TRY(save_result<ModelType>(rki_data, region, path_join(results_dir, "Results_rki.h5")));

    auto rki_data_sum = epi::sum_nodes(std::vector<std::vector<TimeSeries<double>>>{rki_data});
    BOOST_OUTCOME_TRY(save_result<ModelType>({rki_data_sum[0][0]}, {0}, path_join(results_dir, "Results_rki_sum.h5")));

    return success();
}

/**
 * @brief reads population data from population files for the whole country
 * @param model vector of model in which the data is set
 * @param date Date for which the data should be read
 * @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
 * @param scaling_factor_icu factor by which to scale the icu cases of divi data
 * @param dir directory of files
 */
template <class Model>
IOResult<void> read_population_data_germany(std::vector<Model>& model, Date date,
                                            const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                            const std::string& dir)
{
    std::string id_name;
    std::vector<int> region(1, 0);
    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY(
            details::set_divi_data(model, path_join(dir, "germany_divi.json"), id_name, {0}, date, scaling_factor_icu));
    }
    else {
        log_warning("No DIVI data available for this date");
    }
    BOOST_OUTCOME_TRY(
        details::set_rki_data(model, path_join(dir, "all_age_ma_rki.json"), id_name, {0}, date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(dir, "county_current_population.json"), "ID_County", {0}));
    return success();
}

/**
 * @brief reads population data from population files for the specefied state
 * @param model vector of model in which the data is set
 * @param date Date for which the data should be read
 * @param state vector of region keys of states of interest
 * @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
 * @param scaling_factor_icu factor by which to scale the icu cases of divi data
 * @param dir directory of files
 */
template <class Model>
IOResult<void> read_population_data_state(std::vector<Model>& model, Date date, std::vector<int>& state,
                                          const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                          const std::string& dir)
{
    std::string id_name = "ID_State";
    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY(
            details::set_divi_data(model, path_join(dir, "state_divi.json"), id_name, state, date, scaling_factor_icu));
    }
    else {
        log_warning("No DIVI data available for this date");
    }

    BOOST_OUTCOME_TRY(details::set_rki_data(model, path_join(dir, "all_state_age_ma_rki.json"), id_name, state, date,
                                            scaling_factor_inf));
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(dir, "county_current_population.json"), "ID_County", state));
    return success();
}

/**
 * @brief reads population data from population files for the specefied county
 * @param model vector of model in which the data is set
 * @param date Date for which the data should be read
 * @param county vector of region keys of counties of interest
 * @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
 * @param scaling_factor_icu factor by which to scale the icu cases of divi data
 * @param dir directory of files
 */
template <class Model, class ModelType = InfectionState>
IOResult<void> read_population_data_county(std::vector<Model>& model, Date date, const std::vector<int>& county,
                                           const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                           const std::string& dir)
{
    std::string id_name = "ID_County";

    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY((details::set_divi_data<Model, ModelType>(model, path_join(dir, "county_divi_ma7.json"),
                                                                    id_name, county, date, scaling_factor_icu)));
    }
    else {
        log_warning("No DIVI data available for this date");
    }
    BOOST_OUTCOME_TRY((details::set_rki_data<Model, ModelType>(model, path_join(dir, "all_county_age_ma7_rki.json"),
                                                               id_name, county, date, scaling_factor_inf)));
    BOOST_OUTCOME_TRY((details::set_rki_data_export_timeseries<Model, ModelType>(model, dir, dir, county, date, scaling_factor_inf,
                                                                 scaling_factor_icu, 5)));
    BOOST_OUTCOME_TRY((details::set_population_data<Model, ModelType>(
        model, path_join(dir, "county_current_population.json"), "ID_County", county)));
    return success();
}
template <class Model, class ModelType = InfectionStateV>
IOResult<void> read_population_data_county_vaccmodel(std::vector<Model>& model, Date date, const std::vector<int>& county,
                                            const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                            const std::string& dir, int num_days)
{
    std::string id_name       = "ID_County";
    details::get_new_vaccine_growth(model, path_join(dir, "all_county_ageinf_vacc_ma7.json"), date, id_name, county,
                                    num_days);

    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY((details::set_divi_data<Model, ModelType>(model, path_join(dir, "county_divi_ma7.json"),
                                                                    id_name, county, date, scaling_factor_icu)));
    }
    else {
        log_warning("No DIVI data available for this date");
    }

    BOOST_OUTCOME_TRY((details::set_rki_data_vaccmodel<Model, ModelType>(
        model, path_join(dir, "all_county_age_ma7_rki.json"), id_name, county, date, scaling_factor_inf)));
    BOOST_OUTCOME_TRY((details::set_population_data_vaccmodel<Model>(
        model, path_join(dir, "county_current_population.json"), path_join(dir, "all_county_age_ma7_rki.json"),
        "ID_County", county, date)));
    // details::set_vaccine_data(model, path_join(dir, "vaccine_data.json"),
    // id_name, county);
    return success();
}

/**
 * @brief returns a vector with the ids of all german counties
 * @param path directory to population data
 * @return
 */
IOResult<std::vector<int>> get_county_ids(const std::string& path);

} // namespace epi

#endif
