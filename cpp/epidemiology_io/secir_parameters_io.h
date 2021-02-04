#ifndef SAVE_PARAMETERS_H
#define SAVE_PARAMETERS_H

#include <epidemiology/utils/eigen_util.h>
#include <epidemiology/utils/graph.h>
#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology_io/io.h>
#include <epidemiology_io/secir_result_io.h>

#include <tixi.h>


namespace epi
{
/**
 * @brief read contact frequency matrix and damping distributions from xml file
 * @param handle Tixi Document Handle
 * @param path Path to contact frequency matrix Tree of XML file to read from
 * @param io_mode type of xml input (see epi::write_parameter_study for more details)
 * @return returns an UncertainContactMatrix
 */
UncertainContactMatrix read_contact(TixiDocumentHandle handle, const std::string& path, int io_mode);

/**
 * @brief write contact frequency matrix and damping distributions to xml file
 * @param handle Tixi Document Handle
 * @param path Path to contact frequency matrix Tree of XML file
 * @param contact_pattern Contact frequencies, dampings, and distributions
 * @param io_mode type of xml ouput (see epi::write_parameter_study for more details)
 */
void write_contact(TixiDocumentHandle handle, const std::string& path, const UncertainContactMatrix& contact_pattern,
                   int io_mode);

/**
 * @brief reads parameter distribution and/or value from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter
 * @param io_mode type of xml input (see epi::write_parameter_study for more details)
 * @return returns a unique pointer to an UncertainValue
 */
std::unique_ptr<UncertainValue> read_element(TixiDocumentHandle handle, const std::string& path, int io_mode);

/**
 * @brief read parameter distribution from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter
 * @return returns a unique pointer to a ParameterDistribution
 */
std::unique_ptr<ParameterDistribution> read_distribution(TixiDocumentHandle handle, const std::string& path);

/**
 * @brief write parameter distribution and/or value to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter
 * @param element_name Name of parameter
 * @param element Uncertain Value of parameter
 * @param io_mode type of xml output (see epi::write_parameter_study for more details)
 * @param num_runs Number of runs of parameter study (used for predefinied samples)
 */
void write_element(const TixiDocumentHandle& handle, const std::string& path, const std::string& element_name,
                   const UncertainValue& element, int io_mode, int num_runs);

/**
 * @brief write distribution to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter
 * @param element_name Name of parameter
 * @param distribution distribution of parameter
 */
void write_distribution(const TixiDocumentHandle& handle, const std::string& path, const std::string& element_name,
                        const ParameterDistribution& distribution);

/**
 * @brief write predefined samples to xml file
 * @param handle Tixi Document Handle
 * @param path Path to Parameter of predefined samples
 * @param samples Vector of predefined samples
 */
void write_predef_sample(TixiDocumentHandle handle, const std::string& path, const std::vector<double>& samples);

/**
 * @brief read parameter space from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter space
 * @param io_mode type of xml input (see epi::write_parameter_study for more details)
 * @return returns a SecirParams object
 */
template <class AgeGroup>
SecirModel<AgeGroup> read_parameter_space(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    ReturnCode status;
    unused(status);

    int num_groups;
    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);
    assert(status == SUCCESS && ("Failed to read num_groups at " + path).c_str());

    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);

    if (num_groups != (int)AgeGroup::Count) {
        epi::log_error("Only 1, 2,3, 6 or 8 age groups allowed at the moment.");
    }

    SecirModel<AgeGroup> model;
    double read_buffer;
    status = tixiGetDoubleElement(handle, path_join(path, "StartDay").c_str(), &read_buffer);
    assert(status == SUCCESS && ("Failed to read StartDay at " + path).c_str());

    model.parameters.set_start_day(read_buffer);
    model.parameters.set_seasonality(*read_element(handle, path_join(path, "Seasonality"), io_mode));
    model.parameters.set_icu_capacity(*read_element(handle, path_join(path, "ICUCapacity"), io_mode));
    model.parameters.set_test_and_trace_capacity(*read_element(handle, path_join(path, "TestAndTraceCapacity"), io_mode));
    model.parameters.set_contact_patterns(read_contact(handle, path_join(path, "ContactFreq"), io_mode));

    for (size_t i = 0; i < static_cast<size_t>(num_groups); i++) {
        auto group_name = "Group" + std::to_string(i + 1);
        auto group_path = path_join(path, group_name);

        // populations
        auto population_path = path_join(group_path, "Population");

        status = tixiGetDoubleElement(handle, path_join(population_path, "Dead").c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read number of deaths at " + path).c_str());

        model.populations.set(read_buffer, (AgeGroup)i, InfectionType::D);

        model.populations.set(*read_element(handle, path_join(population_path, "Exposed"), io_mode), (AgeGroup)i,
                              InfectionType::E);
        model.populations.set(*read_element(handle, path_join(population_path, "Carrier"), io_mode), (AgeGroup)i,
                              InfectionType::C);
        model.populations.set(*read_element(handle, path_join(population_path, "Infectious"), io_mode), (AgeGroup)i,
                              InfectionType::I);
        model.populations.set(*read_element(handle, path_join(population_path, "Hospitalized"), io_mode), (AgeGroup)i,
                              InfectionType::H);
        model.populations.set(*read_element(handle, path_join(population_path, "ICU"), io_mode), (AgeGroup)i,
                              InfectionType::U);
        model.populations.set(*read_element(handle, path_join(population_path, "Recovered"), io_mode), (AgeGroup)i,
                              InfectionType::R);

        status = tixiGetDoubleElement(handle, path_join(population_path, "Total").c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read total population at " + path).c_str());

        model.populations.set_difference_from_group_total(read_buffer, (AgeGroup)i, (AgeGroup)i, InfectionType::S);

        // times
        auto times_path = path_join(group_path, "StageTimes");

        model.parameters.times[i].set_incubation(*read_element(handle, path_join(times_path, "Incubation"), io_mode));
        model.parameters.times[i].set_infectious_mild(
            *read_element(handle, path_join(times_path, "InfectiousMild"), io_mode));
        model.parameters.times[i].set_serialinterval(
            *read_element(handle, path_join(times_path, "SerialInterval"), io_mode));
        model.parameters.times[i].set_hospitalized_to_home(
            *read_element(handle, path_join(times_path, "HospitalizedToRecovered"), io_mode));
        model.parameters.times[i].set_home_to_hospitalized(
            *read_element(handle, path_join(times_path, "InfectiousToHospitalized"), io_mode));
        model.parameters.times[i].set_infectious_asymp(
            *read_element(handle, path_join(times_path, "InfectiousAsympt"), io_mode));
        model.parameters.times[i].set_hospitalized_to_icu(
            *read_element(handle, path_join(times_path, "HospitalizedToICU"), io_mode));
        model.parameters.times[i].set_icu_to_home(
            *read_element(handle, path_join(times_path, "ICUToRecovered"), io_mode));
        model.parameters.times[i].set_icu_to_death(*read_element(handle, path_join(times_path, "ICUToDead"), io_mode));

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");

        model.parameters.probabilities[i].set_infection_from_contact(
            *read_element(handle, path_join(probabilities_path, "InfectedFromContact"), io_mode));
        model.parameters.probabilities[i].set_carrier_infectability(
            *read_element(handle, path_join(probabilities_path, "Carrierinfectability"), io_mode));
        model.parameters.probabilities[i].set_asymp_per_infectious(
            *read_element(handle, path_join(probabilities_path, "AsympPerInfectious"), io_mode));
        model.parameters.probabilities[i].set_risk_from_symptomatic(
            *read_element(handle, path_join(probabilities_path, "RiskFromSymptomatic"), io_mode));
        model.parameters.probabilities[i].set_test_and_trace_max_risk_from_symptomatic(
            *read_element(handle, path_join(probabilities_path, "TestAndTraceMaxRiskFromSymptomatic"), io_mode));
        model.parameters.probabilities[i].set_dead_per_icu(
            *read_element(handle, path_join(probabilities_path, "DeadPerICU"), io_mode));
        model.parameters.probabilities[i].set_hospitalized_per_infectious(
            *read_element(handle, path_join(probabilities_path, "HospitalizedPerInfectious"), io_mode));
        model.parameters.probabilities[i].set_icu_per_hospitalized(
            *read_element(handle, path_join(probabilities_path, "ICUPerHospitalized"), io_mode));
    }

    return model;
}

/**
 * @brief write parameter space to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter space
 * @param params SecirParams Parameter Space with distributions of all secir parameters
 * @param num_runs Number of runs of parameter study (used for predefinied samples)
 * @param io_mode type of xml output (see epi::write_parameter_study for more details)
 */
template <class Model>
void write_parameter_space(TixiDocumentHandle handle, const std::string& path, Model const& model, int num_runs,
                           int io_mode)
{
    auto num_groups = model.parameters.get_num_groups();
    tixiAddIntegerElement(handle, path.c_str(), "NumberOfGroups", (int)num_groups, "%d");

    tixiAddDoubleElement(handle, path.c_str(), "StartDay", model.parameters.get_start_day(), "%g");
    write_element(handle, path, "Seasonality", model.parameters.get_seasonality(), io_mode, num_runs);
    write_element(handle, path, "ICUCapacity", model.parameters.get_icu_capacity(), io_mode, num_runs);
    write_element(handle, path, "TestAndTraceCapacity", model.parameters.get_test_and_trace_capacity(), io_mode, num_runs);

    for (size_t i = 0; i < num_groups; i++) {
        auto group_name = "Group" + std::to_string(i + 1);
        auto group_path = path_join(path, group_name);

        tixiCreateElement(handle, path.c_str(), group_name.c_str());

        // populations
        auto population_path = path_join(group_path, "Population");
        tixiCreateElement(handle, group_path.c_str(), "Population");

        tixiAddDoubleElement(handle, population_path.c_str(), "Total",
                             model.populations.get_group_total((typename Model::AgeGroup)i), "%g");
        tixiAddDoubleElement(handle, population_path.c_str(), "Dead",
                             model.populations.get((typename Model::AgeGroup)i, InfectionType::D), "%g");
        write_element(handle, population_path, "Exposed",
                      model.populations.get((typename Model::AgeGroup)i, InfectionType::E), io_mode, num_runs);
        write_element(handle, population_path, "Carrier",
                      model.populations.get((typename Model::AgeGroup)i, InfectionType::C), io_mode, num_runs);
        write_element(handle, population_path, "Infectious",
                      model.populations.get((typename Model::AgeGroup)i, InfectionType::I), io_mode, num_runs);
        write_element(handle, population_path, "Hospitalized",
                      model.populations.get((typename Model::AgeGroup)i, InfectionType::H), io_mode, num_runs);
        write_element(handle, population_path, "ICU",
                      model.populations.get((typename Model::AgeGroup)i, InfectionType::U), io_mode, num_runs);
        write_element(handle, population_path, "Recovered",
                      model.populations.get((typename Model::AgeGroup)i, InfectionType::R), io_mode, num_runs);

        // times
        auto times_path = path_join(group_path, "StageTimes");
        tixiCreateElement(handle, group_path.c_str(), "StageTimes");

        write_element(handle, times_path, "Incubation", model.parameters.times[i].get_incubation(), io_mode, num_runs);
        write_element(handle, times_path, "InfectiousMild", model.parameters.times[i].get_infectious_mild(), io_mode,
                      num_runs);
        write_element(handle, times_path, "SerialInterval", model.parameters.times[i].get_serialinterval(), io_mode,
                      num_runs);
        write_element(handle, times_path, "HospitalizedToRecovered",
                      model.parameters.times[i].get_hospitalized_to_home(), io_mode, num_runs);
        write_element(handle, times_path, "InfectiousToHospitalized",
                      model.parameters.times[i].get_home_to_hospitalized(), io_mode, num_runs);
        write_element(handle, times_path, "InfectiousAsympt", model.parameters.times[i].get_infectious_asymp(), io_mode,
                      num_runs);
        write_element(handle, times_path, "HospitalizedToICU", model.parameters.times[i].get_hospitalized_to_icu(),
                      io_mode, num_runs);
        write_element(handle, times_path, "ICUToRecovered", model.parameters.times[i].get_icu_to_home(), io_mode,
                      num_runs);
        write_element(handle, times_path, "ICUToDead", model.parameters.times[i].get_icu_to_dead(), io_mode, num_runs);

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");
        tixiCreateElement(handle, group_path.c_str(), "Probabilities");

        write_element(handle, probabilities_path, "InfectedFromContact",
                      model.parameters.probabilities[i].get_infection_from_contact(), io_mode, num_runs);
        write_element(handle, probabilities_path, "Carrierinfectability",
                      model.parameters.probabilities[i].get_carrier_infectability(), io_mode, num_runs);
        write_element(handle, probabilities_path, "AsympPerInfectious",
                      model.parameters.probabilities[i].get_asymp_per_infectious(), io_mode, num_runs);
        write_element(handle, probabilities_path, "RiskFromSymptomatic",
                      model.parameters.probabilities[i].get_risk_from_symptomatic(), io_mode, num_runs);
        write_element(handle, probabilities_path, "TestAndTraceMaxRiskFromSymptomatic",
                      model.parameters.probabilities[i].get_test_and_trace_max_risk_from_symptomatic(), io_mode, num_runs);
        write_element(handle, probabilities_path, "DeadPerICU", model.parameters.probabilities[i].get_dead_per_icu(),
                      io_mode, num_runs);
        write_element(handle, probabilities_path, "HospitalizedPerInfectious",
                      model.parameters.probabilities[i].get_hospitalized_per_infectious(), io_mode, num_runs);
        write_element(handle, probabilities_path, "ICUPerHospitalized",
                      model.parameters.probabilities[i].get_icu_per_hospitalized(), io_mode, num_runs);
    }

    write_contact(handle, path, model.parameters.get_contact_patterns(), io_mode);
}

/**
 * @brief read parameter study from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of parameters of study
 */
template <class AgeGroup>
ParameterStudy<SecirModel<AgeGroup>> read_parameter_study(TixiDocumentHandle handle, const std::string& path)
{
    ReturnCode status;

    int io_mode;
    int num_runs;
    double t0;
    double tmax;

    status = tixiGetIntegerElement(handle, path_join(path, "IOMode").c_str(), &io_mode);
    assert(status == SUCCESS && ("Failed to read io_mode at " + path).c_str());

    status = tixiGetIntegerElement(handle, path_join(path, "Runs").c_str(), &num_runs);
    assert(status == SUCCESS && ("Failed to read num_runs at " + path).c_str());

    status = tixiGetDoubleElement(handle, path_join(path, "T0").c_str(), &t0);
    assert(status == SUCCESS && ("Failed to read t0 at " + path).c_str());

    status = tixiGetDoubleElement(handle, path_join(path, "TMax").c_str(), &tmax);
    assert(status == SUCCESS && ("Failed to read tmax at " + path).c_str());

    unused(status);

    SecirModel<AgeGroup> model = read_parameter_space<AgeGroup>(handle, path, io_mode);
    return ParameterStudy<SecirModel<AgeGroup>>(model, t0, tmax, num_runs);
}

/**
 * @brief write parameter study to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter study
 * @param parameter_study Parameter study
 * @param io_mode type of xml output
 *        io_mode = 0: only double values of parameters are written
 *        io_mode = 1: only distributions of parameters are written
 *        io_mode = 2: both, values and distributions are written
 *        io_mode = 3: distributions are written and values are saved as predefined samples
 */
template <class Model>
void write_parameter_study(TixiDocumentHandle handle, const std::string& path,
                           const ParameterStudy<Model>& parameter_study, int io_mode = 2)
{
    tixiAddIntegerElement(handle, path.c_str(), "IOMode", io_mode, "%d");
    tixiAddIntegerElement(handle, path.c_str(), "Runs", parameter_study.get_num_runs(), "%d");
    tixiAddDoubleElement(handle, path.c_str(), "T0", parameter_study.get_t0(), "%g");
    tixiAddDoubleElement(handle, path.c_str(), "TMax", parameter_study.get_tmax(), "%g");

    write_parameter_space<Model>(handle, path, parameter_study.get_model(), parameter_study.get_num_runs(), io_mode);
}

/**
 * @brief creates xml file with a single run parameter study with std 0 (used to save parameters of individual runs)
 * @param filename Name of file
 * @param params Secir parameters used during run
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
template <class Model>
        void write_single_run_params(const int run, epi::Graph<epi::ModelNode<epi::Simulation<Model>>, epi::MigrationEdge> graph,
                                     double t0, double tmax)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");


    std::string abs_path;
    bool created = create_directory("results", abs_path);

    if (created) {
        log_info("Results are stored in {:s}/results.",
                 epi::get_current_dir_name());
    }
    else if (run == 0) {
        log_info(
            "Directory '{:s}' already exists. Results are stored in {:s}/results. Files from previous runs will be "
            "overwritten",
            epi::get_current_dir_name());
    }

    int node_id = 0;
    for (auto& node : graph.nodes()) {
        int num_runs     = 1;
        std::string path = "/Parameters";
        TixiDocumentHandle handle;
        tixiCreateDocument("Parameters", &handle);
        ParameterStudy<Model> study(node.property.get_simulation().get_model(), t0, tmax, num_runs);

        write_parameter_study(handle, path, study);

        tixiSaveDocument(handle, path_join(abs_path, ("Parameters_run" + std::to_string(run) + "_node" +
                                                          std::to_string(node_id) + ".xml"))
                                     .c_str());
        tixiCloseDocument(handle);

        save_result(node.property.get_result(), path_join(abs_path, ("Results_run" + std::to_string(run) + "_node" +
                                                                std::to_string(node_id) + ".h5")));
        node_id++;
    }
}

/**
 * @brief Creates xml file containing Parameters of one node of a graph
 * @param handle Tixi document handle
 * @param graph Graph which holds the node
 * @param node Node ID
 */
template <class Model>
void write_node(TixiDocumentHandle handle, const Graph<Model, MigrationEdge>& graph, int node)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");
    int num_runs = 1;
    int io_mode  = 2;

    std::string path = "/Parameters";

    auto model = graph.nodes()[node].property;
    int node_id =  static_cast<int>(graph.nodes()[node].id);

    tixiAddIntegerElement(handle, path.c_str(), "NodeID", node_id, "%d");
    write_parameter_space(handle, path, model, num_runs, io_mode);
}

/**
 * @brief reads parameters of a single node and saves it into the graph
 * @param node_handle Tixi document handle
 * @param graph Graph in which the node is saved
 */
template <class AgeGroup>
void read_node(TixiDocumentHandle node_handle, Graph<SecirModel<AgeGroup>, MigrationEdge>& graph)
{
    int node_id;
    tixiGetIntegerElement(node_handle, path_join("/Parameters", "NodeID").c_str(), &node_id);
    graph.add_node(node_id, read_parameter_space<AgeGroup>(node_handle, "/Parameters", 2));
}

/**
 * @brief Writes the information of a single edge into the xml file corresponding to its start node
 * @param edge_handles vecor of Tixi Document Handle. It's length is the number of nodes in the graph
 * @param path Path to document root
 * @param graph Graph which holds the edge
 * @param edge Edge ID
 */
template <class Model>
void write_edge(const std::vector<TixiDocumentHandle>& edge_handles, const std::string& path,
                const Graph<Model, MigrationEdge>& graph, int edge)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");
    int num_groups  = static_cast<int>(graph.nodes()[0].property.parameters.get_num_groups());
    int num_compart = static_cast<int>(graph.nodes()[0].property.populations.get_num_compartments()) / num_groups;

    auto start_node = static_cast<int>(graph.edges()[edge].start_node_idx);
    auto end_node   = static_cast<int>(graph.edges()[edge].end_node_idx);
    auto handle     = edge_handles[start_node];

    std::string edge_path = path_join(path, "EdgeTo" + std::to_string(end_node));
    tixiCreateElement(handle, path.c_str(), ("EdgeTo" + std::to_string(end_node)).c_str());
    tixiAddIntegerElement(handle, edge_path.c_str(), "StartNode", static_cast<int>(graph.edges()[edge].start_node_idx),
                          "%d");
    tixiAddIntegerElement(handle, edge_path.c_str(), "EndNode", static_cast<int>(graph.edges()[edge].end_node_idx),
                          "%d");
    for (int group = 0; group < num_groups; group++) {
        std::vector<double> weights;
        for (int compart = 0; compart < num_compart; compart++) {
            weights.push_back(graph.edges()[edge].property.get_coefficients()[compart + group * num_compart]);
        }
        tixiAddFloatVector(handle, edge_path.c_str(), ("Group" + std::to_string(group + 1)).c_str(), weights.data(),
                           num_compart, "%g");
    }
}

/**
 * @brief Reads information of a single edge and saves it into the graph
 * @param edge_handles vecor of Tixi Document Handle. It's length is the number of nodes in the graph
 * @param path Path to document root
 * @param graph Graph to which the edge is added
 * @param stort_node origin of the edge
 * @param end_node destination of the edge
 */
template <class Model>
void read_edge(const std::vector<TixiDocumentHandle>& edge_handles, const std::string& path,
               Graph<Model, MigrationEdge>& graph, int start_node, int end_node)
{
    ReturnCode status;
    unused(status);

    auto handle           = edge_handles[start_node];
    std::string edge_path = path_join(path, "EdgeTo" + std::to_string(end_node));
    int num_groups;
    int num_compart;

    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);
    assert(status == SUCCESS && ("Failed to read num_groups at " + path).c_str());

    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfCompartiments").c_str(), &num_compart);
    assert(status == SUCCESS && ("Failed to read num_compart at " + path).c_str());

    auto all_weights = Eigen::VectorXd(num_compart * num_groups);
    for (int group = 0; group < num_groups; group++) {
        double* weights = nullptr;
        status = tixiGetFloatVector(handle, path_join(edge_path, "Group" + std::to_string(group + 1)).c_str(), &weights,
                                    num_compart);
        if (status == SUCCESS) {
            for (int compart = 0; compart < num_compart; compart++) {
                all_weights(compart + group * num_compart) = weights[compart];
            }
            graph.add_edge(start_node, end_node, all_weights);
        }
    }
}


/**
 * @brief creates xml files for each node of a Secir simulation graph and one xml file for its edges for each node
 * @param graph Graph which should be written
 * @param dir_string directory, where graph should be stored
 */
template <class Model>
void write_graph(const Graph<Model, MigrationEdge>& graph, const std::string& dir_string)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");

    std::string abs_path;
    bool created = create_directory(dir_string, abs_path);

    if (created) {
        log_info("Results are stored in {:s}/results.",
                 epi::get_current_dir_name());
    }
    else {
        log_info("Results are stored in {:s}/results. Files from previous "
                 "graph will be "
                 "overwritten",
                 epi::get_current_dir_name());
    }
    int num_nodes   = static_cast<int>(graph.nodes().size());
    int num_edges   = static_cast<int>(graph.edges().size());
    int num_groups  = static_cast<int>(graph.nodes()[0].property.parameters.get_contact_patterns().get_cont_freq_mat().get_num_groups());
    int num_compart = static_cast<int>(graph.nodes()[0].property.populations.get_num_compartments()) / num_groups;

    std::vector<TixiDocumentHandle> edge_handles(num_nodes);
    std::string edges_path = "/Edges";
    for (auto& current_handle : edge_handles) {
        tixiCreateDocument("Edges", &current_handle);

        tixiAddIntegerElement(current_handle, edges_path.c_str(), "NumberOfNodes", num_nodes, "%d");
        tixiAddIntegerElement(current_handle, edges_path.c_str(), "NumberOfEdges", num_edges, "%d");
        tixiAddIntegerElement(current_handle, edges_path.c_str(), "NumberOfGroups", num_groups, "%d");
        tixiAddIntegerElement(current_handle, edges_path.c_str(), "NumberOfCompartiments", num_compart, "%d");
    }

    for (int edge = 0; edge < num_edges; edge++) {
        write_edge(edge_handles, edges_path, graph, edge);
    }

    for (int node = 0; node < num_nodes; node++) {
        tixiSaveDocument(edge_handles[node],
                         (path_join(abs_path.c_str(), ("GraphEdges_node" + std::to_string(node) + ".xml")).c_str()));
        tixiCloseDocument(edge_handles[node]);
    }

    for (int node = 0; node < num_nodes; node++) {
        TixiDocumentHandle node_handle;
        tixiCreateDocument("Parameters", &node_handle);
        write_node(node_handle, graph, node);
        tixiSaveDocument(node_handle, path_join(abs_path, ("GraphNode" + std::to_string(node) + ".xml")).c_str());
        tixiCloseDocument(node_handle);
    }
}

/*
 * @brief reads graph xml files and returns a Secir simulation graph
 * @param dir_string directory from where graph should be read
 */
template <class Model>
Graph<Model, MigrationEdge> read_graph(const std::string& dir_string)
{
    std::string abs_path;
    if (!directory_exists(dir_string, abs_path)) {
        log_error("Directory" + dir_string + " does not exist.");
    }

    ReturnCode status;
    unused(status);
    TixiDocumentHandle handle;
    tixiOpenDocument(path_join(abs_path, "GraphEdges_node0.xml").c_str(), &handle);

    std::string edges_path = "/Edges";

    int num_nodes;
    int num_edges;

    status = tixiGetIntegerElement(handle, path_join(edges_path, "NumberOfNodes").c_str(), &num_nodes);
    assert(status == SUCCESS && ("Failed to read num_nodes at " + edges_path).c_str());

    status = tixiGetIntegerElement(handle, path_join(edges_path, "NumberOfEdges").c_str(), &num_edges);
    assert(status == SUCCESS && ("Failed to read num_edges at " + edges_path).c_str());

    std::vector<TixiDocumentHandle> edge_handles(num_nodes);

    Graph<Model, MigrationEdge> graph;

    for (int node = 0; node < num_nodes; node++) {
        TixiDocumentHandle node_handle;
        tixiOpenDocument(path_join(abs_path,("GraphNode" + std::to_string(node) + ".xml")).c_str(), &node_handle);
        read_node(node_handle, graph);
        tixiCloseDocument(node_handle);
    }

    for (int start_node = 0; start_node < num_nodes; start_node++) {
        tixiOpenDocument(path_join(abs_path, ("GraphEdges_node" + std::to_string(start_node) + ".xml")).c_str(),
                         &edge_handles[start_node]);
        for (int end_node = 0; end_node < num_nodes; end_node++) {
            read_edge(edge_handles, edges_path, graph, start_node, end_node);
        }

        tixiCloseDocument(edge_handles[start_node]);
    }
    return graph;
}

/**
 * @brief interpolates age_ranges to param_ranges and saves ratios in interpolation
 * @param age_ranges original age ranges of the data
 * @param param_ranges age ranges to which the data should be fitted
 * @param interpolation vector of ratios that are aplied to the data of age_ranges
 * @param carry_over boolean vector which indicates whether there is an overflow from one age group to the next while interpolating data
 */
void interpolate_ages(const std::vector<double>& age_ranges, const std::vector<double>& param_ranges,
                      std::vector<std::vector<double>>& interpolation, std::vector<bool>& carry_over);

/**
 * @brief reads populations data from RKI
 * @param path Path to RKI file
 * @param id_name Name of region key column
 * @param region Key of the region of interest
 * @param month Specifies month at which the data is read
 * @param day Specifies day at which the data is read
 * @param num_inf output vector for number of infected
 * @param num_death output vector for number of dead
 * @param num_rec output vector for number of recovered
 */
void read_rki_data(std::string const& path,
                   const std::string& id_name,
                   int region, int month, int day,
                   std::vector<double>& num_inf,
                   std::vector<double>& num_death,
                   std::vector<double>& num_rec);



/**
 * @brief sets populations data from RKI into a SecirModel
 * @param model Object in which the data is set
 * @param param_ranges Age ranges of params
 * @param path Path to RKI file
 * @param id_name Name of region key column
 * @param region Key of the region of interest
 * @param month Specifies month at which the data is read
 * @param day Specifies day at which the data is read
 */
template <class Model>
void set_rki_data(Model& model, const std::vector<double>& param_ranges, const std::string& path,
                  const std::string& id_name, int region, int month, int day)
{
    std::vector<double> age_ranges     = {5., 10., 20., 25., 20., 20.};

    std::vector<std::vector<double>> interpolation(age_ranges.size()+1);
    std::vector<bool> carry_over;

    interpolate_ages(age_ranges, param_ranges, interpolation, carry_over);

    std::vector<double> num_inf;
    std::vector<double> num_death;
    std::vector<double> num_rec;

    read_rki_data(path, id_name, region, month, day, num_inf, num_death, num_rec);

    std::vector<double> interpol_inf(model.parameters.get_num_groups() + 1, 0.0);
    std::vector<double> interpol_death(model.parameters.get_num_groups() + 1, 0.0);
    std::vector<double> interpol_rec(model.parameters.get_num_groups() + 1, 0.0);

    int counter = 0;
    for (size_t i = 0; i < interpolation.size() - 1; i++) {
        for (size_t j = 0; j < interpolation[i].size(); j++) {
            interpol_inf[counter] += interpolation[i][j] * num_inf[i];
            interpol_death[counter] += interpolation[i][j] * num_death[i];
            interpol_rec[counter] += interpolation[i][j] * num_rec[i];
            if (j < interpolation[i].size() - 1 || !carry_over[i]) {
                counter++;
            }
        }
    }

    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {
        interpol_inf[i] += (double)num_inf[num_inf.size() - 1] / (double)model.parameters.get_num_groups();
        interpol_death[i] += (double)num_death[num_death.size() - 1] / (double)model.parameters.get_num_groups();
        interpol_rec[i] += (double)num_rec[num_rec.size() - 1] / (double)model.parameters.get_num_groups();
    }

    if (std::accumulate(num_inf.begin(), num_inf.end(), 0.0) > 0) {
        size_t num_groups = model.parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model.populations.set(interpol_inf[i] - interpol_death[i] - interpol_rec[i], (typename Model::AgeGroup)i, epi::InfectionType::I);
            model.populations.set(interpol_death[i], (typename Model::AgeGroup)i, epi::InfectionType::D);
            model.populations.set(interpol_rec[i], (typename Model::AgeGroup)i, epi::InfectionType::R);
        }
    }
    else {
        log_warning("No infections reported on date " + std::to_string(day) + "-" + std::to_string(month) +
                    " for region " + std::to_string(region) + ". Population data has not been set.");
    }
}

/**
 * @brief reads number of ICU patients from DIVI register into SecirParams
 * @param path Path to DIVI file
 * @param id_name Name of region key column
 * @param region Key of the region of interest
 * @param month Specifies month at which the data is read
 * @param day Specifies day at which the data is read
 * @return number of ICU patients
 */
double read_divi_data(const std::string& path,
                      const std::string& id_name,
                      int region, int month, int day);


/**
 * @brief sets populations data from DIVI register into SecirParams
 * @param params Object in which the data is set
 * @param path Path to DIVI file
 * @param id_name Name of region key column
 * @param region Key of the region of interest
 * @param month Specifies month at which the data is read
 * @param day Specifies day at which the data is read
 */
template <class Model>
void set_divi_data(Model& model, const std::string& path, const std::string& id_name, int region, int month,
                   int day)
{
    double num_icu = read_divi_data(path, id_name, region, month, day);

    if (num_icu > 0) {
        size_t num_groups = model.parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model.populations.set(num_icu / (double)num_groups, (typename Model::AgeGroup)i, epi::InfectionType::U);
        }
    }
    else {
        log_warning("No ICU patients reported on date " + std::to_string(day) + "-" + std::to_string(month) +
                    " for region " + std::to_string(region) + ".");
    }
}

/**
 * @brief reads population data from census data
 * @param path Path to RKI file
 * @param id_name Name of region key column
 * @param region Key of the region of interest
 */
std::vector<double> read_population_data(const std::string& path,
                                         const std::string& id_name,
                                         int region);

/**
 * @brief sets population data from census data
 * @param model Object in which the data is set
 * @param param_ranges Age ranges of params
 * @param path Path to RKI file
 * @param id_name Name of region key column
 * @param region Key of the region of interest
 */
template <class Model>
void set_population_data(Model& model, const std::vector<double>& param_ranges, const std::string& path,
                         const std::string& id_name, int region)
{
    std::vector<double> age_ranges     = {3., 3., 9., 3., 7., 5., 10., 10., 15., 10., 25.};

    std::vector<std::vector<double>> interpolation(age_ranges.size());
    std::vector<bool> carry_over;

    interpolate_ages(age_ranges, param_ranges, interpolation, carry_over);

    std::vector<double> num_population = read_population_data(path, id_name, region);

    std::vector<double> interpol_population(model.parameters.get_num_groups() + 1, 0.0);

    int counter = 0;
    for (size_t i = 0; i < interpolation.size() - 1; i++) {
        for (size_t j = 0; j < interpolation[i].size(); j++) {
            interpol_population[counter] += interpolation[i][j] * num_population[i];
            if (j < interpolation[i].size() - 1 || !carry_over[i]) {
                counter++;
            }
        }
    }

    if (std::accumulate(num_population.begin(), num_population.end(), 0.0) > 0) {
        size_t num_groups = model.parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model.populations.set_difference_from_group_total(interpol_population[i], (typename Model::AgeGroup)i, (typename Model::AgeGroup)i, epi::InfectionType::S);
        }
    }
    else {
        log_warning("No population data available for region " + std::to_string(region) + ". Population data has not been set.");
    }

}


/**
 * @brief reads population data from population files for the whole country
 * @param params Parameters in which the data is set
 * @param param_ranges Vector which specifies the age ranges of params. Needs to add up to 100
 * @param month specifies month at which the data is read
 * @param day specifies day at which the data is read
 * @param dir directory of files
 */
template <class Model>
void read_population_data_germany(Model& model, const std::vector<double>& param_ranges, int month, int day,
                                  const std::string& dir)
{
    assert(param_ranges.size() == model.parameters.get_num_groups() &&
           "size of param_ranges needs to be the same size as the number of groups in model.parameters");
    assert(std::accumulate(param_ranges.begin(), param_ranges.end(), 0.0) == 100. && "param_ranges must add up to 100");

    std::string id_name;

    set_rki_data(model, param_ranges, path_join(dir, "all_age_rki.json"), id_name, 0, month, day);
    set_divi_data(model, path_join(dir, "germany_divi.json"), id_name, 0, month, day);
    set_population_data(model, param_ranges, path_join(dir, "county_current_population.json"), "ID_County", 0);
}

/**
 * @brief reads population data from population files for the specefied state
 * @param params Parameters in which the data is set
 * @param param_ranges Vector which specifies the age ranges of params. Needs to add up to 100
 * @param month specifies month at which the data is read
 * @param day specifies day at which the data is read
 * @param state region key of state of interest
 * @param dir directory of files
 */
template <class Model>
void read_population_data_state(Model& model, const std::vector<double>& param_ranges, int month, int day,
                                int state, const std::string& dir)
{
    assert(state > 0 && state <= 16 && "State must be between 1 and 16");
    assert(param_ranges.size() == model.parameters.get_num_groups() &&
           "size of param_ranges needs to be the same size as the number of groups in params");
    assert(std::accumulate(param_ranges.begin(), param_ranges.end(), 0.0) == 100. && "param_ranges must add up to 100");

    std::string id_name = "ID_State";

    set_rki_data(model, param_ranges, path_join(dir, "all_state_age_rki.json"), id_name, state, month, day);
    set_divi_data(model, path_join(dir, "state_divi.json"), id_name, state, month, day);
    set_population_data(model, param_ranges, path_join(dir, "county_current_population.json"), "ID_County", state);
}

/**
 * @brief reads population data from population files for the specefied county
 * @param params Parameters in which the data is set
 * @param param_ranges Vector which specifies the age ranges of params. Needs to add up to 100
 * @param month specifies month at which the data is read
 * @param day specifies day at which the data is read
 * @param county region key of county of interest
 * @param dir directory of files
 */
template <class Model>
void read_population_data_county(Model& model, const std::vector<double>& param_ranges, int month, int day,
                                 int county, const std::string& dir)
{
    assert(county > 999 && "State must be between 1 and 16");
    assert(param_ranges.size() == model.parameters.get_num_groups() &&
           "size of param_ranges needs to be the same size as the number of groups in params");
    assert(std::accumulate(param_ranges.begin(), param_ranges.end(), 0.0) == 100. && "param_ranges must add up to 100");

    std::string id_name = "ID_County";

    set_rki_data(model, param_ranges, path_join(dir, "all_county_age_rki.json"), id_name, county, month, day);
    set_divi_data(model, path_join(dir, "county_divi.json"), id_name, county, month, day);
    set_population_data(model, param_ranges, path_join(dir, "county_current_population.json"), "ID_County", county);
}

} // namespace epi

#endif
