#ifndef SAVE_PARAMETERS_H
#define SAVE_PARAMETERS_H

#include <epidemiology/utils/eigen_util.h>
#include <epidemiology/utils/graph.h>
#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/parameter_studies.h>
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
auto read_parameter_space(TixiDocumentHandle handle, const std::string& path, int io_mode);

/**
 * @brief write parameter space to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter space
 * @param params SecirParams Parameter Space with distributions of all secir parameters
 * @param num_runs Number of runs of parameter study (used for predefinied samples)
 * @param io_mode type of xml output (see epi::write_parameter_study for more details)
 */
template <class AgeGroup>
void write_parameter_space(
    TixiDocumentHandle handle, const std::string& path,
    CompartmentalModel<Populations<AgeGroup, InfectionType>, SecirParams<(size_t)AgeGroup::Count>>& model, int num_runs,
    int io_mode)
{
    auto num_groups = model.parameters.get_num_groups();
    tixiAddIntegerElement(handle, path.c_str(), "NumberOfGroups", num_groups, "%d");

    tixiAddDoubleElement(handle, path.c_str(), "StartDay", parameters.get_start_day(), "%g");
    write_element(handle, path, "Seasonality", parameters.get_seasonality(), io_mode, num_runs);
    write_element(handle, path, "ICUCapacity", parameters.get_icu_capacity(), io_mode, num_runs);

    for (size_t i = 0; i < num_groups; i++) {
        auto group_name = "Group" + std::to_string(i + 1);
        auto group_path = path_join(path, group_name);

        tixiCreateElement(handle, path.c_str(), group_name.c_str());

        // populations
        auto population_path = path_join(group_path, "Population");
        tixiCreateElement(handle, group_path.c_str(), "Population");

        tixiAddDoubleElement(handle, population_path.c_str(), "Total", model.populations.get_group_total((AgeGroup)i),
                             "%g");
        tixiAddDoubleElement(handle, population_path.c_str(), "Dead", model.populations.get(i, InfectionType::D), "%g");
        write_element(handle, population_path, "Exposed", model.populations.get((AgeGroup)i, InfectionType::E), io_mode,
                      num_runs);
        write_element(handle, population_path, "Carrier", model.populations.get((AgeGroup)i, InfectionType::C), io_mode,
                      num_runs);
        write_element(handle, population_path, "Infectious", model.populations.get((AgeGroup)i, InfectionType::I),
                      io_mode, num_runs);
        write_element(handle, population_path, "Hospitalized", model.populations.get((AgeGroup)i, InfectionType::H),
                      io_mode, num_runs);
        write_element(handle, population_path, "ICU", model.populations.get((AgeGroup)i, InfectionType::U), io_mode,
                      num_runs);
        write_element(handle, population_path, "Recovered", model.populations.get((AgeGroup)i, InfectionType::R),
                      io_mode, num_runs);

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
        write_element(handle, probabilities_path, "DeadPerICU", model.parameters.probabilities[i].get_dead_per_icu(),
                      io_mode, num_runs);
        write_element(handle, probabilities_path, "HospitalizedPerInfectious",
                      model.parameters.probabilities[i].get_hospitalized_per_infectious(), io_mode, num_runs);
        write_element(handle, probabilities_path, "ICUPerHospitalized",
                      model.parameters.probabilities[i].get_icu_per_hospitalized(), io_mode, num_runs);
    }

    write_contact(handle, path, model.parameters.get_contact_patterns(), io_mode, num_runs);
}

/**
 * @brief read parameter study from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of parameters of study
 */
auto read_parameter_study(TixiDocumentHandle handle, const std::string& path);

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

    write_parameter_space(handle, path, parameter_study.get_secir_params(), parameter_study.get_num_runs(), io_mode);
}

/**
 * @brief creates xml file with a single run parameter study with std 0 (used to save parameters of individual runs)
 * @param filename Name of file
 * @param params Secir parameters used during run
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
template <class Model>
void write_single_run_params(const int run, const Model& model, double t0, double tmax,
                             const TimeSeries<double>& result, int node)
{

    int num_runs     = 1;
    std::string path = "/Parameters";
    TixiDocumentHandle handle;
    tixiCreateDocument("Parameters", &handle);
    ParameterStudy<Model> study(model, t0, tmax, num_runs);

    boost::filesystem::path dir("results");

    bool created = boost::filesystem::create_directory(dir);

    if (created) {
        log_info("Directory '{:s}' was created. Results are stored in {:s}/results.", dir.string(),
                 epi::get_current_dir_name());
    }
    else {
        log_info(
            "Directory '{:s}' already exists. Results are stored in {:s}/ results. Files from previous runs will be "
            "overwritten",
            dir.string(), epi::get_current_dir_name());
    }

    write_parameter_study(handle, path, study);
    tixiSaveDocument(
        handle,
        (dir / ("Parameters_run" + std::to_string(run) + "_node" + std::to_string(node) + ".xml")).string().c_str());
    tixiCloseDocument(handle);

    save_result(result,
                (dir / ("Results_run" + std::to_string(run) + "_node" + std::to_string(node) + ".h5")).string());
}

/**
 * @brief Creates xml file containing Parameters of one node of a graph
 * @param graph Graph which holds the node
 * @param node Node ID
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
template <class Model>
void write_node(const Graph<Model, MigrationEdge>& graph, int node, double t0, double tmax)
{
    int num_runs = 1;
    int io_mode  = 2;

    std::string path = "/Parameters";
    TixiDocumentHandle handle;
    tixiCreateDocument("Parameters", &handle);

    tixiAddIntegerElement(handle, path.c_str(), "NodeID", node, "%d");

    auto model = graph.nodes()[node];

    write_parameter_space(handle, path, model, num_runs, io_mode);
    tixiSaveDocument(handle, ("GraphNode" + std::to_string(node) + ".xml").c_str());
    tixiCloseDocument(handle);
}

/**
 * @brief reads parameters of a single node and saves it into the graph
 * @param graph Graph in which the node is saved
 * @param node Node ID
 */
template <class Model>
void read_node(Graph<Model, MigrationEdge>& graph, int node)
{
    TixiDocumentHandle node_handle;
    tixiOpenDocument(("GraphNode" + std::to_string(node) + ".xml").c_str(), &node_handle);

    graph.add_node(read_parameter_space(node_handle, "/Parameters", 2));

    tixiCloseDocument(node_handle);
}

/**
 * @brief Writes the information of a single edge into a xml file
 * @param handle Tixi Document Handle
 * @param path Path to document root
 * @param graph Graph which holds the edge
 * @param edge Edge ID
 */
template <class Model>
void write_edge(TixiDocumentHandle handle, const std::string& path, const Graph<Model, MigrationEdge>& graph, int edge)
{

    int num_groups  = static_cast<int>(graph.nodes()[0].get_num_groups());
    int num_compart = static_cast<int>(graph.nodes()[0].populations.get_num_compartments()) / num_groups;

    std::string edge_path = path_join(path, "Edge" + std::to_string(edge));
    tixiCreateElement(handle, path.c_str(), ("Edge" + std::to_string(edge)).c_str());
    tixiAddIntegerElement(handle, edge_path.c_str(), "StartNode", static_cast<int>(graph.edges()[edge].start_node_idx),
                          "%d");
    tixiAddIntegerElement(handle, edge_path.c_str(), "EndNode", static_cast<int>(graph.edges()[edge].end_node_idx),
                          "%d");
    for (int group = 0; group < num_groups; group++) {
        std::vector<double> weights;
        for (int compart = 0; compart < num_compart; compart++) {
            weights.push_back(graph.edges()[edge].property.coefficients[compart + group * num_compart]);
        }
        tixiAddFloatVector(handle, edge_path.c_str(), ("Group" + std::to_string(group + 1)).c_str(), weights.data(),
                           num_compart, "%g");
    }
}

/**
 * @brief Reads information of a single edge and saves it into the graph
 * @param handle Tixi document handle
 * @param path Path to document root
 * @param graph Graph to which the edge is added
 * @param edge Edge ID
 */
template <class Model>
void read_edge(TixiDocumentHandle handle, const std::string& path, Graph<Model, MigrationEdge>& graph, int edge)
{

    std::string edge_path = path_join(path, "Edge" + std::to_string(edge));
    int num_groups;
    int num_compart;
    int start_node;
    int end_node;

    tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);
    tixiGetIntegerElement(handle, path_join(path, "NumberOfCompartiments").c_str(), &num_compart);
    tixiGetIntegerElement(handle, path_join(edge_path, "StartNode").c_str(), &start_node);
    tixiGetIntegerElement(handle, path_join(edge_path, "EndNode").c_str(), &end_node);

    auto all_weights = Eigen::VectorXd(num_compart * num_groups);
    for (int group = 0; group < num_groups; group++) {
        double* weights = nullptr;
        tixiGetFloatVector(handle, path_join(edge_path, "Group" + std::to_string(group + 1)).c_str(), &weights,
                           num_compart);
        for (int compart = 0; compart < num_compart; compart++) {
            all_weights(compart + group * num_compart) = weights[compart];
        }
    }
    graph.add_edge(start_node, end_node, all_weights);
}

/**
 * @brief creates xml files for each node of a Secir simulation graph and one xml file for its edges
 * @param graph Graph which should be written
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
//TO-DO: Implement apropriate File System for XML FIles
template <class Model>
void write_graph(const Graph<Model, MigrationEdge>& graph, double t0, double tmax)
{
    std::string edges_path = "/Edges";
    TixiDocumentHandle handle;
    tixiCreateDocument("Edges", &handle);

    int num_nodes   = static_cast<int>(graph.nodes().size());
    int num_edges   = static_cast<int>(graph.edges().size());
    int num_groups  = graph.nodes()[0].parameters.get_contact_patterns().get_cont_freq_mat().get_size();
    int num_compart = static_cast<int>(graph.nodes()[0].populations.get_num_compartments()) / num_groups;

    tixiAddIntegerElement(handle, edges_path.c_str(), "NumberOfNodes", num_nodes, "%d");
    tixiAddIntegerElement(handle, edges_path.c_str(), "NumberOfEdges", num_edges, "%d");
    tixiAddIntegerElement(handle, edges_path.c_str(), "NumberOfGroups", num_groups, "%d");
    tixiAddIntegerElement(handle, edges_path.c_str(), "NumberOfCompartiments", num_compart, "%d");

    for (int edge = 0; edge < num_edges; edge++) {
        write_edge(handle, edges_path, graph, edge);
    }

    tixiSaveDocument(handle, "GraphEdges.xml");
    tixiCloseDocument(handle);

    for (int node = 0; node < num_nodes; node++) {
        write_node(graph, node, t0, tmax);
    }
}

/**
 * @brief reads graph xml files and returns a Secir simulation graph
 */
template <class Model>
Graph<Model, MigrationEdge> read_graph()
{
    TixiDocumentHandle handle;
    tixiOpenDocument("GraphEdges.xml", &handle);

    std::string edges_path = "/Edges";

    int num_nodes;
    int num_edges;

    tixiGetIntegerElement(handle, path_join(edges_path, "NumberOfNodes").c_str(), &num_nodes);
    tixiGetIntegerElement(handle, path_join(edges_path, "NumberOfEdges").c_str(), &num_edges);

    Graph<Model, MigrationEdge> graph;

    for (int node = 0; node < num_nodes; node++) {
        read_node(graph, node);
    }

    for (int edge = 0; edge < num_edges; edge++) {
        read_edge(handle, edges_path, graph, edge);
    }
    tixiCloseDocument(handle);
    return model;
}

} // namespace epi

#endif
