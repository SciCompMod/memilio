#ifndef SAVE_PARAMETERS_H
#define SAVE_PARAMETERS_H

#include <epidemiology/utils/eigen_util.h>
#include <epidemiology/utils/graph.h>
#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology_io/io.h>
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology/utils/date.h>
#include <epidemiology/secir/analyze_result.h>

#include <tixi.h>

namespace epi
{
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
 * @param region vector of keys of the region of interest
 * @param year Specifies year at which the data is read
 * @param month Specifies month at which the data is read
 * @param day Specifies day at which the data is read
 * @param num_* output vector for number of people in the corresponding compartement
 * @param t_* vector average time it takes to get from one compartement to another for each age group
 * @param mu_* vector probabilities to get from one compartement to another for each age group
 */
    void read_rki_data(std::string const& path, const std::string& id_name, std::vector<int> const& region, Date date,
                       std::vector<std::vector<double>>& num_exp, std::vector<std::vector<double>>& num_car,
                       std::vector<std::vector<double>>& num_inf, std::vector<std::vector<double>>& num_hosp,
                       std::vector<std::vector<double>>& num_icu, std::vector<std::vector<double>>& num_death,
                       std::vector<std::vector<double>>& num_rec, const std::vector<std::vector<int>>& t_car_to_rec,
                       const std::vector<std::vector<int>>& t_car_to_inf,
                       const std::vector<std::vector<int>>& t_exp_to_car,
                       const std::vector<std::vector<int>>& t_inf_to_rec,
                       const std::vector<std::vector<int>>& t_inf_to_hosp,
                       const std::vector<std::vector<int>>& t_hosp_to_rec,
                       const std::vector<std::vector<int>>& t_hosp_to_icu,
                       const std::vector<std::vector<int>>& t_icu_to_dead,
                       const std::vector<std::vector<double>>& mu_C_R, const std::vector<std::vector<double>>& mu_I_H,
                       const std::vector<std::vector<double>>& mu_H_U, const std::vector<double>& scaling_factor_inf);

    /**
 * @brief sets populations data from RKI into a SecirModel
 * @param model vector of objects in which the data is set
 * @param path Path to RKI file
 * @param id_name Name of region key column
 * @param region vector of keys of the region of interest
 * @param year Specifies year at which the data is read
 * @param month Specifies month at which the data is read
 * @param day Specifies day at which the data is read
 * @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
 */
    void set_rki_data(std::vector<SecirModel>& model, const std::string& path, const std::string& id_name,
                      std::vector<int> const& region, Date date, const std::vector<double>& scaling_factor_inf);

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
    void read_divi_data(const std::string& path, const std::string& id_name, const std::vector<int>& vregion, Date date,
                        std::vector<double>& vnum_icu);

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
    void set_divi_data(std::vector<SecirModel>& model, const std::string& path, const std::string& id_name,
                       const std::vector<int> vregion, Date date, double scaling_factor_icu);

    /**
 * @brief reads population data from census data
 * @param path Path to RKI file
 * @param id_name Name of region key column
 * @param vregion vector of keys of the regions of interest
 */
    std::vector<std::vector<double>> read_population_data(const std::string& path, const std::string& id_name,
                                                          const std::vector<int>& vregion);

    /**
 * @brief sets population data from census data
 * @param model vector of objects in which the data is set
 * @param path Path to RKI file
 * @param id_name Name of region key column
 * @param vregion vector of keys of the regions of interest
 */
    void set_population_data(std::vector<SecirModel>& model, const std::string& path, const std::string& id_name,
                             const std::vector<int> vregion);
} //namespace details

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
 * write an Eigen matrix to xml file.
 * @param handle tixi document handle.
 * @param path path in the document.
 * @param name name of the element in the document.
 * @param m matrix.
 */
template <class M>
void write_matrix(TixiDocumentHandle handle, const std::string& path, const std::string& name, M&& m)
{
    tixiCreateElement(handle, path.c_str(), name.c_str());
    auto matrix_path = path_join(path, name);
    tixiAddIntegerAttribute(handle, matrix_path.c_str(), "Rows", (int)m.rows(), "%d");
    tixiAddIntegerAttribute(handle, matrix_path.c_str(), "Cols", (int)m.cols(), "%d");
    //Matrix may be column major but we want to output row major
    std::vector<double> coeffs(m.size());
    for (Eigen::Index i = 0; i < m.rows(); ++i) {
        for (Eigen::Index j = 0; j < m.cols(); ++j) {
            coeffs[i * m.cols() + j] = m(i, j);
        }
    }
    tixiAddFloatVector(handle, matrix_path.c_str(), "Coefficients", coeffs.data(), (int)coeffs.size(), "%.18g");
}

/**
 * read an Eigen matrix from an xml file.
 * @param handle tixi document handle.
 * @param path path to the matrix element in the document.
 * @return a matrix.
 */
template <class M = Eigen::MatrixXd>
M read_matrix(TixiDocumentHandle handle, const std::string& path)
{
    auto status = SUCCESS;
    unused(status);
    int rows, cols;
    status = tixiGetIntegerAttribute(handle, path.c_str(), "Rows", &rows);
    assert(status == SUCCESS && "Failed to read matrix rows.");
    status = tixiGetIntegerAttribute(handle, path.c_str(), "Cols", &cols);
    assert(status == SUCCESS && "Failed to read matrix columns.");
    double* coeffs;
    M m{rows, cols};
    status = tixiGetFloatVector(handle, path_join(path, "Coefficients").c_str(), &coeffs, (int)m.size());
    assert(status == SUCCESS && "Failed to read matrix coefficients.");
    //values written as row major, but matrix type might be column major
    //so we can't just copy all coeffs to m.data()
    for (Eigen::Index i = 0; i < m.rows(); ++i) {
        for (Eigen::Index j = 0; j < m.cols(); ++j) {
            m(i, j) = coeffs[i * m.cols() + j];
        }
    }
    return m;
}

/**
 * write an instance of DampingMatrixExpressionGroup to an xml document.
 * @param handle tixi document handle.
 * @param path path in the document.
 * @param cfmc the object to write.
 */
template <class C>
void write_damping_matrix_expression_collection(TixiDocumentHandle handle, const std::string& path, const C& cfmc)
{
    tixiCreateElement(handle, path.c_str(), "ContactMatrixGroup");
    auto collection_path = path_join(path, "ContactMatrixGroup");
    for (size_t i = 0; i < cfmc.get_num_matrices(); ++i) {
        auto& cfm     = cfmc[i];
        auto cfm_name = "ContactMatrix" + std::to_string(i + 1);
        tixiCreateElement(handle, collection_path.c_str(), cfm_name.c_str());
        auto cfm_path = path_join(collection_path, cfm_name);
        write_matrix(handle, cfm_path, "Baseline", cfm.get_baseline());
        write_matrix(handle, cfm_path, "Minimum", cfm.get_minimum());
        tixiCreateElement(handle, cfm_path.c_str(), "Dampings");
        auto dampings_path = path_join(cfm_path, "Dampings");
        for (size_t j = 0; j < cfm.get_dampings().size(); ++j) {
            auto& damping     = cfm.get_dampings()[j];
            auto damping_name = "Damping" + std::to_string(j + 1);
            tixiCreateElement(handle, dampings_path.c_str(), damping_name.c_str());
            auto damping_path = path_join(dampings_path, damping_name);
            tixiAddDoubleAttribute(handle, damping_path.c_str(), "Time", double(damping.get_time()), "%.18g");
            tixiAddIntegerAttribute(handle, damping_path.c_str(), "Type", int(damping.get_type()), "%d");
            tixiAddIntegerAttribute(handle, damping_path.c_str(), "Level", int(damping.get_level()), "%d");
            write_matrix(handle, damping_path, "Values", damping.get_coeffs());
        }
    }
}

/**
 * read an instance of DampingMatrixExpressionGroup from an xml document.
 * @param handle tixi document handle.
 * @param path path in the document.
 */
template <class C = ContactMatrixGroup>
C read_damping_matrix_expression_collection(TixiDocumentHandle handle, const std::string& path)
{
    auto status = SUCCESS;
    unused(status);

    auto collection_path = path_join(path, "ContactMatrixGroup");
    int num_matrices;
    status = tixiGetNumberOfChilds(handle, collection_path.c_str(), &num_matrices);
    assert(status == SUCCESS && "Failed to read ContactMatrixGroup.");
    C cfmc{size_t(num_matrices), 1};
    for (size_t i = 0; i < cfmc.get_num_matrices(); ++i) {
        auto cfm_path = path_join(collection_path, "ContactMatrix" + std::to_string(i + 1));
        cfmc[i]       = typename C::value_type(read_matrix<typename C::Matrix>(handle, path_join(cfm_path, "Baseline")),
                                         read_matrix<typename C::Matrix>(handle, path_join(cfm_path, "Minimum")));
        auto dampings_path = path_join(cfm_path, "Dampings");
        int num_dampings;
        status = tixiGetNumberOfChilds(handle, dampings_path.c_str(), &num_dampings);
        assert(status == SUCCESS && "Failed to read Dampings from ContactMatrix.");
        for (int j = 0; j < num_dampings; ++j) {
            auto damping_path = path_join(dampings_path, "Damping" + std::to_string(j + 1));
            double t;
            status = tixiGetDoubleAttribute(handle, damping_path.c_str(), "Time", &t);
            assert(status == SUCCESS && "Failed to read Damping Time.");
            int type;
            status = tixiGetIntegerAttribute(handle, damping_path.c_str(), "Type", &type);
            assert(status == SUCCESS && "Failed to read Damping Type.");
            int level;
            status = tixiGetIntegerAttribute(handle, damping_path.c_str(), "Level", &level);
            assert(status == SUCCESS && "Failed to read Damping Level.");
            cfmc[i].add_damping(read_matrix<typename C::Matrix>(handle, path_join(damping_path, "Values")),
                                DampingLevel(level), DampingType(type), SimulationTime(t));
        }
    }
    return cfmc;
}

/**
 * @brief read parameter space from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter space
 * @param io_mode type of xml input (see epi::write_parameter_study for more details)
 * @return returns a SecirParams object
 */
SecirModel read_parameter_space(TixiDocumentHandle handle, const std::string& path, int io_mode);

/**
 * @brief write parameter space to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter space
 * @param params SecirParams Parameter Space with distributions of all secir parameters
 * @param num_runs Number of runs of parameter study (used for predefinied samples)
 * @param io_mode type of xml output (see epi::write_parameter_study for more details)
 */
void write_parameter_space(TixiDocumentHandle handle, const std::string& path, SecirModel const& model, int num_runs,
                           int io_mode);

/**
 * @brief read parameter study from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of parameters of study
 */
ParameterStudy<SecirModel> read_parameter_study(TixiDocumentHandle handle, const std::string& path);

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
void write_parameter_study(TixiDocumentHandle handle, const std::string& path,
                           const ParameterStudy<SecirModel>& parameter_study, int io_mode = 2);

/**
 * @brief creates xml file with a single run parameter study with std 0 (used to save parameters of individual runs)
 * @param filename Name of file
 * @param params Secir parameters used during run
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
void write_single_run_params(const int run,
                             epi::Graph<epi::ModelNode<epi::Simulation<SecirModel>>, epi::MigrationEdge> graph,
                             double t0, double tmax);

/**
 * @brief Creates xml file containing Parameters of one node of a graph
 * @param handle Tixi document handle
 * @param graph Graph which holds the node
 * @param node Node ID
 */
template <class Model>
void write_node(TixiDocumentHandle handle, const Graph<Model, MigrationParameters>& graph, int node)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");
    int num_runs = 1;
    int io_mode  = 2;

    std::string path = "/Parameters";

    auto model  = graph.nodes()[node].property;
    int node_id = static_cast<int>(graph.nodes()[node].id);

    tixiAddIntegerElement(handle, path.c_str(), "NodeID", node_id, "%d");
    write_parameter_space(handle, path, model, num_runs, io_mode);
}

/**
 * @brief reads parameters of a single node and saves it into the graph
 * @param node_handle Tixi document handle
 * @param graph Graph in which the node is saved
 */
template <class Model>
void read_node(TixiDocumentHandle node_handle, Graph<SecirModel, MigrationParameters>& graph)
{
    int node_id;
    tixiGetIntegerElement(node_handle, path_join("/Parameters", "NodeID").c_str(), &node_id);
    graph.add_node(node_id, read_parameter_space(node_handle, "/Parameters", 2));
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
                const Graph<Model, MigrationParameters>& graph, int edge)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");

    auto start_node = static_cast<int>(graph.edges()[edge].start_node_idx);
    auto end_node   = static_cast<int>(graph.edges()[edge].end_node_idx);
    auto handle     = edge_handles[start_node];

    std::string edge_path = path_join(path, "EdgeTo" + std::to_string(end_node));
    tixiCreateElement(handle, path.c_str(), ("EdgeTo" + std::to_string(end_node)).c_str());
    tixiAddIntegerElement(handle, edge_path.c_str(), "StartNode", static_cast<int>(graph.edges()[edge].start_node_idx),
                          "%d");
    tixiAddIntegerElement(handle, edge_path.c_str(), "EndNode", static_cast<int>(graph.edges()[edge].end_node_idx),
                          "%d");
    tixiCreateElement(handle, edge_path.c_str(), "Parameters");
    write_damping_matrix_expression_collection(handle, path_join(edge_path, "Parameters").c_str(),
                                               graph.edges()[edge].property.get_coefficients());
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
               Graph<Model, MigrationParameters>& graph, int start_node, int end_node)
{
    ReturnCode status;
    unused(status);

    auto handle           = edge_handles[start_node];
    std::string edge_path = path_join(path, "EdgeTo" + std::to_string(end_node));
    int num_groups;
    int num_compart;

    if (tixiCheckElement(handle, edge_path.c_str()) != SUCCESS) {
        return;
    }

    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);
    assert(status == SUCCESS && ("Failed to read num_groups at " + path).c_str());

    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfCompartiments").c_str(), &num_compart);
    assert(status == SUCCESS && ("Failed to read num_compart at " + path).c_str());

    auto coefficients = read_damping_matrix_expression_collection<MigrationCoefficientGroup>(
        handle, path_join(edge_path, "Parameters").c_str());
    graph.add_edge(start_node, end_node, MigrationParameters(coefficients));
}

/**
 * @brief creates xml files for each node of a Secir simulation graph and one xml file for its edges for each node
 * @param graph Graph which should be written
 * @param dir_string directory, where graph should be stored
 */
template <class Model>
void write_graph(const Graph<Model, MigrationParameters>& graph, const std::string& dir_string)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");

    std::string abs_path;
    bool created = create_directory(dir_string, abs_path);

    if (created) {
        log_info("Results are stored in {:s}/results.", epi::get_current_dir_name());
    }
    else {
        log_info("Results are stored in {:s}/results. Files from previous "
                 "graph will be "
                 "overwritten",
                 epi::get_current_dir_name());
    }
    int num_nodes  = static_cast<int>(graph.nodes().size());
    int num_edges  = static_cast<int>(graph.edges().size());
    int num_groups = static_cast<int>(
        graph.nodes()[0].property.parameters.get_contact_patterns().get_cont_freq_mat().get_num_groups());
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
Graph<Model, MigrationParameters> read_graph(const std::string& dir_string)
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

    Graph<Model, MigrationParameters> graph;

    for (int node = 0; node < num_nodes; node++) {
        TixiDocumentHandle node_handle;
        tixiOpenDocument(path_join(abs_path, ("GraphNode" + std::to_string(node) + ".xml")).c_str(), &node_handle);
        read_node<Model>(node_handle, graph);
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
* @brief sets populations data from RKI into a SecirModel
* @param model vector of objects in which the data is set
* @param path Path to RKI file
* @param id_name Name of region key column
* @param region vector of keys of the region of interest
* @param year Specifies year at which the data is read
* @param month Specifies month at which the data is read
* @param day Specifies day at which the data is read
* @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
*/
template <class Model>
void extrapolate_rki_results(std::vector<Model>& model, const std::string& dir, const std::string& results_dir,
                             std::vector<int> const& region, Date date, const std::vector<double>& scaling_factor_inf,
                             double scaling_factor_icu, int num_days)
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

    std::vector<std::vector<double>> mu_C_R{model.size()};
    std::vector<std::vector<double>> mu_I_H{model.size()};
    std::vector<std::vector<double>> mu_H_U{model.size()};

    std::vector<double> sum_mu_I_U(region.size(), 0);
    std::vector<std::vector<double>> mu_I_U{model.size()};

    for (size_t county = 0; county < model.size(); county++) {
        for (size_t group = 0; group < age_ranges.size(); group++) {

            t_car_to_inf[county].push_back(
                static_cast<int>(2 * (model[county].parameters.times[group].get_incubation() -
                                      model[county].parameters.times[group].get_serialinterval())));
            t_car_to_rec[county].push_back(static_cast<int>(
                t_car_to_inf[county][group] + 0.5 * model[county].parameters.times[group].get_infectious_mild()));
            t_exp_to_car[county].push_back(
                static_cast<int>(2 * model[county].parameters.times[group].get_serialinterval() -
                                 model[county].parameters.times[group].get_incubation()));
            t_inf_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.times[group].get_infectious_mild()));
            t_inf_to_hosp[county].push_back(
                static_cast<int>(model[county].parameters.times[group].get_home_to_hospitalized()));
            t_hosp_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.times[group].get_hospitalized_to_home()));
            t_hosp_to_icu[county].push_back(
                static_cast<int>(model[county].parameters.times[group].get_hospitalized_to_icu()));
            t_icu_to_dead[county].push_back(static_cast<int>(model[county].parameters.times[group].get_icu_to_dead()));

            mu_C_R[county].push_back(model[county].parameters.probabilities[group].get_asymp_per_infectious());
            mu_I_H[county].push_back(model[county].parameters.probabilities[group].get_hospitalized_per_infectious());
            mu_H_U[county].push_back(model[county].parameters.probabilities[group].get_icu_per_hospitalized());

            sum_mu_I_U[county] += model[county].parameters.probabilities[group].get_icu_per_hospitalized() *
                                  model[county].parameters.probabilities[group].get_hospitalized_per_infectious();
            mu_I_U[county].push_back(model[county].parameters.probabilities[group].get_icu_per_hospitalized() *
                                     model[county].parameters.probabilities[group].get_hospitalized_per_infectious());
        }
    }

    std::vector<TimeSeries<double>> rki_data(
        region.size(), TimeSeries<double>::zero(num_days, (size_t)InfectionState::Count * age_ranges.size()));

    for (size_t j = 0; j < static_cast<size_t>(num_days); j++) {
        std::vector<std::vector<double>> num_inf(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_death(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_exp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_car(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_hosp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> dummy_icu(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<double> num_icu(model.size(), 0.0);

        details::read_rki_data(path_join(dir, "all_county_age_rki_ma.json"), id_name, region, date, num_exp, num_car,
                               num_inf, num_hosp, dummy_icu, num_death, num_rec, t_car_to_rec, t_car_to_inf,
                               t_exp_to_car, t_inf_to_rec, t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead,
                               mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf);
        details::read_divi_data(path_join(dir, "county_divi.json"), id_name, region, date, num_icu);
        std::vector<std::vector<double>> num_population =
            details::read_population_data(path_join(dir, "county_current_population.json"), id_name, region);

        for (size_t i = 0; i < region.size(); i++) {
            for (size_t age = 0; age < age_ranges.size(); age++) {
                rki_data[i][j]((size_t)InfectionState::Exposed + (size_t)epi::InfectionState::Count * age) =
                    num_exp[i][age];
                rki_data[i][j]((size_t)InfectionState::Carrier + (size_t)InfectionState::Count * age) = num_car[i][age];
                rki_data[i][j]((size_t)InfectionState::Infected + (size_t)InfectionState::Count * age) =
                    num_inf[i][age];
                rki_data[i][j]((size_t)InfectionState::Hospitalized + (size_t)InfectionState::Count * age) =
                    num_hosp[i][age];
                rki_data[i][j]((size_t)InfectionState::ICU + (size_t)InfectionState::Count * age) =
                    scaling_factor_icu * num_icu[i] * mu_I_U[i][age] / sum_mu_I_U[i];
                rki_data[i][j]((size_t)InfectionState::Recovered + (size_t)InfectionState::Count * age) =
                    num_rec[i][age];
                rki_data[i][j]((size_t)InfectionState::Dead + (size_t)InfectionState::Count * age) = num_death[i][age];
                rki_data[i][j]((size_t)InfectionState::Susceptible + (size_t)InfectionState::Count * age) =
                    num_population[i][age] - num_exp[i][age] - num_car[i][age] - num_inf[i][age] - num_hosp[i][age] -
                    num_rec[i][age] - num_death[i][age] -
                    rki_data[i][j]((size_t)InfectionState::ICU + (size_t)InfectionState::Count * age);
            }
        }
        date = offset_date_by_days(date, 1);
    }
    save_result(rki_data, region, path_join(results_dir, "Results_rki.h5"));

    auto rki_data_sum = epi::sum_nodes(std::vector<std::vector<TimeSeries<double>>>{rki_data});
    save_result({rki_data_sum[0][0]}, {0}, path_join(results_dir, "Results_rki_sum.h5"));
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
void read_population_data_germany(std::vector<Model>& model, Date date, std::vector<double>& scaling_factor_inf,
                                  double scaling_factor_icu, const std::string& dir)
{
    std::string id_name;
    std::vector<int> region(1, 0);
    if (date > Date(2020, 4, 23)) {
        details::set_divi_data(model, path_join(dir, "germany_divi.json"), id_name, {0}, date, scaling_factor_icu);
    }
    else {
        log_warning("No DIVI data available for this date");
    }
    details::set_rki_data(model, path_join(dir, "all_age_rki_ma.json"), id_name, {0}, date, scaling_factor_inf);
    details::set_population_data(model, path_join(dir, "county_current_population.json"), "ID_County", {0});
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
void read_population_data_state(std::vector<Model>& model, Date date, std::vector<int>& state,
                                std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                const std::string& dir)
{
    std::string id_name = "ID_State";
    if (date > Date(2020, 4, 23)) {
        details::set_divi_data(model, path_join(dir, "state_divi.json"), id_name, state, date, scaling_factor_icu);
    }
    else {
        log_warning("No DIVI data available for this date");
    }

    details::set_rki_data(model, path_join(dir, "all_state_age_rki_ma.json"), id_name, state, date, scaling_factor_inf);
    details::set_population_data(model, path_join(dir, "county_current_population.json"), "ID_County", state);
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
template <class Model>
void read_population_data_county(std::vector<Model>& model, Date date, std::vector<int> county,
                                 std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                 const std::string& dir)
{
    std::string id_name = "ID_County";

    if (date > Date(2020, 4, 23)) {
        details::set_divi_data(model, path_join(dir, "county_divi.json"), id_name, county, date, scaling_factor_icu);
    }
    else {
        log_warning("No DIVI data available for this date");
    }
    details::set_rki_data(model, path_join(dir, "all_county_age_rki_ma.json"), id_name, county, date,
                          scaling_factor_inf);
    details::set_population_data(model, path_join(dir, "county_current_population.json"), "ID_County", county);
}

/**
 * @brief returns a vector with the ids of all german counties
 * @param path directory to population data
 * @return
 */
std::vector<int> get_county_ids(const std::string& path);

} // namespace epi

#endif
