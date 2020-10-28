#ifndef SAVE_PARAMETERS_H
#define SAVE_PARAMETERS_H

#include <epidemiology/utils/eigen_util.h>
#include <epidemiology/utils/graph.h>
#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/parameter_studies.h>

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
SecirParams read_parameter_space(TixiDocumentHandle handle, const std::string& path, int io_mode);

/**
 * @brief write parameter space to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter space
 * @param params SecirParams Parameter Space with distributions of all secir parameters
 * @param num_runs Number of runs of parameter study (used for predefinied samples)
 * @param io_mode type of xml output (see epi::write_parameter_study for more details)
 */
void write_parameter_space(TixiDocumentHandle handle, const std::string& path, const SecirParams& params, int num_runs,
                           int io_mode);

/**
 * @brief read parameter study from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of parameters of study
 */
ParameterStudy read_parameter_study(TixiDocumentHandle handle, const std::string& path);

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
void write_parameter_study(TixiDocumentHandle handle, const std::string& path, const ParameterStudy& parameter_study,
                           int io_mode = 2);

/**
 * @brief creates xml file with a single run parameter study with std 0 (used to save parameters of individual runs)
 * @param filename Name of file
 * @param params Secir parameters used during run
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
void write_single_run_params(const int run, epi::Graph<epi::ModelNode<SecirSimulation>, MigrationEdge> graph, double t0,
                             double tmax);

/**
 * @brief Creates xml file containing Parameters of one node of a graph
 * @param handle Tixi document handle
 * @param graph Graph which holds the node
 * @param node Node ID
 */
void write_node(TixiDocumentHandle handle, const Graph<SecirParams, MigrationEdge>& graph, int node);

/**
 * @brief reads parameters of a single node and saves it into the graph
 * @param node_handle Tixi document handle
 * @param graph Graph in which the node is saved
 */
void read_node(TixiDocumentHandle node_handle, Graph<SecirParams, MigrationEdge>& graph);

/**
 * @brief Writes the information of a single edge into the xml file corresponding to its start node
 * @param edge_handles vecor of Tixi Document Handle. It's length is the number of nodes in the graph
 * @param path Path to document root
 * @param graph Graph which holds the edge
 * @param edge Edge ID
 */
void write_edge(const std::vector<TixiDocumentHandle>& edge_handles, const std::string& path,
                const Graph<SecirParams, MigrationEdge>& graph, int edge);

/**
 * @brief Reads information of a single edge and saves it into the graph
 * @param edge_handles vecor of Tixi Document Handle. It's length is the number of nodes in the graph
 * @param path Path to document root
 * @param graph Graph to which the edge is added
 * @param stort_node origin of the edge
 * @param end_node destination of the edge
 */
void read_edge(const std::vector<TixiDocumentHandle>& edge_handles, const std::string& path,
               Graph<SecirParams, MigrationEdge>& graph, int start_node, int end_node);

/**
 * @brief creates xml files for each node of a Secir simulation graph and one xml file for its edges for each node
 * @param graph Graph which should be written
 * @param dir_string directory, where graph should be stored
 */
void write_graph(const Graph<SecirParams, MigrationEdge>& graph, const std::string& dir_string);

/**
 * @brief reads graph xml files and returns a Secir simulation graph
 * @param dir_string directory from where graph should be read
 */
Graph<SecirParams, MigrationEdge> read_graph(const std::string& dir_string);

} // namespace epi

#endif
