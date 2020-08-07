#ifndef SAVE_PARAMETERS_H
#define SAVE_PARAMETERS_H

#include <epidemiology/eigen_util.h>
#include <epidemiology/parameter_studies/parameter_studies.h>
#include <epidemiology/graph.h>
#include <epidemiology/migration.h>
#include <tixi.h>

namespace epi
{
/**
 * @brief read contact frequency matrix and damping distributions from xml file
 * @param handle Tixi Document Handle
 * @param path Path to contact frequency matrix Tree of XML file to read from
 */
ContactFrequencyVariableElement read_contact(TixiDocumentHandle handle, const std::string& path);

/**
 * @brief write contact frequency matrix and damping distributions to xml file
 * @param handle Tixi Document Handle
 * @param path Path to contact frequency matrix Tree of XML file
 * @param contact_freq_matrix Contact frequencies and dampings
 * @param nb_runs Number of runs of parameterstudy (used for predefinied samples of dampings)
 */
void write_contact(TixiDocumentHandle handle, const std::string& path,
                   const ContactFrequencyVariableElement& contact_freq_matrix, int nb_runs);

/**
 * @brief read parameter distribution from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter
 */
std::unique_ptr<ParameterDistribution> read_dist(TixiDocumentHandle handle, const std::string& path);

/**
 * @brief write parameter distribution to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter
 * @param element Name of parameter
 * @param dist Distribution of parameter
 */
void write_dist(const TixiDocumentHandle& handle, const std::string& path, const std::string& element,
                const ParameterDistribution& dist);

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
 */
ParameterSpace read_parameter_space(TixiDocumentHandle handle, const std::string& path);

/**
 * @brief write parameter space to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter space
 * @param parameter_space Parameter Space with distributions of all secir parameters
 * @param nb_runs Number of runs of parameterstudy (used for predefinied samples)
 */
void write_parameter_space(TixiDocumentHandle handle, const std::string& path, const ParameterSpace& parameter_space,
                           int nb_runs);

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
 */
void write_parameter_study(TixiDocumentHandle handle, const std::string& path, const ParameterStudy& parameter_study);

/**
 * @brief creates xml file with a single run parameter study with std 0 (used to save parameters of individual runs)
 * @param filename Name of file
 * @param params Secir parameters used during run
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
void write_single_run_params(const int run, const SecirParams& params, double t0, double tmax, std::vector<double> time,
                             std::vector<Eigen::VectorXd> secir_result);

/**
 * @brief Creates xml file containing Parameters of one node of a graph
 * @param graph Graph which holds the node
 * @param node Node ID
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
void write_node(const Graph<ModelNode<SecirSimulation>, MigrationEdge>& graph, int node, double t0, double tmax);

/**
 * @brief reads parameters of a single node and saves it into the graph
 * @param graph Graph in which the node is saved
 * @param node Node ID
 */
void read_node(Graph<ModelNode<SecirSimulation>, MigrationEdge>& graph, int node);

/**
 * @brief Writes the information of a single edge into a xml file
 * @param handle Tixi Document Handle
 * @param path Path to document root
 * @param graph Graph which holds the edge
 * @param edge Edge ID
 */
void write_edge(TixiDocumentHandle handle, const std::string& path,
                const Graph<ModelNode<SecirSimulation>, MigrationEdge>& graph, int edge);

/**
 * @brief Reads information of a single edge and saves it into the graph
 * @param handle Tixi document handle
 * @param path Path to document root
 * @param graph Graph to which the edge is added
 * @param edge Edge ID
 */
void read_edge(TixiDocumentHandle handle, const std::string& path,
               Graph<ModelNode<SecirSimulation>, MigrationEdge>& graph, int edge);

/**
 * @brief creates xml files for each node of a Secir simulation graph and one xml file for its edges
 * @param graph Graph which should be written
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
//TO-DO: Implement apropriate File System for XML FIles
void write_graph(const Graph<ModelNode<SecirSimulation>, MigrationEdge>& graph, double t0, double tmax);

/**
 * @brief reads graph xml files and returns a Secir simulation graph
 */
Graph<ModelNode<SecirSimulation>, MigrationEdge> read_graph();

} // namespace epi

#endif
