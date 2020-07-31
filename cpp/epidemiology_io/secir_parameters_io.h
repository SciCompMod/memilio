#ifndef SAVE_PARAMETERS_H
#define SAVE_PARAMETERS_H

#include <epidemiology/eigen_util.h>
#include <epidemiology/parameter_studies/parameter_studies.h>
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
 * @param cont_freq Contact frequency Matrix used during run
 * @param params Secir parameters used during run
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
void write_single_run_params(const int run, const ContactFrequencyMatrix& cont_freq, const SecirParams& params,
                             double t0, double tmax, std::vector<double> time,
                             std::vector<Eigen::VectorXd> secir_result);

} // namespace epi

#endif
