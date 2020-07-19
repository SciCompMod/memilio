#ifndef SAVE_PARAMETERS_H
#define SAVE_PARAMETERS_H

#include <epidemiology/eigen_util.h>
#include <epidemiology/parameter_studies/parameter_studies.h>
#include <tixi.h>

namespace epi
{

ContactFrequencyVariableElement read_contact(TixiDocumentHandle handle, const std::string& path);
void write_contact(TixiDocumentHandle handle, const std::string& path,
                   const ContactFrequencyVariableElement& contact_freq_matrix, int nb_runs);

std::unique_ptr<ParameterDistribution> read_dist(TixiDocumentHandle handle, const std::string& path);
void write_dist(const TixiDocumentHandle& handle, const std::string& path, const std::string& element,
                const ParameterDistribution& dist);

void write_predef_sample(TixiDocumentHandle handle, const std::string& path, const std::vector<double>& samples);

ParameterSpace read_parameter_space(TixiDocumentHandle handle, const std::string& path);
void write_parameter_space(TixiDocumentHandle handle, const std::string& path, const ParameterSpace& parameter_space,
                           int nb_runs);

ParameterStudy read_parameter_study(TixiDocumentHandle handle, const std::string& path);
void write_parameter_study(TixiDocumentHandle handle, const std::string& path, const ParameterStudy& parameter_study);

void create_document(const std::string& filename, const ContactFrequencyMatrix& cont_freq,
                     const std::vector<SecirParams>& params, double t0, double tmax);

} // namespace epi

#endif
