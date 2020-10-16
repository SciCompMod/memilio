#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology_io/io.h>
#include <epidemiology/utils/memory.h>
#include <epidemiology/utils/uncertain_value.h>
#include <epidemiology/utils/stl_util.h>
#include <epidemiology/utils/graph.h>
#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/damping.h>
#include <epidemiology/model/populations.h>
#include <epidemiology/secir/uncertain_matrix.h>
#include <epidemiology/secir/secir.h>

#include <tixi.h>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <boost/filesystem.hpp>

namespace epi
{

void write_element(const TixiDocumentHandle& handle, const std::string& path, const std::string& element_name,
                   const UncertainValue& element, int io_mode, int num_runs)
{
    auto element_path = path_join(path, element_name);

    if (io_mode == 0) {
        tixiAddDoubleElement(handle, path.c_str(), element_name.c_str(), (double)element, "%g");
    }
    else if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        auto distribution = element.get_distribution().get();
        write_distribution(handle, path, element_name, *distribution);
        if (io_mode == 2) {
            tixiAddDoubleElement(handle, element_path.c_str(), "Value", (double)element, "%g");
        }
    }
    else {
        // TODO error handling
        epi::log_error("Wrong input for io_mode.");
    }

    if (io_mode == 3) {
        std::vector<double> predef_sample;
        for (int run = 0; run < num_runs; run++) {
            predef_sample.push_back((double)element);
        }
        write_predef_sample(handle, element_path, predef_sample);
    }
}

void write_distribution(const TixiDocumentHandle& handle, const std::string& path, const std::string& element_name,
                        const ParameterDistribution& distribution)
{

    struct WriteDistVisitor : public ConstParameterDistributionVisitor {
        WriteDistVisitor(const std::string& xml_path, TixiDocumentHandle tixi_handle)
            : handle(tixi_handle)
            , element_path(xml_path)
        {
        }

        void visit(const ParameterDistributionNormal& normal_distribution) override
        {
            tixiAddTextElement(handle, element_path.c_str(), "Distribution", "Normal");
            tixiAddDoubleElement(handle, element_path.c_str(), "Mean", normal_distribution.get_mean(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Deviation", normal_distribution.get_standard_dev(),
                                 "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Min", normal_distribution.get_lower_bound(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Max", normal_distribution.get_upper_bound(), "%g");
        }

        void visit(const ParameterDistributionUniform& uniform_distribution) override
        {
            tixiAddTextElement(handle, element_path.c_str(), "Distribution", "Uniform");
            tixiAddDoubleElement(handle, element_path.c_str(), "Min", uniform_distribution.get_lower_bound(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Max", uniform_distribution.get_upper_bound(), "%g");
        }

        TixiDocumentHandle handle;
        std::string element_path;
    };

    tixiCreateElement(handle, path.c_str(), element_name.c_str());
    auto element_path = path_join(path, element_name);

    WriteDistVisitor visitor(element_path, handle);
    distribution.accept(visitor);

    tixiAddFloatVector(handle, element_path.c_str(), "PredefinedSamples", distribution.get_predefined_samples().data(),
                       static_cast<int>(distribution.get_predefined_samples().size()), "%g");
}

std::unique_ptr<UncertainValue> read_element(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    std::unique_ptr<UncertainValue> value;

    if (io_mode == 0) {
        double read_buffer;
        tixiGetDoubleElement(handle, path.c_str(), &read_buffer);
        value = std::make_unique<UncertainValue>(read_buffer);
    }
    else if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        std::unique_ptr<ParameterDistribution> distribution = read_distribution(handle, path);

        if (io_mode == 2) {
            double read_buffer;
            tixiGetDoubleElement(handle, path_join(path, "Value").c_str(), &read_buffer);
            value = std::make_unique<UncertainValue>(read_buffer);
        }
        value->set_distribution(*distribution.get());
    }
    else {
        // TODO error handling
        epi::log_error("Wrong input for io_mode.");
    }
    return value;
}

std::unique_ptr<ParameterDistribution> read_distribution(TixiDocumentHandle handle, const std::string& path)
{
    std::unique_ptr<ParameterDistribution> distribution;

    char* distri_str;
    tixiGetTextElement(handle, path_join(path, "Distribution").c_str(), &distri_str);
    if (strcmp("Normal", distri_str) == 0) {
        double mean;
        double dev;
        double min;
        double max;
        tixiGetDoubleElement(handle, path_join(path, "Mean").c_str(), &mean);
        tixiGetDoubleElement(handle, path_join(path, "Deviation").c_str(), &dev);
        tixiGetDoubleElement(handle, path_join(path, "Min").c_str(), &min);
        tixiGetDoubleElement(handle, path_join(path, "Max").c_str(), &max);

        distribution = std::make_unique<ParameterDistributionNormal>(min, max, mean, dev);
    }
    else if (strcmp("Uniform", distri_str) == 0) {
        double min;
        double max;
        tixiGetDoubleElement(handle, path_join(path, "Min").c_str(), &min);
        tixiGetDoubleElement(handle, path_join(path, "Max").c_str(), &max);
        distribution = std::make_unique<ParameterDistributionUniform>(min, max);
    }
    else {
        // TODO: true error handling
        epi::log_error("Unknown distribution.");
        assert(false && "Unknown distribution.");
    }

    auto predef_path = path_join(path, "PredefinedSamples");
    int n_predef;
    tixiGetVectorSize(handle, predef_path.c_str(), &n_predef);

    double* predef = nullptr;
    tixiGetFloatVector(handle, predef_path.c_str(), &predef, n_predef);

    for (int i = 0; i < n_predef; i++) {
        distribution->add_predefined_sample(predef[i]);
    }

    return distribution;
}

void write_predef_sample(TixiDocumentHandle handle, const std::string& path, const std::vector<double>& samples)
{
    tixiRemoveElement(handle, path_join(path, "PredefinedSamples").c_str());
    tixiAddFloatVector(handle, path.c_str(), "PredefinedSamples", samples.data(), static_cast<int>(samples.size()),
                       "%g");
}

void write_contact(TixiDocumentHandle handle, const std::string& path, const UncertainContactMatrix& contact_pattern,
                   int io_mode)
{
    ContactFrequencyMatrix const& contact_freq_matrix = contact_pattern.get_cont_freq_mat();
    int num_groups                                    = contact_freq_matrix.get_size();
    tixiCreateElement(handle, path.c_str(), "ContactFreq");
    auto contact_path = path_join(path, "ContactFreq");
    for (int i = 0; i < num_groups; i++) {
        std::vector<double> row = {};
        for (int j = 0; j < num_groups; j++) {
            row.emplace_back(contact_freq_matrix.get_cont_freq(i, j));
        }
        tixiAddFloatVector(handle, contact_path.c_str(), ("ContactRateGroup_" + std::to_string(i + 1)).c_str(),
                           row.data(), num_groups, "%g");
    }
    for (int i = 0; i < num_groups; i++) {
        for (int j = 0; j < num_groups; j++) {
            int num_damp = static_cast<int>(contact_freq_matrix.get_dampings(i, j).get_dampings_vector().size());
            std::vector<double> row = {};
            for (int k = 0; k < num_damp; k++) {
                row.emplace_back(contact_freq_matrix.get_dampings(i, j).get_dampings_vector()[k].day);
                row.emplace_back(contact_freq_matrix.get_dampings(i, j).get_dampings_vector()[k].factor);
            }
            tixiAddFloatVector(handle, contact_path.c_str(),
                               ("DampingsGroups_" + std::to_string(i + 1) + "_" + std::to_string(j + 1)).c_str(),
                               row.data(), 2 * num_damp, "%g");
        }
    }

    if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        write_distribution(handle, contact_path, "NumDampings", *contact_pattern.get_distribution_damp_nb().get());
        write_distribution(handle, contact_path, "DampingDay", *contact_pattern.get_distribution_damp_days().get());
        write_distribution(handle, contact_path, "DampingDiagBase",
                           *contact_pattern.get_distribution_damp_diag_base().get());
        write_distribution(handle, contact_path, "DampingDiagRel",
                           *contact_pattern.get_distribution_damp_diag_rel().get());
        write_distribution(handle, contact_path, "DampingOffdiagRel",
                           *contact_pattern.get_distribution_damp_offdiag_rel().get());
    }
}

UncertainContactMatrix read_contact(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    int num_groups;
    tixiGetIntegerElement(handle, path_join("/Parameters", "NumberOfGroups").c_str(), &num_groups);
    UncertainContactMatrix contact_patterns{ContactFrequencyMatrix{(size_t)num_groups}};
    for (int i = 0; i < num_groups; i++) {
        double* row = nullptr;
        tixiGetFloatVector(handle, path_join(path, "ContactRateGroup_" + std::to_string(i + 1)).c_str(), &row,
                           num_groups);

        for (int j = 0; j < num_groups; ++j) {
            contact_patterns.get_cont_freq_mat().set_cont_freq(row[j], i, j);
        }
    }

    for (int i = 0; i < num_groups; i++) {
        for (int j = 0; j < num_groups; j++) {
            int num_dampings;
            tixiGetVectorSize(
                handle,
                path_join(path, ("DampingsGroups_" + std::to_string(i + 1) + "_" + std::to_string(j + 1))).c_str(),
                &num_dampings);
            double* dampings = nullptr;
            tixiGetFloatVector(
                handle,
                path_join(path, ("DampingsGroups_" + std::to_string(i + 1) + "_" + std::to_string(j + 1))).c_str(),
                &dampings, num_dampings);
            for (int k = 0; k < num_dampings / 2; k++) {
                contact_patterns.get_cont_freq_mat().add_damping(Damping{dampings[2 * k], dampings[2 * k + 1]}, i, j);
            }
        }
    }

    if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        contact_patterns.set_distribution_damp_nb(*read_distribution(handle, path_join(path, "NumDampings")));
        contact_patterns.set_distribution_damp_days(*read_distribution(handle, path_join(path, "DampingDay")));
        contact_patterns.set_distribution_damp_diag_base(
            *read_distribution(handle, path_join(path, "DampingDiagBase")));
        contact_patterns.set_distribution_damp_diag_rel(*read_distribution(handle, path_join(path, "DampingDiagRel")));
        contact_patterns.set_distribution_damp_offdiag_rel(
            *read_distribution(handle, path_join(path, "DampingOffdiagRel")));
    }
    return contact_patterns;
}

} // namespace epi
