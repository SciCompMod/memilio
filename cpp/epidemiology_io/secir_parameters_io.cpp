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
#include <epidemiology/utils/compiler_diagnostics.h>

#include <tixi.h>

#include <json/json.h>
#include <json/value.h>

#include <boost/filesystem.hpp>

#include <numeric>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <fstream>

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
        assert(element.get_distribution().get() && ("No Distribution detected for " + element_name +
                                                    ". Either define a distribution or choose a different io_mode.")
                                                       .c_str());
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
    ReturnCode status; 
    unused(status);

    if (io_mode == 0) {
        double read_buffer;
        status = tixiGetDoubleElement(handle, path.c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read value at " + path).c_str());
        value = std::make_unique<UncertainValue>(read_buffer);
    }
    else if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        std::unique_ptr<ParameterDistribution> distribution = read_distribution(handle, path);

        if (io_mode == 2) {
            double read_buffer;
            status = tixiGetDoubleElement(handle, path_join(path, "Value").c_str(), &read_buffer);
            assert(status == SUCCESS && ("Failed to read value at " + path).c_str());
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
    ReturnCode status; 
    unused(status);
    std::unique_ptr<ParameterDistribution> distribution;

    char* distri_str;
    status = tixiGetTextElement(handle, path_join(path, "Distribution").c_str(), &distri_str);
    assert(status == SUCCESS && ("Failed to read distribution type at " + path).c_str());
    if (strcmp("Normal", distri_str) == 0) {
        double mean;
        double dev;
        double min;
        double max;
        status = tixiGetDoubleElement(handle, path_join(path, "Mean").c_str(), &mean);
        assert(status == SUCCESS && ("Failed to read mean at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Deviation").c_str(), &dev);
        assert(status == SUCCESS && ("Failed to read deviation at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Min").c_str(), &min);
        assert(status == SUCCESS && ("Failed to read min value at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Max").c_str(), &max);
        assert(status == SUCCESS && ("Failed to read max value at " + path).c_str());

        distribution = std::make_unique<ParameterDistributionNormal>(min, max, mean, dev);
    }
    else if (strcmp("Uniform", distri_str) == 0) {
        double min;
        double max;
        status = tixiGetDoubleElement(handle, path_join(path, "Min").c_str(), &min);
        assert(status == SUCCESS && ("Failed to read min value at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Max").c_str(), &max);
        assert(status == SUCCESS && ("Failed to read max value at " + path).c_str());

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

void write_contact_frequency_matrix_collection(TixiDocumentHandle handle, const std::string& path,
                                               const ContactMatrixGroup& cfmc)
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

ContactMatrixGroup read_contact_frequency_matrix_collection(TixiDocumentHandle handle,
                                                                          const std::string& path)
{
    auto status = SUCCESS;
    unused(status);

    auto collection_path = path_join(path, "ContactMatrixGroup");
    int num_matrices;
    status = tixiGetNumberOfChilds(handle, collection_path.c_str(), &num_matrices);
    assert(status == SUCCESS && "Failed to read ContactMatrixGroup.");
    ContactMatrixGroup cfmc{1, size_t(num_matrices)};
    for (size_t i = 0; i < cfmc.get_num_matrices(); ++i) {
        auto cfm_path      = path_join(collection_path, "ContactMatrix" + std::to_string(i + 1));
        cfmc[i]            = ContactMatrix(read_matrix<>(handle, path_join(cfm_path, "Baseline")),
                                         read_matrix<>(handle, path_join(cfm_path, "Minimum")));
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
            cfmc[i].add_damping(read_matrix<>(handle, path_join(damping_path, "Values")), DampingLevel(level),
                                DampingType(type), SimulationTime(t));
        }
    }
    return cfmc;
}

void write_contact(TixiDocumentHandle handle, const std::string& path, const UncertainContactMatrix& contact_pattern,
                   int io_mode)
{
    tixiCreateElement(handle, path.c_str(), "ContactFreq");
    auto contact_path = path_join(path, "ContactFreq");

    write_contact_frequency_matrix_collection(handle, contact_path, contact_pattern.get_cont_freq_mat());

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
    UncertainContactMatrix contact_patterns{read_contact_frequency_matrix_collection(handle, path)};

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

void interpolate_ages(const std::vector<double>& age_ranges, const std::vector<double>& param_ranges,
                      std::vector<std::vector<double>>& interpolation, std::vector<bool>& carry_over)
{
    //counter for parameter age groups
    size_t counter = 0;

    //residual of param age groups
    double res = 0.0;
    for (size_t i = 0; i < age_ranges.size(); i++) {

        // if current param age group didn't fit into previous rki age group, transfer residual to current age group
        if (res < 0) {
            interpolation[i].push_back(std::min(-res / age_ranges[i], 1.0));
        }

        if (counter < param_ranges.size() - 1) {
            res += age_ranges[i];
            if (std::abs(res) < age_ranges[i]) {
                counter++;
            }
            // iterate over param age groups while there is still room in the current rki age group
            while (res > 0) {
                res -= param_ranges[counter];
                interpolation[i].push_back((param_ranges[counter] + std::min(res, 0.0)) / age_ranges[i]);
                if (res >= 0) {
                    counter++;
                }
            }
            if (res < 0) {
                carry_over.push_back(true);
            }
            else if (res == 0) {
                carry_over.push_back(false);
            }
        }
        // if last param age group is reached
        else {
            interpolation[i].push_back((age_ranges[i] + res) / age_ranges[i]);
            if (res < 0 || counter == 0) {
                carry_over.push_back(true);
            }
            else if (res == 0) {
                carry_over.push_back(false);
            }
            res = 0;
        }
    }
    // last entries for "unknown" age group
    interpolation.push_back({1.0});
    carry_over.push_back(true);
}

void read_rki_data(std::string const& path,
                   const std::string& id_name,
                   int region, int month, int day,
                   std::vector<double>& num_inf,
                   std::vector<double>& num_death,
                   std::vector<double>& num_rec)
{
    Json::Reader reader;
    Json::Value root;

    std::ifstream rki(path);
    reader.parse(rki, root);

    std::vector<std::string> age_names = {"A00-A04", "A05-A14", "A15-A34", "A35-A59", "A60-A79", "A80+", "unknown"};

    num_inf.resize(age_names.size());
    num_death.resize(age_names.size());
    num_rec.resize(age_names.size());

    for (size_t age = 0; age < age_names.size(); age++) {
        for (unsigned int i = 0; i < root.size(); i++) {
            bool correct_region = region == 0 || root[i][id_name] == region;
            std::string date    = root[i]["Date"].asString();
            if (month == std::stoi(date.substr(5, 2)) && day == std::stoi(date.substr(8, 2)) && correct_region) {
                if (root[i]["Age_RKI"].asString() == age_names[age]) {
                    num_inf[age]   = root[i]["Confirmed"].asDouble();
                    num_death[age] = root[i]["Deaths"].asDouble();
                    num_rec[age]   = root[i]["Recovered"].asDouble();
                    break;
                }
            }
        }
    }
}

double read_divi_data(const std::string& path,
                      const std::string& id_name,
                      int region, int month, int day)
{
    Json::Reader reader;
    Json::Value root;

    std::ifstream divi(path);
    reader.parse(divi, root);

    double num_icu = 0;
    for (unsigned int i = 0; i < root.size(); i++) {
        bool correct_region = region == 0 || root[i][id_name] == region;
        std::string date    = root[i]["Date"].asString();
        if (month == std::stoi(date.substr(5, 2)) && day == std::stoi(date.substr(8, 2)) && correct_region) {
            num_icu = root[i]["ICU"].asDouble();
        }
    }
    return num_icu;
}

std::vector<double> read_population_data(const std::string& path,
                                         const std::string& id_name,
                                         int region)
{
    Json::Reader reader;
    Json::Value root;

    std::ifstream census(path);
    reader.parse(census, root);

    std::vector<std::string> age_names = {"<3 years", "3-5 years", "6-14 years", "15-17 years", "18-24 years",
                                          "25-29 years", "30-39 years", "40-49 years", "50-64 years",
                                          "65-74 years", ">74 years"};

    std::vector<double> num_population(age_names.size(), 0.);

    for (size_t age = 0; age < age_names.size(); age++) {
        for (unsigned int i = 0; i < root.size(); i++) {
            bool correct_region = region == 0 || (int) root[i][id_name].asDouble()/1000 == region || root[i][id_name] == region ;
            if (correct_region) {
                num_population[age]   += root[i][age_names[age]].asDouble();
            }
        }
    }
    return num_population;
}


} // namespace epi
