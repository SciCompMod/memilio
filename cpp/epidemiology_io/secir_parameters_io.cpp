#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology_io/io.h>
#include <epidemiology/utils/memory.h>
#include <epidemiology/utils/uncertain_value.h>
#include <epidemiology/utils/stl_util.h>
#include <epidemiology/utils/graph.h>
#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/damping.h>
#include <epidemiology/secir/populations.h>
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

ParameterStudy read_parameter_study(TixiDocumentHandle handle, const std::string& path)
{
    ReturnCode status; 
    unused(status);

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

    return ParameterStudy(read_parameter_space(handle, path, io_mode), t0, tmax, num_runs);
}

SecirParams read_parameter_space(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    ReturnCode status; 
    unused(status);

    int num_groups;
    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);
    assert(status == SUCCESS && ("Failed to read num_groups at " + path).c_str());

    SecirParams params{(size_t)num_groups};

    double read_buffer;
    status = tixiGetDoubleElement(handle, path_join(path, "StartDay").c_str(), &read_buffer);
    assert(status == SUCCESS && ("Failed to read StartDay at " + path).c_str());

    params.set_start_day(read_buffer);
    params.set_seasonality(*read_element(handle, path_join(path, "Seasonality"), io_mode));
    params.set_icu_capacity(*read_element(handle, path_join(path, "ICUCapacity"), io_mode));

    params.set_contact_patterns(read_contact(handle, path_join(path, "ContactFreq"), io_mode));

    for (size_t i = 0; i < static_cast<size_t>(num_groups); i++) {
        auto group_name = "Group" + std::to_string(i + 1);
        auto group_path = path_join(path, group_name);

        // populations
        auto population_path = path_join(group_path, "Population");

        status = tixiGetDoubleElement(handle, path_join(population_path, "Dead").c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read number of deaths at " + path).c_str());

        params.populations.set({i, SecirCompartments::D}, read_buffer);

        params.populations.set({i, SecirCompartments::E},
                               *read_element(handle, path_join(population_path, "Exposed"), io_mode));
        params.populations.set({i, SecirCompartments::C},
                               *read_element(handle, path_join(population_path, "Carrier"), io_mode));
        params.populations.set({i, SecirCompartments::I},
                               *read_element(handle, path_join(population_path, "Infectious"), io_mode));
        params.populations.set({i, SecirCompartments::H},
                               *read_element(handle, path_join(population_path, "Hospitalized"), io_mode));
        params.populations.set({i, SecirCompartments::U},
                               *read_element(handle, path_join(population_path, "ICU"), io_mode));
        params.populations.set({i, SecirCompartments::R},
                               *read_element(handle, path_join(population_path, "Recovered"), io_mode));

        status = tixiGetDoubleElement(handle, path_join(population_path, "Total").c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read total population at " + path).c_str());

        params.populations.set_difference_from_group_total({i, SecirCompartments::S}, epi::SecirCategory::AgeGroup, i,
                                                           read_buffer);

        // times
        auto times_path = path_join(group_path, "StageTimes");

        params.times[i].set_incubation(*read_element(handle, path_join(times_path, "Incubation"), io_mode));
        params.times[i].set_infectious_mild(*read_element(handle, path_join(times_path, "InfectiousMild"), io_mode));
        params.times[i].set_serialinterval(*read_element(handle, path_join(times_path, "SerialInterval"), io_mode));
        params.times[i].set_hospitalized_to_home(
            *read_element(handle, path_join(times_path, "HospitalizedToRecovered"), io_mode));
        params.times[i].set_home_to_hospitalized(
            *read_element(handle, path_join(times_path, "InfectiousToHospitalized"), io_mode));
        params.times[i].set_infectious_asymp(*read_element(handle, path_join(times_path, "InfectiousAsympt"), io_mode));
        params.times[i].set_hospitalized_to_icu(
            *read_element(handle, path_join(times_path, "HospitalizedToICU"), io_mode));
        params.times[i].set_icu_to_home(*read_element(handle, path_join(times_path, "ICUToRecovered"), io_mode));
        params.times[i].set_icu_to_death(*read_element(handle, path_join(times_path, "ICUToDead"), io_mode));

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");

        params.probabilities[i].set_infection_from_contact(
            *read_element(handle, path_join(probabilities_path, "InfectedFromContact"), io_mode));
        params.probabilities[i].set_carrier_infectability(
            *read_element(handle, path_join(probabilities_path, "Carrierinfectability"), io_mode));
        params.probabilities[i].set_asymp_per_infectious(
            *read_element(handle, path_join(probabilities_path, "AsympPerInfectious"), io_mode));
        params.probabilities[i].set_risk_from_symptomatic(
            *read_element(handle, path_join(probabilities_path, "RiskFromSymptomatic"), io_mode));
        params.probabilities[i].set_dead_per_icu(
            *read_element(handle, path_join(probabilities_path, "DeadPerICU"), io_mode));
        params.probabilities[i].set_hospitalized_per_infectious(
            *read_element(handle, path_join(probabilities_path, "HospitalizedPerInfectious"), io_mode));
        params.probabilities[i].set_icu_per_hospitalized(
            *read_element(handle, path_join(probabilities_path, "ICUPerHospitalized"), io_mode));
    }

    return params;
}

void write_parameter_space(TixiDocumentHandle handle, const std::string& path, const SecirParams& parameters,
                           int num_runs, int io_mode)
{
    auto num_groups = parameters.get_num_groups();
    tixiAddIntegerElement(handle, path.c_str(), "NumberOfGroups", static_cast<int>(num_groups), "%d");

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

        tixiAddDoubleElement(handle, population_path.c_str(), "Total",
                             parameters.populations.get_group_total(AgeGroup, i), "%g");
        tixiAddDoubleElement(handle, population_path.c_str(), "Dead",
                             parameters.populations.get({i, SecirCompartments::D}), "%g");
        write_element(handle, population_path, "Exposed", parameters.populations.get({i, SecirCompartments::E}),
                      io_mode, num_runs);
        write_element(handle, population_path, "Carrier", parameters.populations.get({i, SecirCompartments::C}),
                      io_mode, num_runs);
        write_element(handle, population_path, "Infectious", parameters.populations.get({i, SecirCompartments::I}),
                      io_mode, num_runs);
        write_element(handle, population_path, "Hospitalized", parameters.populations.get({i, SecirCompartments::H}),
                      io_mode, num_runs);
        write_element(handle, population_path, "ICU", parameters.populations.get({i, SecirCompartments::U}), io_mode,
                      num_runs);
        write_element(handle, population_path, "Recovered", parameters.populations.get({i, SecirCompartments::R}),
                      io_mode, num_runs);

        // times
        auto times_path = path_join(group_path, "StageTimes");
        tixiCreateElement(handle, group_path.c_str(), "StageTimes");

        write_element(handle, times_path, "Incubation", parameters.times[i].get_incubation(), io_mode, num_runs);
        write_element(handle, times_path, "InfectiousMild", parameters.times[i].get_infectious_mild(), io_mode,
                      num_runs);
        write_element(handle, times_path, "SerialInterval", parameters.times[i].get_serialinterval(), io_mode,
                      num_runs);
        write_element(handle, times_path, "HospitalizedToRecovered", parameters.times[i].get_hospitalized_to_home(),
                      io_mode, num_runs);
        write_element(handle, times_path, "InfectiousToHospitalized", parameters.times[i].get_home_to_hospitalized(),
                      io_mode, num_runs);
        write_element(handle, times_path, "InfectiousAsympt", parameters.times[i].get_infectious_asymp(), io_mode,
                      num_runs);
        write_element(handle, times_path, "HospitalizedToICU", parameters.times[i].get_hospitalized_to_icu(), io_mode,
                      num_runs);
        write_element(handle, times_path, "ICUToRecovered", parameters.times[i].get_icu_to_home(), io_mode, num_runs);
        write_element(handle, times_path, "ICUToDead", parameters.times[i].get_icu_to_dead(), io_mode, num_runs);

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");
        tixiCreateElement(handle, group_path.c_str(), "Probabilities");

        write_element(handle, probabilities_path, "InfectedFromContact",
                      parameters.probabilities[i].get_infection_from_contact(), io_mode, num_runs);
        write_element(handle, probabilities_path, "Carrierinfectability",
                      parameters.probabilities[i].get_carrier_infectability(), io_mode, num_runs);
        write_element(handle, probabilities_path, "AsympPerInfectious",
                      parameters.probabilities[i].get_asymp_per_infectious(), io_mode, num_runs);
        write_element(handle, probabilities_path, "RiskFromSymptomatic",
                      parameters.probabilities[i].get_risk_from_symptomatic(), io_mode, num_runs);
        write_element(handle, probabilities_path, "DeadPerICU", parameters.probabilities[i].get_dead_per_icu(), io_mode,
                      num_runs);
        write_element(handle, probabilities_path, "HospitalizedPerInfectious",
                      parameters.probabilities[i].get_hospitalized_per_infectious(), io_mode, num_runs);
        write_element(handle, probabilities_path, "ICUPerHospitalized",
                      parameters.probabilities[i].get_icu_per_hospitalized(), io_mode, num_runs);
    }

    write_contact(handle, path, parameters.get_contact_patterns(), io_mode);
}

void write_parameter_study(TixiDocumentHandle handle, const std::string& path, const ParameterStudy& parameter_study,
                           int io_mode)
{
    tixiAddIntegerElement(handle, path.c_str(), "IOMode", io_mode, "%d");
    tixiAddIntegerElement(handle, path.c_str(), "Runs", parameter_study.get_num_runs(), "%d");
    tixiAddDoubleElement(handle, path.c_str(), "T0", parameter_study.get_t0(), "%g");
    tixiAddDoubleElement(handle, path.c_str(), "TMax", parameter_study.get_tmax(), "%g");

    write_parameter_space(handle, path, parameter_study.get_secir_params(), parameter_study.get_num_runs(), io_mode);
}

void write_single_run_params(const int run, epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge> graph,
                             double t0, double tmax)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");
    boost::filesystem::path dir("results");

    bool created = boost::filesystem::create_directory(dir);

    if (created) {
        log_info("Directory '{:s}' was created. Results are stored in {:s}/results.", dir.string(),
                 epi::get_current_dir_name());
    }
    else if (run == 0) {
        log_info(
            "Directory '{:s}' already exists. Results are stored in {:s}/ results. Files from previous runs will be "
            "overwritten",
            dir.string(), epi::get_current_dir_name());
    }
    int node_id = 0;
    for (auto& node : graph.nodes()) {
        int num_runs     = 1;
        std::string path = "/Parameters";
        TixiDocumentHandle handle;
        tixiCreateDocument("Parameters", &handle);
        ParameterStudy study(node.get_params(), t0, tmax, num_runs);

        write_parameter_study(handle, path, study);

        tixiSaveDocument(handle, path_join(dir.string(), ("Parameters_run" + std::to_string(run) + "_node" +
                                                          std::to_string(node_id) + ".xml"))
                                     .c_str());
        tixiCloseDocument(handle);

        save_result(node.get_result(), path_join(dir.string(), ("Results_run" + std::to_string(run) + "_node" +
                                                                std::to_string(node_id) + ".h5")));
        node_id++;
    }
}

void write_node(TixiDocumentHandle handle, const Graph<SecirParams, MigrationEdge>& graph, int node)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");
    int num_runs = 1;
    int io_mode  = 2;

    std::string path = "/Parameters";

    tixiAddIntegerElement(handle, path.c_str(), "NodeID", node, "%d");

    auto params = graph.nodes()[node];

    write_parameter_space(handle, path, params, num_runs, io_mode);
}

void read_node(TixiDocumentHandle node_handle, Graph<SecirParams, MigrationEdge>& graph)
{

    graph.add_node(read_parameter_space(node_handle, "/Parameters", 2));
}

void write_edge(const std::vector<TixiDocumentHandle>& edge_handles, const std::string& path,
                const Graph<SecirParams, MigrationEdge>& graph, int edge)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");
    int num_groups  = static_cast<int>(graph.nodes()[0].get_num_groups());
    int num_compart = static_cast<int>(graph.nodes()[0].populations.get_num_compartments()) / num_groups;

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

void read_edge(const std::vector<TixiDocumentHandle>& edge_handles, const std::string& path,
               Graph<SecirParams, MigrationEdge>& graph, int start_node, int end_node)
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

void write_graph(const Graph<SecirParams, MigrationEdge>& graph, const std::string& dir_string)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");

    boost::filesystem::path dir(dir_string);
    bool created = boost::filesystem::create_directory(dir);

    if (created) {
        log_info("Directory '{:s}' was created. Results are stored in {:s}.", dir.string(),
                 path_join(epi::get_current_dir_name(), dir.string()));
    }
    else {
        log_info("Directory '{:s}' already exists. Results are stored in {:s}. Files from previous "
                 "graph will be "
                 "overwritten",
                 dir.string(), path_join(epi::get_current_dir_name(), dir.string()));
    }
    int num_nodes   = static_cast<int>(graph.nodes().size());
    int num_edges   = static_cast<int>(graph.edges().size());
    int num_groups  = static_cast<int>(graph.nodes()[0].get_contact_patterns().get_cont_freq_mat().get_num_groups());
    int num_compart = static_cast<int>(graph.nodes()[0].populations.get_num_compartments()) / num_groups;

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
                         (dir / ("GraphEdges_node" + std::to_string(node) + ".xml")).string().c_str());
        tixiCloseDocument(edge_handles[node]);
    }

    for (int node = 0; node < num_nodes; node++) {
        TixiDocumentHandle node_handle;
        tixiCreateDocument("Parameters", &node_handle);
        write_node(node_handle, graph, node);
        tixiSaveDocument(node_handle, path_join(dir.string(), ("GraphNode" + std::to_string(node) + ".xml")).c_str());
        tixiCloseDocument(node_handle);
    }
}

Graph<SecirParams, MigrationEdge> read_graph(const std::string& dir_string)
{
    boost::filesystem::path dir(dir_string);
    assert(boost::filesystem::exists(dir) && ("Directory " + dir_string + " does not exist.").c_str());

    ReturnCode status; 
    unused(status);
    TixiDocumentHandle handle;
    tixiOpenDocument(path_join(dir.string(), "GraphEdges_node0.xml").c_str(), &handle);

    std::string edges_path = "/Edges";

    int num_nodes;
    int num_edges;

    status = tixiGetIntegerElement(handle, path_join(edges_path, "NumberOfNodes").c_str(), &num_nodes);
    assert(status == SUCCESS && ("Failed to read num_nodes at " + edges_path).c_str());

    status = tixiGetIntegerElement(handle, path_join(edges_path, "NumberOfEdges").c_str(), &num_edges);
    assert(status == SUCCESS && ("Failed to read num_edges at " + edges_path).c_str());

    std::vector<TixiDocumentHandle> edge_handles(num_nodes);

    Graph<SecirParams, MigrationEdge> graph;

    for (int node = 0; node < num_nodes; node++) {
        TixiDocumentHandle node_handle;
        tixiOpenDocument((dir / ("GraphNode" + std::to_string(node) + ".xml")).string().c_str(), &node_handle);
        read_node(node_handle, graph);
        tixiCloseDocument(node_handle);
    }

    for (int start_node = 0; start_node < num_nodes; start_node++) {
        tixiOpenDocument((dir / ("GraphEdges_node" + std::to_string(start_node) + ".xml")).string().c_str(),
                         &edge_handles[start_node]);
        for (int end_node = 0; end_node < num_nodes; end_node++) {
            read_edge(edge_handles, edges_path, graph, start_node, end_node);
        }

        tixiCloseDocument(edge_handles[start_node]);
    }
    return graph;
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

void set_rki_data(epi::SecirParams& params, const std::vector<double>& param_ranges, const std::string& path,
                  const std::string& id_name, int region, int month, int day)
{
    Json::Reader reader;
    Json::Value root;

    std::ifstream rki(path);
    reader.parse(rki, root);

    std::vector<std::string> age_names = {"A00-A04", "A05-A14", "A15-A34", "A35-A59", "A60-A79", "A80+", "unknown"};
    std::vector<double> age_ranges     = {5., 10., 20., 25., 20., 20.};

    std::vector<std::vector<double>> interpolation(age_names.size());
    std::vector<bool> carry_over;

    interpolate_ages(age_ranges, param_ranges, interpolation, carry_over);

    std::vector<double> num_inf(age_names.size(), 0.0);
    std::vector<double> num_death(age_names.size(), 0.0);
    std::vector<double> num_rec(age_names.size(), 0.0);

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

    std::vector<double> interpol_inf(params.get_num_groups() + 1, 0.0);
    std::vector<double> interpol_death(params.get_num_groups() + 1, 0.0);
    std::vector<double> interpol_rec(params.get_num_groups() + 1, 0.0);

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

    for (size_t i = 0; i < params.get_num_groups(); i++) {
        interpol_inf[i] += (double)num_inf[num_inf.size() - 1] / (double)params.get_num_groups();
        interpol_death[i] += (double)num_death[num_death.size() - 1] / (double)params.get_num_groups();
        interpol_rec[i] += (double)num_rec[num_rec.size() - 1] / (double)params.get_num_groups();
    }

    if (std::accumulate(num_inf.begin(), num_inf.end(), 0.0) > 0) {
        size_t num_groups = params.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            params.populations.set({i, epi::SecirCompartments::I},
                                   interpol_inf[i] - interpol_death[i] - interpol_rec[i]);
            params.populations.set({i, epi::SecirCompartments::D}, interpol_death[i]);
            params.populations.set({i, epi::SecirCompartments::R}, interpol_rec[i]);
        }
    }
    else {
        log_warning("No infections reported on date " + std::to_string(day) + "-" + std::to_string(month) +
                    " for region " + std::to_string(region) + ". Population data has not been set.");
    }
}

void set_divi_data(epi::SecirParams& params, const std::string& path, const std::string& id_name, int region, int month,
                   int day)
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

    if (num_icu > 0) {
        size_t num_groups = params.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            params.populations.set({i, epi::SecirCompartments::U}, num_icu / (double)num_groups);
        }
    }
    else {
        log_warning("No ICU patients reported on date " + std::to_string(day) + "-" + std::to_string(month) +
                    " for region " + std::to_string(region) + ".");
    }
}

void set_population_data(epi::SecirParams& params, const std::vector<double>& param_ranges, const std::string& path,
                         const std::string& id_name, int region)
{
    Json::Reader reader;
    Json::Value root;

    std::ifstream census(path);
    reader.parse(census, root);

    std::vector<std::string> age_names = {"<3 years", "3-5 years", "6-14 years", "15-17 years", "18-24 years",
                                          "25-29 years", "30-39 years", "40-49 years", "50-64 years",
                                          "65-74 years", ">74 years"};
    std::vector<double> age_ranges     = {3., 3., 9., 3., 7., 5., 10., 10., 15., 10., 25.};

    std::vector<std::vector<double>> interpolation(age_names.size());
    std::vector<bool> carry_over;

    interpolate_ages(age_ranges, param_ranges, interpolation, carry_over);


    std::vector<double> num_population(age_names.size(), 0.0);

    for (size_t age = 0; age < age_names.size(); age++) {
        for (unsigned int i = 0; i < root.size(); i++) {
            bool correct_region = region == 0 || (int) root[i][id_name].asDouble()/1000 == region || root[i][id_name] == region ;
            if (correct_region) {
                num_population[age]   += root[i][age_names[age]].asDouble();
            }
        }
    }

    std::vector<double> interpol_population(params.get_num_groups() + 1, 0.0);

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
        size_t num_groups = params.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            params.populations.set_difference_from_group_total({i, epi::SecirCompartments::S}, epi::SecirCategory::AgeGroup,
                                                                i, interpol_population[i]);
        }
    }
    else {
        log_warning("No population data available for region " + std::to_string(region) + ". Population data has not been set.");
    }

}

void read_population_data_germany(epi::SecirParams& params, const std::vector<double>& param_ranges, int month, int day,
                                  const std::string& dir)
{
    assert(param_ranges.size() == params.get_num_groups() &&
           "size of param_ranges needs to be the same size as the number of groups in params");
    assert(std::accumulate(param_ranges.begin(), param_ranges.end(), 0.0) == 100. && "param_ranges must add up to 100");

    std::string id_name;

    set_rki_data(params, param_ranges, path_join(dir, "all_age_rki.json"), id_name, 0, month, day);
    set_divi_data(params, path_join(dir, "germany_divi.json"), id_name, 0, month, day);
    set_population_data(params, param_ranges, path_join(dir, "county_current_population.json"), "ID_County", 0);
}

void read_population_data_state(epi::SecirParams& params, const std::vector<double>& param_ranges, int month, int day,
                                int state, const std::string& dir)
{
    assert(state > 0 && state <= 16 && "State must be between 1 and 16");
    assert(param_ranges.size() == params.get_num_groups() &&
           "size of param_ranges needs to be the same size as the number of groups in params");
    assert(std::accumulate(param_ranges.begin(), param_ranges.end(), 0.0) == 100. && "param_ranges must add up to 100");

    std::string id_name = "ID_State";

    set_rki_data(params, param_ranges, path_join(dir, "all_state_age_rki.json"), id_name, state, month, day);
    set_divi_data(params, path_join(dir, "state_divi.json"), id_name, state, month, day);
    set_population_data(params, param_ranges, path_join(dir, "county_current_population.json"), "ID_County", state);
}

void read_population_data_county(epi::SecirParams& params, const std::vector<double>& param_ranges, int month, int day,
                                 int county, const std::string& dir)
{
    assert(county > 999 && "State must be between 1 and 16");
    assert(param_ranges.size() == params.get_num_groups() &&
           "size of param_ranges needs to be the same size as the number of groups in params");
    assert(std::accumulate(param_ranges.begin(), param_ranges.end(), 0.0) == 100. && "param_ranges must add up to 100");

    std::string id_name = "ID_County";

    set_rki_data(params, param_ranges, path_join(dir, "all_county_age_rki_concated_berlin.json"), id_name, county, month, day);
    set_divi_data(params, path_join(dir, "county_divi.json"), id_name, county, month, day);
    set_population_data(params, param_ranges, path_join(dir, "county_current_population.json"), "ID_County", county);
}

} // namespace epi
