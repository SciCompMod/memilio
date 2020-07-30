#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology/secir.h>
#include <epidemiology/damping.h>
#include <epidemiology/stl_util.h>
#include <vector>

#include <iostream>
#include <string>
#include <random>

#include <tixi.h>

namespace epi
{


void write_dist(const TixiDocumentHandle& handle, const std::string& path, const std::string& element,
                const ParameterDistribution& dist)
{

    struct WriteDistVisitor : public ParameterDistributionVisitor
    {
        WriteDistVisitor(const std::string& xml_path, TixiDocumentHandle tixi_handle)
            : handle(tixi_handle), element_path(xml_path)
        {}


        void visit(ParameterDistributionNormal & normal_dist) override
        {
            tixiAddTextElement(handle, element_path.c_str(), "Distribution", "Normal");
            tixiAddDoubleElement(handle, element_path.c_str(), "Mean", normal_dist.get_mean(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Deviation", normal_dist.get_standard_dev(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Min", normal_dist.get_lower_bound(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Max", normal_dist.get_upper_bound(), "%g");
        }

        void visit(ParameterDistributionUniform & uniform_dist) override
        {
            tixiAddTextElement(handle, element_path.c_str(), "Distribution", "Uniform");
            tixiAddDoubleElement(handle, element_path.c_str(), "Min", uniform_dist.get_lower_bound(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Max", uniform_dist.get_upper_bound(), "%g");
        }

        TixiDocumentHandle handle;
        std::string element_path;
    };


    tixiCreateElement(handle, path.c_str(), element.c_str());
    auto element_path = path_join(path, element);

    WriteDistVisitor visitor(element_path, handle);
    const_cast<ParameterDistribution&>(dist).accept(visitor);

    tixiAddFloatVector(handle, element_path.c_str(), "PredefinedSamples", dist.get_predefined_samples().data(),
                       dist.get_predefined_samples().size(), "%g");
}

std::unique_ptr<ParameterDistribution> read_dist(TixiDocumentHandle handle, const std::string& path)
{
    std::unique_ptr<ParameterDistribution> distribution;

    char* dist;
    tixiGetTextElement(handle, path_join(path, "Distribution").c_str(), &dist);
    if (strcmp("Normal", dist) == 0) {
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
    else if (strcmp("Uniform", dist) == 0) {
        double min;
        double max;
        tixiGetDoubleElement(handle, path_join(path, "Min").c_str(), &min);
        tixiGetDoubleElement(handle, path_join(path, "Max").c_str(), &max);
        distribution = std::make_unique<ParameterDistributionUniform>(min, max);
    }
    else {
        //TODO: true error handling
        assert(false && "Unknown distribution.");
    }

    auto predef_path = path_join(path, "PredefinedSamples");
    int n_predef;
    tixiGetVectorSize(handle, predef_path.c_str(), &n_predef);
    auto free = [](double* p) {
        std::free(p);
    };
    std::unique_ptr<double, decltype(free)> predef(
        [&]() {
            double* predef = nullptr;
            tixiGetFloatVector(handle, predef_path.c_str(), &predef, n_predef);
            return predef;
        }(),
        free);
    for (int i = 0; i < n_predef; i++) {
        distribution->add_predefined_sample(predef.get()[i]);
    }

    return distribution;
}

void write_predef_sample(TixiDocumentHandle handle, const std::string& path, const std::vector<double>& samples)
{
    tixiRemoveElement(handle, path_join(path, "PredefinedSamples").c_str());
    tixiAddFloatVector(handle, path.c_str(), "PredefinedSamples", samples.data(), samples.size(), "%g");
}

void write_contact(TixiDocumentHandle handle, const std::string& path,
                   const ContactFrequencyVariableElement& contact_freq_matrix, int nb_runs)
{
    ContactFrequencyMatrix cont_freq = contact_freq_matrix.get_cont_freq();
    int nb_groups                    = cont_freq.get_size();
    tixiCreateElement(handle, path.c_str(), "ContactFreq");
    auto contact_path = path_join(path, "ContactFreq");
    for (int i = 0; i < nb_groups; i++) {
        std::vector<double> row = {};
        for (int j = 0; j < nb_groups; j++) {
            row.emplace_back(cont_freq.get_cont_freq(i, j));
        }
        tixiAddFloatVector(handle, contact_path.c_str(), ("ContactRateGroup" + std::to_string(i + 1)).c_str(),
                           row.data(), nb_groups, "%g");
    }
    write_dist(handle, contact_path, "NumDampings", contact_freq_matrix.get_dist_num_dampings());
    write_dist(handle, contact_path, "DampingDay", contact_freq_matrix.get_dist_day());
    write_dist(handle, contact_path, "DampingDiagBase", contact_freq_matrix.get_dist_damp_diag_base());
    write_dist(handle, contact_path, "DampingDiagRel", contact_freq_matrix.get_dist_damp_diag_rel());
    write_dist(handle, contact_path, "DampingOffdiag", contact_freq_matrix.get_dist_damp_offdiag_rel());

    std::vector<double> predef_num_damp;
    std::vector<double> predef_day;
    std::vector<double> predef_diag_base;
    std::vector<double> predef_diag_rel;
    std::vector<double> predef_offdiag_rel;

    for (int i = 0; i < nb_runs; i++) {
        int nb_damp = cont_freq.get_dampings(0, 0).get_dampings_vector().size();
        predef_num_damp.push_back(nb_damp);
        for (int j = 0; j < nb_damp; j++) {
            predef_day.push_back(cont_freq.get_dampings(0, 0).get_dampings_vector().at(j).day);
            predef_diag_base.push_back(1.0);
            for (int k = 0; k < nb_groups; k++) {
                predef_diag_rel.push_back(cont_freq.get_dampings(k, k).get_dampings_vector().at(j).factor);
                for (int l = k + 1; l < nb_groups; l++) {
                    predef_offdiag_rel.push_back(cont_freq.get_dampings(k, l).get_dampings_vector().at(j).factor /
                                                 cont_freq.get_dampings(k, k).get_dampings_vector().at(j).factor);
                    predef_offdiag_rel.push_back(cont_freq.get_dampings(k, l).get_dampings_vector().at(j).factor /
                                                 cont_freq.get_dampings(l, l).get_dampings_vector().at(j).factor);
                }
            }
        }
    }

    write_predef_sample(handle, path_join(contact_path, "NumDampings"), predef_num_damp);
    write_predef_sample(handle, path_join(contact_path, "DampingDay"), predef_day);
    write_predef_sample(handle, path_join(contact_path, "DampingDiagBase"), predef_diag_base);
    write_predef_sample(handle, path_join(contact_path, "DampingDiagRel"), predef_diag_rel);
    write_predef_sample(handle, path_join(contact_path, "DampingOffdiag"), predef_offdiag_rel);
}

ContactFrequencyVariableElement read_contact(TixiDocumentHandle handle, const std::string& path)
{
    int nb_groups;
    tixiGetIntegerElement(handle, path_join("/Parameters", "NumberOfGroups").c_str(), &nb_groups);
    epi::ContactFrequencyMatrix contact_freq_matrix{(size_t)nb_groups};
    for (size_t i = 0; i < nb_groups; i++) {
        auto free = [](double* p) {
            std::free(p);
        };
        std::unique_ptr<double, decltype(free)> row(
            [&]() {
                double* row = nullptr;
                tixiGetFloatVector(handle, path_join(path, "ContactRateGroup" + std::to_string(i + 1)).c_str(), &row,
                                   nb_groups);
                return row;
            }(),
            free);
        for (int j = 0; j < nb_groups; ++j) {
            double temp = row.get()[j];
            contact_freq_matrix.set_cont_freq(row.get()[j], i, j);
        }
    }

    auto dist_n_damp =
        dynamic_unique_ptr_cast<ParameterDistributionUniform>(read_dist(handle, path_join(path, "NumDampings")));
    auto dist_damp_day =
        dynamic_unique_ptr_cast<ParameterDistributionUniform>(read_dist(handle, path_join(path, "DampingDay")));
    auto dist_damp_diag_base =
        dynamic_unique_ptr_cast<ParameterDistributionUniform>(read_dist(handle, path_join(path, "DampingDiagBase")));
    auto dist_damp_diag_rel =
        dynamic_unique_ptr_cast<ParameterDistributionUniform>(read_dist(handle, path_join(path, "DampingDiagRel")));
    auto dist_damp_offdiag =
        dynamic_unique_ptr_cast<ParameterDistributionUniform>(read_dist(handle, path_join(path, "DampingOffdiag")));

    return {contact_freq_matrix,           std::move(dist_n_damp),
            std::move(dist_damp_day),      std::move(dist_damp_diag_base),
            std::move(dist_damp_diag_rel), std::move(dist_damp_offdiag)};
}

ParameterStudy read_parameter_study(TixiDocumentHandle handle, const std::string& path)
{
    int n_runs;
    double t0;
    double tmax;

    tixiGetIntegerElement(handle, path_join(path, "Runs").c_str(), &n_runs);
    tixiGetDoubleElement(handle, path_join(path, "T0").c_str(), &t0);
    tixiGetDoubleElement(handle, path_join(path, "TMax").c_str(), &tmax);

    return ParameterStudy(&simulate, read_parameter_space(handle, path), n_runs, t0, tmax);
}

ParameterSpace read_parameter_space(TixiDocumentHandle handle, const std::string& path)
{
    int nb_groups;
    tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &nb_groups);

    auto cont_matrix =
        std::make_unique<ContactFrequencyVariableElement>(read_contact(handle, path_join(path, "ContactFreq")));

    //population
    std::vector<double> total_t0;
    std::vector<double> hosp_t0;
    std::vector<double> icu_t0;
    std::vector<double> dead_t0;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_exposed;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_carrier;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_infectious;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_recovered;

    //times
    std::vector<std::unique_ptr<ParameterDistribution>> dist_incubation;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_infectious_mild;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_diff_serial_incubation;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_hosp_to_recovered;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_infectious_to_hosp;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_infectious_asympt;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_hosp_to_icu;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_icu_to_recov;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_icu_to_dead;

    //probabilities
    std::vector<std::unique_ptr<ParameterDistribution>> dist_infect;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_asympt;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_risk_from_sympt;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_dead_per_icu;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_hosp_per_infectious;
    std::vector<std::unique_ptr<ParameterDistribution>> dist_icu_per_hospitalized;

    for (size_t i = 0; i < nb_groups; i++) {
        auto group_name = "Group" + std::to_string(i + 1);
        auto group_path = path_join(path, group_name);

        // populations
        auto population_path = path_join(group_path, "Population");

        total_t0.emplace_back();
        tixiGetDoubleElement(handle, path_join(population_path, "Total").c_str(), &total_t0.back());
        hosp_t0.emplace_back();
        tixiGetDoubleElement(handle, path_join(population_path, "Hospitalized").c_str(), &hosp_t0.back());
        icu_t0.emplace_back();
        tixiGetDoubleElement(handle, path_join(population_path, "ICU").c_str(), &icu_t0.back());
        dead_t0.emplace_back();
        tixiGetDoubleElement(handle, path_join(population_path, "Dead").c_str(), &dead_t0.back());
        dist_exposed.emplace_back(read_dist(handle, path_join(population_path, "Exposed")));
        dist_carrier.emplace_back(read_dist(handle, path_join(population_path, "Carrier")));
        dist_infectious.emplace_back(read_dist(handle, path_join(population_path, "Infectious")));
        dist_recovered.emplace_back(read_dist(handle, path_join(population_path, "Recovered")));

        // times
        auto times_path = path_join(group_path, "StageTimes");

        dist_incubation.emplace_back(read_dist(handle, path_join(times_path, "Incubation")));
        dist_infectious_mild.emplace_back(read_dist(handle, path_join(times_path, "InfectiousMild")));
        dist_diff_serial_incubation.emplace_back(read_dist(handle, path_join(times_path, "DiffIncubationSerial")));
        dist_hosp_to_recovered.emplace_back(read_dist(handle, path_join(times_path, "HospitalizedToRecovered")));
        dist_infectious_to_hosp.emplace_back(read_dist(handle, path_join(times_path, "InfectiousToHospitalized")));
        dist_infectious_asympt.emplace_back(read_dist(handle, path_join(times_path, "InfectiousAsympt")));
        dist_hosp_to_icu.emplace_back(read_dist(handle, path_join(times_path, "HospitalizedToICU")));
        dist_icu_to_recov.emplace_back(read_dist(handle, path_join(times_path, "ICUToRecovered")));
        dist_icu_to_dead.emplace_back(read_dist(handle, path_join(times_path, "ICUToDead")));

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");

        dist_infect.emplace_back(read_dist(handle, path_join(probabilities_path, "InfectedFromContact")));
        dist_asympt.emplace_back(read_dist(handle, path_join(probabilities_path, "AsympPerInfectious")));
        dist_risk_from_sympt.emplace_back(read_dist(handle, path_join(probabilities_path, "RiskFromSymptomatic")));
        dist_dead_per_icu.emplace_back(read_dist(handle, path_join(probabilities_path, "DeadPerICU")));
        dist_hosp_per_infectious.emplace_back(
            read_dist(handle, path_join(probabilities_path, "HospitalizedPerInfectious")));
        dist_icu_per_hospitalized.emplace_back(read_dist(handle, path_join(probabilities_path, "ICUPerHospitalized")));
    }

    return ParameterSpace(
        std::move(cont_matrix),
        //populations
        std::move(total_t0), std::move(dist_exposed), std::move(dist_carrier), std::move(dist_infectious),
        std::move(hosp_t0), std::move(icu_t0), std::move(dist_recovered), std::move(dead_t0),

        //times
        std::move(dist_incubation), std::move(dist_diff_serial_incubation), std::move(dist_infectious_mild),
        std::move(dist_hosp_to_recovered), std::move(dist_infectious_to_hosp), std::move(dist_infectious_asympt),
        std::move(dist_hosp_to_icu), std::move(dist_icu_to_recov), std::move(dist_icu_to_dead),

        //probabilities
        std::move(dist_infect), std::move(dist_asympt), std::move(dist_risk_from_sympt), std::move(dist_dead_per_icu),
        std::move(dist_hosp_per_infectious), std::move(dist_icu_per_hospitalized));
}

void write_parameter_space(TixiDocumentHandle handle, const std::string& path, const ParameterSpace& parameter_space,
                           int nb_runs)
{
    int nb_groups = parameter_space.get_total().size();
    tixiAddIntegerElement(handle, path.c_str(), "NumberOfGroups", nb_groups, "%d");

    for (int i = 0; i < nb_groups; i++) {
        auto group_name = "Group" + std::to_string(i + 1);
        auto group_path = path_join(path, group_name);

        tixiCreateElement(handle, path.c_str(), group_name.c_str());

        // populations
        auto population_path = path_join(group_path, "Population");
        tixiCreateElement(handle, group_path.c_str(), "Population");

        tixiAddDoubleElement(handle, population_path.c_str(), "Total", parameter_space.get_total()[i], "%g");
        tixiAddDoubleElement(handle, population_path.c_str(), "Hospitalized", parameter_space.get_hospitalized()[i],
                             "%g");
        tixiAddDoubleElement(handle, population_path.c_str(), "ICU", parameter_space.get_icu()[i], "%g");
        tixiAddDoubleElement(handle, population_path.c_str(), "Dead", parameter_space.get_dead()[i], "%g");
        write_dist(handle, population_path, "Exposed", parameter_space.get_dist_exposed(i));
        write_dist(handle, population_path, "Carrier", parameter_space.get_dist_carrier(i));
        write_dist(handle, population_path, "Infectious", parameter_space.get_dist_infectious(i));
        write_dist(handle, population_path, "Recovered", parameter_space.get_dist_recovered(i));

        // times
        auto times_path = path_join(group_path, "StageTimes");
        tixiCreateElement(handle, group_path.c_str(), "StageTimes");

        write_dist(handle, times_path, "Incubation", parameter_space.get_dist_incubation(i));
        write_dist(handle, times_path, "InfectiousMild", parameter_space.get_dist_inf_mild(i));
        write_dist(handle, times_path, "DiffIncubationSerial", parameter_space.get_dist_serial_int_incub_diff(i));
        write_dist(handle, times_path, "HospitalizedToRecovered", parameter_space.get_dist_hosp_to_rec(i));
        write_dist(handle, times_path, "InfectiousToHospitalized", parameter_space.get_dist_inf_to_hosp(i));
        write_dist(handle, times_path, "InfectiousAsympt", parameter_space.get_dist_inf_asymp(i));
        write_dist(handle, times_path, "HospitalizedToICU", parameter_space.get_dist_hosp_to_icu(i));
        write_dist(handle, times_path, "ICUToRecovered", parameter_space.get_dist_icu_to_rec(i));
        write_dist(handle, times_path, "ICUToDead", parameter_space.get_dist_icu_to_death(i));

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");
        tixiCreateElement(handle, group_path.c_str(), "Probabilities");

        write_dist(handle, probabilities_path, "InfectedFromContact", parameter_space.get_dist_inf_from_cont(i));
        write_dist(handle, probabilities_path, "AsympPerInfectious", parameter_space.get_dist_asymp_per_inf(i));
        write_dist(handle, probabilities_path, "RiskFromSymptomatic", parameter_space.get_dist_risk_from_symp(i));
        write_dist(handle, probabilities_path, "DeadPerICU", parameter_space.get_dist_death_per_icu(i));
        write_dist(handle, probabilities_path, "HospitalizedPerInfectious", parameter_space.get_dist_hosp_per_inf(i));
        write_dist(handle, probabilities_path, "ICUPerHospitalized", parameter_space.get_dist_icu_per_hosp(i));
    }

    write_contact(handle, path, parameter_space.get_cont_freq_matrix_variable(), nb_runs);
}

void write_parameter_study(TixiDocumentHandle handle, const std::string& path, const ParameterStudy& parameter_study)
{
    tixiAddIntegerElement(handle, path.c_str(), "Runs", parameter_study.get_nb_runs(), "%d");
    tixiAddDoubleElement(handle, path.c_str(), "T0", parameter_study.get_t0(), "%g");
    tixiAddDoubleElement(handle, path.c_str(), "TMax", parameter_study.get_tmax(), "%g");

    write_parameter_space(handle, path, parameter_study.get_parameter_space(), parameter_study.get_nb_runs());
}

void write_single_run_params(const int run, const ContactFrequencyMatrix& cont_freq,
                             const std::vector<SecirParams>& params, double t0, double tmax, std::vector<double> time,
                             std::vector<Eigen::VectorXd> secir_result)
{

    int nb_runs      = 1;
    std::string path = "/Parameters";
    TixiDocumentHandle handle;
    tixiCreateDocument("Parameters", &handle);

    ParameterStudy study(simulate, cont_freq, params, t0, tmax, 0.0, nb_runs);

    write_parameter_study(handle, path, study);
    tixiSaveDocument(handle, ("Parameters_run" + std::to_string(run) + ".xml").c_str());
    tixiCloseDocument(handle);

    save_result(time, secir_result, ("Results_run" + std::to_string(run) + ".h5"));
}

} // namespace epi
