#include <epidemiology/save_parameters.h>
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

void write_dist(const TixiDocumentHandle& handle, const std::string& path, const ParameterDistribution& dist)
{
    switch (dist.get_distribution()) {
    case DIST_NORMAL: {
        tixiAddTextElement(handle, path.c_str(), "Distribution", "Normal");
        auto& normal_dist = static_cast<const ParameterDistributionNormal&>(dist);
        tixiAddDoubleElement(handle, path.c_str(), "Mean", normal_dist.get_mean(), "%g");
        tixiAddDoubleElement(handle, path.c_str(), "Deviation", normal_dist.get_standard_dev(), "%g");
        tixiAddDoubleElement(handle, path.c_str(), "Min", normal_dist.get_lower_bound(), "%g");
        tixiAddDoubleElement(handle, path.c_str(), "Max", normal_dist.get_upper_bound(), "%g");
        break;
    }
    case DIST_UNIFORM: {
        tixiAddTextElement(handle, path.c_str(), "Distribution", "Uniform");
        auto& uniform_dist = static_cast<const ParameterDistributionUniform&>(dist);
        tixiAddDoubleElement(handle, path.c_str(), "Min", uniform_dist.get_lower_bound(), "%g");
        tixiAddDoubleElement(handle, path.c_str(), "Max", uniform_dist.get_upper_bound(), "%g");
        break;
    }
    default:
        //TODO: true error handling
        assert(false && "Unknown distribution.");
        break;
    }
    tixiAddFloatVector(handle, path.c_str(), "PredefinedSamples", dist.get_predefined_samples().data(), dist.get_predefined_samples().size(), "%g");
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
    auto free = [](double* p) { std::free(p); };
    std::unique_ptr<double, decltype(free)> predef([&](){
        double* predef;
        tixiGetFloatVector(handle, predef_path.c_str(), &predef, n_predef);
        return predef;
    }(), free);
    for (int i = 0; i < n_predef; i++)
    {
        distribution->add_predefined_sample(predef.get()[i]);
    }

    return distribution;
}

ContactFrequencyVariableElement read_contact(TixiDocumentHandle handle, const std::string& path)
{
    int nb_groups;
	tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &nb_groups);

	epi::ContactFrequencyMatrix contact_freq_matrix{ (size_t)nb_groups };    
    for (size_t i = 0; i < nb_groups; i++) {
        auto free = [](double* p) { std::free(p); };
        std::unique_ptr<double, decltype(free)> row([&]() {
                double* row;
                tixiGetFloatVector(handle, path_join(path, "ContactRateGroup" + std::to_string(i + 1)).c_str(),
                                    &row, nb_groups);
                return row;
            }(), free);
        for (int j = 0; j < nb_groups; ++j) {
            contact_freq_matrix.set_cont_freq(row.get()[j], i, j);
        }
    }

    auto dist_n_damp =
        dynamic_unique_ptr_cast<ParameterDistributionUniform>(read_dist(handle, path_join(path, "NumDampings")));
    auto dist_damp_day =
        dynamic_unique_ptr_cast<ParameterDistributionUniform>(read_dist(handle, path_join(path, "DampingDay")));
    auto dist_damp_diag_base = dynamic_unique_ptr_cast<ParameterDistributionUniform>(
        read_dist(handle, path_join(path, "DampingDiagBase")));
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

    return ParameterStudy(&simulate, read_parameter_space(handle, path_join(path, "Parameters")), n_runs);
}

ParameterSpace read_parameter_space(TixiDocumentHandle handle, const std::string& path)
{
	int nb_groups;
	tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &nb_groups);
	
    auto cont_matrix = std::make_unique<ContactFrequencyVariableElement>(read_contact(handle, path_join(path, "ContactFreq")));

    //population
    auto total_t0 = std::vector<double>(nb_groups);
    auto hosp_t0 = std::vector<double>(nb_groups);
    auto icu_t0 = std::vector<double>(nb_groups);
    auto dead_t0 = std::vector<double>(nb_groups);
    auto dist_exposed = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_carrier = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_infectious = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_recovered = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);

    //times
    auto dist_incubation = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_infectious_mild = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_diff_serial_incubation = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_hosp_to_recovered = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_infectious_to_hosp = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups); 
    auto dist_infectious_asympt = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups); 
    auto dist_hosp_to_icu = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_icu_to_recov = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_icu_to_dead = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);

    //probabilities
    auto dist_infect = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_asympt = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_risk_from_sympt = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_dead_per_icu = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_hosp_per_infectious = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);
    auto dist_icu_per_hospitalized = std::vector<std::unique_ptr<ParameterDistribution>>(nb_groups);

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

		dist_infect.emplace_back(read_dist(handle, path_join(probabilities_path, "Infected")));
        dist_asympt.emplace_back(read_dist(handle, path_join(probabilities_path, "Asympt")));
        dist_risk_from_sympt.emplace_back(read_dist(handle, path_join(probabilities_path, "RiskFromSymptomatic")));
        dist_dead_per_icu.emplace_back(read_dist(handle, path_join(probabilities_path, "DeadPerICU")));
        dist_hosp_per_infectious.emplace_back(read_dist(handle, path_join(probabilities_path, "HospitalizedPerInfectious")));
        dist_icu_per_hospitalized.emplace_back(read_dist(handle, path_join(probabilities_path, "ICUPerHospitalized")));
    }

    return ParameterSpace(
        std::move(cont_matrix),
        //populations
        std::move(total_t0), std::move(dist_exposed), std::move(dist_carrier), std::move(dist_infectious), 
        std::move(hosp_t0), std::move(icu_t0), std::move(dist_recovered), std::move(dead_t0),

        //times
        std::move(dist_incubation), std::move(dist_infectious_mild), std::move(dist_diff_serial_incubation),
        std::move(dist_hosp_to_recovered), std::move(dist_infectious_to_hosp), std::move(dist_infectious_asympt),
        std::move(dist_hosp_to_icu), std::move(dist_icu_to_recov), std::move(dist_icu_to_dead),

        //probabilities
        std::move(dist_infect), std::move(dist_asympt), std::move(dist_risk_from_sympt), std::move(dist_dead_per_icu),
        std::move(dist_hosp_per_infectious), std::move(dist_icu_per_hospitalized));
}

}
