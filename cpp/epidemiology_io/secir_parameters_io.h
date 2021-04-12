#ifndef SAVE_PARAMETERS_H
#define SAVE_PARAMETERS_H

#include <epidemiology/utils/eigen_util.h>
#include <epidemiology/utils/graph.h>
#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology_io/io.h>
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology/utils/date.h>

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
    template <class Model>
    void set_rki_data(std::vector<Model>& model, const std::string& path, const std::string& id_name,
                      std::vector<int> const& region, Date date, const std::vector<double>& scaling_factor_inf)
    {

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
                t_icu_to_dead[county].push_back(
                    static_cast<int>(model[county].parameters.times[group].get_icu_to_dead()));

                mu_C_R[county].push_back(model[county].parameters.probabilities[group].get_asymp_per_infectious());
                mu_I_H[county].push_back(
                    model[county].parameters.probabilities[group].get_hospitalized_per_infectious());
                mu_H_U[county].push_back(model[county].parameters.probabilities[group].get_icu_per_hospitalized());
            }
        }
        std::vector<std::vector<double>> num_inf(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_death(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_exp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_car(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_hosp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_icu(model.size(), std::vector<double>(age_ranges.size(), 0.0));

        read_rki_data(path, id_name, region, date, num_exp, num_car, num_inf, num_hosp, num_icu, num_death, num_rec,
                      t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec, t_inf_to_hosp, t_hosp_to_rec,
                      t_hosp_to_icu, t_icu_to_dead, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf);

        for (size_t county = 0; county < model.size(); county++) {
            if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
                size_t num_groups = model[county].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[county].populations[{(typename SecirModel<AgeGroup6>::AgeGroup)i, epi::InfectionType::E}] =
                        num_exp[county][i];
                    model[county].populations[{(typename SecirModel<AgeGroup6>::AgeGroup)i, epi::InfectionType::C}] =
                        num_car[county][i];
                    model[county].populations[{(typename SecirModel<AgeGroup6>::AgeGroup)i, epi::InfectionType::I}] =
                        num_inf[county][i];
                    model[county].populations[{(typename SecirModel<AgeGroup6>::AgeGroup)i, epi::InfectionType::H}] =
                        num_hosp[county][i];
                    model[county].populations[{(typename SecirModel<AgeGroup6>::AgeGroup)i, epi::InfectionType::D}] =
                        num_death[county][i];
                    model[county].populations[{(typename SecirModel<AgeGroup6>::AgeGroup)i, epi::InfectionType::R}] =
                        num_rec[county][i];
                }
            }
            else {
                log_warning("No infections reported on date " + std::to_string(date.year) + "-" +
                            std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                            std::to_string(region[county]) + ". Population data has not been set.");
            }
        }
    }

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
    template <class Model>
    void set_divi_data(std::vector<Model>& model, const std::string& path, const std::string& id_name,
                       const std::vector<int> vregion, Date date, double scaling_factor_icu)
    {
        std::vector<double> sum_mu_I_U(vregion.size(), 0);
        std::vector<std::vector<double>> mu_I_U{model.size()};
        for (size_t region = 0; region < vregion.size(); region++) {
            size_t num_groups = model[region].parameters.get_num_groups();
            for (size_t i = 0; i < num_groups; i++) {
                sum_mu_I_U[region] += model[region].parameters.probabilities[i].get_icu_per_hospitalized() *
                                      model[region].parameters.probabilities[i].get_hospitalized_per_infectious();
                mu_I_U[region].push_back(model[region].parameters.probabilities[i].get_icu_per_hospitalized() *
                                         model[region].parameters.probabilities[i].get_hospitalized_per_infectious());
            }
        }
        std::vector<double> num_icu(model.size(), 0.0);
        read_divi_data(path, id_name, vregion, date, num_icu);

        for (size_t region = 0; region < vregion.size(); region++) {
            size_t num_groups = model[region].parameters.get_num_groups();
            for (size_t i = 0; i < num_groups; i++) {
                model[region].populations[{(typename SecirModel<AgeGroup6>::AgeGroup)i, epi::InfectionType::U}] =
                    scaling_factor_icu * num_icu[region] * mu_I_U[region][i] / sum_mu_I_U[region];
            }
        }
    }

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
    template <class Model>
    void set_population_data(std::vector<Model>& model, const std::string& path, const std::string& id_name,
                             const std::vector<int> vregion)
    {
        std::vector<std::vector<double>> num_population = read_population_data(path, id_name, vregion);

        for (size_t region = 0; region < vregion.size(); region++) {
            if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
                size_t num_groups = model[region].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[region].populations.set_difference_from_group_total(
                        num_population[region][i], (typename Model::AgeGroup)i, (typename Model::AgeGroup)i,
                        epi::InfectionType::S);
                }
            }
            else {
                log_warning("No population data available for region " + std::to_string(region) +
                            ". Population data has not been set.");
            }
        }
    }
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
 * @brief read parameter space from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter space
 * @param io_mode type of xml input (see epi::write_parameter_study for more details)
 * @return returns a SecirParams object
 */
template <class AgeGroup>
SecirModel<AgeGroup> read_parameter_space(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    ReturnCode status;
    unused(status);

    int num_groups;
    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);
    assert(status == SUCCESS && ("Failed to read num_groups at " + path).c_str());

    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);

    if (num_groups != (int)AgeGroup::Count) {
        epi::log_error("Only 1, 2,3, 6 or 8 age groups allowed at the moment.");
    }

    SecirModel<AgeGroup> model;
    double read_buffer;
    status = tixiGetDoubleElement(handle, path_join(path, "StartDay").c_str(), &read_buffer);
    assert(status == SUCCESS && ("Failed to read StartDay at " + path).c_str());

    model.parameters.set_start_day(read_buffer);
    model.parameters.set_seasonality(*read_element(handle, path_join(path, "Seasonality"), io_mode));
    model.parameters.set_icu_capacity(*read_element(handle, path_join(path, "ICUCapacity"), io_mode));
    model.parameters.set_test_and_trace_capacity(
        *read_element(handle, path_join(path, "TestAndTraceCapacity"), io_mode));
    model.parameters.set_contact_patterns(read_contact(handle, path_join(path, "ContactFreq"), io_mode));

    for (size_t i = 0; i < static_cast<size_t>(num_groups); i++) {
        auto group_name = "Group" + std::to_string(i + 1);
        auto group_path = path_join(path, group_name);

        // populations
        auto population_path = path_join(group_path, "Population");

        status = tixiGetDoubleElement(handle, path_join(population_path, "Dead").c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read number of deaths at " + path).c_str());

        model.populations[{(AgeGroup)i, InfectionType::D}] = read_buffer;

        model.populations[{(AgeGroup)i, InfectionType::E}] =
            *read_element(handle, path_join(population_path, "Exposed"), io_mode);
        model.populations[{(AgeGroup)i, InfectionType::C}] =
            *read_element(handle, path_join(population_path, "Carrier"), io_mode);
        model.populations[{(AgeGroup)i, InfectionType::I}] =
            *read_element(handle, path_join(population_path, "Infectious"), io_mode);
        model.populations[{(AgeGroup)i, InfectionType::H}] =
            *read_element(handle, path_join(population_path, "Hospitalized"), io_mode);
        model.populations[{(AgeGroup)i, InfectionType::U}] =
            *read_element(handle, path_join(population_path, "ICU"), io_mode);
        model.populations[{(AgeGroup)i, InfectionType::R}] =
            *read_element(handle, path_join(population_path, "Recovered"), io_mode);

        status = tixiGetDoubleElement(handle, path_join(population_path, "Total").c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read total population at " + path).c_str());

        model.populations.set_difference_from_group_total(read_buffer, (AgeGroup)i, (AgeGroup)i, InfectionType::S);

        // times
        auto times_path = path_join(group_path, "StageTimes");

        model.parameters.times[i].set_incubation(*read_element(handle, path_join(times_path, "Incubation"), io_mode));
        model.parameters.times[i].set_infectious_mild(
            *read_element(handle, path_join(times_path, "InfectiousMild"), io_mode));
        model.parameters.times[i].set_serialinterval(
            *read_element(handle, path_join(times_path, "SerialInterval"), io_mode));
        model.parameters.times[i].set_hospitalized_to_home(
            *read_element(handle, path_join(times_path, "HospitalizedToRecovered"), io_mode));
        model.parameters.times[i].set_home_to_hospitalized(
            *read_element(handle, path_join(times_path, "InfectiousToHospitalized"), io_mode));
        model.parameters.times[i].set_infectious_asymp(
            *read_element(handle, path_join(times_path, "InfectiousAsympt"), io_mode));
        model.parameters.times[i].set_hospitalized_to_icu(
            *read_element(handle, path_join(times_path, "HospitalizedToICU"), io_mode));
        model.parameters.times[i].set_icu_to_home(
            *read_element(handle, path_join(times_path, "ICUToRecovered"), io_mode));
        model.parameters.times[i].set_icu_to_death(*read_element(handle, path_join(times_path, "ICUToDead"), io_mode));

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");

        model.parameters.probabilities[i].set_infection_from_contact(
            *read_element(handle, path_join(probabilities_path, "InfectedFromContact"), io_mode));
        model.parameters.probabilities[i].set_carrier_infectability(
            *read_element(handle, path_join(probabilities_path, "Carrierinfectability"), io_mode));
        model.parameters.probabilities[i].set_asymp_per_infectious(
            *read_element(handle, path_join(probabilities_path, "AsympPerInfectious"), io_mode));
        model.parameters.probabilities[i].set_risk_from_symptomatic(
            *read_element(handle, path_join(probabilities_path, "RiskFromSymptomatic"), io_mode));
        model.parameters.probabilities[i].set_test_and_trace_max_risk_from_symptomatic(
            *read_element(handle, path_join(probabilities_path, "TestAndTraceMaxRiskFromSymptomatic"), io_mode));
        model.parameters.probabilities[i].set_dead_per_icu(
            *read_element(handle, path_join(probabilities_path, "DeadPerICU"), io_mode));
        model.parameters.probabilities[i].set_hospitalized_per_infectious(
            *read_element(handle, path_join(probabilities_path, "HospitalizedPerInfectious"), io_mode));
        model.parameters.probabilities[i].set_icu_per_hospitalized(
            *read_element(handle, path_join(probabilities_path, "ICUPerHospitalized"), io_mode));
    }

    return model;
}

/**
 * @brief write parameter space to xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of the parameter space
 * @param params SecirParams Parameter Space with distributions of all secir parameters
 * @param num_runs Number of runs of parameter study (used for predefinied samples)
 * @param io_mode type of xml output (see epi::write_parameter_study for more details)
 */
template <class Model>
void write_parameter_space(TixiDocumentHandle handle, const std::string& path, Model const& model, int num_runs,
                           int io_mode)
{
    auto num_groups = model.parameters.get_num_groups();
    tixiAddIntegerElement(handle, path.c_str(), "NumberOfGroups", (int)num_groups, "%d");

    tixiAddDoubleElement(handle, path.c_str(), "StartDay", model.parameters.get_start_day(), "%g");
    write_element(handle, path, "Seasonality", model.parameters.get_seasonality(), io_mode, num_runs);
    write_element(handle, path, "ICUCapacity", model.parameters.get_icu_capacity(), io_mode, num_runs);
    write_element(handle, path, "TestAndTraceCapacity", model.parameters.get_test_and_trace_capacity(), io_mode,
                  num_runs);

    for (size_t i = 0; i < num_groups; i++) {
        auto group_name = "Group" + std::to_string(i + 1);
        auto group_path = path_join(path, group_name);

        tixiCreateElement(handle, path.c_str(), group_name.c_str());

        // populations
        auto population_path = path_join(group_path, "Population");
        tixiCreateElement(handle, group_path.c_str(), "Population");

        tixiAddDoubleElement(handle, population_path.c_str(), "Total",
                             model.populations.get_group_total((typename Model::AgeGroup)i), "%g");
        tixiAddDoubleElement(handle, population_path.c_str(), "Dead",
                             model.populations[{(typename Model::AgeGroup)i, InfectionType::D}], "%g");
        write_element(handle, population_path, "Exposed",
                      model.populations[{(typename Model::AgeGroup)i, InfectionType::E}], io_mode, num_runs);
        write_element(handle, population_path, "Carrier",
                      model.populations[{(typename Model::AgeGroup)i, InfectionType::C}], io_mode, num_runs);
        write_element(handle, population_path, "Infectious",
                      model.populations[{(typename Model::AgeGroup)i, InfectionType::I}], io_mode, num_runs);
        write_element(handle, population_path, "Hospitalized",
                      model.populations[{(typename Model::AgeGroup)i, InfectionType::H}], io_mode, num_runs);
        write_element(handle, population_path, "ICU",
                      model.populations[{(typename Model::AgeGroup)i, InfectionType::U}], io_mode, num_runs);
        write_element(handle, population_path, "Recovered",
                      model.populations[{(typename Model::AgeGroup)i, InfectionType::R}], io_mode, num_runs);

        // times
        auto times_path = path_join(group_path, "StageTimes");
        tixiCreateElement(handle, group_path.c_str(), "StageTimes");

        write_element(handle, times_path, "Incubation", model.parameters.times[i].get_incubation(), io_mode, num_runs);
        write_element(handle, times_path, "InfectiousMild", model.parameters.times[i].get_infectious_mild(), io_mode,
                      num_runs);
        write_element(handle, times_path, "SerialInterval", model.parameters.times[i].get_serialinterval(), io_mode,
                      num_runs);
        write_element(handle, times_path, "HospitalizedToRecovered",
                      model.parameters.times[i].get_hospitalized_to_home(), io_mode, num_runs);
        write_element(handle, times_path, "InfectiousToHospitalized",
                      model.parameters.times[i].get_home_to_hospitalized(), io_mode, num_runs);
        write_element(handle, times_path, "InfectiousAsympt", model.parameters.times[i].get_infectious_asymp(), io_mode,
                      num_runs);
        write_element(handle, times_path, "HospitalizedToICU", model.parameters.times[i].get_hospitalized_to_icu(),
                      io_mode, num_runs);
        write_element(handle, times_path, "ICUToRecovered", model.parameters.times[i].get_icu_to_home(), io_mode,
                      num_runs);
        write_element(handle, times_path, "ICUToDead", model.parameters.times[i].get_icu_to_dead(), io_mode, num_runs);

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");
        tixiCreateElement(handle, group_path.c_str(), "Probabilities");

        write_element(handle, probabilities_path, "InfectedFromContact",
                      model.parameters.probabilities[i].get_infection_from_contact(), io_mode, num_runs);
        write_element(handle, probabilities_path, "Carrierinfectability",
                      model.parameters.probabilities[i].get_carrier_infectability(), io_mode, num_runs);
        write_element(handle, probabilities_path, "AsympPerInfectious",
                      model.parameters.probabilities[i].get_asymp_per_infectious(), io_mode, num_runs);
        write_element(handle, probabilities_path, "RiskFromSymptomatic",
                      model.parameters.probabilities[i].get_risk_from_symptomatic(), io_mode, num_runs);
        write_element(handle, probabilities_path, "TestAndTraceMaxRiskFromSymptomatic",
                      model.parameters.probabilities[i].get_test_and_trace_max_risk_from_symptomatic(), io_mode,
                      num_runs);
        write_element(handle, probabilities_path, "DeadPerICU", model.parameters.probabilities[i].get_dead_per_icu(),
                      io_mode, num_runs);
        write_element(handle, probabilities_path, "HospitalizedPerInfectious",
                      model.parameters.probabilities[i].get_hospitalized_per_infectious(), io_mode, num_runs);
        write_element(handle, probabilities_path, "ICUPerHospitalized",
                      model.parameters.probabilities[i].get_icu_per_hospitalized(), io_mode, num_runs);
    }

    write_contact(handle, path, model.parameters.get_contact_patterns(), io_mode);
}

/**
 * @brief read parameter study from xml file
 * @param handle Tixi Document Handle
 * @param path Path to XML Tree of parameters of study
 */
template <class AgeGroup>
ParameterStudy<SecirModel<AgeGroup>> read_parameter_study(TixiDocumentHandle handle, const std::string& path)
{
    ReturnCode status;

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

    unused(status);

    SecirModel<AgeGroup> model = read_parameter_space<AgeGroup>(handle, path, io_mode);
    return ParameterStudy<SecirModel<AgeGroup>>(model, t0, tmax, num_runs);
}

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
template <class Model>
void write_parameter_study(TixiDocumentHandle handle, const std::string& path,
                           const ParameterStudy<Model>& parameter_study, int io_mode = 2)
{
    tixiAddIntegerElement(handle, path.c_str(), "IOMode", io_mode, "%d");
    tixiAddIntegerElement(handle, path.c_str(), "Runs", parameter_study.get_num_runs(), "%d");
    tixiAddDoubleElement(handle, path.c_str(), "T0", parameter_study.get_t0(), "%g");
    tixiAddDoubleElement(handle, path.c_str(), "TMax", parameter_study.get_tmax(), "%g");

    write_parameter_space<Model>(handle, path, parameter_study.get_model(), parameter_study.get_num_runs(), io_mode);
}

/**
 * @brief creates xml file with a single run parameter study with std 0 (used to save parameters of individual runs)
 * @param filename Name of file
 * @param params Secir parameters used during run
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
template <class Model>
void write_single_run_params(const int run,
                             epi::Graph<epi::ModelNode<epi::Simulation<Model>>, epi::MigrationEdge> graph, double t0,
                             double tmax)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");

    std::string abs_path;
    bool created = create_directory("results", abs_path);

    if (created) {
        log_info("Results are stored in {:s}/results.", epi::get_current_dir_name());
    }
    else if (run == 0) {
        log_info(
            "Directory '{:s}' already exists. Results are stored in {:s}/results. Files from previous runs will be "
            "overwritten",
            epi::get_current_dir_name());
    }
    std::vector<TimeSeries<double>> all_results;
    std::vector<int> ids;

    ids.reserve(graph.nodes().size());
    all_results.reserve(graph.nodes().size());
    std::transform(graph.nodes().begin(), graph.nodes().end(), std::back_inserter(all_results), [](auto& node) {
        return node.property.get_result();
    });
    std::transform(graph.nodes().begin(), graph.nodes().end(), std::back_inserter(ids), [](auto& node) {
        return node.id;
    });

    int node_id = 0;
    for (auto& node : graph.nodes()) {
        int num_runs     = 1;
        std::string path = "/Parameters";
        TixiDocumentHandle handle;
        tixiCreateDocument("Parameters", &handle);
        ParameterStudy<Model> study(node.property.get_simulation().get_model(), t0, tmax, num_runs);

        write_parameter_study(handle, path, study);

        tixiSaveDocument(handle, path_join(abs_path, ("Parameters_run" + std::to_string(run) + "_node" +
                                                      std::to_string(node_id) + ".xml"))
                                     .c_str());
        tixiCloseDocument(handle);
        node_id++;
    }

    save_result(all_results, ids,
                path_join(abs_path, ("Results_run" + std::to_string(run) + std::to_string(node_id) + ".h5")));
}

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
template <class AgeGroup>
void read_node(TixiDocumentHandle node_handle, Graph<SecirModel<AgeGroup>, MigrationParameters>& graph)
{
    int node_id;
    tixiGetIntegerElement(node_handle, path_join("/Parameters", "NodeID").c_str(), &node_id);
    graph.add_node(node_id, read_parameter_space<AgeGroup>(node_handle, "/Parameters", 2));
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
    int num_groups  = static_cast<int>(graph.nodes()[0].property.parameters.get_num_groups());
    int num_compart = static_cast<int>(graph.nodes()[0].property.populations.get_num_compartments()) / num_groups;

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

    auto all_weights = Eigen::VectorXd(num_compart * num_groups);
    for (int group = 0; group < num_groups; group++) {
        double* weights = nullptr;
        status = tixiGetFloatVector(handle, path_join(edge_path, "Group" + std::to_string(group + 1)).c_str(), &weights,
                                    num_compart);
        assert(status == SUCCESS && "Failed to read coefficients.");
        for (int compart = 0; compart < num_compart; compart++) {
            all_weights(compart + group * num_compart) = weights[compart];
        }
    }
    graph.add_edge(start_node, end_node, all_weights);
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
        read_node(node_handle, graph);
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
