/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kuehn, Lena Ploetzke
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef ODESECIR_PARAMETERS_IO_H
#define ODESECIR_PARAMETERS_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secir/model.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"

namespace mio
{

namespace osecir
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
     * @param[in] path Path to RKI file
     * @param[in] rki_data Vector of ConfirmedCasesDataEntry%s  
     * @param[in] region vector of keys of the region of interest
     * @param[in] year Specifies year at which the data is read
     * @param[in] month Specifies month at which the data is read
     * @param[in] day Specifies day at which the data is read
     * @param[in, out] num_* output vector for number of people in the corresponding compartement
     * @param[in] t_* vector average time it takes to get from one compartement to another for each age group
     * @param[in] mu_* vector probabilities to get from one compartement to another for each age group
     */
IOResult<void> read_confirmed_cases_data(
    std::string const& path, std::vector<ConfirmedCasesDataEntry>& rki_data, std::vector<int> const& vregion, Date date,
    std::vector<std::vector<double>>& vnum_Exposed, std::vector<std::vector<double>>& vnum_InfectedNoSymptoms,
    std::vector<std::vector<double>>& vnum_InfectedSymptoms, std::vector<std::vector<double>>& vnum_InfectedSevere,
    std::vector<std::vector<double>>& vnum_icu, std::vector<std::vector<double>>& vnum_death,
    std::vector<std::vector<double>>& vnum_rec, const std::vector<std::vector<int>>& vt_Exposed,
    const std::vector<std::vector<int>>& vt_InfectedNoSymptoms,
    const std::vector<std::vector<int>>& vt_InfectedSymptoms, const std::vector<std::vector<int>>& vt_InfectedSevere,
    const std::vector<std::vector<int>>& vt_InfectedCritical, const std::vector<std::vector<double>>& vmu_C_R,
    const std::vector<std::vector<double>>& vmu_I_H, const std::vector<std::vector<double>>& vmu_H_U,
    const std::vector<double>& scaling_factor_inf);

/**
     * @brief sets populations data from json file with multiple age groups into a Model with one age group
     * @tparam FP floating point data type, e.g., double
     * @param[in, out] model vector of objects in which the data is set
     * @param[in] path Path to confirmed cases file
     * @param[in] region vector of keys of the region of interest
     * @param[in] year Specifies year at which the data is read
     * @param[in] month Specifies month at which the data is read
     * @param[in] day Specifies day at which the data is read
     * @param[in] scaling_factor_inf factors by which to scale the confirmed cases of rki data
     */
template <typename FP = double>
IOResult<void> set_confirmed_cases_data(std::vector<Model<FP>>& model, const std::string& path,
                                        std::vector<int> const& region, Date date,
                                        const std::vector<double>& scaling_factor_inf)
{
    const size_t num_age_groups = ConfirmedCasesDataEntry::age_group_names.size();
    assert(scaling_factor_inf.size() == num_age_groups);

    std::vector<std::vector<int>> t_InfectedNoSymptoms{model.size()};
    std::vector<std::vector<int>> t_Exposed{model.size()};
    std::vector<std::vector<int>> t_InfectedSymptoms{model.size()};
    std::vector<std::vector<int>> t_InfectedSevere{model.size()};
    std::vector<std::vector<int>> t_InfectedCritical{model.size()};

    std::vector<std::vector<double>> mu_C_R{model.size()};
    std::vector<std::vector<double>> mu_I_H{model.size()};
    std::vector<std::vector<double>> mu_H_U{model.size()};
    std::vector<std::vector<double>> mu_U_D{model.size()};

    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path));

    for (size_t node = 0; node < model.size(); ++node) {
        for (size_t group = 0; group < num_age_groups; group++) {

            t_Exposed[node].push_back(
                static_cast<int>(std::round(model[node].parameters.template get<TimeExposed<FP>>()[(AgeGroup)group])));
            t_InfectedNoSymptoms[node].push_back(static_cast<int>(
                std::round(model[node].parameters.template get<TimeInfectedNoSymptoms<FP>>()[(AgeGroup)group])));
            t_InfectedSymptoms[node].push_back(static_cast<int>(
                std::round(model[node].parameters.template get<TimeInfectedSymptoms<FP>>()[(AgeGroup)group])));
            t_InfectedSevere[node].push_back(static_cast<int>(
                std::round(model[node].parameters.template get<TimeInfectedSevere<FP>>()[(AgeGroup)group])));
            t_InfectedCritical[node].push_back(static_cast<int>(
                std::round(model[node].parameters.template get<TimeInfectedCritical<FP>>()[(AgeGroup)group])));

            mu_C_R[node].push_back(
                model[node].parameters.template get<RecoveredPerInfectedNoSymptoms<FP>>()[(AgeGroup)group]);
            mu_I_H[node].push_back(
                model[node].parameters.template get<SeverePerInfectedSymptoms<FP>>()[(AgeGroup)group]);
            mu_H_U[node].push_back(model[node].parameters.template get<CriticalPerSevere<FP>>()[(AgeGroup)group]);
            mu_U_D[node].push_back(model[node].parameters.template get<DeathsPerCritical<FP>>()[(AgeGroup)group]);
        }
    }
    std::vector<std::vector<double>> num_InfectedSymptoms(model.size(), std::vector<double>(num_age_groups, 0.0));
    std::vector<std::vector<double>> num_death(model.size(), std::vector<double>(num_age_groups, 0.0));
    std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(num_age_groups, 0.0));
    std::vector<std::vector<double>> num_Exposed(model.size(), std::vector<double>(num_age_groups, 0.0));
    std::vector<std::vector<double>> num_InfectedNoSymptoms(model.size(), std::vector<double>(num_age_groups, 0.0));
    std::vector<std::vector<double>> num_InfectedSevere(model.size(), std::vector<double>(num_age_groups, 0.0));
    std::vector<std::vector<double>> num_icu(model.size(), std::vector<double>(num_age_groups, 0.0));

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(path, case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));

    for (size_t node = 0; node < model.size(); node++) {
        if (std::accumulate(num_InfectedSymptoms[node].begin(), num_InfectedSymptoms[node].end(), 0.0) > 0) {
            size_t num_groups = (size_t)model[node].parameters.get_num_groups();
            for (size_t i = 0; i < num_groups; i++) {
                model[node].populations[{AgeGroup(i), InfectionState::Exposed}] = num_Exposed[node][i];
                model[node].populations[{AgeGroup(i), InfectionState::InfectedNoSymptoms}] =
                    num_InfectedNoSymptoms[node][i];
                model[node].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsConfirmed}] = 0;
                model[node].populations[{AgeGroup(i), InfectionState::InfectedSymptoms}] =
                    num_InfectedSymptoms[node][i];
                model[node].populations[{AgeGroup(i), InfectionState::InfectedSymptomsConfirmed}] = 0;
                model[node].populations[{AgeGroup(i), InfectionState::InfectedSevere}] = num_InfectedSevere[node][i];
                model[node].populations[{AgeGroup(i), InfectionState::Dead}]           = num_death[node][i];
                model[node].populations[{AgeGroup(i), InfectionState::Recovered}]      = num_rec[node][i];
            }
        }
        else {
            log_warning("No infections reported on date " + std::to_string(date.year) + "-" +
                        std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                        std::to_string(region[node]) + ". Population data has not been set.");
        }
    }
    return success();
}

/**
     * @brief reads number of ICU patients from DIVI register into Parameters
     * @param path Path to DIVI file
     * @param vregion Keys of the region of interest
     * @param year Specifies year at which the data is read
     * @param month Specifies month at which the data is read
     * @param day Specifies day at which the data is read
     * @param vnum_icu number of ICU patients
     */
IOResult<void> read_divi_data(const std::string& path, const std::vector<int>& vregion, Date date,
                              std::vector<double>& vnum_icu);

/**
     * @brief sets populations data from DIVI register into Model
     * @tparam FP floating point data type, e.g., double
     * @param model vector of objects in which the data is set
     * @param path Path to DIVI file
     * @param vregion vector of keys of the regions of interest
     * @param year Specifies year at which the data is read
     * @param month Specifies month at which the data is read
     * @param day Specifies day at which the data is read
     * @param scaling_factor_icu factor by which to scale the icu cases of divi data
     */
template <typename FP = double>
IOResult<void> set_divi_data(std::vector<Model<FP>>& model, const std::string& path, const std::vector<int>& vregion,
                             Date date, double scaling_factor_icu)
{
    std::vector<double> sum_mu_I_U(vregion.size(), 0);
    std::vector<std::vector<double>> mu_I_U{model.size()};
    for (size_t region = 0; region < vregion.size(); region++) {
        auto num_groups = model[region].parameters.get_num_groups();
        for (auto i = AgeGroup(0); i < num_groups; i++) {
            sum_mu_I_U[region] += model[region].parameters.template get<CriticalPerSevere<FP>>()[i] *
                                  model[region].parameters.template get<SeverePerInfectedSymptoms<FP>>()[i];
            mu_I_U[region].push_back(model[region].parameters.template get<CriticalPerSevere<FP>>()[i] *
                                     model[region].parameters.template get<SeverePerInfectedSymptoms<FP>>()[i]);
        }
    }
    std::vector<double> num_icu(model.size(), 0.0);
    BOOST_OUTCOME_TRY(read_divi_data(path, vregion, date, num_icu));

    for (size_t region = 0; region < vregion.size(); region++) {
        auto num_groups = model[region].parameters.get_num_groups();
        for (auto i = AgeGroup(0); i < num_groups; i++) {
            model[region].populations[{i, InfectionState::InfectedCritical}] =
                scaling_factor_icu * num_icu[region] * mu_I_U[region][(size_t)i] / sum_mu_I_U[region];
        }
    }

    return success();
}

/**
     * @brief Reads population data from census data
     * @tparam FP floating point data type, e.g., double
     * @param[in] path Path to RKI file
     * @param[in] vregion Vector of keys of the regions of interest
     * @param[in] accumulate_age_groups Specifies whether population data sould be accumulated to one age group
     */
IOResult<std::vector<std::vector<double>>>
read_population_data(const std::string& path, const std::vector<int>& vregion, bool accumulate_age_groups = false);

/**
     * @brief sets population data from census data
     * @tparam FP floating point data type, e.g., double
     * @param[in, out] model vector of objects in which the data is set
     * @param[in] path Path to RKI file
     * @param[in] vregion vector of keys of the regions of interest
     * @param[in] accumulate_age_groups specifies whether population data sould be accumulated to one age group
     */
template <typename FP = double>
IOResult<void> set_population_data(std::vector<Model<FP>>& model, const std::string& path,
                                   const std::vector<int>& vregion, bool accumulate_age_groups = false)
{
    BOOST_OUTCOME_TRY(auto&& num_population, read_population_data(path, vregion, accumulate_age_groups));

    for (size_t region = 0; region < vregion.size(); region++) {
        if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
            auto num_groups = model[region].parameters.get_num_groups();
            for (auto i = AgeGroup(0); i < num_groups; i++) {
                model[region].populations.template set_difference_from_group_total<AgeGroup>(
                    {i, InfectionState::Susceptible}, num_population[region][size_t(i)]);
            }
        }
        else {
            log_warning("No population data available for region " + std::to_string(region) +
                        ". Population data has not been set.");
        }
    }

    return success();
}

} //namespace details

#ifdef MEMILIO_HAS_HDF5

/**
* @brief sets populations data from RKI into a Model with one age group
* @param model vector of models in which the data is set. Copy is made to avoid changing the original model
* @param data_dir Path to RKI files
* @param results_dir Path to result files
* @param region vector of keys of the region of interest
* @param year Specifies year at which the data is read
* @param month Specifies month at which the data is read
* @param day Specifies day at which the data is read
* @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
*/
template <class Model>
IOResult<void>
export_input_data_county_timeseries(std::vector<Model> models, const std::string& dir, std::vector<int> const& region,
                                    Date date, const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                    int num_days, const std::string& divi_data_path,
                                    const std::string& confirmed_cases_path, const std::string& population_data_path)
{
    const size_t num_age_groups = scaling_factor_inf.size();
    std::vector<TimeSeries<double>> rki_data(
        region.size(), TimeSeries<double>::zero(num_days, (size_t)InfectionState::Count * num_age_groups));

    BOOST_OUTCOME_TRY(auto&& num_population, details::read_population_data(population_data_path, region));

    for (size_t t = 0; t < static_cast<size_t>(num_days); ++t) {
        auto offset_day = offset_date_by_days(date, t);

        if (offset_day > Date(2020, 4, 23)) {
            BOOST_OUTCOME_TRY(details::set_divi_data(models, divi_data_path, region, offset_day, scaling_factor_icu));
        }
        else {
            log_warning("No DIVI data available for this date");
            // TODO: print specific date
        }

        BOOST_OUTCOME_TRY(
            details::set_confirmed_cases_data(models, confirmed_cases_path, region, offset_day, scaling_factor_inf));

        // set population data
        for (size_t r = 0; r < region.size(); r++) {
            if (std::accumulate(num_population[r].begin(), num_population[r].end(), 0.0) > 0) {
                auto num_groups = models[r].parameters.get_num_groups();
                for (auto i = AgeGroup(0); i < num_groups; i++) {
                    models[r].populations.template set_difference_from_group_total<AgeGroup>(
                        {i, InfectionState::Susceptible}, num_population[r][size_t(i)]);
                }
            }
            else {
                log_warning("No population data available for region " + std::to_string(r) +
                            ". Population data has not been set.");
            }
        }
        for (size_t r = 0; r < region.size(); r++) {
            rki_data[r][t] = models[r].get_initial_values();
        }
    }
    auto num_groups = (int)(size_t)models[0].parameters.get_num_groups();
    BOOST_OUTCOME_TRY(save_result(rki_data, region, num_groups, path_join(dir, "Results_rki.h5")));

    auto rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<double>>>{rki_data});
    BOOST_OUTCOME_TRY(save_result({rki_data_sum[0][0]}, {0}, num_groups, path_join(dir, "Results_rki_sum.h5")));

    return success();
}
#else
template <class Model>
IOResult<void>
export_input_data_county_timeseries(std::vector<Model> models, const std::string& dir, std::vector<int> const& region,
                                    Date date, const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                    int num_days, const std::string& divi_data_path,
                                    const std::string& confirmed_cases_path, const std::string& population_data_path)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}
#endif // MEMILIO_HAS_HDF5

/**
 * @brief reads population data from population files for the whole country
 * @param model vector of model in which the data is set
 * @param date Date for which the data should be read
 * @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
 * @param scaling_factor_icu factor by which to scale the icu cases of divi data
 * @param dir directory of files
 */
template <class Model>
IOResult<void> read_input_data_germany(std::vector<Model>& model, Date date,
                                       const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                       const std::string& dir)
{
    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY(
            details::set_divi_data(model, path_join(dir, "germany_divi.json"), {0}, date, scaling_factor_icu));
    }
    else {
        log_warning("No DIVI data available for this date");
    }
    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(dir, "cases_all_age_ma7.json"), {0}, date,
                                                        scaling_factor_inf));
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(dir, "county_current_population.json"), {0}, false));
    return success();
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
IOResult<void> read_input_data_state(std::vector<Model>& model, Date date, std::vector<int>& state,
                                     const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                     const std::string& dir)
{
    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY(
            details::set_divi_data(model, path_join(dir, "state_divi.json"), state, date, scaling_factor_icu));
    }
    else {
        log_warning("No DIVI data available for this date");
    }

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(dir, "cases_all_state_age_ma7.json"), state,
                                                        date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(dir, "county_current_population.json"), state, false));
    return success();
}

/**
 * @brief reads population data from population files for the specefied county
 * @param model vector of model in which the data is set
 * @param date Date for which the data should be read
 * @param county vector of region keys of counties of interest
 * @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
 * @param scaling_factor_icu factor by which to scale the icu cases of divi data
 * @param dir directory of files
 * @param num_days [Default: 0] Number of days to be simulated; required to extrapolate real data
 * @param export_time_series [Default: false] If true, reads data for each day of simulation and writes it in the same directory as the input files.
 */
template <class Model>
IOResult<void> read_input_data_county(std::vector<Model>& model, Date date, const std::vector<int>& county,
                                      const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                      const std::string& dir, int num_days = 0, bool export_time_series = false)
{
    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY(details::set_divi_data(model, path_join(dir, "pydata/Germany", "county_divi_ma7.json"),
                                                 county, date, scaling_factor_icu));
    }
    else {
        log_warning("No DIVI data available for this date");
    }
    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(
        model, path_join(dir, "pydata/Germany", "cases_all_county_age_ma7.json"), county, date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(details::set_population_data(
        model, path_join(dir, "pydata/Germany", "county_current_population.json"), county, false));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        log_warning("Exporting time series of extrapolated real data. This may take some minutes. "
                    "For simulation runs over the same time period, deactivate it.");
        BOOST_OUTCOME_TRY(
            export_input_data_county_timeseries(model, dir, county, date, scaling_factor_inf, scaling_factor_icu,
                                                num_days, path_join(dir, "pydata/Germany", "county_divi_ma7.json"),
                                                path_join(dir, "pydata/Germany", "cases_all_county_age_ma7.json"),
                                                path_join(dir, "pydata/Germany", "county_current_population.json")));
    }
    return success();
}

/**
 * @brief reads population data from population files for the specified nodes
 * @param model vector of model in which the data is set
 * @param date Date for which the data should be read
 * @param county vector of region keys of interest
 * @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
 * @param scaling_factor_icu factor by which to scale the icu cases of divi data
 * @param dir directory of files
 * @param age_group_names strings specifying age group names
 */
template <class Model>
IOResult<void> read_input_data(std::vector<Model>& model, Date date, const std::vector<int>& node_ids,
                               const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                               const std::string& data_dir, int num_days = 0, bool export_time_series = false)
{
    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY(details::set_divi_data(model, path_join(data_dir, "critical_cases.json"), node_ids, date,
                                                 scaling_factor_icu));
    }
    else {
        log_warning("No DIVI data available for this date");
    }
    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(data_dir, "confirmed_cases.json"), node_ids,
                                                        date, scaling_factor_inf));
    bool single_age_group = scaling_factor_inf.size() == 1 ? true : false;
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(data_dir, "population_data.json"), node_ids, single_age_group));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        log_warning("Exporting time series of extrapolated real data. This may take some minutes. "
                    "For simulation runs over the same time period, deactivate it.");
        BOOST_OUTCOME_TRY(export_input_data_county_timeseries(
            model, data_dir, node_ids, date, scaling_factor_inf, scaling_factor_icu, num_days,
            path_join(data_dir, "critical_cases.json"), path_join(data_dir, "confirmed_cases.json"),
            path_join(data_dir, "population_data.json")));
    }
    return success();
}

} // namespace osecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // ODESECIR_PARAMETERS_IO_H
