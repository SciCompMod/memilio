/* 
* Copyright (C) 2020-2025 MEmilio
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
#include <cassert>

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secir/model.h"
#include "memilio/mobility/graph.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/parameters_io.h"
#include "memilio/io/result_io.h"
#include "memilio/math/math_utils.h"

namespace mio
{

namespace osecir
{

namespace details
{
/**
     * @brief reads populations data from RKI.
     * @param[in] rki_data Vector of ConfirmedCasesDataEntry%s.
     * @param[in] vregion Vector of keys of the region of interest.
     * @param[in] date Date at which the data is read.
     * @param[in, out] vnum_* Output vector for number of people in the corresponding compartement.
     * @param[in] vt_* vector Average time it takes to get from one compartement to another for each age group.
     * @param[in] vmu_* vector Probabilities to get from one compartement to another for each age group.
     * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
     */
IOResult<void> read_confirmed_cases_data(
    std::vector<ConfirmedCasesDataEntry>& rki_data, std::vector<int> const& vregion, Date date,
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
 * @brief Sets populations data from already read case data with multiple age groups into a Model with one age group.
 * @tparam FP Floating point data type, e.g., double.
 * @param[in, out] model Vector of models in which the data is set.
 * @param[in] case_data List of confirmed cases data entries.
 * @param[in] region Vector of keys of the region of interest.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 */
template <typename FP>
IOResult<void> set_confirmed_cases_data(std::vector<Model<FP>>& model, std::vector<ConfirmedCasesDataEntry>& case_data,
                                        const std::vector<int>& region, Date date,
                                        const std::vector<double>& scaling_factor_inf)
{
    const size_t num_age_groups = ConfirmedCasesDataEntry::age_group_names.size();
    // allow single scalar scaling that is broadcast to all age groups
    assert(scaling_factor_inf.size() == 1 || scaling_factor_inf.size() == num_age_groups);

    // Set scaling factors to match num age groups
    std::vector<double> scaling_factor_inf_full;
    if (scaling_factor_inf.size() == 1) {
        scaling_factor_inf_full.assign(num_age_groups, scaling_factor_inf[0]);
    }
    else {
        scaling_factor_inf_full = scaling_factor_inf;
    }

    std::vector<std::vector<int>> t_InfectedNoSymptoms{model.size()};
    std::vector<std::vector<int>> t_Exposed{model.size()};
    std::vector<std::vector<int>> t_InfectedSymptoms{model.size()};
    std::vector<std::vector<int>> t_InfectedSevere{model.size()};
    std::vector<std::vector<int>> t_InfectedCritical{model.size()};

    std::vector<std::vector<double>> mu_C_R{model.size()};
    std::vector<std::vector<double>> mu_I_H{model.size()};
    std::vector<std::vector<double>> mu_H_U{model.size()};
    std::vector<std::vector<double>> mu_U_D{model.size()};

    for (size_t node = 0; node < model.size(); ++node) {
        const size_t model_groups = (size_t)model[node].parameters.get_num_groups();
        assert(model_groups == 1 || model_groups == num_age_groups);
        for (size_t ag = 0; ag < num_age_groups; ag++) {
            // If the model has fewer groups than casedata entries available,
            // reuse group 0 parameters for all RKI age groups
            const size_t group = (model_groups == num_age_groups) ? ag : 0;

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

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf_full));

    for (size_t node = 0; node < model.size(); node++) {

        if (std::accumulate(num_InfectedSymptoms[node].begin(), num_InfectedSymptoms[node].end(), FP(0.0),
                            [](const FP& a, const FP& b) {
                                return evaluate_intermediate<FP>(a + b);
                            }) > 0.0) {
            size_t num_groups = (size_t)model[node].parameters.get_num_groups();
            if (num_groups == num_age_groups) {
                for (size_t i = 0; i < num_groups; i++) {
                    model[node].populations[{AgeGroup(i), InfectionState::Exposed}] = num_Exposed[node][i];
                    model[node].populations[{AgeGroup(i), InfectionState::InfectedNoSymptoms}] =
                        num_InfectedNoSymptoms[node][i];
                    model[node].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsConfirmed}] = 0;
                    model[node].populations[{AgeGroup(i), InfectionState::InfectedSymptoms}] =
                        num_InfectedSymptoms[node][i];
                    model[node].populations[{AgeGroup(i), InfectionState::InfectedSymptomsConfirmed}] = 0;
                    model[node].populations[{AgeGroup(i), InfectionState::InfectedSevere}] =
                        num_InfectedSevere[node][i];
                    // Only set the number of ICU patients here, if the date is not available in the data.
                    if (!is_divi_data_available(date)) {
                        model[node].populations[{AgeGroup(i), InfectionState::InfectedCritical}] = num_icu[node][i];
                    }
                    model[node].populations[{AgeGroup(i), InfectionState::Dead}]      = num_death[node][i];
                    model[node].populations[{AgeGroup(i), InfectionState::Recovered}] = num_rec[node][i];
                }
            }
            else {
                const auto sum_vec = [](const std::vector<FP>& v) {
                    return std::accumulate(v.begin(), v.end(), FP(0.0), [](const FP& a, const FP& b) {
                        return evaluate_intermediate<FP>(a + b);
                    });
                };
                const size_t i0                                                  = 0;
                model[node].populations[{AgeGroup(i0), InfectionState::Exposed}] = sum_vec(num_Exposed[node]);
                model[node].populations[{AgeGroup(i0), InfectionState::InfectedNoSymptoms}] =
                    sum_vec(num_InfectedNoSymptoms[node]);
                model[node].populations[{AgeGroup(i0), InfectionState::InfectedNoSymptomsConfirmed}] = 0;
                model[node].populations[{AgeGroup(i0), InfectionState::InfectedSymptoms}] =
                    sum_vec(num_InfectedSymptoms[node]);
                model[node].populations[{AgeGroup(i0), InfectionState::InfectedSymptomsConfirmed}] = 0;
                model[node].populations[{AgeGroup(i0), InfectionState::InfectedSevere}] =
                    sum_vec(num_InfectedSevere[node]);
                if (!is_divi_data_available(date)) {
                    model[node].populations[{AgeGroup(i0), InfectionState::InfectedCritical}] = sum_vec(num_icu[node]);
                }
                model[node].populations[{AgeGroup(i0), InfectionState::Dead}]      = sum_vec(num_death[node]);
                model[node].populations[{AgeGroup(i0), InfectionState::Recovered}] = sum_vec(num_rec[node]);
            }
        }
        else {
            log_warning("No infections reported on date {} for region {}. Population data has not been set.", date,
                        region[node]);
        }
    }
    return success();
}

/**
 * @brief Sets the infected population for a given model based on confirmed cases data. Here, we
 * read the case data from a file.
 * @tparam FP Floating point data type, e.g., double.
 * @param[in, out] Vector of models for which the confirmed cases data will be set.
 * @param[in] path Path to the confirmed cases data file.
 * @param[in] region Vector of keys of the region of interest.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 */
template <typename FP>
IOResult<void> set_confirmed_cases_data(std::vector<Model<FP>>& model, const std::string& path,
                                        std::vector<int> const& region, Date date,
                                        const std::vector<double>& scaling_factor_inf)
{
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path));
    BOOST_OUTCOME_TRY(set_confirmed_cases_data(model, case_data, region, date, scaling_factor_inf));
    return success();
}

/**
     * @brief Sets populations data from DIVI register into Model.
     * @tparam FP floating point data type, e.g., double.
     * @param[in, out] model Vector of models in which the data is set.
     * @param[in] path Path to DIVI file.
     * @param[in] vregion Vector of keys of the regions of interest.
     * @param[in] date Date for which the arrays are initialized.
     * @param[in] scaling_factor_icu factor by which to scale the icu cases of divi data.
     */
template <typename FP>
IOResult<void> set_divi_data(std::vector<Model<FP>>& model, const std::string& path, const std::vector<int>& vregion,
                             Date date, double scaling_factor_icu)
{
    // DIVI dataset will no longer be updated from CW29 2024 on.
    if (!is_divi_data_available(date)) {
        log_warning("No DIVI data available for date: {}. "
                    "ICU compartment will be set based on Case data.",
                    date);
        return success();
    }
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
* @brief Sets population data from census data which has been read into num_population.
* @tparam FP floating point data type, e.g., double.
* @param[in, out] model Vector of models in which the data is set. There should be one model per region.
* @param[in] num_population Vector of population data. The size should be the same as vregion and model.
* @param[in] vregion Vector of keys of the regions of interest.
*/
template <typename FP>
IOResult<void> set_population_data(std::vector<Model<FP>>& model, const std::vector<std::vector<FP>>& num_population,
                                   const std::vector<int>& vregion)
{
    assert(num_population.size() == vregion.size());
    assert(model.size() == vregion.size());
    for (size_t region = 0; region < vregion.size(); region++) {
        const auto model_groups = (size_t)model[region].parameters.get_num_groups();
        const auto data_groups  = num_population[region].size();
        assert(data_groups == model_groups || (model_groups == 1 && data_groups >= 1));

        if (data_groups == model_groups) {
            for (auto i = AgeGroup(0); i < model[region].parameters.get_num_groups(); i++) {
                model[region].populations.template set_difference_from_group_total<AgeGroup>(
                    {i, InfectionState::Susceptible}, num_population[region][(size_t)i]);
            }
        }
        else if (model_groups == 1 && data_groups >= 1) {
            const FP total = std::accumulate(num_population[region].begin(), num_population[region].end(), FP(0.0),
                                             [](const FP& a, const FP& b) {
                                                 return evaluate_intermediate<FP>(a + b);
                                             });
            model[region].populations.template set_difference_from_group_total<AgeGroup>(
                {AgeGroup(0), InfectionState::Susceptible}, total);
        }
    }
    return success();
}

/**
* @brief Sets population data from census data into a Model.
* @tparam FP Floating point data type, e.g., double.
* @param[in, out] model Vector of models in which the data is set.
* @param[in] path Path to RKI file containing population data.
* @param[in] vregion Vector of keys of the regions of interest.
*/
template <typename FP>
IOResult<void> set_population_data(std::vector<Model<FP>>& model, const std::string& path,
                                   const std::vector<int>& vregion)
{
    BOOST_OUTCOME_TRY(const auto&& num_population, read_population_data(path, vregion));
    BOOST_OUTCOME_TRY(set_population_data(model, num_population, vregion));
    return success();
}

} //namespace details

#ifdef MEMILIO_HAS_HDF5

/**
* @brief Uses the initialisation method, which uses the reported data to set the initial conditions for the model for a given day. 
* The initialisation is applied for a predefined number of days and finally saved in a timeseries for each region. In the end,
* we save the files "Results_rki.h5" and "Results_rki_sum.h5" in the results_dir.
* Results_rki.h5 contains a time series for each region and Results_rki_sum.h5 contains the sum of all regions.
* @param[in] models Vector of models in which the data is set. Copy is made to avoid changing the original model.
* @param[in] results_dir Path to result files.
* @param[in] region Vector of keys of the region of interest.
* @param[in] date Date for which the data should be read.
* @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
* @param[in] scaling_factor_icu Factor by which to scale the icu cases of divi data.
* @param[in] num_days Number of days to be simulated/initialized.
* @param[in] divi_data_path Path to DIVI file.
* @param[in] confirmed_cases_path Path to confirmed cases file.
* @param[in] population_data_path Path to population data file.
*/
template <class Model>
IOResult<void> export_input_data_county_timeseries(
    std::vector<Model> models, const std::string& results_dir, std::vector<int> const& region, Date date,
    const std::vector<double>& scaling_factor_inf, double scaling_factor_icu, int num_days,
    const std::string& divi_data_path, const std::string& confirmed_cases_path, const std::string& population_data_path)
{
    const auto num_age_groups = (size_t)models[0].parameters.get_num_groups();
    // allow scalar scaling factor as convenience for 1-group models
    assert(scaling_factor_inf.size() == 1 || scaling_factor_inf.size() == num_age_groups);
    assert(models.size() == region.size());
    std::vector<TimeSeries<double>> extrapolated_data(
        region.size(), TimeSeries<double>::zero(num_days + 1, (size_t)InfectionState::Count * num_age_groups));

    BOOST_OUTCOME_TRY(auto&& num_population, mio::read_population_data(population_data_path, region));
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(confirmed_cases_path));

    for (int t = 0; t <= num_days; ++t) {
        auto offset_day = offset_date_by_days(date, t);

        BOOST_OUTCOME_TRY(details::set_divi_data(models, divi_data_path, region, offset_day, scaling_factor_icu));
        BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(models, case_data, region, offset_day, scaling_factor_inf));
        BOOST_OUTCOME_TRY(details::set_population_data(models, num_population, region));
        for (size_t r = 0; r < region.size(); r++) {
            extrapolated_data[r][t] = models[r].get_initial_values();
        }
    }

    BOOST_OUTCOME_TRY(
        save_result(extrapolated_data, region, (int)num_age_groups, path_join(results_dir, "Results_rki.h5")));

    auto rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<double>>>{extrapolated_data});
    BOOST_OUTCOME_TRY(
        save_result({rki_data_sum[0][0]}, {0}, (int)num_age_groups, path_join(results_dir, "Results_rki_sum.h5")));

    return success();
}
#else
template <class Model>
IOResult<void> export_input_data_county_timeseries(
    std::vector<Model> models, const std::string& results_dir, std::vector<int> const& region, Date date,
    const std::vector<double>& scaling_factor_inf, double scaling_factor_icu, int num_days,
    const std::string& divi_data_path, const std::string& confirmed_cases_path, const std::string& population_data_path)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}
#endif // MEMILIO_HAS_HDF5

/**
 * @brief Reads population data from population files for the whole country.
 * @param[in, out] model Vector of model in which the data is set.
 * @param[in] date Date for which the data should be read.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 * @param[in] scaling_factor_icu Factor by which to scale the icu cases of divi data.
 * @param[in] pydata_dir Directory of files.
 */
template <class Model>
IOResult<void> read_input_data_germany(std::vector<Model>& model, Date date,
                                       const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                       const std::string& pydata_dir)
{
    BOOST_OUTCOME_TRY(
        details::set_divi_data(model, path_join(pydata_dir, "germany_divi.json"), {0}, date, scaling_factor_icu));
    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(pydata_dir, "cases_all_age_ma7.json"), {0},
                                                        date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(pydata_dir, "county_current_population.json"), {0}));
    return success();
}

/**
 * @brief Reads population data from population files for the specefied state.
 * @param[in, out] model Vector of model in which the data is set.
 * @param[in] date Date for which the data should be read.
 * @param[in] state Vector of region keys of states of interest.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 * @param[in] scaling_factor_icu Factor by which to scale the icu cases of divi data.
 * @param[in] pydata_dir Directory of files.
 */
template <class Model>
IOResult<void> read_input_data_state(std::vector<Model>& model, Date date, std::vector<int>& state,
                                     const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                     const std::string& pydata_dir)
{

    BOOST_OUTCOME_TRY(
        details::set_divi_data(model, path_join(pydata_dir, "state_divi.json"), state, date, scaling_factor_icu));
    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(pydata_dir, "cases_all_state_age_ma7.json"),
                                                        state, date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(pydata_dir, "county_current_population.json"), state));
    return success();
}

/**
 * @brief Reads population data from population files for the specefied county.
 * @param[in, out] model Vector of model in which the data is set.
 * @param[in] date Date for which the data should be read.
 * @param[in] county Vector of region keys of counties of interest.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 * @param[in] scaling_factor_icu Factor by which to scale the icu cases of divi data.
 * @param[in] pydata_dir Directory of files.
 * @param[in] num_days [Default: 0] Number of days to be simulated; required to extrapolate real data.
 * @param[in] export_time_series [Default: false] If true, reads data for each day of simulation and writes it in the same directory as the input files.
 */
template <class Model>
IOResult<void> read_input_data_county(std::vector<Model>& model, Date date, const std::vector<int>& county,
                                      const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                      const std::string& pydata_dir, int num_days = 0, bool export_time_series = false)
{
    BOOST_OUTCOME_TRY(
        details::set_divi_data(model, path_join(pydata_dir, "county_divi_ma7.json"), county, date, scaling_factor_icu));
    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(pydata_dir, "cases_all_county_age_ma7.json"),
                                                        county, date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(pydata_dir, "county_current_population.json"), county));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        log_warning("Exporting time series of extrapolated real data. This may take some minutes. "
                    "For simulation runs over the same time period, deactivate it.");
        BOOST_OUTCOME_TRY(export_input_data_county_timeseries(
            model, pydata_dir, county, date, scaling_factor_inf, scaling_factor_icu, num_days,
            path_join(pydata_dir, "county_divi_ma7.json"), path_join(pydata_dir, "cases_all_county_age_ma7.json"),
            path_join(pydata_dir, "county_current_population.json")));
    }
    return success();
}

/**
 * @brief reads population data from population files for the specified nodes
 * @param[in, out] model vector of model in which the data is set
 * @param[in] date Date for which the data should be read
 * @param[in] county vector of region keys of interest
 * @param[in] scaling_factor_inf factors by which to scale the confirmed cases of rki data
 * @param[in] scaling_factor_icu factor by which to scale the icu cases of divi data
 * @param[in] pydata_dir directory of files
 * @param[in] age_group_names strings specifying age group names
 */
template <class Model>
IOResult<void> read_input_data(std::vector<Model>& model, Date date, const std::vector<int>& node_ids,
                               const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                               const std::string& pydata_dir, int num_days = 0, bool export_time_series = false)
{
    BOOST_OUTCOME_TRY(details::set_divi_data(model, path_join(pydata_dir, "critical_cases.json"), node_ids, date,
                                             scaling_factor_icu));
    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(pydata_dir, "confirmed_cases.json"), node_ids,
                                                        date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(details::set_population_data(model, path_join(pydata_dir, "population_data.json"), node_ids));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        log_warning("Exporting time series of extrapolated real data. This may take some minutes. "
                    "For simulation runs over the same time period, deactivate it.");
        BOOST_OUTCOME_TRY(export_input_data_county_timeseries(
            model, pydata_dir, node_ids, date, scaling_factor_inf, scaling_factor_icu, num_days,
            path_join(pydata_dir, "critical_cases.json"), path_join(pydata_dir, "confirmed_cases.json"),
            path_join(pydata_dir, "population_data.json")));
    }
    return success();
}

/**
 * @brief Converts input data from one range of models to another with different type.
 * 
 * @tparam FP Floating point type.
 * @param[in] model_from VectorRange of Node%s each containing a Model with the input data.
 * @param[in,out] model_to VectorRange of Node%s each containing a Model to be initialized with data.
 */
template<class FP>
void convert_model_data_type(const mio::VectorRange<Node<Model<ScalarType>>>& model_from, const mio::VectorRange<Node<Model<FP>>>& model_to)
{
    assert(model_from.size() == model_to.size());
    assert((size_t)model_from[0].property.parameters.get_num_groups() == (size_t)model_to[0].property.parameters.get_num_groups());
    // Todo: add conversion of ParameterSet and then re-use code from other model parameters io 
    
    for (size_t region_idx = 0; region_idx < model_from.size(); ++region_idx) {
        // convert populations to mio::UncertainValue<FP>
        // needs 2 converts as mio::UncertainValue<ScalarType> -> mio::UncertainValue<FP> does not work
        model_to[region_idx].property.populations = model_from[region_idx].property.populations.template convert<FP>();
    }
}

} // namespace osecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // ODESECIR_PARAMETERS_IO_H
