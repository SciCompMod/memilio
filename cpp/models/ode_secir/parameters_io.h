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
#include "memilio/geography/regions.h"
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
     * @param[in] case_data Vector of ConfirmedCasesDataEntry%s.
     * @param[in] region Vector of keys of the region of interest.
     * @param[in] date Date at which the data is read.
     * @param[in, out] num_* Output vector for number of people in the corresponding compartement.
     * @param[in] t_* vector Average time it takes to get from one compartement to another for each age group.
     * @param[in] mu_* vector Probabilities to get from one compartement to another for each age group.
     * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
     */
IOResult<void> compute_confirmed_cases_data(
    std::vector<ConfirmedCasesDataEntry>& case_data, const int region, Date date,
    std::vector<double>& num_Exposed, std::vector<double>& num_InfectedNoSymptoms,
    std::vector<double>& num_InfectedSymptoms, std::vector<double>& num_InfectedSevere,
    std::vector<double>& num_icu, std::vector<double>& num_death,
    std::vector<double>& num_rec, const std::vector<int>& t_Exposed,
    const std::vector<int>& t_InfectedNoSymptoms,
    const std::vector<int>& t_InfectedSymptoms, const std::vector<int>& t_InfectedSevere,
    const std::vector<int>& t_InfectedCritical, const std::vector<double>& mu_C_R,
    const std::vector<double>& mu_I_H, const std::vector<double>& mu_H_U,
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
IOResult<void> set_confirmed_cases_data(Model<FP>& model, std::vector<ConfirmedCasesDataEntry>& case_data,
                                        const int region, Date date, const std::vector<double>& scaling_factor_inf)
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

    std::vector<int> t_InfectedNoSymptoms;
    std::vector<int> t_Exposed;
    std::vector<int> t_InfectedSymptoms;
    std::vector<int> t_InfectedSevere;
    std::vector<int> t_InfectedCritical;

    std::vector<double> mu_C_R;
    std::vector<double> mu_I_H;
    std::vector<double> mu_H_U;
    std::vector<double> mu_U_D;

    std::vector<double> num_InfectedSymptoms(num_age_groups, 0.0);
    std::vector<double> num_death(num_age_groups, 0.0);
    std::vector<double> num_rec(num_age_groups, 0.0);
    std::vector<double> num_Exposed(num_age_groups, 0.0);
    std::vector<double> num_InfectedNoSymptoms(num_age_groups, 0.0);
    std::vector<double> num_InfectedSevere(num_age_groups, 0.0);
    std::vector<double> num_icu(num_age_groups, 0.0);

    const size_t model_groups = (size_t)model.parameters.get_num_groups();
    assert(model_groups == 1 || model_groups == num_age_groups);
    for (size_t ag = 0; ag < num_age_groups; ag++) {
        // If the model has fewer groups than casedata entries available,
        // reuse group 0 parameters for all RKI age groups
        const size_t group = (model_groups == num_age_groups) ? ag : 0;

        t_Exposed.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeExposed<FP>>()[(AgeGroup)group])));
        t_InfectedNoSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedNoSymptoms<FP>>()[(AgeGroup)group])));
        t_InfectedSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSymptoms<FP>>()[(AgeGroup)group])));
        t_InfectedSevere.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSevere<FP>>()[(AgeGroup)group])));
        t_InfectedCritical.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedCritical<FP>>()[(AgeGroup)group])));

        mu_C_R.push_back(
            model.parameters.template get<RecoveredPerInfectedNoSymptoms<FP>>()[(AgeGroup)group]);
        mu_I_H.push_back(
            model.parameters.template get<SeverePerInfectedSymptoms<FP>>()[(AgeGroup)group]);
        mu_H_U.push_back(
            model.parameters.template get<CriticalPerSevere<FP>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf_full));

    if (model_groups == num_age_groups) {
        for (size_t i = 0; i < model_groups; i++) {
            model.populations[{AgeGroup(i), InfectionState::Exposed}] = num_Exposed[i];
            model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptoms}] =
                num_InfectedNoSymptoms[i];
            model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsConfirmed}] = 0;
            model.populations[{AgeGroup(i), InfectionState::InfectedSymptoms}] =
                num_InfectedSymptoms[i];
            model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsConfirmed}] = 0;
            model.populations[{AgeGroup(i), InfectionState::InfectedSevere}] =
                num_InfectedSevere[i];
            // Only set the number of ICU patients here, if the date is not available in the data.
            if (!is_divi_data_available(date)) {
                model.populations[{AgeGroup(i), InfectionState::InfectedCritical}] = num_icu[i];
            }
            model.populations[{AgeGroup(i), InfectionState::Recovered}] = num_rec[i];
            model.populations[{AgeGroup(i), InfectionState::Dead}]      = num_death[i];
        }
    }
    else {
        const auto sum_vec = [](const std::vector<FP>& v) {
            return std::accumulate(v.begin(), v.end(), FP(0.0), [](const FP& a, const FP& b) {
                return evaluate_intermediate<FP>(a + b);
            });
        };
        model.populations[{AgeGroup(0), InfectionState::Exposed}] = sum_vec(num_Exposed);
        model.populations[{AgeGroup(0), InfectionState::InfectedNoSymptoms}] =
            sum_vec(num_InfectedNoSymptoms);
        model.populations[{AgeGroup(0), InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model.populations[{AgeGroup(0), InfectionState::InfectedSymptoms}] =
            sum_vec(num_InfectedSymptoms);
        model.populations[{AgeGroup(0), InfectionState::InfectedSymptomsConfirmed}] = 0;
        model.populations[{AgeGroup(0), InfectionState::InfectedSevere}] =
            sum_vec(num_InfectedSevere);
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(0), InfectionState::InfectedCritical}] = sum_vec(num_icu);
        }
        model.populations[{AgeGroup(0), InfectionState::Dead}]      = sum_vec(num_death);
        model.populations[{AgeGroup(0), InfectionState::Recovered}] = sum_vec(num_rec);
    }
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), double(0.0),
                        [](const double& a, const double& b) {
                            return evaluate_intermediate<double>(a + b);
                        }) == 0.0) {
        log_warning(
            "No infections for unvaccinated reported on date {} for region {}. Population data has not been set.",
            date, region);
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
IOResult<void> set_confirmed_cases_data(mio::VectorRange<Node<Model<FP>>>& model, const std::string& path,
                                        std::vector<int> const& region, Date date,
                                        const std::vector<double>& scaling_factor_inf)
{
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path));

    // sort case_data into regions and ignore once with no region associated
    std::vector<std::vector<ConfirmedCasesDataEntry>> vcase_data{model.size()};
    for (auto&& entry : case_data) {
        auto it = std::find_if(vregion.begin(), vregion.end(), [&entry](auto r) {
            return r == 0 || get_region_id(entry) == r;
        });
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());
            vcase_data[region_idx].pushback(entry);
        }
    }

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(
            set_confirmed_cases_data(model[region_idx].property, vcase_data[region_idx], vregion[region_idx], date, scaling_factor_inf));
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
IOResult<void> set_population_data(Model<FP>& model, const std::vector<double>& num_population,
                                   const int region)
{
    if (std::accumulate(num_population.begin(), num_population.end(), double(0.0),
                        [](const double& a, const double& b) {
                            return evaluate_intermediate<double>(a + b);
                        }) <= 0)
    {    
        log_warning("No population data available for region " + std::to_string(region) +
                    ". Population data has not been set.");
        return success();
    }

    auto num_groups = model.parameters.get_num_groups();
    auto data_groups  = num_population.size();

    if (data_groups == model_groups) {
        for (auto i = AgeGroup(0); i < num_groups; i++) {
            model.populations.template set_difference_from_group_total<AgeGroup>(
                {i, InfectionState::Susceptible}, num_population[(size_t)i]);
        }
    }
    else if (model_groups == 1 && data_groups >= 1) {
        const FP total = std::accumulate(num_population.begin(), num_population.end(), FP(0.0),
                                            [](const FP& a, const FP& b) {
                                                return evaluate_intermediate<FP>(a + b);
                                            });
        model.populations.template set_difference_from_group_total<AgeGroup>(
            {AgeGroup(0), InfectionState::Susceptible}, total);
    }

    for (auto i = AgeGroup(0); i < AgeGroup(6); i++) {
        for (auto j = Index<InfectionState>(0); j < InfectionState::Count; ++j) {
            if (model.populations[{i, j}] < 0) {
                log_warning("Compartment at age group {}, infection state {}, is negative: {}", size_t(i),
                            size_t(j), model.populations[{i, j}]);
            }
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
IOResult<void> set_population_data(mio::VectorRange<Node<Model<FP>>>& model, const std::string& path,
                                   const std::vector<int>& vregion)
{
    BOOST_OUTCOME_TRY(const auto&& num_population, read_population_data(path, vregion));

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_population_data(model[region_idx].property, num_population[region_idx], vregion[region_idx]));
    }
    return success();
}

/**
 * @brief Sets ICU data from DIVI data into the a vector of models, distributed across age groups.
 *
 * This function reads DIVI data from a file, computes the number of individuals in critical condition (ICU)
 * for each region, and sets these values in the model. The ICU cases are distributed across age groups
 * using the transition probabilities from severe to critical.
 * @tparam FP Floating point type (default: double).
 *
 * @param[in,out] model Vector of models, each representing a region, where the ICU population is updated.
 * @param[in] num_icu icu data
 * @param[in] scaling_factor_icu Scaling factor for reported ICU cases.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> set_divi_data(Model<FP>& model, const double num_icu, double scaling_factor_icu)
{
    FP sum_mu_I_U = 0;
    std::vector<FP> mu_I_U;
    auto num_groups = model.parameters.get_num_groups();
    for (auto i = AgeGroup(0); i < num_groups; i++) {
        sum_mu_I_U += model.parameters.template get<CriticalPerSevere<FP>>()[i] *
                                model.parameters.template get<SeverePerInfectedSymptoms<FP>>()[i];
        mu_I_U.push_back(model.parameters.template get<CriticalPerSevere<FP>>()[i] *
                                    model.parameters.template get<SeverePerInfectedSymptoms<FP>>()[i]);
    }

    for (auto i = AgeGroup(0); i < num_groups; i++) {
        model.populations[{i, InfectionState::InfectedCriticalNaive}] =
            scaling_factor_icu * num_icu * mu_I_U[(size_t)i] / sum_mu_I_U;
    }

    return success();
}

/**
 * @brief sets populations data from DIVI register into Model
 * @param[in, out] model vector of objects in which the data is set
 * @param[in] path Path to transformed DIVI file
 * @param[in] vregion vector of keys of the regions of interest
 * @param[in] date Date for which the arrays are initialized
 * @param[in] scaling_factor_icu factor by which to scale the icu cases of divi data
 */
template <class FP>
IOResult<void> set_divi_data(mio::VectorRange<Node<Model<FP>>>& model, const std::string& path, const std::vector<int>& vregion,
                             Date date, double scaling_factor_icu)
{
    // DIVI dataset will no longer be updated from CW29 2024 on.
    if (!is_divi_data_available(date)) {
        log_warning("No DIVI data available for date: {}. "
                    "ICU compartment will be set based on Case data.",
                    date);
        return success();
    }
    BOOST_OUTCOME_TRY(auto&& num_icu, read_divi_data(path, vregion, date));

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_divi_data(model[region_idx].property, num_icu[region_idx], vregion[region_idx], scaling_factor_icu));
    }

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
* @param[in] date Date for which the data should be read.
* @param[in] node_ids Vector of keys of the node_ids of interest.
* @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
* @param[in] scaling_factor_icu Factor by which to scale the icu cases of divi data.
* @param[in] num_days Number of days to be simulated/initialized.
* @param[in] pydata_dir Directory that contains the data files.
*/
template <class FP>
IOResult<void> export_input_data_timeseries(
    mio::VectorRange<Node<Model<FP>>> models, const std::string& results_dir, Date date, const std::vector<int>& node_ids,
    const std::vector<double>& scaling_factor_inf, double scaling_factor_icu, int num_days,
    const mio::regions::de::EpidataFilenames& epidata_filenames)
{
    const auto num_age_groups = (size_t)models[0].property.parameters.get_num_groups();
    // allow scalar scaling factor as convenience for 1-group models
    assert(scaling_factor_inf.size() == 1 || scaling_factor_inf.size() == num_age_groups);
    assert(models.size() == node_ids.size());
    std::vector<TimeSeries<double>> extrapolated_data(
        node_ids.size(), TimeSeries<double>::zero(num_days + 1, (size_t)InfectionState::Count * num_age_groups));

    for (int t = 0; t <= num_days; ++t) {
        auto offset_day = offset_date_by_days(date, t);

        BOOST_OUTCOME_TRY(read_input_data(model, date, county, scaling_factor_inf, scaling_factor_icu,
                                      num_days, epidata_filenames));

        for (size_t r = 0; r < node_ids.size(); r++) {
            extrapolated_data[r][t] = models[r].property.get_initial_values();
        }
    }

    BOOST_OUTCOME_TRY(save_result(extrapolated_data, node_ids, static_cast<int>(num_age_groups),
                                  path_join(results_dir, "Results_rki.h5")));

    auto extrapolated_rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<double>>>{extrapolated_data});
    BOOST_OUTCOME_TRY(save_result({extrapolated_rki_data_sum[0][0]}, {0}, static_cast<int>(num_age_groups),
                                  path_join(results_dir, "Results_rki_sum.h5")));

    return success();
}
#else
template <class FP>
IOResult<void> export_input_data_county_timeseries(mio::VectorRange<Node<Model<FP>>>, const std::string&, Date, const std::vector<int>&,
                                                   const std::vector<double>&, const double, const int,
                                                   const mio::regions::de::EpidataFilenames&)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}
#endif // MEMILIO_HAS_HDF5

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
template <typename FP>
IOResult<void> read_input_data(mio::VectorRange<Node<Model<FP>>>& model, Date date, const std::vector<int>& node_ids,
                               const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                               const mio::regions::de::EpidataFilenames& epidata_filenames)
{
    BOOST_OUTCOME_TRY(details::set_divi_data(model, epidata_filenames.divi_data_path, node_ids, date,
                                             scaling_factor_icu));

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, epidata_filenames.case_data_path, node_ids,
                                                        date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(details::set_population_data(model, epidata_filenames.population_data_path, node_ids));
    return success();
}

} // namespace osecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // ODESECIR_PARAMETERS_IO_H
