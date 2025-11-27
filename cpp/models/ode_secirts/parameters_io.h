/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker, Wadim Koslow, Daniel Abele, Martin J. Kühn
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
#ifndef MIO_ODE_SECIRTS_PARAMETERS_IO_H
#define MIO_ODE_SECIRTS_PARAMETERS_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secirts/model.h"
#include "ode_secirts/analyze_result.h"
#include "memilio/geography/regions.h"
#include "memilio/math/eigen_util.h"
#include "memilio/math/math_utils.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/parameters_io.h"
#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/date.h"

namespace mio
{
namespace osecirts
{

namespace details
{

/**
 * @brief Reads populations data from rki data.
 * @param[in] case_data Vector of ConfirmedCasesDataEntry%s.
 * @param[in] region Key of the region of interest.
 * @param[in] date Date at which the data is read.
 * @param[in, out] num_* Output vector for number of people in the corresponding compartement.
 * @param[in] t_* vector Average time it takes to get from one compartement to another for each age group.
 * @param[in] mu_* vector Probabilities to get from one compartement to another for each age group.
 * @param[in] reduc_* vector Factors by which to adjust the times and transition rates based on immunity layer.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 * @param[in] layer Specifies the immunity layer: 0 (Naive), 1 (Partial Immunity), 2 (Improved Immunity).
 */
IOResult<void> compute_confirmed_cases_data(
    const std::vector<ConfirmedCasesDataEntry>& case_data, const int region, Date date,
    std::vector<ScalarType>& num_Exposed, std::vector<ScalarType>& num_InfectedNoSymptoms,
    std::vector<ScalarType>& num_InfectedSymptoms, std::vector<ScalarType>& num_InfectedSevere,
    std::vector<ScalarType>& num_icu, std::vector<ScalarType>& num_death, std::vector<ScalarType>& num_imm,
    const std::vector<int>& t_Exposed, const std::vector<int>& t_InfectedNoSymptoms,
    const std::vector<int>& t_InfectedSymptoms, const std::vector<int>& t_InfectedSevere,
    const std::vector<int>& t_InfectedCritical, const std::vector<int>& t_imm_interval_i,
    const std::vector<ScalarType>& mu_C_R, const std::vector<ScalarType>& mu_I_H, const std::vector<ScalarType>& mu_H_U,
    const std::vector<ScalarType>& reduc_t_Infected, const std::vector<ScalarType>& reduc_Exposed,
    const std::vector<ScalarType>& reduc_InfectedSymptoms, const std::vector<ScalarType>& reduc_icu_death,
    const std::vector<ScalarType>& scaling_factor_inf, const size_t layer);

/**
 * @brief Sets the confirmed cases data in the model considering different immunity layers.
 *
 * This function distributes confirmed case data across infection states for regions and age groups
 * in the model. It considers different levels of immunity (naive, partial, and improved).
 *
 * @param[in, out] model Model in which the data is set.
 * @param[in] case_data Vector of ConfirmedCasesDataEntry%s.
 * @param[in] region Key of the region of interest.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] immunity_population Vector containing the immunity distribution for naive, partial, and improved immunity layers.
 *
 * @return An IOResult indicating success or failure.
 */
IOResult<void> set_confirmed_cases_data(Model<ScalarType>& model, const std::vector<ConfirmedCasesDataEntry>& case_data,
                                        const int region, Date date, const std::vector<ScalarType>& scaling_factor_inf,
                                        const std::vector<std::vector<ScalarType>>& immunity_population);

/**
 * @brief Sets confirmed case data from a file into the given models.
 *
 * This function reads transformed rki data from the specified file and calls a subfunction that distributes the 
 * confirmed case data across different infection states for age groups in each model. 
 * It considers naive, partial, and improved immunity layers.
 *
 * @param[in, out] model VectorRange of Node%s each containing a Model for which the confirmed cases data will be set.
 * @param[in] path Path to the confirmed cases data file.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] immunity_population Vector containing the immunity distribution for naive, partial, and improved immunity layers.
 *
 * @return An IOResult indicating success or failure.
 */
IOResult<void> set_confirmed_cases_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path,
                                        Date date, const std::vector<ScalarType>& scaling_factor_inf,
                                        const std::vector<std::vector<ScalarType>>& immunity_population);

/**
 * @brief Sets the population data for the given model based on the provided population distribution and immunity levels.
 * @param[in, out] model Model in which the data is set.
 * @param[in] num_population Vector of population data for the region of interest.
 * @param[in] region Key of the region of interest.
 * @param[in] immunity_population A 2D vector where each row represents the immunity distribution for a specific region
 *                                 across different levels of immunity (e.g., naive, partial, improved immunity).
 *
 * @return An IOResult indicating success or failure.
 */
IOResult<void> set_population_data(Model<ScalarType>& model, const std::vector<ScalarType>& num_population,
                                   const int region, const std::vector<std::vector<ScalarType>>& immunity_population);

/**
 * @brief Reads population data from a file and sets it for the each given model.
 * @param[in, out] model VectorRange of Node%s containing a Model for which the confirmed cases data will be set.
 * @param[in] path Path to population data file.
 * @param[in] immunity_population A 2D vector specifying immunity for each age group and immunity layer.
 *
 * @return An IOResult indicating success or failure.
 */
IOResult<void> set_population_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path,
                                   const std::vector<std::vector<ScalarType>>& immunity_population);

/**
 * @brief Sets ICU data into a Model, distributed across age groups.
 * @param[in, out] model Model in which the data is set.
 * @param[in] num_icu ICU data for the region of interest.
 * @param[in] scaling_factor_icu Scaling factor for reported ICU cases.
 */
IOResult<void> set_divi_data(Model<ScalarType>& model, const ScalarType num_icu, ScalarType scaling_factor_icu);

/**
 * @brief Sets ICU data from DIVI data into the a vector of models, distributed across age groups.
 * 
 * This function reads DIVI data from a file, computes the number of individuals in critical condition (ICU)
 * for each region, and sets these values in the model. The ICU cases are distributed across age groups
 * using the transition probabilities from severe to critical.
 * 
 * @param[in, out] model VectorRange of Node%s each containing a Model for which the confirmed cases data will be set.
 * @param[in] path Path to transformed DIVI file.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_icu Scaling factor for reported ICU cases.
 */
IOResult<void> set_divi_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path, Date date,
                             ScalarType scaling_factor_icu);

/**
 * @brief Sets vaccination data into a Model using provided vaccination (partial, full, and booster) data.
 * @tparam FP Floating point type (default: double).
 * @param[in, out] model Model in which the data is set.
 * @param[in] vacc_data Vector of VaccinationDataEntry objects containing the vaccination data.
 * @param[in] date Start date for the simulation.
 * @param[in] num_days Number of days for which the simulation is run.
 * @param[in] max_date
 * @param[in] days_until_effective_n
 * @param[in] days_until_effective_pi
 * @param[in] days_until_effective_ii
 * 
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> set_vaccination_data(Model<FP>& model, const VaccinationDataEntry& vacc_data_entry, Date date,
                                    int num_days, Date max_date, const int days_until_effective_n,
                                    const int days_until_effective_pi, const int days_until_effective_ii)
{
    auto date_df = vacc_data_entry.date;
    AgeGroup age = vacc_data_entry.age_group;

    // get daily vaccinations for each layer
    for (size_t d = 0; d < (size_t)num_days + 1; ++d) {
        auto offset_first_date = offset_date_by_days(date, (int)d - days_until_effective_n);
        if (max_date >= offset_first_date) {
            if (date_df == offset_first_date) {
                model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] =
                    vacc_data_entry.num_vaccinations_partial;
            }
        }
        else {
            if (date_df == offset_first_date) {
                model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] = 0;
            }
        }

        auto offset_full_date = offset_date_by_days(date, (int)d - days_until_effective_pi);
        if (max_date >= offset_full_date) {
            if (date_df == offset_full_date) {
                model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] =
                    vacc_data_entry.num_vaccinations_completed;
            }
        }
        else {
            if (date_df == offset_first_date) {
                model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] = 0;
            }
        }

        auto offset_booster_date = offset_date_by_days(date, (int)d - days_until_effective_ii);
        if (max_date >= offset_booster_date) {
            if (date_df == offset_booster_date) {
                model.parameters.template get<DailyBoosterVaccinations<FP>>()[{age, SimulationDay(d)}] =
                    vacc_data_entry.num_vaccinations_refreshed_first +
                    vacc_data_entry.num_vaccinations_refreshed_additional;
            }
        }
        else {
            if (date_df == offset_first_date) {
                model.parameters.template get<DailyBoosterVaccinations<FP>>()[{age, SimulationDay(d)}] = 0;
            }
        }
    }
    return success();
}

/**
 * @brief Sets vaccination data for the given models using provided vaccination (partial, full, and booster) data.
 * 
 * @tparam FP Floating point type (default: double).
 * @param[in,out] model A vector of models for which vaccination data will be set.
 * @param[in] vacc_data A vector of VaccinationDataEntry objects containing the vaccination data.
 * @param[in] date The starting date for the simulation.
 * @param[in] num_days The number of days for which the simulation runs.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> set_vaccination_data(Model<FP>& model, const std::vector<VaccinationDataEntry>& vacc_data, Date date,
                                    int num_days)
{
    using std::floor;

    auto max_date_entry = std::max_element(vacc_data.begin(), vacc_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == vacc_data.end()) {
        return failure(StatusCode::InvalidFileFormat, "Vaccination data file is empty.");
    }
    auto max_date       = max_date_entry->date;
    auto min_date_entry = std::min_element(vacc_data.begin(), vacc_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    auto min_date       = min_date_entry->date;
    if (min_date > date || max_date < offset_date_by_days(date, num_days)) {
        log_warning("Vaccination data only available from {} to {}. "
                    "For days before and after, vaccinations will be set to 0.",
                    min_date, max_date);
    }

    auto days_until_effective_n = static_cast<int>(
        floor(static_cast<FP>(model.parameters.template get<DaysUntilEffectivePartialVaccination<FP>>()[AgeGroup(0)])));
    auto days_until_effective_pi = static_cast<int>(
        floor(static_cast<FP>(model.parameters.template get<DaysUntilEffectiveImprovedVaccination<FP>>()[AgeGroup(0)])));
    auto days_until_effective_ii = static_cast<int>(
        floor(static_cast<FP>(model.parameters.template get<DaysUntilEffectiveBoosterImmunity<FP>>()[AgeGroup(0)])));

    for (auto&& vacc_data_entry : vacc_data) {
        BOOST_OUTCOME_TRY(set_vaccination_data<FP>(model, vacc_data_entry, date, num_days, max_date,
                                                   days_until_effective_n, days_until_effective_pi,
                                                   days_until_effective_ii));
    }
    return success();
}

/**
 * @brief Sets vaccination data for the given models using vaccination data from a file.
 *
 * This function reads vaccination data from a specified file, and assigns daily vaccination numbers
 * (partial, full, and booster) to each region and age group in the models.
 *
 * @tparam FP Floating point type (default: double).
 * @param[in,out] model VectorRange of Node%s each containing a Model for which vaccination data will be set.
 * @param[in] path Path to vaccination data file.
 * @param[in] date The starting date for the simulation.
 * @param[in] num_days The number of days for which the simulation runs.
 * 
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> set_vaccination_data(const mio::VectorRange<Node<Model<FP>>>& model, const std::string& path, Date date,
                                    int num_days)
{
    // Set vaccination data to 0 for all models
    for (auto& m : model) {
        m.property.parameters.template get<DailyPartialVaccinations<FP>>().resize(SimulationDay(num_days + 1));
        m.property.parameters.template get<DailyFullVaccinations<FP>>().resize(SimulationDay(num_days + 1));
        m.property.parameters.template get<DailyBoosterVaccinations<FP>>().resize(SimulationDay(num_days + 1));
        for (auto d = SimulationDay(0); d < SimulationDay(num_days + 1); ++d) {
            for (auto a = AgeGroup(0); a < m.property.parameters.get_num_groups(); ++a) {
                m.property.parameters.template get<DailyPartialVaccinations<FP>>()[{a, d}] = 0.0;
                m.property.parameters.template get<DailyFullVaccinations<FP>>()[{a, d}]    = 0.0;
                m.property.parameters.template get<DailyBoosterVaccinations<FP>>()[{a, d}] = 0.0;
            }
        }
    }

    // Check if vaccination data is available for the given date range
    auto end_date = offset_date_by_days(date, num_days);
    if (!is_vaccination_data_available(date, end_date)) {
        log_warning("No vaccination data available in range from {} to {}. "
                    "Vaccination data will be set to 0.",
                    date, end_date);
        return success();
    }

    BOOST_OUTCOME_TRY(auto&& vacc_data, read_vaccination_data(path));

    // Sort vacc_data into regions and ignore once with no region associated
    std::vector<std::vector<VaccinationDataEntry>> vvacc_data{model.size()};
    for (auto&& vacc_data_entry : vacc_data) {
        auto it = std::find_if(model.begin(), model.end(), [&vacc_data_entry](auto&& m) {
            return m.id == 0 ||
                   (vacc_data_entry.county_id && vacc_data_entry.county_id == regions::de::CountyId(m.id)) ||
                   (vacc_data_entry.state_id && vacc_data_entry.state_id == regions::de::StateId(m.id)) ||
                   (vacc_data_entry.district_id && vacc_data_entry.district_id == regions::de::DistrictId(m.id));
        });
        if (it != model.end()) {
            auto region_idx = size_t(it - model.begin());
            vvacc_data[region_idx].push_back(vacc_data_entry);
        }
    }

    for (size_t region_idx = 0; region_idx < model.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_vaccination_data<FP>(model[region_idx].property, vvacc_data[region_idx], date, num_days));
    }

    return success();
}

} // namespace details

/**
 * @brief Reads compartments for geographic units at a specified date from data files.
 *
 * This function estimates all compartments from available data using the provided model parameters.
 *
 * @param[in,out] model VectorRange of Node%s each containing a Model to be initialized with data.
 * @param[in] date Date for which the data should be read.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] scaling_factor_icu Scaling factor for ICU cases.
 * @param[in] num_days Number of days to simulate.
 * @param[in] immunity_population Matrix containing immunity proportions for each age group and immunity layer.
 * @param[in] epidata_filenames Object containing the input data file paths.
 *
 * @return An IOResult indicating success or failure.
 */
IOResult<void> read_input_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, Date date,
                               const std::vector<ScalarType>& scaling_factor_inf, ScalarType scaling_factor_icu,
                               int num_days, const std::vector<std::vector<ScalarType>>& immunity_population,
                               const mio::regions::de::EpidataFilenames& epidata_filenames);

/**
 * @brief Converts input data from one range of models to another with different type.
 * 
 * @tparam FP Floating point type (default: double).
 * @param[in] model_from VectorRange of Node%s each containing a Model with the input data.
 * @param[in,out] model_to VectorRange of Node%s each containing a Model to be initialized with data.
 * @param[in] date Date for which the data should be read.
 * @param[in] num_days Number of days to simulate.
 * @param[in] epidata_filenames Object containing the input data file paths.
 *
 * @return An IOResult indicating success or failure.
 */
template <class FP>
IOResult<void> convert_input_data_type(const mio::VectorRange<Node<Model<ScalarType>>>& model_from,
                                       const mio::VectorRange<Node<Model<FP>>>& model_to, Date date, int num_days,
                                       const mio::regions::de::EpidataFilenames& epidata_filenames)
{
    assert(model_from.size() == model_to.size());
    assert((size_t)model_from[0].property.parameters.get_num_groups() ==
           (size_t)model_to[0].property.parameters.get_num_groups());
    // Todo: add conversion of ParameterSet and then re-use code from all model parameters io
    // For now call set_vacination_data with FP to set correct parameters
    BOOST_OUTCOME_TRY(
        details::set_vaccination_data<FP>(model_to, epidata_filenames.vaccination_data_path, date, num_days));

    for (size_t region_idx = 0; region_idx < model_from.size(); ++region_idx) {
        // convert populations to mio::UncertainValue<FP>
        // needs 2 converts as mio::UncertainValue<ScalarType> -> mio::UncertainValue<FP> does not work
        model_to[region_idx].property.populations = model_from[region_idx]
                                                        .property.populations.template convert<FP>();
    }
    return success();
}

#ifdef MEMILIO_HAS_HDF5

/**
 * @brief Uses the initialisation method, which uses the reported data to set the initial conditions for the model for a given day. 
 * 
 * The initialisation is applied for a predefined number of days and finally saved in a timeseries for each region. In the end,
 * we save the files "Results_rki.h5" and "Results_rki_sum.h5" in the results_dir.
 * Results_rki.h5 contains a time series for each region and Results_rki_sum.h5 contains the sum of all regions.
 * 
 * @param[in] model VectorRange of Node%s each containing a Model in which the data is set.
 * @param[in] results_dir Path to result files.
 * @param[in] date Date for which the data should be read.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] scaling_factor_icu Scaling factor for ICU cases.
 * @param[in] num_days Number of days to be simulated/initialized.
 * @param[in] immunity_population Matrix containing immunity proportions for each age group and immunity layer.
 * @param[in] epidata_filenames Object containing the input data file paths.
 * 
 * @return An IOResult indicating success or failure.
 */
IOResult<void> export_input_data_timeseries(const mio::VectorRange<Node<Model<ScalarType>>> model,
                                            const std::string& results_dir, Date date,
                                            const std::vector<ScalarType>& scaling_factor_inf,
                                            const ScalarType scaling_factor_icu, const int num_days,
                                            const std::vector<std::vector<ScalarType>>& immunity_population,
                                            const mio::regions::de::EpidataFilenames& epidata_filenames);
#else
IOResult<void> export_input_data_timeseries(const mio::VectorRange<Node<Model<ScalarType>>>, const std::string&, Date,
                                            const std::vector<int>&, const std::vector<ScalarType>&, const ScalarType,
                                            const int, const std::vector<std::vector<ScalarType>>,
                                            const mio::regions::de::EpidataFilenames&);

#endif //MEMILIO_HAS_HDF5

} // namespace osecirts

} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // MIO_ODE_SECIRTS_PARAMETERS_IO_H
