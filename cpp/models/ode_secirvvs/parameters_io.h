/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kühn
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
#ifndef MIO_ODE_SECIRVVS_PARAMETERS_IO_H
#define MIO_ODE_SECIRVVS_PARAMETERS_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secirvvs/model.h"
#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/parameters_io.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/logging.h"
#include "memilio/math/math_utils.h"

namespace mio
{
namespace osecirvvs
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
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 */
IOResult<void>
compute_confirmed_cases_data(const std::vector<ConfirmedCasesDataEntry>& case_data, const int region, Date date,
                             std::vector<ScalarType>& num_Exposed, std::vector<ScalarType>& num_InfectedNoSymptoms,
                             std::vector<ScalarType>& num_InfectedSymptoms, std::vector<ScalarType>& num_InfectedSevere,
                             std::vector<ScalarType>& num_icu, std::vector<ScalarType>& num_death,
                             std::vector<ScalarType>& num_rec, const std::vector<int>& t_Exposed,
                             const std::vector<int>& t_InfectedNoSymptoms, const std::vector<int>& t_InfectedSymptoms,
                             const std::vector<int>& t_InfectedSevere, const std::vector<int>& t_InfectedCritical,
                             const std::vector<ScalarType>& mu_C_R, const std::vector<ScalarType>& mu_I_H,
                             const std::vector<ScalarType>& mu_H_U, const std::vector<ScalarType>& scaling_factor_inf);

/**
 * @brief Reads confirmed cases data and translates data of day t0-delay to recovered compartment.
 * @param[in] case_data Vector of ConfirmedCasesDataEntry%s.
 * @param[in] region Key of the region of interest.
 * @param[in] date Date for which the arrays are initialized.
 * @param[in] delay Number of days in the past that are used to set recovered compartment.
 * @return Vector for number of people in the compartment recovered.
 */
IOResult<std::vector<ScalarType>>
compute_confirmed_cases_data_fix_recovered(const std::vector<ConfirmedCasesDataEntry>& case_data, const int region,
                                           Date date, ScalarType delay = 14.);

/**
 * @brief Sets populations data from already read case data with multiple age groups into a Model.
 * @param[in, out] model Model in which the data is set.
 * @param[in] case_data Vector of ConfirmedCasesDataEntry%s.
 * @param[in] region Key of the region of interest.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 * @param[in] set_death If true, set the number of deaths.
 */
IOResult<void> set_confirmed_cases_data(Model<ScalarType>& model, const std::vector<ConfirmedCasesDataEntry>& case_data,
                                        const int region, Date date, const std::vector<ScalarType>& scaling_factor_inf,
                                        bool set_death = false);

/**
 * @brief Sets the infected populations for the given models based on confirmed cases data. 
 * Reads the case data from a file and then calls a subfunction that sets the infected population 
 * for each model.
 * @param[in, out] model VectorRange of Node%s containing a Model for which the confirmed cases data will be set.
 * @param[in] path Path to the confirmed cases data file.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 * @param set_death[in] If true, set the number of deaths.
 */
IOResult<void> set_confirmed_cases_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path,
                                        Date date, const std::vector<ScalarType>& scaling_factor_inf,
                                        bool set_death = false);

/**
 * @brief Sets the population data for the given model based on the provided population distribution.
 * @param[in, out] model Model in which the data is set.
 * @param[in] num_population Vector of population data for the region of interest.
 * @param[in] region Key of the region of interest.
 * @param[in] date Date for which the arrays are initialized.
 */
IOResult<void> set_population_data(Model<ScalarType>& model, const std::vector<ScalarType>& num_population,
                                   const int region, const std::vector<ScalarType>& num_rec);

/**
 * @brief Reads population data from a file and sets it for the each given model.
 * @param[in, out] model VectorRange of Node%s containing a Model for which the confirmed cases data will be set.
 * @param[in] path Path to population data file.
 * @param[in] path_rki Path to rki cases data file.
 * @param[in] date Date for which the arrays are initialized.
 */
IOResult<void> set_population_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path,
                                   const std::string& path_rki, Date date);

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
 * @brief Sets vaccination data into a Model.
 * @tparam FP Floating point type (default: double).
 * @param[in, out] model Model in which the data is set.
 * @param[in] vacc_data Vector of VaccinationDataEntry objects containing the vaccination data.
 * @param[in] date Start date for the simulation.
 * @param[in] num_days Number of days for which the simulation is run.
 * 
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> set_vaccination_data(Model<FP>& model, const std::vector<VaccinationDataEntry>& vacc_data, Date date,
                                    int num_days)
{
    auto max_date_entry = std::max_element(vacc_data.begin(), vacc_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == vacc_data.end()) {
        return failure(StatusCode::InvalidFileFormat, "Vaccination data file is empty.");
    }
    auto max_date = max_date_entry->date;

    auto min_date_entry = std::min_element(vacc_data.begin(), vacc_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    auto min_date       = min_date_entry->date;
    if (min_date > date || max_date < offset_date_by_days(date, num_days)) {
        log_warning("Vaccination data only available from {} to {}. "
                    "For days before, vaccination data will be set to 0. For days after, "
                    "vaccination data will be set to the last available date.",
                    min_date, max_date);
    }

    // type conversion from UncertainValue -> FP -> int
    auto days_until_effective1 = static_cast<int>(
        static_cast<FP>(model.parameters.template get<DaysUntilEffectivePartialImmunity<FP>>()[AgeGroup(0)]));
    auto days_until_effective2 = static_cast<int>(
        static_cast<FP>(model.parameters.template get<DaysUntilEffectiveImprovedImmunity<FP>>()[AgeGroup(0)]));
    auto vaccination_distance =
        static_cast<int>(static_cast<FP>(model.parameters.template get<VaccinationGap<FP>>()[AgeGroup(0)]));

    for (auto&& vacc_data_entry : vacc_data) {
        auto date_df = vacc_data_entry.date;
        AgeGroup age = vacc_data_entry.age_group;

        for (size_t d = 0; d < (size_t)num_days + 1; ++d) {
            int days_plus;
            // In the following, second dose means previous 'full immunization', now 'Grundimmunisierung'.
            // ---
            // date: start_date of the simulation (Input from IO call read_input_data_county_vaccmodel())
            // d: day of simulation, counted from 0 to num_days (for which we need (approximated) vaccination numbers)
            // root[i]["Vacc_completed"]: accumulated number of total second doses up to day date_df;
            //                               taken from input dataframe, single value, per county and age group
            // ----
            // An averaged distance between first and second doses (vaccination_distance) is assumed in the following
            // and the first doses are computed based on the second doses given 'vaccination_distance' days later.
            // ----
            // a person whose second dose is reported at start_date + simulation_day - days_until_effective1 + vaccination_distance
            // had the first dose on start_date + simulation_day - days_until_effective1. Furthermore, he/she has the full protection
            // of the first dose at day X = start_date + simulation_day
            // Storing its value in get<DailyPartialVaccinations>() will eventually (in the simulation)
            // transfer the difference (between get<DailyPartialVaccinations>() at d and d-1) of
            // N susceptible individuals to 'Susceptible Partially Vaccinated' state at day d; see secir_vaccinated.h
            auto offset_first_date = offset_date_by_days(date, (int)d - days_until_effective1 + vaccination_distance);
            if (max_date >= offset_first_date) {
                // Option 1: considered offset_first_date is available in input data frame
                if (date_df == offset_first_date) {
                    model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] =
                        vacc_data_entry.num_vaccinations_completed;
                }
            }
            else { // offset_first_date > max_date
                // Option 2: considered offset_first_date is NOT available in input data frame
                // Here, a constant number of first and second doses is assumed, i.e.,
                // the the number of vaccinationes at day d (N days after max_date) will be:
                // total number of vaccinations up to day max_date + N * number of vaccinations ON max_date
                // (where the latter is computed as the difference between the total number at max_date and max_date-1)
                days_plus = get_offset_in_days(offset_first_date, max_date);
                if (date_df == offset_date_by_days(max_date, -1)) {
                    model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] -=
                        days_plus * vacc_data_entry.num_vaccinations_completed;
                }
                else if (date_df == max_date) {
                    model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] +=
                        (days_plus + 1) * vacc_data_entry.num_vaccinations_completed;
                }
            }

            // a person whose second dose is reported at start_date + simulation_day - days_until_effective2
            // has the full protection of the second dose at day X = start_date + simulation_day
            // Storing its value in get<DailyFullVaccinations>() will eventually (in the simulation)
            // transfer the difference (between get<DailyFullVaccinations>() at d and d-1) of
            // N susceptible, partially vaccinated individuals to 'SusceptibleImprovedImmunity' state at day d; see secir_vaccinated.h
            auto offset_full_date = offset_date_by_days(date, (int)d - days_until_effective2);
            if (max_date >= offset_full_date) {
                // Option 1: considered offset_full_date is available in input data frame
                if (date_df == offset_full_date) {
                    model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] =
                        vacc_data_entry.num_vaccinations_completed;
                }
            }
            else { // offset_full_date > max_full_date
                // Option 2: considered offset_full_date is NOT available in input data frame
                days_plus = get_offset_in_days(offset_full_date, max_date);
                if (date_df == offset_date_by_days(max_date, -1)) {
                    model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] -=
                        days_plus * vacc_data_entry.num_vaccinations_completed;
                }
                else if (date_df == max_date) {
                    model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] +=
                        (days_plus + 1) * vacc_data_entry.num_vaccinations_completed;
                }
            }
        }
    }
    return success();
}

/**
 * @brief Reads vaccination data from a file and sets it for each model.
 * @tparam FP Floating point type (default: double).
 * @param[in, out] model VectorRange of Node%s each containing a Model in which the vaccination data is set.
 * @param[in] path Path to vaccination data file.
 * @param[in] date Start date for the simulation.
 * @param[in] num_days Number of days for which the simulation is run.
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
        for (auto d = SimulationDay(0); d < SimulationDay(num_days + 1); ++d) {
            for (auto a = AgeGroup(0); a < m.property.parameters.get_num_groups(); ++a) {
                m.property.parameters.template get<DailyPartialVaccinations<FP>>()[{a, d}] = 0.0;
                m.property.parameters.template get<DailyFullVaccinations<FP>>()[{a, d}]    = 0.0;
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
 * @param[in] epidata_filenames Object containing the input data file paths.
 *
 * @return An IOResult indicating success or failure.
 */
IOResult<void> read_input_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, Date date,
                               const std::vector<ScalarType>& scaling_factor_inf, ScalarType scaling_factor_icu,
                               int num_days, const mio::regions::de::EpidataFilenames& epidata_filenames);

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
        model_to[region_idx].property.populations = model_to[region_idx]
                                                        .property.populations.template convert<ScalarType>()
                                                        .template convert<mio::UncertainValue<FP>>();
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
 * @param[in] epidata_filenames Object containing the input data file paths.
 * 
 * @return An IOResult indicating success or failure.
 */
IOResult<void> export_input_data_timeseries(const mio::VectorRange<Node<Model<ScalarType>>> model,
                                            const std::string& results_dir, Date date,
                                            const std::vector<ScalarType>& scaling_factor_inf,
                                            const ScalarType scaling_factor_icu, const int num_days,
                                            const mio::regions::de::EpidataFilenames& epidata_filenames);
#else
IOResult<void> export_input_data_county_timeseries(const mio::VectorRange<Node<Model<ScalarType>>>, const std::string&,
                                                   Date, const std::vector<int>&, const std::vector<ScalarType>&,
                                                   const ScalarType, const int,
                                                   const mio::regions::de::EpidataFilenames&);

#endif //MEMILIO_HAS_HDF5

} // namespace osecirvvs

} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // MIO_ODE_SECIRVVS_PARAMETERS_IO_H
