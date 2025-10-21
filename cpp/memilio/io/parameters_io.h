/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef MEMILIO_IO_PARAMETER_H
#define MEMILIO_IO_PARAMETER_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "json/value.h"
#include <string>
#include <vector>

namespace mio
{
/**
 * @brief Gets the region ID (county, state, or district) of an EpiDataEntry. If none are available, 
 * it defaults to 0 which is representing the whole country.
 * 
 * If none are available, it defaults to 0 which is representing the whole country.
 * 
 * @tparam EpiDataEntry The type of the data entry.
 * @param data_entry The (RKI) data entry to extract the region ID from.
 * @return The region ID as integer, or 0 if no specific region information is available.
 */
template <class EpiDataEntry>
int get_region_id(const EpiDataEntry& data_entry)
{
    return data_entry.county_id ? data_entry.county_id->get()
                                : (data_entry.state_id ? data_entry.state_id->get()
                                                       : (data_entry.district_id ? data_entry.district_id->get() : 0));
}

/**
 * @brief Extracts the number of individuals in critical condition (ICU) for each region 
 * on a specified date from the provided DIVI data-
 *
 * @param[in] divi_data Vector of DIVI data entries containing date, region, and ICU information.
 * @param[in] vregion Vector of region IDs for which the data is computed.
 * @param[in] date Date for which the ICU data is computed.
 * @return An IOResult containing a vector with the number of ICU cases for each region, or an 
 *         error if the function fails.
 */
IOResult<std::vector<ScalarType>> compute_divi_data(const std::vector<DiviEntry>& divi_data, const std::vector<int>& vregion, Date date);

/**
 * @brief Reads DIVI data from a file and computes the ICU data for specified regions and date.
 *
 * @param[in] path Path to the file containing DIVI data.
 * @param[in] vregion Vector of region IDs for which the data is computed.
 * @param[in] date Date for which the ICU data is computed.
 * @return An IOResult containing a vector with the number of ICU cases for each region, or an 
 *         error if the function fails.
 */
IOResult<std::vector<ScalarType>> read_divi_data(const std::string& path, const std::vector<int>& vregion, Date date);

/**
 * @brief Sets ICU data from DIVI data into the a vector of models, distributed across age groups.
 *
 * This function reads DIVI data from a file, computes the number of individuals in critical condition (ICU)
 * for each region, and sets these values in the model. The ICU cases are distributed across age groups
 * using the transition probabilities from severe to critical.
 *
 * @tparam Model The type of the model used.
 * @tparam FP Floating point type (default: double).
 *
 * @param[in,out] model Vector of models, each representing a region, where the ICU population is updated.
 * @param[in] num_icu vector of icu data
 * @param[in] vregion Vector of region IDs for which the data is computed.
 * @param[in] date Date for which the ICU data is computed.
 * @param[in] scaling_factor_icu Scaling factor for reported ICU cases.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP, class Model>
IOResult<void> set_divi_data(mio::VectorRange<Model>& model, const std::vector<double>& num_icu, const std::vector<int>& vregion,
                             Date date, FP scaling_factor_icu)
{
    std::vector<FP> sum_mu_I_U(vregion.size(), 0);
    std::vector<std::vector<FP>> mu_I_U{model.size()};
    for (size_t region = 0; region < vregion.size(); region++) {
        auto num_groups = model[region].parameters.get_num_groups();
        for (auto i = AgeGroup(0); i < num_groups; i++) {
            sum_mu_I_U[region] += model[region].parameters.template get<CriticalPerSevere<FP>>()[i] *
                                  model[region].parameters.template get<SeverePerInfectedSymptoms<FP>>()[i];
            mu_I_U[region].push_back(model[region].parameters.template get<CriticalPerSevere<FP>>()[i] *
                                     model[region].parameters.template get<SeverePerInfectedSymptoms<FP>>()[i]);
        }
    }

    for (size_t region = 0; region < vregion.size(); region++) {
        auto num_groups = model[region].parameters.get_num_groups();
        for (auto i = AgeGroup(0); i < num_groups; i++) {
            model[region].populations[{i, InfectionState::InfectedCriticalNaive}] =
                scaling_factor_icu * num_icu[region] * mu_I_U[region][(size_t)i] / sum_mu_I_U[region];
        }
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
template <class Model>
IOResult<void> set_divi_data(mio::VectorRange<Model>& model, const std::string& path, const std::vector<int>& vregion,
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
    BOOST_OUTCOME_TRY(set_divi_data(model, num_icu, vregion, date, scaling_factor_icu));
    return success();
}

/**
 * @brief Reads population data from a vector of population data entries.
 * 
 * @param[in] population_data Vector of population data entries.
 * @param[in] vregion Vector of keys representing the regions of interest.
 * @return An IOResult containing a vector of vectors, where each inner vector represents the population
 *         distribution across age groups for a specific region, or an error if the function fails.
 */
IOResult<std::vector<std::vector<ScalarType>>>
read_population_data(const std::vector<PopulationDataEntry>& population_data, const std::vector<int>& vregion);

/**
 * @brief Reads population data from census data.
 * 
 * @param[in] path Path to the population data file.
 * @param[in] vregion Vector of keys representing the regions of interest.
 * @return An IOResult containing a vector of vectors, where each inner vector represents the population
 *         distribution across age groups for a specific region, or an error if the function fails.
 */
IOResult<std::vector<std::vector<ScalarType>>> read_population_data(const std::string& path,
                                                                    const std::vector<int>& vregion);

} // namespace mio

#endif //MEMILIO_HAS_JSONCPP

#endif //MEMILIO_IO_PARAMETER_H
