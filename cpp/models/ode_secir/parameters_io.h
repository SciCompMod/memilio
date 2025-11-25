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
 * @brief Reads populations data from rki data.
 * @param[in] case_data Vector of ConfirmedCasesDataEntry%s.
 * @param[in] region Key of the region of interest.
 * @param[in] date Date at which the data is read.
 * @param[in, out] num_* Output vector for number of people in the corresponding compartement.
 * @param[in] t_* vector Average time it takes to get from one compartement to another for each age group.
 * @param[in] mu_* vector Probabilities to get from one compartement to another for each age group.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 */
IOResult<void> compute_confirmed_cases_data(
    std::vector<ConfirmedCasesDataEntry>& case_data, const int region, Date date,
    std::vector<ScalarType>& num_Exposed, std::vector<ScalarType>& num_InfectedNoSymptoms,
    std::vector<ScalarType>& num_InfectedSymptoms, std::vector<ScalarType>& num_InfectedSevere,
    std::vector<ScalarType>& num_icu, std::vector<ScalarType>& num_death,
    std::vector<ScalarType>& num_rec, const std::vector<int>& t_Exposed,
    const std::vector<int>& t_InfectedNoSymptoms,
    const std::vector<int>& t_InfectedSymptoms, const std::vector<int>& t_InfectedSevere,
    const std::vector<int>& t_InfectedCritical, const std::vector<ScalarType>& mu_C_R,
    const std::vector<ScalarType>& mu_I_H, const std::vector<ScalarType>& mu_H_U,
    const std::vector<ScalarType>& scaling_factor_inf);

/**
 * @brief Sets populations data from already read case data with multiple age groups into a Model.
 * @param[in, out] model Model in which the data is set.
 * @param[in] case_data Vector of ConfirmedCasesDataEntry%s.
 * @param[in] region Key of the region of interest.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 */
IOResult<void> set_confirmed_cases_data(Model<ScalarType>& model, std::vector<ConfirmedCasesDataEntry>& case_data,
                                        const int region, Date date, const std::vector<ScalarType>& scaling_factor_inf);

/**
 * @brief Sets the infected populations for the given models based on confirmed cases data. 
 * Reads the case data from a file and then calls a subfunction that sets the infected population 
 * for each model.
 * @param[in, out] model VectorRange of Node%s each containing a Model for which the confirmed cases data will be set.
 * @param[in] path Path to the confirmed cases data file.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
 */
IOResult<void> set_confirmed_cases_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path,
                                        Date date, const std::vector<ScalarType>& scaling_factor_inf);

/**
 * @brief Sets population data from census data which has been read into num_population.
 * @param[in, out] model Model in which the data is set.
 * @param[in] num_population Vector of population data for the region of interest.
 * @param[in] region Key of the region of interest.
 */
IOResult<void> set_population_data(Model<ScalarType>& model, const std::vector<ScalarType>& num_population,
                                   const int region);

/**
 * @brief Reads population data from a file and sets it for the each given model.
 * @param[in, out] model VectorRange of Node%s containing a Model for which the confirmed cases data will be set.
 * @param[in] path Path to population data file.
 */
IOResult<void> set_population_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path);

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
 * @param[in, out] model VectorRange of Node%s each containing a Model.
 * @param[in] path Path to transformed DIVI file.
 * @param[in] date Date at which the data is read.
 * @param[in] scaling_factor_icu Scaling factor for reported ICU cases.
 */
IOResult<void> set_divi_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path, Date date, 
                             ScalarType scaling_factor_icu);

} //namespace details


/**
 * @brief Reads compartments for geographic units at a specified date from data files.
 *
 * This function estimates all compartments from available data using the provided model parameters.
 * 
 * @param[in,out] model VectorRange of Node%s each containing a Model to be initialized with data.
 * @param[in] date Date for which the data should be read.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] scaling_factor_icu Scaling factor for ICU cases.
 * @param[in] epidata_filenames Object containing the input data file paths.
 * 
 * @return An IOResult indicating success or failure.
 */
IOResult<void> read_input_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, Date date,
                               const std::vector<ScalarType>& scaling_factor_inf, ScalarType scaling_factor_icu,
                               const mio::regions::de::EpidataFilenames& epidata_filenames);

/**
 * @brief Converts input data from one range of models to another with different type.
 * 
 * @tparam FP Floating point type (default: double).
 * @param[in] model_from VectorRange of Node%s each containing a Model with the input data.
 * @param[in,out] model_to VectorRange of Node%s each containing a Model to be initialized with data.
 *
 * @return An IOResult indicating success or failure.
 */
template<class FP>
void convert_input_data_type(const mio::VectorRange<Node<Model<ScalarType>>>& model_from, const mio::VectorRange<Node<Model<FP>>>& model_to)
{
    assert(model_from.size() == model_to.size());
    assert((size_t)model_from[0].property.parameters.get_num_groups() == (size_t)model_to[0].property.parameters.get_num_groups());
    // Todo: add conversion of ParameterSet and then re-use code from other model parameters io 

    for (size_t region_idx = 0; region_idx < model_from.size(); ++region_idx) {
        // convert populations to mio::UncertainValue<FP>
        // needs 2 converts as mio::UncertainValue<ScalarType> -> mio::UncertainValue<FP> does not work
        model_to[region_idx].property.populations = model_to[region_idx].property.populations.template convert<ScalarType>().template convert<mio::UncertainValue<FP>>();
    }
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
IOResult<void> export_input_data_timeseries(
    const mio::VectorRange<Node<Model<ScalarType>>> model, const std::string& results_dir, Date date,
    const std::vector<ScalarType>& scaling_factor_inf, ScalarType scaling_factor_icu, int num_days,
    const mio::regions::de::EpidataFilenames& epidata_filenames);
#else
IOResult<void> export_input_data_county_timeseries(const mio::VectorRange<Node<Model<ScalarType>>>, const std::string&, Date, const std::vector<int>&,
                                                   const std::vector<ScalarType>&, const ScalarType, const int,
                                                   const mio::regions::de::EpidataFilenames&);
#endif // MEMILIO_HAS_HDF5

} // namespace osecir

} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // ODESECIR_PARAMETERS_IO_H
