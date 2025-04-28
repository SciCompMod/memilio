/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, Khoa Nguyen, David Kerkmann, Rene Schmieding
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

#ifndef MIO_ABM_MODEL_FUNCTIONS_H
#define MIO_ABM_MODEL_FUNCTIONS_H

#include "abm/location.h"
#include "abm/person.h"

#include <vector>

namespace mio
{
namespace abm
{

/**
 * @brief Compute the number of daily transmissions for contact transmission of a virus in a cell.
 * @param[in] rates The local exposure rates.
 * @param[in] cell_index Cell index of the Cell.
 * @param[in] virus VirusVariant of interest.
 * @param[in] age_receiver AgeGroup of the receiving Person.
 * @param[in] params The local infection parameters.
 * @return Average amount of Infection%s with the virus from the AgeGroup of the transmitter per day.
 */
ScalarType daily_transmissions_by_contacts(const ContactExposureRates& rates, const CellIndex cell_index,
                                           const VirusVariant virus, const AgeGroup age_receiver,
                                           const LocalInfectionParameters& params);

/**
 * @brief Compute the number of daily transmissions for aerosol transmission of a virus in a cell.
 * @param[in] rates The local exposure rates.
 * @param[in] cell_index Cell index of the Cell.
 * @param[in] virus VirusVariant of interest.
 * @param[in] global_params The parameter set of the Model.
 * @return Average amount of Infection%s with the virus per day.
 */
ScalarType daily_transmissions_by_air(const AirExposureRates& rates, const CellIndex cell_index,
                                      const VirusVariant virus, const Parameters& global_params);

/**
 * @brief Add the contribution of a person to the local exposure rates.
 * @param[in, out] local_air_exposure Exposure rates by aerosols for the local population.
 * @param[in, out] local_contact_exposure Exposure by rates contacts for the local population.
 * @param[in] person A person from the local population.
 * @param[in] location The person's current location.
 * @param[in] t Current Simulation time.
 * @param[in] dt Length of the current Simulation time step.
 */
void add_exposure_contribution(AirExposureRates& local_air_exposure, ContactExposureRates& local_contact_exposure,
                               const Person& person, const Location& location, const TimePoint t, const TimeSpan dt);

/**
 * @brief Let a Person interact with the population at its current Location, possibly getting infected.
 * @param[in, out] rng PersonalRandomNumberGenerator for this Person.
 * @param[in, out] person The person to interact with the local population.
 * @param[in] location The person's current location.
 * @param[in] local_air_exposure Precomputed exposure rates by aerosols for the local population.
 * @param[in] local_contact_exposure Precomputed exposure by rates contacts for the local population.
 * @param[in] t Current Simulation time.
 * @param[in] dt Length of the current Simulation time step.
 * @param[in] global_parameters Parameters of the Model.
 */
void interact(PersonalRandomNumberGenerator& personal_rng, Person& person, const Location& location,
              const AirExposureRates& local_air_exposure, const ContactExposureRates& local_contact_exposure,
              const TimePoint t, const TimeSpan dt, const Parameters& global_parameters);
/**
 * @brief Change a persons location to another location.
 * If the person already is at the destination, neither mode nor cells are set.
 * @param[in, out] person The person to change location.
 * @param[in] destination The destination to change location to.
 * @param[in] mode The transport mode the person uses to change location.
 * @param[in] cells The cells within the destination the person should be in.
 * @return Returns false if the person already is at the given destination, true otherwise.
 */
bool change_location(Person& person, const Location& destination, const TransportMode mode = TransportMode::Unknown,
                     const std::vector<uint32_t>& cells = {0});

/**
 * @brief Adjust ContactRates of location by MaximumContacts.
 * Every ContactRate is adjusted by the proportion MaximumContacts of the location has on the total 
 * number of contacts according to the ContactRates.
 * @param[in, out] location The location whose ContactRates are adjusted.
 * @param[in] num_agegroup The number of AgeGroups in the model.
 */
void adjust_contact_rates(Location& location, size_t num_agegroups);

} // namespace abm
} // namespace mio

#endif // MIO_ABM_MODEL_FUNCTIONS_H
