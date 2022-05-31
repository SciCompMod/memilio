/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele
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
#ifndef SECIRV_INFECTIONSTATE_H
#define SECIRV_INFECTIONSTATE_H
namespace mio
{
namespace osecirvvs
{

    /**
    * @brief The InfectionState enum describes the possible
    * categories for the infectious state of persons.
    * Enum is usable as an index, e.g. in CustomIndexArray.
    * Suffix `Naive` means no immunity through vaccination or infection.
    * Suffix `PartiaImmunity` means vaccinated once.
    * Suffix `ImprovedImmunity` means vaccinated twice or recovered from infection.
    * Suffix `Confirmed` means infection has been confirmed by a test, e.g. during commute.
    */
    enum class InfectionState
    {
        SusceptibleNaive = 0,
        SusceptiblePartialImmunity, //ImprovedImmunity is included in Recovered
        ExposedNaive,
        ExposedPartialImmunity,
        ExposedImprovedImmunity,
        CarrierNaive,
        CarrierPartialImmunity,
        CarrierImprovedImmunity ,
        CarrierNaiveConfirmed,
        CarrierPartialImmunityConfirmed,
        CarrierImprovedImmunityConfirmed,
        InfectedNaive,
        InfectedPartialImmunity,
        InfectedImprovedImmunity,
        InfectedNaiveConfirmed,
        InfectedPartialImmunityConfirmed,
        InfectedImprovedImmunityConfirmed,
        HospitalizedNaive,
        HospitalizedPartialImmunity,
        HospitalizedImprovedImmunity,
        ICUNaive,
        ICUPartialImmunity,
        ICUImprovedImmunity,
        Recovered, //includes all with improved immunity, either through infection or at least two vaccinations
        Dead, //no division by immunity
        TotalInfections, //total number of infections during the simulation, for data tracking only, does not contribute to simulation
        
        Count
    };

} // namespace osecirvvs
} // namespace mio

#endif //SECIRV_INFECTIONSTATE_H
