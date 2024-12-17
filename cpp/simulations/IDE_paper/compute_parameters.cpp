/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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

#include <iostream>

int main()
{
    /** With this file parameters without age distribution can be calculated to match those specified in the
     covasim paper (https://doi.org/10.1371/journal.pcbi.1009149).
    First, we calculate a weighted average time across the age groups.
    If other probabilites than required are given, we calculate the right probabilities.
    */

    // Age group sizes are calculated using table number 12411-04-02-4-B from www.regionalstatistik.de for the date 31.12.2020.
    const double age_group_sizes[] = {7752706.0, 7581868,  9483430, 10871964, 10070748,
                                      13304542,  10717241, 7436098, 5092743,  843691};
    const int total                = 83155031;
    const int numagegroups         = 10;

    // Calculate value for probability InfectedSymptomsPerInfectedNoSymptoms.
    const double InfectedSymptomsPerInfectedNoSymptoms[] = {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.90};
    double resultInfectedSymptomsPerInfectedNoSymptoms   = 0;
    for (int i = 0; i < numagegroups; i++) {
        resultInfectedSymptomsPerInfectedNoSymptoms += age_group_sizes[i] * InfectedSymptomsPerInfectedNoSymptoms[i];
    }
    resultInfectedSymptomsPerInfectedNoSymptoms = resultInfectedSymptomsPerInfectedNoSymptoms / total;

    std::cout << "InfectedSymptomsPerInfectedNoSymptoms: " << resultInfectedSymptomsPerInfectedNoSymptoms << std::endl;

    // Calculate value for probability SeverePerInfectedSymptoms.
    const double SeverePerInfectedNoSymptoms[] = {0.00050, 0.00165, 0.00720, 0.02080, 0.03430,
                                                  0.07650, 0.13280, 0.20655, 0.24570, 0.24570};
    double average_SeverePerInfectedNoSymptoms = 0;
    for (int i = 0; i < numagegroups; i++) {
        average_SeverePerInfectedNoSymptoms += age_group_sizes[i] * SeverePerInfectedNoSymptoms[i];
    }
    average_SeverePerInfectedNoSymptoms = average_SeverePerInfectedNoSymptoms / total;
    double resultSeverePerInfectedSymptoms =
        average_SeverePerInfectedNoSymptoms / resultInfectedSymptomsPerInfectedNoSymptoms;

    std::cout << "SeverePerInfectedSymptoms: " << resultSeverePerInfectedSymptoms << std::endl;

    // Calculate value for probability CriticalPerSevere.
    const double CriticalPerInfectedNoSymptoms[] = {0.00003, 0.00008, 0.00036, 0.00104, 0.00216,
                                                    0.00933, 0.03639, 0.08923, 0.17420, 0.17420};
    double average_CriticalPerInfectedNoSymptoms = 0;
    for (int i = 0; i < numagegroups; i++) {
        average_CriticalPerInfectedNoSymptoms += age_group_sizes[i] * CriticalPerInfectedNoSymptoms[i];
    }
    average_CriticalPerInfectedNoSymptoms = average_CriticalPerInfectedNoSymptoms / total;
    double resultCriticalPerSevere        = average_CriticalPerInfectedNoSymptoms / average_SeverePerInfectedNoSymptoms;

    std::cout << "CriticalPerSevere: " << resultCriticalPerSevere << std::endl;

    // Calculate value for probability DeathsPerCritical.
    const double DeathsPerInfectedNoSymptoms[] = {0.00002, 0.00002, 0.00010, 0.00032, 0.00098,
                                                  0.00265, 0.00766, 0.02439, 0.08292, 0.16190};
    double average_DeathsPerInfectedNoSymptoms = 0;
    for (int i = 0; i < numagegroups; i++) {
        average_DeathsPerInfectedNoSymptoms += age_group_sizes[i] * DeathsPerInfectedNoSymptoms[i];
    }
    average_DeathsPerInfectedNoSymptoms = average_DeathsPerInfectedNoSymptoms / total;
    double resultDeathsPerCritical      = average_DeathsPerInfectedNoSymptoms / average_CriticalPerInfectedNoSymptoms;

    std::cout << "DeathsPerCritical: " << resultDeathsPerCritical << std::endl;

    // ---- Stay times for the ODE model. ----
    std::cout << "\nAveraged stay times for the ODE model: " << std::endl;
    std::cout << "TimeExposed: " << 4.5 << std::endl;
    std::cout << "TimeInfectedNoSymptoms: "
              << 1.1 * resultInfectedSymptomsPerInfectedNoSymptoms +
                     8. * (1. - resultInfectedSymptomsPerInfectedNoSymptoms)
              << std::endl;
    std::cout << "TimeInfectedSymptoms: "
              << 6.6 * resultSeverePerInfectedSymptoms + 8. * (1. - resultSeverePerInfectedSymptoms) << std::endl;
    std::cout << "TimeInfectedSevere: " << 1.5 * resultCriticalPerSevere + 18.1 * (1. - resultCriticalPerSevere)
              << std::endl;
    std::cout << "TimeInfectedCritical: " << 10.7 * resultDeathsPerCritical + 18.1 * (1. - resultDeathsPerCritical)
              << std::endl;
}
