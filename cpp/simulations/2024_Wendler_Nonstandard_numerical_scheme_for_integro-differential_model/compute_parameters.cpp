/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Lena Ploetzke, Anna Wendler
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
    /** Here, the parameters for the probability of transitioning from InfectedNoSymptoms to InfectedSymptoms, from 
    InfectedSymptoms to InfectedSevere and from InfectedSevere to InfectedCritical without age resolution are 
    calculated based on the values in the paper "Assessment of effective mitigation and prediction of the spread of 
    SARS-CoV-2 in Germany using demographic information and spatial resolution" 
    (https://doi.org/10.1016/j.mbs.2021.108648).
    We do this by calculating a weighted average across the age groups. Here, we use 6 age groups according to the 
    paper.
    */

    // Age group sizes are calculated using table number 12411-04-02-4-B from www.regionalstatistik.de for the date 31.12.2020.
    double age_group_sizes[] = {3969138, 7508662, 18921292, 28666166, 18153339, 5936434};
    int total                = 83155031;
    int numagegroups         = 6;

    const double RecoveredPerInfectedNoSymptoms[]    = {0.25, 0.25, 0.20, 0.20, 0.20, 0.20};
    const double InfectedSeverePerInfectedSymptoms[] = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
    const double InfectedCriticalPerInfectedSevere[] = {0.075, 0.075, 0.075, 0.15, 0.30, 0.40};

    double resultRecoveredPerInfectedNoSymptoms    = 0.;
    double resultInfectedSeverePerInfectedSymptoms = 0.;
    double resultInfectedCriticalPerInfectedSevere = 0.;
    for (int i = 0; i < numagegroups; i++) {
        resultRecoveredPerInfectedNoSymptoms += age_group_sizes[i] * RecoveredPerInfectedNoSymptoms[i];
        resultInfectedSeverePerInfectedSymptoms += age_group_sizes[i] * InfectedSeverePerInfectedSymptoms[i];
        resultInfectedCriticalPerInfectedSevere += age_group_sizes[i] * InfectedCriticalPerInfectedSevere[i];
    }
    double resultInfectedSymptomsPerInfectedNoSymptoms = 1. - resultRecoveredPerInfectedNoSymptoms / total;
    resultInfectedSeverePerInfectedSymptoms            = resultInfectedSeverePerInfectedSymptoms / total;
    resultInfectedCriticalPerInfectedSevere            = resultInfectedCriticalPerInfectedSevere / total;

    std::cout << "InfectedSymptomsPerInfectedNoSymptoms: " << resultInfectedSymptomsPerInfectedNoSymptoms << std::endl;
    std::cout << "InfectedSeverePerInfectedSymptoms: " << resultInfectedSeverePerInfectedSymptoms << std::endl;
    std::cout << "InfectedCriticalPerInfectedSevere: " << resultInfectedCriticalPerInfectedSevere << std::endl;

    /** Here, the parameter for the probability of transitioning from the ICU compartment to the Dead compartment 
    without age resolution is calculated based on the values in the Covasim paper (https://doi.org/10.1371/journal.pcbi.1009149).
    We do this by calculating a weighted average across the age groups. In the following, we use 10 age groups 
    according to the Covasim paper.
    In the paper, probabilities for other transitions than required for our model are given, hence we calculate the 
    probabilities for the necessary transitions as described below.
    */

    // Age group sizes are calculated using table number 12411-04-02-4-B from www.regionalstatistik.de for the date 31.12.2020.
    double age_group_sizes_covasim[] = {7752706.0, 7581868,  9483430, 10871964, 10070748,
                                        13304542,  10717241, 7436098, 5092743,  843691};
    int numagegroups_covasim         = 10;

    // Calculate value for probability DeathsPerCritical.
    // This is done by first computing the average CriticalPerInfectedNoSymptoms and the average
    // probability DeathsPerInfectedNoSymptoms based on the values from the Covasim paper. With this, we can compute
    // the probability DeathsPerCritical that we need for our simulation.

    // Calculate value for probability CriticalPerInfectedNoSymptoms.
    const double CriticalPerInfectedNoSymptoms[] = {0.00003, 0.00008, 0.00036, 0.00104, 0.00216,
                                                    0.00933, 0.03639, 0.08923, 0.17420, 0.17420};
    double average_CriticalPerInfectedNoSymptoms = 0;
    for (int i = 0; i < numagegroups_covasim; i++) {
        average_CriticalPerInfectedNoSymptoms += age_group_sizes_covasim[i] * CriticalPerInfectedNoSymptoms[i];
    }
    average_CriticalPerInfectedNoSymptoms = average_CriticalPerInfectedNoSymptoms / total;

    // Calculate value for probability DeathsPerInfectedNoSymptoms.
    const double DeathsPerInfectedNoSymptoms[] = {0.00002, 0.00002, 0.00010, 0.00032, 0.00098,
                                                  0.00265, 0.00766, 0.02439, 0.08292, 0.16190};
    double average_DeathsPerInfectedNoSymptoms = 0;
    for (int i = 0; i < numagegroups_covasim; i++) {
        average_DeathsPerInfectedNoSymptoms += age_group_sizes_covasim[i] * DeathsPerInfectedNoSymptoms[i];
    }
    average_DeathsPerInfectedNoSymptoms = average_DeathsPerInfectedNoSymptoms / total;

    // Calculate value for probability DeathsPerCritical.
    double resultDeathsPerCritical = average_DeathsPerInfectedNoSymptoms / average_CriticalPerInfectedNoSymptoms;

    std::cout << "DeathsPerCritical: " << resultDeathsPerCritical << std::endl;
}
