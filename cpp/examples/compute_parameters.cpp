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
    /** With this file the parameters for a ODE-SECIR model without age distribution can be calculated to match those specified in the first published paper
     (https://doi.org/10.1016/j.mbs.2021.108648). 
     
    For this, the parameters from the paper are used first to derive just one transitiontime if 2 are specified, eg T_C=\mu_C^R*T_C^R+(1-\mu_C^R)*T_C^I.
    For each age group, average values are calculated if lower and upper bounds are given. 
    For this we assume a uniform distribution so that T=(T_min+T_max)/2.
    Finally we calculate a weighted average time across the age groups.
    */
    bool printResult = true;

    // Age group sizes are calculated using table number 12411-04-02-4-B from www.regionalstatistik.de for the date 31.12.2020
    const double age_group_sizes[] = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
    const int total                = 83155031.0;
    const int numagegroups         = 6;

    // transmission parameters
    const double transmissionProbabilityOnContactMin[] = {0.02, 0.05, 0.05, 0.05, 0.08, 0.15};
    const double transmissionProbabilityOnContactMax[] = {0.04, 0.07, 0.07, 0.07, 0.10, 0.20};

    double transmissionProbabilityOnContact = 0;
    for (int i = 0; i < numagegroups; i++) {
        transmissionProbabilityOnContact +=
            (age_group_sizes[i] / total) * 0.5 *
            (transmissionProbabilityOnContactMin[i] + transmissionProbabilityOnContactMax[i]);
    }

    if (printResult) {
        std::cout << "transmissionProbabilityOnContact: " << transmissionProbabilityOnContact << std::endl;
    }

    const double relativeTransmissionNoSymptoms = 1;
    const double riskOfInfectionFromSymptomatic = 0.3;

    if (printResult) {
        std::cout << "relativeTransmissionNoSymptoms: " << relativeTransmissionNoSymptoms << std::endl;
        std::cout << "riskOfInfectionFromSymptomatic: " << riskOfInfectionFromSymptomatic << std::endl;
    }

    // E
    const double timeExposedMin = 2.67;
    const double timeExposedMax = 4.00;

    double timeExposed = 0.5 * (timeExposedMin + timeExposedMax);

    if (printResult) {
        std::cout << "timeExposed: " << timeExposed << std::endl;
    }

    // Calculate parameters for I first because a value of I is needed for C
    // I
    const double timeInfectedSymptomstoRecoveredMin        = 5.6;
    const double timeInfectedSymptomstoRecoveredMax        = 8.4;
    const double timeInfectedSymptomstoInfectedSevereMin[] = {9, 9, 9, 5, 5, 5};
    const double timeInfectedSymptomstoInfectedSevereMax[] = {12, 12, 12, 7, 7, 7};
    const double severePerInfectedSymptomsMin[]            = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20};
    const double severePerInfectedSymptomsMax[]            = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25};

    double timeInfectedSymptomsMindummy;
    double timeInfectedSymptomsMaxdummy;
    double severePerInfectedSymptomsdummy;
    double timeInfectedSymptoms      = 0;
    double severePerInfectedSymptoms = 0;

    for (int i = 0; i < numagegroups; i++) {
        severePerInfectedSymptomsdummy = 0.5 * (severePerInfectedSymptomsMin[i] + severePerInfectedSymptomsMax[i]);
        timeInfectedSymptomsMindummy   = (1 - severePerInfectedSymptomsdummy) * timeInfectedSymptomstoRecoveredMin +
                                       severePerInfectedSymptomsdummy * timeInfectedSymptomstoInfectedSevereMin[i];
        timeInfectedSymptomsMaxdummy = (1 - severePerInfectedSymptomsdummy) * timeInfectedSymptomstoRecoveredMax +
                                       severePerInfectedSymptomsdummy * timeInfectedSymptomstoInfectedSevereMax[i];

        timeInfectedSymptoms +=
            (age_group_sizes[i] / total) * 0.5 * (timeInfectedSymptomsMindummy + timeInfectedSymptomsMaxdummy);
        severePerInfectedSymptoms += (age_group_sizes[i] / total) * severePerInfectedSymptomsdummy;
    }

    // C
    const double timeInfectedNoSymptomstoInfectedSymptoms = 5.2 - timeExposed;
    const double timeInfectedNoSymptomstoRecovered =
        timeInfectedNoSymptomstoInfectedSymptoms +
        0.5 * (timeInfectedSymptomstoRecoveredMin + timeInfectedSymptomstoRecoveredMax);
    const double recoveredPerInfectedNoSymptomsMin[] = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15};
    const double recoveredPerInfectedNoSymptomsMax[] = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25};

    double timeInfectedNoSymptoms         = 0;
    double recoveredPerInfectedNoSymptoms = 0;

    for (int i = 0; i < numagegroups; i++) {
        recoveredPerInfectedNoSymptoms += (age_group_sizes[i] / total) * 0.5 *
                                          (recoveredPerInfectedNoSymptomsMin[i] + recoveredPerInfectedNoSymptomsMax[i]);
    }

    timeInfectedNoSymptoms = recoveredPerInfectedNoSymptoms * timeInfectedNoSymptomstoRecovered +
                             (1 - recoveredPerInfectedNoSymptoms) * timeInfectedNoSymptomstoInfectedSymptoms;

    if (printResult) {
        std::cout << "timeInfectedNoSymptoms: " << timeInfectedNoSymptoms << std::endl;
        std::cout << "recoveredPerInfectedNoSymptoms: " << recoveredPerInfectedNoSymptoms << std::endl;
        std::cout << "timeInfectedSymptoms: " << timeInfectedSymptoms << std::endl;
        std::cout << "severePerInfectedSymptoms: " << severePerInfectedSymptoms << std::endl;
    }

    // H
    const double timeInfectedSeveretoRecoveredMin[]      = {4, 4, 5, 7, 9, 13};
    const double timeInfectedSeveretoRecoveredMax[]      = {6, 6, 7, 9, 11, 17};
    const double timeInfectedSeveretoInfectedCriticalMin = 3;
    const double timeInfectedSeveretoInfectedCriticalMax = 7;
    const double criticalPerSevereMin[]                  = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35};
    const double criticalPerSevereMax[]                  = {0.10, 0.10, 0.10, 0.20, 0.35, 0.45};

    double timeInfectedSevereMindummy;
    double timeInfectedSevereMaxdummy;
    double criticalPerSeveredummy;
    double timeInfectedSevere = 0;
    double criticalPerSevere  = 0;

    for (int i = 0; i < numagegroups; i++) {
        criticalPerSeveredummy     = 0.5 * (criticalPerSevereMin[i] + criticalPerSevereMax[i]);
        timeInfectedSevereMindummy = (1 - criticalPerSeveredummy) * timeInfectedSeveretoRecoveredMin[i] +
                                     criticalPerSeveredummy * timeInfectedSeveretoInfectedCriticalMin;
        timeInfectedSevereMaxdummy = (1 - criticalPerSeveredummy) * timeInfectedSeveretoRecoveredMax[i] +
                                     criticalPerSeveredummy * timeInfectedSeveretoInfectedCriticalMax;

        timeInfectedSevere +=
            (age_group_sizes[i] / total) * 0.5 * (timeInfectedSevereMindummy + timeInfectedSevereMaxdummy);
        criticalPerSevere += (age_group_sizes[i] / total) * criticalPerSeveredummy;
    }

    if (printResult) {
        std::cout << "timeInfectedSevere: " << timeInfectedSevere << std::endl;
        std::cout << "criticalPerSevere: " << criticalPerSevere << std::endl;
    }

    // U
    const double timeInfectedCriticaltoRecoveredMin[] = {5, 5, 5, 14, 14, 10};
    const double timeInfectedCriticaltoRecoveredMax[] = {9, 9, 9, 21, 21, 15};
    const double timeInfectedCriticaltoDeadMin[]      = {4, 4, 4, 15, 15, 10};
    const double timeInfectedCriticaltoDeadMax[]      = {8, 8, 8, 18, 18, 12};
    const double deathsPerCriticalMin[]               = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5};
    const double deathsPerCriticalMax[]               = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7};

    double timeInfectedCriticalMindummy;
    double timeInfectedCriticalMaxdummy;
    double deathsPerCriticaldummy;
    double timeInfectedCritical = 0;
    double deathsPerCritical    = 0;

    for (int i = 0; i < numagegroups; i++) {
        deathsPerCriticaldummy       = 0.5 * (deathsPerCriticalMin[i] + deathsPerCriticalMax[i]);
        timeInfectedCriticalMindummy = (1 - deathsPerCriticaldummy) * timeInfectedCriticaltoRecoveredMin[i] +
                                       deathsPerCriticaldummy * timeInfectedCriticaltoDeadMin[i];
        timeInfectedCriticalMaxdummy = (1 - deathsPerCriticaldummy) * timeInfectedCriticaltoRecoveredMax[i] +
                                       deathsPerCriticaldummy * timeInfectedCriticaltoDeadMax[i];

        timeInfectedCritical +=
            (age_group_sizes[i] / total) * 0.5 * (timeInfectedCriticalMindummy + timeInfectedCriticalMaxdummy);
        deathsPerCritical += (age_group_sizes[i] / total) * deathsPerCriticaldummy;
    }

    if (printResult) {
        std::cout << "timeInfectedCritical: " << timeInfectedCritical << std::endl;
        std::cout << "deathsPerCritical: " << deathsPerCritical << std::endl;
    }
}