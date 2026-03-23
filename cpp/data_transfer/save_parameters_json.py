import os
import pandas as pd
import json

variantFactor = 1.4

parameter_dict_default = {
    "TimeExposed": [[2.67, 2.67, 2.67, 2.67, 2.67, 2.67], [4., 4., 4., 4., 4., 4.]],
    "TimeInfectedNoSymptoms": [[1.2, 1.2, 1.2, 1.2, 1.2, 1.2], [2.53, 2.53, 2.53, 2.53, 2.53, 2.53]],
    "TimeInfectedSymptoms": [[5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465], [8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085]],
    "TimeInfectedSevere": [[
        3.925, 3.925, 4.85, 6.4, 7.2, 9.], [6.075, 6.075, 7., 8.7, 9.8, 13.]],
    "TimeInfectedCritical": [[4.95, 4.95, 4.86, 14.14, 14.4, 10.], [8.95, 8.95, 8.86, 20.58, 19.8, 13.2]],
    "TransmissionProbabilityOnContact": [[
        0.02 * variantFactor, 0.05 * variantFactor, 0.05 * variantFactor,
        0.05 * variantFactor, 0.08 * variantFactor, 0.1 * variantFactor
    ], [
        0.04 * variantFactor, 0.07 * variantFactor, 0.07 * variantFactor,
        0.07 * variantFactor, 0.10 * variantFactor, 0.15 * variantFactor
    ]],
    "RelativeTransmissionNoSymptoms": [[0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]],
    "RiskOfInfectionFromSymptomatic": [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.2, 0.2, 0.2, 0.2, 0.2, 0.2]],
    "MaxRiskOfInfectionFromSymptomatic": [[0.4, 0.4, 0.4, 0.4, 0.4, 0.4], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]],
    "RecoveredPerInfectedNoSymptoms": [[
        0.2, 0.2, 0.15, 0.15, 0.15, 0.15], [
        0.3, 0.3, 0.25, 0.25, 0.25, 0.25]],
    "SeverePerInfectedSymptoms": [[
        0.006, 0.006, 0.015, 0.049, 0.15, 0.20], [
        0.009, 0.009, 0.023, 0.074, 0.18, 0.25]],
    "CriticalPerSevere": [[0.05, 0.05, 0.05, 0.10, 0.25, 0.35], [0.10, 0.10, 0.10, 0.20, 0.35, 0.45]],
    "DeathsPerCritical": [[0.00, 0.00, 0.10, 0.10, 0.30, 0.5], [0.10, 0.10, 0.18, 0.18, 0.50, 0.7]],
    "ReducedExposedPartialImmunity": [[0.75, 0.75, 0.75, 0.75, 0.75, 0.75], [0.85, 0.85, 0.85, 0.85, 0.85, 0.85]],
    "ReducedExposedImprovedImmunity": [[0.281, 0.281, 0.281, 0.281, 0.281, 0.281], [0.381, 0.381, 0.381, 0.381, 0.381, 0.381]],
    "ReducedInfectedSymptomsPartialImmunity": [[0.6, 0.6, 0.6, 0.6, 0.6, 0.6], [0.7, 0.7, 0.7, 0.7, 0.7, 0.7]],
    "ReducedInfectedSymptomsImprovedImmunity": [[0.193, 0.193, 0.193, 0.193, 0.193, 0.193], [0.293, 0.293, 0.293, 0.293, 0.293, 0.293]],
    "ReducedInfectedSevereCriticalDeadPartialImmunity": [[0.05, 0.05, 0.05, 0.05, 0.05, 0.05], [0.15, 0.15, 0.15, 0.15, 0.15, 0.15]],
    "ReducedInfectedSevereCriticalDeadImprovedImmunity": [[0.041, 0.041, 0.041, 0.041, 0.041, 0.041], [0.141, 0.141, 0.141, 0.141, 0.141, 0.141]],
    "ReducedTimeInfectedMild": [[1., 1., 1., 1., 1., 1.], [1., 1., 1., 1., 1., 1.]],
    "Seasonality": [[0.1, 0.1, 0.1, 0.1, 0.1, 0.1], [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]]
}

if __name__ == "__main__":

    with open(os.path.join("temp_dir", "parameter_list_lha_direct.json"), 'w') as f:
        json.dump(parameter_dict_default, f, ensure_ascii=False)
