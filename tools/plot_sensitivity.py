import matplotlib.pyplot as plt
import numpy as np

# Liste der Elementary Effects für die Parameter
elementary_effects_daily_infections = [-418.244, 613.054, 416.906, -65.8956, -145.059, -116.708, -45.795, 195.713, -38.5292,
                                       21533.7, -288.471, 305.075, 2808, -462.085, -19.1413, 284.15, 115.067, 592.906,
                                       3656.1, 296.339, 540.716, 201.785, 545.794, 1326.06]

# Liste der Parameternamen
parameter_names = [
    "TimeExposed", "TimeInfectedNoSymptoms", "TimeInfectedSymptoms", "TimeInfectedSevere",
    "TimeInfectedCritical", "TimeTemporaryImmunityPI", "TimeTemporaryImmunityII", "TimeWaningPartialImmunity",
    "TimeWaningImprovedImmunity", "TransmissionProbabilityOnContact", "RelativeTransmissionNoSymptoms",
    "RiskOfInfectionFromSymptomatic", "MaxRiskOfInfectionFromSymptomatic", "RecoveredPerInfectedNoSymptoms",
    "SeverePerInfectedSymptoms", "CriticalPerSevere", "DeathsPerCritical", "ReducExposedPartialImmunity",
    "ReducExposedImprovedImmunity", "ReducInfectedSymptomsPartialImmunity", "ReducInfectedSymptomsImprovedImmunity",
    "ReducInfectedSevereCriticalDeadPartialImmunity", "ReducInfectedSevereCriticalDeadImprovedImmunity",
    "ReducTimeInfectedMild"
]

# Plotting
plt.figure(figsize=(10, 8))
plt.barh(parameter_names, elementary_effects_daily_infections, color='skyblue')
plt.xlabel('Elementary Effect')
plt.ylabel('Parameter')
plt.title('Elementary Effects of Parameters on New Infections')
plt.grid(True)
plt.tight_layout()
plt.show()
