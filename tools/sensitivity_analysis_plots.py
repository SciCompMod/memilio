import matplotlib.pyplot as plt
import numpy as np

# Liste der Elementary Effects für neue Infektionen (kopiere die Werte aus deinem Ergebnis
elementary_effects_infections = [-176.201, 408.631, 94.0238, -4.84166, -0.607096, -7.88622, -36.1108, -0.728007,
                                 -0.781813, 20479.9, 1.26494, 59.3305, -471.647, -275.052, -121.037, -64.3411,
                                 548.504, 4205.08, 222.71, 1073.79, -51.827, -61.827, 1485.26]
# reverse
elementary_effects_infections = elementary_effects_infections[::-1]

# Liste der Elementary Effects für ICU-Belegung (kopiere die Werte aus deinem Ergebnis)
elementary_effects_icu = [-26.6798, 3.19165, 0.417507, -3.17338, 64.0284, -0.0942643, -1.44256, -0.0992241, -0.0205701, 1201.36,
                          0.913348, -2.22022, -329.531, 2093.55, 2366, -11.6818, -
                          46.2347, -191.932, 30.6898, 75.0195, 286.659,
                          336.659, 7.19287]
# reverse
elementary_effects_icu = elementary_effects_icu[::-1]

# Liste der Parameternamen
# parameter_names = [
#     "TimeExposed", "TimeInfectedNoSymptoms", "TimeInfectedSymptoms", "TimeInfectedSevere",
#     "TimeInfectedCritical", "TimeTemporaryImmunityPI", "TimeTemporaryImmunityII", "TimeWaningPartialImmunity",
#     "TimeWaningImprovedImmunity", "TransmissionProbabilityOnContact", "RelativeTransmissionNoSymptoms",
#     "RiskOfInfectionFromSymptomatic", "MaxRiskOfInfectionFromSymptomatic", "RecoveredPerInfectedNoSymptoms",
#     "SeverePerInfectedSymptoms", "CriticalPerSevere", "DeathsPerCritical", "ReducExposedPartialImmunity",
#     "ReducExposedImprovedImmunity", "ReducInfectedSymptomsPartialImmunity", "ReducInfectedSymptomsImprovedImmunity",
#     "ReducInfectedSevereCriticalDeadPartialImmunity", "ReducInfectedSevereCriticalDeadImprovedImmunity",
#     "ReducTimeInfectedMild"
# ]

parameter_names = [
    r"$T_{E}$", r"$T_{I_{NS}}$", r"$T_{I_{Sy}}$", r"$T_{I_{Sev}}$",
    r"$T_{I_{Cr}}$", r"$T_{\mathcal{I}_{PI}}$", r"$T_{\mathcal{I}_{II}}$", r"$T_{W_{PI}}$",
    r"$T_{W_{II}}$", r"$\rho$", r"${\xi_{I_{NS}}}$",
    r"$\xi_{I_{Sy}}$", r"$\mu_{I_{NS}}^{I_{Sy}}$",
    r"$\mu_{I_{Sy, N}}^{I_{Sev, N}}$", r"${\mu}_{I_{Sev, N}}^{I_{Cr, N}}$", r"$\mu_{I_{Cr, N}}^{D_{N}}$", r"$p_{{E_{PI}}}$",
    r"$p_{{E_{II}}}$", r"$p_{{I_{Sy,PI}}}$", r"$p_{{I_{Sy,II}}}$",
    r"$p_{I_{Sev, PI}}$", r"$p_{I_{Sev,II}}$",
    r"$\kappa$"
]
# reverse
parameter_names = parameter_names[::-1]

# Erstellen von zwei Subplots
fig, ax = plt.subplots(1, 2, figsize=(18, 10))
fontsize_ticks = 15

# Plot für neue Infektionen
ax[0].barh(parameter_names, elementary_effects_infections, color='skyblue')
ax[0].set_xlabel('Elementary Effect on New Infections')
ax[0].set_ylabel('Parameter')
ax[0].set_title('Elementary Effects of Parameters on New Infections')
ax[0].set_yticklabels(parameter_names, fontsize=fontsize_ticks)
ax[0].grid(True)
ax[0].set_xscale('symlog')

# Plot für ICU-Belegung
ax[1].barh(parameter_names, elementary_effects_icu, color='lightcoral')
ax[1].set_xlabel('Elementary Effect on ICU Occupancy')
ax[1].set_ylabel('Parameter')
ax[1].set_title('Elementary Effects of Parameters on ICU Occupancy')
ax[1].set_yticklabels(parameter_names, fontsize=fontsize_ticks)
ax[1].grid(True)
ax[1].set_xscale('symlog')


plt.tight_layout()
# plt.show()

plt.savefig('sensitivity_analysis.png')
