
from plot_realistic_scenario import *
import os
import subprocess


def run_simulation(start_date, simulation_parameter, data_dir):
    RelativeTransmissionNoSymptoms = simulation_parameter["RelativeTransmissionNoSymptoms"]
    RiskOfInfectionFromSymptomatic = simulation_parameter["RiskOfInfectionFromSymptomatic"]
    scale_confirmed_cases = simulation_parameter["scale_confirmed_cases"]
    scale_contacts = simulation_parameter["scale_contacts"]
    lockdown_hard = simulation_parameter["lockdown_hard"]
    subprocess.call([f"./../build/bin/lct_realistic_scenario", f"{start_date.year}", f"{start_date.month}", f"{start_date.day}", data_dir,
                    f"{RelativeTransmissionNoSymptoms}", f"{RiskOfInfectionFromSymptomatic}",
                     f"{scale_confirmed_cases}", f"{scale_contacts}", f"{lockdown_hard}"])


def get_file_name(start_date, subcompartment, data_dir, boolagedistributed=False, boolsubcomp=False):
    filename = "real_" + start_date+"_" + f"{subcompartment}"+"_"
    if boolagedistributed:
        filename += "ageres"
    else:
        filename += "accumulated"
    if boolsubcomp:
        filename += "_subcompartments"
    return os.path.join(data_dir, filename)


def main_october():
    folder = "../data/simulation_lct_real/"
    datafile = "../data/pydata/Germany/cases_all_age.json"
    start_date = '2020-10-1'
    start_date_timestamp = pd.Timestamp(start_date)

    cases = [0, 1, 2]
    for case in cases:

        if case == 0:
            simulation_parameter = {"RelativeTransmissionNoSymptoms": 1.,
                                    "RiskOfInfectionFromSymptomatic": 0.3,
                                    "scale_confirmed_cases": 1.,
                                    "scale_contacts": 5747.075/8799.08555,
                                    "lockdown_hard": 0.3
                                    }
            run_simulation(start_date_timestamp,
                           simulation_parameter, "../data")
        if case == 1:
            plot_new_infections_real([get_file_name(start_date, 1, folder, False),
                                      get_file_name(
                                          start_date, 3, folder, False),
                                      get_file_name(
                                          start_date, 10, folder, False),
                                      get_file_name(
                                          start_date, 50, folder, False),
                                      get_file_name(
                                          start_date, 0, folder, False)],
                                     -1, datafile, start_date_timestamp, 45, 1.0,
                                     legendplot=list(
                ["Extrapolated RKI Data", "ODE", "LCT3", "LCT10", "LCT50", "LCTvar"]),
                filename_plot="real_new_infections_"+start_date+"_allage")
        if case == 2:
            for age in range(6):
                plot_new_infections_real([get_file_name(start_date, 1, folder, True),
                                          get_file_name(
                                              start_date, 3, folder, True),
                                          get_file_name(
                                              start_date, 10, folder, True),
                                          get_file_name(
                                              start_date, 50, folder, True),
                                          get_file_name(start_date, 0, folder, True)],
                                         age, datafile, start_date_timestamp, 45, 1.0,
                                         legendplot=list(
                    ["Extrapolated RKI Data", "ODE", "LCT3", "LCT10", "LCT50", "LCTvar"]),
                    filename_plot="real_new_infections_"+start_date+"_age"+f"{age}")


def main_july():
    folder = "../data/simulation_lct_real/"
    datafile = "../data/pydata/Germany/cases_all_age.json"
    start_date = '2020-7-1'
    start_date_timestamp = pd.Timestamp(start_date)

    cases = [0, 1, 2]
    for case in cases:

        if case == 0:
            simulation_parameter = {"RelativeTransmissionNoSymptoms": 1.,
                                    "RiskOfInfectionFromSymptomatic": 0.3,
                                    "scale_confirmed_cases": 1.,
                                    "scale_contacts": 492.3667 / 1200.0,
                                    "lockdown_hard": -0.3
                                    }
            run_simulation(start_date_timestamp,
                           simulation_parameter, "../data")
        if case == 1:
            plot_new_infections_real([get_file_name(start_date, 1, folder, False), get_file_name(
                start_date, 3, folder, False), get_file_name(
                start_date, 10, folder, False), get_file_name(
                start_date, 50, folder, False)],
                -1, datafile, start_date_timestamp, 45, 1.0,
                legendplot=list(
                ["Extrapolated RKI Data", "ODE", "LCT3", "LCT10", "LCT50"]),
                filename_plot="real_new_infections_"+start_date+"_allage")
        if case == 2:
            for age in range(6):
                plot_new_infections_real([get_file_name(start_date, 1, folder, True),
                                          get_file_name(
                                              start_date, 3, folder, True),
                                          get_file_name(
                                              start_date, 10, folder, True),
                                          get_file_name(
                                              start_date, 50, folder, True)
                                          ],
                                         age, datafile, start_date_timestamp, 45, 1.0,
                                         legendplot=list(
                    ["Extrapolated RKI Data", "ODE", "LCT3", "LCT10", "LCT50"]),
                    filename_plot="real_new_infections_"+start_date+"_age"+f"{age}")


if __name__ == "__main__":
    main_july()
