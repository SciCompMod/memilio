/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/initializer_flows.h"

#include "memilio/config.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <iostream>
#include <vector>

namespace params
{
// num_subcompartments is used as a template argument and has to be a constexpr.
constexpr size_t num_subcompartments = NUM_SUBCOMPARTMENTS;

// Define (non-age-resolved) parameters.
const ScalarType dt                             = 0.01;
const ScalarType seasonality                    = 0.;
const ScalarType relativeTransmissionNoSymptoms = 1.;
const ScalarType riskOfInfectionFromSymptomatic = 0.3;
const ScalarType total_population               = 83155031.0;

ScalarType timeExposed                            = 3.335;
const ScalarType timeInfectedNoSymptoms           = 2.58916;
const ScalarType timeInfectedSymptoms             = 6.94547;
const ScalarType timeInfectedSevere               = 7.28196;
const ScalarType timeInfectedCritical             = 13.066;
const ScalarType transmissionProbabilityOnContact = 0.07333;
const ScalarType recoveredPerInfectedNoSymptoms   = 0.206901;
const ScalarType severePerInfectedSymptoms        = 0.07864;
const ScalarType criticalPerSevere                = 0.17318;
const ScalarType deathsPerCritical                = 0.21718;
} // namespace params

/** 
* @brief Constructs an initial value vector such that the initial infection dynamic is constant.
*   
*   The initial value vector is constructed based on a value of approximately 4050 for the daily new transmissions. 
*   This value is based on some official reporting numbers for Germany. The vector is constructed such that the 
*   daily new transmissions remain constant if the effective reproduction number is equal to 1. 
*   Derived numbers for compartments are distributed equally to the subcompartments.
*   To define the vector, we assume that individuals behave exactly as defined by the epidemiological parameters 
*   and that the new transmissions are constant over time. 
*   The value for the Recovered compartment is also set according to the reported numbers 
*   and the number of death is set to zero.
*   
* @returns The initial value vector with subcompartments constructed.
*/
std::vector<ScalarType> get_initial_values()
{
    using namespace params;
    using InfState                                 = mio::lsecir::InfectionState;
    const ScalarType dailyNewTransmissionsReported = (34.1 / 7) * total_population / 100000;

    // Firstly, we calculate an initial value vector without division in subcompartments.
    std::vector<ScalarType> init_compartments((size_t)InfState::Count);

    init_compartments[(size_t)InfState::Exposed]            = timeExposed * dailyNewTransmissionsReported;
    init_compartments[(size_t)InfState::InfectedNoSymptoms] = timeInfectedNoSymptoms * dailyNewTransmissionsReported;
    init_compartments[(size_t)InfState::InfectedSymptoms] =
        (1 - recoveredPerInfectedNoSymptoms) * timeInfectedSymptoms * dailyNewTransmissionsReported;
    init_compartments[(size_t)InfState::InfectedSevere] = (1 - recoveredPerInfectedNoSymptoms) *
                                                          severePerInfectedSymptoms * timeInfectedSevere *
                                                          dailyNewTransmissionsReported;
    init_compartments[(size_t)InfState::InfectedCritical] = (1 - recoveredPerInfectedNoSymptoms) *
                                                            severePerInfectedSymptoms * criticalPerSevere *
                                                            timeInfectedCritical * dailyNewTransmissionsReported;
    init_compartments[(size_t)InfState::Recovered]   = 275292.;
    init_compartments[(size_t)InfState::Dead]        = 0.;
    init_compartments[(size_t)InfState::Susceptible] = total_population;
    for (size_t i = (size_t)InfState::Exposed; i < (size_t)InfState::Count; i++) {
        init_compartments[(size_t)InfState::Susceptible] -= init_compartments[i];
    }

    // Now, we construct an initial value vector with division in subcompartments.
    // Compartment sizes are distributed equally to the subcompartments.
    std::vector<ScalarType> initial_value_vector;
    initial_value_vector.push_back(init_compartments[(size_t)InfState::Susceptible]);
    for (size_t i = 0; i < num_subcompartments; i++) {
        initial_value_vector.push_back(init_compartments[(size_t)InfState::Exposed] / num_subcompartments);
    }
    for (size_t i = 0; i < num_subcompartments; i++) {
        initial_value_vector.push_back(init_compartments[(size_t)InfState::InfectedNoSymptoms] / num_subcompartments);
    }
    for (size_t i = 0; i < num_subcompartments; i++) {
        initial_value_vector.push_back(init_compartments[(size_t)InfState::InfectedSymptoms] / num_subcompartments);
    }
    for (size_t i = 0; i < num_subcompartments; i++) {
        initial_value_vector.push_back(init_compartments[(size_t)InfState::InfectedSevere] / num_subcompartments);
    }
    for (size_t i = 0; i < num_subcompartments; i++) {
        initial_value_vector.push_back(init_compartments[(size_t)InfState::InfectedCritical] / num_subcompartments);
    }
    initial_value_vector.push_back(init_compartments[(size_t)InfState::Recovered]);
    initial_value_vector.push_back(init_compartments[(size_t)InfState::Dead]);

    return initial_value_vector;
}

/** 
* @brief Perform simulation to examine the impact of the distribution assumption. 
*
*   The simulation uses LCT models with Covid-19 inspired parameters and an initial contact rate such that the 
*   effective reproduction number is initially approximately equal to one. The initial values are chosen such that the
*   daily new transmissions remain constant. The contact rate is changed at simulation day 2 such that the effective 
*   reproduction number at simulation day 2 is equal to the input Reff2. 
*   Therefore, we simulate a change point at day 2.
*   
* @param[in] Reff2 The effective reproduction number to be set at simulation time 2. Please use a number greater zero.
* @param[in] tmax Time horizon of the simulation.
* @param[in] save_dir Specifies the directory where the results should be stored.
*    Provide an empty string if results should not be saved.
* @param[in] save_subcompartments If true, the result will be saved with division in subcompartments. Default is false.
* @param[in] scale_TimeExposed The value for TimeExposed (=3.335) is scaled by this value before the simulation.
*   Default is 1 (no scaling).
* @param[in] print_final_size If true, the final size will be printed. Default is false. 
*   The printed values refer to the final size only if tmax is chosen large enough.
* @returns Any IO errors that occur when saving the results.
*/
mio::IOResult<void> simulate(ScalarType Reff2, ScalarType tmax, std::string save_dir = "",
                             bool save_subcompartments = false, ScalarType scale_TimeExposed = 1.,
                             bool print_final_size = false)
{
    using namespace params;
    std::cout << "Simulation with " << num_subcompartments << " subcompartments and reproduction number " << Reff2
              << "." << std::endl;

    // Initialize LCT model.
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Set parameters.
    // Scale TimeExposed for some numerical experiments.
    timeExposed                                                              = scale_TimeExposed * timeExposed;
    model.parameters.get<mio::lsecir::TimeExposed>()[0]                      = timeExposed;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0]           = timeInfectedNoSymptoms;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]             = timeInfectedSymptoms;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]               = timeInfectedSevere;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]             = timeInfectedCritical;
    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = transmissionProbabilityOnContact;

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = relativeTransmissionNoSymptoms;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = riskOfInfectionFromSymptomatic;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = recoveredPerInfectedNoSymptoms;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = severePerInfectedSymptoms;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = criticalPerSevere;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = deathsPerCritical;

    // Determine the contact rate such that the effective reproduction number is initially approximately equal to one.
    ScalarType contacts_R1 =
        1. / (transmissionProbabilityOnContact *
              (timeInfectedNoSymptoms * relativeTransmissionNoSymptoms +
               (1 - recoveredPerInfectedNoSymptoms) * timeInfectedSymptoms * riskOfInfectionFromSymptomatic));

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    if (Reff2 <= 1.) {
        // Perform a simulation with a decrease in the effective reproduction number on day 2.
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, contacts_R1));
        contact_matrix[0].add_damping(0., mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(Reff2, mio::SimulationTime(2.));
    }
    else {
        // Perform a simulation with an increase in the effective reproduction number on day 2.
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, Reff2 * contacts_R1));
        contact_matrix[0].add_damping(1 - 1. / Reff2, mio::SimulationTime(-1.));
        contact_matrix[0].add_damping(1 - 1. / Reff2, mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(0., mio::SimulationTime(2.));
    }

    model.parameters.get<mio::lsecir::ContactPatterns>() = mio::UncertainContactMatrix<ScalarType>(contact_matrix);
    model.parameters.get<mio::lsecir::Seasonality>()     = seasonality;

    // Set initial values.
    std::vector<ScalarType> initial_values = get_initial_values();
    for (size_t i = 0; i < model.populations.get_num_compartments(); i++) {
        model.populations[i] = initial_values[i];
    }

    // Set integrator of fifth order with fixed step size and perform simulation.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max to get a fixed step size.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);

    // Save results and print desired information.
    if (!save_dir.empty()) {
        std::string Reff2string = std::to_string(Reff2);
        std::string filename = save_dir + "lct_Reff" + Reff2string.substr(0, Reff2string.find(".") + 2) + "_subcomp" +
                               std::to_string(num_subcompartments);
        if (save_subcompartments) {
            filename = filename + "_subcompartments.h5";
            // Just store daily results in this case.
            auto result_interpolated   = mio::interpolate_simulation_result(result, dt / 2);
            mio::IOResult<void> status = mio::save_result({result_interpolated}, {0}, 1, filename);
            if (!status) {
                return status;
            }
        }
        else {
            // Calculate result without division in subcompartments.
            mio::TimeSeries<ScalarType> populations = model.calculate_compartments(result);
            filename                                = filename + ".h5";
            mio::IOResult<void> status              = mio::save_result({populations}, {0}, 1, filename);
            if (!status) {
                return status;
            }
        }
    }
    if (print_final_size) {
        std::cout << "Final size: " << std::fixed << std::setprecision(6)
                  << total_population - result.get_last_value()[0] << std::endl;
        std::cout << std::endl;
    }
    return mio::success();
}

/** 
* Usage: lct_impact_distribution_assumption <Reff2> <tmax> <save_dir> <save_subcompartments> <scale_TimeExposed> 
*           <print_final_size>
*   All command line arguments are optional. Simple default values are provided if not specified.
*   All parameters are passed to the simulation() function. See the documentation for a description of the parameters.
*
*   The numbers of subcompartments used in the LCT model is determined by the preprocessor macro NUM_SUBCOMPARTMENTS.
*   You can set the number via the flag -DNUM_SUBCOMPARTMENTS=... . 
*/
int main(int argc, char** argv)
{
    ScalarType Reff2             = 1.;
    ScalarType tmax              = 40;
    std::string save_dir         = "";
    bool save_subcompartments    = false;
    ScalarType scale_TimeExposed = 1.;
    bool print_final_size        = false;

    switch (argc) {
    case 7:
        print_final_size = std::stoi(argv[6]);
        [[fallthrough]];
    case 6:
        scale_TimeExposed = std::stod(argv[5]);
        [[fallthrough]];
    case 5:
        save_subcompartments = std::stoi(argv[4]);
        [[fallthrough]];
    case 4:
        save_dir = argv[3];
        [[fallthrough]];
    case 3:
        tmax = std::stod(argv[2]);
        [[fallthrough]];
    case 2:
        Reff2 = std::stod(argv[1]);
    }
    auto result = simulate(Reff2, tmax, save_dir, save_subcompartments, scale_TimeExposed, print_final_size);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
