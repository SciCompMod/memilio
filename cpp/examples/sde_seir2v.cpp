/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Nils Wassmuth, Rene Schmieding, Martin J. Kuehn
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
#include <fstream>
#include <vector>
#include <iostream>

#include "memilio/utils/logging.h"
#include "sde_seir2v/model.h"
#include "sde_seir2v/simulation.h"
#include "memilio/utils/random_number_generator.h"

template <typename T>
void print_to_file(const T& history, const std::vector<std::string>& column_labels = {}, const std::string filepath = "seir2v_output.txt")
{
    std::ofstream myfile(filepath);
    for (size_t k = 0; k < static_cast<size_t>(history.get_num_elements()); k++) {
        if (k < column_labels.size()) {
            myfile << " ";
            myfile <<  column_labels[k];
        }
        else {
            myfile << " ";
            myfile << "#" + std::to_string(k + 1);
        }
    }
    // print values as table
    auto num_points = static_cast<size_t>(history.get_num_time_points());
    for (size_t i = 0; i < num_points; i++) {
        myfile << "\n";
        myfile << history.get_time(i);
        auto res_i = history.get_value(i);
        for (size_t j = 0; j < static_cast<size_t>(res_i.size()); j++) {
            myfile << " ";
            myfile << res_i[j];
        }
    }
    myfile << "\n";
    myfile.close();
}

/*template <typename T>
void print_to_file2(const T& history_pre, const T& history_post, const double tmid,  const std::vector<std::string>& column_labels = {}, const std::string filepath = "seir2v_output2.txt")
{   
    
    std::vector<double> N;
    std::vector<double> S_prob;
    std::vector<double> E_prob;
    std::vector<double> I_prob;
    std::vector<double> R_prob;

    std::vector<double> S;

    //Set timepoints
    std::vector<size_t> timepoints {20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380}
    auto num_points = static_cast<size_t>(history.get_num_time_points());

    for (size_t it = 0; it < timepoints.size(); ++it)
    {   
        if (timepoints[it] < tmid) 
        {
            S.append(history_pre.get_value(timepoints[it])[0])
        }
    }
    std::ofstream myfile(filepath);
    for (size_t it = 0; it < S.size(); ++it){
        myfile << " ";
        myfile << S[it];
    }
    myfile << "\n";
    myfile.close();
}*/

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmid = 100.;
    double tmax = 400.;
    double dt   = 0.1;

    double total_population = 1000;

    mio::RandomNumberGenerator rng;
    std::initializer_list<uint32_t> seeds = {14159265u, 35897932u};
    rng.seed(seeds);

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::sseir2v::Model model;

    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::ExposedV1)}]  = 0;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::ExposedV2)}]  = 0;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::InfectedV1)}]  = 100;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::InfectedV2)}]  = 0.01;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::RecoveredV1)}] = 0;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::RecoveredV2)}] = 0;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::ExposedV1V2)}]  = 0;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::InfectedV1V2)}]  = 0;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::RecoveredV1V2)}] = 0;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::ExposedV1)}] -
        model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::ExposedV2)}] -
        model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::InfectedV1)}] -
        model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::InfectedV2)}] -
        model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::RecoveredV1)}] -
        model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::RecoveredV2)}] -
        model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::ExposedV1V2)}] -
        model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::InfectedV1V2)}] -
        model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::RecoveredV1V2)}];

    double beta = mio::DistributionAdapter<std::lognormal_distribution<double>>::get_instance()(rng, std::log(0.08), 0.5);
    double kappa = mio::DistributionAdapter<std::lognormal_distribution<double>>::get_instance()(rng, std::log(8), 0.2);
    double gamma = mio::DistributionAdapter<std::lognormal_distribution<double>>::get_instance()(rng, std::log(16), 0.2);

    model.parameters.set<mio::sseir2v::TransmissionProbabilityOnContactV1>(beta);
    model.parameters.set<mio::sseir2v::TransmissionProbabilityOnContactV2>(beta);
    model.parameters.set<mio::sseir2v::TimeExposedV1>(kappa);
    model.parameters.set<mio::sseir2v::TimeExposedV2>(kappa);      
    model.parameters.set<mio::sseir2v::TimeInfectedV1>(gamma);
    model.parameters.set<mio::sseir2v::TimeInfectedV2>(gamma);

    model.check_constraints();

    auto ssirs = mio::sseir2v::simulate(t0, 23, dt, model);
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::Susceptible)}] = ssirs.get_value(static_cast<size_t>(ssirs.get_num_time_points()) - 1)[0];
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::ExposedV1)}] = ssirs.get_value(static_cast<size_t>(ssirs.get_num_time_points()) - 1)[1];
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::InfectedV1)}]  = ssirs.get_value(static_cast<size_t>(ssirs.get_num_time_points()) - 1)[2];
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::RecoveredV1)}] = ssirs.get_value(static_cast<size_t>(ssirs.get_num_time_points()) - 1)[3];    
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::ExposedV2)}]  = ssirs.get_value(static_cast<size_t>(ssirs.get_num_time_points()) - 1)[4];
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::InfectedV2)}]  = 10;
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::RecoveredV2)}] = ssirs.get_value(static_cast<size_t>(ssirs.get_num_time_points()) - 1)[6];
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::ExposedV1V2)}]  = ssirs.get_value(static_cast<size_t>(ssirs.get_num_time_points()) - 1)[7];
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::InfectedV1V2)}]  = ssirs.get_value(static_cast<size_t>(ssirs.get_num_time_points()) - 1)[8];
    model.populations[{mio::Index<mio::sseir2v::InfectionState>(mio::sseir2v::InfectionState::RecoveredV1V2)}] = ssirs.get_value(static_cast<size_t>(ssirs.get_num_time_points()) - 1)[9];
    auto ssirs2 = mio::sseir2v::simulate(tmid, tmax, dt, model);

    //print_to_file(ssirs, {"Susceptible", "ExposedV1", "InfectedV1", "RecoveredV1", "ExposedV2", "InfectedV2", "RecoveredV2", "ExposedV1V2", "InfectedV1V2", "RecoveredV1V2"}, "seir2v_1.txt");
    ssirs.print_table({"Susceptible", "ExposedV1", "InfectedV1", "RecoveredV1", "ExposedV2", "InfectedV2", "RecoveredV2", "ExposedV1V2", "InfectedV1V2", "RecoveredV1V2"});
    //ssirs2.print_table({"Susceptible", "ExposedV1", "InfectedV1", "RecoveredV1", "ExposedV2", "InfectedV2", "RecoveredV2", "ExposedV1V2", "InfectedV1V2", "RecoveredV1V2"});
    getchar();
}
