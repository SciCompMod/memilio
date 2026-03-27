/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Carlotta Gerstein
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

#include "memilio/compartments/simulation.h"
#include "memilio/utils/custom_index_array.h"
#include "models/ode_seir_metapop/model.h"
#include "models/ode_seir_metapop/parameters.h"
#include "memilio/geography/regions.h"
#include "memilio/data/analyze_result.h"

int main()
{
    const ScalarType t0   = 0.;
    const ScalarType tmax = 10;
    ScalarType dt         = 0.1;

    using namespace mio::oseirmetapop;

    Model<ScalarType> model(3, 1);

    for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
        model.populations[{mio::regions::Region(i), mio::AgeGroup(0), InfectionState::Susceptible}] = 10000;
    }

    model.populations[{mio::regions::Region(0), mio::AgeGroup(0), InfectionState::Infected}] += 100;
    model.populations[{mio::regions::Region(0), mio::AgeGroup(0), InfectionState::Susceptible}] -= 100;

    Eigen::MatrixXd mobility_data_commuter(3, 3);
    mobility_data_commuter << 0.4, 0.3, 0.3, 0.2, 0.7, 0.1, 0.4, 0.1, 0.5;

    model.set_commuting_strengths(mobility_data_commuter);

    model.parameters.template get<ContactPatterns<>>().get_cont_freq_mat()[0].get_baseline().setConstant(2.7);

    model.parameters.get<mio::oseirmetapop::TimeExposed<>>()[{mio::regions::Region(0), mio::AgeGroup(0)}]  = 3.;
    model.parameters.get<mio::oseirmetapop::TimeExposed<>>()[{mio::regions::Region(1), mio::AgeGroup(0)}]  = 4.;
    model.parameters.get<mio::oseirmetapop::TimeExposed<>>()[{mio::regions::Region(2), mio::AgeGroup(0)}]  = 5.;
    model.parameters.get<mio::oseirmetapop::TimeInfected<>>()[{mio::regions::Region(0), mio::AgeGroup(0)}] = 7.;
    model.parameters.get<mio::oseirmetapop::TimeInfected<>>()[{mio::regions::Region(1), mio::AgeGroup(0)}] = 8.;
    model.parameters.get<mio::oseirmetapop::TimeInfected<>>()[{mio::regions::Region(2), mio::AgeGroup(0)}] = 9.;
    model.parameters.set<TransmissionProbabilityOnContact<>>(0.07333);

    auto result              = simulate(t0, tmax, dt, model);
    auto interpolated_result = mio::interpolate_simulation_result(result);

    std::vector<std::string> vars = {"S", "E", "I", "R"};
    printf("Infected individuals per Region over time [%%]:\n");
    for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
        printf("\t Region %ld", i);
    }
    for (size_t t : interpolated_result.get_times()) {
        printf("\n %ld", t);
        // for (size_t k = 0; k < (size_t)mio::oseir::InfectionState::Count; k++) {
        for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {

            printf("\t %.5f ",
                   interpolated_result.get_value(t)[(size_t)model.parameters.get_num_regions() *
                                                        ((size_t)mio::oseirmetapop::InfectionState::Count - 1) +
                                                    i] /
                       model.populations.get_group_total(mio::regions::Region(i)) * 100);
        }
    }
    return 0;
}
