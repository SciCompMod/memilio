/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#include "ode_secir/model.h"
#include "ode_secir/parameters_io.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/graph.h"
#include "memilio/utils/logging.h"
#include <iostream>
#include <vector>
#include <utility>

int main()
{
    mio::set_log_level(mio::LogLevel::info);
    using FP = double;

    // params for 1 age group
    mio::osecir::Parameters<FP> params(1);
    params.get<mio::osecir::TimeExposed<FP>>()[mio::AgeGroup(0)]                      = 3.2;
    params.get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[mio::AgeGroup(0)]           = 2.0;
    params.get<mio::osecir::TimeInfectedSymptoms<FP>>()[mio::AgeGroup(0)]             = 5.8;
    params.get<mio::osecir::TimeInfectedSevere<FP>>()[mio::AgeGroup(0)]               = 9.5;
    params.get<mio::osecir::TimeInfectedCritical<FP>>()[mio::AgeGroup(0)]             = 7.1;
    params.get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[mio::AgeGroup(0)] = 0.05;
    params.get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[mio::AgeGroup(0)]   = 0.7;
    params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[mio::AgeGroup(0)]   = 0.09;
    params.get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[mio::AgeGroup(0)]        = 0.2;
    params.get<mio::osecir::CriticalPerSevere<FP>>()[mio::AgeGroup(0)]                = 0.25;
    params.get<mio::osecir::DeathsPerCritical<FP>>()[mio::AgeGroup(0)]                = 0.3;

    // input data
    const mio::Date date{2020, 12, 1};
    const auto& data_dir                   = "/localdata1/code_2025/memilio/data";
    const std::string pydata_dir           = mio::path_join(data_dir, "Germany", "pydata");
    const std::string population_data_path = mio::path_join(pydata_dir, "county_current_population.json");

    // scaling factors
    std::vector<double> scaling_factor_inf(static_cast<size_t>(params.get_num_groups()), 1.0);
    const double scaling_factor_icu  = 1.0;
    const double tnt_capacity_factor = 7.5 / 100000.0;

    // grraph
    mio::Graph<mio::osecir::Model<FP>, mio::MobilityParameters<FP>> graph;

    const auto& read_function_nodes = mio::osecir::read_input_data_county<mio::osecir::Model<FP>>;
    auto node_id_function           = [](const std::string&, bool, bool) -> mio::IOResult<std::vector<int>> {
        return mio::success(std::vector<int>{1002});
    };

    const auto& set_node_function =
        mio::set_nodes<mio::osecir::TestAndTraceCapacity<FP>, mio::osecir::ContactPatterns<FP>, mio::osecir::Model<FP>,
                       mio::MobilityParameters<FP>, mio::osecir::Parameters<FP>, decltype(read_function_nodes),
                       decltype(node_id_function), FP>;

    auto io = set_node_function(params, date, date, pydata_dir, population_data_path, /*is_county=*/true, graph,
                                read_function_nodes, std::move(node_id_function), scaling_factor_inf,
                                scaling_factor_icu, tnt_capacity_factor, /*num_days=*/0,
                                /*export_time_series=*/false, /*rki_age_groups=*/true);
    if (!io) {
        std::cerr << io.error().formatted_message() << std::endl;
        return 1;
    }
    // check output
    const auto& m = graph.nodes()[0].property;
    const auto ag = mio::AgeGroup(0);
    std::cout << "Initialized via set_nodes for county 1002 on " << date << "\n";
    std::cout << "S =" << m.populations[{ag, mio::osecir::InfectionState::Susceptible}].value() << ", ";
    std::cout << "E =" << m.populations[{ag, mio::osecir::InfectionState::Exposed}].value() << ", ";
    std::cout << "C =" << m.populations[{ag, mio::osecir::InfectionState::InfectedNoSymptoms}].value() << ", ";
    std::cout << "I =" << m.populations[{ag, mio::osecir::InfectionState::InfectedSymptoms}].value() << ", ";
    std::cout << "H =" << m.populations[{ag, mio::osecir::InfectionState::InfectedSevere}].value() << ", ";
    std::cout << "U =" << m.populations[{ag, mio::osecir::InfectionState::InfectedCritical}].value() << ", ";
    std::cout << "R =" << m.populations[{ag, mio::osecir::InfectionState::Recovered}].value() << ", ";
    std::cout << "D =" << m.populations[{ag, mio::osecir::InfectionState::Dead}].value() << std::endl;

    return 0;
}
