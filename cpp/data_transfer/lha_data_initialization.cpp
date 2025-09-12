/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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

#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/parameters_io.cpp"
#include <vector>

int main()
{
    // County ID. Use Cologne for now.
    std::vector<int> lha_ids = {5314};

    size_t num_age_groups = 6;

    // Set up Graph-SECIRVVS model.

    mio::Date current_date = mio::Date(2021, 2, 1);

    mio::osecirvvs::Parameters<double> params(num_age_groups);

    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> graph_model;

    const std::string data_dir = "../../data";

    mio::IOResult<void> result =
        mio::osecirvvs::set_lha_data<double>(params, graph_model, data_dir, current_date, lha_ids);

    return 0;
}