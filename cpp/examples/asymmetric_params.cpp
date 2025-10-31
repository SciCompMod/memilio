/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Kilian Volmer
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
#include "memilio/config.h"
#include "memilio/geography/geolocation.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/geography/rtree.h"
#include "memilio/geography/distance.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_asymmetric.h"
#include "memilio/mobility/graph.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/logging.h"
#include "memilio/timer/auto_timer.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/base_dir.h"
#include "smm/simulation.h"
#include "smm/parameters.h"
#include "abm/time.h"
#include "thirdparty/csv.h"
#include <ranges>
#include <omp.h>

enum class InfectionState
{
    S,
    E,
    I,
    INS,
    ICS,
    R,
    D,
    Count
};

template <class T>
MPI_Datatype mpi_type();

template <>
inline MPI_Datatype mpi_type<int>()
{
    return MPI_INT;
}
template <>
inline MPI_Datatype mpi_type<double>()
{
    return MPI_DOUBLE;
}
template <>
inline MPI_Datatype mpi_type<float>()
{
    return MPI_FLOAT;
}
template <>
inline MPI_Datatype mpi_type<long long>()
{
    return MPI_LONG_LONG;
}
// If your MPI supports it (MPI-3+), for fixed-width:
#include <cstdint>
template <>
inline MPI_Datatype mpi_type<std::int64_t>()
{
    return MPI_INT64_T;
}
template <>
inline MPI_Datatype mpi_type<std::uint64_t>()
{
    return MPI_UINT64_T;
}

template <class T>
void bcast_vector(std::vector<T>& v, int root, MPI_Comm comm)
{
    int n = static_cast<int>(v.size());
    MPI_Bcast(&n, 1, MPI_INT, root, comm);
    if (n < 0)
        std::abort();
    if (n == 0) {
        v.clear();
        return;
    }
    if (v.size() != static_cast<size_t>(n))
        v.resize(n);
    MPI_Bcast(v.data(), n, mpi_type<T>(), root, comm);
}

int main(int /*argc*/, char** /*argv*/)
{
    mio::mpi::init();
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::string result_dir = mio::path_join(mio::base_dir(), "example_results");

    const auto t0   = 0.;
    const auto tmax = 100.;
    const auto dt   = 1.; //initial time step

    //total compartment sizes

    using Model = mio::smm::Model<ScalarType, 1, InfectionState>;
    auto home   = mio::regions::Region(0);
    auto S      = InfectionState::S;
    auto E      = InfectionState::E;
    auto I      = InfectionState::I;
    auto INS    = InfectionState::INS;
    auto ICS    = InfectionState::ICS;
    auto R      = InfectionState::R;
    auto D      = InfectionState::D;

    std::vector<mio::AdoptionRate<ScalarType, InfectionState>> adoption_rates;
    // Adoption rates corresponding to our model, paramters are arbitrary
    adoption_rates.push_back({S, E, home, 0.2, {{I, 0.8}, {INS, 0.1}, {ICS, 0.5}}});
    adoption_rates.push_back({S, E, home, 0.0, {}});
    adoption_rates.push_back({E, I, home, 0.2, {}});
    adoption_rates.push_back({I, INS, home, 0.1, {}});
    adoption_rates.push_back({I, ICS, home, 0.1, {}});
    adoption_rates.push_back({ICS, D, home, 0.6, {}});
    adoption_rates.push_back({ICS, R, home, 0.4, {}});
    adoption_rates.push_back({INS, R, home, 0.5, {}});

    mio::Graph<mio::LocationNode<ScalarType, mio::smm::Simulation<ScalarType, 1, InfectionState>>,
               mio::MobilityEdgeDirected<ScalarType>>
        graph;

    // usage:
    // if (rank==0) v = {1,2,3,4};
    // bcast_vector(v, 0, MPI_COMM_WORLD);
    // now every rank has v
    std::vector<ScalarType> latitudes, longitudes;
    std::vector<int> farm_ids, num_cows_vec, dates, num_animals_exchanges;
    std::vector<int> froms, tos, from_exchanges, to_exchanges;
    if (rank == 0) {
        io::CSVReader<4> farms("../../farms10000.csv");
        farms.read_header(io::ignore_extra_column, "farms", "num_cows", "latitude", "longitude");
        int farm_id, num_cows;
        double latitude, longitude;
        while (farms.read_row(farm_id, num_cows, latitude, longitude)) {
            farm_ids.push_back(farm_id);
            num_cows_vec.push_back(num_cows);
            latitudes.push_back(latitude);
            longitudes.push_back(longitude);
        }
        io::CSVReader<2> edges("../../edges10000.csv");
        edges.read_header(io::ignore_extra_column, "from", "to");
        size_t from, to;
        while (edges.read_row(from, to)) {
            froms.push_back(from);
            tos.push_back(to);
        }
        io::CSVReader<4> exchanges("../../trade10000.csv");
        exchanges.read_header(io::ignore_extra_column, "date", "num_animals", "from", "to");

        int date, num_animals;
        while (exchanges.read_row(date, num_animals, from, to)) {
            dates.push_back(date);
            num_animals_exchanges.push_back(num_animals);
            from_exchanges.push_back(from);
            to_exchanges.push_back(to);
        }
        if (!mio::create_directory(result_dir)) {
            mio::log_error("Could not create result directory \"{}\".", result_dir);
            return 1;
        }
    }
    bcast_vector(farm_ids, 0, MPI_COMM_WORLD);
    bcast_vector(num_cows_vec, 0, MPI_COMM_WORLD);
    bcast_vector(latitudes, 0, MPI_COMM_WORLD);
    bcast_vector(longitudes, 0, MPI_COMM_WORLD);
    bcast_vector(froms, 0, MPI_COMM_WORLD);
    bcast_vector(tos, 0, MPI_COMM_WORLD);
    bcast_vector(dates, 0, MPI_COMM_WORLD);
    bcast_vector(num_animals_exchanges, 0, MPI_COMM_WORLD);
    bcast_vector(from_exchanges, 0, MPI_COMM_WORLD);
    bcast_vector(to_exchanges, 0, MPI_COMM_WORLD);

    for (size_t i = 0; i < farm_ids.size(); ++i) {
        Model curr_model;
        curr_model.populations[{home, InfectionState::S}]                                = num_cows_vec[i];
        curr_model.populations[{home, InfectionState::E}]                                = 0;
        curr_model.populations[{home, InfectionState::I}]                                = 0;
        curr_model.populations[{home, InfectionState::INS}]                              = 0;
        curr_model.populations[{home, InfectionState::ICS}]                              = 0;
        curr_model.populations[{home, InfectionState::R}]                                = 0;
        curr_model.populations[{home, InfectionState::D}]                                = 0;
        curr_model.parameters.get<mio::smm::AdoptionRates<ScalarType, InfectionState>>() = adoption_rates;
        graph.add_node(farm_ids[i], longitudes[i], latitudes[i], curr_model, t0);
    }
    auto rng = mio::RandomNumberGenerator();

    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back({Model().populations.get_flat_index({home, InfectionState::I})});

    for (size_t i = 0; i < froms.size(); ++i) {
        graph.add_edge(froms[i], tos[i], interesting_indices);
    }

    auto nodes = graph.nodes() | std::views::transform([](const auto& node) {
                     return &node.property;
                 });
    auto tree  = mio::geo::RTree(nodes.begin(), nodes.end());

    for (auto& node : graph.nodes()) {
        node.property.set_regional_neighbors(
            tree.in_range_indices_query(node.property.get_location(), {mio::geo::kilometers(2.0)}));
    }

    auto sim = mio::make_mobility_sim(t0, dt, std::move(graph));

    for (size_t i = 0; i < dates.size(); ++i) {
        sim.add_exchange(dates[i], num_animals_exchanges[i], from_exchanges[i], to_exchanges[i]);
    }

    const size_t num_runs = 4000;

    mio::ParameterStudy study(sim, t0, tmax, dt, num_runs);

    auto get_simulation = [](const auto& sim, ScalarType, ScalarType, size_t run_id) {
        auto sim2      = sim;
        sim2.get_rng() = mio::thread_local_rng();
        if (run_id % 2 == 0) {
            sim2.infectionrisk = 0.1;
        }
        auto index = sim2.get_graph().nodes()[0].property.get_simulation().get_model().populations.get_flat_index(
            {mio::regions::Region(0), InfectionState::E});
        sim2.get_graph().nodes()[0].property.get_result().get_last_value()[index] = 10;
        return sim2;
    };
    auto handle_result = [](auto&& sim, auto&& run) {
        auto abs_path = mio::path_join(mio::base_dir(), "example_results");
        auto result   = sim.statistics_per_timestep().export_csv(
            mio::path_join(abs_path, "AsymmetricParams_run" + std::to_string(run) + ".csv"));
        return 0;
    };
    auto result = study.run(get_simulation, handle_result);
    mio::mpi::finalize();
    return 0;
}
