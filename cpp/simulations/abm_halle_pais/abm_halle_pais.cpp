/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Sascha Korf
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
/**
 * @file
 * @brief Halle ABM with PAIS, and the simulator side of a simulation based inference fit.
 *
 * Three modes:
 * - "simulate" runs the model once for one parameter vector and writes the observable time series.
 * - "ensemble" draws parameter vectors from the prior, simulates each, and writes one row per run of
 *   (parameters, observable). That table is the training set for the neural posterior estimator in
 *   fit_npe.py; the estimator itself is not implemented here, since it needs a neural density estimator.
 * - "priors" writes the prior bounds, so the Python side does not duplicate them.
 *
 * The ensemble mode distributes the runs over MPI ranks, which is where essentially all of the compute
 * of a fit sits.
 */
#include "halle_model.h"

#include "abm/simulation.h"
#include "memilio/io/cli.h"
#include "memilio/io/directories.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/miompi.h"

#include <fstream>
#include <iostream>
#include <numeric>

#ifdef MEMILIO_ENABLE_MPI
#include <mpi.h>
#endif

namespace
{

/// @brief Number of observable channels, see mio::halle::observable_channels().
constexpr size_t num_channels = 3;

/**
 * @brief Logger for the channels the fit and the PAIS projection need.
 *
 * Records cumulative deaths since the start of the simulation, and the number of Person%s currently in
 * each of the two PAIS severity states.
 */
struct LogChannels : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, std::array<double, num_channels>>;

    static Type log(const mio::abm::Simulation<>& sim)
    {
        std::array<double, num_channels> values{0.0, 0.0, 0.0};
        const auto t = sim.get_time();
        for (auto& person : sim.get_model().get_persons()) {
            if (person.get_infection_state(t) == mio::abm::InfectionState::Dead) {
                values[0] += 1.0;
            }
            switch (person.get_pais_state(t)) {
            case mio::abm::PAISState::Medium:
                values[1] += 1.0;
                break;
            case mio::abm::PAISState::Severe:
                values[2] += 1.0;
                break;
            default:
                break;
            }
        }
        return {t, values};
    }
};

/// @brief One simulation result: the daily value of every observable channel.
struct Observable {
    std::vector<double> days{}; ///< Day offset from the start of the simulation.
    std::vector<std::array<double, num_channels>> values{}; ///< One entry per day.

    /// @brief Flatten into the vector the fit consumes: all days of channel 0, then all days of channel 1, ...
    std::vector<double> flatten() const
    {
        std::vector<double> flat;
        flat.reserve(values.size() * num_channels);
        for (size_t channel = 0; channel < num_channels; ++channel) {
            for (const auto& value : values) {
                flat.push_back(value[channel]);
            }
        }
        return flat;
    }
};

/**
 * @brief Run the model once and return the observable.
 * @param[in] setup The non-fitted parts of the model.
 * @param[in] theta The fitted parameters.
 * @param[in] start_date Calendar date of the start of the simulation.
 * @param[in] num_days Length of the simulation in days.
 * @param[in] seed Seed of this run, so that a run is reproducible.
 */
mio::IOResult<Observable> run_once(const mio::halle::ModelSetup& setup, const std::vector<double>& theta,
                                   mio::Date start_date, int num_days, uint64_t seed)
{
    const auto t0   = mio::abm::TimePoint(0);
    const auto tmax = t0 + mio::abm::days(num_days);

    auto rng = mio::RandomNumberGenerator();
    rng.seed({static_cast<uint32_t>(seed & 0xFFFFFFFF), static_cast<uint32_t>(seed >> 32)});

    BOOST_OUTCOME_TRY(auto&& model, mio::halle::make_model(setup, theta, start_date, t0, tmax, rng));
    auto sim = mio::abm::Simulation(t0, std::move(model));

    mio::History<mio::DataWriterToMemory, LogChannels> history;
    sim.advance(tmax, history);

    Observable observable;
    // The simulation logs at every internal step, but the fit compares daily values, so only the last
    // entry of each day is kept.
    const auto& log = std::get<0>(history.get_log());
    for (size_t i = 0; i < log.size(); ++i) {
        const double day = log[i].first.days();
        if (i + 1 == log.size() || static_cast<int>(log[i + 1].first.days()) > static_cast<int>(day)) {
            observable.days.push_back(std::floor(day));
            observable.values.push_back(log[i].second);
        }
    }

    // Rebase the deaths channel onto the start of the simulation. The seeded infection history kills a
    // part of the population before t0, so the raw count is dominated by deaths that already happened;
    // the fit compares deaths that occur during the simulated window. The real data must be rebased the
    // same way, which fit_npe.py does. PAIS is not rebased: its level at t0 is a meaningful initial state.
    if (!observable.values.empty()) {
        const double deaths_at_t0 = observable.values.front()[0];
        for (auto& value : observable.values) {
            value[0] -= deaths_at_t0;
        }
    }
    return mio::success(observable);
}

/// @brief Write one observable as a CSV with a Day column and one column per channel.
mio::IOResult<void> write_observable(const std::string& path, const Observable& observable)
{
    std::ofstream file(path);
    if (!file.is_open()) {
        return mio::failure(mio::StatusCode::FileNotFound, "Could not open " + path + " for writing.");
    }
    file << "day";
    for (const auto& name : mio::halle::observable_channels()) {
        file << ',' << name;
    }
    file << '\n';
    for (size_t i = 0; i < observable.days.size(); ++i) {
        file << observable.days[i];
        for (size_t channel = 0; channel < num_channels; ++channel) {
            file << ',' << observable.values[i][channel];
        }
        file << '\n';
    }
    return mio::success();
}

/**
 * @brief Write the design matrix of an ensemble: one row per run, holding its parameters and observable.
 *
 * The header names every column, so that fit_npe.py does not need to know the parameter or channel
 * layout, and adding a fitted parameter needs no change on the Python side.
 */
mio::IOResult<void> write_ensemble(const std::string& path, const std::vector<std::vector<double>>& thetas,
                                   const std::vector<std::vector<double>>& observables, int num_days)
{
    std::ofstream file(path);
    if (!file.is_open()) {
        return mio::failure(mio::StatusCode::FileNotFound, "Could not open " + path + " for writing.");
    }
    for (const auto& parameter : mio::halle::fit_parameters()) {
        file << parameter.name << ',';
    }
    const auto& channels = mio::halle::observable_channels();
    for (size_t channel = 0; channel < channels.size(); ++channel) {
        for (int day = 0; day <= num_days; ++day) {
            file << channels[channel] << '_' << day;
            if (channel + 1 < channels.size() || day < num_days) {
                file << ',';
            }
        }
    }
    file << '\n';
    for (size_t run = 0; run < thetas.size(); ++run) {
        for (const auto& value : thetas[run]) {
            file << value << ',';
        }
        for (size_t i = 0; i < observables[run].size(); ++i) {
            file << observables[run][i];
            if (i + 1 < observables[run].size()) {
                file << ',';
            }
        }
        file << '\n';
    }
    return mio::success();
}

} // namespace

int main(int argc, char** argv)
{
    mio::mpi::init();

    auto parameters =
        mio::cli::ParameterSetBuilder()
            .add<"mode", std::string>("simulate", {"m", "One of \"simulate\" (a single run), \"ensemble\" (draw runs from "
                                                        "the prior for the fit) or \"priors\" (write the prior bounds)."})
            .add<"population_file", std::string>(mio::path_join(mio::data_dir().string(),
                                                                "Germany/halle_population_data.csv"),
                                                 {.description = "CSV holding the Halle population."})
            .add<"contact_dir", std::string>(mio::path_join(mio::data_dir().string(), "Germany/contacts"),
                                             {.description = "Directory with the German baseline contact matrices."})
            .add<"cases_file", std::string>("", {.description = "CSV holding the reported cases, for the "
                                                                "infection history. See prepare_data.py."})
            .add<"vaccinations_file", std::string>("", {.description = "CSV holding the vaccinations, for the "
                                                                       "vaccination history. See prepare_data.py."})
            .add<"start_date", std::string>("2022-07-01", {.description = "First day of the simulation."})
            .add<"num_days">(90, {"d", "Length of the simulation in days."})
            .add<"history_lookback_days">(90, {.description = "Length of the window before the start date that "
                                                               "the histories are seeded from."})
            .add<"allow_missing_history">(false, {.description = "Run without infection and vaccination history. "
                                                                  "Only for smoke tests, never for a fit."})
            .add<"num_runs">(1, {"n", "Number of runs. In ensemble mode, this many parameter vectors are drawn "
                                      "from the prior, distributed over the MPI ranks."})
            .add<"seed">(0, {"s", "Base seed. Every run derives its own seed from this and its index."})
            .add<"theta", std::vector<double>>({}, {.description = "Parameter vector for simulate mode. Empty means "
                                                                   "the centre of the prior."})
            .add<"output_dir", std::string>(mio::path_join(mio::base_dir().string(), "halle_pais_results"),
                                            {"o", "Directory to write results to."})
            .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, parameters);
    if (!cli_result) {
        std::cout << cli_result.error().message();
        mio::mpi::finalize();
        return cli_result.error().code().value();
    }

    mio::set_log_level(mio::LogLevel::warn);

    int rank = 0, num_procs = 1;
#ifdef MEMILIO_ENABLE_MPI
    MPI_Comm_rank(mio::mpi::get_world(), &rank);
    MPI_Comm_size(mio::mpi::get_world(), &num_procs);
#endif

    mio::halle::ModelSetup setup;
    setup.person_file           = parameters.get<"population_file">();
    setup.contact_dir           = parameters.get<"contact_dir">();
    setup.cases_file            = parameters.get<"cases_file">();
    setup.vaccinations_file     = parameters.get<"vaccinations_file">();
    setup.history_lookback_days = parameters.get<"history_lookback_days">();
    setup.allow_missing_history = parameters.get<"allow_missing_history">();

    const auto date_result = mio::parse_date(parameters.get<"start_date">());
    if (!date_result) {
        if (rank == 0) {
            std::cout << "Could not parse start_date: " << date_result.error().formatted_message() << "\n";
        }
        mio::mpi::finalize();
        return 1;
    }
    const auto start_date = date_result.value();
    const int num_days    = parameters.get<"num_days">();
    const auto output_dir = parameters.get<"output_dir">();
    const auto mode       = parameters.get<"mode">();
    const uint64_t seed   = static_cast<uint64_t>(parameters.get<"seed">());

    if (rank == 0 && !mio::create_directory(output_dir)) {
        std::cout << "Could not create output directory " << output_dir << "\n";
        mio::mpi::finalize();
        return 1;
    }

    if (mode == "simulate") {
        auto theta = parameters.get<"theta">();
        if (theta.empty()) {
            for (const auto& parameter : mio::halle::fit_parameters()) {
                theta.push_back(0.5 * (parameter.lower + parameter.upper));
            }
        }
        auto result = run_once(setup, theta, start_date, num_days, seed);
        if (!result) {
            std::cout << result.error().formatted_message() << "\n";
            mio::mpi::finalize();
            return 1;
        }
        if (rank == 0) {
            auto written = write_observable(mio::path_join(output_dir, "observable.csv"), result.value());
            if (!written) {
                std::cout << written.error().formatted_message() << "\n";
                mio::mpi::finalize();
                return 1;
            }
            std::cout << "Wrote " << mio::path_join(output_dir, "observable.csv") << "\n";
            const auto& last = result.value().values.back();
            std::cout << "Final cumulative deaths: " << last[0] << ", PAIS medium: " << last[1]
                      << ", PAIS severe: " << last[2] << "\n";
        }
    }
    else if (mode == "ensemble") {
        const int num_runs = parameters.get<"num_runs">();
        // Every rank takes the runs whose index is congruent to its own rank, so no communication is
        // needed to hand out work and the run index alone determines the seed and the parameter vector.
        std::vector<std::vector<double>> thetas, observables;
        for (int run = rank; run < num_runs; run += num_procs) {
            // Draw this run's parameters from a generator seeded by the run index, so the drawn vector
            // does not depend on how many ranks the ensemble is spread over.
            auto prior_rng = mio::RandomNumberGenerator();
            prior_rng.seed({static_cast<uint32_t>(seed), static_cast<uint32_t>(run)});
            const auto theta = mio::halle::sample_prior(prior_rng);

            auto result = run_once(setup, theta, start_date, num_days, seed + 1000000ULL + static_cast<uint64_t>(run));
            if (!result) {
                std::cout << "Rank " << rank << " run " << run << " failed: "
                          << result.error().formatted_message() << "\n";
                continue;
            }
            thetas.push_back(theta);
            observables.push_back(result.value().flatten());
            std::cout << "Rank " << rank << " finished run " << run << " of " << num_runs << "\n" << std::flush;
        }

        // Each rank writes its own shard. Concatenating shards is cheaper and less fragile than gathering
        // large observable vectors through MPI, and fit_npe.py reads a directory of shards.
        const auto shard = mio::path_join(output_dir, "ensemble_rank" + std::to_string(rank) + ".csv");
        auto written     = write_ensemble(shard, thetas, observables, num_days);
        if (!written) {
            std::cout << written.error().formatted_message() << "\n";
            mio::mpi::finalize();
            return 1;
        }
        std::cout << "Rank " << rank << " wrote " << thetas.size() << " runs to " << shard << "\n";
    }
    else if (mode == "priors") {
        // Emit the prior so that fit_npe.py does not have to duplicate the bounds. fit_parameters() stays
        // the single place where the fitted parameters and their priors are defined.
        if (rank == 0) {
            const auto path = mio::path_join(output_dir, "priors.csv");
            std::ofstream file(path);
            if (!file.is_open()) {
                std::cout << "Could not open " << path << " for writing.\n";
                mio::mpi::finalize();
                return 1;
            }
            file << "name,lower,upper\n";
            for (const auto& parameter : mio::halle::fit_parameters()) {
                file << parameter.name << ',' << parameter.lower << ',' << parameter.upper << '\n';
            }
            std::cout << "Wrote " << path << "\n";
        }
    }
    else {
        if (rank == 0) {
            std::cout << "Unknown mode \"" << mode << "\". Use \"simulate\", \"ensemble\" or \"priors\".\n";
        }
        mio::mpi::finalize();
        return 1;
    }

    mio::mpi::finalize();
    return 0;
}
