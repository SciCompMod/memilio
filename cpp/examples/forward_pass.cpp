/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Nils Wassmuth
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
#include "abm/model.h"
#include "abm/common_abm_loggers.h"
#include "abm/intervention_type.h"
#include "city_builder.h"
#include "forward_pass.h"

#include <array>
#include <cmath>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
// forward_pass() may be called concurrently from multiple GIL-released Python threads
// set_log_level sets global state, so it must only ever run once, not once per call.
std::once_flag g_log_level_once;

// Per-day probability a symptomatic, non-isolated Person seeks a diagnostic PCR test.
constexpr ScalarType kCareSeekPerDay = 0.19;
} // namespace

struct ABMPopulation::Impl {
    mio::abm::Model model;

    explicit Impl(int total_population)
        : model(CityBuilder::build_world(CityConfig{total_population}, mio::thread_local_rng()))
    {
        // Stochastic viral-load kinetics.
        model.parameters.get<mio::abm::ViralLoadDistributions>() = mio::abm::ViralLoadDistributionsParameters{
            mio::AbstractParameterDistribution(mio::ParameterDistributionUniform(7.0, 9.2)),      // peak (log10)
            mio::AbstractParameterDistribution(mio::ParameterDistributionUniform(1.4, 2.6)),      // incline
            mio::AbstractParameterDistribution(mio::ParameterDistributionUniform(-0.24, -0.10)),  // decline
        };
    }
};

ABMPopulation::ABMPopulation(int total_population)
    : impl(std::make_shared<Impl>(total_population))
{
}

namespace
{
// Apply the inference parameters in `params` to `model`
void apply_params(mio::abm::Model& model, const std::map<std::string, ScalarType>& params)
{
    for (const auto& [name, value] : params) {
        if (name == "beta") {
            model.parameters.get<mio::abm::InfectionRateFromViralShed>() = value;
        }
        else if (name == "t_exposed") {
            model.parameters.get<mio::abm::TimeExposedToNoSymptoms>() =
                mio::ParameterDistributionLogNormal(log(value), 1.);
        }
        else if (name == "time_presymptomatic") {
            model.parameters.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>() =
                mio::ParameterDistributionLogNormal(log(value), 1.);
        }
        else if (name == "time_asymptomatic_recovery") {
            model.parameters.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() =
                mio::ParameterDistributionLogNormal(log(value), 1.);
        }
        else if (name == "symptom_prob") {
            model.parameters.get<mio::abm::SymptomsPerInfectedNoSymptoms>() =
                mio::UncertainValue<ScalarType>(value);
        }
        else if (name == "quarantine_compliance") {
            for (auto& person : model.get_persons()) {
                person.set_compliance(mio::abm::InterventionType::Isolation, value);
            }
        }
        else {
            throw std::invalid_argument("forward_pass: unknown parameter '" + name + "'");
        }
    }
}
} // namespace

std::map<std::string, Eigen::MatrixXd> forward_pass(const ABMPopulation& population,
                                                    const std::map<std::string, ScalarType>& params,
                                                    const mio::abm::TestingBudget& design)
{
    std::call_once(g_log_level_once, [] {
        mio::set_log_level(mio::LogLevel::warn);
    });

    // Copy the pre-built population and set inference parameters on the copy.
    auto model = population.impl->model;
    apply_params(model, params);

    // Configure the surveillance design
    model.get_surveillance_testing().configure(design, mio::abm::PcrAssayParameters{}, kCareSeekPerDay);

    auto start_date = mio::abm::TimePoint(0);
    std::vector<ScalarType> infection_distribution{0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (auto& person : model.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, start_date, infection_state));
        }
    }

    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(30);
    auto sim  = mio::abm::Simulation(t0, std::move(model));

    mio::History<mio::DataWriterToMemory, mio::abm::LogSurveillanceRecords> history;
    sim.advance(tmax, history);

    // Flatten: positives -> a table (one row per positive test),
    //          negatives -> counts aggregated by (day, source, age, location).
    const auto& record_log = std::get<0>(history.get_log()); // vector<vector<PcrTestRecord>>, one entry per logged step

    std::vector<std::array<double, 6>> positives;             // [day, person_id, age, location, ct, source]
    std::map<std::array<int, 4>, int> negative_counts;        // key [day, source, age, location] -> count
    for (const auto& step_records : record_log) {
        for (const auto& r : step_records) {
            const int day = static_cast<int>(std::round(r.test_time.days()));
            const int age = static_cast<int>(r.age.get());
            const int loc = static_cast<int>(r.location_type);
            const int src = static_cast<int>(r.source);
            if (r.positive) {
                positives.push_back({static_cast<double>(day), static_cast<double>(r.person_id.get()),
                                     static_cast<double>(age), static_cast<double>(loc), r.ct_value,
                                     static_cast<double>(src)});
            }
            else {
                negative_counts[{day, src, age, loc}] += 1;
            }
        }
    }

    Eigen::MatrixXd pos_mat(static_cast<Eigen::Index>(positives.size()), 6);
    for (size_t i = 0; i < positives.size(); ++i) {
        for (int j = 0; j < 6; ++j) {
            pos_mat(static_cast<Eigen::Index>(i), j) = positives[i][static_cast<size_t>(j)];
        }
    }

    Eigen::MatrixXd neg_mat(static_cast<Eigen::Index>(negative_counts.size()), 5);
    {
        Eigen::Index i = 0;
        for (const auto& [key, count] : negative_counts) {
            neg_mat(i, 0) = key[0];
            neg_mat(i, 1) = key[1];
            neg_mat(i, 2) = key[2];
            neg_mat(i, 3) = key[3];
            neg_mat(i, 4) = count;
            ++i;
        }
    }

    // non-susceptible at tmax
    Eigen::Index n_ever_infected = 0;
    for (auto& person : sim.get_model().get_persons()) {
        if (person.get_infection_state(tmax) != mio::abm::InfectionState::Susceptible) {
            ++n_ever_infected;
        }
    }
    Eigen::MatrixXd meta(1, 1);
    meta(0, 0) = static_cast<double>(n_ever_infected);

    return {{"positives", pos_mat}, {"negatives", neg_mat}, {"n_ever_infected", meta}};
}
