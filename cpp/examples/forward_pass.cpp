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
#include "city_builder.h"
#include "forward_pass.h"

namespace
{
const auto AGE_GROUP_5_TO_14 = mio::AgeGroup(1);
} // namespace

struct ABMPopulation::Impl {
    mio::abm::Model model;

    explicit Impl(int total_population)
        : model(CityBuilder::build_world(CityConfig{total_population}, mio::thread_local_rng()))
    {
        // People can get tested at work (and do this with 0.5 probability) from time point 0 to day 10.
        auto testing_criteria_work = mio::abm::TestingCriteria();
        auto testing_scheme_work   = mio::abm::TestingScheme(
            testing_criteria_work, mio::abm::days(1), mio::abm::TimePoint(0),
            mio::abm::TimePoint(0) + mio::abm::days(10),
            model.parameters.get<mio::abm::TestData>()[mio::abm::TestType::Antigen], 0.5);
        model.get_testing_strategy().add_scheme(mio::abm::LocationType::Work, testing_scheme_work);
    }
};

ABMPopulation::ABMPopulation(int total_population)
    : impl(std::make_shared<Impl>(total_population))
{
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> forward_pass(const ABMPopulation& population,
                                                          ScalarType beta, ScalarType kappa,
                                                          int cohort_budget)
{
    mio::set_log_level(mio::LogLevel::warn);

    // Copy the pre-built population, then set inference parameters on the copy.
    auto model = population.impl->model;
    model.parameters.get<mio::abm::TimeExposedToNoSymptoms>() =
        mio::ParameterDistributionLogNormal(log(kappa), 1.);
    model.parameters.check_constraints();
    model.parameters.get<mio::abm::AerosolTransmissionRates>() = 1 / beta;

    // Randomise initial infection states on the copied model.
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

    mio::History<mio::DataWriterToMemory, mio::abm::LogCTCluster, mio::abm::LogSchoolCohort>
        historyTimeSeries{
            mio::abm::LogCTCluster{AGE_GROUP_5_TO_14, mio::abm::LocationType::School},
            mio::abm::LogSchoolCohort{cohort_budget}};
    sim.advance(tmax, historyTimeSeries);

    const auto& hist_log = std::get<0>(historyTimeSeries.get_log());
    Eigen::MatrixXd hist_result(static_cast<Eigen::Index>(hist_log.size()), 42);
    for (size_t i = 0; i < hist_log.size(); ++i) {
        const auto& [time, hist] = hist_log[i];
        hist_result(i, 0)           = std::round(time.days());
        hist_result.row(i).tail(41) = hist.cast<double>().transpose();
    }

    const auto& cohort_log = std::get<1>(historyTimeSeries.get_log());
    Eigen::MatrixXd cohort_result(static_cast<Eigen::Index>(cohort_log.size()), cohort_budget + 1);
    for (size_t i = 0; i < cohort_log.size(); ++i) {
        const auto& [time, cts] = cohort_log[i];
        cohort_result(i, 0) = std::round(time.days());
        for (int j = 0; j < cohort_budget; ++j)
            cohort_result(i, j + 1) = static_cast<double>(cts[j]);
    }

    return {hist_result, cohort_result};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> forward_pass(ScalarType beta, ScalarType kappa, int cohort_budget)
{
    return forward_pass(ABMPopulation{}, beta, kappa, cohort_budget);
}
