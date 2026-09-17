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
#include "lsmm/simulation.h"
#include "lsmm/tau_leaping.h"
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <limits>

namespace
{
using namespace mio::lsmm;
using Dependencies = std::vector<Eigen::Index>;

// Control integrated reaction thresholds without coupling these tests to another model.
struct ClockSequence {
    using Distribution = mio::ExponentialDistribution<double>;
    Distribution::GeneratorFunction previous = Distribution::get_instance().get_generator();
    explicit ClockSequence(std::vector<double> values)
    {
        Distribution::get_instance().set_generator([values = std::move(values), i = size_t(0)](const auto&) mutable {
            return values.at(i++);
        });
    }
    ~ClockSequence()
    {
        Distribution::get_instance().set_generator(previous);
    }
};

void check_margins(const Matrix& initial, const Matrix& history, const Matrix& endpoints)
{
    EXPECT_TRUE((endpoints.array() >= 0).all());
    EXPECT_TRUE(endpoints.colwise().sum() == initial.colwise().sum());
    EXPECT_TRUE(endpoints.rowwise().sum() == history.rowwise().sum());
    EXPECT_TRUE(history.colwise().sum().transpose() == initial.rowwise().sum());
}

TEST(LSMM, rejectsInvalidInputs)
{
    mio::RandomNumberGenerator rng;
    rng.seed({123});
    auto rate = [](const State&) { return 1.; };
    EXPECT_THROW(Model(0, {}), std::invalid_argument);
    EXPECT_THROW(Model(2, {{0, 2, rate}}), std::invalid_argument);
    EXPECT_THROW(Model(2, {{0, 0, rate}}), std::invalid_argument);
    EXPECT_THROW(Model(2, {{0, 1, {}}}), std::invalid_argument);
    EXPECT_THROW(Model(2, {{0, 1, rate, Dependencies{-1}}}), std::invalid_argument);
    EXPECT_THROW(Model(2, {{0, 1, rate, Dependencies{2}}}), std::invalid_argument);
    const Model model(2, {{0, 1, rate}});
    Matrix initial = Matrix::Zero(2, 1);
    EXPECT_THROW(Simulation(model, Matrix::Zero(3, 1), rng), std::invalid_argument);
    EXPECT_THROW(Simulation(model, Matrix::Zero(2, 0), rng), std::invalid_argument);
    initial(0, 0) = -1;
    EXPECT_THROW(Simulation(model, initial, rng), std::invalid_argument);
    initial(0, 0) = Count(std::numeric_limits<int>::max()) + 1;
    EXPECT_THROW(Simulation(model, initial, rng), std::invalid_argument);
    initial(0, 0) = 2;
    for (double bad : {-1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::max()}) {
        const Model invalid(2, {{0, 1, [bad](const State&) { return bad; }}});
        Simulation simulation(invalid, initial, rng);
        EXPECT_THROW(simulation.advance(1.), std::domain_error);
    }
    EXPECT_THROW(Simulation(model, initial, rng, std::numeric_limits<double>::infinity()), std::invalid_argument);
    Simulation simulation(model, initial, rng, 2.);
    EXPECT_THROW(simulation.advance(1.), std::invalid_argument);
    EXPECT_THROW(simulation.advance(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    EXPECT_THROW(reconstruct(initial, Matrix::Zero(3, 3), rng), std::invalid_argument);
    EXPECT_THROW(reconstruct(initial, Matrix::Zero(2, 2), rng), std::invalid_argument);
}

TEST(LSMM, emptySourcesAndZeroRates)
{
    mio::RandomNumberGenerator rng;
    rng.seed({124});
    const Model model(2, {{0, 1, [](const State&) -> double { throw std::logic_error("Empty source evaluated"); }},
                          {1, 0, [](const State&) { return 0.; }}});
    for (Count population : {0, 4}) {
        Matrix initial = Matrix::Zero(2, 1);
        initial(1, 0)  = population;
        Simulation simulation(model, initial, rng);
        simulation.advance(10.);
        EXPECT_EQ(simulation.get_time(), 10.);
        EXPECT_EQ(simulation.get_num_events(), 0u);
        EXPECT_TRUE(simulation.sample_endpoints() == initial);
    }
    Simulation no_channels(Model(2, {}), Matrix::Zero(2, 3), rng);
    no_channels.advance(1.);
    EXPECT_EQ(no_channels.get_num_events(), 0u);
}

TEST(LSMM, retainsClocksAcrossAdvanceCalls)
{
    ClockSequence clocks({2., 10.});
    mio::RandomNumberGenerator rng;
    rng.seed({125});
    Matrix initial(2, 1);
    initial << 1, 0;
    Simulation simulation(Model(2, {{0, 1, [](const State&) { return 1.; }}}), initial, rng, 5.);
    simulation.advance(5.);
    simulation.advance(6.);
    simulation.advance(6.5);
    EXPECT_EQ(simulation.get_num_events(), 0u);
    EXPECT_TRUE(simulation.get_state() == initial.col(0));
    simulation.advance(7.); // The pending event is included exactly at the requested endpoint.
    EXPECT_EQ(simulation.get_num_events(), 1u);
    EXPECT_EQ(simulation.get_state()[1], 1);
    simulation.advance(9.);
    EXPECT_EQ(simulation.get_num_events(), 1u);
    EXPECT_TRUE(simulation.sample_endpoints().col(0) == simulation.get_state());
}

TEST(LSMM, selectiveUpdatesPreserveTheSameRealization)
{
    // Dependencies can be unrelated to a channel's source or target, including nonlinear rates.
    const std::vector<Transition> selective{
        {0, 1, [](const State& z) { return 0.2 + 0.03 * z[2] * z[2]; }, Dependencies{2}},
        {0, 2, [](const State&) { return 0.1; }, Dependencies{}},
        {1, 2, [](const State& z) { return 0.2 + 0.01 * z[3]; }, Dependencies{3}},
        {2, 0, [](const State&) { return 0.4; }, Dependencies{}},
        {3, 1, [](const State& z) { return 0.5 + 0.02 * z[0]; }, Dependencies{0}},
        {2, 3, [](const State&) { return 0.2; }, Dependencies{}}};
    auto full = selective;
    for (auto& transition : full) {
        transition.dependencies.reset();
    }
    mio::RandomNumberGenerator full_rng, selective_rng;
    full_rng.seed({135});
    selective_rng.seed({135});
    Matrix initial(4, 3);
    initial << 3, 2, 0, 1, 2, 0, 0, 1, 3, 1, 0, 2;
    Simulation reference(Model(4, full), initial, full_rng);
    Simulation optimized(Model(4, selective), initial, selective_rng);
    for (int step = 1; step <= 80; ++step) {
        reference.advance(0.25 * step);
        optimized.advance(0.25 * step);
        EXPECT_EQ(optimized.get_num_events(), reference.get_num_events());
        EXPECT_TRUE(optimized.get_state() == reference.get_state());
        EXPECT_TRUE(optimized.get_history() == reference.get_history());
        EXPECT_TRUE(optimized.sample_endpoints() == reference.sample_endpoints());
    }
    EXPECT_GT(optimized.get_num_events(), 0u);
}

TEST(LSMM, unaffectedAndConstantRatesRemainCached)
{
    ClockSequence clocks({1., 100., 100., 100.});
    mio::RandomNumberGenerator rng;
    rng.seed({136});
    std::array<int, 3> calls{};
    const Model model(4, {{0, 1, [&calls](const State&) {
                              ++calls[0];
                              return 1.;
                          }, Dependencies{}},
                          {2, 3, [&calls](const State& z) {
                              ++calls[1];
                              return 1. + z[3];
                          }, Dependencies{3}},
                          {1, 0, [&calls](const State&) {
                              ++calls[2];
                              return 1.;
                          }, Dependencies{}}});
    Matrix initial(4, 1);
    initial << 2, 0, 1, 1;
    Simulation simulation(model, initial, rng);
    simulation.advance(0.25);
    EXPECT_EQ(calls, (std::array<int, 3>{1, 1, 0}));
    simulation.advance(0.6); // Only 0->1 fires, leaving the second channel untouched.
    ASSERT_EQ(simulation.get_num_events(), 1u);
    EXPECT_EQ(calls, (std::array<int, 3>{1, 1, 1}));
    simulation.advance(0.7);
    simulation.advance(0.8);
    EXPECT_EQ(calls, (std::array<int, 3>{1, 1, 1}));
}

TEST(LSMM, dependencyChangesWhileSourceIsEmpty)
{
    // 0->1 empties its source at t=0.1. Its rate changes at t=0.2 while empty,
    // and 3->0 reactivates it at t=0.3. The next 0->1 event must occur at t=0.5.
    ClockSequence clocks({0.3, 0.2, 0.1, 0.4, 100., 100., 100.});
    mio::RandomNumberGenerator rng;
    rng.seed({137});
    int calls = 0;
    const Model model(4, {{0, 1, [&calls](const State& z) {
                              if (z[0] == 0) {
                                  throw std::logic_error("Empty source evaluated");
                              }
                              ++calls;
                              return 2. + z[2];
                          }, Dependencies{2}},
                          {2, 3, [](const State&) { return 1.; }, Dependencies{}},
                          {3, 0, [](const State&) { return 1.; }, Dependencies{}}});
    Matrix initial(4, 1);
    initial << 1, 0, 1, 0;
    Simulation simulation(model, initial, rng);
    EXPECT_NO_THROW(simulation.advance(0.49));
    ASSERT_EQ(simulation.get_num_events(), 3u);
    EXPECT_EQ(simulation.get_state()[0], 1);
    EXPECT_EQ(calls, 2);
    EXPECT_NO_THROW(simulation.advance(0.51));
    EXPECT_EQ(simulation.get_num_events(), 4u);
    EXPECT_EQ(simulation.get_state()[1], 2);
    EXPECT_EQ(calls, 2);
}

TEST(LSMM, cachedRatesStillCheckTotalRateOverflow)
{
    const double largest = std::numeric_limits<double>::max();
    ClockSequence clocks({largest, 0.25, 100.});
    mio::RandomNumberGenerator rng;
    rng.seed({138});
    int calls = 0;
    const Model model(3, {{0, 1, [&calls, largest](const State&) {
                              ++calls;
                              return 0.75 * largest;
                          }, Dependencies{}},
                          {2, 0, [](const State&) { return 1.; }, Dependencies{}}});
    Matrix initial(3, 1);
    initial << 1, 0, 1;
    Simulation simulation(model, initial, rng);
    // The rate is initially finite, but 2->0 doubles its source count at t=0.25.
    EXPECT_THROW(simulation.advance(0.5), std::domain_error);
    EXPECT_EQ(calls, 1);
}

TEST(LSMM, historyPreservesRecoveryBeforeInfection)
{
    // Group A starts with one S, B with two I. Recovery at 0.25 precedes infection at 0.75.
    // Distributing summed transitions in S->I->R order would incorrectly let A recover.
    ClockSequence clocks({0.5, 1., 100., 100.});
    mio::RandomNumberGenerator rng;
    rng.seed({126});
    const Model model(3, {{1, 2, [](const State&) { return 1.; }},
                          {0, 1, [](const State& z) { return double(z[1]); }}});
    Matrix initial(3, 2), expected(3, 2), history(3, 3);
    initial << 1, 0, 0, 2, 0, 0;
    expected << 0, 0, 1, 1, 0, 1;
    history << 0, 0, 0, 1, 1, 0, 0, 1, 0;
    Simulation simulation(model, initial, rng);
    simulation.advance(1.);
    ASSERT_EQ(simulation.get_num_events(), 2u);
    EXPECT_TRUE(simulation.get_history() == history);
    EXPECT_TRUE(simulation.sample_endpoints() == expected);
}

TEST(LSMM, branchesAndCyclesConserveBothMargins)
{
    mio::RandomNumberGenerator rng;
    rng.seed({127});
    const Model model(3, {{0, 1, [](const State&) { return 0.8; }},
                          {0, 2, [](const State&) { return 0.3; }},
                          {1, 2, [](const State&) { return 0.5; }},
                          {2, 0, [](const State&) { return 0.7; }}});
    Matrix initial(3, 3);
    initial << 3, 1, 0, 1, 0, 3, 0, 2, 1;
    Simulation simulation(model, initial, rng);
    for (int t = 1; t <= 20; ++t) {
        simulation.advance(t);
        EXPECT_TRUE((simulation.get_history().array() >= 0).all());
        EXPECT_TRUE(simulation.get_history().rowwise().sum() == simulation.get_state());
        check_margins(initial, simulation.get_history(), simulation.sample_endpoints());
    }
    EXPECT_GT(simulation.get_num_events(), 0u);
}

TEST(LSMM, aggregateDeathLaw)
{
    mio::RandomNumberGenerator rng;
    rng.seed({128});
    const Model model(2, {{0, 1, [](const State&) { return 1.; }}});
    Matrix initial(2, 2);
    initial << 1, 2, 0, 0;
    constexpr int samples = 16000;
    std::array<int, 4> histogram{};
    for (int n = 0; n < samples; ++n) {
        Simulation simulation(model, initial, rng);
        ++histogram[simulation.advance(std::log(2.))[0]];
    }
    const std::array<double, 4> probabilities{0.125, 0.375, 0.375, 0.125};
    for (size_t k = 0; k < histogram.size(); ++k) {
        const auto p = probabilities[k]; // Three independent survivors, each with probability 1/2.
        EXPECT_NEAR(histogram[k], samples * p, 7 * std::sqrt(samples * p * (1 - p)) + 2);
    }
}

TEST(LSMM, conditionalLawWithMixedInitialCompartments)
{
    mio::RandomNumberGenerator rng;
    rng.seed({129});
    Matrix initial(2, 2), history(2, 2);
    initial << 2, 1, 1, 2;
    history << 1, 2, 2, 1;
    constexpr int samples = 18000;
    std::array<int, 3> histogram{};
    for (int n = 0; n < samples; ++n) {
        const Matrix endpoints = reconstruct(initial, history, rng);
        check_margins(initial, history, endpoints);
        ASSERT_GE(endpoints(0, 0), 0);
        ASSERT_LE(endpoints(0, 0), 2);
        ++histogram[endpoints(0, 0)];
    }
    // Independent initial-state allocations give P(X_A,0 = 0,1,2) = (1,4,4)/9.
    const std::array<double, 3> probabilities{1. / 9, 4. / 9, 4. / 9};
    for (size_t k = 0; k < histogram.size(); ++k) {
        const auto p = probabilities[k];
        EXPECT_NEAR(histogram[k], samples * p, 7 * std::sqrt(samples * p * (1 - p)) + 2);
    }
}

TEST(LSMM, nonlinearSIRLaw)
{
    mio::RandomNumberGenerator rng;
    rng.seed({133});
    const Model model(3, {{0, 1, [](const State& z) { return 2. * z[1] / z.sum(); }},
                          {1, 2, [](const State&) { return 1.; }}});
    Matrix initial(3, 1);
    initial << 1, 1, 0;
    constexpr int samples = 16000;
    std::array<int, 5> histogram{};
    for (int n = 0; n < samples; ++n) {
        Simulation simulation(model, initial, rng);
        const State& z = simulation.advance(1.);
        const size_t outcome = z[0] == 1 ? (z[1] == 1 ? 0 : 1) : 2 + (2 - z[1]);
        ASSERT_LT(outcome, histogram.size());
        ++histogram[outcome];
    }
    // (S,I) = (1,1), (1,0), (0,2), (0,1), (0,0): analytic forward-equation solution.
    const double e = std::exp(-1.);
    const std::array<double, 5> probabilities{
        e * e, (1. - e * e) / 2., e * e, 2 * e - 4 * e * e, 0.5 - 2 * e + 2.5 * e * e};
    for (size_t k = 0; k < histogram.size(); ++k) {
        const auto p = probabilities[k];
        EXPECT_NEAR(histogram[k], samples * p, 7 * std::sqrt(samples * p * (1 - p)) + 2);
    }
}

TEST(LSMM, onlineHistorySamplesMixedOrigins)
{
    mio::RandomNumberGenerator rng;
    rng.seed({134});
    const Model model(2, {{0, 1, [](const State&) { return 1.; }},
                          {1, 0, [](const State&) { return 1.; }}});
    const Matrix initial = Matrix::Identity(2, 2);
    constexpr int samples = 16000;
    std::array<int, 4> histogram{};
    for (int n = 0; n < samples; ++n) {
        Simulation simulation(model, initial, rng);
        simulation.advance(std::log(2.) / 2.);
        const Matrix& h = simulation.get_history();
        ++histogram[2 * h(0, 0) + h(1, 1)];
    }
    // Each original individual independently stays in its initial state with probability 3/4.
    const std::array<double, 4> probabilities{1. / 16, 3. / 16, 3. / 16, 9. / 16};
    for (size_t k = 0; k < histogram.size(); ++k) {
        const auto p = probabilities[k];
        EXPECT_NEAR(histogram[k], samples * p, 7 * std::sqrt(samples * p * (1 - p)) + 2);
    }
}

TEST(LSMM, hypergeometricCombinatorialLaw)
{
    mio::RandomNumberGenerator rng;
    rng.seed({132});
    constexpr int samples = 24000;
    std::array<int, 10> histogram{};
    for (int n = 0; n < samples; ++n) {
        ++histogram[details::hypergeometric(30, 12, 9, rng)];
    }
    const auto choose = [](int n, int k) {
        double value = 1.;
        for (int i = 1; i <= k; ++i) {
            value *= double(n - i + 1) / i;
        }
        return value;
    };
    for (int k = 0; k <= 9; ++k) {
        const double p = choose(12, k) * choose(18, 9 - k) / choose(30, 9);
        EXPECT_NEAR(histogram[k], samples * p, 7 * std::sqrt(samples * p * (1 - p)) + 2);
    }
}

TEST(LSMM, hypergeometricBoundariesAndLargeCounts)
{
    mio::RandomNumberGenerator rng;
    rng.seed({130});
    EXPECT_EQ(details::hypergeometric(0, 0, 0, rng), 0);
    EXPECT_EQ(details::hypergeometric(8, 0, 3, rng), 0);
    EXPECT_EQ(details::hypergeometric(8, 3, 8, rng), 3);
    EXPECT_EQ(details::hypergeometric(8, 8, 3, rng), 3);
    const Count balanced = details::hypergeometric(1000000, 500000, 500000, rng);
    EXPECT_GE(balanced, 245000);
    EXPECT_LE(balanced, 255000);
    const Count total = std::numeric_limits<int>::max();
    for (int n = 0; n < 20; ++n) {
        const Count x = details::hypergeometric(total, total - 2, total - 3, rng);
        EXPECT_GE(x, total - 5);
        EXPECT_LE(x, total - 3);
    }
    Matrix initial(2, 2), history(2, 2);
    initial << total - 3, 3, 0, 0;
    history << total - 2, 0, 2, 0;
    check_margins(initial, history, reconstruct(initial, history, rng));
}

TEST(LSMM, tauLeapingUsesOldHistoryAndPreservesMargins)
{
    mio::RandomNumberGenerator rng;
    rng.seed({131});
    const Model cycle(2, {{0, 1, [](const State&) { return 1000.; }},
                          {1, 0, [](const State&) { return 1000.; }}});
    Matrix initial(2, 2), expected(2, 2);
    initial << 3, 1, 1, 2;
    expected << 1, 2, 3, 1;
    const auto result = simulate_tau_leaping(cycle, initial, rng, 0., 1., 1.);
    EXPECT_TRUE(result.endpoints == expected); // Everyone moves once; arrivals cannot move again in the step.
    EXPECT_EQ(result.history.diagonal().sum(), 0);
    check_margins(initial, result.history, result.endpoints);
    EXPECT_THROW(simulate_tau_leaping(cycle, initial, rng, 0., 1., 0.), std::invalid_argument);
    EXPECT_TRUE(simulate_tau_leaping(cycle, initial, rng, 2., 2., 1.).endpoints == initial);
    EXPECT_EQ(simulate_tau_leaping(cycle, Matrix::Zero(2, 2), rng, 0., 2., 1.).state.sum(), 0);

    const Model branch(3, {{0, 1, [](const State&) { return 1000.; }},
                           {0, 2, [](const State&) { return 2000.; }},
                           {1, 0, [](const State&) { return 1000.; }},
                           {2, 0, [](const State&) { return 1000.; }}});
    initial = Matrix::Constant(3, 3, 3);
    const auto branching = simulate_tau_leaping(branch, initial, rng, 0., 1., 1.);
    EXPECT_EQ(branching.history.diagonal().sum(), 0);
    check_margins(initial, branching.history, branching.endpoints);
}

TEST(LSMM, tauLeapingOneTransitionAndChainConvergence)
{
    mio::RandomNumberGenerator rng;
    rng.seed({132});
    constexpr Count population = 1000000;
    Matrix initial = Matrix::Zero(3, 1);
    initial(0, 0) = population;
    const Model death(3, {{0, 1, [](const State&) { return 1.; }}});
    const auto once = simulate_tau_leaping(death, initial, rng, 2., 2. + std::log(2.), 1.);
    EXPECT_NEAR(once.state[0], population / 2., 7 * std::sqrt(population / 4.));

    const Model chain(3, {{0, 1, [](const State&) { return 1.; }},
                          {1, 2, [](const State&) { return 1.; }}});
    const double exact = 1. - 2. * std::exp(-1.);
    double coarse_error = 0.;
    for (int steps : {2, 100}) {
        const auto result = simulate_tau_leaping(chain, initial, rng, 0., 1., 1. / steps);
        const double stay = std::exp(-1. / steps);
        const double expected = 1. - std::pow(stay, steps) - steps * (1. - stay) * std::pow(stay, steps - 1);
        EXPECT_NEAR(result.state[2], population * expected,
                    7 * std::sqrt(population * expected * (1. - expected)) + 2);
        const double error = std::abs(double(result.state[2]) / population - exact);
        if (steps == 2) {
            coarse_error = error;
        }
        else {
            EXPECT_LT(error, 0.01);
            EXPECT_LT(error, coarse_error / 5.);
        }
    }
}
} // namespace
