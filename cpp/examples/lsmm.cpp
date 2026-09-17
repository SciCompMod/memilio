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
#include <iostream>

int main()
{
    using namespace mio::lsmm;
    Matrix initial(3, 2); // rows: S, I, R; columns: two groups at one location.
    initial << 90, 45, 10, 5, 0, 0;
    const double population = initial.sum();
    // Infection depends on I; recovery per person is constant. Source counts are handled automatically.
    const Model sir(3, {{0, 1, [population](const State& z) { return 0.3 * z[1] / population; },
                        std::vector<Eigen::Index>{1}},
                       {1, 2, [](const State&) { return 0.1; }, std::vector<Eigen::Index>{}}});
    mio::RandomNumberGenerator rng;
    rng.seed({1234});

    Simulation exact(sir, initial, rng);
    exact.advance(10.);
    std::cout << "Exact aggregate: " << exact.get_state().transpose() << '\n'
              << "History (current state, initial state):\n" << exact.get_history() << '\n'
              << "Group endpoints:\n" << exact.sample_endpoints() << '\n';

    const auto approximate = simulate_tau_leaping(sir, initial, rng, 0., 10., 0.05);
    std::cout << "Approximate binomial tau-leaping aggregate: " << approximate.state.transpose() << '\n';
}
