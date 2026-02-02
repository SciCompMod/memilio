/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Julia Bicker, René Schmieding, Kilian Volmer
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
#include "memilio/epidemiology/age_group.h"
#include "smm/simulation.h"
#include "smm/model.h"
#include "smm/parameters.h"
#include "memilio/data/analyze_result.h"
#include <cstdlib>
#include "memilio/timer/auto_timer.h"

#pragma GCC target("avx2")

enum class InfectionState
{
    S,
    E,
    I,
    Sus,
    Lab,
    C,
    Desinfected,
    Count
};

struct FarmType : public mio::Index<FarmType> {
    FarmType(size_t val)
        : Index<FarmType>(val)
    {
    }
};

int main()
{
    using Status = mio::Index<InfectionState, FarmType>;
    using mio::regions::Region;
    using enum InfectionState;

    //Example how to run the stochastic metapopulation models with four regions
    const size_t num_regions    = 100; //557;
    const size_t num_farm_types = 6;
    using Model                 = mio::smm::Model<ScalarType, InfectionState, Status, Region>;

    Model model(Status{Count, FarmType(num_farm_types)}, Region(num_regions));

    for (size_t r = 0; r < num_regions; ++r) {
        mio::timing::AutoTimer<"Population"> timer;

        for (size_t f = 0; f < num_farm_types; ++f) {
            model.populations[{Region(r), S, FarmType(f)}]           = 10.0;
            model.populations[{Region(r), E, FarmType(f)}]           = 0.0;
            model.populations[{Region(r), I, FarmType(f)}]           = 0.0;
            model.populations[{Region(r), Sus, FarmType(f)}]         = 0.0;
            model.populations[{Region(r), Lab, FarmType(f)}]         = 0.0;
            model.populations[{Region(r), C, FarmType(f)}]           = 0.0;
            model.populations[{Region(r), Desinfected, FarmType(f)}] = 0.0;
        }
    }

    using AR = mio::smm::AdoptionRates<ScalarType, Status, Region>;
    using TR = mio::smm::TransitionRates<ScalarType, Status, Region>;

    //Set infection state adoption and spatial transition rates
    AR::Type adoption_rates;
    TR::Type transition_rates;
    for (size_t r = 0; r < num_regions; ++r) {
        mio::timing::AutoTimer<"AdoptionRates"> timer;

        for (size_t f = 0; f < num_farm_types - 1; ++f) {
            adoption_rates.push_back({{S, FarmType(f)}, {E, FarmType(f)}, Region(r), 0.1, {{{I, FarmType(5)}, 1}}});
            adoption_rates.push_back({{E, FarmType(f)}, {I, FarmType(f)}, Region(r), 1. / 3., {}});
            adoption_rates.push_back({{I, FarmType(f)}, {Sus, FarmType(f)}, Region(r), 1. / 3., {}});
            adoption_rates.push_back({{Sus, FarmType(f)}, {Lab, FarmType(f)}, Region(r), 1. / 2.6, {}});
            adoption_rates.push_back({{Lab, FarmType(f)}, {C, FarmType(f)}, Region(r), 1. / 2.6, {}});
            adoption_rates.push_back({{C, FarmType(f)}, {Desinfected, FarmType(f)}, Region(r), 1. / 2., {}});
            adoption_rates.push_back({{Desinfected, FarmType(f)}, {S, FarmType(f)}, Region(r), 1. / 22., {}});
            adoption_rates.push_back({{S, FarmType(f)}, {Desinfected, FarmType(f)}, Region(r), 1. / 21., {}});
        }
    }

    for (size_t r = 0; r < num_regions; ++r) {
        adoption_rates.push_back({{S, FarmType(5)}, {E, FarmType(5)}, Region(r), 0.1, {{{I, FarmType(5)}, 1}}});
        adoption_rates.push_back({{E, FarmType(5)}, {I, FarmType(5)}, Region(r), 1. / 3., {}});
    }

    for (size_t i = 0; i < num_regions; ++i) {
        for (size_t j = 0; j < num_regions; ++j) {
            if (i != j && (std::abs(int(i) - int(j)) <= 6)) {
                transition_rates.push_back({{S, FarmType(5)}, Region(i), Region(j), 0.01});
                transition_rates.push_back({{E, FarmType(5)}, Region(j), Region(i), 0.01});
                transition_rates.push_back({{I, FarmType(5)}, Region(j), Region(i), 0.01});
            }
        }
    }

    model.parameters.get<AR>() = adoption_rates;
    model.parameters.get<TR>() = transition_rates;

    ScalarType dt   = 10.0;
    ScalarType tmax = 90.0;

    auto sim = mio::smm::Simulation(model, 0.0, dt);
    mio::log_info("Starting simulation...");
    sim.advance(tmax);
    mio::log_info("Simulation finished.");
    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());
    // interpolated_results.print_table({"S", "E", "C", "I", "R", "D "});
}
