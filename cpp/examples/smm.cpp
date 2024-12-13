/* 
* Copyright (C) 2020-2024 German Aerospace Center (DLR-SC)
*
* Authors: Julia Bicker, René Schmieding
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

#include "smm/simulation.h"
#include "smm/parameters.h"
#include "memilio/data/analyze_result.h"

enum class InfectionState
{
    S,
    E,
    C,
    I,
    R,
    D,
    Count

};

int main()
{

    //Example how to run the stochastic metapopulation models with four regions
    const size_t num_regions = 4;
    using Model              = mio::smm::Model<num_regions, InfectionState>;

    double numE = 12, numC = 4, numI = 12, numR = 0, numD = 0;

    Model model;
    //Population are distributed uniformly to the four regions
    for (size_t r = 0; r < num_regions; ++r) {
        model.populations[{mio::smm::Region(r), InfectionState::S}] =
            (1000 - numE - numC - numI - numR - numD) / num_regions;
        model.populations[{mio::smm::Region(r), InfectionState::E}] = numE / num_regions;
        model.populations[{mio::smm::Region(r), InfectionState::C}] = numC / num_regions;
        model.populations[{mio::smm::Region(r), InfectionState::I}] = numI / num_regions;
        model.populations[{mio::smm::Region(r), InfectionState::R}] = numR / num_regions;
        model.populations[{mio::smm::Region(r), InfectionState::D}] = numD / num_regions;
    }

    //Set infection state adoption and spatial transition rates
    std::vector<mio::smm::AdoptionRate<InfectionState>> adoption_rates;
    std::vector<mio::smm::TransitionRate<InfectionState>> transition_rates;
    for (size_t r = 0; r < num_regions; ++r) {
        adoption_rates.push_back({InfectionState::S,
                                  InfectionState::E,
                                  mio::smm::Region(r),
                                  0.1,
                                  {InfectionState::C, InfectionState::I},
                                  {1, 0.5}});
        adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::smm::Region(r), 1.0 / 5., {}, {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::R, mio::smm::Region(r), 0.2 / 3., {}, {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::I, mio::smm::Region(r), 0.8 / 3., {}, {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::R, mio::smm::Region(r), 0.99 / 5., {}, {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::D, mio::smm::Region(r), 0.01 / 5., {}, {}});
    }

    //Agents in infection state D do not transition
    for (size_t s = 0; s < static_cast<size_t>(InfectionState::D); ++s) {
        for (size_t i = 0; i < num_regions; ++i) {
            for (size_t j = 0; j < num_regions; ++j)
                if (i != j) {
                    transition_rates.push_back({InfectionState(s), mio::smm::Region(i), mio::smm::Region(j), 0.01});
                    transition_rates.push_back({InfectionState(s), mio::smm::Region(j), mio::smm::Region(i), 0.01});
                }
        }
    }

    model.parameters.get<mio::smm::AdoptionRates<InfectionState>>()   = adoption_rates;
    model.parameters.get<mio::smm::TransitionRates<InfectionState>>() = transition_rates;

    double dt   = 0.1;
    double tmax = 30.;

    auto sim = mio::Simulation<double, Model>(model, 0.0, dt);
    sim.advance(tmax);

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());
    interpolated_results.print_table({"S", "E", "C", "I", "R", "D "});
}
