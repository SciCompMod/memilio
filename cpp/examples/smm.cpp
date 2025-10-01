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

#include "memilio/epidemiology/age_group.h"
#include "smm/simulation.h"
#include "smm/model.h"
#include "smm/parameters.h"
#include "memilio/data/analyze_result.h"
#include "memilio/epidemiology/adoption_rate.h"

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

using Age     = mio::AgeGroup;
using Species = mio::AgeGroup;

int main()
{

    //Example how to run the stochastic metapopulation models with four regions
    const size_t num_regions    = 4;
    const size_t num_age_groups = 1;
    const size_t num_groups     = 1;
    using Model                 = mio::smm::Model<num_regions, InfectionState, num_age_groups, num_groups>;

    double numE = 12, numC = 4, numI = 12, numR = 0, numD = 0;

    Model model;
    //Population are distributed uniformly to the four regions
    for (size_t r = 0; r < num_regions; ++r) {
        model.populations[{mio::regions::Region(r), InfectionState::S, Age(0), Species(0)}] =
            (1000 - numE - numC - numI - numR - numD) / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::E, Age(0), Species(0)}] = numE / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::C, Age(0), Species(0)}] = numC / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::I, Age(0), Species(0)}] = numI / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::R, Age(0), Species(0)}] = numR / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::D, Age(0), Species(0)}] = numD / num_regions;
    }

    //Set infection state adoption and spatial transition rates
    std::vector<mio::AdoptionRate<InfectionState, Age, Species>> adoption_rates;
    std::vector<mio::smm::TransitionRate<InfectionState, Age, Species>> transition_rates;
    for (size_t r = 0; r < num_regions; ++r) {
        adoption_rates.push_back({InfectionState::S,
                                  InfectionState::E,
                                  mio::regions::Region(r),
                                  0.1,
                                  {{InfectionState::C, 1, mio::regions::Region(3), {Age(0), Species(0)}},
                                   {InfectionState::I, 0.5, mio::regions::Region(1), {Age(0), Species(0)}}},
                                  {Age(0), Species(0)}});
        adoption_rates.push_back(
            {InfectionState::C, InfectionState::R, mio::regions::Region(r), 0.2 / 3., {}, {Age(0), Species(0)}});
        adoption_rates.push_back(
            {InfectionState::E, InfectionState::C, mio::regions::Region(r), 1.0 / 5., {}, {Age(0), Species(0)}});
        adoption_rates.push_back(
            {InfectionState::C, InfectionState::I, mio::regions::Region(r), 0.8 / 3., {}, {Age(0), Species(0)}});
        adoption_rates.push_back(
            {InfectionState::I, InfectionState::R, mio::regions::Region(r), 0.99 / 5., {}, {Age(0), Species(0)}});
        adoption_rates.push_back(
            {InfectionState::I, InfectionState::D, mio::regions::Region(r), 0.01 / 5., {}, {Age(0), Species(0)}});
    }

    //Agents in infection state D do not transition
    for (size_t s = 0; s < static_cast<size_t>(InfectionState::D); ++s) {
        for (size_t i = 0; i < num_regions; ++i) {
            for (size_t j = 0; j < num_regions; ++j)
                if (i != j) {
                    transition_rates.push_back({InfectionState(s),
                                                mio::regions::Region(i),
                                                mio::regions::Region(j),
                                                0.01,
                                                {Age(0), Species(0)},
                                                {Age(0), Species(0)}});
                    transition_rates.push_back({InfectionState(s),
                                                mio::regions::Region(j),
                                                mio::regions::Region(i),
                                                0.01,
                                                {Age(0), Species(0)},
                                                {Age(0), Species(0)}});
                }
        }
    }

    model.parameters.get<mio::smm::AdoptionRates<InfectionState, Age, Species>>()   = adoption_rates;
    model.parameters.get<mio::smm::TransitionRates<InfectionState, Age, Species>>() = transition_rates;

    double dt   = 0.1;
    double tmax = 30.;

    auto sim = mio::smm::Simulation(model, 0.0, dt);
    sim.advance(tmax);

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());
    interpolated_results.print_table({"S", "E", "C", "I", "R", "D "});
}
