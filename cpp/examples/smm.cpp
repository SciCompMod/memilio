/*
* Copyright (C) 2020-2025 MEmilio
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

struct Species : public mio::Index<Species> {
    Species(size_t val)
        : Index<Species>(val)
    {
    }
};

int main()
{
    using Age    = mio::AgeGroup;
    using Status = mio::Index<InfectionState, Age, Species>;
    using mio::regions::Region;
    using enum InfectionState;

    /* Example how to run the stochastic metapopulation models with four regions. Within each region we differentiate by
       age groups, species and infection states. The infection states are S, E, C, I, R, D. For the number of age groups
       and species we choose: */
    const size_t num_regions    = 4;
    const size_t num_age_groups = 1;
    const size_t num_species    = 1;
    using Model                 = mio::smm::Model<ScalarType, InfectionState, Status, Region>;

    ScalarType numE = 12, numC = 4, numI = 12, numR = 0, numD = 0;

    Model model(Status{Count, Age(num_age_groups), Species(num_species)}, Region(num_regions));
    //Population are distributed uniformly to the four regions
    for (size_t r = 0; r < num_regions; ++r) {
        model.populations[{Region(r), S, Age(0), Species(0)}] = (1000 - numE - numC - numI - numR - numD) / num_regions;
        model.populations[{Region(r), E, Age(0), Species(0)}] = numE / num_regions;
        model.populations[{Region(r), C, Age(0), Species(0)}] = numC / num_regions;
        model.populations[{Region(r), I, Age(0), Species(0)}] = numI / num_regions;
        model.populations[{Region(r), R, Age(0), Species(0)}] = numR / num_regions;
        model.populations[{Region(r), D, Age(0), Species(0)}] = numD / num_regions;
    }

    using AR = mio::smm::AdoptionRates<ScalarType, Status, Region>;
    using TR = mio::smm::TransitionRates<ScalarType, Status, Region>;

    //Set infection state adoption rates. Adoptions only happen within a region.
    AR::Type adoption_rates;
    for (size_t r = 0; r < num_regions; ++r) {
        adoption_rates.push_back({{S, Age(0), Species(0)},
                                  {E, Age(0), Species(0)},
                                  Region(r),
                                  0.1,
                                  {{{C, Age(0), Species(0)}, 1}, {{I, Age(0), Species(0)}, 0.5}}});
        adoption_rates.push_back({{C, Age(0), Species(0)}, {R, Age(0), Species(0)}, Region(r), 0.2 / 3., {}});
        adoption_rates.push_back({{E, Age(0), Species(0)}, {C, Age(0), Species(0)}, Region(r), 1.0 / 5., {}});
        adoption_rates.push_back({{C, Age(0), Species(0)}, {I, Age(0), Species(0)}, Region(r), 0.8 / 3., {}});
        adoption_rates.push_back({{I, Age(0), Species(0)}, {R, Age(0), Species(0)}, Region(r), 0.99 / 5., {}});
        adoption_rates.push_back({{I, Age(0), Species(0)}, {D, Age(0), Species(0)}, Region(r), 0.01 / 5., {}});
    }

    //Set transition rates between regions. Agents in infection state D do not transition.
    TR::Type transition_rates;
    for (size_t s = 0; s < static_cast<size_t>(D); ++s) {
        for (size_t i = 0; i < num_regions; ++i) {
            for (size_t j = 0; j < num_regions; ++j)
                if (i != j) {
                    transition_rates.push_back({{InfectionState(s), Age(0), Species(0)}, Region(i), Region(j), 0.01});
                    transition_rates.push_back({{InfectionState(s), Age(0), Species(0)}, Region(j), Region(i), 0.01});
                }
        }
    }

    model.parameters.get<AR>() = adoption_rates;
    model.parameters.get<TR>() = transition_rates;

    ScalarType dt   = 0.1;
    ScalarType tmax = 30.0;

    auto sim = mio::smm::Simulation(model, 0.0, dt);
    sim.advance(tmax);

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());
    interpolated_results.print_table({"S", "E", "C", "I", "R", "D "});
}
