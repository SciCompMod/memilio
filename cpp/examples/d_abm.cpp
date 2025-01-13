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

#include "d_abm/quad_well.h"
#include "d_abm/simulation.h"
#include "d_abm/parameters.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/data/analyze_result.h"
#include <vector>

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
    //Example how to run a simulation of the diffusive ABM using the quadwell potential
    using Model = mio::dabm::Model<QuadWellModel<InfectionState>>;
    std::vector<Model::Agent> agents(1000);
    //Random variables for initialization of agents' position and infection state
    auto& pos_rng = mio::UniformDistribution<double>::get_instance();
    auto& sta_rng = mio::DiscreteDistribution<size_t>::get_instance();
    //Infection state distribution
    std::vector<double> pop_dist{0.98, 0.01, 0.005, 0.005, 0., 0.};
    for (auto& a : agents) {
        //Agents are equally distributed in [-2,2]x[-2,2] at the beginning
        a.position =
            Eigen::Vector2d{pos_rng(mio::thread_local_rng(), -2., 2.), pos_rng(mio::thread_local_rng(), -2., 2.)};
        a.status = static_cast<InfectionState>(sta_rng(mio::thread_local_rng(), pop_dist));
    }

    //Set adoption rates
    std::vector<mio::dabm::AdoptionRate<InfectionState>> adoption_rates;
    for (size_t region = 0; region < 4; ++region) {
        adoption_rates.push_back({InfectionState::S,
                                  InfectionState::E,
                                  mio::dabm::Region(region),
                                  0.1,
                                  {InfectionState::C, InfectionState::I},
                                  {1, 0.5}});
        adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::dabm::Region(region), 1.0 / 5., {}, {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::R, mio::dabm::Region(region), 0.2 / 3., {}, {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::I, mio::dabm::Region(region), 0.8 / 3., {}, {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::R, mio::dabm::Region(region), 0.99 / 5., {}, {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::D, mio::dabm::Region(region), 0.01 / 5., {}, {}});
    }

    //Set interaction radius and noise term of the diffusion process
    double interaction_radius = 0.5;
    double noise              = 0.4;

    double dt   = 0.1;
    double tmax = 30.;

    Model model(agents, adoption_rates, interaction_radius, noise, {InfectionState::D});
    auto sim = mio::Simulation<double, Model>(model, 0.0, dt);
    sim.advance(tmax);

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());
    interpolated_results.print_table({"S", "E", "C", "I", "R", "D "});
}
