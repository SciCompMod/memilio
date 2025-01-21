/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Jan Kleinert
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
#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/parameterset.h"
#include "memilio/compartments/simulation.h"
#include "gtest/gtest.h"

TEST(TestCompartmentalModel, secir)
{
    /********************************
     * Define population categories *
     ********************************/

    // Note: There can be an arbitrary numober of categories
    // Each category is a dimension of a multidimensional array
    // defining the populatons of each individual compartment

    enum class InfectionType
    {
        S,
        E,
        C,
        I,
        H,
        U,
        R,
        D,
        Count = 8

    };

    enum class AgeGroup
    {
        leq20,
        leq40,
        leq60,
        leq80,
        o80,
        Count = 5
    };

    enum class Gender
    {
        Female,
        Male,
        Divers,
        Count = 3
    };

    enum class Income
    {
        poor,
        rich,
        Count = 2
    };

    using Po = mio::Populations<InfectionType, AgeGroup, Gender, Income>;

    /***********************************************
     * Define parameters and instantiate the model *
     ***********************************************/

    // Note: There can be an arbitrary amount of parameters
    // These are internally added to a compile time map for fast lookup

    struct IncubationTime {
        using Type = ScalarType;
        static constexpr Type get_default()
        {
            return 1.0;
        }
        static constexpr const char* name()
        {
            return "IncubationTime";
        }
    };

    struct SerialInterval {
        using Type = ScalarType;
        static constexpr Type get_default()
        {
            return 1.0;
        }
        static constexpr const char* name()
        {
            return "SerialInterval";
        }
    };

    //ADD MORE PARAMETERS HERE

    using Pa = mio::ParameterSet<IncubationTime, SerialInterval>;

    mio::CompartmentalModel<InfectionType, Po, Pa> model;

    /********************
     * Define the flows *
     ********************/

    for (size_t i = 0; i < static_cast<size_t>(AgeGroup::Count); ++i) {
        AgeGroup ai = static_cast<AgeGroup>(i);
        for (size_t j = 0; j < static_cast<size_t>(Gender::Count); ++j) {
            Gender gj = static_cast<Gender>(j);
            for (size_t k = 0; k < static_cast<size_t>(Income::Count); ++k) {
                Income ik = static_cast<Income>(k);

                //Ei to Ci
                model.add_flow(std::make_tuple(InfectionType::S, ai, gj, ik),
                               std::make_tuple(InfectionType::E, ai, gj, ik),
                               [ai, gj, ik](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                                   return Po::get_from(y, InfectionType::E, ai, gj, ik) /
                                          (2 * p.get<SerialInterval>() - p.get<IncubationTime>());
                               });
            }
        }
    }

    /****************************
     * Define initial conditios *
     ****************************/

    for (size_t i = 0; i < static_cast<size_t>(AgeGroup::Count); ++i) {
        AgeGroup ai = static_cast<AgeGroup>(i);
        for (size_t j = 0; j < static_cast<size_t>(Gender::Count); ++j) {
            Gender gj = static_cast<Gender>(j);
            for (size_t k = 0; k < static_cast<size_t>(Income::Count); ++k) {
                Income ik = static_cast<Income>(k);

                model.populations.set(9750, InfectionType::S, ai, gj, ik);
                model.populations.set(100, InfectionType::E, ai, gj, ik);
                model.populations.set(100, InfectionType::C, ai, gj, ik);
                model.populations.set(50, InfectionType::I, ai, gj, ik);
                // all other populations are zero initialized
            }
        }
    }

    /********************************************************
     *  Do some simulations with different incubation times *
     ********************************************************/

    double t0 = 0, tmax = 10, dt = 0.01;

    std::vector<double> inc_times{2., 3., 4.};
    std::vector<mio::TimeSeries<double>> results;

    for (auto inc_time : inc_times) {
        model.parameters.set<IncubationTime>(inc_time);
        results.push_back(simulate(t0, tmax, dt, model));
    }
}
