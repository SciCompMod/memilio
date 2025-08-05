/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Jan Kleinert, Daniel Abele, Rene Schmieding
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

#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/stochastic_simulation.h"
#include "memilio/config.h"
#include "memilio/math/integrator.h"
#include "memilio/utils/time_series.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

TEST(TestCompartmentSimulation, integrator_uses_model_reference)
{
    struct MockModel {
        Eigen::VectorXd get_initial_values() const
        {
            return Eigen::VectorXd::Zero(1);
        }
        void eval_right_hand_side(const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&,
                                  double, Eigen::Ref<Eigen::VectorXd> dydt) const
        {
            dydt[0] = this->m_dydt;
        }
        double m_dydt = 1.0;
    };

    auto sim = mio::Simulation<double, MockModel>(MockModel(), 0.0);
    sim.advance(1.0);

    ASSERT_NEAR(sim.get_result().get_last_value()[0], 1.0, 1e-5);

    //modifying the model from the outside should affect the integration result
    sim.get_model().m_dydt = 2.0;
    sim.advance(2.0);

    ASSERT_NEAR(sim.get_result().get_last_value()[0], 3.0, 1e-5);
}

struct MockSimulateSim { // looks just enough like a simulation for the simulate functions not to notice

    // this "model" converts to and from int implicitly, exposing its value after calling chech_constraints
    // this enables us to check whether check_constraints is called before the simulation is constructed
    struct Model {
        Model(int val_in)
            : hidden_val(val_in)
            , val(0)
        {
        }

        void check_constraints() const
        {
            val = hidden_val;
        }

        operator int() const
        {
            return val;
        }

        int hidden_val;
        mutable int val;
    };

    template <class ...Integrands>
    struct Core: public mio::IntegratorCore<double, Integrands...>
    {
        Core(int val_in)
            : mio::IntegratorCore<double, Integrands...>(val_in, 0)
        {
        }

        bool step(const Integrands&..., Eigen::Ref<const Eigen::VectorX<double>>, double&, double&,
                      Eigen::Ref<Eigen::VectorX<double>>) const override
        {
            return true;
        }
    };

    MockSimulateSim(int model_in, double t0_in, double dt_in)
    {
        model = model_in;
        t0    = t0_in;
        dt    = dt_in;
    }

    template <class... Integrands>
    void set_integrator(std::unique_ptr<mio::IntegratorCore<double, Integrands...>> integrator_in)
    {
        integrator = (int)integrator_in->get_dt_min();
    }

    auto get_result()
    {
        // basically, return 17
        mio::TimeSeries<double> ts(0);
        ts.add_time_point(17);
        return ts;
    }

    auto get_flows()
    {
        return get_result();
    }

    void advance(double tmax_in) // wrong return type, but simulate functions do not use it anyways
    {
        tmax = tmax_in;
    }

    static void clear() // reset variables
    {
        t0 = dt = tmax = model = integrator = 0;
    }

    // static public members used to check whether a simulate function works as expected.
    inline static double t0, dt, tmax;
    inline static int model, integrator;
};

TEST(TestCompartmentSimulation, simulate_functions)
{
    // this checks that the (not model-specific) simulate functions make all calls as expected, like "advance(tmax)"

    // this works by misusing a simulate function on the MockSimulateSim to write out 1,2,3,4,5 and return 17

    // the lambdas in this test help deal with the integrator pointer and TimeSeries used and returned by a simulate
    // function, by misusing these types to store simple values

    const double t0 = 1, dt = 2, tmax = 3;
    const int model = 4, integrator = 5;

    using Sim = MockSimulateSim;

    // evaluate the timeseries
    // we expect 17 as result stored as the first time, but just in case we also check for the number of time points
    const auto eval_ts = [](auto&& ts) {
        return ts.get_num_time_points() == 0 ? 0 : ts.get_num_time_points() * ts.get_time(0);
    };
    // this handles flows (or any Vector of TimeSeries) as well, using an average to immitate a logical "and"
    const auto eval_vts = [](auto&& vts) {
        double sum = 0;
        for (auto& ts : vts) {
            sum += ts.get_num_time_points() == 0 ? 0 : ts.get_num_time_points() * ts.get_time(0);
        }
        return sum / vts.size();
    };

    // helpers to deal with different orders of cores. do not reuse this or something similar in actual code
    const auto evil_pointer_cast_1 = [](int i) {
        return std::make_unique<Sim::Core<mio::DerivFunction<double>>>(i);
    };
    const auto evil_pointer_cast_2 = [](int i) {
        return std::make_unique<Sim::Core<mio::DerivFunction<double>, mio::DerivFunction<double>>>(i);
    };

    // helper function to compose a "simulate" function
    const auto compose = [](auto&& f, auto&& g, auto&& h) {
        return [f, g, h](double t0_, double tmax_, double dt_, int model_, int i_) {
            return f(g(t0_, tmax_, dt_, model_, h(i_)));
        };
    };

    // list of functions to test, with helpers attached
    // note that std::functions only work with lambdas that do not use captures
    std::vector<std::function<double(double, double, double, int, int)>> simulate_fcts = {
        compose(eval_ts, mio::simulate<double, Sim::Model, Sim>, evil_pointer_cast_1),
        compose(eval_vts, mio::simulate_flows<double, Sim::Model, Sim>, evil_pointer_cast_1),
        compose(eval_ts, mio::simulate_stochastic<double, Sim::Model, Sim>, evil_pointer_cast_2)};

    // test all simulate functions
    for (auto&& simulate : simulate_fcts) {
        // reset static values in the mock, then call simulate to set them
        Sim::clear();
        auto result = simulate(t0, tmax, dt, model, integrator);
        // check that all values are set correctly
        EXPECT_NEAR(result, 17.0, mio::Limits<double>::zero_tolerance());
        // use equal (instead of near) since there is no math happening in this test
        EXPECT_EQ(Sim::t0, t0);
        EXPECT_EQ(Sim::dt, dt);
        EXPECT_EQ(Sim::tmax, tmax);
        EXPECT_EQ(Sim::model, model);
        EXPECT_EQ(Sim::integrator, integrator);
    }
}
