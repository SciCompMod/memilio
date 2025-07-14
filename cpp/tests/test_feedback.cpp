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
#include "load_test_data.h"
#include "memilio/compartments/compartmental_model.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/feedback_simulation.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/model.h"
#include "ode_seir/parameters.h"

#include "gtest/gtest.h"

class TestFeedbackSimulation : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create an ODE SEIR model with one age group.
        auto model                                                                   = mio::oseir::Model<double>(1);
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 1000;
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 1000;
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 1000;
        model.populations.set_difference_from_total({mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}, 10000);
        model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(1.0);
        model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
        model.parameters.set<mio::oseir::TimeInfected<double>>(2);

        mio::ContactMatrixGroup<double>& contact_matrix =
            model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
        contact_matrix[0].get_baseline().setConstant(2.7);
        contact_matrix[0].add_damping(0.6, mio::SimulationTime<double>(12.5));

        // ICU compartment index
        icu_indices = {2};

        // Create a simulation based on the model.
        auto sim = mio::Simulation<double, mio::oseir::Model<double>>(model);

        // Create the FeedbackSimulation by moving the simulation and providing the ICU indices.
        feedback_sim =
            std::make_unique<mio::FeedbackSimulation<double, mio::Simulation<double, mio::oseir::Model<double>>,
                                                     mio::oseir::ContactPatterns<double>>>(std::move(sim), icu_indices);
    }

    std::vector<size_t> icu_indices;
    std::unique_ptr<mio::FeedbackSimulation<double, mio::Simulation<double, mio::oseir::Model<double>>,
                                            mio::oseir::ContactPatterns<double>>>
        feedback_sim;
};

// Test that advancing the simulation adds ICU occupancy and perceived risk entries.
TEST_F(TestFeedbackSimulation, AdvanceAddICUAndRisk)
{
    // initially, the perceived risk time series should be empty
    EXPECT_EQ(feedback_sim->get_perceived_risk().get_num_time_points(), 0);

    // advance the simulation to t = 2.0 in two steps (dt_feedback = 1.0).
    feedback_sim->advance(2.0, 1.0);

    // ICU occupancy time series should now have 2 time points.
    const auto& icu_occ = feedback_sim->get_parameters().template get<mio::ICUOccupancyHistory<double>>();
    EXPECT_EQ(icu_occ.get_num_time_points(), 2);

    // similarly, the perceived risk time series should also have 2 entries.
    EXPECT_EQ(feedback_sim->get_perceived_risk().get_num_time_points(), 2);
}

TEST_F(TestFeedbackSimulation, CalcPerceivedRisk)
{
    // set GammaShapeParameter to 1 and GammaScaleParameter to 1 so that gamma becomes exp(-day)
    auto& fb_params                                            = feedback_sim->get_parameters();
    fb_params.template get<mio::GammaShapeParameter<double>>() = 1;
    fb_params.template get<mio::GammaScaleParameter<double>>() = 1;
    fb_params.template get<mio::NominalICUCapacity<double>>()  = 10;

    // add a single ICU occupancy time points with value 2.
    Eigen::VectorXd icu_value(1);
    icu_value << 2;
    auto& icu_occ = fb_params.template get<mio::ICUOccupancyHistory<double>>();
    icu_occ.add_time_point(0.0, icu_value);
    double risk = feedback_sim->calc_risk_perceived();
    // For day 0, we have gamma = exp(0) = 1 and therefore perc_risk = 2 / 10 = 0.2.
    EXPECT_NEAR(risk, 0.2, 1e-10);

    // Add another time point with value 2.
    icu_occ.add_time_point(1.0, icu_value);

    // Here, we have gamma = exp(-1) = 0.367879441 and therefore expect a risk of 0.2 + 0.07357588823428847.
    risk = feedback_sim->calc_risk_perceived();
    EXPECT_NEAR(risk, 0.27357588823428847, 1e-10);
}

TEST_F(TestFeedbackSimulation, ApplyFeedback)
{
    // bounds for contact reduction measures
    auto& feedback_params                                            = feedback_sim->get_parameters();
    feedback_params.template get<mio::ContactReductionMax<double>>() = std::vector<double>{0.8};
    feedback_params.template get<mio::ContactReductionMin<double>>() = std::vector<double>{0.2};

    // get initial number of damping samples.
    auto& contact_patterns  = feedback_sim->get_model().parameters.get<mio::oseir::ContactPatterns<double>>();
    size_t initial_dampings = contact_patterns.get_dampings().size();

    // set all historical ICU occupancy values to 100.0, to have maximum perceived risk.
    auto& icu_occ = feedback_params.template get<mio::ICUOccupancyHistory<double>>();
    Eigen::VectorXd icu_value(1);
    icu_value << 100.0;
    for (int t = -45; t <= 0; ++t) {
        icu_occ.add_time_point(t, icu_value);
    }

    // apply_feedback at time t = 0.
    feedback_sim->apply_feedback(0.0);

    // number of new damping samples added should equal the number of locations,
    // as given by the size of the ContactReductionMax vector.
    size_t num_locations = feedback_sim->get_parameters().template get<mio::ContactReductionMax<double>>().size();
    EXPECT_EQ(contact_patterns.get_dampings().size(), initial_dampings + num_locations);

    // contact reduction should be near the maximum allowed contact reduction.
    EXPECT_NEAR(contact_patterns.get_dampings().back().get_value().value(),
                feedback_params.template get<mio::ContactReductionMax<double>>()[0], 1e-3);
}

TEST_F(TestFeedbackSimulation, AddICUOccupancy)
{
    // check that the ICU occupancy time series is empty
    auto& feedback_params = feedback_sim->get_parameters();
    EXPECT_EQ(feedback_params.template get<mio::ICUOccupancyHistory<double>>().get_num_time_points(), 0);

    // add ICU occupancy for t = 1.0.
    feedback_sim->add_icu_occupancy(1.0);

    // check that the ICU occupancy time series now has one time point.
    EXPECT_EQ(feedback_params.template get<mio::ICUOccupancyHistory<double>>().get_num_time_points(), 1);
    // check that stored value is equal to the inital value relative to 100,000
    auto stored_val = feedback_params.template get<mio::ICUOccupancyHistory<double>>().get_value(0);
    auto expected_val =
        feedback_sim->get_model().populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] /
        feedback_sim->get_model().populations.get_total() * 100000;
    EXPECT_EQ(stored_val[0], expected_val);
}
