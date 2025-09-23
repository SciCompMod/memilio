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
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/utils/type_list.h"
#include "state_estimators.h"
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace mio::examples;
using FlowSim = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>;

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    std::vector<ScalarType> dts = {1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5};
    std::vector<ScalarType> pcs = {0.95, 0.96, 0.97, 0.98, 0.99, 1.00};

    const std::string save_dir         = "/localdata1/code/memilio/saves";
    const std::string metrics_csv_path = save_dir + "/commuter_metrics.csv";
    std::ofstream metrics_csv(metrics_csv_path, std::ios::out);
    metrics_csv << std::setprecision(15);

    // Aggregierte Metriken je (p_c, dt)
    metrics_csv << "p_c,dt,n_steps,"
                << "viol_frac,neg_frac,over_S_frac,over_E_frac,over_I_frac,over_R_frac,"
                << "Linf_E,Linf_I,Linf_R,Linf_all,"
                << "rmax_dt_max,rmax_dt_mean\n";

    for (auto p_c : pcs) {
        std::cout << "Starting p_c = " << p_c << " ..." << std::endl;
        for (auto curr_dt : dts) {

            // --- Zeithorizont & Modelle ---
            ScalarType t0   = 0.;
            ScalarType dt   = curr_dt; // Euler-Schrittweite und Ausgabeabstand
            ScalarType tmax = 150.;

            // Flow/Total SEIR (1 Altersgruppe)
            mio::oseir::Model<ScalarType> model_flow(1);
            const double sus = 7700, exp = 1500, inf = 1500, rec = 300;
            model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = sus;
            model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = exp;
            model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = inf;
            model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}]   = rec;

            model_flow.parameters.set<mio::oseir::TimeExposed<ScalarType>>(1.0);
            model_flow.parameters.set<mio::oseir::TimeInfected<ScalarType>>(1.0);
            model_flow.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);

            mio::ContactMatrixGroup& contact_matrix =
                model_flow.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
            contact_matrix[0].get_baseline().setConstant(2.7);
            model_flow.check_constraints();

            // Integrator & Simulation (Totals)
            auto integrator_flow =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            FlowSim sim_flow(model_flow, t0, dt); // Ausgabeabstand = dt
            sim_flow.set_integrator(integrator_flow);
            sim_flow.advance(tmax);
            const auto& seir_res = sim_flow.get_result();

            // Parameter (für rmax_dt Monitoring)
            const double rho =
                model_flow.parameters.get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(0)];
            const auto phi_ct = 2.7;
            const double TE   = model_flow.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(0)];
            const double TI   = model_flow.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(0)];

            // Initiale Pendler-Population (hier ohne Multiplikatoren < 1, sauberer für Heatmap)
            const double S_c = sus * p_c;
            const double E_c = exp * p_c;
            const double I_c = inf * p_c;
            const double R_c = rec * p_c;
            Eigen::VectorXd initial_mobile_pop(4);
            initial_mobile_pop << S_c, E_c, I_c, R_c;

            // --- Explizites Kommuter-Modell (Referenz) ---
            ExplicitModel model_explicit(1, 1);

            model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
                mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::S}] = sus - S_c;
            model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
                mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::E}] = exp - E_c;
            model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
                mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::I}] = inf - I_c;
            model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
                mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::R}] = rec - R_c;

            model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
                mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::S}] = S_c;
            model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
                mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::E}] = E_c;
            model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
                mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::I}] = I_c;
            model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
                mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::R}] = R_c;

            // Parameter übernehmen
            model_explicit.parameters.set<mio::oseir::TimeExposed<ScalarType>>(
                model_flow.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(0)]);
            model_explicit.parameters.set<mio::oseir::TimeInfected<ScalarType>>(
                model_flow.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(0)]);
            model_explicit.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(
                model_flow.parameters
                    .get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(0)]);
            model_explicit.parameters.get<mio::oseir::ContactPatterns<ScalarType>>() =
                model_flow.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
            model_explicit.check_constraints();

            auto integrator_explicit =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            ExplicitSim sim_explicit(model_explicit, t0, dt);
            sim_explicit.set_integrator(integrator_explicit);
            sim_explicit.advance(tmax);
            const auto& explicit_results_full = sim_explicit.get_result();

            // Pendler aus explizitem Modell extrahieren (letzte 4 Spalten)
            mio::TimeSeries<ScalarType> explicit_commuter_results(4);
            for (Eigen::Index i = 0; i < explicit_results_full.get_num_time_points(); ++i) {
                explicit_commuter_results.add_time_point(explicit_results_full.get_time(i),
                                                         explicit_results_full.get_value(i).tail(4));
            }

            // --- Euler-Auxiliary-Step (Test-Methode) ---
            const auto euler_commuter_results = calculate_mobile_population(
                t0, tmax, dt, integrate_mobile_population_euler, sim_flow, seir_res, initial_mobile_pop);

            // --- Aggregation der Metriken über die Zeitpunkte ---
            const auto n_points = static_cast<size_t>(std::min<ptrdiff_t>(
                euler_commuter_results.get_num_time_points(), explicit_commuter_results.get_num_time_points()));

            size_t steps = 0, cnt_violation = 0, cnt_neg = 0, cnt_overS = 0, cnt_overE = 0, cnt_overI = 0,
                   cnt_overR = 0;
            double LinfE = 0.0, LinfI = 0.0, LinfR = 0.0;

            double rmax_dt_max = 0.0, rmax_dt_sum = 0.0;

            for (size_t i = 0; i < n_points; ++i) {
                // const auto time   = euler_commuter_results.get_time(static_cast<Eigen::Index>(i));
                const auto eul_v  = euler_commuter_results.get_value(static_cast<Eigen::Index>(i));
                const auto expl_v = explicit_commuter_results.get_value(static_cast<Eigen::Index>(i));
                const auto tot_v  = seir_res.get_value(static_cast<Eigen::Index>(i)); // [S,E,I,R] totals

                // Totals + Rateparameter (für rmax_dt)
                const double S_tot = tot_v[0], E_tot = tot_v[1], I_tot = tot_v[2], R_tot = tot_v[3];
                const double N_tot  = S_tot + E_tot + I_tot + R_tot;
                const double lambda = (N_tot > 0.) ? (rho * phi_ct * (I_tot / N_tot)) : 0.0;
                const double rS = lambda, rE = 1.0 / TE, rI = 1.0 / TI;
                const double rmax_dt = std::max({rS, rE, rI}) * static_cast<double>(dt);
                rmax_dt_max          = std::max(rmax_dt_max, rmax_dt);
                rmax_dt_sum += rmax_dt;

                // Linf-Fehler (E, I, R)
                LinfE = std::max(LinfE, std::abs(eul_v[1] - expl_v[1]));
                LinfI = std::max(LinfI, std::abs(eul_v[2] - expl_v[2]));
                LinfR = std::max(LinfR, std::abs(eul_v[3] - expl_v[3]));

                // Feasibility-Checks
                const bool neg_euler = (eul_v[0] < 0.) || (eul_v[1] < 0.) || (eul_v[2] < 0.) || (eul_v[3] < 0.);
                const bool over_S    = (eul_v[0] > S_tot);
                const bool over_E    = (eul_v[1] > E_tot);
                const bool over_I    = (eul_v[2] > I_tot);
                const bool over_R    = (eul_v[3] > R_tot);
                const bool any_bad   = neg_euler || over_S || over_E || over_I || over_R;

                if (neg_euler)
                    ++cnt_neg;
                if (over_S)
                    ++cnt_overS;
                if (over_E)
                    ++cnt_overE;
                if (over_I)
                    ++cnt_overI;
                if (over_R)
                    ++cnt_overR;
                if (any_bad)
                    ++cnt_violation;
                ++steps;
            }

            steps--;

            const double viol_frac    = (steps > 0) ? static_cast<double>(cnt_violation) / steps : 0.0;
            const double neg_frac     = (steps > 0) ? static_cast<double>(cnt_neg) / steps : 0.0;
            const double overS_frac   = (steps > 0) ? static_cast<double>(cnt_overS) / steps : 0.0;
            const double overE_frac   = (steps > 0) ? static_cast<double>(cnt_overE) / steps : 0.0;
            const double overI_frac   = (steps > 0) ? static_cast<double>(cnt_overI) / steps : 0.0;
            const double overR_frac   = (steps > 0) ? static_cast<double>(cnt_overR) / steps : 0.0;
            const double Linf_all     = std::max({LinfE, LinfI, LinfR});
            const double rmax_dt_mean = (steps > 0) ? (rmax_dt_sum / steps) : 0.0;

            metrics_csv << p_c << "," << curr_dt << "," << steps << "," << viol_frac << "," << neg_frac << ","
                        << overS_frac << "," << overE_frac << "," << overI_frac << "," << overR_frac << "," << LinfE
                        << "," << LinfI << "," << LinfR << "," << Linf_all << "," << rmax_dt_max << "," << rmax_dt_mean
                        << "\n";
        }
    }
    return 0;
}
