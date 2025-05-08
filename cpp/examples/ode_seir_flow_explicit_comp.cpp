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

using namespace mio::examples;
using FlowSim = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>;

int main()
{
    mio::set_log_level(mio::LogLevel::info);

    // Time parameters
    ScalarType t0   = 0.;
    ScalarType dt   = 0.5;
    ScalarType tmax = 100.;

    // Initialize flow-based SEIR model
    mio::oseir::Model<ScalarType> model_flow(1); // one age group
    const auto sus = 9700, exp = 100, inf = 100, rec = 100;
    model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = sus;
    model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = exp;
    model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = inf;
    model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}]   = rec;

    model_flow.parameters.set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model_flow.parameters.set<mio::oseir::TimeInfected<ScalarType>>(6.0);
    model_flow.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);

    mio::ContactMatrixGroup& contact_matrix = model_flow.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);
    // contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.0));
    model_flow.check_constraints();

    // Integrator for flow simulation
    auto integrator = std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
    auto integrator2 =
        std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
    // auto integrator  = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
    // auto integrator2 = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
    // const ScalarType max_dt_ref = 1e-1;
    // integrator->set_dt_max(max_dt_ref);

    // Run flow-based simulation
    FlowSim sim_flow(model_flow, t0, dt);
    sim_flow.set_integrator(integrator);
    sim_flow.advance(tmax);
    const auto& seir_res = sim_flow.get_result();

    // commuter state estimations
    const auto p_c   = 0.3;
    const double S_c = sus * p_c;
    const double E_c = exp * p_c;
    const double I_c = inf * p_c;
    const double R_c = rec * p_c;
    Eigen::VectorXd initial_mobile_pop(4);
    initial_mobile_pop << S_c, E_c, I_c, R_c;

    // high-precision reference solution for mobile population
    mio::log_info("Calculating high-precision reference solution for mobile population...");
    mio::TimeSeries<ScalarType> flow_commuter_results(4); // S_c, E_c, I_c, R_c
    flow_commuter_results.add_time_point(t0, initial_mobile_pop);
    Eigen::VectorXd mobile_pop_ref = initial_mobile_pop;

    for (ScalarType t = t0; t < tmax; t += dt) {
        if (t + dt > tmax) {
            break;
        }
        const auto closest_idx_total = static_cast<size_t>(t / dt);
        const auto& total_pop_at_t   = seir_res.get_value(closest_idx_total);
        flow_based_mobility_returns(mobile_pop_ref, sim_flow, total_pop_at_t, t, dt);
        flow_commuter_results.add_time_point(t + dt, mobile_pop_ref);
    }
    // auto flow_end = mobile_pop_ref;

    // Initialize explicit compartmental model
    ExplicitModel model_explicit(1);
    const ScalarType p_nc                                                                 = 1.0 - p_c;
    model_explicit.populations[{mio::AgeGroup(0), InfectionStateExplicit::S_NonCommuter}] = sus * p_nc;
    model_explicit.populations[{mio::AgeGroup(0), InfectionStateExplicit::E_NonCommuter}] = exp * p_nc;
    model_explicit.populations[{mio::AgeGroup(0), InfectionStateExplicit::I_NonCommuter}] = inf * p_nc;
    model_explicit.populations[{mio::AgeGroup(0), InfectionStateExplicit::R_NonCommuter}] = rec * p_nc;
    model_explicit.populations[{mio::AgeGroup(0), InfectionStateExplicit::S_Commuter}]    = sus * p_c;
    model_explicit.populations[{mio::AgeGroup(0), InfectionStateExplicit::E_Commuter}]    = exp * p_c;
    model_explicit.populations[{mio::AgeGroup(0), InfectionStateExplicit::I_Commuter}]    = inf * p_c;
    model_explicit.populations[{mio::AgeGroup(0), InfectionStateExplicit::R_Commuter}]    = rec * p_c;

    // Copy parameters from flow model
    model_explicit.parameters.set<mio::oseir::TimeExposed<ScalarType>>(
        model_flow.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(0)]);
    model_explicit.parameters.set<mio::oseir::TimeInfected<ScalarType>>(
        model_flow.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(0)]);
    model_explicit.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(
        model_flow.parameters.get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(0)]);
    model_explicit.parameters.get<mio::oseir::ContactPatterns<ScalarType>>() =
        model_flow.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();

    model_explicit.check_constraints();

    // Run explicit simulation
    ExplicitSim sim_explicit(model_explicit, t0, dt);
    sim_explicit.set_integrator(integrator2);
    sim_explicit.advance(tmax);
    const auto& explicit_results_full = sim_explicit.get_result();

    // Extract commuter compartments (last 4) from explicit results
    mio::TimeSeries<ScalarType> explicit_commuter_results(4);
    for (Eigen::Index i = 0; i < explicit_results_full.get_num_time_points(); ++i) {
        explicit_commuter_results.add_time_point(explicit_results_full.get_time(i),
                                                 explicit_results_full.get_value(i).tail(4));
    }

    // --- Write results to CSV file ---
    std::ofstream outfile("/localdata1/code/memilio/saves/comparison_data.csv");
    outfile << std::setprecision(15);
    outfile << "time,flow_S_c,flow_E_c,flow_I_c,flow_R_c,explicit_S_c,explicit_E_c,explicit_I_c,explicit_R_c\n";

    // Assuming both time series have the same time points (due to same dt and tmax)
    for (Eigen::Index i = 0; i < flow_commuter_results.get_num_time_points(); ++i) {
        outfile << flow_commuter_results.get_time(i);
        for (Eigen::Index j = 0; j < flow_commuter_results.get_num_elements(); ++j) {
            outfile << "," << flow_commuter_results.get_value(i)[j];
        }
        if (i < explicit_commuter_results.get_num_time_points() &&
            std::abs(explicit_commuter_results.get_time(i) - flow_commuter_results.get_time(i)) < 1e-10) {
            for (Eigen::Index j = 0; j < explicit_commuter_results.get_num_elements(); ++j) {
                outfile << "," << explicit_commuter_results.get_value(i)[j];
            }
        }
        else {
            outfile << ",,,,";
            mio::log_warning("Time point mismatch between flow ({}) and explicit ({}) results at index {}.",
                             flow_commuter_results.get_time(i), explicit_commuter_results.get_time(i), i);
        }
        outfile << "\n";
    }
    outfile.close();

    return 0;
}
