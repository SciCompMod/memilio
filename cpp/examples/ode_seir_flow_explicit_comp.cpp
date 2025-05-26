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

using namespace mio::examples;
using FlowSim = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>;

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    // Time parameters
    ScalarType t0   = 0.;
    ScalarType dt   = 0.5;
    ScalarType tmax = 10.;

    // Initialize flow-based SEIR model
    mio::oseir::Model<ScalarType> model_flow(1); // one age group
    const double sus = 9700, exp = 100, inf = 100, rec = 100;
    model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = sus;
    model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = exp;
    model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = inf;
    model_flow.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}]   = rec;

    model_flow.parameters.set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model_flow.parameters.set<mio::oseir::TimeInfected<ScalarType>>(6.0);
    model_flow.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);

    mio::ContactMatrixGroup& contact_matrix = model_flow.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);
    model_flow.check_constraints();

    // Integrator for flow simulation
    auto integrator_flow =
        std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();

    // Run flow-based simulation (aggregate model)
    FlowSim sim_flow(model_flow, t0, dt); // dt here is the desired output step, integrator might use smaller steps
    sim_flow.set_integrator(integrator_flow);
    sim_flow.advance(tmax);
    const auto& seir_res = sim_flow.get_result();

    // Commuter state estimations (sub-population for flow-based approach)
    const auto p_c   = 0.3; // Fraction of commuters
    const double S_c = sus * p_c;
    const double E_c = exp * p_c;
    const double I_c = inf * p_c;
    const double R_c = rec * p_c;
    Eigen::VectorXd initial_mobile_pop(4);
    initial_mobile_pop << S_c, E_c, I_c, R_c;

    // Calculate flow-based commuter results
    mio::log_info("Calculating flow-based commuter results...");
    mio::TimeSeries<ScalarType> flow_commuter_results(4);
    flow_commuter_results.add_time_point(t0, initial_mobile_pop);
    Eigen::VectorXd mobile_pop_current = initial_mobile_pop;

    for (size_t i_res = 0; i_res < static_cast<size_t>(seir_res.get_num_time_points()) - 1; ++i_res) {
        ScalarType current_t_res  = seir_res.get_time(i_res);
        ScalarType next_t_res     = seir_res.get_time(i_res + 1);
        ScalarType current_dt_res = next_t_res - current_t_res;

        if (current_t_res >= tmax) {
            break;
        }

        const auto& total_pop_at_t = seir_res.get_value(i_res); // Population at the start of the interval

        // Use flow_based_mobility_returns
        // It modifies mobile_pop_current in place
        mio::examples::flow_based_mobility_returns(mobile_pop_current, sim_flow, total_pop_at_t, current_t_res,
                                                   current_dt_res);

        flow_commuter_results.add_time_point(next_t_res, mobile_pop_current);
    }

    // Initialize explicit compartmental model
    // One age group, one commuter group (CommuterBase)
    ExplicitModel model_explicit(1, 1);
    const ScalarType p_nc = 1.0 - p_c;

    // NonCommuter
    model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::S}] = sus * p_nc;
    model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::E}] = exp * p_nc;
    model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::I}] = inf * p_nc;
    model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::R}] = rec * p_nc;

    // Commuter (using CommuterBase for the single commuter group)
    model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::S}] = sus * p_c;
    model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::E}] = exp * p_c;
    model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::I}] = inf * p_c;
    model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::R}] = rec * p_c;

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

    // Integrator for explicit simulation
    auto integrator_explicit =
        std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();

    // Run explicit simulation
    ExplicitSim sim_explicit(model_explicit, t0, dt); // dt is the desired output step
    sim_explicit.set_integrator(integrator_explicit);
    sim_explicit.advance(tmax);
    const auto& explicit_results_full = sim_explicit.get_result();

    // Extract commuter compartments from explicit results
    // For 1 age group, NonCommuter (4 states), CommuterBase (4 states) -> CommuterBase is tail(4)
    mio::TimeSeries<ScalarType> explicit_commuter_results(4);
    for (Eigen::Index i = 0; i < explicit_results_full.get_num_time_points(); ++i) {
        explicit_commuter_results.add_time_point(explicit_results_full.get_time(i),
                                                 explicit_results_full.get_value(i).tail(4));
    }

    // --- Print results to console ---
    std::cout << std::fixed << std::setprecision(6); // Default precision for most output
    mio::log_info("Comparing Flow-based Commuter vs Explicit Commuter results:");
    std::cout << "Time      | Flow S_c  E_c  I_c  R_c             | Explicit S_c  E_c  I_c  R_c         | Diff S_c E_c "
                 "I_c R_c (Precision 13)"
              << std::endl;
    std::cout << "----------|-----------------------------------|-----------------------------------|------------------"
                 "------------------------"
              << std::endl;

    ScalarType max_abs_diff_norm = 0.0;

    for (Eigen::Index i = 0; i < flow_commuter_results.get_num_time_points(); ++i) {
        ScalarType current_time       = flow_commuter_results.get_time(i);
        Eigen::VectorXd flow_vals     = flow_commuter_results.get_value(i);
        Eigen::VectorXd explicit_vals = explicit_commuter_results.get_value(i);

        std::cout << std::setw(9) << current_time << " |";
        for (Eigen::Index k = 0; k < flow_vals.size(); ++k) {
            std::cout << std::setw(7) << flow_vals[k] << (k == flow_vals.size() - 1 ? "" : " ");
        }
        std::cout << " |";
        for (Eigen::Index k = 0; k < explicit_vals.size(); ++k) {
            std::cout << std::setw(7) << explicit_vals[k] << (k == explicit_vals.size() - 1 ? "" : " ");
        }
        std::cout << " | ";

        Eigen::VectorXd diff = flow_vals - explicit_vals;
        std::cout << std::fixed << std::setprecision(13);
        for (Eigen::Index k = 0; k < diff.size(); ++k) {
            std::cout << std::scientific << std::setw(20) << diff[k] << (k == diff.size() - 1 ? "" : " ");
        }
        std::cout << std::fixed << std::setprecision(6);
        max_abs_diff_norm = std::max(max_abs_diff_norm, diff.lpNorm<Eigen::Infinity>());
        std::cout << std::endl;
    }
    return 0;
}
