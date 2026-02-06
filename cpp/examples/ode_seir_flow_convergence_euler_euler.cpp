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
#include "memilio/data/analyze_result.h"
#include <iomanip>
#include <vector>
#include <cmath>

using namespace mio::examples;
using FlowSim = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>;

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    std::vector<ScalarType> dts = {2.0, 1.0, 0.5, 0.25, 0.125, 0.0001}; // 0.001 for reference solution

    for (auto curr_dt : dts) {

        // Time parameters
        ScalarType t0   = 0.;
        ScalarType dt   = curr_dt; // use the current dt from the loop
        ScalarType tmax = 100.;

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

        // Integrator for flow simulation (Euler)
        auto integrator_flow = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();

        // Run flow-based simulation (aggregate model)
        FlowSim sim_flow(model_flow, t0, dt);
        sim_flow.set_integrator(integrator_flow);
        sim_flow.advance(tmax);
        const auto& seir_res = sim_flow.get_result();

        // Commuter state estimations (sub-population for flow-based approach)
        const auto p_c   = 0.3; // Fraction of commuters
        const double S_c = sus * p_c * 1.0;
        const double E_c = exp * p_c * 0.5;
        const double I_c = inf * p_c * 0.8;
        const double R_c = rec * p_c * 0.1;
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

            // flow_based_mobility_returns_exact(commuter, totals, model, t, dt)
            Eigen::VectorXd totals_step = total_pop_at_t; // mutable copy for in-place update of totals
            // mio::examples::flow_based_mobility_returns_rk4(mobile_pop_current, totals_step, model_flow,
            //                                                     current_t_res, current_dt_res);
            mio::examples::flow_based_mobility_returns(mobile_pop_current, sim_flow, total_pop_at_t, current_t_res,
                                                       current_dt_res);

            flow_commuter_results.add_time_point(next_t_res, mobile_pop_current);
        }

        // Initialize explicit compartmental model
        // One age group, one commuter group (CommuterBase)
        ExplicitModel model_explicit(1, 1);

        // NonCommuter
        model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
            mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::S}] = sus - S_c;
        model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
            mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::E}] = exp - E_c;
        model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
            mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::I}] = inf - I_c;
        model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
            mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::R}] = rec - R_c;

        // Commuter (using CommuterBase for the single commuter group)
        model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
            mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::S}] = S_c;
        model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
            mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::E}] = E_c;
        model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
            mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::I}] = I_c;
        model_explicit.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
            mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::R}] = R_c;

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

        // Integrator for explicit simulation (Euler)
        auto integrator_explicit = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();

        // Run explicit simulation
        ExplicitSim sim_explicit(model_explicit, t0, dt); // dt is the desired output step
        sim_explicit.set_integrator(integrator_explicit);
        sim_explicit.advance(tmax);
        const auto& explicit_results_full = mio::interpolate_simulation_result(sim_explicit.get_result());

        // print results to CSV file
        const std::string save_dir                        = "/localdata1/code/memilio/saves";
        const std::vector<std::string> compartment_labels = {"S", "E", "I", "R"};
        const auto precision_csv                          = 20;

        // add current dt to filenames
        std::ostringstream dt_ss;
        dt_ss.setf(std::ios::fixed);
        dt_ss << std::setprecision(3) << dt;
        const std::string dt_suffix = "_dt_" + dt_ss.str();

        explicit_results_full.export_csv(save_dir + "/explicit_commuter_compare_sol_euler" + dt_suffix + ".csv",
                                         compartment_labels, ',', precision_csv);
        flow_commuter_results.export_csv(save_dir + "/flow_commuter_compare_sol_euler" + dt_suffix + ".csv",
                                         compartment_labels, ',', precision_csv);
        // flow_commuter_results.export_csv(save_dir + "/flow_exact_commuter_compare_sol" + dt_suffix + ".csv",
        //                                  compartment_labels, ',', precision_csv);
    }
    return 0;
}
