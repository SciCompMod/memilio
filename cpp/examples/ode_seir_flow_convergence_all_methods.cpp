#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/math/euler.h"
#include "memilio/math/rk2.h"
#include "memilio/math/rk3.h"
#include "memilio/utils/type_list.h"
#include "state_estimators.h"
#include "memilio/data/analyze_result.h"
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

using namespace mio::examples;
using FlowSim = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>;

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    // Test different step sizes
    std::vector<ScalarType> dts = {2.0, 1.0, 0.5, 0.25, 0.125, 0.0625};

    // Reference solution with very small fixed dt
    const ScalarType dt_ref = 1e-5;

    // Methods to test: Euler, RK2, RK3, RK4
    std::vector<std::string> method_names = {"euler", "rk2", "rk3", "rk4"};

    // Common setup
    ScalarType t0   = 0.;
    ScalarType tmax = 100.;

    // Initialize model
    mio::oseir::Model<ScalarType> model(1);
    const double sus = 9700, exp = 100, inf = 100, rec = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = sus;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = exp;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = inf;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}]   = rec;

    model.parameters.set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<ScalarType>>(6.0);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);
    model.check_constraints();

    // Commuter initial conditions
    const auto p_c   = 0.3;
    const double S_c = sus * p_c * 1.0;
    const double E_c = exp * p_c * 0.5;
    const double I_c = inf * p_c * 0.8;
    const double R_c = rec * p_c * 0.1;
    Eigen::VectorXd initial_commuter(4);
    initial_commuter << S_c, E_c, I_c, R_c;

    // ========== Compute Reference Solution (with RK4 at dt_ref) ==========
    std::cout << "Computing reference solution with dt=" << dt_ref << "..." << std::endl;

    auto integrator_ref =
        std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_fehlberg78>>();
    FlowSim sim_ref(model, t0, dt_ref);
    sim_ref.set_integrator(integrator_ref);
    sim_ref.advance(tmax);
    const auto& ref_totals = sim_ref.get_result();

    // Compute reference commuter trajectory
    // RK4 integrates totals internally, so we track totals independently
    mio::TimeSeries<ScalarType> ref_commuter_results(4);
    ref_commuter_results.add_time_point(t0, initial_commuter);
    Eigen::VectorXd commuter_ref = initial_commuter;
    Eigen::VectorXd totals_ref   = model.populations.get_compartments();

    for (size_t i = 0; i < static_cast<size_t>(ref_totals.get_num_time_points()) - 1; ++i) {
        ScalarType t_curr  = ref_totals.get_time(i);
        ScalarType t_next  = ref_totals.get_time(i + 1);
        ScalarType dt_step = t_next - t_curr;

        if (t_curr >= tmax)
            break;

        // RK4 updates both commuter_ref and totals_ref
        mio::examples::flow_based_mobility_returns_rk4(commuter_ref, totals_ref, model, t_curr, dt_step);
        ref_commuter_results.add_time_point(t_next, commuter_ref);
    }

    Eigen::VectorXd ref_final = ref_commuter_results.get_last_value();
    std::cout << "Reference final state: S=" << ref_final[0] << ", E=" << ref_final[1] << ", I=" << ref_final[2]
              << ", R=" << ref_final[3] << std::endl;

    // ========== Test Each Method at Different Step Sizes ==========
    const std::string save_dir = "/localdata1/code/memilio/saves";

    for (const auto& method_name : method_names) {
        std::cout << "\n========== Testing Method: " << method_name << " ==========" << std::endl;

        for (auto dt : dts) {
            std::cout << "  dt = " << dt << "..." << std::endl;

            // Setup integrator for totals simulation - use same order as commuter method
            std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator;
            if (method_name == "euler") {
                integrator = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
            }
            else if (method_name == "rk2") {
                integrator = std::make_shared<mio::RK2IntegratorCore<ScalarType>>();
            }
            else if (method_name == "rk3") {
                integrator = std::make_shared<mio::RK3IntegratorCore<ScalarType>>();
            }
            else if (method_name == "rk4") {
                integrator =
                    std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            }

            FlowSim sim(model, t0, dt);
            sim.set_integrator(integrator);
            sim.advance(tmax);
            const auto& totals_result = sim.get_result();

            // Compute commuter trajectory with selected method
            mio::TimeSeries<ScalarType> commuter_results(4);
            commuter_results.add_time_point(t0, initial_commuter);
            Eigen::VectorXd commuter_current = initial_commuter;

            // For RK2/RK3/RK4: track totals independently (integrated by the RK function)
            // For Euler: use simulator totals
            Eigen::VectorXd totals_current = model.populations.get_compartments();

            for (size_t i = 0; i < static_cast<size_t>(totals_result.get_num_time_points()) - 1; ++i) {
                ScalarType t_curr  = totals_result.get_time(i);
                ScalarType t_next  = totals_result.get_time(i + 1);
                ScalarType dt_step = t_next - t_curr;

                if (t_curr >= tmax)
                    break;

                // Apply the appropriate method
                if (method_name == "euler") {
                    // Euler uses flows from simulator
                    mio::examples::flow_based_mobility_returns(commuter_current, sim, totals_result.get_value(i),
                                                               t_curr, dt_step);
                }
                else if (method_name == "rk2") {
                    // RK2 integrates totals internally - use and update totals_current
                    mio::examples::flow_based_mobility_returns_rk2(commuter_current, totals_current, model, t_curr,
                                                                   dt_step);
                }
                else if (method_name == "rk3") {
                    // RK3 integrates totals internally - use and update totals_current
                    mio::examples::flow_based_mobility_returns_rk3(commuter_current, totals_current, model, t_curr,
                                                                   dt_step);
                }
                else if (method_name == "rk4") {
                    // RK4 integrates totals internally - use and update totals_current
                    mio::examples::flow_based_mobility_returns_rk4(commuter_current, totals_current, model, t_curr,
                                                                   dt_step);
                }

                commuter_results.add_time_point(t_next, commuter_current);
            }

            // Compute error against reference
            Eigen::VectorXd final_state = commuter_results.get_last_value();
            Eigen::VectorXd error       = final_state - ref_final;

            ScalarType max_abs_error = error.cwiseAbs().maxCoeff();
            ScalarType max_rel_error = 0.0;
            for (int i = 0; i < 4; ++i) {
                if (std::abs(ref_final[i]) > 1e-10) {
                    ScalarType rel_err = std::abs(error[i]) / std::abs(ref_final[i]);
                    max_rel_error      = std::max(max_rel_error, rel_err);
                }
            }

            std::cout << "    Max abs error: " << max_abs_error << ", Max rel error: " << max_rel_error << std::endl;

            // Save to CSV
            std::ostringstream dt_ss;
            dt_ss.setf(std::ios::fixed);
            dt_ss << std::setprecision(4) << dt;
            std::string dt_str = dt_ss.str();

            std::string filename = save_dir + "/convergence_" + method_name + "_dt_" + dt_str + ".csv";
            std::ofstream out(filename);
            out << std::setprecision(20);
            out << "dt,max_abs_error,max_rel_error,err_S,err_E,err_I,err_R\n";
            out << dt << "," << max_abs_error << "," << max_rel_error << "," << std::abs(error[0]) << ","
                << std::abs(error[1]) << "," << std::abs(error[2]) << "," << std::abs(error[3]) << "\n";
            out.close();
        }
    }

    std::cout << "\n========== All methods completed ==========" << std::endl;
    return 0;
}
