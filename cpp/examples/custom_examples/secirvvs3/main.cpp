// #include <iostream>

// #include "IpTNLP.hpp"
// #include "IpIpoptApplication.hpp"

// #include "settings.h"
// #include "secirvvs_ipopt.h"

// #include "ad/ad.hpp"

// #include "memilio/utils/logging.h"

// // ============================================================================
// // --- Main Function ---------------------------------------------------------
// // ============================================================================

// int main()
// {
//     // Switch off logging for mio
//     mio::set_log_level(mio::LogLevel::off);

//     // // Set up problem settings
//     // ProblemSettings problem;

//     // // Create NLP and solver
//     // Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new Seair_NLP(problem);
//     // Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

//     // // Configure solver
//     // int verbose = 5;
//     // bool print_timings = true;
//     // int max_iter = 500;
//     // double tol = 1e-6;
//     // bool useHessianApproximation = false; 
//     // assert(!useHessianApproximation);
//     // std::string hessianApproximation = useHessianApproximation ? "exact" : "limited-memory";
//     // std::string linearSolver = "mumps";
//     // std::string muStrategy = "adaptive";
//     // std::string memoryUpdateType = "bfgs";

//     // app->Options()->SetIntegerValue("print_level", verbose);
//     // app->Options()->SetIntegerValue("print_frequency_iter", 10);
//     // app->Options()->SetIntegerValue("max_iter", max_iter);
//     // app->Options()->SetNumericValue("tol", tol);
//     // app->Options()->SetStringValue("linear_solver", linearSolver);
//     // app->Options()->SetStringValue("hessian_approximation", hessianApproximation);
//     // app->Options()->SetStringValue("mu_strategy", muStrategy);
//     // app->Options()->SetStringValue("limited_memory_update_type", memoryUpdateType);
//     // if(print_timings){
//     //     app->Options()->SetStringValue("timing_statistics", "yes");
//     //     app->Options()->SetStringValue("print_timing_statistics", "yes");
//     // }

//     // // Initialize and solve
//     // Ipopt::ApplicationReturnStatus status = app->Initialize();
//     // if (status != Ipopt::Solve_Succeeded) {
//     //     std::cout << "\n*** Error during initialization!\n";
//     //     return (int)status;
//     // }

//     // status = app->OptimizeTNLP(mynlp);

//     // if (status == Ipopt::Solve_Succeeded) {
//     //     std::cout << "\n*** The problem was solved successfully!\n";
//     // }
//     // else {
//     //     std::cout << "\n*** The problem failed to solve!\n";
//     // }

//     return 0;
// }
























// /*
// * Copyright (C) 2020-2024 MEmilio
// *
// * Authors: Ralf Hannemann-Tamas
// *
// * Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
// *
// * Licensed under the Apache License, Version 2.0 (the "License");
// * you may not use this file except in compliance with the License.
// * You may obtain a copy of the License at
// *
// *     http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS,
// * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// */


// #include "ad/ad.hpp"
// #include "ad/ad_spdlog_formatter.h" // IWYU pragma: keep

// #include "ode_seair/model.h"
// #include "ode_seair/infection_state.h"
// #include "ode_seair/parameters.h"
// #include "ode_secirvvs/model.h"
// #include "ode_secirvvs/infection_state.h"
// #include "ode_secirvvs/parameters.h"
// #include "memilio/compartments/simulation.h"
// #include "memilio/utils/logging.h"
// #include "memilio/utils/time_series.h"
// #include "IpTNLP.hpp"
// #include "IpIpoptApplication.hpp"
// #include <fstream>

// class Secirvvs_NLP : public Ipopt::TNLP
// {
// public:
//     Secirvvs_NLP()                               = default;
//     Secirvvs_NLP(const Secirvvs_NLP&)            = delete;
//     Secirvvs_NLP(Secirvvs_NLP&&)                 = delete;
//     Secirvvs_NLP& operator=(const Secirvvs_NLP&) = delete;
//     Secirvvs_NLP& operator=(Secirvvs_NLP&&)      = delete;
//     ~Secirvvs_NLP()                              = default;
//     bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
//                       IndexStyleEnum& index_style) override;
//     bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l,
//                          Ipopt::Number* g_u) override;

//     bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
//                             Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;

//     /** Method to return the objective value */
//     virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) override;

//     /** Method to return the gradient of the objective */
//     virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) override;

//     /** Method to return the constraint residuals */
//     virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) override;

//     /** Method to return:
//     *   1) The structure of the jacobian (if "values" is NULL)
//     *   2) The values of the jacobian (if "values" is not NULL)
//     */
//     virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
//                             Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) override;

//     /** Method to return:
//     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
//     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
//     */
//     virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
//                         const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
//                         Ipopt::Index* jCol, Ipopt::Number* values) override;

//     /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
//     virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
//                                    const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
//                                    const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
//                                    const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) override;
//     //@}

//     bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
//                                Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm,
//                                Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr,
//                                Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
//                                Ipopt::IpoptCalculatedQuantities* ip_cq) override;

//     template <typename FP = double>
//     void eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective);

// public:
//     int getN()
//     {
//         return n_;
//     }
//     int getM()
//     {
//         return m_;
//     }
//     static constexpr int tmax = 140;

// private:
//     const int numControlIntervals_ = 20;
//     const int numControls_         = 1;
//     const int numPathConstraints_  = 1;
//     const int pcresolution_ =
//         7; // the resolution of path constraints is by this factor higher than the control discretization
//     const int numIntervals_        = pcresolution_ * numControlIntervals_;
//     const int n_                   = numControlIntervals_ * numControls_;
//     const int m_                   = numIntervals_ * numPathConstraints_;
// };

// template <typename FP>
// void set_initial_values(mio::osecirvvs::Model<FP>& model)
// {

//     for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
//         model.populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}]                                = 10;
//         model.populations[{i, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}]                     = 11;
//         model.populations[{i, mio::osecirvvs::InfectionState::ExposedPartialImmunity}]                      = 12;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]                     = 13;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = 13;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]           = 14;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = 14;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = 15;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 15;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]                       = 5;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]              = 5;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]             = 6;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = 6;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]            = 7;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = 7;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereNaive}]                         = 8;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity}]              = 1;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity}]               = 2;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalNaive}]                       = 3;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}]             = 4;
//         model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}]            = 5;
//         model.populations[{i, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}]                 = 6;
//         model.populations[{i, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}]                  = 7;
//         model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadNaive}]                    = 0;
//         model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadPartialImmunity}]          = 0;
//         model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]         = 0;
//         model.populations.template set_difference_from_group_total<mio::AgeGroup>(
//             {i, mio::osecirvvs::InfectionState::SusceptibleNaive}, FP(1000));
//     }

//     model.parameters.template get<mio::osecirvvs::ICUCapacity<FP>>()          = 100;
//     model.parameters.template get<mio::osecirvvs::TestAndTraceCapacity<FP>>() = 0.0143;
//     const size_t daily_vaccinations                                           = 10;
//     model.parameters.template get<mio::osecirvvs::DailyPartialVaccinations<FP>>().resize(
//         mio::SimulationDay((size_t)ad::value(Secirvvs_NLP::tmax + 1)));
//     model.parameters.template get<mio::osecirvvs::DailyFullVaccinations<FP>>().resize(
//         mio::SimulationDay((size_t)ad::value(Secirvvs_NLP::tmax) + 1));
//     for (size_t i = 0; i < Secirvvs_NLP::tmax + 1; ++i) {
//         auto num_vaccinations = static_cast<FP>(i * daily_vaccinations);
//         model.parameters
//             .template get<mio::osecirvvs::DailyPartialVaccinations<FP>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
//             num_vaccinations;
//         model.parameters
//             .template get<mio::osecirvvs::DailyFullVaccinations<FP>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
//             num_vaccinations;
//     }
//     auto& contacts       = model.parameters.template get<mio::osecirvvs::ContactPatterns<FP>>();
//     auto& contact_matrix = contacts.get_cont_freq_mat();
//     contact_matrix[0].get_baseline().setConstant(0.5);
//     contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
//     //contact_matrix[0].add_damping(0.3, mio::SimulationTime<FP>(5.0));

//     //times
//     model.parameters.template get<mio::osecirvvs::TimeExposed<FP>>()[mio::AgeGroup(0)]            = 3.33;
//     model.parameters.template get<mio::osecirvvs::TimeInfectedNoSymptoms<FP>>()[mio::AgeGroup(0)] = 1.87;
//     model.parameters.template get<mio::osecirvvs::TimeInfectedSymptoms<FP>>()[mio::AgeGroup(0)]   = 7;
//     model.parameters.template get<mio::osecirvvs::TimeInfectedSevere<FP>>()[mio::AgeGroup(0)]     = 6;
//     model.parameters.template get<mio::osecirvvs::TimeInfectedCritical<FP>>()[mio::AgeGroup(0)]   = 7;

//     //probabilities
//     model.parameters.template get<mio::osecirvvs::TransmissionProbabilityOnContact<FP>>()[mio::AgeGroup(0)] = 0.15;
//     model.parameters.template get<mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>>()[mio::AgeGroup(0)]   = 0.5;
//     // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
//     // depends on incidence and test and trace capacity
//     model.parameters.template get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>>()[mio::AgeGroup(0)]    = 0.0;
//     model.parameters.template get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>>()[mio::AgeGroup(0)] = 0.4;
//     model.parameters.template get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>>()[mio::AgeGroup(0)]    = 0.2;
//     model.parameters.template get<mio::osecirvvs::SeverePerInfectedSymptoms<FP>>()[mio::AgeGroup(0)]         = 0.1;
//     model.parameters.template get<mio::osecirvvs::CriticalPerSevere<FP>>()[mio::AgeGroup(0)]                 = 0.1;
//     model.parameters.template get<mio::osecirvvs::DeathsPerCritical<FP>>()[mio::AgeGroup(0)]                 = 0.1;

//     model.parameters.template get<mio::osecirvvs::ReducExposedPartialImmunity<FP>>()[mio::AgeGroup(0)]          = 0.8;
//     model.parameters.template get<mio::osecirvvs::ReducExposedImprovedImmunity<FP>>()[mio::AgeGroup(0)]         = 0.331;
//     model.parameters.template get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>>()[mio::AgeGroup(0)] = 0.65;
//     model.parameters.template get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>>()[mio::AgeGroup(0)] =
//         0.243;
//     model.parameters
//         .template get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[mio::AgeGroup(0)] = 0.1;
//     model.parameters
//         .template get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[mio::AgeGroup(0)] = 0.091;
//     model.parameters.template get<mio::osecirvvs::ReducTimeInfectedMild<FP>>()[mio::AgeGroup(0)]               = 0.9;

//     model.parameters.template get<mio::osecirvvs::Seasonality<FP>>() = 0.2;
// }

// template <typename FP>
// void Secirvvs_NLP::eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective)
// {

//     FP t0   = 0;
//     FP tmax = this->tmax;
//     FP dt   = 0.2;
//     std::vector<FP> grid(numIntervals_ + 1);
//     for (int i = 0; i < numIntervals_ + 1; ++i) {
//         grid[i] = (tmax / numIntervals_) * i + (t0 / numIntervals_) * (numIntervals_ - i);
//     }
//     mio::osecirvvs::Model<FP> model(1);
//     auto& params = model.parameters;

//     set_initial_values(model);
//     int gridindex = 0;
//     objective     = FP(0.0);
//     for (int controlIndex = 0; controlIndex < numControlIntervals_; ++controlIndex) {
//         auto& contactControl = model.parameters.template get<mio::osecirvvs::ContactControl<FP>>();

//         // contactControl = x[controlIndex];

//         auto& contacts       = model.parameters.template get<mio::osecirvvs::ContactPatterns<FP>>();
//         auto& contact_matrix = contacts.get_cont_freq_mat();
//         //contact_matrix[0].get_baseline().setConstant(0.5);
//         //contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
//         if constexpr (std::is_same_v<FP, double>) {
//             contact_matrix[0].add_damping(FP(FP(1.) - x[controlIndex]),
//                                           mio::SimulationTime<FP>(controlIndex * pcresolution_));
//         }
//         else {
//             mio::SquareMatrixShape<FP> shape(1);
//             contact_matrix[0].add_damping(FP(FP(1.) - x[controlIndex]),
//                                           mio::SimulationTime<FP>(controlIndex * pcresolution_), std::move(shape));
//         }

//         objective -= x[controlIndex];

//         for (int i = 0; i < pcresolution_; ++i, ++gridindex) {

//             auto result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(grid[gridindex], grid[gridindex + 1], dt, model);

//             for (int j = 0; j < (int)mio::osecirvvs::InfectionState::Count; ++j) {
//                 model.populations[{mio::AgeGroup(0), mio::osecirvvs::InfectionState(j)}] = result.get_last_value()[j];
//             }

//             constraints[gridindex] =
//                 result.get_last_value()[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive];
//         }
//     }

//     return;
// }

// int main()
// {
//     std::cout<<"A"<<std::endl;

//     //switch of logging for mio
//     mio::set_log_level(mio::LogLevel::off);

//     // Create a new instance of your nlp
//     //  (use a SmartPtr, not raw)
//     Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new Secirvvs_NLP();

//     // Create a new instance of IpoptApplication
//     //  (use a SmartPtr, not raw)
//     // We are using the factory, since this allows us to compile this
//     // example with an Ipopt Windows DLL
//     Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

//     std::cout<<"B"<<std::endl;

//     // Change some options
//     // Note: The following choices are only examples, they might not be
//     //       suitable for your optimization problem.
//     // gut: 0,1,2,3,4
//     // schlecht: 5,6,7
//     app->Options()->SetIntegerValue("print_level", 4);
//     app->Options()->SetNumericValue("tol", 1e-6);
//     app->Options()->SetStringValue("mu_strategy", "adaptive");
//     app->Options()->SetStringValue("output_file", "ipopt.out");
//     app->Options()->SetStringValue("hessian_approximation", "limited-memory");
//     app->Options()->SetStringValue("limited_memory_update_type", "bfgs");

//     std::cout<<"C"<<std::endl;

//     // Initialize the IpoptApplication and process the options
//     Ipopt::ApplicationReturnStatus status;
//     status = app->Initialize();
//     if (status != Ipopt::Solve_Succeeded) {
//         std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
//         return (int)status;
//     }

//     std::cout<<"D"<<std::endl;

//     // Ask Ipopt to solve the problem
//     status = app->OptimizeTNLP(mynlp);

//     std::cout<<"E"<<std::endl;

//     if (status == Ipopt::Solve_Succeeded) {
//         std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
//     }
//     else {
//         std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
//     }

//     // As the SmartPtrs go out of scope, the reference count
//     // will be decremented and the objects will automatically
//     // be deleted.

//     std::cout<<"F"<<std::endl;

//     return (int)status;
// }

// bool Secirvvs_NLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
//                                 IndexStyleEnum& index_style)
// {
//     n = n_;

//     m = m_;

//     // in this example the jacobian is dense
//     nnz_jac_g = m_ * n_;

//     // the Hessian is also dense
//     nnz_h_lag = n_ * n_;

//     // use the C style indexing (0-based)
//     index_style = Ipopt::TNLP::C_STYLE;

//     return true;
// }

// bool Secirvvs_NLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m,
//                                    Ipopt::Number* g_l, Ipopt::Number* g_u)
// {
//     // controls order: 1. alpha_a, 2. alpha_i, 3. kappa
//     for (int i = 0; i < numControlIntervals_; ++i) {
//         x_l[i] = 0.2; // lower bound of ContactControl
//         x_u[i] = 1.0; // upper bound of ContactControl
//     }

//     // path constraints
//     for (int i = 0; i < m_; ++i) {
//         g_l[i] = 0.0;
//         g_u[i] = 15.;
//     }
//     return true;
// }

// bool Secirvvs_NLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
//                                       Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda)
// {
//     assert(init_z == false);
//     assert(init_lambda == false);

//     for (int i = 0; i < n; ++i) {
//         x[i] = 1.0;
//     }
//     return true;
// }

// bool Secirvvs_NLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
// {
//     std::vector<double> xx(getN());
//     std::vector<double> constraints(getM());
//     for (int i = 0; i < n; ++i)
//         xx[i] = x[i];
//     eval_objective_constraints(xx, constraints, obj_value);
//     return true;
// }

// bool Secirvvs_NLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
// {
//     using FP = ad::gt1s<double>::type;
//     std::vector<FP> xx(getN());
//     std::vector<FP> constraints(getM());
//     FP objective;
//     for (int i = 0; i < n; ++i)
//         ad::value(xx[i]) = x[i];
//     for (int i = 0; i < n; ++i) {
//         ad::derivative(xx[i]) = 1.0;
//         eval_objective_constraints(xx, constraints, objective);
//         grad_f[i]             = ad::derivative(objective);
//         ad::derivative(xx[i]) = 0.0;
//     }
//     return true;
// }

// bool Secirvvs_NLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
// {
//     std::vector<double> xx(getN());
//     std::vector<double> constraints(getM());
//     double obj_value = 0;
//     for (int i = 0; i < n; ++i)
//         xx[i] = x[i];
//     eval_objective_constraints(xx, constraints, obj_value);
//     for (int i = 0; i < m; ++i)
//         g[i] = constraints[i];
//     return true;
// }

// bool Secirvvs_NLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
//                               Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
// {

//     if (values == nullptr) {
//         int jac_index = 0;
//         for (int i = 0; i < n; ++i) {
//             for (int j = 0; j < m; ++j) {
//                 iRow[jac_index] = j;
//                 jCol[jac_index] = i;
//                 ++jac_index;
//             }
//         }
//     }
//     else {
//         using FP = ad::gt1s<double>::type;
//         std::vector<FP> xx(getN());
//         std::vector<FP> constraints(getM());
//         FP objective;
//         int jac_index = 0;
//         for (int i = 0; i < n; ++i)
//             ad::value(xx[i]) = x[i];
//         for (int i = 0; i < n; ++i) {
//             ad::derivative(xx[i]) = 1.0;
//             eval_objective_constraints(xx, constraints, objective);
//             for (int j = 0; j < m; ++j) {
//                 values[jac_index] = ad::derivative(constraints[j]);
//                 ++jac_index;
//             }
//             ad::derivative(xx[i]) = 0.0;
//         }
//     }
//     return true;
// }

// bool Secirvvs_NLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
//                           const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
//                           Ipopt::Index* jCol, Ipopt::Number* values)
// {
//     return true;
// }

// void Secirvvs_NLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
//                                      const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
//                                      const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
//                                      const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
// {
//     std::cout << "optimal solution is\n";
//     for (Ipopt::Index i = 0; i < n; ++i) {
//         std::cout << x[i] << std::endl;
//     }
//     std::cout << "Constraints are" << std::endl;
//     for (Ipopt::Index i = 0; i < m; ++i) {
//         std::cout << g[i] << std::endl;
//     }
//     std::cout << "Writing output to text files" << std::endl;

//     double t0   = 0;
//     double tmax = this->tmax;
//     double dt   = 0.2;

//     mio::osecirvvs::Model<double> model(1);
//     set_initial_values(model);

//     std::vector<double> grid(numIntervals_ + 1);
//     for (int i = 0; i < numIntervals_ + 1; ++i) {
//         grid[i] = (tmax / numIntervals_) * i + (t0 / numIntervals_) * (numIntervals_ - i);
//     }
//     std::ofstream outFileContactControl("contactControl.txt");

//     int gridindex = 0;
//     for (int controlIndex = 0; controlIndex < numControlIntervals_; ++controlIndex) {
//         model.parameters.template get<mio::osecirvvs::ContactControl<double>>() = x[controlIndex];
//         outFileContactControl << grid[gridindex] << " "
//                               << model.parameters.template get<mio::osecirvvs::ContactControl<double>>() << "\n";
//         for (int i = 0; i < pcresolution_; ++i, ++gridindex) {
//             auto result =
//                 mio::simulate<double, mio::osecirvvs::Model<double>>(grid[gridindex], grid[gridindex + 1], dt, model);

//             for (int j = 0; j < (int)mio::oseair::InfectionState::Count; ++j) {
//                 model.populations[{mio::AgeGroup(0), mio::oseair::InfectionState(j)}] = result.get_last_value()[j];
//             }
//         }
//         outFileContactControl << grid[gridindex] << " "
//                               << model.parameters.template get<mio::osecirvvs::ContactControl<double>>() << "\n";
//     }

//     outFileContactControl.close();

//     return;
// }

// bool Secirvvs_NLP::intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
//                                          Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
//                                          Ipopt::Number d_norm, Ipopt::Number regularization_size,
//                                          Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
//                                          const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
// {
//     return true;
// }
























































/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Ralf Hannemann-Tamas
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


#include "ad/ad.hpp"

#include "ode_seair/model.h"
#include "ode_seair/infection_state.h"
#include "ode_seair/parameters.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include <fstream>

class Secirvvs_NLP : public Ipopt::TNLP
{
public:
    Secirvvs_NLP()                               = default;
    Secirvvs_NLP(const Secirvvs_NLP&)            = delete;
    Secirvvs_NLP(Secirvvs_NLP&&)                 = delete;
    Secirvvs_NLP& operator=(const Secirvvs_NLP&) = delete;
    Secirvvs_NLP& operator=(Secirvvs_NLP&&)      = delete;
    ~Secirvvs_NLP()                              = default;
    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      IndexStyleEnum& index_style) override;
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l,
                         Ipopt::Number* g_u) override;

    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
                            Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;

    /** Method to return the objective value */
    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) override;

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) override;

    /** Method to return the constraint residuals */
    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) override;

    /** Method to return:
    *   1) The structure of the jacobian (if "values" is NULL)
    *   2) The values of the jacobian (if "values" is not NULL)
    */
    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                            Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) override;

    /** Method to return:
    *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
    *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
    */
    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
                        const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                        Ipopt::Index* jCol, Ipopt::Number* values) override;

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                   const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                   const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                   const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) override;
    //@}

    bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
                               Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm,
                               Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr,
                               Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
                               Ipopt::IpoptCalculatedQuantities* ip_cq) override;

    template <typename FP = double>
    void eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective);

public:
    int getN()
    {
        return n_;
    }
    int getM()
    {
        return m_;
    }
    static constexpr int tmax = 140;

private:
    const int numControlIntervals_ = 20;
    const int numControls_         = 1;
    const int numPathConstraints_  = 1;
    const int pcresolution_ =
        7; // the resolution of path constraints is by this factor higher than the control discretization
    const int numIntervals_        = pcresolution_ * numControlIntervals_;
    const int n_                   = numControlIntervals_ * numControls_;
    const int m_                   = numIntervals_ * numPathConstraints_;
};

template <typename FP>
void set_initial_values(mio::osecirvvs::Model<FP>& model)
{

    for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}]                                = 10;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}]                     = 11;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedPartialImmunity}]                      = 12;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]                     = 13;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = 13;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]           = 14;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = 14;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = 15;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 15;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]                       = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]              = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]             = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]            = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereNaive}]                         = 8;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity}]              = 1;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity}]               = 2;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalNaive}]                       = 3;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}]             = 4;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}]            = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}]                 = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}]                  = 7;
        model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadNaive}]                    = 0;
        model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadPartialImmunity}]          = 0;
        model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]         = 0;
        model.populations.template set_difference_from_group_total<mio::AgeGroup>(
            {i, mio::osecirvvs::InfectionState::SusceptibleNaive}, FP(1000));
    }

    model.parameters.template get<mio::osecirvvs::ICUCapacity<FP>>()          = 100;
    model.parameters.template get<mio::osecirvvs::TestAndTraceCapacity<FP>>() = 0.0143;
    const size_t daily_vaccinations                                           = 10;
    model.parameters.template get<mio::osecirvvs::DailyPartialVaccinations<FP>>().resize(
        mio::SimulationDay((size_t)ad::value(Secirvvs_NLP::tmax + 1)));
    model.parameters.template get<mio::osecirvvs::DailyFullVaccinations<FP>>().resize(
        mio::SimulationDay((size_t)ad::value(Secirvvs_NLP::tmax) + 1));
    for (size_t i = 0; i < Secirvvs_NLP::tmax + 1; ++i) {
        auto num_vaccinations = static_cast<FP>(i * daily_vaccinations);
        model.parameters
            .template get<mio::osecirvvs::DailyPartialVaccinations<FP>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters
            .template get<mio::osecirvvs::DailyFullVaccinations<FP>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
    }
    auto& contacts       = model.parameters.template get<mio::osecirvvs::ContactPatterns<FP>>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    //contact_matrix[0].add_damping(0.3, mio::SimulationTime<FP>(5.0));

    //times
    model.parameters.template get<mio::osecirvvs::TimeExposed<FP>>()[mio::AgeGroup(0)]            = 3.33;
    model.parameters.template get<mio::osecirvvs::TimeInfectedNoSymptoms<FP>>()[mio::AgeGroup(0)] = 1.87;
    model.parameters.template get<mio::osecirvvs::TimeInfectedSymptoms<FP>>()[mio::AgeGroup(0)]   = 7;
    model.parameters.template get<mio::osecirvvs::TimeInfectedSevere<FP>>()[mio::AgeGroup(0)]     = 6;
    model.parameters.template get<mio::osecirvvs::TimeInfectedCritical<FP>>()[mio::AgeGroup(0)]   = 7;

    //probabilities
    model.parameters.template get<mio::osecirvvs::TransmissionProbabilityOnContact<FP>>()[mio::AgeGroup(0)] = 0.15;
    model.parameters.template get<mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>>()[mio::AgeGroup(0)]   = 0.5;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    model.parameters.template get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>>()[mio::AgeGroup(0)]    = 0.0;
    model.parameters.template get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>>()[mio::AgeGroup(0)] = 0.4;
    model.parameters.template get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>>()[mio::AgeGroup(0)]    = 0.2;
    model.parameters.template get<mio::osecirvvs::SeverePerInfectedSymptoms<FP>>()[mio::AgeGroup(0)]         = 0.1;
    model.parameters.template get<mio::osecirvvs::CriticalPerSevere<FP>>()[mio::AgeGroup(0)]                 = 0.1;
    model.parameters.template get<mio::osecirvvs::DeathsPerCritical<FP>>()[mio::AgeGroup(0)]                 = 0.1;

    model.parameters.template get<mio::osecirvvs::ReducExposedPartialImmunity<FP>>()[mio::AgeGroup(0)]          = 0.8;
    model.parameters.template get<mio::osecirvvs::ReducExposedImprovedImmunity<FP>>()[mio::AgeGroup(0)]         = 0.331;
    model.parameters.template get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>>()[mio::AgeGroup(0)] = 0.65;
    model.parameters.template get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>>()[mio::AgeGroup(0)] =
        0.243;
    model.parameters
        .template get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[mio::AgeGroup(0)] = 0.1;
    model.parameters
        .template get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[mio::AgeGroup(0)] = 0.091;
    model.parameters.template get<mio::osecirvvs::ReducTimeInfectedMild<FP>>()[mio::AgeGroup(0)]               = 0.9;

    model.parameters.template get<mio::osecirvvs::Seasonality<FP>>() = 0.2;
}

template <typename FP>
void Secirvvs_NLP::eval_objective_constraints(const std::vector<FP>& x, std::vector<FP>& constraints, FP& objective)
{

    FP t0   = 0;
    FP tmax = this->tmax;
    FP dt   = 0.2;
    std::vector<FP> grid(numIntervals_ + 1);
    for (int i = 0; i < numIntervals_ + 1; ++i) {
        grid[i] = (tmax / numIntervals_) * i + (t0 / numIntervals_) * (numIntervals_ - i);
    }
    mio::osecirvvs::Model<FP> model(1);
    auto& params = model.parameters;

    set_initial_values(model);
    int gridindex = 0;
    objective     = FP(0.0);
    for (int controlIndex = 0; controlIndex < numControlIntervals_; ++controlIndex) {
        // auto& contactControl = model.parameters.template get<mio::osecirvvs::ContactControl<FP>>();

        auto contactControl = x[controlIndex];

        auto& contacts       = model.parameters.template get<mio::osecirvvs::ContactPatterns<FP>>();
        auto& contact_matrix = contacts.get_cont_freq_mat();
        //contact_matrix[0].get_baseline().setConstant(0.5);
        //contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
        if constexpr (std::is_same_v<FP, double>) {
            contact_matrix[0].add_damping(FP(FP(1.) - x[controlIndex]),
                                          mio::SimulationTime<FP>(controlIndex * pcresolution_));
        }
        else {
            mio::SquareMatrixShape<FP> shape(1);
            contact_matrix[0].add_damping(FP(FP(1.) - x[controlIndex]),
                                          mio::SimulationTime<FP>(controlIndex * pcresolution_), std::move(shape));
        }

        objective -= x[controlIndex];

        for (int i = 0; i < pcresolution_; ++i, ++gridindex) {

            auto result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(grid[gridindex], grid[gridindex + 1], dt, model);

            for (int j = 0; j < (int)mio::osecirvvs::InfectionState::Count; ++j) {
                model.populations[{mio::AgeGroup(0), mio::osecirvvs::InfectionState(j)}] = result.get_last_value()[j];
            }

            constraints[gridindex] =
                result.get_last_value()[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive];
        }
    }

    return;
}

int main()
{
    //switch of logging for mio
    mio::set_log_level(mio::LogLevel::off);

    // Create a new instance of your nlp
    //  (use a SmartPtr, not raw)
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new Secirvvs_NLP();

    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // Change some options
    // Note: The following choices are only examples, they might not be
    //       suitable for your optimization problem.
    app->Options()->SetNumericValue("tol", 1e-6);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("limited_memory_update_type", "bfgs");

    // Initialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int)status;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if (status == Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else {
        std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    }

    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.

    return (int)status;
}

bool Secirvvs_NLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                                IndexStyleEnum& index_style)
{
    n = n_;

    m = m_;

    // in this example the jacobian is dense
    nnz_jac_g = m_ * n_;

    // the Hessian is also dense
    nnz_h_lag = n_ * n_;

    // use the C style indexing (0-based)
    index_style = Ipopt::TNLP::C_STYLE;

    return true;
}

bool Secirvvs_NLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m,
                                   Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    // controls order: 1. alpha_a, 2. alpha_i, 3. kappa
    for (int i = 0; i < numControlIntervals_; ++i) {
        x_l[i] = 0.1; // lower bound of ContactControl
        x_u[i] = 1.0; // upper bound of ContactControl
    }

    // path constraints
    for (int i = 0; i < m_; ++i) {
        g_l[i] = 0.0;
        g_u[i] = 150.;
    }
    return true;
}

bool Secirvvs_NLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
                                      Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda)
{
    assert(init_z == false);
    assert(init_lambda == false);

    for (int i = 0; i < n; ++i) {
        x[i] = 1.0;
    }
    return true;
}

bool Secirvvs_NLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
    std::vector<double> xx(getN());
    std::vector<double> constraints(getM());
    for (int i = 0; i < n; ++i)
        xx[i] = x[i];
    eval_objective_constraints(xx, constraints, obj_value);
    return true;
}

bool Secirvvs_NLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    using FP = ad::gt1s<double>::type;
    std::vector<FP> xx(getN());
    std::vector<FP> constraints(getM());
    FP objective;
    for (int i = 0; i < n; ++i)
        ad::value(xx[i]) = x[i];
    for (int i = 0; i < n; ++i) {
        ad::derivative(xx[i]) = 1.0;
        eval_objective_constraints(xx, constraints, objective);
        grad_f[i]             = ad::derivative(objective);
        ad::derivative(xx[i]) = 0.0;
    }
    return true;
}

bool Secirvvs_NLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
    std::vector<double> xx(getN());
    std::vector<double> constraints(getM());
    double obj_value = 0;
    for (int i = 0; i < n; ++i)
        xx[i] = x[i];
    eval_objective_constraints(xx, constraints, obj_value);
    for (int i = 0; i < m; ++i)
        g[i] = constraints[i];
    return true;
}

bool Secirvvs_NLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                              Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{

    if (values == nullptr) {
        int jac_index = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                iRow[jac_index] = j;
                jCol[jac_index] = i;
                ++jac_index;
            }
        }
    }
    else {
        using FP = ad::gt1s<double>::type;
        std::vector<FP> xx(getN());
        std::vector<FP> constraints(getM());
        FP objective;
        int jac_index = 0;
        for (int i = 0; i < n; ++i)
            ad::value(xx[i]) = x[i];
        for (int i = 0; i < n; ++i) {
            ad::derivative(xx[i]) = 1.0;
            eval_objective_constraints(xx, constraints, objective);
            for (int j = 0; j < m; ++j) {
                values[jac_index] = ad::derivative(constraints[j]);
                ++jac_index;
            }
            ad::derivative(xx[i]) = 0.0;
        }
    }
    return true;
}

bool Secirvvs_NLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
                          const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                          Ipopt::Index* jCol, Ipopt::Number* values)
{
    return true;
}

void Secirvvs_NLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
                                     const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
                                     const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                     const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    std::cout << "optimal solution is\n";
    for (Ipopt::Index i = 0; i < n; ++i) {
        std::cout << x[i] << std::endl;
    }
    std::cout << "Constraints are" << std::endl;
    for (Ipopt::Index i = 0; i < m; ++i) {
        std::cout << g[i] << std::endl;
    }
    std::cout << "Writing output to text files" << std::endl;

    double t0   = 0;
    double tmax = this->tmax;
    double dt   = 0.2;

    mio::osecirvvs::Model<double> model(1);
    set_initial_values(model);

    std::vector<double> grid(numIntervals_ + 1);
    for (int i = 0; i < numIntervals_ + 1; ++i) {
        grid[i] = (tmax / numIntervals_) * i + (t0 / numIntervals_) * (numIntervals_ - i);
    }
    // std::ofstream outFileContactControl("contactControl.txt");

    // int gridindex = 0;
    // for (int controlIndex = 0; controlIndex < numControlIntervals_; ++controlIndex) {
    //     model.parameters.template get<mio::osecirvvs::ContactControl<double>>() = x[controlIndex];
    //     outFileContactControl << grid[gridindex] << " "
    //                           << model.parameters.template get<mio::osecirvvs::ContactControl<double>>() << "\n";
    //     for (int i = 0; i < pcresolution_; ++i, ++gridindex) {
    //         auto result =
    //             mio::simulate<double, mio::osecirvvs::Model<double>>(grid[gridindex], grid[gridindex + 1], dt, model);

    //         for (int j = 0; j < (int)mio::oseair::InfectionState::Count; ++j) {
    //             model.populations[{mio::AgeGroup(0), mio::oseair::InfectionState(j)}] = result.get_last_value()[j];
    //         }
    //     }
    //     outFileContactControl << grid[gridindex] << " "
    //                           << model.parameters.template get<mio::osecirvvs::ContactControl<double>>() << "\n";
    // }

    // outFileContactControl.close();

    return;
}

bool Secirvvs_NLP::intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
                                         Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
                                         Ipopt::Number d_norm, Ipopt::Number regularization_size,
                                         Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
                                         const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    return true;
}