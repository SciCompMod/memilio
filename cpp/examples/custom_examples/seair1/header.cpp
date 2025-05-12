#include "header.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <array>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "settings.h"
#include "IpIpoptData.hpp"

#include "models/ode_seair/model.h"
#include "memilio/compartments/simulation.h"

template <typename FP>
FP Seair_NLP::objective_function(const FP* parameters, std::size_t n) {
    FP objective = 0.0;
    for (std::size_t controlIndex = 0; controlIndex < settings_.numControlIntervals; controlIndex++){
        FP socialDistancing = parameters[controlIndex * settings_.numControls + 0];
        FP quarantined      = parameters[controlIndex * settings_.numControls + 1];
        FP testingRate      = parameters[controlIndex * settings_.numControls + 2];

        objective += settings_.pcResolution * (-socialDistancing - quarantined + 0.1 * testingRate);
    }
    return objective;
}

template <typename FP>
void Seair_NLP::constraint_functions(const FP* parameters, std::size_t n, FP* constraints, std::size_t m) {

    std::vector<FP> grid(settings_.numIntervals + 1);
    FP grid_spacing = (settings_.tmax - settings_.t0) / settings_.numIntervals;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        grid[i] = settings_.t0 + i * grid_spacing;
    }

    mio::oseair::Model<FP> model;
    set_initial_values(model, settings_);

    for (std::size_t controlIndex = 0; controlIndex < settings_.numControlIntervals; controlIndex++){

        FP socialDistancing = parameters[controlIndex * settings_.numControls + 0];
        FP quarantined      = parameters[controlIndex * settings_.numControls + 1];
        FP testingRate      = parameters[controlIndex * settings_.numControls + 2];
        model.parameters.template get<mio::oseair::SocialDistancing<FP>>() = socialDistancing;
        model.parameters.template get<mio::oseair::Quarantined<FP>>()      = quarantined;
        model.parameters.template get<mio::oseair::TestingRate<FP>>()      = testingRate;

        for (std::size_t i = 0; i < settings_.pcResolution ; i++)
        {
            std::size_t gridindex = controlIndex * settings_.pcResolution + i;
            auto result = mio::simulate<FP, mio::oseair::Model<FP>>(grid[gridindex], grid[gridindex + 1], settings_.dt, model);

            for (int j = 0; j < (int)mio::oseair::InfectionState::Count; ++j) {
                model.populations[mio::oseair::InfectionState(j)] = result.get_last_value()[j];
            }
            constraints[gridindex] = result.get_last_value()[(int)mio::oseair::InfectionState::Infected];
        }
    }
}

template <typename FP>
void set_initial_values(mio::oseair::Model<FP>& model, const ProblemSettings& settings)
{
    using IS = mio::oseair::InfectionState;

    const std::array<std::pair<IS, FP>, 6> init = {{
        {IS::Susceptible,  0.9977558755803503 * settings.N},
        {IS::Exposed,      0.0003451395725394549 * settings.N},
        {IS::Asymptomatic, 0.00037846880968213874 * settings.N},
        {IS::Infected,     337072.0},
        {IS::Recovered,    17448.0},
        {IS::Dead,         9619.0}
    }};

    for (const auto& [state, value] : init) {
        model.populations[{mio::Index<IS>(state)}] = value;
    }
}

constexpr Ipopt::Number INF = 1e19;

Seair_NLP::Seair_NLP(const ProblemSettings& settings) : 
    settings_(settings)
{
    n_ = settings_.numControlIntervals * settings_.numControls;
    m_ = settings_.numIntervals * settings_.numPathConstraints;

    nnz_jac_g_ = n_ * m_;
    bool use_hessian_approximation = false;
    nnz_h_lag_ = use_hessian_approximation ? n_ * (n_ + 1) / 2 : 0;

    // Resize vectors
    x_l_.resize(n_);
    x_u_.resize(n_);
    g_l_.resize(m_);
    g_u_.resize(m_);

    // Fill x_l_ and x_u_ from control bounds
    for (int controlInterval = 0; controlInterval < settings_.numControlIntervals; ++controlInterval) {
        for (int controlIndex = 0; controlIndex < settings_.numControls; ++controlIndex) {
            int idx = controlInterval * settings_.numControls + controlIndex;

            const auto& controlTuple = settings_.controlBounds[controlIndex];
            const auto& bounds = std::get<1>(controlTuple);
            x_l_[idx] = bounds.first;
            x_u_[idx] = bounds.second;
        }
    }

    // Fill g_l_ and g_u_ from path constraint bounds
    for (int interval = 0; interval < settings_.numIntervals; interval++){
        for (int constraintIndex = 0; constraintIndex < settings_.numPathConstraints; ++constraintIndex) {
            int idx = interval * settings_.numPathConstraints + constraintIndex;

            const auto& bounds = settings_.pathConstraints[constraintIndex].second;
            g_l_[idx] = bounds.first;
            g_u_[idx] = bounds.second;
        }
    }
    
    if (!ad::ga1s<double>::global_tape) {
        ad::ga1s<double>::global_tape = tape_t::create();
    }
    tape_ = ad::ga1s<double>::global_tape;
}

Seair_NLP::~Seair_NLP() {
    if (tape_) {
        tape_t::remove(tape_);
        tape_ = nullptr;
    }
}

bool Seair_NLP::get_nlp_info(
    Ipopt::Index& n,
    Ipopt::Index& m,
    Ipopt::Index& nnz_jac_g,
    Ipopt::Index& nnz_h_lag,
    Ipopt::TNLP::IndexStyleEnum& index_style
) {
    n = n_;
    m = m_;
    nnz_jac_g = nnz_jac_g_;
    nnz_h_lag = nnz_h_lag_;
    index_style = Ipopt::TNLP::C_STYLE;
    return true;
}

bool Seair_NLP::get_bounds_info(
    Ipopt::Index n,
    Ipopt::Number* x_l,
    Ipopt::Number* x_u,
    Ipopt::Index m,
    Ipopt::Number* g_l,
    Ipopt::Number* g_u
) {
    assert(n == n_ && m == m_);
    for (Ipopt::Index i = 0; i < n_; ++i) {
        x_l[i] = x_l_[i];
        x_u[i] = x_u_[i];
    }

    for (Ipopt::Index i = 0; i < m_; ++i) {
        g_l[i] = g_l_[i];
        g_u[i] = g_u_[i];
    }

    return true;
}

bool Seair_NLP::get_starting_point(
    Ipopt::Index n,
    bool init_x,
    Ipopt::Number* x,
    bool init_z,
    Ipopt::Number* /*z_L*/,
    Ipopt::Number* /*z_U*/,
    Ipopt::Index m,
    bool init_lambda,
    Ipopt::Number* /*lambda*/
) {
    assert(init_x && !init_z && !init_lambda);
    assert(n == n_ && m == m_);

    for (int controlInterval = 0; controlInterval < settings_.numControlIntervals; ++controlInterval) {
        for (int controlIndex = 0; controlIndex < settings_.numControls; ++controlIndex) {
            int idx = controlInterval * settings_.numControls + controlIndex;

            const auto& controlTuple = settings_.controlBounds[controlIndex];
            double initial = std::get<2>(controlTuple);
            x[idx] = initial;
        }
    }

    return true;
}

bool Seair_NLP::eval_f(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool /*new_x*/,
    Ipopt::Number& obj_value
) {
    assert(n == n_);

    obj_value = objective_function(x, static_cast<std::size_t>(n));

    return true;
}

bool Seair_NLP::eval_g(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool /*new_x*/,
    Ipopt::Index m,
    Ipopt::Number* g
) {
    assert(n == n_ && m == m_);

    constraint_functions(x, static_cast<std::size_t>(n), g, static_cast<std::size_t>(m));

    return true;
}

bool Seair_NLP::eval_grad_f(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool new_x,
    Ipopt::Number* grad_f
) {
    bool use_forward_mode = false;

    if (use_forward_mode) {
        // Use forward mode for Jacobian evaluation
        return eval_grad_f_forward(n, x, new_x, grad_f);
    } else {
        // Use reverse mode for Jacobian evaluation
        return eval_grad_f_reverse(n, x, new_x, grad_f);
    }
}

bool Seair_NLP::eval_jac_g(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool new_x,
    Ipopt::Index m,
    Ipopt::Index nele_jac,
    Ipopt::Index* iRow,
    Ipopt::Index* jCol,
    Ipopt::Number* values
) {
    assert(n == n_ && m == m_ && nele_jac == nnz_jac_g_);

    bool use_forward_mode = true;

    if (use_forward_mode) {
        // Use forward mode for Jacobian evaluation
        return eval_jac_g_forward(n, x, new_x, m, nele_jac, iRow, jCol, values);
    } else {
        // Use reverse mode for Jacobian evaluation
        return eval_jac_g_reverse(n, x, new_x, m, nele_jac, iRow, jCol, values);
    }
}

bool Seair_NLP::eval_h(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool /*new_x*/,
    Ipopt::Number obj_factor,
    Ipopt::Index m,
    const Ipopt::Number* lambda,
    bool /*new_lambda*/,
    Ipopt::Index nele_hess,
    Ipopt::Index* iRow,
    Ipopt::Index* jCol,
    Ipopt::Number* values
) {
    assert(!use_hessian_approximation_);
    return false;
}

bool Seair_NLP::eval_grad_f_forward(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool /*new_x*/,
    Ipopt::Number* grad_f
) {
    assert(n == n_);

    using ad_t = ad::gt1s<double>::type;

    #pragma omp parallel num_threads(1)
    {
        std::vector<ad_t> x_ad(n);
        for (Ipopt::Index i = 0; i < n; ++i) {
            x_ad[i] = x[i];
            ad::derivative(x_ad[i]) = 0.0;
        } 
        
        #pragma omp for
        for (Ipopt::Index column = 0; column < n; ++column) {
            ad::derivative(x_ad[column]) = 1.0;
            ad_t obj_ad = objective_function(x_ad.data(), static_cast<std::size_t>(n));
            grad_f[column] = ad::derivative(obj_ad);
            ad::derivative(x_ad[column]) = 0.0;
        }
    }

    return true;
}

bool Seair_NLP::eval_grad_f_reverse(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool /*new_x*/,
    Ipopt::Number* grad_f
) {
    assert(n == n_);
    using ad_t = ad::ga1s<double>::type;

    tape_->reset();

    std::vector<ad_t> x_ad(n);

    for (Ipopt::Index i = 0; i < n; ++i) {
        ad::value(x_ad[i]) = x[i];
        ad::derivative(x_ad[i]) = 0.0;
        tape_->register_variable(x_ad[i]);
    }

    ad_t obj_ad = objective_function(x_ad.data(), static_cast<std::size_t>(n));
    tape_->register_output_variable(obj_ad);
    ad::derivative(obj_ad) = 1.0;

    tape_->interpret_adjoint();

    for (Ipopt::Index i = 0; i < n; ++i) {
        grad_f[i] = ad::derivative(x_ad[i]);
    }

    return true;
}

bool Seair_NLP::eval_jac_g_forward(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool /*new_x*/,
    Ipopt::Index m,
    Ipopt::Index nele_jac,
    Ipopt::Index* iRow,
    Ipopt::Index* jCol,
    Ipopt::Number* values
) {
    assert(n == n_ && m == m_ && nele_jac == nnz_jac_g_);

    if (values == nullptr) {
        Ipopt::Index idx = 0;
        for (Ipopt::Index column = 0; column < n; ++column) {
            for (Ipopt::Index row = 0; row < m; ++row) {
                iRow[idx] = row;
                jCol[idx] = column;
                ++idx;
            }
        }
    }
    else {
        using ad_t = ad::gt1s<double>::type;

        #pragma omp parallel num_threads(1)
        {
            std::vector<ad_t> x_ad(n);
            std::vector<ad_t> g_ad(m);

            for (Ipopt::Index i = 0; i < n; ++i) {
                x_ad[i] = x[i];
                ad::derivative(x_ad[i]) = 0.0;
            }

            #pragma omp for
            for (Ipopt::Index column = 0; column < n; ++column) {
                ad::derivative(x_ad[column]) = 1.0;
                constraint_functions(x_ad.data(), static_cast<std::size_t>(n),
                                     g_ad.data(), static_cast<std::size_t>(m));

                for (Ipopt::Index row = 0; row < m; ++row) {
                    values[column * m + row] = ad::derivative(g_ad[row]);
                }

                ad::derivative(x_ad[column]) = 0.0;
            }
        }
    }

    return true;
}

bool Seair_NLP::eval_jac_g_reverse(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool /*new_x*/,
    Ipopt::Index m,
    Ipopt::Index nele_jac,
    Ipopt::Index* iRow,
    Ipopt::Index* jCol,
    Ipopt::Number* values
) {
    assert(n == n_ && m == m_ && nele_jac == nnz_jac_g_);

    if (values == nullptr) {
        // Sparsity pattern setup
        Ipopt::Index idx = 0;
        for (Ipopt::Index row = 0; row < m; ++row) {
            for (Ipopt::Index column = 0; column < n; ++column) {
                iRow[idx] = row;
                jCol[idx] = column;
                ++idx;
            }
        }
    } else {
        using ad_t = ad::ga1s<double>::type;

        tape_->reset();
        
        std::vector<ad_t> x_ad(n);
        for (Ipopt::Index j = 0; j < n; ++j) {
            ad::value(x_ad[j]) = x[j];
            ad::derivative(x_ad[j]) = 0.0;
            tape_->register_variable(x_ad[j]);
        }

        std::vector<ad_t> g_ad(m);
        constraint_functions(x_ad.data(), n, g_ad.data(), m);

        for (Ipopt::Index i = 0; i < m; ++i) {
            tape_->register_output_variable(g_ad[i]);
        }

        for (Ipopt::Index i = 0; i < m; ++i) {
            tape_->zero_adjoints();
            ad::derivative(g_ad[i]) = 1.0;
            tape_->interpret_adjoint(); // Takes up 97% of the time
            for (Ipopt::Index j = 0; j < n; ++j) {
                values[i * n + j] = ad::derivative(x_ad[j]);
            }
        }
    }

    return true;
}

void Seair_NLP::finalize_solution(
    Ipopt::SolverReturn status,
    Ipopt::Index n,
    const Ipopt::Number* x,
    const Ipopt::Number* z_L,
    const Ipopt::Number* z_U,
    Ipopt::Index m,
    const Ipopt::Number* g,
    const Ipopt::Number* lambda,
    Ipopt::Number obj_value,
    const Ipopt::IpoptData* ip_data,
    Ipopt::IpoptCalculatedQuantities* /*ip_cq*/
) {
    // Print the final objective value
    std::cout << "Final Objective Value: " << obj_value << "\n";

    // // Print the final solution for x
    // std::cout << "Optimal Solution (x): \n";
    // for (Ipopt::Index i = 0; i < n; ++i) {
    //     std::cout << "x[" << i << "] = " << x[i] << "\n";
    // }

    // // Optionally, print the constraints values and multipliers
    // std::cout << "Constraints (g): \n";
    // for (Ipopt::Index i = 0; i < m; ++i) {
    //     std::cout << "g[" << i << "] = " << g[i] << "\n";
    // }

    if (ip_data) {
        std::cout << "Number of iterations: " << ip_data->iter_count() << "\n";
    }

    // Load the Data in a .csv file

    std::ofstream outputFile("output.csv");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // YOu can use settings_. to acces the settings
    
    outputFile << "Time,Infected,Recovered,Dead,Susceptible,Exposed,Asymptomatic\n";
    mio::oseair::Model<double> model;
    set_initial_values(model, settings_);

    std::vector<double> grid(settings_.numIntervals + 1);
    double grid_spacing = (settings_.tmax - settings_.t0) / settings_.numIntervals;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        grid[i] = settings_.t0 + i * grid_spacing;
    }

    for (std::size_t controlIndex = 0; controlIndex < settings_.numControlIntervals; controlIndex++){
        double socialDistancing = x[controlIndex * settings_.numControls + 0];
        double quarantined      = x[controlIndex * settings_.numControls + 1];
        double testingRate      = x[controlIndex * settings_.numControls + 2];
        model.parameters.template get<mio::oseair::SocialDistancing<double>>() = socialDistancing;
        model.parameters.template get<mio::oseair::Quarantined<double>>()      = quarantined;
        model.parameters.template get<mio::oseair::TestingRate<double>>()      = testingRate;

        for (std::size_t i = 0; i < settings_.pcResolution ; i++)
        {
            std::size_t gridindex = controlIndex * settings_.pcResolution + i;
            auto result = mio::simulate<double, mio::oseair::Model<double>>(grid[gridindex], grid[gridindex + 1], settings_.dt, model);

            for (int j = 0; j < (int)mio::oseair::InfectionState::Count; ++j) {
                model.populations[mio::oseair::InfectionState(j)] = result.get_last_value()[j];
            }
            outputFile << grid[gridindex] << "," << result.get_last_value()[(int)mio::oseair::InfectionState::Infected] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Recovered] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Dead] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Susceptible] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Exposed] << ","
                       << result.get_last_value()[(int)mio::oseair::InfectionState::Asymptomatic] << "\n";
        }
    }
    outputFile.close();
    std::cout << "Solution saved to output.csv\n";


    // Next save the optimal control values
    std::ofstream controlFile("control.csv");
    if (!controlFile.is_open()) {
        std::cerr << "Error opening control output file!" << std::endl;
        return;
    }
    controlFile << "ControlInterval,SocialDistancing,Quarantined,TestingRate\n";
    controlFile << std::scientific << std::setprecision(6);

    for (std::size_t controlIndex = 0; controlIndex < settings_.numControlIntervals; controlIndex++){
        controlFile << controlIndex << "," 
                    << x[controlIndex * settings_.numControls + 0] << ","
                    << x[controlIndex * settings_.numControls + 1] << ","
                    << x[controlIndex * settings_.numControls + 2] << "\n";
    }
    controlFile.close();
    std::cout << "Control values saved to control.csv\n";


}


// struct ProblemSettings {
//     int numControlIntervals = 20;
//     int pcResolution = 5;
//     double t0 = 0.0;
//     double tmax = 100.0;
//     int numIntervals = numControlIntervals * pcResolution;

//     int integratorResolution = 5;
//     double dt = (tmax - t0) / (numIntervals * integratorResolution);

//     double N = 327'167'434; // total US population

//     // Vector of control bounds: { name, {lower, upper} }
//     std::vector<std::pair<std::string, std::pair<double, double>>> controlBounds = {
//         {"SocialDistancing", {0.05, 0.5}},
//         {"Quarantined",      {0.01, 0.3}},
//         {"TestingRate",      {0.15, 0.3}}
//     };
    
//     // Vector of path constraint bounds: { name, {lower, upper} }
//     std::vector<std::pair<std::string, std::pair<double, double>>> pathConstraints = {
//         {"Infections", {0.0, 1'000'000}}
//     };

//     int numControls = static_cast<int>(controlBounds.size());
//     int numPathConstraints = static_cast<int>(pathConstraints.size());
// };

