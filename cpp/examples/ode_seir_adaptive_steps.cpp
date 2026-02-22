#include "memilio/compartments/flow_simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"
#include "ode_seir/model.h"
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>

namespace mio
{
namespace examples
{

inline void build_seir_matrix_A_val(Eigen::MatrixXd& A, const mio::oseir::Model<double>& model,
                                    const Eigen::VectorXd& z, double t)
{
    using namespace mio::oseir;
    const size_t num_groups = (size_t)model.parameters.get_num_groups();

    Eigen::VectorXd force_of_infection = Eigen::VectorXd::Zero(num_groups);
    for (size_t j = 0; j < num_groups; ++j) {
        const size_t Sj =
            model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), InfectionState::Susceptible});
        const size_t Ej =
            model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), InfectionState::Exposed});
        const size_t Ij =
            model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), InfectionState::Infected});
        const size_t Rj =
            model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), InfectionState::Recovered});

        const double Nj    = z[Sj] + z[Ej] + z[Ij] + z[Rj];
        const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;

        for (size_t i = 0; i < num_groups; ++i) {
            const double coeffStoE =
                model.parameters.template get<ContactPatterns<double>>().get_cont_freq_mat().get_matrix_at(t)(
                    static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) *
                model.parameters
                    .template get<TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(static_cast<int>(i))] *
                divNj;
            force_of_infection[i] += coeffStoE * z[Ij];
        }
    }

    A.setZero();
    for (size_t g = 0; g < num_groups; ++g) {
        double time_exposed  = model.parameters.get<TimeExposed<double>>()[mio::AgeGroup(g)];
        double time_infected = model.parameters.get<TimeInfected<double>>()[mio::AgeGroup(g)];

        double sigma  = (time_exposed > 1e-10) ? (1.0 / time_exposed) : 0.0;
        double gamma  = (time_infected > 1e-10) ? (1.0 / time_infected) : 0.0;
        double lambda = force_of_infection[g];

        size_t iS = model.populations.get_flat_index({mio::AgeGroup(g), InfectionState::Susceptible});
        size_t iE = model.populations.get_flat_index({mio::AgeGroup(g), InfectionState::Exposed});
        size_t iI = model.populations.get_flat_index({mio::AgeGroup(g), InfectionState::Infected});
        size_t iR = model.populations.get_flat_index({mio::AgeGroup(g), InfectionState::Recovered});

        // S -> E
        A(iS, iS) = -lambda;
        A(iE, iS) = lambda;

        // E -> I
        A(iE, iE) = -sigma;
        A(iI, iE) = sigma;

        // I -> R
        A(iI, iI) = -gamma;
        A(iR, iI) = gamma;
    }
}

struct AugmentedPhiSystemVal {
    const mio::oseir::Model<double>& model;
    size_t NC;
    Eigen::MatrixXd A_mem;

    AugmentedPhiSystemVal(const mio::oseir::Model<double>& m)
        : model(m)
        , NC((size_t)m.populations.get_num_compartments())
    {
        A_mem = Eigen::MatrixXd::Zero(NC, NC);
    }

    void operator()(const Eigen::VectorXd& y, Eigen::VectorXd& dydt, double t)
    {
        // 1. z and dz
        const auto z = y.head(NC);
        Eigen::VectorXd dz(NC);
        model.get_derivatives(z, z, t, dz);
        dydt.head(NC) = dz;

        // 2. A(z)
        build_seir_matrix_A_val(A_mem, model, z, t);

        // 3. Phi Dynamik
        const auto Phi       = Eigen::Map<const Eigen::MatrixXd>(y.tail(NC * NC).data(), NC, NC);
        Eigen::MatrixXd dPhi = A_mem * Phi;
        dydt.tail(NC * NC)   = Eigen::Map<const Eigen::VectorXd>(dPhi.data(), NC * NC);
    }
};

struct StandardSystemVal {
    const mio::oseir::Model<double>& model;

    StandardSystemVal(const mio::oseir::Model<double>& m)
        : model(m)
    {
    }

    void operator()(const Eigen::VectorXd& y, Eigen::VectorXd& dydt, double t)
    {
        model.get_derivatives(y, y, t, dydt);
    }
};

void setup_model_validation(mio::oseir::Model<double>& model)
{
    // Same as validation file
    double total_population = 10000.0;
    size_t num_age_groups   = (size_t)model.parameters.get_num_groups();

    for (size_t i = 0; i < num_age_groups; ++i) {
        model.parameters.get<mio::oseir::TimeExposed<double>>()[mio::AgeGroup(i)]                      = 5.2;
        model.parameters.get<mio::oseir::TimeInfected<double>>()[mio::AgeGroup(i)]                     = 6.0;
        model.parameters.get<mio::oseir::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(i)] = 0.05;

        // Initial pops
        model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Susceptible}] =
            0.9 * total_population / num_age_groups;
        model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Exposed}] = 0.0;
        model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Infected}] =
            0.1 * total_population / num_age_groups;
        model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Recovered}] = 0.0;
    }

    // Contact matrix
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_age_groups, num_age_groups, 10.0));
}

} // namespace examples
} // namespace mio

// Observer to track steps
struct StepObserver {
    double last_t;
    int step_count;
    std::vector<double> dts;

    StepObserver()
        : last_t(0.0)
        , step_count(0)
    {
    }

    void operator()(const Eigen::VectorXd& /*x*/, double t)
    {
        if (step_count > 0) {
            double dt = t - last_t;
            // The step observer is called after each step.
            if (dt > 1e-12) {
                dts.push_back(dt);
            }
        }
        last_t = t;
        step_count++;
    }

    void print_stats()
    {
        if (dts.empty())
            return;

        double min_dt = dts[0];
        double max_dt = dts[0];
        double sum_dt = 0;

        for (double dt : dts) {
            if (dt < min_dt)
                min_dt = dt;
            if (dt > max_dt)
                max_dt = dt;
            sum_dt += dt;
        }
        double avg = sum_dt / dts.size();

        std::cout << "Steps: " << dts.size() << "\n"
                  << "Min dt: " << min_dt << "\n"
                  << "Max dt: " << max_dt << "\n"
                  << "Avg dt: " << avg << std::endl;

        // Print first 10 and last 10 steps
        std::cout << "First 10 steps: ";
        for (size_t i = 0; i < std::min((size_t)10, dts.size()); ++i)
            std::cout << dts[i] << " ";
        std::cout << "\nLast 10 steps: ";
        for (size_t i = (dts.size() > 10 ? dts.size() - 10 : 0); i < dts.size(); ++i)
            std::cout << dts[i] << " ";
        std::cout << std::endl;
    }
};

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    const size_t num_age_groups = 10;
    const double t0             = 0.0;
    const double t_max          = 10.0;

    std::cout << "Comparison of Adaptive Step Sizes" << std::endl
              << "AgeGroups: " << num_age_groups << ", T_max: " << t_max << std::endl;

    mio::oseir::Model<double> model(num_age_groups);
    mio::examples::setup_model_validation(model);

    // Initial state setup
    size_t NC = (size_t)model.populations.get_num_compartments();

    // 1. Matrix-Phi System (NC + NC*NC variables)
    size_t system_size_phi = NC + NC * NC;
    Eigen::VectorXd y0_phi(system_size_phi);
    y0_phi.head(NC)      = model.populations.get_compartments();
    Eigen::MatrixXd Id   = Eigen::MatrixXd::Identity(NC, NC);
    y0_phi.tail(NC * NC) = Eigen::Map<Eigen::VectorXd>(Id.data(), NC * NC);

    // 2. Standard System (NC variables)
    Eigen::VectorXd y0_standard = model.populations.get_compartments();

    typedef boost::numeric::odeint::runge_kutta_cash_karp54<Eigen::VectorXd, double, Eigen::VectorXd, double,
                                                            boost::numeric::odeint::vector_space_algebra>
        error_stepper_type;

    // 1. Matrix-Phi Adaptive Simulation
    {
        std::cout << "\n--- Phi ---" << std::endl;

        auto stepper = boost::numeric::odeint::make_controlled<error_stepper_type>(1e-10, 1e-5 // abs, rel tol
        );

        mio::examples::AugmentedPhiSystemVal sys(model);
        StepObserver obs;
        Eigen::VectorXd y = y0_phi;

        boost::numeric::odeint::integrate_adaptive(stepper, sys, y, t0, t_max, 10.0, std::ref(obs));

        obs.print_stats();
    }

    // 2. Standard SEIR Adaptive Simulation for the totals
    {
        std::cout << "\n--- Stage Aligned ---" << std::endl;

        auto stepper = boost::numeric::odeint::make_controlled<error_stepper_type>(1e-10, 1e-5);

        mio::examples::StandardSystemVal sys(model);
        StepObserver obs;
        Eigen::VectorXd y = y0_standard;

        boost::numeric::odeint::integrate_adaptive(stepper, sys, y, t0, t_max, 10., std::ref(obs));

        obs.print_stats();
    }

    return 0;
}
