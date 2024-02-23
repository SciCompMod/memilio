// #include "Eigen/src/Core/util/Meta.h"
#include "boost/numeric/odeint/stepper/controlled_runge_kutta.hpp"
// #include "load_test_data.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/model.h"
#include "memilio/math/adapt_rk.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/analyze_result.h"
#include "ode_secir/parameters.h"

#include "boost/fusion/functional/invocation/invoke.hpp"
#include "memilio/io/result_io.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>
#include <string>

int main()
{
    // General set up.
    ScalarType t0   = 0;
    ScalarType tmax = 10.00;
    ScalarType dt   = 1.;

    /**********************************
    *         ODE simulation          *
    **********************************/

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 20, nb_car_t0 = 20, nb_inf_t0 = 3, nb_hosp_t0 = 1, nb_icu_t0 = 1,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model model_ode(1);

    // Set parameters
    ScalarType cont_freq = 1.0;

    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix
    model_ode.parameters.set<mio::osecir::Seasonality>(0.0);
    mio::ContactMatrixGroup& contact_matrix = model_ode.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model_ode.parameters.get<mio::osecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // Parameters needed to determine transition rates
    model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] = 2.6;
    model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] = 2.0;

    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 0.3;
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 0.3;
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 0.3;

    // Set initial values for compartments
    model_ode.populations.set_total(nb_total_t0);
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]                     = nb_exp_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]          = nb_car_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = nb_inf_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = nb_hosp_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = nb_icu_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = nb_rec_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = nb_dead_t0;
    model_ode.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                    nb_total_t0);

    // Set probabilities that determine proportion between compartments
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.5;
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 0.5;
    model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 0.5;
    model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.5;

    // Further model parameters
    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 1.0;
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 1.0;
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 1.0;
    // Choose TestAndTraceCapacity very large so that riskFromInfectedSymptomatic = RiskOfInfectionFromSymptomatic
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity>() = std::numeric_limits<ScalarType>::max();
    // Choose ICUCapacity very large so that CriticalPerSevereAdjusted = CriticalPerSevere and deathsPerSevereAdjusted = 0
    model_ode.parameters.get<mio::osecir::ICUCapacity>() = std::numeric_limits<ScalarType>::max();

    model_ode.check_constraints();

    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // choose dt_min = dt_max so that we have a fixed time step and can compare to IDE
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);

    integrator->set_rel_tolerance(1e-5);
    integrator->set_abs_tolerance(1e-5);

    mio::TimeSeries<ScalarType> secihurd_ode = simulate(t0, tmax, dt, model_ode, integrator);

    printf("\n # t");

    auto num_points = static_cast<size_t>(secihurd_ode.get_num_time_points());
    for (size_t i = 0; i < num_points; i++) {
        printf("\n%.14f ", secihurd_ode.get_time(i));
    }
    std::cout << "\n";
}