#include "gtest/gtest.h"
#include <iomanip>
#include "../../cpp/models/ode_secir/reproduction_nb_secir.h"
#include "Eigen/src/Core/Matrix.h"
#include "load_test_data.h"
#include "ode_secir/model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/math/stepper_wrapper.h"

TEST(SecirReproductionNumberTest, SecirAllNumbersCalculation)
{ //Tests the function getReproductionNumbers()
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model model(3);
    mio::AgeGroup nb_groups = model.parameters.get_num_groups();
    double fact             = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::StartDay>(60);
    params.set<mio::osecir::Seasonality>(0.2);
    params.get<mio::osecir::TestAndTraceCapacity>() = 35;

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.8;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 9.5;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 7.1;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i]  = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]    = 0.7;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]    = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]    = 0.25;
        params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[i] = 0.45;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]         = 0.2;
        params.get<mio::osecir::CriticalPerSevere>()[i]                 = 0.25;
        params.get<mio::osecir::DeathsPerCritical>()[i]                 = 0.3;
    }

    params.apply_constraints();

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    auto integrator = std::make_shared<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    mio::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model, integrator);

    Eigen::VectorXd check_reproduction_numbers;

    check_reproduction_numbers << 2.1638667159763307, 2.1632954864514815, 2.1571571556004301, 2.1504039811502995,
        2.1431719654877202, 2.1355248862531422, 2.1274913816514456, 2.1190832031240379, 2.1103042347730474,
        2.1011549831023335, 2.0916348409435193, 2.0817432543239063, 2.0714803376012867, 2.0608465992296456,
        2.0498425176270807, 2.038469375669048, 2.0267293901902863, 2.0146257521044664, 2.0021626671156745,
        1.9893453973851463, 1.9761803044253157, 1.9626748932796978, 1.9488378577795786, 1.9346791263776075,
        1.9202099077584158, 1.9054427351259804, 1.8903915077733227, 1.8750715282602988, 1.8594995332682944,
        1.8436937159777183, 0.54828792974181884, 0.54544311652706057, 0.54401554130951646, 0.54266867306450506,
        0.54142368536764063, 0.54028399823827622, 0.53924527706280578, 0.53830112000268104, 0.53744545595568516,
        0.53667164779084897, 0.5359728831634043, 0.53534252888392309, 0.53477430966292916, 0.53426238824239447,
        0.53380139075112232, 0.53338640224931155, 0.5330129466704745, 0.53267695922386515, 0.53237475580277338,
        0.53210300192645732, 0.5318586825849374, 0.53165999117327567;

    for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
        EXPECT_NEAR(get_reproduction_numbers(secihurd, params)[i], check_reproduction_numbers[i], 3e-12);
    }
}