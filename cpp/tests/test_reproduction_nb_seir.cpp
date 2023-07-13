#include "gtest/gtest.h"
#include <iomanip>
#include "../../cpp/models/ode_seir/reproduction_nb.h"
#include "Eigen/src/Core/Matrix.h"
#include "memilio/compartments/simulation.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/model.h"
#include "ode_seir/parameters.h"
#include "load_test_data.h"
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/math/euler.h"
#include "memilio/compartments/simulation.h"

TEST(ReproductionNumberTest, SingleNumberCalculation){

    double t0   = 0;
    double tmax = 1;
    double dt   =  0.001;

    mio::oseir::Model model;

    double total_population                                                                            = 10000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]   = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}]  = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}] = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.set<mio::oseir::TimeExposed>(5.2);
    model.parameters.set<mio::oseir::TimeInfected>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 10;

    model.check_constraints();

    double checkReproductionNumber = 2.328; //Correct result for these input values


    double coeffStoE =  4e-05;       
    Eigen::Index timept = 0;
    double TimeInfected = 6;
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);
    EXPECT_NEAR(getReproductionNumber(timept, coeffStoE, TimeInfected, result), checkReproductionNumber, 1e-12)<<"Reproduction number at time 0 is wrong";
}

TEST(ReproductionNumberTest, AllNumbersCalculation){
    double t0   = 0;
    double tmax = 1;
    double dt   =  0.001;

    mio::oseir::Model model;

     double total_population                                                                            = 10000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]   = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}]  = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}] = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.set<mio::oseir::TimeExposed>(5.2);
    model.parameters.set<mio::oseir::TimeInfected>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 10;

    model.check_constraints();

    double coeffStoE = 4e-05;
    double TimeInfected = 6;

    Eigen::VectorXd checkReproductionNumbers(7); 
    checkReproductionNumbers << 2.328,
2.3279906878991881,
2.3279487809434576,
2.3277601483151549,
2.3269102025388899,
2.323058005241374,
2.318540062468307;

    model.parameters.get<mio::oseir::TimeInfected>()=6;
    
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    for(Eigen::Index i = 0; i < getReproductionNumbers(coeffStoE, TimeInfected, result).size(); i++){
        EXPECT_NEAR(getReproductionNumbers(coeffStoE, TimeInfected, result)[i], checkReproductionNumbers[i], 1e-12)<<"Reproduction numbers are wrong at time: "<<i;
    }

}
