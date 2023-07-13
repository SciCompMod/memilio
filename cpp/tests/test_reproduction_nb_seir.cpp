#include <gtest/gtest.h>
#include "../../cpp/models/ode_seir/reproduction_nb.h"
#include "Eigen/src/Core/Matrix.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/model.h"


TEST(ReproductionNumberTest, SingleNumberCalculation){
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::oseir::Model model;

    double checkReproductionNumber = 2.328; //Correct result for these input values

    double coeffStoE =  0;       
    Eigen::Index timept = 0;
    double TimeInfected = 6;
    mio::TimeSeries<double> y = mio::simulate(t0, tmax, dt, model);
    EXPECT_NEAR(getReproductionNumber(timept, coeffStoE, TimeInfected, y), checkReproductionNumber, 1e-12);
}

TEST(ReproductionNumberTest, AllNumbersCalculation){
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::oseir::Model model;
    double coeffStoE = 1;
    double TimeInfected = 1;

    Eigen::VectorXd checkReproductionNumbers(7); 
    checkReproductionNumbers << 2.328,
2.32799,
2.32795,
2.32776,
2.32691,
2.32306,
2.31854;

    mio::TimeSeries<double> y = mio::simulate(t0, tmax, dt, model);

    Eigen::VectorXd temp = getReproductionNumbers(coeffStoE, TimeInfected, y);

    for(Eigen::Index i = 0; i < temp.size(); i++){
        EXPECT_NEAR(temp[i], checkReproductionNumbers[i], 1e-12);
    }

}