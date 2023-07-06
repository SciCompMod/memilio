#include "reproduction_nb.h"
#include "ode_seir/parameters.h"
#include "parameters.h"

double ReproductionNumber::getReproductionNumber(Eigen::Index timept, double coeffStoE, double TimeInfected, mio::TimeSeries<ScalarType> y){//Computes the reproduction number at a certain time (actually only needs number of susceptibles from the TimeSeries)
            return y.get_value(timept)[(Eigen::Index)mio::oseir::InfectionState::Susceptible]*TimeInfected*coeffStoE;
    }

    Eigen::VectorXd ReproductionNumber::getReproductionNumbers(double coeffStoE, double TimeInfected, mio::TimeSeries<ScalarType> y){//Computes the reproduction numbers at all times
            Eigen::VectorXd temp(y.get_num_time_points());
            for(int i = 0; i < y.get_num_time_points(); i++){
                temp[i] = coeffStoE*TimeInfected*y.get_value((Eigen::Index)i)[(Eigen::Index)mio::oseir::InfectionState::Susceptible];
            }
            return temp;
    }