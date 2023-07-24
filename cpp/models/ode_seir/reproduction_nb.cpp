#include "reproduction_nb.h"
#include "ode_seir/parameters.h"
#include "parameters.h"

ScalarType get_reproduction_number(Eigen::Index timept, ScalarType coeffStoE, ScalarType TimeInfected,
                                   mio::TimeSeries<ScalarType> y)
{ //Computes the reproduction number at a certain time (actually only needs number of susceptibles from the TimeSeries)
    return y.get_value(timept)[(Eigen::Index)mio::oseir::InfectionState::Susceptible] * TimeInfected * coeffStoE;
}

Eigen::VectorXd get_reproduction_numbers(ScalarType coeffStoE, ScalarType TimeInfected, mio::TimeSeries<ScalarType> y)
{
    auto num_time_points = y.get_num_time_points();
    Eigen::VectorXd temp(num_time_points);
    for (int i = 0; i < num_time_points; i++) {
        temp[i] = coeffStoE * TimeInfected *
                  y.get_value((Eigen::Index)i)[(Eigen::Index)mio::oseir::InfectionState::Susceptible];
    }
    return temp;
}