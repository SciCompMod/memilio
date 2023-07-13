#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/config.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "parameters.h"


double getReproductionNumber(Eigen::Index timept, double coeffStoE, double TimeInfected, mio::TimeSeries<ScalarType> y);
Eigen::VectorXd getReproductionNumbers(double coeffStoE, double TimeInfected, mio::TimeSeries<ScalarType> y);
