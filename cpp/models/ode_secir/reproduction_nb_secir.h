#include "memilio/epidemiology/age_group.h"
#include "parameters.h"
#include "memilio/utils/time_series.h"
#include "Eigen/src/Eigenvalues/ComplexEigenSolver.h"
#include "ode_secir/infection_state.h"
#include "memilio/data/analyze_result.h"




double getReproductionNumber(Eigen::Index timept, mio::TimeSeries<double> y, mio::osecir::Parameters &params);
Eigen::VectorXd getReproductionNumbers(mio::TimeSeries<double>y, mio::osecir::Parameters &params);  
