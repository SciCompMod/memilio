#include "memilio/epidemiology/age_group.h"
#include "parameters.h"
#include "memilio/utils/time_series.h"
#include "Eigen/src/Eigenvalues/ComplexEigenSolver.h"
#include "ode_secir/infection_state.h"
#include "memilio/data/analyze_result.h"

namespace mio
{

double get_reproduction_number(Eigen::Index timept, mio::TimeSeries<double> y, mio::osecir::Parameters& params);
Eigen::VectorXd get_reproduction_numbers(mio::TimeSeries<double> y, mio::osecir::Parameters& params);

} // namespace mio