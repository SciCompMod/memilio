#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include <assert.h>
#include <epidemiology/secir.h>
#include <string>
#include <vector>

namespace epi {

/* TODO: Add more distributions here. */
typedef enum { DIST_UNIFORM } parameter_distribution;

struct parameter_info {
  std::string name; /*< The name of this parameter */
  double min_value; /*< The minumum value of this parameter */
  double max_value; /*< The maximum value of this parameter */
  parameter_distribution
      dist; /*< The statistical distribution of this parameter */
};

/* The class parameter_space_t stores ranges of parameters
 * together with information on step sizes,
 * a start and end time as well as an initial time step.
 * The class provides an iterator that iterates over all
 * generated parameter combinations.
 *
 * Currently all parameters are of type double.
 */
class parameter_space_t {
public:
  /* Constructor
   * \param [in] paramter_filename filename of a file storing ranges of input
   * parameters. Reads parameter names and values from an input file.
   */
  parameter_space_t(std::string &parameter_filename);

  /* Constructor from given SecirParams. Mainly used for testing.
   * \param [in] seir Input parameters
   * \param [in] eps  0 <= \a eps < 1
   * This will construct a parameter space with all parameters in \a seir
   * and for each parameter p in \a seir values in the range
   *  (1-\a eps) * p to (1 + \a eps) * p
   */
  parameter_space_t(const SecirParams &seir, double eps);

private:
  // A vector of all parameters with names and min/max values
  std::vector<struct parameter_info> parameters;
};

parameter_space_t::parameter_space_t(std::string &parameter_filename) {
  // TODO: implement
  assert(0 && "This function is not implemented yet.");
}

parameter_space_t::parameter_space_t(const SecirParams &seir,
                                     double eps) {
  assert(0 <= eps && eps < 1);
  const double min_factor = 1 - eps;
  const double max_factor = 1 + eps;

  /* Read all the parameters from seir and store them in our parameters list.
   * Many parameters are stored inverse in seir, so we need to reinvert them. */
  // TODO: Currently we use UNIFORM distribution for all. Change this later when
  // we know more about distributions.
  // times
  parameters.push_back({"T_inc", min_factor * 1. / seir.times.get_incubation_inv(),
                        max_factor * 1. / seir.times.get_incubation_inv(), DIST_UNIFORM});
  parameters.push_back({"T_serint", min_factor * 1. / seir.times.get_serialinterval_inv(),
                        max_factor * 1. / seir.times.get_serialinterval_inv(), DIST_UNIFORM});
  parameters.push_back({"T_infmild", min_factor * 1. / seir.times.get_infectious_mild_inv(),
                        max_factor * 1. / seir.times.get_infectious_mild_inv(), DIST_UNIFORM});
  parameters.push_back({"T_hosp2home", min_factor * 1. / seir.times.get_hospitalized_to_home_inv(),
                        max_factor * 1. / seir.times.get_hospitalized_to_home_inv(), DIST_UNIFORM});
  parameters.push_back({"T_home2hosp", min_factor * 1. / seir.times.get_home_to_hospitalized_inv(),
                        max_factor * 1. / seir.times.get_home_to_hospitalized_inv(), DIST_UNIFORM});
  parameters.push_back({"T_hosp2icu", min_factor * 1. / seir.times.get_hospitalized_to_icu_inv(),
                        max_factor * 1. / seir.times.get_hospitalized_to_icu_inv(), DIST_UNIFORM});
  parameters.push_back({"T_icu2home", min_factor * 1. / seir.times.get_icu_to_home_inv(),
                        max_factor * 1. / seir.times.get_icu_to_home_inv(), DIST_UNIFORM});
  parameters.push_back({"T_infasy", min_factor * 1. / seir.times.get_infectious_asymp_inv(),
                        max_factor * 1. / seir.times.get_infectious_asymp_inv(), DIST_UNIFORM});
  // probabilities
  parameters.push_back({"alpha", min_factor * seir.probabilities.get_asymp_per_infectious(),
                        max_factor * seir.probabilities.get_asymp_per_infectious(), DIST_UNIFORM});
  parameters.push_back({"beta", min_factor * seir.probabilities.get_risk_from_symptomatic(),
                        max_factor * seir.probabilities.get_risk_from_symptomatic(), DIST_UNIFORM});
  parameters.push_back({"rho", min_factor * seir.probabilities.get_hospitalized_per_infectious(),
                        max_factor * seir.probabilities.get_hospitalized_per_infectious(), DIST_UNIFORM});
  parameters.push_back({"theta", min_factor * seir.probabilities.get_icu_per_hospitalized(),
                        max_factor * seir.probabilities.get_icu_per_hospitalized(), DIST_UNIFORM});
  parameters.push_back({"delta", min_factor * seir.probabilities.get_dead_per_icu(),
                        max_factor * seir.probabilities.get_dead_per_icu(), DIST_UNIFORM});
}

} // namespace epi

#endif // PARAMETER_SPACE_H