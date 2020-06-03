#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include <assert.h>
#include <epidemiology/secir.h>
#include <string>
#include <vector>
#include <random>
#include <epidemiology/logging.h>

namespace epi
{

/* TODO: Add more distributions here. */
typedef enum
{
    DIST_UNIFORM,
    DIST_NORMAL
} parameter_distribution;

class ParameterDistribution
{
public:
    ParameterDistribution()
    {
        std::random_device random_device;
        m_random_generator = std::mt19937{random_device()};
    }

    ParameterDistribution(std::string name, double lower_bound, double upper_bound, parameter_distribution dist)
    {
        std::random_device random_device;
        m_random_generator = std::mt19937{random_device()};
        m_name             = name;
        m_lower_bound      = lower_bound;
        m_upper_bound      = upper_bound;
        m_dist             = dist;
    }

    void set_name(std::string name)
    {
        m_name = name;
    }

    void set_lower_bound(double lower_bound)
    {
        m_lower_bound = lower_bound;
    }

    void set_upper_bound(double upper_bound)
    {
        m_upper_bound = upper_bound;
    }

    void set_distribution(parameter_distribution dist)
    {
        m_dist = dist;
    }

    std::string const& get_name() const
    {
        return m_name;
    }

    double get_lower_bound() const
    {
        return m_lower_bound;
    }

    double get_upper_bound() const
    {
        return m_upper_bound;
    }

    parameter_distribution const& get_distribution() const
    {
        return m_dist;
    }

    virtual double get_sample_point()
    {
    }

protected:
    std::string m_name; /*< The name of this parameter */
    double m_lower_bound; /*< A realistic lower bound on the given parameter */
    double m_upper_bound; /*< A realistic upper bound on the given parameter */
    parameter_distribution m_dist; /*< The statistical distribution of this parameter */
    std::mt19937 m_random_generator;
};

class ParameterDistributionNormal : ParameterDistribution
{
public:
    void set_mean(double mean)
    {
        m_mean = mean;
    }

    void set_standard_dev(double standard_dev)
    {
        m_standard_dev = standard_dev;
    }

    double get_mean() const
    {
        return m_mean;
    }

    double get_standard_dev() const
    {
        return m_standard_dev;
    }

    /*
     * @brief gets a sample of a normally distributed variable
     */
    double get_sample_point()
    {
        if (m_distribution.mean() != m_mean || m_distribution.stddev() != m_standard_dev) {
            m_distribution = std::normal_distribution<double>{m_mean, m_standard_dev};
        }

        int i        = 0;
        int retries  = 10;
        double rnumb = m_distribution(m_random_generator);
        while ((rnumb > m_upper_bound || rnumb < m_lower_bound) && i < retries) {
            rnumb = m_distribution(m_random_generator);
            i++;
            if (i == retries) {
                log_warning("Not successfully sampled within [min,max].");
                if (rnumb > m_upper_bound) {
                    rnumb = m_upper_bound;
                }
                else {
                    rnumb = m_lower_bound;
                }
            }
        }
        return rnumb;
    }

private:
    double m_mean; // the mean value of the normal distribution
    double m_standard_dev; // the standard deviation of the normal distribution
    std::normal_distribution<double> m_distribution;
};

class ParameterDistributionUniform : ParameterDistribution
{
public:
    /*
     * @brief gets a sample of a uniformly distributed variable
     */
    double get_sample_point()
    {
        if (m_distribution.max() != m_upper_bound || m_distribution.min() != m_lower_bound) {
            m_distribution = std::uniform_real_distribution<double>{m_lower_bound, m_upper_bound};
        }

        return m_distribution(m_random_generator);
    }

private:
    std::uniform_real_distribution<double> m_distribution;
};

/* The class parameter_space_t stores ranges of parameters
 * together with information on step sizes,
 * a start and end time as well as an initial time step.
 * The class provides an iterator that iterates over all
 * generated parameter combinations.
 *
 * Currently all parameters are of type double.
 */
class parameter_space_t
{
public:
    /* Constructor
   * \param [in] paramter_filename filename of a file storing ranges of input
   * parameters. Reads parameter names and values from an input file.
   */
    parameter_space_t(std::string& parameter_filename);

    /* Constructor from given SecirParams. Mainly used for testing.
   * \param [in] seir Input parameters
   * \param [in] eps  0 <= \a eps < 1
   * This will construct a parameter space with all parameters in \a seir
   * and for each parameter p in \a seir values in the range
   *  (1-\a eps) * p to (1 + \a eps) * p
   */
    parameter_space_t(const SecirParams& seir, double eps);

private:
    // A vector of all parameters with names and min/max values
    std::vector<ParameterDistribution> parameters;
};

parameter_space_t::parameter_space_t(std::string& parameter_filename)
{
    // TODO: implement
    assert(0 && "This function is not implemented yet.");
}

parameter_space_t::parameter_space_t(const SecirParams& seir, double eps)
{
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
    parameters.push_back({"infprob", min_factor * seir.probabilities.get_infection_from_contact(),
                          max_factor * seir.probabilities.get_infection_from_contact(), DIST_UNIFORM});
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