#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include <assert.h>
#include <epidemiology/secir.h>
#include <string>
#include <vector>
#include <random>
#include <epidemiology/logging.h>
#include <memory>

namespace epi
{

/* TODO: Add more distributions here. */
typedef enum
{
    DIST_UNIFORM,
    DIST_NORMAL
} distribution;

/*
 * Parameter Distribution class which contains the name of a variable as string
 * the lower bound and the upper bound as maximum admissible values and an enum
 * item with the name of the distribution  
 */
class ParameterDistribution
{
public:
    ParameterDistribution()
        : m_lower_bound{0}
        , m_upper_bound{0}
        , m_dist{DIST_UNIFORM}
    {
        std::random_device random_device;
        m_random_generator = std::mt19937{random_device()};
    }

    ParameterDistribution(double lower_bound, double upper_bound, distribution dist)
    {
        std::random_device random_device;
        m_random_generator = std::mt19937{random_device()};
        m_lower_bound      = lower_bound;
        m_upper_bound      = upper_bound;
        m_dist             = dist;
    }
    void set_lower_bound(double lower_bound)
    {
        m_lower_bound = lower_bound;
    }

    void set_upper_bound(double upper_bound)
    {
        m_upper_bound = upper_bound;
    }

    void set_distribution(distribution dist)
    {
        m_dist = dist;
    }

    void add_predefined_sample(double sample)
    {
        m_predefined_samples.push_back(sample);
    }

    void remove_predefined_samples()
    {
        m_predefined_samples.resize(0);
    }

    double get_lower_bound() const
    {
        return m_lower_bound;
    }

    double get_upper_bound() const
    {
        return m_upper_bound;
    }

    distribution const& get_distribution() const
    {
        return m_dist;
    }

    /*
     * @brief returns a value for the given parameter distribution
     * in case some predefined samples are set, these values are taken
     * first, in case the vector of predefined values is empty, a 'real' 
     * random sample is taken
     */
    double get_sample()
    {
        if (m_predefined_samples.size() > 0) {
            double rnumb = m_predefined_samples[0];
            m_predefined_samples.erase(m_predefined_samples.begin());
            return rnumb;
        }
        else {
            return get_rand_sample();
        }
    }

    virtual double get_rand_sample()
    {
        printf("\nMaster func");
        return 0.0;
    }

protected:
    double m_lower_bound; /*< A realistic lower bound on the given parameter */
    double m_upper_bound; /*< A realistic upper bound on the given parameter */
    distribution m_dist; /*< The statistical distribution of this parameter */
    std::mt19937 m_random_generator;
    std::vector<double>
        m_predefined_samples; // if these values are set; no real sample will occur but these values will be taken
};

/*
 * Child class of Parameter Distribution class which additionally contains
 * the mean value and the standard deviation for a normal distribution 
 */
class ParameterDistributionNormal : public ParameterDistribution
{
public:
    ParameterDistributionNormal()
        : ParameterDistribution()
    {
        m_dist         = DIST_NORMAL;
        m_mean         = 0;
        m_standard_dev = 1;
    }

    ParameterDistributionNormal(double mean, double standard_dev)
        : ParameterDistribution()
    {
        m_dist         = DIST_NORMAL;
        m_standard_dev = standard_dev;
        check_quantiles();
    }

    ParameterDistributionNormal(double lower_bound, double upper_bound, double mean, double standard_dev)
        : ParameterDistribution(lower_bound, upper_bound, DIST_NORMAL)
    {
        m_mean         = mean;
        m_standard_dev = standard_dev;
        check_quantiles();
    }

    void set_mean(double mean)
    {
        m_mean = mean;
    }

    /*
     * @brief verification that at least 99% of the density
     * function lie in the interval defined by the boundaries
     */
    bool check_quantiles()
    {
        bool changed = false;
        if (m_mean < m_lower_bound || m_mean > m_upper_bound) {
            m_mean  = 0.5 * (m_upper_bound - m_lower_bound);
            changed = true;
        }

        // ensure that 0.99 % of the distribution are within lower bound and upper bound
        if (m_mean + m_standard_dev * m_quantile > m_upper_bound) {
            m_standard_dev = (m_upper_bound - m_mean) / m_quantile;
            changed        = true;
        }
        if (m_mean - m_standard_dev * m_quantile < m_lower_bound) {
            m_standard_dev = (m_mean - m_lower_bound) / m_quantile;
            changed        = true;
        }

        if (changed && m_log_stddev_change) {
            log_warning("Standard deviation reduced to fit 95% of the distribution within [lowerbound,upperbound].");
        }

        return changed;
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

    void set_log(bool log_stddev_change)
    {
        m_log_stddev_change = log_stddev_change;
    }

    /*
     * @brief gets a sample of a normally distributed variable
     * before sampling, it is verified that at least 99% of the
     * density function lie in the interval defined by the boundaries
     * otherwise the normal distribution is adapted
     */
    double get_rand_sample() override
    {
        printf("\n child func");
        if (check_quantiles() || m_distribution.mean() != m_mean || m_distribution.stddev() != m_standard_dev) {
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
    constexpr static double m_quantile = 2.5758; // 0.995 quartile
    std::normal_distribution<double> m_distribution;
    bool m_log_stddev_change = true;
};

/*
 * Child class of Parameter Distribution class which represents an uniform distribution 
 */
class ParameterDistributionUniform : public ParameterDistribution
{
public:
    ParameterDistributionUniform()
        : ParameterDistribution()
    {
    }

    ParameterDistributionUniform(double lower_bound, double upper_bound)
        : ParameterDistribution(lower_bound, upper_bound, DIST_UNIFORM)
    {
    }

    /*
     * @brief gets a sample of a uniformly distributed variable
     */
    double get_rand_sample() override
    {
        printf("\n child func");
        if (m_distribution.max() != m_upper_bound || m_distribution.min() != m_lower_bound) {
            m_distribution = std::uniform_real_distribution<double>{m_lower_bound, m_upper_bound};
        }

        return m_distribution(m_random_generator);
    }

private:
    std::uniform_real_distribution<double> m_distribution;
};

/*
 * base class for variable elements
 * an element can be a real-valued parameter
 * or a dampings set 
 */
class VariableElement
{
public:
    /*
     * @brief initializes a VariableElement corresponding to a variable real parameter
     * @param[in] name the name of the element
     */
    VariableElement(std::string name)
        : m_nb_attr{1}
    {
        m_name = name;
    }

    /*
     * @brief initializes a VariableElement corresponding to a variable real parameter
     * @param[in] name the name of the element
     * @param[in] nb_attr number of attributes
     */
    VariableElement(std::string name, size_t nb_attr)
    {
        m_name    = name;
        m_nb_attr = nb_attr;
    }

    /*
     * @brief sets the number of attributes of the VariableElement
     */
    void set_nb_attributes(int nb_attr)
    {
        m_nb_attr = nb_attr;
    }

    /*
     * @brief returns the name of the VariableElement
     */
    std::string get_name() const
    {
        return m_name;
    }

    /*
     * @brief returns the number of attributes of the VariableElement
     */
    int get_nb_attributes() const
    {
        return m_nb_attr;
    }

protected:
    std::string m_name; // name of the model
    int m_nb_attr; // number of attributes
};

class RealVariableElement : public VariableElement
{
public:
    RealVariableElement(std::string name, std::unique_ptr<ParameterDistribution>&& distribution)
        : VariableElement(name)
    {
        m_distribution = std::move(distribution);
    }

    double get_sample()
    {
        return m_distribution->get_sample();
    }

private:
    std::unique_ptr<ParameterDistribution> m_distribution;
};

class ContactFrequencyVariableElement : public VariableElement
{
public:
    /*
     * @brief Construction of a DampingsVariableElement
     * @param[in] cont_freq contact frequency matrix
     * @param[in] nb_dampings uniform distribution on number of dampings
     * @param[in] day uniform distribution on day where one damping is implemented
     * @param[in] damp_diag_base uniform distribution on diagonal base value for one damping matrix
     * @param[in] damp_diag_rel uniform distribution for variation between diagonal values, based on the diagonal base value
     * @param[in] damp_offdiag_rel uniform distribution for variation between offdiagonal values of one line, based on the diagonal value
     */
    ContactFrequencyVariableElement(ContactFrequencyMatrix cont_freq, ParameterDistributionUniform nb_dampings,
                                    ParameterDistributionUniform day, ParameterDistributionUniform damp_diag_base,
                                    ParameterDistributionUniform damp_diag_rel,
                                    ParameterDistributionUniform damp_offdiag_rel)
        : VariableElement("ContactFrequencyMatrix", 2)
    {
        m_cont_freq        = cont_freq;
        m_nb_dampings      = nb_dampings;
        m_day              = day;
        m_damp_diag_base   = damp_diag_base;
        m_damp_diag_rel    = damp_diag_rel;
        m_damp_offdiag_rel = damp_offdiag_rel;
    }

    ContactFrequencyMatrix get_sample()
    {
        int nb_dampings = (int)(m_nb_dampings.get_sample() + 0.5);
        for (int i = 0; i < nb_dampings; i++) {

            double day            = m_day.get_sample();
            double damp_diag_base = m_damp_diag_base.get_sample();

            // diagonal entries
            std::vector<double> damp_diag_val(m_cont_freq.get_size(), 0);
            for (int j = 0; j < m_cont_freq.get_size(); j++) {
                damp_diag_val[j] = damp_diag_base * m_damp_diag_rel.get_sample();
                m_cont_freq.add_damping(Damping(day, damp_diag_val[j]), j, j);
            }

            // offdiagonal entries
            for (int j = 0; j < m_cont_freq.get_size(); j++) {

                for (int k = j + 1; k < m_cont_freq.get_size(); k++) {
                    double damp_offdiag_val = 0.5 * damp_diag_val[j] * m_damp_offdiag_rel.get_sample() +
                                              0.5 * damp_diag_val[k] * m_damp_offdiag_rel.get_sample();
                    m_cont_freq.add_damping(Damping(day, damp_diag_val[j]), j, k);
                }
            }
        }

        return m_cont_freq;
    }

private:
    ContactFrequencyMatrix m_cont_freq;
    ParameterDistributionUniform
        m_nb_dampings; // random number of dampings (one damping is understood as nb_groups^2 many dampings at the same day)
    ParameterDistributionUniform m_day; // random number of day where to implement damping
    ParameterDistributionUniform m_damp_diag_base; // random number of base value for the diagonal of the damping matrix
    ParameterDistributionUniform m_damp_diag_rel; // random number of variation from base value for diagonal
    ParameterDistributionUniform m_damp_offdiag_rel; // random number of variation from diagonal value for offdiagonal
};

// /*
//  * base class of model element variability
//  * a model element can be a real-valued parameter,
//  * a dampings of related parameters, a
//  * set of both or a combination of all of that
//  */
// class VariableElement
// {
// public:
//     /*
//      * @brief initializes a VariableElement corresponding to a variable real parameter
//      * @param[in] name the name of the element
//      * @param[in] distribution the ParameterDistribution of the single attribute of dimension one (i.e., one parameter)
//      */
//     VariableElement(std::string name, ParameterDistribution distribution)
//         : m_attr{1}
//         , m_attr_expr_eval{1}
//         , m_dimensions{1}
//         , m_max_varibility{1}
//     {
//         m_name          = name;
//         m_distributions = {1, distribution};
//     }

//     /*
//      * @brief initializes a VariableElement that can have any number of attributes, attribute samples etc.
//      * @param[in] name the name of the element
//      * @param[in] attr the number of the attributes
//      * @param[in] dimensions the dimensions per attribute
//      * @param[in] distributions the vector of ParameterDistributions of the attributes
//      * @param[in] attr_min_indep_expr minimum number of independent expression of an attribute
//      * @param[in] attr_max_indep_expr maximum number of independent expression of an attribute
//      * @param[in] max_variability the maximum variability between the parameters within a higher dimensional attribute
//      */
//     VariableElement(std::string name, int attr, std::vector<int> dimensions, int attr_min_indep_expr,
//                     int attr_max_indep_expr, std::vector<ParameterDistribution> distributions,
//                     std::vector<double> max_varibility)
//     {
//         m_name = name;
//         m_attr = attr;
//         if (attr_max_indep_expr <= 1) {
//             m_attr_expr_eval = 1;
//         }
//         else {
//             m_attr_expr.set_lower_bound(attr_min_indep_expr);
//             m_attr_expr.set_upper_bound(attr_max_indep_expr);
//         }

//         m_distributions  = distributions;
//         m_dimensions     = dimensions;
//         m_max_varibility = max_varibility;
//     }

//     /*
//      * @brief sets the number of attributes of the VariableElement
//      */
//     void set_independent_attributes(int attr)
//     {
//         m_attr = attr;
//     }

//     /*
//      * @brief sets the distributions of the independent attributes of the VariableElement
//      */
//     void set_distributions(std::vector<ParameterDistribution> distributions)
//     {
//         m_distributions = distributions;
//     }

//     /*
//      * @brief sets the dimensions of the independent attributes of the VariableElement
//      */
//     void get_dimensions(std::vector<int> dimensions)
//     {
//         m_dimensions = dimensions;
//     }
//     /*
//      * @brief returns the maximum variability in one of the independent attributes of larger dimensions of the VariableElement
//      */
//     void get_max_variability(std::vector<double> max_variability)
//     {
//         m_max_varibility = max_variability;
//     }

//     /*
//      * @brief returns the name of the VariableElement
//      */
//     std::string get_name() const
//     {
//         return m_name;
//     }

//     /*
//      * @brief returns the number of attributes of the VariableElement
//      */
//     int get_nb_attributes() const
//     {
//         return m_attr;
//     }

//     /*
//      * @brief returns the distributions of the independent attributes of the VariableElement
//      */
//     std::vector<ParameterDistribution> get_distributions() const
//     {
//         return m_distributions;
//     }

//     /*
//      * @brief returns the dimensions of the independent attributes of the VariableElement
//      */
//     std::vector<int> get_dimensions() const
//     {
//         return m_dimensions;
//     }

//     /*
//      * @brief returns the maximum variability in one of the independent attributes of larger dimensions of the VariableElement
//      */
//     std::vector<double> get_max_variability() const
//     {
//         return m_max_varibility;
//     }

// private:
//     std::string m_name; // name of the model
//     int m_attr; // number of attributes
//     ParameterDistributionUniform m_attr_expr;
//     int m_attr_expr_eval; // number of expressions per attribute; if an attribute is a matrix, two "expressions" means two matrices
//     std::vector<int> m_dimensions; // number of dimensions per attribute
//     std::vector<ParameterDistribution> m_distributions; // distributions per attribute
//     std::vector<double>
//         m_max_varibility; // maximum variability in independent attribute (only applies for multivalued attributes)
// };

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
    parameter_space_t(ContactFrequencyMatrix const& cont_freq_matrix, std::vector<SecirParams> const& params, double t0,
                      double tmax);

private:
    // A vector of all parameters with names and min/max values
    std::unordered_map<std::string, VariableElement> parameters;
};

parameter_space_t::parameter_space_t(std::string& parameter_filename)
{
    // TODO: implement
    assert(0 && "This function is not implemented yet.");
}

parameter_space_t::parameter_space_t(ContactFrequencyMatrix const& cont_freq_matrix,
                                     std::vector<SecirParams> const& params, double t0, double tmax)
{

    RealVariableElement a{"incubation time",
                          std::unique_ptr<ParameterDistribution>(new ParameterDistributionNormal(1, 14, 5.2, 3))};

    a.get_sample();

    // maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    // damping base values are between 0.1 and 1; diagonal values vary lie in the range of 0.6 to 1.4 times the base value
    // off diagonal values vary between 0.7 to 1.1 of the corresponding diagonal value (symmetrization is conducted)
    ContactFrequencyVariableElement bbb{cont_freq_matrix,
                                        ParameterDistributionUniform(1, (tmax - t0) / 10),
                                        ParameterDistributionUniform(t0, tmax),
                                        ParameterDistributionUniform(0.1, 1),
                                        ParameterDistributionUniform(0.6, 1.4),
                                        ParameterDistributionUniform(0.7, 1.1)};

    bbb.get_sample();

    /* Read all the parameters from seir and store them in our parameters list.
   * Many parameters are stored inverse in seir, so we need to reinvert them. */
    // TODO: Currently we use UNIFORM distribution for all. Change this later when
    // we know more about distributions.
    // times
    // parameters.push_back({"T_inc", min_factor * 1. / seir.times.get_incubation_inv(),
    //                       max_factor * 1. / seir.times.get_incubation_inv(), DIST_UNIFORM});
    // parameters.push_back({"T_serint", min_factor * 1. / seir.times.get_serialinterval_inv(),
    //                       max_factor * 1. / seir.times.get_serialinterval_inv(), DIST_UNIFORM});
    // parameters.push_back({"T_infmild", min_factor * 1. / seir.times.get_infectious_mild_inv(),
    //                       max_factor * 1. / seir.times.get_infectious_mild_inv(), DIST_UNIFORM});
    // parameters.push_back({"T_hosp2home", min_factor * 1. / seir.times.get_hospitalized_to_home_inv(),
    //                       max_factor * 1. / seir.times.get_hospitalized_to_home_inv(), DIST_UNIFORM});
    // parameters.push_back({"T_home2hosp", min_factor * 1. / seir.times.get_home_to_hospitalized_inv(),
    //                       max_factor * 1. / seir.times.get_home_to_hospitalized_inv(), DIST_UNIFORM});
    // parameters.push_back({"T_hosp2icu", min_factor * 1. / seir.times.get_hospitalized_to_icu_inv(),
    //                       max_factor * 1. / seir.times.get_hospitalized_to_icu_inv(), DIST_UNIFORM});
    // parameters.push_back({"T_icu2home", min_factor * 1. / seir.times.get_icu_to_home_inv(),
    //                       max_factor * 1. / seir.times.get_icu_to_home_inv(), DIST_UNIFORM});
    // parameters.push_back({"T_infasy", min_factor * 1. / seir.times.get_infectious_asymp_inv(),
    //                       max_factor * 1. / seir.times.get_infectious_asymp_inv(), DIST_UNIFORM});
    // // probabilities
    // parameters.push_back({"infprob", min_factor * seir.probabilities.get_infection_from_contact(),
    //                       max_factor * seir.probabilities.get_infection_from_contact(), DIST_UNIFORM});
    // parameters.push_back({"alpha", min_factor * seir.probabilities.get_asymp_per_infectious(),
    //                       max_factor * seir.probabilities.get_asymp_per_infectious(), DIST_UNIFORM});
    // parameters.push_back({"beta", min_factor * seir.probabilities.get_risk_from_symptomatic(),
    //                       max_factor * seir.probabilities.get_risk_from_symptomatic(), DIST_UNIFORM});
    // parameters.push_back({"rho", min_factor * seir.probabilities.get_hospitalized_per_infectious(),
    //                       max_factor * seir.probabilities.get_hospitalized_per_infectious(), DIST_UNIFORM});
    // parameters.push_back({"theta", min_factor * seir.probabilities.get_icu_per_hospitalized(),
    //                       max_factor * seir.probabilities.get_icu_per_hospitalized(), DIST_UNIFORM});
    // parameters.push_back({"delta", min_factor * seir.probabilities.get_dead_per_icu(),
    //                       max_factor * seir.probabilities.get_dead_per_icu(), DIST_UNIFORM});
}

} // namespace epi

#endif // PARAMETER_SPACE_H