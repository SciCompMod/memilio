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
        m_mean         = mean;
        m_standard_dev = standard_dev;
        check_quantiles(m_mean, m_standard_dev);
    }

    ParameterDistributionNormal(double lower_bound, double upper_bound, double mean)
        : ParameterDistribution(lower_bound, upper_bound, DIST_NORMAL)
    {
        m_mean         = mean;
        m_standard_dev = upper_bound; // set as to high and adapt then
        adapt_standard_dev(m_standard_dev);
    }

    ParameterDistributionNormal(double lower_bound, double upper_bound, double mean, double standard_dev)
        : ParameterDistribution(lower_bound, upper_bound, DIST_NORMAL)
    {
        m_mean         = mean;
        m_standard_dev = standard_dev;
        check_quantiles(m_mean, m_standard_dev);
    }

    void set_mean(double mean)
    {
        m_mean = mean;
    }

    bool check_quantiles()
    {
        return check_quantiles(m_mean, m_standard_dev);
    }

    /*
     * @brief verification that at least 99% of the density
     * function lie in the interval defined by the boundaries
     */
    bool check_quantiles(double& mean, double& standard_dev)
    {
        bool changed = false;

        changed = adapt_mean(mean);

        changed = adapt_standard_dev(standard_dev);

        if (changed && m_log_stddev_change) {
            log_warning("Standard deviation reduced to fit 99% of the distribution within [lowerbound,upperbound].");
        }

        return changed;
    }

    bool adapt_mean(double& mean)
    {
        bool changed = false;
        if (mean < m_lower_bound || mean > m_upper_bound) {
            mean    = 0.5 * (m_upper_bound - m_lower_bound);
            changed = true;
        }
        return changed;
    }

    // ensure that 0.99 % of the distribution are within lower bound and upper bound
    bool adapt_standard_dev(double& standard_dev)
    {
        bool changed = false;
        if (m_mean + standard_dev * m_quantile > m_upper_bound) {
            standard_dev = (m_upper_bound - m_mean) / m_quantile;
            changed      = true;
        }
        if (m_mean - standard_dev * m_quantile < m_lower_bound) {
            standard_dev = (m_mean - m_lower_bound) / m_quantile;
            changed      = true;
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
        if (check_quantiles(m_mean, m_standard_dev) || m_distribution.mean() != m_mean ||
            m_distribution.stddev() != m_standard_dev) {
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
        if (m_distribution.max() != m_upper_bound || m_distribution.min() != m_lower_bound) {
            m_distribution = std::uniform_real_distribution<double>{m_lower_bound, m_upper_bound};
        }

        return m_distribution(m_random_generator);
    }

private:
    std::uniform_real_distribution<double> m_distribution;
};

class ParamsVariableElement
{
public:
    /*
     * @brief creates a ParamsVariableElement from a string name and a distribution via unique_ptr
     * @param[in] distribution unique pointer to a distribution
     */
    ParamsVariableElement(std::vector<std::unique_ptr<ParameterDistribution>>&& distribution)
    {
        m_distribution = std::move(distribution);
    }

    std::vector<double> get_sample()
    {
        std::vector<double> samples(m_distribution.size(), 0);
        for (size_t i = 0; i < m_distribution.size(); i++) {
            samples[i] = m_distribution[i]->get_sample();
        }
        return samples;
    }

private:
    std::vector<std::unique_ptr<ParameterDistribution>> m_distribution;
};

class ContactFrequencyVariableElement
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
    ContactFrequencyVariableElement(ContactFrequencyMatrix cont_freq,
                                    std::unique_ptr<ParameterDistributionUniform>&& nb_dampings,
                                    std::unique_ptr<ParameterDistributionUniform>&& day,
                                    std::unique_ptr<ParameterDistributionUniform>&& damp_diag_base,
                                    std::unique_ptr<ParameterDistributionUniform>&& damp_diag_rel,
                                    std::unique_ptr<ParameterDistributionUniform>&& damp_offdiag_rel)
    {
        m_cont_freq        = cont_freq;
        m_nb_dampings      = std::move(nb_dampings);
        m_day              = std::move(day);
        m_damp_diag_base   = std::move(damp_diag_base);
        m_damp_diag_rel    = std::move(damp_diag_rel);
        m_damp_offdiag_rel = std::move(damp_offdiag_rel);
    }

    ContactFrequencyMatrix get_sample()
    {
        int nb_dampings = (int)(m_nb_dampings->get_sample() + 0.5);
        for (int i = 0; i < nb_dampings; i++) {

            double day            = m_day->get_sample();
            double damp_diag_base = m_damp_diag_base->get_sample();

            // diagonal entries
            std::vector<double> damp_diag_val(m_cont_freq.get_size(), 0);
            for (int j = 0; j < m_cont_freq.get_size(); j++) {
                damp_diag_val[j] = damp_diag_base * m_damp_diag_rel->get_sample();
                m_cont_freq.add_damping(Damping(day, damp_diag_val[j]), j, j);
            }

            // offdiagonal entries
            for (int j = 0; j < m_cont_freq.get_size(); j++) {

                for (int k = j + 1; k < m_cont_freq.get_size(); k++) {
                    double damp_offdiag_val = 0.5 * damp_diag_val[j] * m_damp_offdiag_rel->get_sample() +
                                              0.5 * damp_diag_val[k] * m_damp_offdiag_rel->get_sample();
                    m_cont_freq.add_damping(Damping(day, damp_offdiag_val), j, k);
                }
            }
        }

        return m_cont_freq;
    }

private:
    ContactFrequencyMatrix m_cont_freq;
    std::unique_ptr<ParameterDistributionUniform>
        m_nb_dampings; // random number of dampings (one damping is understood as nb_groups^2 many dampings at the same day)
    std::unique_ptr<ParameterDistributionUniform> m_day; // random number of day where to implement damping
    std::unique_ptr<ParameterDistributionUniform>
        m_damp_diag_base; // random number of base value for the diagonal of the damping matrix
    std::unique_ptr<ParameterDistributionUniform>
        m_damp_diag_rel; // random number of variation from base value for diagonal
    std::unique_ptr<ParameterDistributionUniform>
        m_damp_offdiag_rel; // random number of variation from diagonal value for offdiagonal
};

/* The class ParameterSpace stores ranges of parameters
 * together with information on step sizes,
 * a start and end time as well as an initial time step.
 * The class provides an iterator that iterates over all
 * generated parameter combinations.
 *
 * Currently all parameters are of type double.
 */
class ParameterSpace
{
public:
    /* Constructor
     * @param [in] paramter_filename filename of a file storing ranges of input
     * parameters. Reads parameter names and values from an input file.
     */
    ParameterSpace(std::string& parameter_filename);

    /* Constructor from given ContactFrequencyMatrix and SecirParams. Mainly used for testing.
     * @param[in] cont_freq_matrix basic contact frequency matrix
     * @param[in] params SecirParams for alle age groups
     * @param[in] t0 start time
     * @param[in] tmax end time
     * @param[in] dev_rel maximum relative deviation from particular value(s) given in params
     */
    ParameterSpace(ContactFrequencyMatrix const& cont_freq_matrix, std::vector<SecirParams> const& params, double t0,
                   double tmax, double dev_rel);

    ContactFrequencyMatrix get_cont_freq_matrix_sample()
    {
        return m_cont_freq_matrix_variable->get_sample();
    }

    std::vector<SecirParams> get_secir_params_sample()
    {
        std::vector<SecirParams> secir_params_vec_sample;
        for (size_t i = 0; i < m_nb_age_groups; i++) {
            SecirParams secir_params_sample;

            secir_params_sample.populations.set_total_t0(m_total[i]);
            secir_params_sample.populations.set_hospital_t0(m_hospitalized[i]);
            secir_params_sample.populations.set_icu_t0(m_icu[i]);
            secir_params_sample.populations.set_dead_t0(m_dead[i]);

            secir_params_sample.populations.set_exposed_t0(m_exposed[i]->get_sample());
            secir_params_sample.populations.set_carrier_t0(m_carrier[i]->get_sample());
            secir_params_sample.populations.set_infectious_t0(m_infectious[i]->get_sample());
            secir_params_sample.populations.set_recovered_t0(m_recovered[i]->get_sample());

            secir_params_sample.populations.set_suscetible_t0(); // calculated as remaing from total minus other groups

            double inc_dummy    = m_incubation[i]->get_sample();
            double serint_dummy = inc_dummy - m_serial_int_incub_diff[i]->get_sample();

            secir_params_sample.times.set_incubation(inc_dummy);
            secir_params_sample.times.set_infectious_mild(m_inf_mild[i]->get_sample());
            secir_params_sample.times.set_serialinterval(serint_dummy);
            secir_params_sample.times.set_hospitalized_to_home(m_hosp_to_rec[i]->get_sample()); // here: home=recovered
            secir_params_sample.times.set_home_to_hospitalized(m_inf_to_hosp[i]->get_sample()); // here: home=infectious
            secir_params_sample.times.set_infectious_asymp(m_inf_asymp[i]->get_sample());
            secir_params_sample.times.set_hospitalized_to_icu(m_hosp_to_icu[i]->get_sample());
            secir_params_sample.times.set_icu_to_death(m_icu_to_death[i]->get_sample());
            secir_params_sample.times.set_icu_to_home(m_icu_to_rec[i]->get_sample());

            secir_params_sample.probabilities.set_infection_from_contact(m_inf_from_cont[i]->get_sample());
            secir_params_sample.probabilities.set_asymp_per_infectious(m_asymp_per_inf[i]->get_sample());
            secir_params_sample.probabilities.set_risk_from_symptomatic(m_risk_from_symp[i]->get_sample());
            secir_params_sample.probabilities.set_dead_per_icu(m_death_per_icu[i]->get_sample());
            secir_params_sample.probabilities.set_hospitalized_per_infectious(m_hosp_per_inf[i]->get_sample());
            secir_params_sample.probabilities.set_icu_per_hospitalized(m_icu_per_hosp[i]->get_sample());

            secir_params_vec_sample.push_back(std::move(secir_params_sample));
        }

        return secir_params_vec_sample;
    }

private:
    size_t m_nb_age_groups;

    // contact frequency matrix
    std::unique_ptr<ContactFrequencyVariableElement> m_cont_freq_matrix_variable;

    // populations
    std::vector<double> m_total;
    std::vector<std::unique_ptr<ParameterDistribution>> m_exposed;
    std::vector<std::unique_ptr<ParameterDistribution>> m_carrier;
    std::vector<std::unique_ptr<ParameterDistribution>> m_infectious;
    std::vector<double> m_hospitalized;
    std::vector<double> m_icu;
    std::vector<std::unique_ptr<ParameterDistribution>> m_recovered;
    std::vector<double> m_dead;

    // times
    std::vector<std::unique_ptr<ParameterDistribution>> m_incubation;
    std::vector<std::unique_ptr<ParameterDistribution>> m_serial_int_incub_diff;
    std::vector<std::unique_ptr<ParameterDistribution>> m_inf_mild;
    std::vector<std::unique_ptr<ParameterDistribution>> m_hosp_to_rec;
    std::vector<std::unique_ptr<ParameterDistribution>> m_inf_to_hosp;
    std::vector<std::unique_ptr<ParameterDistribution>> m_inf_asymp;
    std::vector<std::unique_ptr<ParameterDistribution>> m_hosp_to_icu;
    std::vector<std::unique_ptr<ParameterDistribution>> m_icu_to_rec;
    std::vector<std::unique_ptr<ParameterDistribution>> m_icu_to_death;

    // probabilities
    std::vector<std::unique_ptr<ParameterDistribution>> m_inf_from_cont;
    std::vector<std::unique_ptr<ParameterDistribution>> m_asymp_per_inf;
    std::vector<std::unique_ptr<ParameterDistribution>> m_risk_from_symp;
    std::vector<std::unique_ptr<ParameterDistribution>> m_death_per_icu;
    std::vector<std::unique_ptr<ParameterDistribution>> m_hosp_per_inf;
    std::vector<std::unique_ptr<ParameterDistribution>> m_icu_per_hosp;
};

ParameterSpace::ParameterSpace(std::string& parameter_filename)
{
    // TODO: implement
    assert(0 && "This function not implemented yet and needs a file read method.");
}

ParameterSpace::ParameterSpace(ContactFrequencyMatrix const& cont_freq_matrix, std::vector<SecirParams> const& params,
                               double t0, double tmax, double dev_rel)
    : m_nb_age_groups{params.size()}
{
    double min_val = 0.001;

    // populations
    for (size_t i = 0; i < m_nb_age_groups; i++) {

        // fixed size groups
        // total
        m_total.push_back(params[i].populations.get_total_t0());
        m_hospitalized.push_back(params[i].populations.get_hospitalized_t0());
        m_icu.push_back(params[i].populations.get_icu_t0());
        m_dead.push_back(params[i].populations.get_dead_t0());

        // variably sized groups
        // exposed
        double value_params = params[i].populations.get_exposed_t0();
        m_exposed.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // carrier
        value_params = params[i].populations.get_carrier_t0();
        m_carrier.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // infectious
        value_params = params[i].populations.get_infectious_t0();
        m_infectious.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // recovered
        value_params = params[i].populations.get_recovered_t0();
        m_recovered.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));
    }

    // times
    for (size_t i = 0; i < m_nb_age_groups; i++) {
        // incubation time
        double value_params = 1.0 / params[i].times.get_incubation_inv();
        m_incubation.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // serial interval
        value_params = 1.0 / params[i].times.get_incubation_inv() - 1.0 / params[i].times.get_serialinterval_inv();
        m_serial_int_incub_diff.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // infectious time (mild)
        value_params = 1.0 / params[i].times.get_infectious_mild_inv();
        m_inf_mild.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // infectious to recovered (hospitalized to home)
        value_params = 1.0 / params[i].times.get_hospitalized_to_home_inv();
        m_hosp_to_rec.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // infectious to hospitalized (home to hospitalized)
        value_params = 1.0 / params[i].times.get_home_to_hospitalized_inv();
        m_inf_to_hosp.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // infectious (asymptomatic)
        value_params = 1.0 / params[i].times.get_infectious_asymp_inv();
        m_inf_asymp.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // hospitalized to ICU
        value_params = 1.0 / params[i].times.get_hospitalized_to_icu_inv();
        m_hosp_to_icu.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // ICU to recovered
        value_params = 1.0 / params[i].times.get_icu_to_home_inv();
        m_icu_to_rec.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // ICU to death
        value_params = 1.0 / params[i].times.get_icu_to_dead_inv();
        m_icu_to_death.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));
    }

    // probabilities
    for (size_t i = 0; i < m_nb_age_groups; i++) {
        // infection from contact
        double value_params = params[i].probabilities.get_infection_from_contact();
        m_inf_from_cont.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // asymptomatic per infectious
        value_params = params[i].probabilities.get_asymp_per_infectious();
        m_asymp_per_inf.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // risk of infection from infectious
        value_params = params[i].probabilities.get_risk_from_symptomatic();
        m_risk_from_symp.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // deaths per icu treatments
        value_params = params[i].probabilities.get_dead_per_icu();
        m_death_per_icu.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // hospitalized per infections
        value_params = params[i].probabilities.get_hospitalized_per_infectious();
        m_hosp_per_inf.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // icu treatments per hospitalized
        value_params = params[i].probabilities.get_icu_per_hospitalized();
        m_icu_per_hosp.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));
    }

    // maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    // damping base values are between 0.1 and 1; diagonal values vary lie in the range of 0.6 to 1.4 times the base value
    // off diagonal values vary between 0.7 to 1.1 of the corresponding diagonal value (symmetrization is conducted)
    m_cont_freq_matrix_variable = std::move(std::make_unique<ContactFrequencyVariableElement>(
        cont_freq_matrix, std::make_unique<ParameterDistributionUniform>(1, (tmax - t0) / 10),
        std::make_unique<ParameterDistributionUniform>(t0, tmax),
        std::make_unique<ParameterDistributionUniform>(0.1, 1),
        std::make_unique<ParameterDistributionUniform>(0.6, 1.4),
        std::make_unique<ParameterDistributionUniform>(0.7, 1.1)));
}

} // namespace epi

#endif // PARAMETER_SPACE_H