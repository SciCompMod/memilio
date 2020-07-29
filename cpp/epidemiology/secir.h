#ifndef SECIR_H
#define SECIR_H

#include <epidemiology/populations.h>
#include <epidemiology/damping.h>
#include <epidemiology/adapt_rk.h>
#include <epidemiology/uncertain_value.h>

#include <vector>
#include <Eigen/Core>

namespace epi
{

enum SecirCategory
{
    AgeGroup,
    InfectionType,
    CategoryCount
};

enum SecirCompartments
{
    S,
    E,
    C,
    I,
    H,
    U,
    R,
    D,
    SecirCount
};

/**
 * @brief Initializes a Contact Frequency matrix for the SECIR/SECIHURD model
 *
 * @todo parameter description
 *
 * @param m_cont_freq
 **/
class ContactFrequencyMatrix
{
public:
    /**
     * @brief Standard constructor of contact frequencies 1x1-matrix in the SECIR model
     */
    ContactFrequencyMatrix();

    /**
     * @brief Constructor of contact frequencies nb_groups x nb_groups-matrix in the SECIR model
     * @param[in] nb_groups number of groups in the model
     */
    ContactFrequencyMatrix(size_t const nb_groups);

    /**
     * @brief returns the size of the contact frequency matrix
     */
    int get_size() const;

    /**
     * @brief sets the contact frequency in the SECIR model; in case of multiple groups, set the contact rate cr_ij=cr_ji=cont_freq
     * @param cont_freq contact rate/frequency in 1/day unit
     * @param self_group own group
     * @param contact_group group which gets in contact with own group
     */
    void set_cont_freq(double cont_freq, int self_group, int contact_group);

    /**
     * @brief returns the contact frequency set for the SECIR model in 1/day unit; in case of multiple groups, returns the contact rate cr_ij=cr_ji
     */
    double get_cont_freq(int self_group, int contact_group) const;

    /**
     * @brief sets the damping in the SECIR model; in case of multiple groups, set the contact rate d_ij=d_ji=cont_freq
     * @param damping dampings over the whole time line in day unit
     * @param self_group own group
     * @param contact_group group which gets in contact with own group
     */
    void set_dampings(Dampings const& damping, int self_group, int contact_group);

    /**
     * @brief returns the dampings set for the SECIR model in 1/day unit; in case of multiple groups, returns the damping d_ij=d_ji
     */
    const Dampings& get_dampings(int self_group, int contact_group) const;

    /**
     * @brief add damping to the dampings object specified by self_ and contact_group
     * @param damping one damping in day unit
     * @param self_group own group
     * @param contact_group group which gets in contact with own group
     */
    void add_damping(Damping const& damping, int self_group, int contact_group);

private:
    std::vector<std::vector<double>> m_cont_freq;
    // This defines a damping factor for a mitigation strategy for different points in time.
    std::vector<std::vector<Dampings>> m_dampings;
}; // namespace epi

/**
 * Paramters of the SECIR/SECIHURD model:
 * T_inc (also sigma^(-1) or R_2^(-1)+R_3^(-1)): mean incubation period (default: 5.2);
 *          R_2^(-1) is the first part of the incubation time where the person is not yet infectioous
 *          R_3 is the exchange between asymptomatic carriers and infectious people; R_3^(-1) is the second part of the incubation time where the person is infectious WITHOUT showing symptoms
 * T_serint (also R_2^(-1)+0.5*R_3^(-1)): serial interval (default: 4.2);
 * T_infmild (also gamma^(-1) or R_4^(-1)): time a person remains infective after disease (if 'hospitalized' is considered a state, it does not apply to them but only to 'mildly infected' people in SECIR)
 * T_hosp2home (also R_5^(1)): duration for which the hospitalized patients not requiring further intensive care remain under general hospital care (=INF or R_5=0 in standard SEIR to waive influence of this parameter)
 * T_home2hosp (also R_6^(-1)): mean time a patient with mild symptoms spends at home before hospital admission due to worsening of the disease condition  (=INF or R_6=0 in standard SEIR to waive influence of this parameter)
 * T_hosp2icu (also R_7^(-1)): mean time a patient who entered the hospital will be hopistalized without ICU before being connected to an ICU  (=INF or R_7=0 in standard SEIR to waive influence of this parameter)
 * T_icu2home (also R_8^(-1)): mean time a patient is connected to an ICU before returning home (=INF or R_8=0 in standard SEIR to waive influence of this parameter)
 * T_infasy (also R_9^(-1)): mean time an asymptomatic person remains infective (=INF or R_9=0 in standard SEIR to waive influence of this parameter)
 * T_icu2death (also d; better would be R_10^(-1)): mean time a person needs ICU support before dying (=INF or R_10=0 in standard SEIR to waive influence of this parameter)
 * cont_freq (also R_1: contact frequency/rate; called beta in the standard SEIR model)
 * alpha: share of asymptomatic cases
 * beta (Not the beta in SEIR model): risk of infection from the infected symptomatic patients
 * rho: H/I; hospitalized per infected (=0 in standard SEIR)
 * theta: U/H; intensive care units per hospitalized
 * delta: D/U; deaths per intensive care units
**/
/**
     * @brief Initializes a SECIR/SECIHURD model
     *
     * @todo parameter description
     *
     * @param tinc
     * @param tinfmild
     * @param tserint
     * @param thosp2home
     * @param thome2hosp
     * @param thosp2icu
     * @param ticu2home
     * @param tinfasy
     * @param ticu2death
     * @param alpha_in
     * @param beta_in
     * @param delta_in
     * @param rho_in
     * @param theta_in
     * @param nb_total_t0_in
     * @param nb_exp_t0_in
     * @param nb_car_t0_in
     * @param nb_inf_t0_in
     * @param nb_hosp_t0_in
     * @param nb_icu_t0_in
     * @param nb_rec_t0_in
     * @param nb_dead_t0_in
     */
class SecirParams
{
public:
    SecirParams(size_t n_groups = 1)
        : populations(Populations({n_groups, SecirCount}))
    {
        times         = std::vector<StageTimes>(n_groups, StageTimes());
        probabilities = std::vector<Probabilities>(n_groups, Probabilities());
    }

    size_t size() const
    {
        return times.size();
    }

    double base_reprod;

    // time parameters for the different 'stages' of the disease of scale day or 1/day
    // 'stages' does not refer to the 'states' of the SECIR model but also includes incubation time or contact frequency
    class StageTimes
    {
    public:
        /**
         * @brief Standard constructor of a time parameters' class in the SECIR model
         */
        StageTimes();
        StageTimes(StageTimes&&)      = default;
        StageTimes(const StageTimes&) = default;
        StageTimes& operator=(const StageTimes&) = default;

        /**
         * @brief sets the incubation time in the SECIR model
         * @param tinc incubation time in day unit
         */
        void set_incubation(double const& tinc);

        /**
         * @brief sets the infectious time for symptomatic cases that are infected but who do not need to be hsopitalized in the SECIR model
         * @param tinfmild infectious time for symptomatic cases (if not hospitalized) in day unit 
         */
        void set_infectious_mild(double const& tinfmild);

        /**
         * @brief sets the serial interval in the SECIR model
         * @param tserint serial interval in day unit 
         */
        void set_serialinterval(double const& tserint);

        /**
         * @brief sets the time people are 'simply' hospitalized before returning home in the SECIR model
         * @param thosp2home time people are 'simply' hospitalized before returning home in day unit 
         */
        void set_hospitalized_to_home(double const& thosp2home);

        /**
         * @brief sets the time people are infectious at home before 'simply' hospitalized in the SECIR model
         * @param thome2hosp time people are infectious at home before 'simply' hospitalized in day unit 
         */
        void set_home_to_hospitalized(double const& thome2hosp);

        /**
         * @brief sets the time people are 'simply' hospitalized before being treated by ICU in the SECIR model
         * @param thosp2icu time people are 'simply' hospitalized before being treated by ICU in day unit 
         */
        void set_hospitalized_to_icu(double const& thosp2icu);

        /**
         * @brief sets the time people are treated by ICU before returning home in the SECIR model
         * @param ticu2home time people are treated by ICU before returning home in day unit 
         */
        void set_icu_to_home(double const& ticu2home);

        /**
         * @brief sets the infectious time for asymptomatic cases in the SECIR model
         * @param tinfasy infectious time for asymptomatic cases in day unit 
         */
        void set_infectious_asymp(double const& tinfasy);

        /**
         * @brief sets the time people are treated by ICU before dying in the SECIR model
         * @param ticu2death time people are treated by ICU before dying in day unit 
         */
        void set_icu_to_death(double const& ticu2death);

        /**
         * @brief returns 1.0 over the incubation time set for the SECIR model in day unit
         */
        const UncertainValue& get_incubation_inv() const;
        UncertainValue& get_incubation_inv();

        /**
         * @brief returns 1.0 over the infectious time set for the SECIR model in day unit
         */
        const UncertainValue& get_infectious_mild_inv() const;
        UncertainValue& get_infectious_mild_inv();

        /**
         * @brief returns 1.0 over the serial interval in the SECIR model
         */
        const UncertainValue& get_serialinterval_inv() const;
        UncertainValue& get_serialinterval_inv();

        /**
         * @brief returns 1.0 over the time people are 'simply' hospitalized before returning home in the SECIR model 
         */
        const UncertainValue& get_hospitalized_to_home_inv() const;
        UncertainValue& get_hospitalized_to_home_inv();

        /**
         * @brief returns 1.0 over the time people are infectious at home before 'simply' hospitalized in the SECIR model 
         */
        const UncertainValue& get_home_to_hospitalized_inv() const;
        UncertainValue& get_home_to_hospitalized_inv();

        /**
         * @brief returns 1.0 over the time people are 'simply' hospitalized before being treated by ICU in the SECIR model
         */
        const UncertainValue& get_hospitalized_to_icu_inv() const;
        UncertainValue& get_hospitalized_to_icu_inv();

        /**
         * @brief returns 1.0 over the time people are treated by ICU before returning home in the SECIR model
         */
        const UncertainValue& get_icu_to_home_inv() const;
        UncertainValue& get_icu_to_home_inv();

        /**
         * @brief returns 1.0 over the infectious time for asymptomatic cases in the SECIR model
         */
        const UncertainValue& get_infectious_asymp_inv() const;
        UncertainValue& get_infectious_asymp_inv();

        /**
         * @brief returns 1.0 over the time people are treated by ICU before dying in the SECIR model
         */
        const UncertainValue& get_icu_to_dead_inv() const;
        UncertainValue& get_icu_to_dead_inv();

    private:
        UncertainValue m_tinc_inv, m_tinfmild_inv; // parameters also available in SEIR
        UncertainValue m_tserint_inv, m_thosp2home_inv, m_thome2hosp_inv, m_thosp2icu_inv, m_ticu2home_inv,
            m_tinfasy_inv,
            m_ticu2death_inv; // new SECIR params
    };

    class Probabilities
    {
    public:
        /**
         * @brief Standard constructor of probabilites parameters' class in the SECIR model
         */
        Probabilities();
        Probabilities(const Probabilities&) = default;
        Probabilities(Probabilities&&)      = default;
        Probabilities& operator=(const Probabilities&) = default;

        /**
        * @brief sets probability of getting infected from a contact
        * @param infprob the probability of getting infected from a contact
        */
        void set_infection_from_contact(double const& infprob);

        /**
        * @brief sets the percentage of asymptomatic cases in the SECIR model
        * @param alpha the percentage of asymptomatic cases
        */
        void set_asymp_per_infectious(double const& alpha);

        /**
        * @brief sets the risk of infection from symptomatic cases in the SECIR model
        * @param beta the risk of infection from symptomatic cases 
        */
        void set_risk_from_symptomatic(double const& beta);

        /**
        * @brief sets the percentage of hospitalized patients per infected patients in the SECIR model
        * @param rho percentage of hospitalized patients per infected patients
        */
        void set_hospitalized_per_infectious(double const& rho);

        /**
        * @brief sets the percentage of ICU patients per hospitalized patients in the SECIR model
        * @param theta percentage of ICU patients per hospitalized patients
        */
        void set_icu_per_hospitalized(double const& theta);

        /**
        * @brief sets the percentage of dead patients per ICU patients in the SECIR model
        * @param delta percentage of dead patients per ICU patients 
        */
        void set_dead_per_icu(double const& delta);

        /**
        * @brief gets probability of getting infected from a contact
        */
        double get_infection_from_contact() const;

        /**
        * @brief returns the percentage of asymptomatic cases in the SECIR model
        */
        double get_asymp_per_infectious() const;

        /**
        * @brief returns the risk of infection from symptomatic cases in the SECIR model
        */
        double get_risk_from_symptomatic() const;

        /**
        * @brief returns the percentage of hospitalized patients per infected patients in the SECIR model
        */
        double get_hospitalized_per_infectious() const;

        /**
        * @brief returns the percentage of ICU patients per hospitalized patients in the SECIR model
        */
        double get_icu_per_hospitalized() const;

        /**
        * @brief returns the percentage of dead patients per ICU patients in the SECIR model
        */
        double get_dead_per_icu() const;

    private:
        double m_infprob, m_alpha, m_beta, m_rho, m_theta, m_delta; // probabilities
    };

    Populations populations;
    std::vector<StageTimes> times;
    std::vector<Probabilities> probabilities;
};

/**
 * @brief returns the actual, approximated reproduction rate 
 */
double get_reprod_rate(ContactFrequencyMatrix const& cont_freq_matrix, SecirParams const& params, double t,
                       std::vector<double> const& yt);

/**
 * prints given parameters
 * @param[in] params the SecirParams parameter object
 */
void print_secir_params(ContactFrequencyMatrix const& cont_freq, SecirParams const& params);

/**
 * Computes the current time-derivative of S, E, C, I, (H, U,) R, (D) in the SECIR/SECIHURD model
 * @param[in] params SECIR/SECIHURD Model parameters, created by secir_param
 * @tparam T the datatype of the cases
 * @param[in] y current  S, E, C, I, (H, U,) R, (D) values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, C, I, (H, U,) R, (D)
 */
void secir_get_derivatives(ContactFrequencyMatrix const& cont_freq_matrix, SecirParams const& params,
                           Eigen::VectorXd const& y, double t, Eigen::VectorXd& dydt);

/**
 * @brief simulate SECIR compartment model.
 * The simulation supports multiple groups with different parameters interacting with each other.
 */
class SecirSimulation
{
public:
    /**
     * @brief setup the SECIR simulation
     * @param cont_freq_matrix contact frequencies between groups
     * @param params parameters of each group
     * @param t0 start time
     * @param dt_init initial step size of integration
     */
    SecirSimulation(const ContactFrequencyMatrix& cont_freq_matrix, const SecirParams& params, double t0 = 0.,
                    double dt_init = 0.1);

    /**
     * @brief advance simulation to tmax
     * must be greater than get_t().back()
     */
    Eigen::VectorXd& advance(double tmax);

    /**
     * @brief the integration time points
     */
    const std::vector<double>& get_t() const
    {
        return m_integrator.get_t();
    }

    /**
     * @brief values of compartments at each time point
     */
    const std::vector<Eigen::VectorXd>& get_y() const
    {
        return m_integrator.get_y();
    }
    std::vector<Eigen::VectorXd>& get_y()
    {
        return m_integrator.get_y();
    }

    const ContactFrequencyMatrix& get_cont_freq() const
    {
        return m_cont_freq;
    }

    const SecirParams& get_params() const
    {
        return m_params;
    }

private:
    std::shared_ptr<RKIntegratorCore> m_integratorCore;
    OdeIntegrator m_integrator;
    ContactFrequencyMatrix m_cont_freq;
    SecirParams m_params;
};

/**
 * @brief run secir simulation over fixed time
 * @param cont_freq_matrix contact frequencies between groups
 * @param params parameters of each group
 * @param t0 start time
 * @param tmax end time
 * @param dt_init initial step size of integration
 * @param[out] secir value of compartments at each integration time point
 * @returns integration time points
 */
std::vector<double> simulate(double t0, double tmax, double dt, ContactFrequencyMatrix const& cont_freq_matrix,
                             SecirParams const& params, std::vector<Eigen::VectorXd>& secir);

} // namespace epi

#endif // SECIR_H
