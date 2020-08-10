#ifndef SECIR_H
#define SECIR_H

#include <epidemiology/populations.h>
#include <epidemiology/adapt_rk.h>
#include <epidemiology/uncertain_value.h>
#include <epidemiology/uncertain_matrix.h>

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
    SecirParams(size_t nb_groups = 1)
        : contact_patterns(ContactFrequencyMatrix{nb_groups})
        , populations(Populations({nb_groups, SecirCount}))
    {
        times         = std::vector<StageTimes>(nb_groups, StageTimes());
        probabilities = std::vector<Probabilities>(nb_groups, Probabilities());
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
        const UncertainValue& get_incubation() const;
        UncertainValue& get_incubation();

        /**
         * @brief returns 1.0 over the infectious time set for the SECIR model in day unit
         */
        const UncertainValue& get_infectious_mild() const;
        UncertainValue& get_infectious_mild();

        /**
         * @brief returns 1.0 over the serial interval in the SECIR model
         */
        const UncertainValue& get_serialinterval() const;
        UncertainValue& get_serialinterval();

        /**
         * @brief returns 1.0 over the time people are 'simply' hospitalized before returning home in the SECIR model 
         */
        const UncertainValue& get_hospitalized_to_home() const;
        UncertainValue& get_hospitalized_to_home();

        /**
         * @brief returns 1.0 over the time people are infectious at home before 'simply' hospitalized in the SECIR model 
         */
        const UncertainValue& get_home_to_hospitalized() const;
        UncertainValue& get_home_to_hospitalized();

        /**
         * @brief returns 1.0 over the time people are 'simply' hospitalized before being treated by ICU in the SECIR model
         */
        const UncertainValue& get_hospitalized_to_icu() const;
        UncertainValue& get_hospitalized_to_icu();

        /**
         * @brief returns 1.0 over the time people are treated by ICU before returning home in the SECIR model
         */
        const UncertainValue& get_icu_to_home() const;
        UncertainValue& get_icu_to_home();

        /**
         * @brief returns 1.0 over the infectious time for asymptomatic cases in the SECIR model
         */
        const UncertainValue& get_infectious_asymp() const;
        UncertainValue& get_infectious_asymp();

        /**
         * @brief returns 1.0 over the time people are treated by ICU before dying in the SECIR model
         */
        const UncertainValue& get_icu_to_dead() const;
        UncertainValue& get_icu_to_dead();

    private:
        UncertainValue m_tinc, m_tinfmild; // parameters also available in SEIR
        UncertainValue m_tserint, m_thosp2home, m_thome2hosp, m_thosp2icu, m_ticu2home, m_tinfasy,
            m_ticu2death; // new SECIR params
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
        const UncertainValue& get_infection_from_contact() const;
        UncertainValue& get_infection_from_contact();

        /**
        * @brief returns the percentage of asymptomatic cases in the SECIR model
        */
        const UncertainValue& get_asymp_per_infectious() const;
        UncertainValue& get_asymp_per_infectious();

        /**
        * @brief returns the risk of infection from symptomatic cases in the SECIR model
        */
        const UncertainValue& get_risk_from_symptomatic() const;
        UncertainValue& get_risk_from_symptomatic();

        /**
        * @brief returns the percentage of hospitalized patients per infected patients in the SECIR model
        */
        const UncertainValue& get_hospitalized_per_infectious() const;
        UncertainValue& get_hospitalized_per_infectious();

        /**
        * @brief returns the percentage of ICU patients per hospitalized patients in the SECIR model
        */
        const UncertainValue& get_icu_per_hospitalized() const;
        UncertainValue& get_icu_per_hospitalized();

        /**
        * @brief returns the percentage of dead patients per ICU patients in the SECIR model
        */
        const UncertainValue& get_dead_per_icu() const;
        UncertainValue& get_dead_per_icu();

    private:
        UncertainValue m_infprob, m_asympinf, m_risksymp, m_hospinf, m_icuhosp, m_deathicu; // probabilities
    };

    /**
     * @brief returns the contact frequency matrix set in UncertainContactMatrix
     */
    ContactFrequencyMatrix& get_cont_freq_matrix();

    /**
     * @brief returns the contact frequency matrix set in UncertainContactMatrix
     */
    ContactFrequencyMatrix const& get_cont_freq_matrix() const;

    UncertainContactMatrix contact_patterns;
    Populations populations;
    std::vector<StageTimes> times;
    std::vector<Probabilities> probabilities;
};

/**
 * @brief WIP !! TO DO: returns the actual, approximated reproduction rate 
 */
double get_reprod_rate(SecirParams const& params, double t, std::vector<double> const& yt);

/**
 * Computes the current time-derivative of the SECIR compartment populations (e.g., S, E, C, I, H, U, R, D)
 * @param[in]  params SECIR params: contact frequencies, population sizes and epidemiological parameters
 * @param[in] y current  SECIR compartments' values at t; (e.g., y: [0:S, 1:E, ...])
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatives of the SECIR compartments (e.g., S, E, C, I, H, U, R, D)
 */
void secir_get_derivatives(SecirParams const& params, Eigen::VectorXd const& y, double t, Eigen::VectorXd& dydt);

/**
 * @brief simulate SECIR compartment model.
 * The simulation supports multiple groups with different parameters interacting with each other.
 */
class SecirSimulation
{
public:
    /**
     * @brief setup the SECIR simulation 
     * @param[in] params SECIR params: contact frequencies, population sizes and epidemiological parameters
     * @param[in] t0 start time
     * @param[in] dt initial step size of integration
     */
    SecirSimulation(const SecirParams& params, double t0 = 0., double dt = 0.1);

    /**
     * @brief advance simulation to tmax
     * must be greater than get_t().back()
     * @param tmax next stopping point of simulation
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

    /**
     * @brief returns the SECIR params used in simulation
     */
    const SecirParams& get_params() const
    {
        return m_params;
    }

private:
    std::shared_ptr<RKIntegratorCore> m_integratorCore;
    OdeIntegrator m_integrator;
    SecirParams m_params;
};

/**
 * @brief run secir simulation over fixed time
 * @param[in] t0 start time
 * @param[in] tmax end time
 * @param[in] dt initial step size of integration
 * @param[in] params SECIR params: contact frequencies, population sizes and epidemiological parameters
 * @param[out] secir value of compartments at each integration time point
 * @returns integration time points
 */
std::vector<double> simulate(double t0, double tmax, double dt, SecirParams const& params,
                             std::vector<Eigen::VectorXd>& secir);

} // namespace epi

#endif // SECIR_H
