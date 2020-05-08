#ifndef SECIR_H
#define SECIR_H

#include <epidemiology/damping.h>

#include <vector>

namespace epi
{

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
class SecirParams
{
public:
    double base_reprod;
    double cont_freq, tinc_inv, tinfmild_inv; // parameters of the standard SEIR model
    double tserint_inv, thosp2home_inv, thome2hosp_inv, thosp2icu_inv, ticu2home_inv, tinfasy_inv, ticu2death_inv;
    double alpha, beta, rho, theta, delta; // probabilities

    // population parameters of unit scale
    class Populations
    {
    public:
        /**
         * @brief Initializes a time parameters' struct of the SECIR model
         */
        Populations();

        /**
         * @brief sets the number of total people at t0 for the SECIR model
         * automatically calls set_suscetible_t0() to subtract from the total number
         * @param nb_total_t0 total number of people at t0
         */
        void set_total_t0(double nb_total_t0);

        /**
         * @brief sets the number of exposed people at t0 for the SECIR model
         * automatically calls set_suscetible_t0() to subtract from the total number
         * @param nb_exp_t0 number of exposed people at t0
         */
        void set_exposed_t0(double nb_exp_t0);

        /**
         * @brief sets the number of carrier people at t0 for the SECIR model
         * automatically calls set_suscetible_t0() to subtract from the total number
         * @param nb_car_t0 number of recovered people at t0
         */
        void set_carrier_t0(double nb_car_t0);

        /**
         * @brief sets the number of infectious people at t0 for the SECIR model
         * automatically calls set_suscetible_t0() to subtract from the total number
         * @param nb_inf_t0 number of infectious people at t0
         */
        void set_infectious_t0(double nb_inf_t0);

        /**
         * @brief sets the number of hospitalized people at t0 for the SECIR model
         * automatically calls set_suscetible_t0() to subtract from the total number
         * @param nb_hosp_t0 number of recovered people at t0
         */
        void set_hospital_t0(double nb_hosp_t0);

        /**
         * @brief sets the number of ICU-treated people at t0 for the SECIR model
         * automatically calls set_suscetible_t0() to subtract from the total number
         * @param nb_icu_t0 number of recovered people at t0
         */
        void set_icu_t0(double nb_icu_t0);

        /**
         * @brief sets the number of recovered people at t0 for the SECIR model
         * automatically calls set_suscetible_t0() to subtract from the total number
         * @param nb_rec_t0 number of recovered people at t0
         */
        void set_recovered_t0(double nb_rec_t0);

        /**
         * @brief sets the number of dead people at t0 for the SECIR model
         * automatically calls set_suscetible_t0() to subtract from the total number
         * @param nb_dead_t0 number of recovered people at t0
         */
        void set_dead_t0(double nb_dead_t0);

        /**
         * @brief sets the number of suscetible people at t0 for the SECIR model
         * only to be called after all other populations have been called
         */
        void set_suscetible_t0();

        /**
         * @brief returns the number of total people at t0 for the SECIR model
         */
        double get_total_t0() const;

        /**
         * @brief returns the number of exposed people at t0 for the SECIR model
         */
        double get_exposed_t0() const;

        /**
         * @brief returns the number of carrier people at t0 for the SECIR model
         */
        double get_carrier_t0() const;

        /**
         * @brief returns the number of infectious people at t0 for the SECIR model
         */
        double get_infectious_t0() const;

        /**
         * @brief returns the number of hospitalized people at t0 for the SECIR model
         */
        double get_hospitalized_t0() const;

        /**
         * @brief returns the number of ICU-treated people at t0 for the SECIR model
         */
        double get_icu_t0() const;

        /**
         * @brief returns the number of recovered people at t0 for the SECIR model
         */
        double get_recovered_t0() const;

        /**
         * @brief returns the number of dead people at t0 for the SECIR model
         */
        double get_dead_t0() const;

        /**
         * @brief returns the number of suscetible people at t0 for the SECIR model
         */
        double get_suscetible_t0() const;

    private:
        double m_nb_total_t0, m_nb_sus_t0, m_nb_exp_t0, m_nb_car_t0, m_nb_inf_t0, m_nb_hosp_t0, m_nb_icu_t0,
            m_nb_rec_t0, m_nb_dead_t0;
    };

    Populations populations;

    // This defines a damping factor for a mitigation strategy for different points in time.
    Dampings dampings;

    /**
     * @brief Initializes a SECIR/SECIHURD model without default parameters 
     */
    SecirParams();

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
     * @param cont_freq_in
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
    SecirParams(double tinc, double tinfmild, double tserint, double thosp2home, double thome2hosp, double thosp2icu,
                double ticu2home, double tinfasy, double ticu2death, double cont_freq_in, double alpha_in,
                double beta_in, double delta_in, double rho_in, double theta_in);
};

/**
 * prints given parameters
 * @param[in] params the SecirParams parameter object
 */
void print_secir_params(SecirParams const& params);

/**
 * Computes the current time-derivative of S, E, C, I, (H, U,) R, (D) in the SECIR/SECIHURD model
 * @param[in] params SECIR/SECIHURD Model parameters, created by secir_param
 * @tparam T the datatype of the cases
 * @param[in] y current  S, E, C, I, (H, U,) R, (D) values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, C, I, (H, U,) R, (D)
 */
void secir_get_derivatives(SecirParams const& params, std::vector<double> const& y, double t,
                           std::vector<double>& dydt);

/**
 * Computes the SECIR curve by integration
 * @param[in] secir_0 Initial S, E, C, I, (H, U,) R, (D) values at t0
 * @param[in] t0 start time of simulation
 * @param[in] tmax end time of simulation
 * @param[in] dt initial time step
 * @param[in] params SECIR/SECIHURD model parameters
 *
 * @returns Vector of times t
 */
std::vector<double> simulate(double t0, double tmax, double dt, SecirParams const& params,
                             std::vector<std::vector<double>>& secir);

} // namespace epi

#endif // SECIR_H
