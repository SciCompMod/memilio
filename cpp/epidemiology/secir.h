#ifndef SECIR_H
#define SECIR_H

#include <vector>

/**
 * This defined a damping factor for a
 * mitigation strategy for one point in time.
 */
class Damping
{
public:
    double day;
    double factor;

    Damping(double day_in, double factor_in);
};

/**
 * Paramters of the model(s):
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
 * cont_freq (also R_1: contact frequency
 * alpha: share of asymptomatic cases
 * beta (Not the beta in SEIR model): risk of infection from the infected symptomatic patients
 * rho: H/I; hospitalized per infected (=0 in standard SEIR)
 * theta: U/H; intensive care units per hospitalized
 * delta: D/U; deaths per intensive care units
**/
class SecirParams
{
public:
    int model;

    double base_reprod, b, cont_freq;
    double tinc_inv, tinfmild_inv;
    double tserint_inv, thosp2home_inv, thome2hosp_inv, thosp2icu_inv, ticu2home_inv, tinfasy_inv, ticu2death_inv;
    double alpha, beta, rho, theta, delta;

    // double nb_total, nb_exp, nb_car, nb_inf, nb_hosp, nb_icu, nb_rec, nb_dead;
    double nb_total_t0, nb_sus_t0, nb_exp_t0, nb_car_t0, nb_inf_t0, nb_hosp_t0, nb_icu_t0, nb_rec_t0, nb_dead_t0;

    // This defines a damping factor for a mitigation strategy for different points in time.
    std::vector<Damping> dampings;

    /**
     * @brief Initializes a SEIR model with some default parameters
     */
    SecirParams();

    /**
     * @brief Initializes a SEIR model with given parameters
     */
    SecirParams(double tinc, double tinfmild, double base_reprod_in, double nb_total_t0_in, double nb_exp_t0_in,
                double nb_inf_t0_in, double nb_rec_t0_in);

    /**
     * @brief Initializes a SECIR model
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
                double beta_in, double delta_in, double rho_in, double theta_in, double nb_total_t0_in,
                double nb_exp_t0_in, double nb_car_t0_in, double nb_inf_t0_in, double nb_hosp_t0_in,
                double nb_icu_t0_in, double nb_rec_t0_in, double nb_dead_t0_in);

    /**
     * @brief Adds a damping to the current model
     * @param d The damping, which is a factor and day from which the mitigation acts
     */
    void add_damping(const Damping& d);
};

/**
 * prints given parameters
 * @param[in] params the SeirParams parameter object
 */
void printSecirParams(SecirParams const& params);

/**
 * Returns the damping factor
 *
 * @param[in] damping_array Array of dampings
 * @param[in] day Current day
 */
double getDampingFactor(std::vector<Damping> const& damping_array, double day);

/**
 * Computes the current time-derivative of S, E, I, and R in the SEIR model
 * @tparam T the datatype of the cases
 * @param[in] params SEIR Model parameters, created by seir_param
 * @param[in] y current  S, E, I, and R values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, I, and R
 */
void seir_getDerivatives(SecirParams const& params, std::vector<double> const& y, double t, std::vector<double>& dydt);

/**
 * Computes the current time-derivative of S, E, I, and R in the SEIR model
 * @param[in] params SEIR Model parameters, created by seir_param
 * @tparam T the datatype of the cases
 * @param[in] y current  S, E, I, and R values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, I, and R
 */
void secir_getDerivatives(SecirParams const& params, std::vector<double> const& y, double t, std::vector<double>& dydt);

/**
 * Computes the seir curve by integration
 * @param[in] seir_0 Initial S, E, I, and R values at t0
 * @param[in] t0 start time of simulation
 * @param[in] tmax end time of simulation
 * @param[in] dt initial time step
 * @param[in] params SEIR model parameters
 *
 * @returns Vector of times t
 */
std::vector<double> simulate(double t0, double tmax, double dt, SecirParams const& params,
                             std::vector<std::vector<double>>& seir);
#endif // SECIR_H
