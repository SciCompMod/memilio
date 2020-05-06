#ifndef SEIR_H
#define SEIR_H

#include <epidemiology/damping.h>

#include <vector>

#include "Eigen/Core"

namespace epi
{

/**
 * Paramters of the SEIR model:
 * T_inc (also sigma^(-1) or, in the SECIR modile, R_2^(-1)+R_3^(-1)): mean incubation period (default: 5.2);
 *          in SECIR: R_2^(-1) is the first part of the incubation time where the person is not yet infectioous
 *          in SECIR: R_3 is the exchange between asymptomatic carriers and infectious people; R_3^(-1) is the second part of the incubation time where the person is infectious WITHOUT showing symptoms
 * T_infmild (also gamma^(-1) or R_4^(-1)): time a person remains infective after disease (if 'hospitalized' is considered a state, it does not apply to them but only to 'mildly infected' people in SECIR)
 * cont_freq (contact frequency/rate; called beta in the standard SEIR model, also R_1 in the SECIR model)
 *  NOTE: Here, the contact frequency is not calculated as R_0 * tinf_inv but directly input. This means that the Rt at the beginning (possibly, the R0) is given by Rt=tinf*cont_freq 
**/
class SeirParams
{
public:
    double base_reprod;
    double cont_freq, tinc_inv, tinfmild_inv;

    // double nb_total, nb_exp, nb_car, nb_inf, nb_hosp, nb_icu, nb_rec, nb_dead;
    double nb_total_t0, nb_sus_t0, nb_exp_t0, nb_inf_t0, nb_rec_t0;

    // This defines a damping factor for a mitigation strategy for different points in time.
    Dampings dampings;

    /**
     * @brief Initializes a SEIR model with some default parameters
     */
    SeirParams();

    /**
     * @brief Initializes a SEIR model with given parameters
     *
     * @todo parameter description
     *
     * @param tinc
     * @param tinfmild
     * @param cont_freq_in
     * @param nb_total_t0_in
     * @param nb_exp_t0_in
     * @param nb_inf_t0_in
     * @param nb_rec_t0_in
     */
    SeirParams(double tinc, double tinfmild, double cont_freq_in, double nb_total_t0_in, double nb_exp_t0_in,
               double nb_inf_t0_in, double nb_rec_t0_in);
};

/**
 * prints given parameters
 * @param[in] params the SeirParams parameter object
 */
void print_seir_params(SeirParams const& params);

/**
 * Computes the current time-derivative of S, E, I, and R in the SEIR model
 * @param[in] params SEIR Model parameters, created by seir_param
 * @tparam T the datatype of the cases
 * @param[in] y current  S, E, I, and R values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, I, and R
 */
void seir_get_derivatives(SeirParams const& params, const Eigen::VectorXd& y, double t, Eigen::VectorXd& dydt);

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
std::vector<double> simulate(double t0, double tmax, double dt, SeirParams const& params,
                             std::vector<Eigen::VectorXd>& seir);

} // namespace epi

#endif // SEIR_H
