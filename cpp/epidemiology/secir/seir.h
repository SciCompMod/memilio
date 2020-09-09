#ifndef SEIR_H
#define SEIR_H

#include "epidemiology/math/integrator.h"
#include "epidemiology/secir/populations.h"
#include "epidemiology/secir/damping.h"

#include <Eigen/Core>
#include <vector>

namespace epi
{

enum SeirCompartments
{
    S,
    E,
    I,
    R,
    SeirCount
};

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
    // time parameters for the different 'stages' of the disease of scale day or 1/day
    // 'stages' does not refer to the 'states' of the SEIR model but also includes incubation time or contact frequency
    class StageTimes
    {
    public:
        /**
         * @brief Initializes a time parameters' struct of the SEIR model
         */
        StageTimes();

        /**
         * @brief sets the contact frequency for the SEIR model
         * @param cont_freq contact rate/frequency in 1/day unit
         */
        void set_cont_freq(double const& cont_freq);

        /**
         * @brief sets the incubation time for the SEIR model
         * @param tinc incubation time in day unit
         */
        void set_incubation(double const& tinc);

        /**
         * @brief sets the infectious time for the SEIR model
         * @param tinfmild infectious time in day unit (in a generalized model, only for cases not treated in a hospital)
         */
        void set_infectious(double const& tinfmild);

        /**
         * @brief returns the contact frequency set for the SEIR model in 1/day unit
         */
        double get_cont_freq() const;

        /**
         * @brief returns 1.0 over the incubation time set for the SEIR model in day unit
         */
        double get_incubation_inv() const;

        /**
         * @brief returns 1.0 over the infectious time set for the SEIR model in day unit
         */
        double get_infectious_inv() const;

    private:
        double m_cont_freq, m_tinc_inv, m_tinfmild_inv;
    };

    StageTimes times;

    Populations populations{Populations({SeirCount})};

    // This defines a damping factor for a mitigation strategy for different points in time.
    Dampings dampings;
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
 * @brief simulate SEIR compartment model
 */
class SeirSimulation
{
public:
    SeirSimulation(const SeirParams& params, double t0 = 0., double dt_init = 0.1);
    Eigen::Ref<Eigen::VectorXd> advance(double tmax);
    TimeSeries<double>& get_result()
    {
        return m_integrator.get_result();
    }
    const TimeSeries<double>& get_result() const
    {
        return m_integrator.get_result();
    }

private:
    OdeIntegrator m_integrator;
};

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
