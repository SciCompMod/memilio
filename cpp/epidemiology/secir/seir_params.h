#ifndef SEIR_PARAMS_H
#define SEIR_PARAMS_H

#include "epidemiology/math/integrator.h"
#include "epidemiology/model/populations.h"
#include "epidemiology/secir/damping.h"

#include <Eigen/Core>
#include <vector>

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

    void check_constraints() const {}; //dummy check_constraints function

    StageTimes times;

    // This defines a damping factor for a mitigation strategy for different points in time.
    Dampings dampings;
};

} // namespace epi

#endif // SEIR_H
