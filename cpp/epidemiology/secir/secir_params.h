#ifndef SECIR_PARAMS_H
#define SECIR_PARAMS_H

#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/uncertain_value.h"
#include "epidemiology/math/adapt_rk.h"
#include "epidemiology/secir/uncertain_matrix.h"

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
template <int N>
class SecirParams
{
public:
    SecirParams()
        : times(N, StageTimes())
        , probabilities(N, Probabilities())
        , m_num_groups{N}
        , m_contact_patterns(ContactFrequencyMatrix{N})
        , m_tstart{0}
        , m_seasonality{0}
        , m_icu_capacity{std::numeric_limits<double>::max()}
    {
    }

    SecirParams(ContactFrequencyMatrix cont_freq_matrix)
        : times(cont_freq_matrix.get_size(), StageTimes())
        , probabilities(cont_freq_matrix.get_size(), Probabilities())
        , m_num_groups{(size_t)cont_freq_matrix.get_size()}
        , m_contact_patterns(cont_freq_matrix)
        , m_tstart{0}
        , m_seasonality{0}
        , m_icu_capacity{std::numeric_limits<double>::max()}
    {
    }

    size_t get_num_groups() const
    {
        return m_num_groups;
    }

    double base_reprod;

    /**
     * @brief sets the start day in the SECIR model
     * The start day defines in which season the simulation can be started
     * If the start day is 180 and simulation takes place from t0=0 to
     * tmax=100 the days 180 to 280 of the year are simulated
     * @param tstart start day
     */
    void set_start_day(double tstart)
    {
        m_tstart = tstart;
    }

    /**
     * @brief returns the start day in the SECIR model
     * The start day defines in which season the simulation can be started
     * If the start day is 180 and simulation takes place from t0=0 to
     * tmax=100 the days 180 to 280 of the year are simulated
     */
    double get_start_day() const
    {
        return m_tstart;
    }

    /**
     * @brief sets the seasonality in the SECIR model
     * the seasonality is given as (1+k*sin()) where the sine
     * curve is below one in summer and above one in winter
     *
     * @param seasonality seasonality
     */
    void set_seasonality(UncertainValue const& seasonality)
    {
        m_seasonality = seasonality;
    }

    /**
     * @brief sets the seasonality in the SECIR model
     * the seasonality is given as (1+k*sin()) where the sine
     * curve is below one in summer and above one in winter
     *
     * @param seasonality seasonality
     */
    void set_seasonality(double seasonality)
    {
        m_seasonality = seasonality;
    }

    /**
     * @brief sets the seasonality in the SECIR model
     * the seasonality is given as (1+k*sin()) where the sine
     * curve is below one in summer and above one in winter
     *
     * @param seasonality seasonality
     */
    void set_seasonality(ParameterDistribution const& seasonality)
    {
        m_seasonality.set_distribution(seasonality);
    }

    /**
     * @brief returns the seasonality in the SECIR model
     * the seasonality is given as (1+k*sin()) where the sine
     * curve is below one in summer and above one in winter
     *
     */
    const UncertainValue& get_seasonality() const
    {
        return m_seasonality;
    }

    /**
     * @brief returns the seasonality in the SECIR model
     * the seasonality is given as (1+k*sin()) where the sine
     * curve is below one in summer and above one in winter
     *
     */
    UncertainValue& get_seasonality()
    {
        return m_seasonality;
    }

    /**
     * @brief sets the icu capacity in the SECIR model
     * @param icu_capacity icu capacity
     */
    void set_icu_capacity(UncertainValue const& icu_capacity)
    {
        m_icu_capacity = icu_capacity;
    }

    /**
     * @brief sets the icu capacity in the SECIR model
     * @param icu_capacity icu capacity
     */
    void set_icu_capacity(double icu_capacity)
    {
        m_icu_capacity = icu_capacity;
    }

    /**
     * @brief sets the icu capacity in the SECIR model
     * @param icu_capacity icu capacity
     */
    void set_icu_capacity(ParameterDistribution const& icu_capacity)
    {
        m_icu_capacity.set_distribution(icu_capacity);
    }

    /**
     * @brief returns the icu capacity in the SECIR model
     */
    const UncertainValue& get_icu_capacity() const
    {
        return m_icu_capacity;
    }

    /**
     * @brief returns the icu capacity in the SECIR model
     */
    UncertainValue& get_icu_capacity()
    {
        return m_icu_capacity;
    }

    // time parameters for the different 'stages' of the disease of scale day or 1/day
    // 'stages' does not refer to the 'states' of the SECIR model but also includes incubation time or contact frequency
    class StageTimes
    {
    public:
        /**
         * @brief Standard constructor of a time parameters' class in the SECIR model
         */
        StageTimes()
            : m_tinc{1.0}
            , m_tinfmild{1.0}
            , m_tserint{1.0}
            , m_thosp2home{1.0}
            , m_thome2hosp{1.0}
            , m_thosp2icu{1.0}
            , m_ticu2home{1.0}
            , m_tinfasy{1.0}
            , m_ticu2death{1.0}
        {
        }
        StageTimes(StageTimes&&)      = default;
        StageTimes(const StageTimes&) = default;
        StageTimes& operator=(const StageTimes&) = default;

        /**
         * @brief sets the incubation time in the SECIR model
         * @param tinc incubation time in day unit
         */
        void set_incubation(UncertainValue const& tinc)
        {
            m_tinc = tinc;
        }

        /**
         * @brief sets the incubation time in the SECIR model
         * @param tinc incubation time in day unit
         */
        void set_incubation(double tinc)
        {
            m_tinc = tinc;
        }

        /**
         * @brief sets the incubation time in the SECIR model
         * @param tinc incubation time in day unit
         */
        void set_incubation(ParameterDistribution const& tinc)
        {
            m_tinc.set_distribution(tinc);
        }

        /**
         * @brief sets the infectious time for symptomatic cases that are infected but
         *        who do not need to be hsopitalized in the SECIR model
         * @param tinfmild infectious time for symptomatic cases (if not hospitalized) in day unit
         */
        void set_infectious_mild(UncertainValue const& tinfmild)
        {
            m_tinfmild = tinfmild;
        }

        /**
         * @brief sets the infectious time for symptomatic cases that are infected but 
         *        who do not need to be hsopitalized in the SECIR model
         * @param tinfmild infectious time for symptomatic cases (if not hospitalized) in day unit 
         */
        void set_infectious_mild(double tinfmild)
        {
            m_tinfmild = tinfmild;
        }

        /**
         * @brief sets the infectious time for symptomatic cases that are infected but 
         *        who do not need to be hsopitalized in the SECIR model
         * @param tinfmild infectious time for symptomatic cases (if not hospitalized) in day unit 
         */
        void set_infectious_mild(ParameterDistribution const& tinfmild)
        {
            m_tinfmild.set_distribution(tinfmild);
        }

        /**
         * @brief sets the serial interval in the SECIR model
         * @param tserint serial interval in day unit
         */
        void set_serialinterval(UncertainValue const& tserint)
        {
            m_tserint = tserint;
        }

        /**
         * @brief sets the serial interval in the SECIR model
         * @param tserint serial interval in day unit 
         */
        void set_serialinterval(double tserint)
        {
            m_tserint = tserint;
        }

        /**
         * @brief sets the serial interval in the SECIR model
         * @param tserint serial interval in day unit 
         */
        void set_serialinterval(ParameterDistribution const& tserint)
        {
            m_tserint.set_distribution(tserint);
        }

        /**
         * @brief sets the time people are 'simply' hospitalized before returning home in the SECIR model
         * @param thosp2home time people are 'simply' hospitalized before returning home in day unit
         */
        void set_hospitalized_to_home(UncertainValue const& thosp2home)
        {
            m_thosp2home = thosp2home;
        }

        /**
         * @brief sets the time people are 'simply' hospitalized before returning home in the SECIR model
         * @param thosp2home time people are 'simply' hospitalized before returning home in day unit 
         */
        void set_hospitalized_to_home(double thosp2home)
        {
            m_thosp2home = thosp2home;
        }

        /**
         * @brief sets the time people are 'simply' hospitalized before returning home in the SECIR model
         * @param thosp2home time people are 'simply' hospitalized before returning home in day unit 
         */
        void set_hospitalized_to_home(ParameterDistribution const& thosp2home)
        {
            m_thosp2home.set_distribution(thosp2home);
        }

        /**
         * @brief sets the time people are infectious at home before 'simply' hospitalized in the SECIR model
         * @param thome2hosp time people are infectious at home before 'simply' hospitalized in day unit
         */
        void set_home_to_hospitalized(UncertainValue const& thome2hosp)
        {
            m_thome2hosp = thome2hosp;
        }

        /**
         * @brief sets the time people are infectious at home before 'simply' hospitalized in the SECIR model
         * @param thome2hosp time people are infectious at home before 'simply' hospitalized in day unit 
         */
        void set_home_to_hospitalized(double thome2hosp)
        {
            m_thome2hosp = thome2hosp;
        }

        /**
         * @brief sets the time people are infectious at home before 'simply' hospitalized in the SECIR model
         * @param thome2hosp time people are infectious at home before 'simply' hospitalized in day unit 
         */
        void set_home_to_hospitalized(ParameterDistribution const& thome2hosp)
        {
            m_thome2hosp.set_distribution(thome2hosp);
        }

        /**
         * @brief sets the time people are 'simply' hospitalized before being treated by ICU in the SECIR model
         * @param thosp2icu time people are 'simply' hospitalized before being treated by ICU in day unit
         */
        void set_hospitalized_to_icu(UncertainValue const& thosp2icu)
        {
            m_thosp2icu = thosp2icu;
        }

        /**
         * @brief sets the time people are 'simply' hospitalized before being treated by ICU in the SECIR model
         * @param thosp2icu time people are 'simply' hospitalized before being treated by ICU in day unit 
         */
        void set_hospitalized_to_icu(double thosp2icu)
        {
            m_thosp2icu = thosp2icu;
        }

        /**
         * @brief sets the time people are 'simply' hospitalized before being treated by ICU in the SECIR model
         * @param thosp2icu time people are 'simply' hospitalized before being treated by ICU in day unit 
         */
        void set_hospitalized_to_icu(ParameterDistribution const& thosp2icu)
        {
            m_thosp2icu.set_distribution(thosp2icu);
        }

        /**
         * @brief sets the time people are treated by ICU before returning home in the SECIR model
         * @param ticu2home time people are treated by ICU before returning home in day unit
         */
        void set_icu_to_home(UncertainValue const& ticu2home)
        {
            m_ticu2home = ticu2home;
        }

        /**
         * @brief sets the time people are treated by ICU before returning home in the SECIR model
         * @param ticu2home time people are treated by ICU before returning home in day unit 
         */
        void set_icu_to_home(double ticu2home)
        {
            m_ticu2home = ticu2home;
        }

        /**
         * @brief sets the time people are treated by ICU before returning home in the SECIR model
         * @param ticu2home time people are treated by ICU before returning home in day unit 
         */
        void set_icu_to_home(ParameterDistribution const& ticu2home)
        {
            m_ticu2home.set_distribution(ticu2home);
        }

        /**
         * @brief sets the infectious time for asymptomatic cases in the SECIR model
         * @param tinfasy infectious time for asymptomatic cases in day unit
         */
        void set_infectious_asymp(UncertainValue const& tinfasy)
        {
            log_warning(
                "The parameter of the asymptomatic infectious period is meant to be defined by the other parameters of "
                "the system. Do you really want to set it?");
            m_tinfasy = tinfasy;
        }

        /**
         * @brief sets the infectious time for asymptomatic cases in the SECIR model
         * @param tinfasy infectious time for asymptomatic cases in day unit 
         */
        void set_infectious_asymp(double tinfasy)
        {
            log_warning(
                "The parameter of the asymptomatic infectious period is meant to be defined by the other parameters of "
                "the system. Do you really want to set it?");
            m_tinfasy = tinfasy;
        }

        /**
         * @brief sets the infectious time for asymptomatic cases in the SECIR model
         * @param tinfasy infectious time for asymptomatic cases in day unit 
         */
        void set_infectious_asymp(ParameterDistribution const& tinfasy)
        {
            m_tinfasy.set_distribution(tinfasy);
        }

        /**
         * @brief sets the time people are treated by ICU before dying in the SECIR model
         * @param ticu2death time people are treated by ICU before dying in day unit
         */
        void set_icu_to_death(UncertainValue const& ticu2death)
        {
            m_ticu2death = ticu2death;
        }

        /**
         * @brief sets the time people are treated by ICU before dying in the SECIR model
         * @param ticu2death time people are treated by ICU before dying in day unit 
         */
        void set_icu_to_death(double ticu2death)
        {
            m_ticu2death = ticu2death;
        }

        /**
         * @brief sets the time people are treated by ICU before dying in the SECIR model
         * @param ticu2death time people are treated by ICU before dying in day unit 
         */
        void set_icu_to_death(ParameterDistribution const& ticu2death)
        {
            m_ticu2death.set_distribution(ticu2death);
        }

        /**
         * @brief returns incubation time set for the SECIR model in day unit
         */
        const UncertainValue& get_incubation() const
        {
            return m_tinc;
        }
        UncertainValue& get_incubation()
        {
            return m_tinc;
        }

        /**
         * @brief returns infectious time set for the SECIR model in day unit
         */
        const UncertainValue& get_infectious_mild() const
        {
            return m_tinfmild;
        }
        UncertainValue& get_infectious_mild()
        {
            return m_tinfmild;
        }

        /**
         * @brief returns serial interval in the SECIR model
         */
        const UncertainValue& get_serialinterval() const
        {
            return m_tserint;
        }
        UncertainValue& get_serialinterval()
        {
            return m_tserint;
        }

        /**
         * @brief returns time people are 'simply' hospitalized before returning home in the SECIR model 
         */
        const UncertainValue& get_hospitalized_to_home() const
        {
            return m_thosp2home;
        }
        UncertainValue& get_hospitalized_to_home()
        {
            return m_thosp2home;
        }

        /**
         * @brief returns time people are infectious at home before 'simply' hospitalized in the SECIR model 
         */
        const UncertainValue& get_home_to_hospitalized() const
        {
            return m_thome2hosp;
        }
        UncertainValue& get_home_to_hospitalized()
        {
            return m_thome2hosp;
        }

        /**
         * @brief returns time people are 'simply' hospitalized before being treated by ICU in the SECIR model
         */
        const UncertainValue& get_hospitalized_to_icu() const
        {
            return m_thosp2icu;
        }
        UncertainValue& get_hospitalized_to_icu()
        {
            return m_thosp2icu;
        }

        /**
         * @brief returns time people are treated by ICU before returning home in the SECIR model
         */
        const UncertainValue& get_icu_to_home() const
        {
            return m_ticu2home;
        }
        UncertainValue& get_icu_to_home()
        {
            return m_ticu2home;
        }

        /**
         * @brief returns infectious time for asymptomatic cases in the SECIR model
         */
        const UncertainValue& get_infectious_asymp() const
        {
            return m_tinfasy;
        }
        UncertainValue& get_infectious_asymp()
        {
            return m_tinfasy;
        }

        /**
         * @brief returns time people are treated by ICU before dying in the SECIR model
         */
        const UncertainValue& get_icu_to_dead() const
        {
            return m_ticu2death;
        }
        UncertainValue& get_icu_to_dead()
        {
            return m_ticu2death;
        }

        /**
         * @brief checks whether the stage times Parameters satisfy their corresponding constraints and applies them, if they do not
         * For certain stage time values, the simulation can produce (small) negative population results.
         * To prevent this from happening, mathematical constraints have been derived and implemented here.
         * For further details, see the overleaf discussion chapter.
         * The constraints are step size-dependent, in the current implementation, we have assumed dt_max = 1.
         */
        void apply_constraints()
        {

            if (m_tinc < 2.0) {
                log_warning("Constraint check: Parameter m_tinc changed from {:.4f} to {:.4f}", m_tinc, 2.0);
                m_tinc = 2.0;
            }

            if (2 * m_tserint < m_tinc + 1.0) {
                log_warning("Constraint check: Parameter m_tserint changed from {:.4f} to {:.4f}", m_tserint,
                            0.5 * m_tinc + 0.5);
                m_tserint = 0.5 * m_tinc + 0.5;
            }
            else if (m_tserint > m_tinc - 0.5) {
                log_warning("Constraint check: Parameter m_tserint changed from {:.4f} to {:.4f}", m_tserint,
                            m_tinc - 0.5);
                m_tserint = m_tinc - 0.5;
            }

            if (m_tinfmild < 1.0) {
                log_warning("Constraint check: Parameter m_tinfmild changed from {:.4f} to {:.4f}", m_tinfmild, 1.0);
                m_tinfmild = 1.0;
            }

            if (m_thosp2home < 1.0) {
                log_warning("Constraint check: Parameter m_thosp2home changed from {:.4f} to {:.4f}", m_thosp2home,
                            1.0);
                m_thosp2home = 1.0;
            }

            if (m_thome2hosp < 1.0) {
                log_warning("Constraint check: Parameter m_thome2hosp changed from {:.4f} to {:.4f}", m_thome2hosp,
                            1.0);
                m_thome2hosp = 1.0;
            }

            if (m_thosp2icu < 1.0) {
                log_warning("Constraint check: Parameter m_thosp2icu changed from {:.4f} to {:.4f}", m_thosp2icu, 1.0);
                m_thosp2icu = 1.0;
            }

            if (m_ticu2home < 1.0) {
                log_warning("Constraint check: Parameter m_ticu2home changed from {:.4f} to {:.4f}", m_ticu2home, 1.0);
                m_ticu2home = 1.0;
            }

            if (m_tinfasy != 1.0 / (0.5 / (m_tinc - m_tserint) + 0.5 / m_tinfmild)) {
                log_info("Constraint check: Parameter m_tinfasy set as fully dependent on tinc, tserint and tinfmild. "
                         "See HZI "
                         "paper.");
                m_tinfasy = 1.0 / (0.5 / (m_tinc - m_tserint) + 0.5 / m_tinfmild);
            }

            if (m_ticu2death < 1.0) {
                log_warning("Constraint check: Parameter m_ticu2death changed from {:.4f} to {:.4f}", m_ticu2death,
                            1.0);
                m_ticu2death = 1.0;
            }
        }

        /**
         * @brief checks whether the stage times Parameters satisfy their corresponding constraints and throws errors, if they do not
         * For certain stage time values, the simulation can produce (small) negative population results.
         * To prevent this from happening, mathematical constraints have been derived and implemented here.
         * For further details, see the overleaf discussion chapter.
         * The constraints are step size-dependent, in the current implementation, we have assumed dt_max = 1.
         */
        void check_constraints() const
        {

            if (m_tinc < 2.0) {
                log_error("Constraint check: Parameter m_tinc {:.4f} smaller {:.4f}", m_tinc, 2.0);
            }

            if (2 * m_tserint < m_tinc + 1.0) {
                log_error("Constraint check: Parameter m_tserint {:.4f} smaller {:.4f}", m_tserint, 0.5 * m_tinc + 0.5);
            }
            else if (m_tserint > m_tinc - 0.5) {
                log_error("Constraint check: Parameter m_tserint {:.4f} smaller {:.4f}", m_tserint, m_tinc - 0.5);
            }

            if (m_tinfmild < 1.0) {
                log_error("Constraint check: Parameter m_tinfmild {:.4f} smaller {:.4f}", m_tinfmild, 1.0);
            }

            if (m_thosp2home < 1.0) {
                log_error("Constraint check: Parameter m_thosp2home {:.4f} smaller {:.4f}", m_thosp2home, 1.0);
            }

            if (m_thome2hosp < 1.0) {
                log_error("Constraint check: Parameter m_thome2hosp {:.4f} smaller {:.4f}", m_thome2hosp, 1.0);
            }

            if (m_thosp2icu < 1.0) {
                log_error("Constraint check: Parameter m_thosp2icu {:.4f} smaller {:.4f}", m_thosp2icu, 1.0);
            }

            if (m_ticu2home < 1.0) {
                log_error("Constraint check: Parameter m_ticu2home {:.4f} smaller {:.4f}", m_ticu2home, 1.0);
            }

            if (m_tinfasy != 1.0 / (0.5 / (m_tinc - m_tserint) + 0.5 / m_tinfmild)) {
                log_error("Constraint check: Parameter m_tinfasy not set as fully dependent on tinc, tserint and "
                          "tinfmild. See "
                          "HZI paper.");
            }

            if (m_ticu2death < 1.0) {
                log_error("Constraint check: Parameter m_ticu2death {:.4f} smaller {:.4f}", m_ticu2death, 1.0);
            }
        }

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
        Probabilities()
            : m_infprob{1}
            , m_carrinf{1}
            , m_asympinf{0}
            , m_risksymp{0}
            , m_hospinf{0}
            , m_icuhosp{0}
            , m_deathicu{0}
        {
        }
        Probabilities(const Probabilities&) = default;
        Probabilities(Probabilities&&)      = default;
        Probabilities& operator=(const Probabilities&) = default;

        /**
        * @brief sets probability of getting infected from a contact
        * @param infprob the probability of getting infected from a contact
        */
        void set_infection_from_contact(UncertainValue const& infprob)
        {
            m_infprob = infprob;
        }

        /**
        * @brief sets probability of getting infected from a contact
        * @param infprob the probability of getting infected from a contact
        */
        void set_infection_from_contact(double infprob)
        {
            m_infprob = infprob;
        }

        /**
        * @brief sets probability of getting infected from a contact
        * @param infprob the probability of getting infected from a contact
        */
        void set_infection_from_contact(ParameterDistribution const& infprob)
        {
            m_infprob.set_distribution(infprob);
        }

        /**
        * @brief sets the relative carrier infectability 
        * @param carrinf relative carrier infectability
        */
        void set_carrier_infectability(UncertainValue const& carrinf)
        {
            m_carrinf = carrinf;
        }

        /**
        * @brief sets the relative carrier infectability
        * @param carrinf relative carrier infectability
        */
        void set_carrier_infectability(double carrinf)
        {
            m_carrinf = carrinf;
        }

        /**
        * @brief sets the relative carrier infectability
        * @param carrinf relative carrier infectability
        */
        void set_carrier_infectability(ParameterDistribution const& carrinf)
        {
            m_carrinf.set_distribution(carrinf);
        }

        /**
        * @brief sets the percentage of asymptomatic cases in the SECIR model
        * @param alpha the percentage of asymptomatic cases
        */
        void set_asymp_per_infectious(UncertainValue const& asympinf)
        {
            m_asympinf = asympinf;
        }

        /**
        * @brief sets the percentage of asymptomatic cases in the SECIR model
        * @param alpha the percentage of asymptomatic cases
        */
        void set_asymp_per_infectious(double asympinf)
        {
            m_asympinf = asympinf;
        }

        /**
        * @brief sets the percentage of asymptomatic cases in the SECIR model
        * @param alpha the percentage of asymptomatic cases
        */
        void set_asymp_per_infectious(ParameterDistribution const& asympinf)
        {
            m_asympinf.set_distribution(asympinf);
        }

        /**
        * @brief sets the risk of infection from symptomatic cases in the SECIR model
        * @param beta the risk of infection from symptomatic cases
        */
        void set_risk_from_symptomatic(UncertainValue const& risksymp)
        {
            m_risksymp = risksymp;
        }

        /**
        * @brief sets the risk of infection from symptomatic cases in the SECIR model
        * @param beta the risk of infection from symptomatic cases 
        */
        void set_risk_from_symptomatic(double risksymp)
        {
            m_risksymp = risksymp;
        }

        /**
        * @brief sets the risk of infection from symptomatic cases in the SECIR model
        * @param beta the risk of infection from symptomatic cases 
        */
        void set_risk_from_symptomatic(ParameterDistribution const& risksymp)
        {
            m_risksymp.set_distribution(risksymp);
        }

        /**
        * @brief sets the percentage of hospitalized patients per infected patients in the SECIR model
        * @param rho percentage of hospitalized patients per infected patients
        */
        void set_hospitalized_per_infectious(UncertainValue const& hospinf)
        {
            m_hospinf = hospinf;
        }

        /**
        * @brief sets the percentage of hospitalized patients per infected patients in the SECIR model
        * @param rho percentage of hospitalized patients per infected patients
        */
        void set_hospitalized_per_infectious(double hospinf)
        {
            m_hospinf = hospinf;
        }

        /**
        * @brief sets the percentage of hospitalized patients per infected patients in the SECIR model
        * @param rho percentage of hospitalized patients per infected patients
        */
        void set_hospitalized_per_infectious(ParameterDistribution const& hospinf)
        {
            m_hospinf.set_distribution(hospinf);
        }

        /**
        * @brief sets the percentage of ICU patients per hospitalized patients in the SECIR model
        * @param theta percentage of ICU patients per hospitalized patients
        */
        void set_icu_per_hospitalized(UncertainValue const& icuhosp)
        {
            m_icuhosp = icuhosp;
        }

        /**
        * @brief sets the percentage of ICU patients per hospitalized patients in the SECIR model
        * @param theta percentage of ICU patients per hospitalized patients
        */
        void set_icu_per_hospitalized(double icuhosp)
        {
            m_icuhosp = icuhosp;
        }

        /**
        * @brief sets the percentage of ICU patients per hospitalized patients in the SECIR model
        * @param theta percentage of ICU patients per hospitalized patients
        */
        void set_icu_per_hospitalized(ParameterDistribution const& icuhosp)
        {
            m_icuhosp.set_distribution(icuhosp);
        }

        /**
        * @brief sets the percentage of dead patients per ICU patients in the SECIR model
        * @param delta percentage of dead patients per ICU patients
        */
        void set_dead_per_icu(UncertainValue const& deathicu)
        {
            m_deathicu = deathicu;
        }

        /**
        * @brief sets the percentage of dead patients per ICU patients in the SECIR model
        * @param delta percentage of dead patients per ICU patients 
        */
        void set_dead_per_icu(double deathicu)
        {
            m_deathicu = deathicu;
        }

        /**
        * @brief sets the percentage of dead patients per ICU patients in the SECIR model
        * @param delta percentage of dead patients per ICU patients 
        */
        void set_dead_per_icu(ParameterDistribution const& deathicu)
        {
            m_deathicu.set_distribution(deathicu);
        }

        /**
        * @brief gets probability of getting infected from a contact
        */
        const UncertainValue& get_infection_from_contact() const
        {
            return m_infprob;
        }
        UncertainValue& get_infection_from_contact()
        {
            return m_infprob;
        }

        /**
        * @brief set the relative carrier infectability
        */
        const UncertainValue& get_carrier_infectability() const
        {
            return m_carrinf;
        }
        UncertainValue& get_carrier_infectability()
        {
            return m_carrinf;
        }

        /**
        * @brief returns the percentage of asymptomatic cases in the SECIR model
        */
        const UncertainValue& get_asymp_per_infectious() const
        {
            return m_asympinf;
        }
        UncertainValue& get_asymp_per_infectious()
        {
            return m_asympinf;
        }

        /**
        * @brief returns the risk of infection from symptomatic cases in the SECIR model
        */
        const UncertainValue& get_risk_from_symptomatic() const
        {
            return m_risksymp;
        }
        UncertainValue& get_risk_from_symptomatic()
        {
            return m_risksymp;
        }

        /**
        * @brief returns the percentage of hospitalized patients per infected patients in the SECIR model
        */
        const UncertainValue& get_hospitalized_per_infectious() const
        {
            return m_hospinf;
        }
        UncertainValue& get_hospitalized_per_infectious()
        {
            return m_hospinf;
        }

        /**
        * @brief returns the percentage of ICU patients per hospitalized patients in the SECIR model
        */
        const UncertainValue& get_icu_per_hospitalized() const
        {
            return m_icuhosp;
        }
        UncertainValue& get_icu_per_hospitalized()
        {
            return m_icuhosp;
        }

        /**
        * @brief returns the percentage of dead patients per ICU patients in the SECIR model
        */
        const UncertainValue& get_dead_per_icu() const
        {
            return m_deathicu;
        }
        UncertainValue& get_dead_per_icu()
        {
            return m_deathicu;
        }

        /**
         * @brief checks whether the probability Parameters satisfy their corresponding constraints and applies them, if they do not
         */
        void apply_constraints()
        {
            if (m_infprob < 0.0) {
                log_warning("Constraint check: Parameter m_asympinf changed from {:0.4f} to {:d} ", m_infprob, 0);
                m_infprob = 0;
            }

            if (m_carrinf < 0.0) {
                log_warning("Constraint check: Parameter m_asympinf changed from {:0.4f} to {:d} ", m_carrinf, 0);
                m_carrinf = 0;
            }

            if (m_asympinf < 0.0 || m_asympinf > 1.0) {
                log_warning("Constraint check: Parameter m_asympinf changed from {:0.4f} to {:d} ", m_asympinf, 0);
                m_asympinf = 0;
            }

            if (m_risksymp < 0.0 || m_risksymp > 1.0) {
                log_warning("Constraint check: Parameter m_risksymp changed from {:0.4f} to {:d}", m_risksymp, 0);
                m_risksymp = 0;
            }

            if (m_hospinf < 0.0 || m_hospinf > 1.0) {
                log_warning("Constraint check: Parameter m_risksymp changed from {:0.4f} to {:d}", m_hospinf, 0);
                m_hospinf = 0;
            }

            if (m_icuhosp < 0.0 || m_icuhosp > 1.0) {
                log_warning("Constraint check: Parameter m_risksymp changed from {:0.4f} to {:d}", m_icuhosp, 0);
                m_icuhosp = 0;
            }

            if (m_deathicu < 0.0 || m_deathicu > 1.0) {
                log_warning("Constraint check: Parameter m_risksymp changed from {:0.4f} to {:d}", m_deathicu, 0);
                m_deathicu = 0;
            }
        }

        /**
         * @brief checks whether the probability Parameters satisfy their corresponding constraints and throws errors, if they do not
         */
        void check_constraints() const
        {
            if (m_infprob < 0.0) {
                log_warning("Constraint check: Parameter m_infprob smaller {:d}", 0);
            }

            if (m_carrinf < 0.0) {
                log_warning("Constraint check: Parameter m_carrinf smaller {:d}", 0);
            }

            if (m_asympinf < 0.0 || m_asympinf > 1.0) {
                log_warning("Constraint check: Parameter m_asympinf smaller {:d} or larger {:d}", 0, 1);
            }

            if (m_risksymp < 0.0 || m_risksymp > 1.0) {
                log_warning("Constraint check: Parameter m_risksymp smaller {:d} or larger {:d}", 0, 1);
            }

            if (m_hospinf < 0.0 || m_hospinf > 1.0) {
                log_warning("Constraint check: Parameter m_risksymp smaller {:d} or larger {:d}", 0, 1);
            }

            if (m_icuhosp < 0.0 || m_icuhosp > 1.0) {
                log_warning("Constraint check: Parameter m_risksymp smaller {:d} or larger {:d}", 0, 1);
            }

            if (m_deathicu < 0.0 || m_deathicu > 1.0) {
                log_warning("Constraint check: Parameter m_risksymp smaller {:d} or larger {:d}", 0, 1);
            }
        }

    private:
        UncertainValue m_infprob, m_carrinf, m_asympinf, m_risksymp, m_hospinf, m_icuhosp, m_deathicu; // probabilities
    };

    /**
     * @brief checks whether all Parameters satisfy their corresponding constraints and applies them, if they do not
     */
    void apply_constraints()
    {
        if (m_seasonality < 0.0 || m_seasonality > 0.5) {
            log_warning("Constraint check: Parameter m_seasonality changed from {:0.4f} to {:d}", m_seasonality, 0);
            m_seasonality = 0;
        }

        if (m_icu_capacity < 0.0) {
            log_warning("Constraint check: Parameter m_icu_capacity changed from {:0.4f} to {:d}", m_icu_capacity, 0);
            m_icu_capacity = 0;
        }

        for (size_t i = 0; i < times.size(); i++) {
            times[i].apply_constraints();
            probabilities[i].apply_constraints();
        }
    }

    /**
     * @brief checks whether all Parameters satisfy their corresponding constraints and throws errors, if they do not
     */
    void check_constraints() const
    {
        if (m_seasonality < 0.0 || m_seasonality > 0.5) {
            log_warning("Constraint check: Parameter m_seasonality smaller {:d} or larger {:d}", 0, 0.5);
        }

        if (m_icu_capacity < 0.0) {
            log_warning("Constraint check: Parameter m_icu_capacity smaller {:d}", 0);
        }

        for (size_t i = 0; i < times.size(); i++) {
            times[i].check_constraints();
            probabilities[i].check_constraints();
        }
    }
    /**
     * @brief sets the UncertainContactMatrix
     */
    void set_contact_patterns(UncertainContactMatrix contact_patterns)
    {
        m_contact_patterns = contact_patterns;
    }

    /**
     * @brief returns the UncertainContactMatrix
     */
    UncertainContactMatrix& get_contact_patterns()
    {
        return m_contact_patterns;
    }

    /**
     * @brief returns the UncertainContactMatrix
     */
    UncertainContactMatrix const& get_contact_patterns() const
    {
        return m_contact_patterns;
    }

    std::vector<StageTimes> times;
    std::vector<Probabilities> probabilities;

private:
    size_t m_num_groups;

    UncertainContactMatrix m_contact_patterns;

    double m_tstart;
    UncertainValue m_seasonality;
    UncertainValue m_icu_capacity;
};

/**
 * @brief WIP !! TO DO: returns the actual, approximated reproduction rate 
 */
//double get_reprod_rate(SecirParams const& params, double t, std::vector<double> const& yt);

} // namespace epi

#endif // SECIR_H
