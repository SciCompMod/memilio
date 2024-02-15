/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: David Kerkmann, Sascha Korf, Khoa Nguyen
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef EPI_ABM_INFECTION_H
#define EPI_ABM_INFECTION_H

#include "abm/time.h"
#include "abm/infection_state.h"
#include "abm/virus_variant.h"
#include "abm/parameters.h"
#include "abm/person.h"
#include "ad/ad.hpp"

#include <vector>
#include <cmath>

namespace mio
{
namespace abm
{

/**
 * @brief Models the ViralLoad for an Infection, modelled on a log_10 scale.
 * Based on https://www.science.org/doi/full/10.1126/science.abi5273
 * Examplary ViralLoad courses can be seen for example in Fig. 4 B.
 */
struct ViralLoad {
    TimePoint start_date; ///< Start date of the ViralLoad concentration in the Person.
    TimePoint end_date; ///< End date of the ViralLoad concentration in the Person.
    ScalarType peak; ///< Peak amplitude of the ViralLoad.
    ScalarType incline; ///< Incline of the ViralLoad during incline phase in log_10 scale per day (always positive).
    ScalarType decline; ///< Decline of the ViralLoad during decline phase in log_10 scale per day (always negative).
};

template <typename FP = double>
class Infection
{

public:
    /**
     * @brief Create an Infection for a single Person.
     * Draws a random infection course.
     * @param[inout] rng Person::RandomNumberGenerator for the Person.
     * @param[in] virus Virus type of the Infection.
     * @param[in] age AgeGroup to determine the ViralLoad course.
     * @param[in] params Parameters of the Model.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state [Default: InfectionState::Exposed] #InfectionState at time of initializing the Infection.
     * @param[in] latest_exposure [Default: {ExposureType::NoProtection, TimePoint(0)}] The pair value of last ExposureType (previous Infection/Vaccination) and TimePoint of that protection.
     * @param[in] detected [Default: false] If the Infection is detected.     
     */
    Infection(typename Person<FP>::RandomNumberGenerator& rng, VirusVariant virus, AgeGroup age,
              const Parameters<FP>& params, TimePoint init_date, InfectionState init_state = InfectionState::Exposed,
              std::pair<ExposureType, TimePoint> latest_exposure = {ExposureType::NoProtection, TimePoint(0)},
              bool detected                                      = false)
        : m_virus_variant(virus)
        , m_detected(detected)
    {
        assert(age.get() < params.get_num_groups());
        m_viral_load.start_date = draw_infection_course(rng, age, params, init_date, init_state, latest_exposure);

        auto vl_params                    = params.template get<ViralLoadDistributions>()[{virus, age}];
        ScalarType high_viral_load_factor = 1;
        if (latest_exposure.first != ExposureType::NoProtection) {
            high_viral_load_factor -=
                params.template get<HighViralLoadProtectionFactor>()(init_date.days() - latest_exposure.second.days());
        }
        m_viral_load.peak =
            vl_params.viral_load_peak.get_distribution_instance()(rng, vl_params.viral_load_peak.params) *
            high_viral_load_factor;
        m_viral_load.incline =
            vl_params.viral_load_incline.get_distribution_instance()(rng, vl_params.viral_load_incline.params);
        m_viral_load.decline =
            vl_params.viral_load_decline.get_distribution_instance()(rng, vl_params.viral_load_decline.params);
        m_viral_load.end_date =
            m_viral_load.start_date +
            days(int(m_viral_load.peak / m_viral_load.incline - m_viral_load.peak / m_viral_load.decline));

        auto inf_params = params.template get<InfectivityDistributions>()[{virus, age}];
        m_log_norm_alpha =
            inf_params.infectivity_alpha.get_distribution_instance()(rng, inf_params.infectivity_alpha.params);
        m_log_norm_beta =
            inf_params.infectivity_beta.get_distribution_instance()(rng, inf_params.infectivity_beta.params);
    }

    /**
     * @brief Gets the ViralLoad of the Infection at a given TimePoint.
     * @param[in] t TimePoint of querry.
     */
    ScalarType get_viral_load(TimePoint t) const
    {
        if (t >= m_viral_load.start_date && t <= m_viral_load.end_date) {
            if (t.days() <= m_viral_load.start_date.days() + m_viral_load.peak / m_viral_load.incline) {
                return m_viral_load.incline * (t - m_viral_load.start_date).days();
            }
            else {
                return m_viral_load.peak + m_viral_load.decline * (t.days() - m_viral_load.peak / m_viral_load.incline -
                                                                   m_viral_load.start_date.days());
            }
        }
        else {
            return 0.;
        }
    }

    /**
     * @brief Get infectivity at a given time.
     * Computed depending on current ViralLoad and individual invlogit function of each Person
     * corresponding to https://www.science.org/doi/full/10.1126/science.abi5273
     * The mapping corresponds to Fig. 2 C.
     * Formula of invlogit function can be found here:
     * https://github.com/VirologyCharite/SARS-CoV-2-VL-paper/tree/main
     * in ExtendedMethods.html, Section 3.1.2.1.
     * @param[in] t TimePoint of the querry.
     * @return Infectivity at given TimePoint.
     */
    ScalarType get_infectivity(TimePoint t) const
    {
        if (m_viral_load.start_date >= t || get_infection_state(t) == InfectionState::Exposed)
            return 0;
        return 1 / (1 + exp(-(m_log_norm_alpha + m_log_norm_beta * get_viral_load(t))));
    }

    /**
     * @brief: Get VirusVariant.
     * @return VirusVariant of the Infection.
     */
    VirusVariant get_virus_variant() const
    {
        return m_virus_variant;
    }

    /**
     * @brief Get the #InfectionState of the Infection.
     * @param[in] t TimePoint of the querry.
     * @return #InfectionState at the given TimePoint.
     */
    InfectionState get_infection_state(TimePoint t) const
    {
        if (t < m_infection_course[0].first)
            return InfectionState::Susceptible;

        return (*std::prev(std::upper_bound(m_infection_course.begin(), m_infection_course.end(), t,
                                            [](const TimePoint& s, std::pair<TimePoint, InfectionState> state) {
                                                return state.first > s;
                                            })))
            .second;
    }

    /**
     * @brief Set the Infection to detected.
     */
    void set_detected()
    {
        m_detected = true;
    }

    /**
     * @return Get the detected state.
     */
    bool is_detected() const
    {
        return m_detected;
    }

    /**
     * @returns Get the start date of the infection.
    */
    TimePoint get_start_date() const
    {
        return m_viral_load.start_date;
    }

private:
    /**
     * @brief Determine ViralLoad course and Infection course based on init_state.
     * Calls draw_infection_course_backward for all #InfectionState%s prior and draw_infection_course_forward for all
     * subsequent #InfectionState%s.
     * @param[inout] rng Person::RandomNumberGenerator of the Person.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state #InfectionState at time of initializing the Infection.
     * @return The starting date of the Infection.
     */
    TimePoint draw_infection_course(typename Person<FP>::RandomNumberGenerator& rng, AgeGroup age,
                                    const Parameters<FP>& params, TimePoint init_date, InfectionState init_state,
                                    std::pair<ExposureType, TimePoint> latest_protection)
    {
        assert(age.get() < params.get_num_groups());
        TimePoint start_date = draw_infection_course_backward(rng, age, params, init_date, init_state);
        draw_infection_course_forward(rng, age, params, init_date, init_state, latest_protection);
        return start_date;
    }

    /**
     * @brief Determine ViralLoad course and Infection course prior to the given start_state.
     * @param[inout] rng Person::RandomNumberGenerator of the Person.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state #InfectionState at time of initializing the Infection.
     */
    void draw_infection_course_forward(typename Person<FP>::RandomNumberGenerator& rng, AgeGroup age,
                                       const Parameters<FP>& params, TimePoint init_date, InfectionState start_state,
                                       std::pair<ExposureType, TimePoint> latest_exposure)
    {
        assert(age.get() < params.get_num_groups());
        auto t = init_date;
        TimeSpan time_period{}; // time period for current infection state
        InfectionState next_state{start_state}; // next state to enter
        m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, next_state));
        auto& uniform_dist = UniformDistribution<double>::get_instance();
        ScalarType v; // random draws
        while ((next_state != InfectionState::Recovered && next_state != InfectionState::Dead)) {
            switch (next_state) {
            case InfectionState::Exposed:
                // roll out how long until infected without symptoms
                time_period = days(ad::value(
                    FP(params.template get<IncubationPeriod<FP>>()[{m_virus_variant, age}]))); // subject to change
                next_state  = InfectionState::InfectedNoSymptoms;
                break;
            case InfectionState::InfectedNoSymptoms:
                // roll out next infection step
                v = uniform_dist(rng);
                if (v < 0.5) { // TODO: subject to change
                    time_period = days(ad::value(FP(params.template get<InfectedNoSymptomsToSymptoms<FP>>()[{
                        m_virus_variant, age}]))); // TODO: subject to change
                    next_state  = InfectionState::InfectedSymptoms;
                }
                else {
                    time_period = days(ad::value(FP(params.template get<InfectedNoSymptomsToRecovered<FP>>()[{
                        m_virus_variant, age}]))); // TODO: subject to change
                    next_state  = InfectionState::Recovered;
                }

                break;
            case InfectionState::InfectedSymptoms:
                // roll out next infection step
                {
                    ScalarType severity_protection_factor = 0.5;
                    v                                     = uniform_dist(rng);
                    if (latest_exposure.first != ExposureType::NoProtection) {
                        severity_protection_factor = params.template get<SeverityProtectionFactor>()[{
                            latest_exposure.first, age, m_virus_variant}](t.days() - latest_exposure.second.days());
                    }
                    if (v < (1 - severity_protection_factor) * 0.5) {
                        time_period = days(ad::value(FP(params.template get<InfectedSymptomsToSevere<FP>>()[{
                            m_virus_variant, age}]))); // TODO: subject to change
                        next_state  = InfectionState::InfectedSevere;
                    }
                    else {
                        time_period = days(ad::value(FP(params.template get<InfectedSymptomsToRecovered<FP>>()[{
                            m_virus_variant, age}]))); // TODO: subject to change
                        next_state  = InfectionState::Recovered;
                    }
                    break;
                }
            case InfectionState::InfectedSevere:
                // roll out next infection step
                v = uniform_dist(rng);
                if (v < 0.5) { // TODO: subject to change
                    time_period = days(ad::value(FP(params.template get<SevereToCritical<FP>>()[{
                        m_virus_variant, age}]))); // TODO: subject to change
                    next_state  = InfectionState::InfectedCritical;
                }
                else {
                    time_period = days(ad::value(FP(params.template get<SevereToRecovered<FP>>()[{
                        m_virus_variant, age}]))); // TODO: subject to change
                    next_state  = InfectionState::Recovered;
                }
                break;
            case InfectionState::InfectedCritical:
                // roll out next infection step
                v = uniform_dist(rng);
                if (v < 0.5) { // TODO: subject to change
                    time_period = days(ad::value(FP(
                        params.template get<CriticalToDead<FP>>()[{m_virus_variant, age}]))); // TODO: subject to change
                    next_state  = InfectionState::Dead;
                }
                else {
                    time_period = days(ad::value(FP(params.template get<CriticalToRecovered<FP>>()[{
                        m_virus_variant, age}]))); // TODO: subject to change
                    next_state  = InfectionState::Recovered;
                }
                break;
            default:
                break;
            }
            t = t + time_period;
            m_infection_course.push_back({t, next_state});
        }
    }

    /**
     * @brief Determine ViralLoad course and Infection course subsequent to the given start_state.
     * @param[inout] rng Person::RandomNumberGenerator of the Person.
     * @param[in] age AgeGroup of the person.
     * @param[in] params Parameters of the Model.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state InfectionState at time of initializing the Infection.
     * @return The starting date of the Infection.
     */
    TimePoint draw_infection_course_backward(typename Person<FP>::RandomNumberGenerator& rng, AgeGroup age,
                                             const Parameters<FP>& params, TimePoint init_date,
                                             InfectionState init_state)
    {
        assert(age.get() < params.get_num_groups());
        auto start_date = init_date;
        TimeSpan time_period{}; // time period for current infection state
        InfectionState previous_state{init_state}; // next state to enter
        auto& uniform_dist = UniformDistribution<double>::get_instance();
        ScalarType v; // random draws

        while ((previous_state != InfectionState::Exposed)) {
            switch (previous_state) {

            case InfectionState::InfectedNoSymptoms:
                time_period    = days(ad::value(FP(
                    params.template get<IncubationPeriod<FP>>()[{m_virus_variant, age}]))); // TODO: subject to change
                previous_state = InfectionState::Exposed;
                break;

            case InfectionState::InfectedSymptoms:
                time_period    = days(ad::value(FP(params.template get<InfectedNoSymptomsToSymptoms<FP>>()[{
                    m_virus_variant, age}]))); // TODO: subject to change
                previous_state = InfectionState::InfectedNoSymptoms;
                break;

            case InfectionState::InfectedSevere:
                time_period    = days(ad::value(FP(params.template get<InfectedSymptomsToSevere<FP>>()[{
                    m_virus_variant, age}]))); // TODO: subject to change
                previous_state = InfectionState::InfectedSymptoms;
                break;

            case InfectionState::InfectedCritical:
                time_period    = days(ad::value(FP(
                    params.template get<SevereToCritical<FP>>()[{m_virus_variant, age}]))); // TODO: subject to change
                previous_state = InfectionState::InfectedSevere;
                break;

            case InfectionState::Recovered:
                // roll out next infection step
                v = uniform_dist(rng);
                if (v < 0.25) {
                    time_period    = days(ad::value(FP(params.template get<InfectedNoSymptomsToRecovered<FP>>()[{
                        m_virus_variant, age}]))); // TODO: subject to change
                    previous_state = InfectionState::InfectedNoSymptoms;
                }
                else if (v < 0.5) { // TODO: subject to change
                    time_period    = days(ad::value(FP(params.template get<InfectedSymptomsToRecovered<FP>>()[{
                        m_virus_variant, age}]))); // TODO: subject to change
                    previous_state = InfectionState::InfectedSymptoms;
                }
                else if (v < 0.75) {
                    time_period    = days(ad::value(FP(params.template get<SevereToRecovered<FP>>()[{
                        m_virus_variant, age}]))); // TODO: subject to change
                    previous_state = InfectionState::InfectedSevere;
                }
                else {
                    time_period    = days(ad::value(FP(params.template get<CriticalToRecovered<FP>>()[{
                        m_virus_variant, age}]))); // TODO: subject to change
                    previous_state = InfectionState::InfectedCritical;
                }
                break;

            case InfectionState::Dead:
                time_period    = days(ad::value(
                    FP(params.template get<CriticalToDead<FP>>()[{m_virus_variant, age}]))); // TODO: subject to change
                previous_state = InfectionState::InfectedCritical;
                break;

            default:
                break;
            }
            start_date = start_date - time_period;
            m_infection_course.insert(m_infection_course.begin(), {start_date, previous_state});
        }
        return start_date;
    }

    std::vector<std::pair<TimePoint, InfectionState>> m_infection_course; ///< Start date of each #InfectionState.
    VirusVariant m_virus_variant; ///< Variant of the Infection.
    ViralLoad m_viral_load; ///< ViralLoad of the Infection.
    ScalarType m_log_norm_alpha,
        m_log_norm_beta; ///< Parameters for the infectivity mapping, which is modelled through an invlogit function.
    bool m_detected; ///< Whether an Infection is detected or not.
};

} // namespace abm
} // namespace mio

#endif
