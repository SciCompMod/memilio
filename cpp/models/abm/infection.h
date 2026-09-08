/*
* Copyright (C) 2020-2026 MEmilio
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
#ifndef MIO_ABM_INFECTION_H
#define MIO_ABM_INFECTION_H

#include "abm/personal_rng.h"
#include "memilio/io/default_serialize.h"
#include "abm/time.h"
#include "abm/symptom_state.h"
#include "abm/virus_variant.h"
#include "abm/parameters.h"

#include <vector>

namespace mio
{
namespace abm
{

/**
 * @brief Represents a transition period between two symptom states.
 */
struct StateTransition {
    SymptomState from_state;
    SymptomState to_state;
    TimeSpan duration; // Duration that the infection stays in from_state before transitioning to to_state.
};

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

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("ViralLoad")
            .add("start_date", start_date)
            .add("end_date", end_date)
            .add("peak", peak)
            .add("incline", incline)
            .add("decline", decline);
    }
};

/**
 * @brief Distributions of the relative time that people have been in their initial symptom state at the beginning of the simulation.
 * Values have to be within [0, 1].
 * This makes it possible to draw from a user-defined distribution instead of drawing from a uniform distribution.
 */
using InitialSymptomStateDistribution = CustomIndexArray<AbstractParameterDistribution, VirusVariant, AgeGroup>;

class Infection
{
public:
    /**
     * @brief Create an Infection for a single Person.
     * Draws a random infection course.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] virus Virus type of the Infection.
     * @param[in] age AgeGroup to determine the ViralLoad course.
     * @param[in] params Parameters of the Model.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state [Default: SymptomState::Exposed] #SymptomState at time of initializing the Infection.
    * @param[in] detected [Default: false] If the Infection is detected.  
    * @param[in] recovered [Default: false] If the Person is already recovered at the time of initializing the Infection. Only works if the Person has start_state = SymptomState::None.
    *  @param[in] latest_protection [Default: {ProtectionType::NoProtection, TimePoint(0)}] The pair value of last ProtectionType (previous Infection/Vaccination) and TimePoint of that protection.
       
     */
    Infection(PersonalRandomNumberGenerator& rng, VirusVariant virus, AgeGroup age, const Parameters& params,
              TimePoint start_date, SymptomState start_state = SymptomState::None, bool recovered = false,
              bool detected = false, ProtectionEvent latest_protection = {ProtectionType::NoProtection, TimePoint(0)});

    /**
     * @brief Create an Infection for a single Person with a time spent in the given initial state that is drawn from the given distribution.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] virus Virus type of the Infection.
     * @param[in] age AgeGroup to determine the ViralLoad course.
     * @param[in] params Parameters of the Model.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state #SymptomState at time of initializing the Infection.
     * @param[in] init_state_dist Distribution to draw the relative time spent in the initial state from. Values have to be within [0, 1].
     * @param[in] detected If the Infection is detected.
     * @param[in] latest_protection The pair value of last ProtectionType (previous Infection/Vaccination) and TimePoint of that protection.
     */
    Infection(PersonalRandomNumberGenerator& rng, VirusVariant virus, AgeGroup age, const Parameters& params,
              TimePoint init_date, SymptomState init_state, const InitialSymptomStateDistribution& init_state_dist,
              bool detected = false, ProtectionEvent latest_protection = {ProtectionType::NoProtection, TimePoint(0)});

    /**
     * @brief Gets the ViralLoad of the Infection at a given TimePoint.
     * @param[in] t TimePoint of querry.
     */
    ScalarType get_viral_load(TimePoint t) const;

    /**
     * @brief Get viral shed at a specific time.
     * Computed depending on current ViralLoad and individual invlogit function of each Person
     * corresponding to https://www.science.org/doi/full/10.1126/science.abi5273
     * The mapping corresponds to Fig. 2 C.
     * Formula of invlogit function can be found here:
     * https://github.com/VirologyCharite/SARS-CoV-2-VL-paper/tree/main
     * in ExtendedMethods.html, Section 3.1.2.1.
     * * Also in accordance to Fig. 3d of another publication:
     * https://www.nature.com/articles/s41564-022-01105-z/figures/3
     * The result is in arbitrary units and has to be scaled to the rate "infections per day".
     * @param[in] t TimePoint of the querry.
     * @return Viral shed at given TimePoint.
     */
    ScalarType get_viral_shed(TimePoint t) const;

    /**
     * @brief: Get VirusVariant.
     * @return VirusVariant of the Infection.
     */
    VirusVariant get_virus_variant() const;

    /**
     * @brief Get the #SymptomState of the Infection.
     * @param[in] t TimePoint of the querry.
     * @return #SymptomState at the given TimePoint.
     */
    SymptomState get_symptom_state(TimePoint t) const;

    /**
     * @brief Set the Infection to detected.
     */
    void set_detected();

    /**
     * @return Get the detected state.
     */
    bool is_detected() const;

    /**
     * @returns Get the start date of the infection.
    */
    TimePoint get_start_date() const;

    /**
     * @brief If the Person has recovered from the Infection at a given TimePoint.
     * @param[in] t TimePoint of the querry.
     * @return True if the Person has recovered, false otherwise.
     */
    bool is_recovered(TimePoint t) const;

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("Infection")
            .add("symptom_course", m_symptom_course)
            .add("virus_variant", m_virus_variant)
            .add("viral_load", m_viral_load)
            .add("log_norm_alpha", m_log_norm_alpha)
            .add("log_norm_beta", m_log_norm_beta)
            .add("individual_viral_shed_factor", m_individual_viral_shed_factor)
            .add("detected", m_detected);
    }

private:
    friend DefaultFactory<Infection>;
    Infection() = default;

    /**
     * @brief Determine Infection course subsequent to the given #SymptomState start_state.
     * From the start_state, a random path through the #SymptomState tree is chosen, that is
     * Susceptible -> InfectedNoSymptoms,
     * InfectedNoSymptoms -> InfectedSymptoms or InfectedNoSymptoms -> Recovered,
     * InfectedSymptoms -> Infected_Severe or InfectedSymptoms -> Recovered,
     * InfectedSevere -> InfectedCritical or InfectedSevere -> Recovered or InfectedSevere -> Dead,
     * InfectedCritical -> Recovered or InfectedCritical -> Dead,
     * until either Recoverd or Dead is reached.
     * The duration in each #SymptomState is taken from the respective parameter.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state #SymptomState at time of initializing the Infection.
     * @param[in] recovered If the Person is already recovered at the time of initializing the Infection. Only works if the Person has start_state = SymptomState::None.
     * @param[in] latest_protection Latest protection against Infection, has an influence on transition probabilities.
     */
    void draw_symptom_course_forward(PersonalRandomNumberGenerator& rng, AgeGroup age, const Parameters& params,
                                     TimePoint init_date, SymptomState init_state, bool recovered,
                                     ProtectionEvent latest_protection);

    /**
     * @brief Determine Infection course prior to the given #SymptomState start_state.
     * From the start_state, a random path through the #SymptomState tree is chosen backwards, until Susceptible is reached.
     * For more detailed information, refer to draw_symptom_course_forward.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] age AgeGroup of the person.
     * @param[in] params Parameters of the Model.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state #SymptomState at time of initializing the Infection.
     * @param[in] recovered If the Person is already recovered at the time of initializing the Infection. Only works if the Person has start_state = SymptomState::None.
     * @return The starting date of the Infection.
     */
    TimePoint draw_symptom_course_backward(PersonalRandomNumberGenerator& rng, AgeGroup age, const Parameters& params,
                                           TimePoint init_date, SymptomState init_state, bool recovered);

    /**
     * @brief Initialize the viral load parameters for the infection.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] virus Virus type of the Infection.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     * @param[in] latest_protection Latest protection against Infection.
     */
    void initialize_viral_load(PersonalRandomNumberGenerator& rng, VirusVariant virus, AgeGroup age,
                               const Parameters& params, ProtectionEvent latest_protection);

    /**
     * @brief Initialize the viral shed parameters and individual factor for the infection.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] virus Virus type of the Infection.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     */
    void initialize_viral_shed(PersonalRandomNumberGenerator& rng, VirusVariant virus, AgeGroup age,
                               const Parameters& params);

    /**
     * @brief Get the forward transition from a given symptom state.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     * @param[in] current_state Current symptom state.
     * @param[in] current_time Current time point.
     * @param[in] latest_protection Latest protection against Infection.
     * @return StateTransition representing the next transition.
     */
    StateTransition get_forward_transition(PersonalRandomNumberGenerator& rng, AgeGroup age, const Parameters& params,
                                           SymptomState current_state, TimePoint current_time,
                                           ProtectionEvent latest_protection) const;

    /**
     * @brief Get the backward transition from a given symptom state.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     * @param[in] current_state Current symptom state.
     * @return StateTransition representing the previous transition.
     */
    StateTransition get_backward_transition(PersonalRandomNumberGenerator& rng, AgeGroup age, const Parameters& params,
                                            SymptomState current_state) const;

    /**
     * @brief Get the backward transition from recovered state.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     * @return StateTransition representing the transition that led to recovery.
     */
    StateTransition get_recovered_backward_transition(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                                      const Parameters& params) const;

    /**
     * @brief Get the backward transition from dead state.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     * @return StateTransition representing the transition that led to death.
     */
    StateTransition get_dead_backward_transition(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                                 const Parameters& params) const;

    /**
     * @brief Calculate the overall death probability for the infection.
     * @param[in] age AgeGroup of the Person.
     * @param[in] params Parameters of the Model.
     * @return The probability of death for this infection.
     */
    ScalarType calculate_death_probability(AgeGroup age, const Parameters& params) const;

    /**
     * @brief Get the severity protection factor based on latest protection.
     * @param[in] params Parameters of the Model.
     * @param[in] latest_protection Latest protection against Infection.
     * @param[in] age AgeGroup of the Person.
     * @param[in] current_time Current time point.
     * @return The protection factor against severe outcomes.
     */
    ScalarType get_severity_protection_factor(const Parameters& params, ProtectionEvent latest_protection, AgeGroup age,
                                              TimePoint current_time) const;

    std::vector<std::pair<TimePoint, SymptomState>> m_symptom_course; ///< Start date of each #SymptomState.
    VirusVariant m_virus_variant; ///< Variant of the Infection.
    ViralLoad m_viral_load; ///< ViralLoad of the Infection.
    ScalarType m_log_norm_alpha,
        m_log_norm_beta; ///< Parameters for the viral shed mapping, which is modelled through an invlogit function.
    ScalarType m_individual_viral_shed_factor; ///< Individual viral shed factor.
    bool m_detected; ///< Whether an Infection is detected or not.
};

} // namespace abm
} // namespace mio

#endif
