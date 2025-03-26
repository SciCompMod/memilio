/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kühn
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
#ifndef MIO_ODE_SECIRVVS_PARAMETERS_H
#define MIO_ODE_SECIRVVS_PARAMETERS_H

#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/dynamic_npis.h"
#include "memilio/epidemiology/simulation_day.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/uncertain_value.h"

namespace mio
{
namespace osecirvvs
{

/**
* @brief the start day in the SECIRVVS model
* The start day defines in which season the simulation can be started
* If the start day is 180 and simulation takes place from t0=0 to
* tmax=100 the days 180 to 280 of the year are simulated
*/
struct StartDay {
    using Type = double;
    static Type get_default(AgeGroup)
    {
        return 0.;
    }
    static std::string name()
    {
        return "StartDay";
    }
};

/**
* @brief the start day of a new variant in the SECIRVVS model
* The start day of the new variant defines in which day of the simulation the new variant is introduced.
* Starting on this day, the new variant will impact the transmission probability depending on the
* infectiousness of the new variant in the parameter InfectiousnessNewVariant.
*/
struct StartDayNewVariant {
    using Type = double;
    static Type get_default(AgeGroup)
    {
        return std::numeric_limits<double>::max();
    }
    static std::string name()
    {
        return "StartDayNewVariant";
    }
};

/**
* @brief the seasonality in the SECIR model
* the seasonality is given as (1+k*sin()) where the sine
* curve is below one in summer and above one in winter
*/
template <typename FP = double>
struct Seasonality {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(0.);
    }
    static std::string name()
    {
        return "Seasonality";
    }
};

/**
* @brief the icu capacity in the SECIR model
*/
template <typename FP = double>
struct ICUCapacity {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(std::numeric_limits<FP>::max());
    }
    static std::string name()
    {
        return "ICUCapacity";
    }
};

/**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
template <typename FP = double>
struct TestAndTraceCapacity {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(std::numeric_limits<FP>::max());
    }
    static std::string name()
    {
        return "TestAndTraceCapacity";
    }
};

/**
 * @brief Multiplier for the test and trace capacity to determine when it is considered overloaded from cases without symptoms.
 */
template <typename FP = double>
struct TestAndTraceCapacityMaxRiskNoSymptoms {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(2.0);
    }
    static std::string name()
    {
        return "TestAndTraceCapacityMaxRiskNoSymptoms";
    }
};

/**
 * @brief Multiplier for the test and trace capacity to determine when it is considered overloaded by symptomatic cases.
 */
template <typename FP = double>
struct TestAndTraceCapacityMaxRiskSymptoms {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(15.0);
    }
    static std::string name()
    {
        return "TestAndTraceCapacityMaxRiskSymptoms";
    }
};

/**
 * @brief the contact patterns within the society are modelled using an UncertainContactMatrix
 */
template <typename FP = double>
struct ContactPatterns {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

/**
 * @brief the NPIs that are enforced if certain infection thresholds are exceeded.
 */
template <typename FP = double>
struct DynamicNPIsInfectedSymptoms {
    using Type = DynamicNPIs<FP>;
    static Type get_default(AgeGroup /*size*/)
    {
        return {};
    }
    static std::string name()
    {
        return "DynamicNPIsInfectedSymptoms";
    }
};

/**
 * @brief The delay with which DynamicNPIs are implemented and enforced after exceedance of threshold.
 */
template <typename FP = double>
struct DynamicNPIsImplementationDelay {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup /*size*/)
    {
        return 0.;
    }
    static std::string name()
    {
        return "DynamicNPIsImplementationDelay";
    }
};

/**
 * @brief the (mean) latent time in day unit
 */
template <typename FP = double>
struct TimeExposed {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "TimeExposed";
    }
};

/**
 * @brief the (mean) time in day unit for asymptomatic cases that are infected but
 *        have not yet developed symptoms.
 */
template <typename FP = double>
struct TimeInfectedNoSymptoms {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptoms";
    }
};

/**
* @brief the infectious time for symptomatic cases that are infected but
*        who do not need to be hospitalized in the SECIR model in day unit
*/
template <typename FP = double>
struct TimeInfectedSymptoms {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedSymptoms";
    }
};

/**
 * @brief the time people are 'simply' hospitalized before returning home in the SECIR model
 *        in day unit
 */
template <typename FP = double>
struct TimeInfectedSevere {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedSevere";
    }
};

/**
 * @brief the time people are treated by ICU before returning home in the SECIR model
 *        in day unit
 */
template <typename FP = double>
struct TimeInfectedCritical {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedCritical";
    }
};

/**
* @brief probability of getting infected from a contact
*/
template <typename FP = double>
struct TransmissionProbabilityOnContact {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
* @brief the relative InfectedNoSymptoms infectability
*/
template <typename FP = double>
struct RelativeTransmissionNoSymptoms {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RelativeTransmissionNoSymptoms";
    }
};

/**
* @brief the percentage of asymptomatic cases in the SECIR model
*/
template <typename FP = double>
struct RecoveredPerInfectedNoSymptoms {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "RecoveredPerInfectedNoSymptoms";
    }
};

/**
* @brief the risk of infection from symptomatic cases in the SECIR model
*/
template <typename FP = double>
struct RiskOfInfectionFromSymptomatic {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

/**
* @brief risk of infection from symptomatic cases increases as test and trace capacity is exceeded.
*/
template <typename FP = double>
struct MaxRiskOfInfectionFromSymptomatic {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "MaxRiskOfInfectionFromSymptomatic";
    }
};

/**
* @brief the percentage of hospitalized patients per infected patients in the SECIR model
*/
template <typename FP = double>
struct SeverePerInfectedSymptoms {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "SeverePerInfectedSymptoms";
    }
};

/**
* @brief the percentage of ICU patients per hospitalized patients in the SECIR model
*/
template <typename FP = double>
struct CriticalPerSevere {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "CriticalPerSevere";
    }
};

/**
* @brief the percentage of dead patients per ICU patients in the SECIR model
*/
template <typename FP = double>
struct DeathsPerCritical {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "DeathsPerCritical";
    }
};

/**
 * @brief Time in days between first and second vaccine dose.
 */
template <typename FP = double>
struct VaccinationGap {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 49.0);
    }
    static std::string name()
    {
        return "VaccinationGap";
    }
};

/**
 * @brief Time in days until first vaccine dose takes full effect.
 */
template <typename FP = double>
struct DaysUntilEffectivePartialImmunity {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 14.0);
    }
    static std::string name()
    {
        return "DaysUntilEffectivePartialImmunity";
    }
};

/**
 * @brief Time in days until second vaccine dose takes full effect.
 */
template <typename FP = double>
struct DaysUntilEffectiveImprovedImmunity {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 7.0);
    }
    static std::string name()
    {
        return "DaysUntilEffectiveImprovedImmunity";
    }
};

/**
* @brief Total number of first vaccinations up to the given day.
*/
template <typename FP = double>
struct DailyPartialVaccinations {
    using Type = CustomIndexArray<FP, AgeGroup, SimulationDay>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, SimulationDay(0)});
    }
    static std::string name()
    {
        return "DailyPartialVaccinations";
    }
};

/**
* @brief Total number of full vaccinations up to the given day.
*/
template <typename FP = double>
struct DailyFullVaccinations {
    using Type = CustomIndexArray<FP, AgeGroup, SimulationDay>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, SimulationDay(0)});
    }
    static std::string name()
    {
        return "DailyFullVaccinations";
    }
};

/**
 * @brief Factor to reduce infection risk for persons with partial immunity.
 */
template <typename FP = double>
struct ReducExposedPartialImmunity {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "ReducExposedPartialImmunity";
    }
};

/**
 * @brief Factor to reduce infection risk for persons with improved immunity.
 */
template <typename FP = double>
struct ReducExposedImprovedImmunity {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "ReducExposedImprovedImmunity";
    }
};

/**
 * @brief Factor to reduce risk of developing symptoms for persons with partial immunity.
 */
template <typename FP = double>
struct ReducInfectedSymptomsPartialImmunity {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "ReducInfectedSymptomsPartialImmunity";
    }
};

/**
 * @brief Factor to reduce risk of developing symptoms for persons with improved immunity.
 */
template <typename FP = double>
struct ReducInfectedSymptomsImprovedImmunity {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "ReducInfectedSymptomsImprovedImmunity";
    }
};

/**
 * @brief Factor to reduce risk of hospitalization for persons with partial immunity.
 * Also applies to ICU and Death risk.
 */
template <typename FP = double>
struct ReducInfectedSevereCriticalDeadPartialImmunity {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "ReducInfectedSevereCriticalDeadPartialImmunity";
    }
};

/**
 * @brief Factor to reduce risk of hospitalization for persons with improved immunity.
 */
template <typename FP = double>
struct ReducInfectedSevereCriticalDeadImprovedImmunity {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "ReducInfectedSevereCriticalDeadImprovedImmunity";
    }
};

/**
 * @brief Factor to reduce infectious time of persons with partial or improved immunity.
 */
template <typename FP = double>
struct ReducTimeInfectedMild {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.5);
    }
    static std::string name()
    {
        return "ReducTimeInfectedMild";
    }
};

/**
 * @brief Represents the relative infectiousness of a new variant.
 */
template <typename FP = double>
struct InfectiousnessNewVariant {
    using Type = CustomIndexArray<FP, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "InfectiousnessNewVariant";
    }
};

template <typename FP = double>
using ParametersBase = ParameterSet<
    StartDay, Seasonality<FP>, ICUCapacity<FP>, TestAndTraceCapacity<FP>, TestAndTraceCapacityMaxRiskNoSymptoms<FP>,
    TestAndTraceCapacityMaxRiskSymptoms<FP>, ContactPatterns<FP>, DynamicNPIsImplementationDelay<FP>,
    DynamicNPIsInfectedSymptoms<FP>, TimeExposed<FP>, TimeInfectedNoSymptoms<FP>, TimeInfectedSymptoms<FP>,
    TimeInfectedSevere<FP>, TimeInfectedCritical<FP>, TransmissionProbabilityOnContact<FP>,
    RelativeTransmissionNoSymptoms<FP>, RecoveredPerInfectedNoSymptoms<FP>, RiskOfInfectionFromSymptomatic<FP>,
    MaxRiskOfInfectionFromSymptomatic<FP>, SeverePerInfectedSymptoms<FP>, CriticalPerSevere<FP>, DeathsPerCritical<FP>,
    VaccinationGap<FP>, DaysUntilEffectivePartialImmunity<FP>, DaysUntilEffectiveImprovedImmunity<FP>,
    DailyFullVaccinations<FP>, DailyPartialVaccinations<FP>, ReducExposedPartialImmunity<FP>,
    ReducExposedImprovedImmunity<FP>, ReducInfectedSymptomsPartialImmunity<FP>,
    ReducInfectedSymptomsImprovedImmunity<FP>, ReducInfectedSevereCriticalDeadPartialImmunity<FP>,
    ReducInfectedSevereCriticalDeadImprovedImmunity<FP>, ReducTimeInfectedMild<FP>, InfectiousnessNewVariant<FP>,
    StartDayNewVariant>;

/**
 * @brief Parameters of an age-resolved SECIR/SECIHURD model with paths for partial and improved immunity through vaccination.
 */
template <typename FP = double>
class Parameters : public ParametersBase<FP>
{
public:
    Parameters(AgeGroup num_agegroups)
        : ParametersBase<FP>(num_agegroups)
        , m_num_groups{num_agegroups}
    {
    }

    AgeGroup get_num_groups() const
    {
        return m_num_groups;
    }

    /**
     * Percentage of infected commuters that are not detected.
     */
    double& get_commuter_nondetection()
    {
        return m_commuter_nondetection;
    }
    double get_commuter_nondetection() const
    {
        return m_commuter_nondetection;
    }

    /**
     * Time in simulation before which no infected commuters are detected.
     */
    double& get_start_commuter_detection()
    {
        return m_start_commuter_detection;
    }

    double get_start_commuter_detection() const
    {
        return m_start_commuter_detection;
    }

    /**
     * Time in simulation after which no infected commuters are detected.
     */
    double& get_end_commuter_detection()
    {
        return m_end_commuter_detection;
    }

    double get_end_commuter_detection() const
    {
        return m_end_commuter_detection;
    }

    /**
     * Time in simulation after which no dynamic NPIs are applied.
     */
    double& get_end_dynamic_npis()
    {
        return m_end_dynamic_npis;
    }
    double get_end_dynamic_npis() const
    {
        return m_end_dynamic_npis;
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and applies them, if they do not.
     * Time spans cannot be negative and probabilities can only take values between [0,1]. 
     *
     * Attention: This function should be used with care. It is necessary for some test problems to run through quickly,
     *            but in a manual execution of an example, check_constraints() may be preferred. Note that the apply_constraints()
     *            function can and will not set Parameters to meaningful values in an epidemiological or virological context,
     *            as all models are designed to be transferable to multiple diseases. Consequently, only acceptable
     *            (like 0 or 1 for probabilities or small positive values for time spans) values are set here and a manual adaptation
     *            may often be necessary to have set meaningful values.
     *
     * @return Returns true if one ore more constraint were corrected, false otherwise.  
     */
    bool apply_constraints()
    {
        int corrected = false;
        if (this->template get<Seasonality<FP>>() < 0.0 || this->template get<Seasonality<FP>>() > 0.5) {
            log_warning("Constraint check: Parameter Seasonality changed from {} to {}",
                        this->template get<Seasonality<FP>>(), 0);
            this->template set<Seasonality<FP>>(0);
            corrected = true;
        }

        if (this->template get<ICUCapacity<FP>>() < 0.0) {
            log_warning("Constraint check: Parameter ICUCapacity changed from {} to {}",
                        this->template get<ICUCapacity<FP>>(), 0);
            this->template set<ICUCapacity<FP>>(0);
            corrected = true;
        }

        if (this->template get<TestAndTraceCapacity<FP>>() < 0.0) {
            log_warning("Constraint check: Parameter TestAndTraceCapacity changed from {} to {}",
                        this->template get<TestAndTraceCapacity<FP>>(), 0);
            this->template set<TestAndTraceCapacity<FP>>(0);
            corrected = true;
        }

        if (this->template get<TestAndTraceCapacityMaxRiskSymptoms<FP>>() < 0.0) {
            log_warning("Constraint check: Parameter TestAndTraceCapacityMaxRiskSymptoms changed from {} to {}",
                        this->template get<TestAndTraceCapacityMaxRiskSymptoms<FP>>(), 0);
            this->template set<TestAndTraceCapacityMaxRiskSymptoms<FP>>(0);
            corrected = true;
        }

        if (this->template get<TestAndTraceCapacityMaxRiskNoSymptoms<FP>>() < 0.0) {
            log_warning("Constraint check: Parameter TestAndTraceCapacityMaxRiskNoSymptoms changed from {} to {}",
                        this->template get<TestAndTraceCapacityMaxRiskNoSymptoms<FP>>(), 0);
            this->template set<TestAndTraceCapacityMaxRiskNoSymptoms<FP>>(0);
            corrected = true;
        }

        if (this->template get<DynamicNPIsImplementationDelay<FP>>() < 0.0) {
            log_warning("Constraint check: Parameter DynamicNPIsImplementationDelay changed from {} to {}",
                        this->template get<DynamicNPIsImplementationDelay<FP>>(), 0);
            this->template set<DynamicNPIsImplementationDelay<FP>>(0);
            corrected = true;
        }

        const double tol_times = 1e-1; // accepted tolerance for compartment stays

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {
            if (this->template get<TimeExposed<FP>>()[i] < tol_times) {
                log_warning("Constraint check: Parameter TimeExposed changed from {:.4f} to {:.4f}. Please "
                            "note that unreasonably small compartment stays lead to massively increased run time. "
                            "Consider to cancel and reset parameters.",
                            this->template get<TimeExposed<FP>>()[i], tol_times);
                this->template get<TimeExposed<FP>>()[i] = tol_times;
                corrected                                = true;
            }

            if (this->template get<TimeInfectedNoSymptoms<FP>>()[i] < tol_times) {
                log_warning("Constraint check: Parameter TimeInfectedNoSymptoms changed from {:.4f} to {:.4f}. Please "
                            "note that unreasonably small compartment stays lead to massively increased run time. "
                            "Consider to cancel and reset parameters.",
                            this->template get<TimeInfectedNoSymptoms<FP>>()[i], tol_times);
                this->template get<TimeInfectedNoSymptoms<FP>>()[i] = tol_times;
                corrected                                           = true;
            }

            if (this->template get<TimeInfectedSymptoms<FP>>()[i] < tol_times) {
                log_warning("Constraint check: Parameter TimeInfectedSymptoms changed from {} to {}. Please "
                            "note that unreasonably small compartment stays lead to massively increased run time. "
                            "Consider to cancel and reset parameters.",
                            this->template get<TimeInfectedSymptoms<FP>>()[i], tol_times);
                this->template get<TimeInfectedSymptoms<FP>>()[i] = tol_times;
                corrected                                         = true;
            }

            if (this->template get<TimeInfectedSevere<FP>>()[i] < tol_times) {
                log_warning("Constraint check: Parameter TimeInfectedSevere changed from {} to {}. Please "
                            "note that unreasonably small compartment stays lead to massively increased run time. "
                            "Consider to cancel and reset parameters.",
                            this->template get<TimeInfectedSevere<FP>>()[i], tol_times);
                this->template get<TimeInfectedSevere<FP>>()[i] = tol_times;
                corrected                                       = true;
            }

            if (this->template get<TimeInfectedCritical<FP>>()[i] < tol_times) {
                log_warning("Constraint check: Parameter TimeInfectedCritical changed from {} to {}. Please "
                            "note that unreasonably small compartment stays lead to massively increased run time. "
                            "Consider to cancel and reset parameters.",
                            this->template get<TimeInfectedCritical<FP>>()[i], tol_times);
                this->template get<TimeInfectedCritical<FP>>()[i] = tol_times;
                corrected                                         = true;
            }

            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter TransmissionProbabilityOnContact changed from {} to {} ",
                            this->template get<TransmissionProbabilityOnContact<FP>>()[i], 0.0);
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] = 0.0;
                corrected                                                     = true;
            }

            if (this->template get<RelativeTransmissionNoSymptoms<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RelativeTransmissionNoSymptoms changed from {} to {} ",
                            this->template get<RelativeTransmissionNoSymptoms<FP>>()[i], 0);
                this->template get<RelativeTransmissionNoSymptoms<FP>>()[i] = 0;
                corrected                                                   = true;
            }

            if (this->template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] < 0.0 ||
                this->template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter RecoveredPerInfectedNoSymptoms changed from {} to {} ",
                            this->template get<RecoveredPerInfectedNoSymptoms<FP>>()[i], 0);
                this->template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] = 0;
                corrected                                                   = true;
            }

            if (this->template get<RiskOfInfectionFromSymptomatic<FP>>()[i] < 0.0 ||
                this->template get<RiskOfInfectionFromSymptomatic<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter RiskOfInfectionFromSymptomatic changed from {} to {}",
                            this->template get<RiskOfInfectionFromSymptomatic<FP>>()[i], 0);
                this->template get<RiskOfInfectionFromSymptomatic<FP>>()[i] = 0;
                corrected                                                   = true;
            }

            if (this->template get<SeverePerInfectedSymptoms<FP>>()[i] < 0.0 ||
                this->template get<SeverePerInfectedSymptoms<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter SeverePerInfectedSymptoms changed from {} to {}",
                            this->template get<SeverePerInfectedSymptoms<FP>>()[i], 0);
                this->template get<SeverePerInfectedSymptoms<FP>>()[i] = 0;
                corrected                                              = true;
            }

            if (this->template get<CriticalPerSevere<FP>>()[i] < 0.0 ||
                this->template get<CriticalPerSevere<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter CriticalPerSevere changed from {} to {}",
                            this->template get<CriticalPerSevere<FP>>()[i], 0);
                this->template get<CriticalPerSevere<FP>>()[i] = 0;
                corrected                                      = true;
            }

            if (this->template get<DeathsPerCritical<FP>>()[i] < 0.0 ||
                this->template get<DeathsPerCritical<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter DeathsPerCritical changed from {} to {}",
                            this->template get<DeathsPerCritical<FP>>()[i], 0);
                this->template get<DeathsPerCritical<FP>>()[i] = 0;
                corrected                                      = true;
            }

            if (this->template get<DaysUntilEffectivePartialImmunity<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter DaysUntilEffectivePartialImmunity changed from {} to {}",
                            this->template get<DaysUntilEffectivePartialImmunity<FP>>()[i], 0);
                this->template get<DaysUntilEffectivePartialImmunity<FP>>()[i] = 0;
                corrected                                                      = true;
            }
            if (this->template get<DaysUntilEffectiveImprovedImmunity<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter DaysUntilEffectiveImprovedImmunity changed from {} to {}",
                            this->template get<DaysUntilEffectiveImprovedImmunity<FP>>()[i], 0);
                this->template get<DaysUntilEffectiveImprovedImmunity<FP>>()[i] = 0;
                corrected                                                       = true;
            }

            if (this->template get<ReducExposedPartialImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducExposedPartialImmunity<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter ReducExposedPartialImmunity changed from {} to {}",
                            this->template get<ReducExposedPartialImmunity<FP>>()[i], 1);
                this->template get<ReducExposedPartialImmunity<FP>>()[i] = 1;
                corrected                                                = true;
            }
            if (this->template get<ReducExposedImprovedImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducExposedImprovedImmunity<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter ReducExposedImprovedImmunity changed from {} to {}",
                            this->template get<ReducExposedImprovedImmunity<FP>>()[i], 1);
                this->template get<ReducExposedImprovedImmunity<FP>>()[i] = 1;
                corrected                                                 = true;
            }
            if (this->template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter ReducInfectedSymptomsPartialImmunity changed from {} to {}",
                            this->template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i], 1);
                this->template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i] = 1;
                corrected                                                         = true;
            }
            if (this->template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter ReducInfectedSymptomsImprovedImmunity changed from {} to {}",
                            this->template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i], 1.0);
                this->template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i] = 1.0;
                corrected                                                          = true;
            }
            if (this->template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter ReducInfectedSevereCriticalDeadPartialImmunity changed from "
                            "{} to {}",
                            this->template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i], 1.0);
                this->template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i] = 1.0;
                corrected                                                                   = true;
            }
            if (this->template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter ReducInfectedSevereCriticalDeadImprovedImmunity changed from "
                            "{} to {}",
                            this->template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i], 1.0);
                this->template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i] = 1.0;
                corrected                                                                    = true;
            }
            if (this->template get<ReducTimeInfectedMild<FP>>()[i] <= 0.0 ||
                this->template get<ReducTimeInfectedMild<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter ReducTimeInfectedMild changed from {} to {}",
                            this->template get<ReducTimeInfectedMild<FP>>()[i], 1.0);
                this->template get<ReducTimeInfectedMild<FP>>()[i] = 1.0;
                corrected                                          = true;
            }
            if (this->template get<InfectiousnessNewVariant<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter InfectiousnessNewVariant changed from {} to {}",
                            this->template get<InfectiousnessNewVariant<FP>>()[i], 1.0);
                this->template get<InfectiousnessNewVariant<FP>>()[i] = 1.0;
                corrected                                             = true;
            }
            if (this->template get<VaccinationGap<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter VaccinationGap changed from {} to {}",
                            this->template get<VaccinationGap<FP>>()[i], 0);
                this->template get<VaccinationGap<FP>>()[i] = 0;
                corrected                                   = true;
            }
        }
        return corrected;
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and logs an error 
     * if constraints are not satisfied.
     * @return Returns true if one constraint is not satisfied, otherwise false.    
     */
    bool check_constraints() const
    {
        const double tol_times = 1e-1; // accepted tolerance for compartment stays
        if (this->template get<Seasonality<FP>>() < 0.0 || this->template get<Seasonality<FP>>() > 0.5) {
            log_error("Constraint check: Parameter Seasonality smaller {} or larger {}", 0, 0.5);
            return true;
        }

        if (this->template get<ICUCapacity<FP>>() < 0.0) {
            log_error("Constraint check: Parameter ICUCapacity smaller {}", 0);
            return true;
        }

        if (this->template get<TestAndTraceCapacity<FP>>() < 0.0) {
            log_error("Constraint check: Parameter TestAndTraceCapacity smaller {}", 0);
            return true;
        }

        if (this->template get<TestAndTraceCapacityMaxRiskSymptoms<FP>>() < 0.0) {
            log_error("Constraint check: Parameter TestAndTraceCapacityMaxRiskSymptoms smaller {}", 0);
            return true;
        }

        if (this->template get<TestAndTraceCapacityMaxRiskNoSymptoms<FP>>() < 0.0) {
            log_error("Constraint check: Parameter TestAndTraceCapacityMaxRiskNoSymptoms smaller {}", 0);
            return true;
        }

        if (this->template get<DynamicNPIsImplementationDelay<FP>>() < 0.0) {
            log_error("Constraint check: Parameter DynamicNPIsImplementationDelay smaller {:d}", 0);
            return true;
        }

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {

            if (this->template get<TimeExposed<FP>>()[i] < tol_times) {
                log_error("Constraint check: Parameter TimeExposed {:.4f} smaller {:.4f}. Please "
                          "note that unreasonably small compartment stays lead to massively increased run time. "
                          "Consider to cancel and reset parameters.",
                          this->template get<TimeExposed<FP>>()[i], tol_times);
                return true;
            }

            if (this->template get<TimeInfectedNoSymptoms<FP>>()[i] < tol_times) {
                log_error("Constraint check: Parameter TimeInfectedNoSymptoms {:.4f} smaller {:.4f}. Please "
                          "note that unreasonably small compartment stays lead to massively increased run time. "
                          "Consider to cancel and reset parameters.",
                          this->template get<TimeInfectedNoSymptoms<FP>>()[i], tol_times);
                return true;
            }

            if (this->template get<TimeInfectedSymptoms<FP>>()[i] < tol_times) {
                log_error("Constraint check: Parameter TimeInfectedSymptoms {} smaller {}. Please "
                          "note that unreasonably small compartment stays lead to massively increased run time. "
                          "Consider to cancel and reset parameters.",
                          this->template get<TimeInfectedSymptoms<FP>>()[i], tol_times);
                return true;
            }

            if (this->template get<TimeInfectedSevere<FP>>()[i] < tol_times) {
                log_error("Constraint check: Parameter TimeInfectedSevere {} smaller {}. Please "
                          "note that unreasonably small compartment stays lead to massively increased run time. "
                          "Consider to cancel and reset parameters.",
                          this->template get<TimeInfectedSevere<FP>>()[i], tol_times);
                return true;
            }

            if (this->template get<TimeInfectedCritical<FP>>()[i] < tol_times) {
                log_error("Constraint check: Parameter TimeInfectedCritical {} smaller {}. Please "
                          "note that unreasonably small compartment stays lead to massively increased run time. "
                          "Consider to cancel and reset parameters.",
                          this->template get<TimeInfectedCritical<FP>>()[i], tol_times);
                return true;
            }

            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContact smaller {} or larger {}", 0, 1);
                return true;
            }

            if (this->template get<RelativeTransmissionNoSymptoms<FP>>()[i] < 0.0) {
                log_error("Constraint check: Parameter RelativeTransmissionNoSymptoms smaller {}", 0);
                return true;
            }

            if (this->template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] < 0.0 ||
                this->template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter RecoveredPerInfectedNoSymptoms smaller {} or larger {}", 0, 1);
                return true;
            }

            if (this->template get<RiskOfInfectionFromSymptomatic<FP>>()[i] < 0.0 ||
                this->template get<RiskOfInfectionFromSymptomatic<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter RiskOfInfectionFromSymptomatic smaller {} or larger {}", 0, 1);
                return true;
            }

            if (this->template get<SeverePerInfectedSymptoms<FP>>()[i] < 0.0 ||
                this->template get<SeverePerInfectedSymptoms<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter SeverePerInfectedSymptoms smaller {} or larger {}", 0, 1);
                return true;
            }

            if (this->template get<CriticalPerSevere<FP>>()[i] < 0.0 ||
                this->template get<CriticalPerSevere<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter CriticalPerSevere smaller {} or larger {}", 0, 1);
                return true;
            }

            if (this->template get<DeathsPerCritical<FP>>()[i] < 0.0 ||
                this->template get<DeathsPerCritical<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter DeathsPerCritical smaller {} or larger {}", 0, 1);
                return true;
            }

            if (this->template get<VaccinationGap<FP>>()[i] < 1) {
                log_error("Constraint check: Parameter VaccinationGap smaller {}", 1);
                return true;
            }

            if (this->template get<DaysUntilEffectivePartialImmunity<FP>>()[i] < 0.0) {
                log_error("Constraint check: Parameter DaysUntilEffectivePartialImmunity smaller {}", 0);
                return true;
            }
            if (this->template get<DaysUntilEffectiveImprovedImmunity<FP>>()[i] < 0.0) {
                log_error("Constraint check: Parameter DaysUntilEffectiveImprovedImmunity smaller {}", 0);
                return true;
            }

            if (this->template get<ReducExposedPartialImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducExposedPartialImmunity<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter ReducExposedPartialImmunity smaller {} or larger {}", 0, 1);
                return true;
            }
            if (this->template get<ReducExposedImprovedImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducExposedImprovedImmunity<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter ReducExposedImprovedImmunity smaller {} or larger {}", 0, 1);
                return true;
            }
            if (this->template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter ReducInfectedSymptomsPartialImmunity smaller {} or larger {}", 0,
                          1);
                return true;
            }
            if (this->template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter ReducInfectedSymptomsImprovedImmunity smaller {} or larger {}",
                          0, 1);
                return true;
            }
            if (this->template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter ReducInfectedSevereCriticalDeadPartialImmunity smaller {} or "
                          "larger {}",
                          0, 1);
                return true;
            }
            if (this->template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i] <= 0.0 ||
                this->template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter ReducInfectedSevereCriticalDeadImprovedImmunity smaller {} or "
                          "larger {}",
                          0, 1);
                return true;
            }
            if (this->template get<ReducTimeInfectedMild<FP>>()[i] <= 0.0 ||
                this->template get<ReducTimeInfectedMild<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter ReducTimeInfectedMild smaller {} or larger {}", 0, 1);
                return true;
            }
            if (this->template get<InfectiousnessNewVariant<FP>>()[i] < 0.0) {
                log_error("Constraint check: Parameter InfectiousnessNewVariant smaller {}", 0);
                return true;
            }
        }
        return false;
    }

private:
    Parameters(ParametersBase<FP>&& base)
        : ParametersBase<FP>(std::move(base))
        , m_num_groups(this->template get<ContactPatterns<FP>>().get_cont_freq_mat().get_num_groups())
    {
    }

public:
    /**
     * deserialize an object of this class.
     * @see epi::deserialize
     */
    template <class IOContext>
    static IOResult<Parameters> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& base, ParametersBase<FP>::deserialize(io));
        return success(Parameters(std::move(base)));
    }

private:
    AgeGroup m_num_groups;
    double m_commuter_nondetection    = 0.0;
    double m_start_commuter_detection = 0.0;
    double m_end_commuter_detection   = 0.0;
    double m_end_dynamic_npis         = std::numeric_limits<double>::max();
};

} // namespace osecirvvs
} // namespace mio

#endif // MIO_ODE_SECIRVVS_PARAMETERS_H
