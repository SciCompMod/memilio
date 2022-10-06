/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
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
#ifndef ODESECIRVVS_PARAMETERS_H
#define ODESECIRVVS_PARAMETERS_H

#include "memilio/math/eigen.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/simulation_day.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/dynamic_npis.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/custom_index_array.h"

#include <vector>

namespace mio
{
namespace osecirvvs
{

/**
* @brief the start day in the SECIR model
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
* @brief the seasonality in the SECIR model
* the seasonality is given as (1+k*sin()) where the sine
* curve is below one in summer and above one in winter
*/
struct Seasonality {
    using Type = UncertainValue;
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
struct ICUCapacity {
    using Type = UncertainValue;
    static Type get_default(AgeGroup)
    {
        return Type(std::numeric_limits<double>::max());
    }
    static std::string name()
    {
        return "ICUCapacity";
    }
};

/**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
struct TestAndTraceCapacity {
    using Type = UncertainValue;
    static Type get_default(AgeGroup)
    {
        return Type(std::numeric_limits<double>::max());
    }
    static std::string name()
    {
        return "TestAndTraceCapacity";
    }
};

/**
 * @brief the contact patterns within the society are modelled using an UncertainContactMatrix
 */
struct ContactPatterns {
    using Type = UncertainContactMatrix;
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
 * @brief the NPIs that are enacted if certain infection thresholds are exceeded.
 */
struct DynamicNPIsInfectedSymptoms {
    using Type = DynamicNPIs;
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
* @brief the incubation time in the SECIR model
* @param tinc incubation time in day unit
*/
struct IncubationTime {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "IncubationTime";
    }
};

/**
* @brief the infectious time for symptomatic cases that are infected but
*        who do not need to be hsopitalized in the SECIR model in day unit
*/
struct TimeInfectedSymptoms {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
 * @brief the serial interval in the SECIR model in day unit
 */
struct SerialInterval {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "SerialInterval";
    }
};

/**
 * @brief the time people are 'simply' hospitalized before returning home in the SECIR model
 *        in day unit
 */
struct TimeInfectedSevere {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct TimeInfectedCritical {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct TransmissionProbabilityOnContact {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
* @brief the relative carrier infectability
*/
struct RelativeTransmissionNoSymptoms {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct AsymptoticCasesPerInfectious {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "AsymptoticCasesPerInfectious";
    }
};

/**
* @brief the risk of infection from symptomatic cases in the SECIR model
*/
struct RiskOfInfectionFromSympomatic {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSympomatic";
    }
};

/**
* @brief risk of infection from symptomatic cases increases as test and trace capacity is exceeded.
*/
struct MaxRiskOfInfectionFromSympomatic {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "MaxRiskOfInfectionFromSympomatic";
    }
};

/**
* @brief the percentage of hospitalized patients per infected patients in the SECIR model
*/
struct SeverePerInfectedSymptoms {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct CriticalPerSevere {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct DeathsPerCritical {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct VaccinationGap {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct DaysUntilEffectivePartialImmunity {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct DaysUntilEffectiveImprovedImmunity {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct DailyFirstVaccination {
    using Type = CustomIndexArray<double, AgeGroup, SimulationDay>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, SimulationDay(0)});
    }
    static std::string name()
    {
        return "DailyFirstVaccination";
    }
};

/**
* @brief Total number of full vaccinations up to the given day.
*/
struct DailyFullVaccination {
    using Type = CustomIndexArray<double, AgeGroup, SimulationDay>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, SimulationDay(0)});
    }
    static std::string name()
    {
        return "DailyFullVaccination";
    }
};

/**
 * @brief Factor to reduce infection risk for persons with partial immunity.
 */
struct ExposedFactorPartialImmunity {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "ExposedFactorPartialImmunity";
    }
};

/**
 * @brief Factor to reduce infection risk for persons with improved immunity.
 */
struct ExposedFactorImprovedImmunity {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "ExposedFactorImprovedImmunity";
    }
};

/**
 * @brief Factor to reduce risk of developing symptoms for persons with partial immunity.
 */
struct InfectedFactorPartialImmunity {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "InfectedFactorPartialImmunity";
    }
};

/**
 * @brief Factor to reduce risk of developing symptoms for persons with improved immunity.
 */
struct InfectedFactorImprovedImmunity {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "InfectedFactorImprovedImmunity";
    }
};

/**
 * @brief Factor to reduce risk of hospitalization for persons with partial immunity.
 * Also applies to ICU and Death risk.
 */
struct HospitalizedFactorPartialImmunity {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "HospitalizedFactorPartialImmunity";
    }
};

/**
 * @brief Factor to reduce risk of hospitalization for persons with improved immunity.
 */
struct HospitalizedFactorImprovedImmunity {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "HospitalizedFactorImprovedImmunity";
    }
};

/**
 * @brief Factor to reduce infectious time of persons with partial or improved immunity.
 */
struct InfectiousTimeFactorImmune {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.5);
    }
    static std::string name()
    {
        return "InfectiousTimeFactorImmune";
    }
};

/**
 * @brief Infectiousness of variant B117.
 */
struct BaseInfectiousnessB117 {
    using Type = CustomIndexArray<double, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "BaseInfectiousnessB117";
    }
};

/**
 * @brief Infectiousness of variant B161.
 */
struct BaseInfectiousnessB161 {
    using Type = CustomIndexArray<double, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "BaseInfectiousnessB161";
    }
};

using ParametersBase =
    ParameterSet<StartDay, Seasonality, ICUCapacity, TestAndTraceCapacity, ContactPatterns,
                 DynamicNPIsInfectedSymptoms, IncubationTime, TimeInfectedSymptoms, SerialInterval,
                 TimeInfectedSevere, TimeInfectedCritical,
                 TransmissionProbabilityOnContact, RelativeTransmissionNoSymptoms, AsymptoticCasesPerInfectious,
                 RiskOfInfectionFromSympomatic, MaxRiskOfInfectionFromSympomatic, SeverePerInfectedSymptoms,
                 CriticalPerSevere, DeathsPerCritical, VaccinationGap, DaysUntilEffectivePartialImmunity,
                 DaysUntilEffectiveImprovedImmunity, DailyFullVaccination, DailyFirstVaccination,
                 ExposedFactorPartialImmunity, ExposedFactorImprovedImmunity, InfectedFactorPartialImmunity,
                 InfectedFactorImprovedImmunity, HospitalizedFactorPartialImmunity, HospitalizedFactorImprovedImmunity,
                 InfectiousTimeFactorImmune, BaseInfectiousnessB117, BaseInfectiousnessB161>;

/**
 * @brief Parameters of an age-resolved SECIR/SECIHURD model with paths for partial and improved immunity through vaccination.
 */
class Parameters : public ParametersBase
{
public:
    Parameters(AgeGroup num_agegroups)
        : ParametersBase(num_agegroups)
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
     * @brief checks whether all Parameters satisfy their corresponding constraints and applies them, if they do not
     */
    void apply_constraints()
    {
        if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
            log_warning("Constraint check: Parameter Seasonality changed from {:0.4f} to {:d}",
                        this->get<Seasonality>(), 0);
            this->set<Seasonality>(0);
        }

        if (this->get<ICUCapacity>() < 0.0) {
            log_warning("Constraint check: Parameter ICUCapacity changed from {:0.4f} to {:d}",
                        this->get<ICUCapacity>(), 0);
            this->set<ICUCapacity>(0);
        }

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {

            if (this->get<IncubationTime>()[i] < 2.0) {
                log_warning("Constraint check: Parameter IncubationTime changed from {:.4f} to {:.4f}",
                            this->get<IncubationTime>()[i], 2.0);
                this->get<IncubationTime>()[i] = 2.0;
            }

            if (2 * this->get<SerialInterval>()[i] < this->get<IncubationTime>()[i] + 1.0) {
                log_warning("Constraint check: Parameter SerialInterval changed from {:.4f} to {:.4f}",
                            this->get<SerialInterval>()[i], 0.5 * this->get<IncubationTime>()[i] + 0.5);
                this->get<SerialInterval>()[i] = 0.5 * this->get<IncubationTime>()[i] + 0.5;
            }
            else if (this->get<SerialInterval>()[i] > this->get<IncubationTime>()[i] - 0.5) {
                log_warning("Constraint check: Parameter SerialInterval changed from {:.4f} to {:.4f}",
                            this->get<SerialInterval>()[i], this->get<IncubationTime>()[i] - 0.5);
                this->get<SerialInterval>()[i] = this->get<IncubationTime>()[i] - 0.5;
            }

            if (this->get<TimeInfectedSymptoms>()[i] < 1.0) {
                log_warning("Constraint check: Parameter TimeInfectedSymptoms changed from {:.4f} to {:.4f}",
                            this->get<TimeInfectedSymptoms>()[i], 1.0);
                this->get<TimeInfectedSymptoms>()[i] = 1.0;
            }

            if (this->get<TimeInfectedSevere>()[i] < 1.0) {
                log_warning("Constraint check: Parameter TimeInfectedSevere changed from {:.4f} to {:.4f}",
                            this->get<TimeInfectedSevere>()[i], 1.0);
                this->get<TimeInfectedSevere>()[i] = 1.0;
            }

            if (this->get<TimeInfectedCritical>()[i] < 1.0) {
                log_warning("Constraint check: Parameter TimeInfectedCritical changed from {:.4f} to {:.4f}",
                            this->get<TimeInfectedCritical>()[i], 1.0);
                this->get<TimeInfectedCritical>()[i] = 1.0;
            }

            if (this->get<TransmissionProbabilityOnContact>()[i] < 0.0) {
                log_warning("Constraint check: Parameter TransmissionProbabilityOnContact changed from {:0.4f} to {:d} ",
                            this->get<TransmissionProbabilityOnContact>()[i], 0);
                this->get<TransmissionProbabilityOnContact>()[i] = 0;
            }

            if (this->get<RelativeTransmissionNoSymptoms>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RelativeTransmissionNoSymptoms changed from {:0.4f} to {:d} ",
                            this->get<RelativeTransmissionNoSymptoms>()[i], 0);
                this->get<RelativeTransmissionNoSymptoms>()[i] = 0;
            }

            if (this->get<AsymptoticCasesPerInfectious>()[i] < 0.0 ||
                this->get<AsymptoticCasesPerInfectious>()[i] > 1.0) {
                log_warning("Constraint check: Parameter AsymptoticCasesPerInfectious changed from {:0.4f} to {:d} ",
                            this->get<AsymptoticCasesPerInfectious>()[i], 0);
                this->get<AsymptoticCasesPerInfectious>()[i] = 0;
            }

            if (this->get<RiskOfInfectionFromSympomatic>()[i] < 0.0 ||
                this->get<RiskOfInfectionFromSympomatic>()[i] > 1.0) {
                log_warning("Constraint check: Parameter RiskOfInfectionFromSympomatic changed from {:0.4f} to {:d}",
                            this->get<RiskOfInfectionFromSympomatic>()[i], 0);
                this->get<RiskOfInfectionFromSympomatic>()[i] = 0;
            }

            if (this->get<SeverePerInfectedSymptoms>()[i] < 0.0 ||
                this->get<SeverePerInfectedSymptoms>()[i] > 1.0) {
                log_warning("Constraint check: Parameter SeverePerInfectedSymptoms changed from {:0.4f} to {:d}",
                            this->get<SeverePerInfectedSymptoms>()[i], 0);
                this->get<SeverePerInfectedSymptoms>()[i] = 0;
            }

            if (this->get<CriticalPerSevere>()[i] < 0.0 || this->get<CriticalPerSevere>()[i] > 1.0) {
                log_warning("Constraint check: Parameter CriticalPerSevere changed from {:0.4f} to {:d}",
                            this->get<CriticalPerSevere>()[i], 0);
                this->get<CriticalPerSevere>()[i] = 0;
            }

            if (this->get<DeathsPerCritical>()[i] < 0.0 || this->get<DeathsPerCritical>()[i] > 1.0) {
                log_warning("Constraint check: Parameter DeathsPerCritical changed from {:0.4f} to {:d}",
                            this->get<DeathsPerCritical>()[i], 0);
                this->get<DeathsPerCritical>()[i] = 0;
            }
        }
    }

    /**
     * @brief checks whether all Parameters satisfy their corresponding constraints and throws errors, if they do not
     */
    void check_constraints() const
    {
        if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
            log_warning("Constraint check: Parameter m_seasonality smaller {:d} or larger {:d}", 0, 0.5);
        }

        if (this->get<ICUCapacity>() < 0.0) {
            log_warning("Constraint check: Parameter m_icu_capacity smaller {:d}", 0);
        }

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {

            if (this->get<IncubationTime>()[i] < 2.0) {
                log_error("Constraint check: Parameter IncubationTime {:.4f} smaller {:.4f}",
                          this->get<IncubationTime>()[i], 2.0);
            }

            if (2 * this->get<SerialInterval>()[i] < this->get<IncubationTime>()[i] + 1.0) {
                log_error("Constraint check: Parameter SerialInterval {:.4f} smaller {:.4f}",
                          this->get<SerialInterval>()[i], 0.5 * this->get<IncubationTime>()[i] + 0.5);
            }
            else if (this->get<SerialInterval>()[i] > this->get<IncubationTime>()[i] - 0.5) {
                log_error("Constraint check: Parameter SerialInterval {:.4f} smaller {:.4f}",
                          this->get<SerialInterval>()[i], this->get<IncubationTime>()[i] - 0.5);
            }

            if (this->get<TimeInfectedSymptoms>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedSymptoms {:.4f} smaller {:.4f}",
                          this->get<TimeInfectedSymptoms>()[i], 1.0);
            }

            if (this->get<TimeInfectedSevere>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedSevere {:.4f} smaller {:.4f}",
                          this->get<TimeInfectedSevere>()[i], 1.0);
            }
            
            if (this->get<TimeInfectedCritical>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedCritical {:.4f} smaller {:.4f}",
                          this->get<TimeInfectedCritical>()[i], 1.0);
            }

            if (this->get<TransmissionProbabilityOnContact>()[i] < 0.0) {
                log_warning("Constraint check: Parameter TransmissionProbabilityOnContact smaller {:d}", 0);
            }

            if (this->get<RelativeTransmissionNoSymptoms>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RelativeTransmissionNoSymptoms smaller {:d}", 0);
            }

            if (this->get<AsymptoticCasesPerInfectious>()[i] < 0.0 ||
                this->get<AsymptoticCasesPerInfectious>()[i] > 1.0) {
                log_warning("Constraint check: Parameter AsymptoticCasesPerInfectious smaller {:d} or larger {:d}", 0,
                            1);
            }

            if (this->get<RiskOfInfectionFromSympomatic>()[i] < 0.0 ||
                this->get<RiskOfInfectionFromSympomatic>()[i] > 1.0) {
                log_warning("Constraint check: Parameter RiskOfInfectionFromSympomatic smaller {:d} or larger {:d}", 0,
                            1);
            }

            if (this->get<SeverePerInfectedSymptoms>()[i] < 0.0 ||
                this->get<SeverePerInfectedSymptoms>()[i] > 1.0) {
                log_warning("Constraint check: Parameter SeverePerInfectedSymptoms smaller {:d} or larger {:d}", 0,
                            1);
            }

            if (this->get<CriticalPerSevere>()[i] < 0.0 || this->get<CriticalPerSevere>()[i] > 1.0) {
                log_warning("Constraint check: Parameter CriticalPerSevere smaller {:d} or larger {:d}", 0, 1);
            }

            if (this->get<DeathsPerCritical>()[i] < 0.0 || this->get<DeathsPerCritical>()[i] > 1.0) {
                log_warning("Constraint check: Parameter DeathsPerCritical smaller {:d} or larger {:d}", 0, 1);
            }
        }
    }

private:
    Parameters(ParametersBase&& base)
        : ParametersBase(std::move(base))
        , m_num_groups(get<ContactPatterns>().get_cont_freq_mat().get_num_groups())
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
        BOOST_OUTCOME_TRY(base, ParametersBase::deserialize(io));
        return success(Parameters(std::move(base)));
    }

private:
    AgeGroup m_num_groups;
    double m_commuter_nondetection    = 0.0;
    double m_start_commuter_detection = 0.0;
    double m_end_commuter_detection   = 0.0;
    double m_end_dynamic_npis         = 0.0;
};

} // namespace osecirvvs
} // namespace mio

#endif // ODESECIRVVS_PARAMETERS_H
