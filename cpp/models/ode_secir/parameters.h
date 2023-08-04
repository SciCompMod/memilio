/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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
#ifndef SECIR_PARAMETERS_H
#define SECIR_PARAMETERS_H

#include "memilio/math/eigen.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/dynamic_npis.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/custom_index_array.h"

#include <vector>

namespace mio
{
namespace osecir
{

/*******************************************
 * Define Parameters of the SECIHURD model *
 *******************************************/

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
 * @brief the incubation time in the SECIR model
 * @param tinc incubation time in day unit
 */
struct IncubationTime {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 2.);
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
        return Type(size, 1.5);
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
* @brief the relative InfectedNoSymptoms infectability
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
struct RecoveredPerInfectedNoSymptoms {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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
struct RiskOfInfectionFromSymptomatic {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

/**
* @brief risk of infection from symptomatic cases increases as test and trace capacity is exceeded.
*/
struct MaxRiskOfInfectionFromSymptomatic {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
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

using ParametersBase =
    ParameterSet<StartDay, Seasonality, ICUCapacity, TestAndTraceCapacity, ContactPatterns, DynamicNPIsInfectedSymptoms,
                 IncubationTime, TimeInfectedSymptoms, SerialInterval, TimeInfectedSevere, TimeInfectedCritical,
                 TransmissionProbabilityOnContact, RelativeTransmissionNoSymptoms, RecoveredPerInfectedNoSymptoms,
                 RiskOfInfectionFromSymptomatic, MaxRiskOfInfectionFromSymptomatic, SeverePerInfectedSymptoms,
                 CriticalPerSevere, DeathsPerCritical>;

/**
 * @brief Parameters of an age-resolved SECIR/SECIHURD model.
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
        double tol_times = 1e-4; // accepted tolerance for compartment stays

        int corrected = false;
        if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
            log_warning("Constraint check: Parameter Seasonality changed from {:0.4f} to {:d}",
                        this->get<Seasonality>(), 0);
            this->set<Seasonality>(0);
            corrected = true;
        }

        if (this->get<ICUCapacity>() < 0.0) {
            log_warning("Constraint check: Parameter ICUCapacity changed from {:0.4f} to {:d}",
                        this->get<ICUCapacity>(), 0);
            this->set<ICUCapacity>(0);
            corrected = true;
        }

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {

            if (this->get<IncubationTime>()[i] < 2 * tol_times) {
                log_warning("Constraint check: Parameter IncubationTime changed from {:.4f} to {:.4f}",
                            this->get<IncubationTime>()[i], tol_times);
                this->get<IncubationTime>()[i] = 2 * tol_times;
                corrected                      = true;
            }

            if (2 * this->get<SerialInterval>()[i] < this->get<IncubationTime>()[i] + tol_times) {
                log_warning("Constraint check: Parameter SerialInterval changed from {:.4f} to {:.4f}",
                            this->get<SerialInterval>()[i], 0.5 * this->get<IncubationTime>()[i] + tol_times / 2);
                this->get<SerialInterval>()[i] = 0.5 * this->get<IncubationTime>()[i] + tol_times / 2;
                corrected                      = true;
            }
            else if (this->get<SerialInterval>()[i] > this->get<IncubationTime>()[i] - tol_times / 2) {
                log_warning("Constraint check: Parameter SerialInterval changed from {:.4f} to {:.4f}",
                            this->get<SerialInterval>()[i], this->get<IncubationTime>()[i] - tol_times / 2);
                this->get<SerialInterval>()[i] = this->get<IncubationTime>()[i] - tol_times / 2;
                corrected                      = true;
            }

            if (this->get<TimeInfectedSymptoms>()[i] < tol_times) {
                log_warning("Constraint check: Parameter TimeInfectedSymptoms changed from {:.4f} to {:.4f}. Please "
                            "note that unreasonably small compartment stays lead to massively increased run time. "
                            "Consider to cancel and reset parameters.",
                            this->get<TimeInfectedSymptoms>()[i], tol_times);
                this->get<TimeInfectedSymptoms>()[i] = tol_times;
                corrected                            = true;
            }

            if (this->get<TimeInfectedSevere>()[i] < tol_times) {
                log_warning("Constraint check: Parameter TimeInfectedSevere changed from {:.4f} to {:.4f}. Please note "
                            "that unreasonably small compartment stays lead to massively increased run time. Consider "
                            "to cancel and reset parameters.",
                            this->get<TimeInfectedSevere>()[i], tol_times);
                this->get<TimeInfectedSevere>()[i] = tol_times;
                corrected                          = true;
            }

            if (this->get<TimeInfectedCritical>()[i] < tol_times) {
                log_warning("Constraint check: Parameter TimeInfectedCritical changed from {:.4f} to {:.4f}. Please "
                            "note that unreasonably small compartment stays lead to massively increased run time. "
                            "Consider to cancel and reset parameters.",
                            this->get<TimeInfectedCritical>()[i], tol_times);
                this->get<TimeInfectedCritical>()[i] = tol_times;
                corrected                            = true;
            }

            if (this->get<TransmissionProbabilityOnContact>()[i] < 0.0 ||
                this->get<TransmissionProbabilityOnContact>()[i] > 1.0) {
                log_warning(
                    "Constraint check: Parameter TransmissionProbabilityOnContact changed from {:0.4f} to {:d} ",
                    this->get<TransmissionProbabilityOnContact>()[i], 0.0);
                this->get<TransmissionProbabilityOnContact>()[i] = 0.0;
                corrected                                        = true;
            }

            if (this->get<RelativeTransmissionNoSymptoms>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RelativeTransmissionNoSymptoms changed from {:0.4f} to {:d} ",
                            this->get<RelativeTransmissionNoSymptoms>()[i], 0);
                this->get<RelativeTransmissionNoSymptoms>()[i] = 0;
                corrected                                      = true;
            }

            if (this->get<RecoveredPerInfectedNoSymptoms>()[i] < 0.0 ||
                this->get<RecoveredPerInfectedNoSymptoms>()[i] > 1.0) {
                log_warning("Constraint check: Parameter RecoveredPerInfectedNoSymptoms changed from {:0.4f} to {:d} ",
                            this->get<RecoveredPerInfectedNoSymptoms>()[i], 0);
                this->get<RecoveredPerInfectedNoSymptoms>()[i] = 0;
                corrected                                      = true;
            }

            if (this->get<RiskOfInfectionFromSymptomatic>()[i] < 0.0 ||
                this->get<RiskOfInfectionFromSymptomatic>()[i] > 1.0) {
                log_warning("Constraint check: Parameter RiskOfInfectionFromSymptomatic changed from {:0.4f} to {:d}",
                            this->get<RiskOfInfectionFromSymptomatic>()[i], 0);
                this->get<RiskOfInfectionFromSymptomatic>()[i] = 0;
                corrected                                      = true;
            }

            if (this->get<SeverePerInfectedSymptoms>()[i] < 0.0 || this->get<SeverePerInfectedSymptoms>()[i] > 1.0) {
                log_warning("Constraint check: Parameter SeverePerInfectedSymptoms changed from {:0.4f} to {:d}",
                            this->get<SeverePerInfectedSymptoms>()[i], 0);
                this->get<SeverePerInfectedSymptoms>()[i] = 0;
                corrected                                 = true;
            }

            if (this->get<CriticalPerSevere>()[i] < 0.0 || this->get<CriticalPerSevere>()[i] > 1.0) {
                log_warning("Constraint check: Parameter CriticalPerSevere changed from {:0.4f} to {:d}",
                            this->get<CriticalPerSevere>()[i], 0);
                this->get<CriticalPerSevere>()[i] = 0;
                corrected                         = true;
            }

            if (this->get<DeathsPerCritical>()[i] < 0.0 || this->get<DeathsPerCritical>()[i] > 1.0) {
                log_warning("Constraint check: Parameter DeathsPerCritical changed from {:0.4f} to {:d}",
                            this->get<DeathsPerCritical>()[i], 0);
                this->get<DeathsPerCritical>()[i] = 0;
                corrected                         = true;
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
        if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
            log_error("Constraint check: Parameter Seasonality smaller {:d} or larger {:d}", 0, 0.5);
            return true;
        }

        if (this->get<ICUCapacity>() < 0.0) {
            log_error("Constraint check: Parameter ICUCapacity smaller {:d}", 0);
            return true;
        }

        double tol_times = 1e-4; // accepted tolerance for compartment stays

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {

            if (this->get<IncubationTime>()[i] < 2 * tol_times) {
                log_error("Constraint check: Parameter IncubationTime {:.4f} smaller {:.4f}",
                          this->get<IncubationTime>()[i], 2 * tol_times);
                return true;
            }

            if (2 * this->get<SerialInterval>()[i] < this->get<IncubationTime>()[i] + tol_times) {
                log_error("Constraint check: Parameter SerialInterval {:.4f} smaller {:.4f}",
                          this->get<SerialInterval>()[i], 0.5 * this->get<IncubationTime>()[i] + tol_times / 2);
                return true;
            }
            else if (this->get<SerialInterval>()[i] > this->get<IncubationTime>()[i] - tol_times / 2) {
                log_error("Constraint check: Parameter SerialInterval {:.4f} greater {:.4f}",
                          this->get<SerialInterval>()[i], this->get<IncubationTime>()[i] - tol_times / 2);
                return true;
            }

            if (this->get<TimeInfectedSymptoms>()[i] < tol_times) {
                log_error("Constraint check: Parameter TimeInfectedSymptoms {:.4f} smaller {:.4f}. Please "
                          "note that unreasonably small compartment stays lead to massively increased run time. "
                          "Consider to cancel and reset parameters.",
                          this->get<TimeInfectedSymptoms>()[i], tol_times);
                return true;
            }

            if (this->get<TimeInfectedSevere>()[i] < tol_times) {
                log_error("Constraint check: Parameter TimeInfectedSevere {:.4f} smaller {:.4f}. Please "
                          "note that unreasonably small compartment stays lead to massively increased run time. "
                          "Consider to cancel and reset parameters.",
                          this->get<TimeInfectedSevere>()[i], tol_times);
                return true;
            }

            if (this->get<TimeInfectedCritical>()[i] < tol_times) {
                log_error("Constraint check: Parameter TimeInfectedCritical {:.4f} smaller {:.4f}. Please "
                          "note that unreasonably small compartment stays lead to massively increased run time. "
                          "Consider to cancel and reset parameters.",
                          this->get<TimeInfectedCritical>()[i], tol_times);
                return true;
            }

            if (this->get<TransmissionProbabilityOnContact>()[i] < 0.0 ||
                this->get<TransmissionProbabilityOnContact>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContact smaller {:d} or larger {:d}", 0,
                          1);
                return true;
            }

            if (this->get<RelativeTransmissionNoSymptoms>()[i] < 0.0) {
                log_error("Constraint check: Parameter RelativeTransmissionNoSymptoms smaller {:d}", 0);
                return true;
            }

            if (this->get<RecoveredPerInfectedNoSymptoms>()[i] < 0.0 ||
                this->get<RecoveredPerInfectedNoSymptoms>()[i] > 1.0) {
                log_error("Constraint check: Parameter RecoveredPerInfectedNoSymptoms smaller {:d} or larger {:d}", 0,
                          1);
                return true;
            }

            if (this->get<RiskOfInfectionFromSymptomatic>()[i] < 0.0 ||
                this->get<RiskOfInfectionFromSymptomatic>()[i] > 1.0) {
                log_error("Constraint check: Parameter RiskOfInfectionFromSymptomatic smaller {:d} or larger {:d}", 0,
                          1);
                return true;
            }

            if (this->get<SeverePerInfectedSymptoms>()[i] < 0.0 || this->get<SeverePerInfectedSymptoms>()[i] > 1.0) {
                log_error("Constraint check: Parameter SeverePerInfectedSymptoms smaller {:d} or larger {:d}", 0, 1);
                return true;
            }

            if (this->get<CriticalPerSevere>()[i] < 0.0 || this->get<CriticalPerSevere>()[i] > 1.0) {
                log_error("Constraint check: Parameter CriticalPerSevere smaller {:d} or larger {:d}", 0, 1);
                return true;
            }

            if (this->get<DeathsPerCritical>()[i] < 0.0 || this->get<DeathsPerCritical>()[i] > 1.0) {
                log_error("Constraint check: Parameter DeathsPerCritical smaller {:d} or larger {:d}", 0, 1);
                return true;
            }
        }
        return false;
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
     * @see mio::deserialize
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
};

/**
 * @brief WIP !! TO DO: returns the actual, approximated reproduction rate 
 */
//double get_reprod_rate(Parameters const& params, double t, std::vector<double> const& yt);

} // namespace osecir
} // namespace mio

#endif // SECIR_PARAMETERS_H
