/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Annika Jungklaus, Lena Ploetzke
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

#ifndef LCT_SECIR_2_DISEASES_PARAMS_H
#define LCT_SECIR_2_DISEASES_PARAMS_H

#include "memilio/config.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/uncertain_matrix.h"

namespace mio
{
namespace lsecir2d
{

/*********************************************************
* Define Parameters of the LCT-SECIHURD-2-DISEASES model *
**********************************************************/

/**
 * @brief Average time spent in the Exposed compartment for disease a.
 */
struct TimeExposed_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeExposed_a";
    }
};

/**
 * @brief Average time spent in the TimeInfectedNoSymptoms before developing 
 *  symptoms or recover for disease a in day unit.
 */
struct TimeInfectedNoSymptoms_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptoms_a";
    }
};

/**
 * @brief Average time spent in the TimeInfectedSymptoms before going to hospital 
 *  or recover for disease a in day unit.
 */
struct TimeInfectedSymptoms_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptoms_a";
    }
};

/**
 * @brief Average time being in the Hospital before treated by ICU or recover for disease a in day unit.
 */
struct TimeInfectedSevere_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedSevere_a";
    }
};

/**
 * @brief Average time treated by ICU before dead or recover for disease a in day unit.
 */
struct TimeInfectedCritical_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedCritical_a";
    }
};

/**
 * @brief Probability of getting infected from a contact for disease a.
 */
struct TransmissionProbabilityOnContact_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact_a";
    }
};

/**
 * @brief Average time spent in the Exposed compartment for disease b.
 */
struct TimeExposed_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeExposed_b";
    }
};

/**
 * @brief Average time spent in the TimeInfectedNoSymptoms before developing 
 *  symptoms or recover for disease b in day unit.
 */
struct TimeInfectedNoSymptoms_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptoms_b";
    }
};

/**
 * @brief Average time spent in the TimeInfectedSymptoms before going to hospital 
 *  or recover for disease b in day unit.
 */
struct TimeInfectedSymptoms_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptoms_b";
    }
};

/**
 * @brief Average time being in the Hospital before treated by ICU or recover for disease b in day unit.
 */
struct TimeInfectedSevere_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedSevere_b";
    }
};

/**
 * @brief Average time treated by ICU before dead or recover for disease b in day unit.
 */
struct TimeInfectedCritical_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeInfectedCritical_b";
    }
};

/**
 * @brief Probability of getting infected from a contact for disease b.
 */
struct TransmissionProbabilityOnContact_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact_b";
    }
};

/**
 * @brief The contact patterns within the society are modelled using an UncertainContactMatrix.
 */
struct ContactPatterns {
    using Type = UncertainContactMatrix<ScalarType>;
    static Type get_default(size_t size = 1) // no age groups
    {
        mio::ContactMatrixGroup contact_matrix(1, (Eigen::Index)size);
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant((Eigen::Index)size, (Eigen::Index)size, 10.));
        return Type(contact_matrix);
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

/**
 * @brief The relative InfectedNoSymptoms infectability for disease a.
 */
struct RelativeTransmissionNoSymptoms_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "RelativeTransmissionNoSymptoms_a";
    }
};

/**
 * @brief The risk of infection from symptomatic cases for disease a.
 */
struct RiskOfInfectionFromSymptomatic_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic_a";
    }
};

/**
 * @brief The relative InfectedNoSymptoms infectability for disease b.
 */
struct RelativeTransmissionNoSymptoms_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "RelativeTransmissionNoSymptoms_b";
    }
};

/**
 * @brief The risk of infection from symptomatic cases for disease b.
 */
struct RiskOfInfectionFromSymptomatic_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic_b";
    }
};

/**
 * @brief The percentage of asymptomatic cases for disease a.
 */
struct RecoveredPerInfectedNoSymptoms_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 0.5);
    }
    static std::string name()
    {
        return "RecoveredPerInfectedNoSymptoms_a";
    }
};

/**
 * @brief The percentage of hospitalized patients per infected patients for disease a.
 */
struct SeverePerInfectedSymptoms_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 0.5);
    }
    static std::string name()
    {
        return "SeverePerInfectedSymptoms_a";
    }
};

/**
 * @brief The percentage of ICU patients per hospitalized patients for disease a.
 */
struct CriticalPerSevere_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 0.5);
    }
    static std::string name()
    {
        return "CriticalPerSevere_a";
    }
};

/**
 * @brief The percentage of dead patients per ICU patients for disease a.
 */
struct DeathsPerCritical_a {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 0.1);
    }
    static std::string name()
    {
        return "DeathsPerCritical_a";
    }
};

/**
 * @brief The percentage of asymptomatic cases for disease b.
 */
struct RecoveredPerInfectedNoSymptoms_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 0.5);
    }
    static std::string name()
    {
        return "RecoveredPerInfectedNoSymptoms_b";
    }
};

/**
 * @brief The percentage of hospitalized patients per infected patients for disease b.
 */
struct SeverePerInfectedSymptoms_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 0.5);
    }
    static std::string name()
    {
        return "SeverePerInfectedSymptoms_b";
    }
};

/**
 * @brief The percentage of ICU patients per hospitalized patients for disease b.
 */
struct CriticalPerSevere_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 0.5);
    }
    static std::string name()
    {
        return "CriticalPerSevere_b";
    }
};

/**
 * @brief The percentage of dead patients per ICU patients for disease b.
 */
struct DeathsPerCritical_b {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size = 1)
    {
        return Type::Constant(size, 1, 0.1);
    }
    static std::string name()
    {
        return "DeathsPerCritical_b";
    }
};

/**
 * @brief The start day in the LCT-SECIR-2-DISEASES  model.
 * The start day defines in which season the simulation is started.
 * If the start day is 180 and simulation takes place from t0=0 to
 * tmax=100 the days 180 to 280 of the year are simulated.
 */
struct StartDay {
    using Type = ScalarType;
    static Type get_default(size_t)
    {
        return 0.;
    }
    static std::string name()
    {
        return "StartDay";
    }
};

/**
 * @brief The seasonality in the LCT-SECIR-2-DISEASES model.
 * The seasonality is given as (1+k*sin()) where the sine
 * curve is below one in summer and above one in winter.
 */
struct Seasonality {
    using Type = ScalarType;
    static Type get_default(size_t)
    {
        return 0.;
    }
    static std::string name()
    {
        return "Seasonality";
    }
};

using ParametersBase =
    ParameterSet<TimeExposed_a, TimeInfectedNoSymptoms_a, TimeInfectedSymptoms_a, TimeInfectedSevere_a,
                 TimeInfectedCritical_a, TimeExposed_b, TimeInfectedNoSymptoms_b, TimeInfectedSymptoms_b,
                 TimeInfectedSevere_b, TimeInfectedCritical_b, TransmissionProbabilityOnContact_a,
                 TransmissionProbabilityOnContact_b, ContactPatterns, RelativeTransmissionNoSymptoms_a,
                 RiskOfInfectionFromSymptomatic_a, RecoveredPerInfectedNoSymptoms_a, SeverePerInfectedSymptoms_a,
                 CriticalPerSevere_a, DeathsPerCritical_a, RelativeTransmissionNoSymptoms_b,
                 RiskOfInfectionFromSymptomatic_b, RecoveredPerInfectedNoSymptoms_b, SeverePerInfectedSymptoms_b,
                 CriticalPerSevere_b, DeathsPerCritical_b, StartDay, Seasonality>;

/**
 * @brief Parameters of an LCT-SECIR-2-DISEASES model.
 */
class Parameters : public ParametersBase
{
public:
    /**
     * @brief Constructor.
     * @param num_groups The number of groups considered in the LCT2D model.
     */
    Parameters(size_t num_groups)
        : ParametersBase(num_groups)
        , m_num_groups{num_groups}
    {
    }

    size_t get_num_groups() const
    {
        return m_num_groups;
    }

    /**
     * @brief Checks whether all parameters satisfy their corresponding constraints and throws errors, if they do not.
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false. 
     */
    bool check_constraints() const
    {
        for (size_t i = 0; i < m_num_groups; ++i) {
            if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
                log_warning("Constraint check: Parameter Seasonality should lie between {:0.4f} and {:.4f}", 0.0, 0.5);
                return true;
            }

            if (this->get<TimeExposed_a>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeExposed_a is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TimeExposed_b>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeExposed_b is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TimeInfectedNoSymptoms_a>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedNoSymptoms_a is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TimeInfectedNoSymptoms_b>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedNoSymptoms_b is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TimeInfectedSymptoms_a>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedSymptoms_a is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TimeInfectedSymptoms_b>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedSymptoms_b is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TimeInfectedSevere_a>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedSevere_a is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TimeInfectedSevere_b>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedSevere_b is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TimeInfectedCritical_a>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedCritical_a is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TimeInfectedCritical_b>()[i] < 1.0) {
                log_error("Constraint check: Parameter TimeInfectedCritical_b is smaller than {:.4f}", 1.0);
                return true;
            }

            if (this->get<TransmissionProbabilityOnContact_a>()[i] < 0.0 ||
                this->get<TransmissionProbabilityOnContact_a>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContact_a smaller {:d} or larger {:d}",
                          0, 1);
                return true;
            }

            if (this->get<TransmissionProbabilityOnContact_b>()[i] < 0.0 ||
                this->get<TransmissionProbabilityOnContact_b>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContact_b smaller {:d} or larger {:d}",
                          0, 1);
                return true;
            }

            if (this->get<RelativeTransmissionNoSymptoms_a>()[i] < 0.0 ||
                this->get<RelativeTransmissionNoSymptoms_a>()[i] > 1.0) {
                log_error("Constraint check: Parameter RelativeTransmissionNoSymptoms_a smaller {:d} or larger {:d}", 0,
                          1);
                return true;
            }

            if (this->get<RelativeTransmissionNoSymptoms_b>()[i] < 0.0 ||
                this->get<RelativeTransmissionNoSymptoms_b>()[i] > 1.0) {
                log_error("Constraint check: Parameter RelativeTransmissionNoSymptoms_b smaller {:d} or larger {:d}", 0,
                          1);
                return true;
            }

            if (this->get<RiskOfInfectionFromSymptomatic_a>()[i] < 0.0 ||
                this->get<RiskOfInfectionFromSymptomatic_a>()[i] > 1.0) {
                log_error("Constraint check: Parameter  RiskOfInfectionFromSymptomatic_a smaller {:d} or larger {:d}",
                          0, 1);
                return true;
            }

            if (this->get<RiskOfInfectionFromSymptomatic_b>()[i] < 0.0 ||
                this->get<RiskOfInfectionFromSymptomatic_b>()[i] > 1.0) {
                log_error("Constraint check: Parameter  RiskOfInfectionFromSymptomatic_b smaller {:d} or larger {:d}",
                          0, 1);
                return true;
            }

            if (this->get<RecoveredPerInfectedNoSymptoms_a>()[i] < 0.0 ||
                this->get<RecoveredPerInfectedNoSymptoms_a>()[i] > 1.0) {
                log_error("Constraint check: Parameter RecoveredPerInfectedNoSymptoms_a smaller {:d} or larger {:d}", 0,
                          1);
                return true;
            }

            if (this->get<RecoveredPerInfectedNoSymptoms_b>()[i] < 0.0 ||
                this->get<RecoveredPerInfectedNoSymptoms_b>()[i] > 1.0) {
                log_error("Constraint check: Parameter RecoveredPerInfectedNoSymptoms_b smaller {:d} or larger {:d}", 0,
                          1);
                return true;
            }

            if (this->get<SeverePerInfectedSymptoms_a>()[i] < 0.0 ||
                this->get<SeverePerInfectedSymptoms_a>()[i] > 1.0) {
                log_error("Constraint check: Parameter SeverePerInfectedSymptoms_a smaller {:d} or larger {:d}", 0, 1);
                return true;
            }

            if (this->get<SeverePerInfectedSymptoms_b>()[i] < 0.0 ||
                this->get<SeverePerInfectedSymptoms_b>()[i] > 1.0) {
                log_error("Constraint check: Parameter SeverePerInfectedSymptoms_b smaller {:d} or larger {:d}", 0, 1);
                return true;
            }

            if (this->get<CriticalPerSevere_a>()[i] < 0.0 || this->get<CriticalPerSevere_a>()[i] > 1.0) {
                log_error("Constraint check: Parameter CriticalPerSevere_a smaller {:d} or larger {:d}", 0, 1);
                return true;
            }

            if (this->get<CriticalPerSevere_b>()[i] < 0.0 || this->get<CriticalPerSevere_b>()[i] > 1.0) {
                log_error("Constraint check: Parameter CriticalPerSevere_b smaller {:d} or larger {:d}", 0, 1);
                return true;
            }

            if (this->get<DeathsPerCritical_a>()[i] < 0.0 || this->get<DeathsPerCritical_a>()[i] > 1.0) {
                log_error("Constraint check: Parameter DeathsPerCritical_a smaller {:d} or larger {:d}", 0, 1);
                return true;
            }

            if (this->get<DeathsPerCritical_b>()[i] < 0.0 || this->get<DeathsPerCritical_b>()[i] > 1.0) {
                log_error("Constraint check: Parameter DeathsPerCritical_b smaller {:d} or larger {:d}", 0, 1);
                return true;
            }
        }

        return false;
    }

private:
    Parameters(ParametersBase&& base)
        : ParametersBase(std::move(base))
        , m_num_groups(this->template get<ContactPatterns>().get_cont_freq_mat().get_num_groups())
    {
    }

    size_t m_num_groups;

public:
    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Parameters> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& base, ParametersBase::deserialize(io));
        return success(Parameters(std::move(base)));
    }
};

} // namespace lsecir2d
} // namespace mio

#endif // LCT_SECIR_2_DISEASES_PARAMS_H
