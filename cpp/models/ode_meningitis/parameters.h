/* 
* Copyright (C) 2020-2026 MEmilio
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

#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/dynamic_npis.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/uncertain_value.h"
#include <limits>

namespace mio
{
namespace omeng
{

/*******************************************
 * Define Parameters of the meningitis model *
 *******************************************/

/**
 * @brief the rate with which Carrier individuals get infected.
 */
template <typename FP>
struct RateCarrierToInfected {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RateCarrierToInfected";
    }
};

/**
 * @brief the rate with which Carrier individuals get recovered.
 */
template <typename FP>
struct RateCarrierToRecovered {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RateCarrierToRecovered";
    }
};

/**
 * @brief the rate with which Carrier individuals get infected.
 */
template <typename FP>
struct RateInfectedToDead {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RateInfectedToDead";
    }
};

/**
 * @brief the rate with which Carrier individuals get recovered.
 */
template <typename FP>
struct RateInfectedToRecovered {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RateInfectedToRecovered";
    }
};

/**
 * @brief The natural death rate.
 */
template <typename FP>
struct RateNaturalDeath {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RateNaturalDeath";
    }
};

/**
 * @brief The immunity loss rate.
 */
template <typename FP>
struct RateImmunityLoss {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RateImmunityLoss";
    }
};

/**
 * @brief The probability of going to SusceptibleLow Risk upon immunity loss.
 */
template <typename FP>
struct ProbabilityImmunityLossSusLow {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "ProbabilityImmunityLossSusLow";
    }
};

/**
 * @brief probability of getting infected from a contact
 */
template <typename FP>
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
template <typename FP>
struct RiskOfInfectionFromFromCarrier {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromFromCarrier";
    }
};

/**
 * @brief the percentage of asymptomatic cases in the SECIR model
 */
template <typename FP>
struct RiskOfInfectionFromFromInfected {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromFromInfected";
    }
};

/**
 * @brief modifies the rate of SusceptibleLow to Carrier
 */
template <typename FP>
struct ModificationRate {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "ModificationRate";
    }
};

/**
 * @brief fraction of incoming individuals with low susceptibility.
 */
template <typename FP>
struct IncomeFractionSusLow {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "IncomeFractionSusLow";
    }
};

/**
 * @brief rate of incoming individuals.
 */
template <typename FP>
struct IncomeRate {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "IncomeRate";
    }
};

/**
 * @brief the contact patterns within the society are modelled using an UncertainContactMatrix
 */
template <typename FP>
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
template <typename FP>
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
template <typename FP>
struct DynamicNPIsImplementationDelay {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup /*size*/)
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "DynamicNPIsImplementationDelay";
    }
};

template <typename FP>
using ParametersBase =
    ParameterSet<RateCarrierToInfected<FP>, RateCarrierToRecovered<FP>, RateInfectedToDead<FP>,
                 RateInfectedToRecovered<FP>, RateNaturalDeath<FP>, RateImmunityLoss<FP>,
                 ProbabilityImmunityLossSusLow<FP>, TransmissionProbabilityOnContact<FP>,
                 RiskOfInfectionFromFromCarrier<FP>, RiskOfInfectionFromFromInfected<FP>,
                 ModificationRate<FP>, IncomeFractionSusLow<FP>, IncomeRate<FP>,
                 ContactPatterns<FP>, DynamicNPIsImplementationDelay<FP>, DynamicNPIsInfectedSymptoms<FP>>;

/**
 * @brief Parameters of the meningitis ODE model.
 */
template <typename FP>
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
     * Time in simulation after which no dynamic NPIs are applied.
     */
    FP& get_end_dynamic_npis()
    {
        return m_end_dynamic_npis;
    }
    FP get_end_dynamic_npis() const
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
     * @return Returns true if one or more constraints were corrected, false otherwise.
     */
    bool apply_constraints()
    {
        int corrected = false;

        if (this->template get<DynamicNPIsImplementationDelay<FP>>() < 0.0) {
            log_warning("Constraint check: Parameter DynamicNPIsImplementationDelay changed from {} to {}",
                        this->template get<DynamicNPIsImplementationDelay<FP>>(), 0);
            this->template set<DynamicNPIsImplementationDelay<FP>>(0);
            corrected = true;
        }

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {
            if (this->template get<RateCarrierToInfected<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RateCarrierToInfected changed from {} to {}",
                            this->template get<RateCarrierToInfected<FP>>()[i], 0);
                this->template get<RateCarrierToInfected<FP>>()[i] = 0.0;
                corrected                                          = true;
            }
            if (this->template get<RateCarrierToRecovered<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RateCarrierToRecovered changed from {} to {}",
                            this->template get<RateCarrierToRecovered<FP>>()[i], 0);
                this->template get<RateCarrierToRecovered<FP>>()[i] = 0.0;
                corrected                                           = true;
            }
            if (this->template get<RateInfectedToDead<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RateInfectedToDead changed from {} to {}",
                            this->template get<RateInfectedToDead<FP>>()[i], 0);
                this->template get<RateInfectedToDead<FP>>()[i] = 0.0;
                corrected                                       = true;
            }
            if (this->template get<RateInfectedToRecovered<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RateInfectedToRecovered changed from {} to {}",
                            this->template get<RateInfectedToRecovered<FP>>()[i], 0);
                this->template get<RateInfectedToRecovered<FP>>()[i] = 0.0;
                corrected                                            = true;
            }
            if (this->template get<RateNaturalDeath<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RateNaturalDeath changed from {} to {}",
                            this->template get<RateNaturalDeath<FP>>()[i], 0);
                this->template get<RateNaturalDeath<FP>>()[i] = 0.0;
                corrected                                     = true;
            }
            if (this->template get<RateImmunityLoss<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter RateImmunityLoss changed from {} to {}",
                            this->template get<RateImmunityLoss<FP>>()[i], 0);
                this->template get<RateImmunityLoss<FP>>()[i] = 0.0;
                corrected                                     = true;
            }
            if (this->template get<ProbabilityImmunityLossSusLow<FP>>()[i] < 0.0 ||
                this->template get<ProbabilityImmunityLossSusLow<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter ProbabilityImmunityLossSusLow changed from {} to {}",
                            this->template get<ProbabilityImmunityLossSusLow<FP>>()[i], 0);
                this->template get<ProbabilityImmunityLossSusLow<FP>>()[i] = 0.0;
                corrected                                                   = true;
            }
            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter TransmissionProbabilityOnContact changed from {} to {}",
                            this->template get<TransmissionProbabilityOnContact<FP>>()[i], 0);
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] = 0.0;
                corrected                                                     = true;
            }
            if (this->template get<ModificationRate<FP>>()[i] < 0.0 ||
                this->template get<ModificationRate<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter ModificationRate changed from {} to {}",
                            this->template get<ModificationRate<FP>>()[i], 0);
                this->template get<ModificationRate<FP>>()[i] = 0.0;
                corrected                                     = true;
            }
            if (this->template get<IncomeFractionSusLow<FP>>()[i] < 0.0 ||
                this->template get<IncomeFractionSusLow<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter IncomeFractionSusLow changed from {} to {}",
                            this->template get<IncomeFractionSusLow<FP>>()[i], 0);
                this->template get<IncomeFractionSusLow<FP>>()[i] = 0.0;
                corrected                                         = true;
            }
            if (this->template get<IncomeRate<FP>>()[i] < 0.0) {
                log_warning("Constraint check: Parameter IncomeRate changed from {} to {}",
                            this->template get<IncomeRate<FP>>()[i], 0);
                this->template get<IncomeRate<FP>>()[i] = 0.0;
                corrected                               = true;
            }
        }
        return corrected;
    }

    /**
     * @brief Checks whether all Parameters satisfy their constraints and logs an error if not.
     * @return Returns true if any constraint is violated, false otherwise.
     */
    bool check_constraints() const
    {
        if (this->template get<DynamicNPIsImplementationDelay<FP>>() < 0.0) {
            log_error("Constraint check: DynamicNPIsImplementationDelay smaller 0");
            return true;
        }

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {
            if (this->template get<RateCarrierToInfected<FP>>()[i] < 0.0) {
                log_error("Constraint check: RateCarrierToInfected smaller 0");
                return true;
            }
            if (this->template get<RateCarrierToRecovered<FP>>()[i] < 0.0) {
                log_error("Constraint check: RateCarrierToRecovered smaller 0");
                return true;
            }
            if (this->template get<RateInfectedToDead<FP>>()[i] < 0.0) {
                log_error("Constraint check: RateInfectedToDead smaller 0");
                return true;
            }
            if (this->template get<RateInfectedToRecovered<FP>>()[i] < 0.0) {
                log_error("Constraint check: RateInfectedToRecovered smaller 0");
                return true;
            }
            if (this->template get<RateNaturalDeath<FP>>()[i] < 0.0) {
                log_error("Constraint check: RateNaturalDeath smaller 0");
                return true;
            }
            if (this->template get<RateImmunityLoss<FP>>()[i] < 0.0) {
                log_error("Constraint check: RateImmunityLoss smaller 0");
                return true;
            }
            if (this->template get<ProbabilityImmunityLossSusLow<FP>>()[i] < 0.0 ||
                this->template get<ProbabilityImmunityLossSusLow<FP>>()[i] > 1.0) {
                log_error("Constraint check: ProbabilityImmunityLossSusLow smaller 0 or larger 1");
                return true;
            }
            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_error("Constraint check: TransmissionProbabilityOnContact smaller 0 or larger 1");
                return true;
            }
            if (this->template get<ModificationRate<FP>>()[i] < 0.0 ||
                this->template get<ModificationRate<FP>>()[i] > 1.0) {
                log_error("Constraint check: ModificationRate smaller 0 or larger 1");
                return true;
            }
            if (this->template get<IncomeFractionSusLow<FP>>()[i] < 0.0 ||
                this->template get<IncomeFractionSusLow<FP>>()[i] > 1.0) {
                log_error("Constraint check: IncomeFractionSusLow smaller 0 or larger 1");
                return true;
            }
            if (this->template get<IncomeRate<FP>>()[i] < 0.0) {
                log_error("Constraint check: IncomeRate smaller 0");
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
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Parameters> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& base, ParametersBase<FP>::deserialize(io));
        return success(Parameters(std::move(base)));
    }

private:
    AgeGroup m_num_groups;
    FP m_end_dynamic_npis = std::numeric_limits<FP>::max();
};

} // namespace omeng
} // namespace mio

#endif // SECIR_PARAMETERS_H
