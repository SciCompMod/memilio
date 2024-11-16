/* 
* Copyright (C) 2020-2024 MEmilio
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

#ifndef UI10_PARAMETERS_H
#define UI10_PARAMETERS_H

#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/uncertain_value.h"

namespace mio
{
namespace oui10
{

/**************************************
 * Define Parameters of the UI model *
 **************************************/

/**
 * @brief probability of getting infected with variant 1 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV1 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV1";
    }
};

/**
 * @brief probability of getting infected with variant 2 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV2 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV2";
    }
};

/**
 * @brief probability of getting infected with variant 3 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV3 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV3";
    }
};

/**
 * @brief probability of getting infected with variant 4 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV4 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV4";
    }
};

/**
 * @brief probability of getting infected with variant 5 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV5 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV5";
    }
};

/**
 * @brief probability of getting infected with variant 6 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV6 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV6";
    }
};

/**
 * @brief probability of getting infected with variant 7 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV7 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV7";
    }
};

/**
 * @brief probability of getting infected with variant 8 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV8 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV8";
    }
};

/**
 * @brief probability of getting infected with variant 9 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV9 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV9";
    }
};

/**
 * @brief probability of getting infected with variant 10 from a contact
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContactV10 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV10";
    }
};

/**
     * @brief the infectious time of variant 1 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV1 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV1";
    }
};

/**
     * @brief the infectious time of variant 2 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV2 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV2";
    }
};

/**
     * @brief the infectious time of variant 3 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV3 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV3";
    }
};

/**
     * @brief the infectious time of variant 4 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV4 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV4";
    }
};

/**
     * @brief the infectious time of variant 5 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV5 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV5";
    }
};

/**
     * @brief the infectious time of variant 6 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV6 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV6";
    }
};

/**
     * @brief the infectious time of variant 7 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV7 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV7";
    }
};

/**
     * @brief the infectious time of variant 8 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV8 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV8";
    }
};

/**
     * @brief the infectious time of variant 9 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV9 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV9";
    }
};

/**
     * @brief the infectious time of variant 10 in day unit
     */
template <typename FP = ScalarType>
struct TimeInfectedV10 {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV10";
    }
};

/**
 * @brief the recovery patterns of variant 1 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV1 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV1";
    }
};

/**
 * @brief the recovery patterns of variant 2 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV2 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV2";
    }
};

/**
 * @brief the recovery patterns of variant 3 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV3 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV3";
    }
};

/**
 * @brief the recovery patterns of variant 4 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV4 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV4";
    }
};

/**
 * @brief the recovery patterns of variant 5 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV5 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV5";
    }
};

/**
 * @brief the recovery patterns of variant 6 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV6 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV6";
    }
};

/**
 * @brief the recovery patterns of variant 7 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV7 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV7";
    }
};

/**
 * @brief the recovery patterns of variant 8 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV8 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV8";
    }
};

/**
 * @brief the recovery patterns of variant 9 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV9 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV9";
    }
};

/**
 * @brief the recovery patterns of variant 10 are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
struct RecoveryPatternsV10 {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "RecoveryPatternsV10";
    }
};

/**
 * @brief the contact patterns within the society are modelled using a ContactMatrix
 */
template <typename FP = ScalarType>
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

template <typename FP = ScalarType>
using ParametersBase = ParameterSet<TransmissionProbabilityOnContactV1<FP>, TransmissionProbabilityOnContactV2<FP>, TransmissionProbabilityOnContactV3<FP>,
                                    TransmissionProbabilityOnContactV4<FP>, TransmissionProbabilityOnContactV5<FP>, TransmissionProbabilityOnContactV6<FP>, 
                                    TransmissionProbabilityOnContactV7<FP>, TransmissionProbabilityOnContactV8<FP>, TransmissionProbabilityOnContactV9<FP>, 
                                    TransmissionProbabilityOnContactV10<FP>, TimeInfectedV1<FP>, TimeInfectedV2<FP>, TimeInfectedV3<FP>, TimeInfectedV4<FP>,
                                    TimeInfectedV5<FP>, TimeInfectedV6<FP>, TimeInfectedV7<FP>, TimeInfectedV8<FP>, TimeInfectedV9<FP>, TimeInfectedV10<FP>,
                                    RecoveryPatternsV1<FP>, RecoveryPatternsV2<FP>, RecoveryPatternsV3<FP>, RecoveryPatternsV4<FP>, RecoveryPatternsV5<FP>, 
                                    RecoveryPatternsV6<FP>, RecoveryPatternsV7<FP>, RecoveryPatternsV8<FP>, RecoveryPatternsV9<FP>, RecoveryPatternsV10<FP>, 
                                    ContactPatterns<FP>>;

/**
 * @brief Parameters of UI model.
 */
template <typename FP = ScalarType>
class Parameters : public ParametersBase<FP>
{
public:
    Parameters(AgeGroup num_agegroups)
        : ParametersBase<FP>(num_agegroups)
        , m_num_groups{num_agegroups}
    {
    }

    AgeGroup get_num_groups()
    {
        return m_num_groups;
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
        double tol_times = 1e-1;

        int corrected = false;

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); i++) {
            if (this->template get<TimeInfectedV1<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeInfectedV1 changed from {:.4f} to {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfectedV1<FP>>()[i], tol_times);
                this->template get<TimeInfectedV1<FP>>()[i] = tol_times;
                corrected                                 = true;
            }
            if (this->template get<TimeInfectedV2<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeInfectedV2 changed from {:.4f} to {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfectedV2<FP>>()[i], tol_times);
                this->template get<TimeInfectedV2<FP>>()[i] = tol_times;
                corrected                                 = true;
            }
            if (this->template get<TimeInfectedV3<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeInfectedV3 changed from {:.4f} to {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfectedV3<FP>>()[i], tol_times);
                this->template get<TimeInfectedV3<FP>>()[i] = tol_times;
                corrected                                 = true;
            }
            if (this->template get<TransmissionProbabilityOnContactV1<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContactV1<FP>>()[i] > 1.0) {
                log_warning(
                    "Constraint check: Parameter TransmissionProbabilityOnContactV1 changed from {:0.4f} to {:d} ",
                    this->template get<TransmissionProbabilityOnContactV1<FP>>()[i], 0.0);
                this->template get<TransmissionProbabilityOnContactV1<FP>>() = 0.0;
                corrected                                                  = true;
            }
            if (this->template get<TransmissionProbabilityOnContactV2<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContactV2<FP>>()[i] > 1.0) {
                log_warning(
                    "Constraint check: Parameter TransmissionProbabilityOnContactV2 changed from {:0.4f} to {:d} ",
                    this->template get<TransmissionProbabilityOnContactV2<FP>>()[i], 0.0);
                this->template get<TransmissionProbabilityOnContactV2<FP>>() = 0.0;
                corrected                                                  = true;
            }
            if (this->template get<TransmissionProbabilityOnContactV3<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContactV3<FP>>()[i] > 1.0) {
                log_warning(
                    "Constraint check: Parameter TransmissionProbabilityOnContactV3 changed from {:0.4f} to {:d} ",
                    this->template get<TransmissionProbabilityOnContactV3<FP>>()[i], 0.0);
                this->template get<TransmissionProbabilityOnContactV3<FP>>() = 0.0;
                corrected                                                  = true;
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
        double tol_times = 1e-1;

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); i++) {

            if (this->template get<TimeInfectedV1<FP>>()[i] < tol_times) {
                log_error(
                    "Constraint check: Parameter TimeInfectedV1 {:.4f} smaller or equal {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfectedV1<FP>>()[i], 0.0);
                return true;
            }
            if (this->template get<TimeInfectedV2<FP>>()[i] < tol_times) {
                log_error(
                    "Constraint check: Parameter TimeInfectedV2 {:.4f} smaller or equal {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfectedV2<FP>>()[i], 0.0);
                return true;
            }
            if (this->template get<TimeInfectedV3<FP>>()[i] < tol_times) {
                log_error(
                    "Constraint check: Parameter TimeInfectedV3 {:.4f} smaller or equal {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfectedV3<FP>>()[i], 0.0);
                return true;
            }
            if (this->template get<TransmissionProbabilityOnContactV1<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContactV1<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContactV1 {:.4f} smaller {:.4f} or "
                          "greater {:.4f}",
                          this->template get<TransmissionProbabilityOnContactV1<FP>>()[i], 0.0, 1.0);
                return true;
            }
            if (this->template get<TransmissionProbabilityOnContactV2<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContactV2<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContactV2 {:.4f} smaller {:.4f} or "
                          "greater {:.4f}",
                          this->template get<TransmissionProbabilityOnContactV2<FP>>()[i], 0.0, 1.0);
                return true;
            }
            if (this->template get<TransmissionProbabilityOnContactV3<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContactV3<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContactV3 {:.4f} smaller {:.4f} or "
                          "greater {:.4f}",
                          this->template get<TransmissionProbabilityOnContactV3<FP>>()[i], 0.0, 1.0);
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
};

} // namespace oui10
} // namespace mio

#endif // UI_PARAMETERS_H
