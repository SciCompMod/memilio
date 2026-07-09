/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Annika Jungklaus
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

#ifndef MIO_GLCT_VECTOR_PARAMS_H
#define MIO_GLCT_VECTOR_PARAMS_H

#include "memilio/utils/parameter_set.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/logging.h"

namespace mio
{
namespace glvector
{



/***********************************************
* Define Parameters of the GLCT-SEIR-VECTOR model *
***********************************************/

/// @brief Vector with the probability to start in any of the subcompartments of the ExposedHuman compartment.
template <typename FP>
struct StartingProbabilitiesExposedHuman {
    using Type = Eigen::VectorX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in the ExposedHuman compartment.
     * @param[in] numExposedHuman Number of subcompartments of the ExposedHuman compartment.
     */
    static Type get_default(size_t numExposedHuman)
    {
        Eigen::VectorX<FP> def = Eigen::VectorX<FP>::Zero(numExposedHuman);
        def[0]                 = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesExposedHuman";
    }
};

/// @brief Vector with the probability to start in any of the subcompartments of the ExposedVector compartment.
template <typename FP>
struct StartingProbabilitiesExposedVector {
    using Type = Eigen::VectorX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in the ExposedVector compartment.
     * @param[in] numExposedVector Number of subcompartments of the ExposedVector compartment.
     */
    static Type get_default(size_t numExposedVector)
    {
        Eigen::VectorX<FP> def = Eigen::VectorX<FP>::Zero(numExposedVector);
        def[0]                 = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesExposedVector";
    }
};

/// @brief Transition matrix of the ExposedHuman compartment.
template <typename FP>
struct TransitionMatrixExposedToInfectedHuman {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in the ExposedHuman compartment.
     * @param[in] numExposedHuman Number of subcompartments of the ExposedHuman compartment.
     * @param[in] timeExposedHuman Average time spent in ExposedHuman compartment in day unit.
     */
    static Type get_default(size_t numExposedHuman, FP timeExposedHuman = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(numExposedHuman, -(FP)numExposedHuman / timeExposedHuman).asDiagonal();
        def.diagonal(1).setConstant((FP)numExposedHuman / timeExposedHuman);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixExposedToInfectedHuman";
    }
};

/// @brief Transition matrix of the ExposedVector compartment.
template <typename FP>
struct TransitionMatrixExposedToTransmitterVector {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in the ExposedVector compartment.
     * @param[in] numExposedVector Number of subcompartments of the ExposedVector compartment.
     * @param[in] timeExposedVector Average time spent in ExposedVector compartment in day unit.
     */
    static Type get_default(size_t numExposedVector, FP timeExposedVector = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(numExposedVector, -(FP)numExposedVector / timeExposedVector).asDiagonal();
        def.diagonal(1).setConstant((FP)numExposedVector / timeExposedVector);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixExposedToTransmitterVector";
    }
};

/// @brief Vector with the probability to start in any of the subcompartments of the InfectedHuman compartment.
template <typename FP>
struct StartingProbabilitiesInfectedHuman {
    using Type = Eigen::VectorX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedHuman compartment.
     * @param[in] numInfectedHuman Number of subcompartments of the InfectedHuman compartment.
     */
    static Type get_default(size_t numInfectedHuman)
    {
        Eigen::VectorX<FP> def = Eigen::VectorX<FP>::Zero(numInfectedHuman);
        def[0]                 = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedHuman";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedHuman
 *      compartment before recovering.
 */
template <typename FP>
struct TransitionMatrixInfectedToRecoveredHuman {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedHuman compartment
     *   before recovering.
     * @param[in] numInfectedHuman Number of rows/columns of the transition matrix.
     * @param[in] timeInfectedBeforeRecoveredHuman Average time spent in InfectedHuman before recovering in day unit.
     */
    static Type get_default(size_t numInfectedHuman, FP timeInfectedBeforeRecoveredHuman = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(numInfectedHuman, -(FP)numInfectedHuman / timeInfectedBeforeRecoveredHuman).asDiagonal();
        def.diagonal(1).setConstant((FP)numInfectedHuman / timeInfectedBeforeRecoveredHuman);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedToRecoveredHuman";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedHuman
 *      compartment before dying.
 */
template <typename FP>
struct TransitionMatrixInfectedToDeadHuman {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedHuman compartment
     *   before recovery.
     * @param[in] dimension Number of rows/columns of the transition matrix.
     * @param[in] time Average time spent in InfectedHuman before recovery in day unit.
     */
    static Type get_default(size_t dimension, FP time = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(dimension, -(FP)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((FP)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedToDeadHuman";
    }
};


/// @brief Probability of a person getting infected from a contact with an infectious vector.
template <typename FP>
struct TransmissionProbabilityOnContactVectorToHuman {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactVectorToHuman";
    }
};

/// @brief Probability of a vector getting infected from a contact with an infectious person.
template <typename FP>
struct TransmissionProbabilityOnContactHumanToVector {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactHumanToVector";
    }
};

/// @brief Birth rate per day for the human population.
template <typename FP>
struct BirthRateHuman {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "BirthRateHuman";
    }
};

/// @brief Mortality rate per day for the human population, only counting non-disease related deaths.
template <typename FP>
struct MortalityRateHuman {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "MortalityRateHuman";
    }
};

/// @brief Mortality rate per day for the aquatic stage in the vector population.
template <typename FP>
struct MortalityRateAquaticVector {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "MortalityRateAquaticVector";
    }
};

/// @brief Mortality rate per day for the adult stages (Susceptible, Exposed, Transmitter) in the vector population.
template <typename FP>
struct MortalityRateAdultVector {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "MortalityRateAdultVector";
    }
};

/// @brief Biting Rate (number of humans a vector bites per day).
template <typename FP>
struct BitingRateVector {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "BitingRateVector";
    }
};


/// @brief Oviposition rate per day for the vector population.
template <typename FP>
struct OvipositionRateVector {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "OvipositionRateVector";
    }
};

/// @brief Aquatic-to-Adult transition rate per day for the vector population.
template <typename FP>
struct TransitionRateAquaticToAdultVector {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "TransitionRateAquaticToAdultVector";
    }
};

/**
* @brief Proportion of the vector population that is female (for e.g. mosquitoes only females bite). 
*       If all individuals can carry the disease, the sex ratio can be set to 1.
*/ 
template <typename FP>
struct SexRatioVector {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "SexRatioVector";
    }
};


/// @brief Mean vector-to-human ratio.
template <typename FP>
struct RatioVectorToHuman {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "RatioVectorToHuman";
    }
};



/// @brief Duration of the seasonal cycle in days (e.g. 1 year = 365).
template <typename FP>
struct SeasonalCycleDuration {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "SeasonalCycleDuration";
    }
};

/**
*  @brief Proportion of variation of the carrying capacity for the vector species. 
*       0 equals no variation, so no seasonality, while 1 equals maximum oscillations around the mean.
*/
template <typename FP>
struct SeasonalVariationImpact {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "SeasonalVariationImpact";
    }
};


template <typename FP>
using ParametersBase =
    ParameterSet<StartingProbabilitiesExposedHuman<FP>, StartingProbabilitiesExposedVector<FP>, TransitionMatrixExposedToInfectedHuman<FP>,
                 TransitionMatrixExposedToTransmitterVector<FP>, StartingProbabilitiesInfectedHuman<FP>, TransitionMatrixInfectedToDeadHuman<FP>,
                 TransitionMatrixInfectedToRecoveredHuman<FP>, TransmissionProbabilityOnContactVectorToHuman<FP>, 
                 TransmissionProbabilityOnContactHumanToVector<FP>, BirthRateHuman<FP>, MortalityRateHuman<FP>, MortalityRateAquaticVector<FP>,
                 MortalityRateAdultVector<FP>, TransitionRateAquaticToAdultVector<FP>, OvipositionRateVector<FP>, 
                 SexRatioVector<FP>, RatioVectorToHuman<FP>, SeasonalCycleDuration<FP>, SeasonalVariationImpact<FP>, BitingRateVector<FP>>;

/// @brief Parameters of an GLCT-VECTOR model.
template <typename FP>
class Parameters : public ParametersBase<FP>
{
public:
    /// @brief Default constructor.
    Parameters()
        : ParametersBase<FP>()
    {
    }

    /**
     * @brief Checks that all parameters satisfy their corresponding constraints and logs an error
     *      if constraints are not satisfied.
     *
     * @return Returns true if one or more constraints are not satisfied, false otherwise.
     */
    bool check_constraints() const
    {   

        // --- Parameters affecting population growth/shrink not connected to the disease. ---
        if (this->template get<BirthRateHuman<FP>>() < 0.0) {
            log_error("Constraint check: Parameter BirthRateHuman smaller {}", 0);
            return true;
        }

        if (this->template get<MortalityRateHuman<FP>>() < 0.0) {
            log_error("Constraint check: Parameter MortalityRateHuman smaller {}", 0);
            return true;
        }

        if (this->template get<MortalityRateAquaticVector<FP>>() < 0.0) {
            log_error("Constraint check: Parameter MortalityRateAquaticVector smaller {}", 0);
            return true;
        }

        if (this->template get<MortalityRateAdultVector<FP>>() < 0.0) {
            log_error("Constraint check: Parameter MortalityRateAdultVector smaller {}", 0);
            return true;
        }
        
        if (this->template get<OvipositionRateVector<FP>>() < 0.0) {
            log_error("Constraint check: Parameter OvipostionRateVector smaller {}", 0);
            return true;
        }

        if (this->template get<TransitionRateAquaticToAdultVector<FP>>() < 0.0) {
            log_error("Constraint check: Parameter TransitionRateToAdultVector smaller {}", 0);
            return true;
        }

        if (this->template get<RatioVectorToHuman<FP>>() < 0.0) {
            log_error("Constraint check: Parameter RatioVectorToHuman smaller {}", 0);
            return true;
        }

        if (this->template get<SexRatioVector<FP>>() < 0.0 || this->template get<SexRatioVector<FP>>() > 1.0) {
            log_error("Constraint check: Parameter SexRatioVector smaller {} or larger {}", 0, 1);
            return true;
        }

        if (this->template get<SeasonalCycleDuration<FP>>() < 1e-10) {
            log_error("Constraint check: Parameter SeasonalCycleDuration smaller or euqal to {}", 0);
            return true;
        }

        if (this->template get<SeasonalVariationImpact<FP>>() < 0.0 || this->template get<SeasonalVariationImpact<FP>>() > 1.0) {
            log_error("Constraint check: Parameter SeasonalVariationImpact smaller {} or larger {}", 0, 1);
            return true;
        }


        // --- Parameters affecting the transmission of the virus. ---
        if (this->template get<BitingRateVector<FP>>() < 0.0) {
            log_error("Constraint check: Parameter BitingRateVector smaller {}", 0);
            return true;
        }

        if (this->template get<TransmissionProbabilityOnContactVectorToHuman<FP>>() < 0.0 ||
            this->template get<TransmissionProbabilityOnContactVectorToHuman<FP>>() > 1.0) {
            log_error("Constraint check: Parameter TransmissionProbabilityOnContactVectorToHuman smaller {} or larger {}", 0, 1);
            return true;
        }

        if (this->template get<TransmissionProbabilityOnContactHumanToVector<FP>>() < 0.0 ||
        this->template get<TransmissionProbabilityOnContactHumanToVector<FP>>() > 1.0) {
        log_error("Constraint check: Parameter TransmissionProbabilityOnContactHumanToVector smaller {} or larger {}", 0, 1);
        return true;
        }        

        // --- Parameters affecting the phase-type distributions. ---
        // --- Check that the dimensions are consistent. ---
        if ((this->template get<TransitionMatrixExposedToInfectedHuman<FP>>().cols() !=
             this->template get<TransitionMatrixExposedToInfectedHuman<FP>>().rows()) ||
            (this->template get<TransitionMatrixExposedToTransmitterVector<FP>>().cols() !=
             this->template get<TransitionMatrixExposedToTransmitterVector<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedToDeadHuman<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedToDeadHuman<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedToRecoveredHuman<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedToRecoveredHuman<FP>>().rows())) {
            log_error("Constraint check: At least one of the matrices used for the TransitionMatrix parameters is not "
                      "quadratic.");
            return true;
        }

        if (this->template get<StartingProbabilitiesExposedHuman<FP>>().rows() !=
            this->template get<TransitionMatrixExposedToInfectedHuman<FP>>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesExposedHuman and "
                      "TransitionMatrixExposedToInfectedHuman are not matching.");
            return true;
        }

        if (this->template get<StartingProbabilitiesInfectedHuman<FP>>().rows() !=
            this->template get<TransitionMatrixInfectedToDeadHuman<FP>>().rows() +
                this->template get<TransitionMatrixInfectedToRecoveredHuman<FP>>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesInfectedHuman and "
                      "TransitionMatrices of InfectedHuman compartment are not matching.");
            return true;
        }

        if (this->template get<StartingProbabilitiesExposedVector<FP>>().rows() !=
            this->template get<TransitionMatrixExposedToTransmitterVector<FP>>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesExposedVector and "
                      "TransitionMatrixExposedVector are not matching.");
            return true;
        }

        // --- Check constraints of the starting probability vectors. ---
        if ((!floating_point_equal<FP>(1., this->template get<StartingProbabilitiesExposedHuman<FP>>().sum())) ||
            (!floating_point_equal<FP>(1., this->template get<StartingProbabilitiesInfectedHuman<FP>>().sum())) ||
            (!floating_point_equal<FP>(1., this->template get<StartingProbabilitiesExposedVector<FP>>().sum()))) {
            log_warning(
                "Constraint check: At least one of the vectors for the starting probabilities does not sum to one.");
            return true;
        }

        if ((this->template get<StartingProbabilitiesExposedHuman<FP>>().array() < -1e-10).any() ||
            (this->template get<StartingProbabilitiesInfectedHuman<FP>>().array() < -1e-10).any() ||
            (this->template get<StartingProbabilitiesExposedVector<FP>>().array() < -1e-10).any()) {
            log_warning("Constraint check: At least one of the vectors for the starting probabilities has at least one "
                        "negative entry.");
            return true;
        }

        // --- Check that we have no flows back from one compartment to the previous one
        // (only in between of the subcompartments). ---
        if (((this->template get<TransitionMatrixExposedToInfectedHuman<FP>>() *
              Eigen::VectorX<FP>::Ones(this->template get<TransitionMatrixExposedToInfectedHuman<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixExposedToInfectedHuman lead to a negative "
                "flow ExposedToInfectedHuman.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedToDeadHuman<FP>>() *
              Eigen::VectorX<FP>::Ones(
                  this->template get<TransitionMatrixInfectedToDeadHuman<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning("Constraint check: The entries of TransitionMatrixInfectedToDeadHuman lead to "
                        "a negative "
                        "flow InfectedToDeadHuman.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedToRecoveredHuman<FP>>() *
              Eigen::VectorX<FP>::Ones(this->template get<TransitionMatrixInfectedToRecoveredHuman<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedToRecoveredHuman lead to a negative "
                "flow InfectedToRecoveredHuman.");
            return true;
        }
        if (((this->template get<TransitionMatrixExposedToTransmitterVector<FP>>() *
              Eigen::VectorX<FP>::Ones(
                  this->template get<TransitionMatrixExposedToTransmitterVector<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixExposedToTransmitterVector lead to a negative "
                "flow InfectedSymptomsToInfectedSevere.");
            return true;
        }
        return false;
    }

private:
    Parameters(ParametersBase<FP>&& base)
        : ParametersBase<FP>(std::move(base))
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
};

} // namespace glvector
} // namespace mio

#endif // MIO_GLCT_VECTOR_PARAMS_H
