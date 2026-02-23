/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Carlotta Gerstein
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
#ifndef ODESEIRMETAPOP_PARAMETERS_H
#define ODESEIRMETAPOP_PARAMETERS_H

#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/geography/regions.h"

#include <vector>

namespace mio
{
namespace oseirmetapop
{

using Region = mio::regions::Region;

/****************************************************
* Define Parameters of the SEIR model with mobility *
****************************************************/

/**
 * @brief Probability of getting infected from a contact.
 */
template <typename FP = ScalarType>
struct TransmissionProbabilityOnContact {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(Region, AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
 * @brief The latent time in day unit.
 */
template <typename FP = ScalarType>
struct TimeExposed {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(Region, AgeGroup size)
    {
        return Type(size, 5.2);
    }
    static std::string name()
    {
        return "TimeExposed";
    }
};

/**
 * @brief The infectious time in day unit.
 */
template <typename FP = ScalarType>
struct TimeInfected {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(Region, AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfected";
    }
};

/**
 * @brief The contact patterns within the society are modelled using a 
 * ContactMatrix.
 */
template <typename FP = ScalarType>
struct ContactPatterns {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(Region, AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

/**
 * @brief The commuting patterns between different Region%s are modelled using a ContactMatrix of size n_regions x n_regions.
 * Each entry of the matrix represents the fraction of individuals commuting from one Region to another. The diagonal corresponds 
 * to the fraction of individuals staying in their Region.
 */
template <typename FP = ScalarType>
struct CommutingStrengths {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(Region size, AgeGroup)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "CommutingStrengths";
    }
};

/**
 * @brief The number of individuals in each Region and AgeGroup if commuting was applied.
 * Computed as the sum of the number of individuals staying in their Region and the number of individuals commuting to this Region 
 * minus the number of individuals commuting from this Region.
 */
template <typename FP = ScalarType>
struct PopulationAfterCommuting {
    using Type = Populations<FP, Region, AgeGroup>;
    static Type get_default(Region size_regions, AgeGroup size_agegroups)
    {
        return Type({size_regions, size_agegroups}, 0.);
    }
    static std::string name()
    {
        return "PopulationAfterCommuting";
    }
};

template <typename FP = ScalarType>
using ParametersBase = ParameterSet<TransmissionProbabilityOnContact<FP>, TimeExposed<FP>, TimeInfected<FP>,
                                    ContactPatterns<FP>, CommutingStrengths<FP>, PopulationAfterCommuting<FP>>;

/**
 * @brief Parameters of the SEIR metapopulation model.
 */
template <typename FP = ScalarType>
class Parameters : public ParametersBase<FP>
{
public:
    Parameters(Region num_regions, AgeGroup num_agegroups)
        : ParametersBase<FP>(num_regions, num_agegroups)
        , m_num_regions{num_regions}
        , m_num_agegroups(num_agegroups)
    {
    }

    Region get_num_regions() const
    {
        return m_num_regions;
    }

    AgeGroup get_num_agegroups() const
    {
        return m_num_agegroups;
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and applies them, if they do not.
     *
     * Attention: This function should be used with care. It is necessary for some test problems to run through quickly,
     *            but in a manual execution of an example, check_constraints() may be preferred. Note that the apply_constraints()
     *            function can and will not set Parameters to meaningful values in an epidemiological or virological context,
     *            as all models are designed to be transferable to multiple diseases. Consequently, only acceptable
     *            (like 0 or 1 for probabilities or small positive values for time spans) values are set here and a manual adaptation
     *            may often be necessary to have set meaningful values.
     *
     * @return Returns true if one ore more constraints were corrected, false otherwise.  
     */
    bool apply_constraints()
    {
        const FP tol_times = 1e-1;
        bool corrected     = false;

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); i++) {
            if (this->template get<TimeExposed<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeExposed changed from {} to {}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeExposed<FP>>()[i], tol_times);
                this->template get<TimeExposed<FP>>()[i] = tol_times;
                corrected                                = true;
            }
            if (this->template get<TimeInfected<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeInfected changed from {} to {}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfected<FP>>()[i], tol_times);
                this->template get<TimeInfected<FP>>()[i] = tol_times;
                corrected                                 = true;
            }
            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter TransmissionProbabilityOnContact changed from {} to {} ",
                            this->template get<TransmissionProbabilityOnContact<FP>>()[i], 0.0);
                this->template get<TransmissionProbabilityOnContact<FP>>() = 0.0;
                corrected                                                  = true;
            }
            for (auto j = Region(0); j < Region(m_num_regions); j++) {
                if (this->template get<PopulationAfterCommuting<FP>>()[{j, i}] <= 0.0) {
                    log_warning(
                        "Constraint check: Parameter PopulationAfterCommuting changed from {} to {}. Please "
                        "note that this only prevents division by zero. Consider to cancel and reset parameters.",
                        this->template get<PopulationAfterCommuting<FP>>()[{j, i}], 1.0);
                    this->template get<PopulationAfterCommuting<FP>>()[{j, i}] = 1.0;
                    corrected                                                  = true;
                }
            }
        }
        if ((this->template get<CommutingStrengths<FP>>()
                 .get_cont_freq_mat()
                 .get_matrix_at(SimulationTime<FP>(0))
                 .rowwise()
                 .sum() -
             Eigen::VectorXd::Ones((size_t)this->get_num_regions()))
                    .cwiseAbs()
                    .maxCoeff() > 1e-10 ||
            this->template get<CommutingStrengths<FP>>()
                    .get_cont_freq_mat()
                    .get_matrix_at(SimulationTime<FP>(0))
                    .minCoeff() < 0.0 ||
            this->template get<CommutingStrengths<FP>>()
                    .get_cont_freq_mat()
                    .get_matrix_at(SimulationTime<FP>(0))
                    .maxCoeff() > 1.0) {
            log_warning("Constraint check: Parameter CommutingStrengths does not ensure that the number of people "
                        "staying equals the complement of those leaving. Running without commuting.");
            this->template get<CommutingStrengths<FP>>().get_cont_freq_mat()[0].get_baseline() =
                Eigen::MatrixXd::Identity((size_t)this->get_num_regions(), (size_t)this->get_num_regions());
            corrected = true;
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
        const double tol_times = 1e-1;
        for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); i++) {
            if (this->template get<TimeExposed<FP>>()[i] < tol_times) {
                log_error(
                    "Constraint check: Parameter TimeExposed {} smaller or equal {}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeExposed<FP>>()[i], 0.0);
                return true;
            }
            if (this->template get<TimeInfected<FP>>()[i] < tol_times) {
                log_error(
                    "Constraint check: Parameter TimeInfected {} smaller or equal {}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfected<FP>>()[i], 0.0);
                return true;
            }
            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContact {} smaller {} or "
                          "greater {:.4f}",
                          this->template get<TransmissionProbabilityOnContact<FP>>()[i], 0.0, 1.0);
                return true;
            }
            for (auto j = Region(0); j < Region(m_num_regions); j++) {
                if (this->template get<PopulationAfterCommuting<FP>>()[{j, i}] <= 0.0) {
                    log_error("Constraint check: Parameter PopulationAfterCommuting {} smaller or equal {}",
                              this->template get<PopulationAfterCommuting<FP>>()[{j, i}], 0.0);
                    return true;
                }
            }
        }
        if ((this->template get<CommutingStrengths<FP>>()
                 .get_cont_freq_mat()
                 .get_matrix_at(SimulationTime<FP>(0))
                 .rowwise()
                 .sum() -
             Eigen::VectorXd::Ones((size_t)this->get_num_regions()))
                    .cwiseAbs()
                    .maxCoeff() > 1e-10 ||
            this->template get<CommutingStrengths<FP>>()
                    .get_cont_freq_mat()
                    .get_matrix_at(SimulationTime<FP>(0))
                    .minCoeff() < 0.0 ||
            this->template get<CommutingStrengths<FP>>()
                    .get_cont_freq_mat()
                    .get_matrix_at(SimulationTime<FP>(0))
                    .maxCoeff() > 1.0) {
            log_error("Constraint check: Parameter CommutingStrengths does not ensure that the number of people "
                      "staying equals the complement of those leaving.");
            return true;
        }

        return false;
    }

private:
    Parameters(ParametersBase<FP>&& base)
        : ParametersBase<FP>(std::move(base))
        , m_num_regions(base.get_num_regions())
        , m_num_agegroups(base.get_num_agegroups())
    {
    }

public:
    /**
     * @brief Deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Parameters> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& base, ParametersBase<FP>::deserialize(io));
        return success(Parameters(std::move(base)));
    }

private:
    Region m_num_regions;
    AgeGroup m_num_agegroups;
};

} // namespace oseirmetapop
} // namespace mio

#endif // SEIRMETAPOP_PARAMETERS_H
