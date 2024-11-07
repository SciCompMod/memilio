
#ifndef SEIRMOBILITY_PARAMETERS_H
#define SEIRMOBILITY_PARAMETERS_H

#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/custom_index_array.h"
#include "models/ode_seir_mobility_improved/regions.h"
#include "Eigen/Sparse"

#include <vector>

namespace mio
{
namespace oseirmobilityimproved
{

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
 * @brief the latent time in day unit
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
     * @brief The contact patterns within the society are modelled using a ContactMatrix.
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
 * @brief The contact patterns between different Region%s are modelled using a ContactMatrix.
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
 * @brief The sizes of the populations after commuting.
 */
template <typename FP = ScalarType>
struct PopulationSizes {
    using Type = CustomIndexArray<FP, Region>;
    static Type get_default(Region size, AgeGroup)
    {
        return Type(size, 0.);
    }
    static std::string name()
    {
        return "PopulationSizes";
    }
};

template <typename FP = ScalarType>
using ParametersBase = ParameterSet<TransmissionProbabilityOnContact<FP>, TimeExposed<FP>, TimeInfected<FP>,
                                    ContactPatterns<FP>, CommutingStrengths<FP>, PopulationSizes<FP>>;

/**
 * @brief Parameters of SEIR model.
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

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); i++) {
            if (this->template get<TimeExposed<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeInfected changed from {:.4f} to {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeExposed<FP>>()[i], tol_times);
                this->template get<TimeExposed<FP>>()[i] = tol_times;
                corrected                                = true;
            }
            if (this->template get<TimeInfected<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeInfected changed from {:.4f} to {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfected<FP>>()[i], tol_times);
                this->template get<TimeInfected<FP>>()[i] = tol_times;
                corrected                                 = true;
            }
            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_warning(
                    "Constraint check: Parameter TransmissionProbabilityOnContact changed from {:0.4f} to {:d} ",
                    this->template get<TransmissionProbabilityOnContact<FP>>()[i], 0.0);
                this->template get<TransmissionProbabilityOnContact<FP>>() = 0.0;
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

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); i++) {
            if (this->template get<TimeExposed<FP>>()[i] < tol_times) {
                log_error(
                    "Constraint check: Parameter TimeExposed {:.4f} smaller or equal {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeExposed<FP>>()[i], 0.0);
                return true;
            }
            if (this->template get<TimeInfected<FP>>()[i] < tol_times) {
                log_error(
                    "Constraint check: Parameter TimeInfected {:.4f} smaller or equal {:.4f}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfected<FP>>()[i], 0.0);
                return true;
            }
            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContact {:.4f} smaller {:.4f} or "
                          "greater {:.4f}",
                          this->template get<TransmissionProbabilityOnContact<FP>>()[i], 0.0, 1.0);
                return true;
            }
        }
        return false;
    }

private:
    // Parameters(ParametersBase&& base)
    //     : ParametersBase(std::move(base)) //TODO: Adjust
    // {
    // }

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
    Region m_num_regions;
    AgeGroup m_num_agegroups;
};

} // namespace oseirmobilityimproved
} // namespace mio

#endif // SEIR_PARAMETERS_H
