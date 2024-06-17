
#ifndef SEIRMOBILITY_PARAMETERS_H
#define SEIRMOBILITY_PARAMETERS_H

#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/custom_index_array.h"
#include "models/ode_seir_mobility/regions.h"

#include <vector>

namespace mio
{
namespace oseirmobility
{

/****************************************************
 * Define Parameters of the SEIR model with mobility *
 ****************************************************/

/**
 * @brief Probability of getting infected from a contact.
 */
struct TransmissionProbabilityOnContact {
    using Type = UncertainValue;
    static Type get_default(Region)
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
     * @brief the latent time in day unit
     */
struct TimeExposed {
    using Type = UncertainValue;
    static Type get_default()
    {
        return Type(5.2);
    }
    static std::string name()
    {
        return "TimeExposed";
    }
};

/**
     * @brief The infectious time in day unit.
     */
struct TimeInfected {
    using Type = UncertainValue;
    static Type get_default(Region)
    {
        return Type(6.0);
    }
    static std::string name()
    {
        return "TimeInfected";
    }
};

/**
     * @brief The contact patterns within the society are modelled using a ContactMatrix.
     */
struct ContactPatterns {
    using Type = ContactMatrix;
    static Type get_default(Region)
    {
        return Type{1};
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

/**
 * @brief The mean number of people migrating from one Region to another during a TimeStep.
 */
struct CommutingRatio {
    using Type = std::vector<std::tuple<Region, Region, double>>;
    static Type get_default(Region)
    {
        return Type({{Region(0), Region(0), 0.}});
    }
    static std::string name()
    {
        return "CommutingRatio";
    }
};

/**
 * @brief The ratio that regulates the infections during commuting.
*/
struct ImpactCommuters {
    using Type = UncertainValue;
    static Type get_default(Region)
    {
        return Type(0.);
    }
    static std::string name()
    {
        return "ImpactCommuters";
    }
};

/**
 * @brief The Region%s that a person crosses when travelling from one Region to another. 
*/
struct PathIntersections {
    using Type = CustomIndexArray<std::vector<Region>, Region, Region>;
    static Type get_default(Region size)
    {
        return Type({size, size});
    }
    static std::string name()
    {
        return "PathIntersections";
    }
};

using ParametersBase = ParameterSet<TransmissionProbabilityOnContact, TimeExposed, TimeInfected, ContactPatterns,
                                    CommutingRatio, ImpactCommuters, PathIntersections>;

/**
 * @brief Parameters of SEIR model.
 */
class Parameters : public ParametersBase
{
public:
    Parameters(Region num_regions, AgeGroup num_agegroups)
        : ParametersBase(num_regions)
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
        if (this->get<TimeExposed>() < tol_times) {
            log_warning("Constraint check: Parameter TimeExposed changed from {:.4f} to {:.4f}. Please note that "
                        "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                        "and reset parameters.",
                        this->get<TimeExposed>(), tol_times);
            this->get<TimeExposed>() = tol_times;
            corrected                = true;
        }
        if (this->get<TimeInfected>() < tol_times) {
            log_warning("Constraint check: Parameter TimeInfected changed from {:.4f} to {:.4f}. Please note that "
                        "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                        "and reset parameters.",
                        this->get<TimeInfected>(), tol_times);
            this->get<TimeInfected>() = tol_times;
            corrected                 = true;
        }
        if (this->get<TransmissionProbabilityOnContact>() < 0.0 ||
            this->get<TransmissionProbabilityOnContact>() > 1.0) {
            log_warning("Constraint check: Parameter TransmissionProbabilityOnContact changed from {:0.4f} to {:d} ",
                        this->get<TransmissionProbabilityOnContact>(), 0.0);
            this->get<TransmissionProbabilityOnContact>() = 0.0;
            corrected                                     = true;
        }
        if (this->get<ImpactCommuters>() < 0.0 || this->get<ImpactCommuters>() > 1.0) {
            log_warning("Constraint check: Parameter ImpactCommuters changed from {:.4f} to {:.4f}.",
                        this->get<ImpactCommuters>(), 0.0);
            this->get<ImpactCommuters>() = 0.0;
            corrected                    = true;
        }
        for (auto& i : this->get<CommutingRatio>()) {
            if (std::get<double>(i) < 0.0 || std::get<double>(i) > 1.0) {
                log_warning("Constraint check: Parameter CommutingRatio changed from {:.4f} to {:.4f}.",
                            std::get<double>(i), 0.0);
                std::get<double>(i) = 0.0;
                corrected           = true;
            }
            if (std::get<0>(i) < Region(0) || std::get<1>(i) < Region(0) || std::get<0>(i) >= m_num_regions ||
                std::get<1>(i) >= m_num_regions) {
                log_warning(
                    "Constraint check: Removed entry of Parameter CommutingRatio because of non-existing Regions.");
                auto it = std::find(this->get<CommutingRatio>().begin(), this->get<CommutingRatio>().end(), i);
                this->get<CommutingRatio>().erase(it);
                corrected = true;
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

        if (this->get<TimeExposed>() < tol_times) {
            log_error("Constraint check: Parameter TimeExposed {:.4f} smaller or equal {:.4f}. Please note that "
                      "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                      "and reset parameters.",
                      this->get<TimeExposed>(), 0.0);
            return true;
        }
        if (this->get<TimeInfected>() < tol_times) {
            log_error("Constraint check: Parameter TimeInfected {:.4f} smaller or equal {:.4f}. Please note that "
                      "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                      "and reset parameters.",
                      this->get<TimeInfected>(), 0.0);
            return true;
        }
        if (this->get<TransmissionProbabilityOnContact>() < 0.0 ||
            this->get<TransmissionProbabilityOnContact>() > 1.0) {
            log_error(
                "Constraint check: Parameter TransmissionProbabilityOnContact {:.4f} smaller {:.4f} or greater {:.4f}",
                this->get<TransmissionProbabilityOnContact>(), 0.0, 1.0);
            return true;
        }
        if (this->get<ImpactCommuters>() < 0.0 || this->get<ImpactCommuters>() > 1.0) {
            log_error("Constraint check: Parameter ImpactCommuters {:.4f} smaller {:.4f} or greater {:.4f}",
                      this->get<ImpactCommuters>(), 0.0, 1.0);
            return true;
        }
        for (auto i : this->get<CommutingRatio>()) {
            if (std::get<double>(i) < 0.0 || std::get<double>(i) > 1.0) {
                log_error("Constraint check: Parameter CommutingRatio entry {:.4f} smaller {:.4f} or greater {:.4f}",
                          std::get<double>(i), 0.0, 1.0);
                return true;
            }
            if (std::get<0>(i) < Region(0) || std::get<1>(i) < Region(0) || std::get<0>(i) > m_num_regions ||
                std::get<1>(i) > m_num_regions) {
                log_error("Constraint check: Parameter CommutingRatio has an entry with start or end Region that does "
                          "not appear in the model.");
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
        BOOST_OUTCOME_TRY(base, ParametersBase::deserialize(io));
        return success(Parameters(std::move(base)));
    }

private:
    Region m_num_regions;
    AgeGroup m_num_agegroups;
};

} // namespace oseirmobility
} // namespace mio

#endif // SEIR_PARAMETERS_H
