
#ifndef SIRMOBILITY_PARAMETERS_H
#define SIRMOBILITY_PARAMETERS_H

#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/custom_index_array.h"
#include <iostream>
#include "regions.h"

#include <vector>

namespace mio
{
namespace osirmobility
{

/*******************************************
      * Define Parameters of the SIR model *
    *******************************************/

/**
     * @brief probability of getting infected from a contact
     */
struct TransmissionProbabilityOnContact {
    using Type = UncertainValue;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
     * @brief the infectious time in day unit
     */
struct TimeInfected {
    using Type = UncertainValue;
    static Type get_default()
    {
        return Type(6.0);
    }
    static std::string name()
    {
        return "TimeInfected";
    }
};

/**
     * @brief the contact patterns within the society are modelled using a ContactMatrix
     */
struct ContactPatterns {
    using Type = ContactMatrix;
    static Type get_default()
    {
        return Type{1};
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

/**
     * @brief The mean number of Persons migrating from one Region to another during a Time interval
     */
struct CommutingRatio {
    using Type = std::vector<std::tuple<Region, Region, double>>;
    static Type get_default()
    {
        return Type({{Region(0), Region(0), 0.}});
    }
    static std::string name()
    {
        return "CommutingRatio";
    }
};

struct ImpactCommuters {
    using Type = UncertainValue;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "ImpactCommuters";
    }
};

using ParametersBase =
    ParameterSet<TransmissionProbabilityOnContact, TimeInfected, ContactPatterns, CommutingRatio, ImpactCommuters>;

/**
 * @brief Parameters of SIR model.
 */
class Parameters : public ParametersBase
{
public:
    Parameters(Region num_regions)
        : ParametersBase() //TODO: Is this fine?
        , m_num_regions{num_regions}
    {
    }

    Region get_num_regions() const
    {
        return m_num_regions;
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
};

} // namespace osirmobility
} // namespace mio

#endif // SIR_PARAMETERS_H
