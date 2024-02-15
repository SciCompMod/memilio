/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC),
* Forschungszentrum Juelich (FZJ)
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn, Ralf Hannemann-Tamas
*
* Contact: Ralf Hannemann-Tamas <r.hannemann-tamas@fz-juelich.de>
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
#ifndef SEAIR_PARAMETERS_H
#define SEAIR_PARAMETERS_H

//#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/parameter_set.h"


namespace mio
{
namespace oseair
{

/*******************************************
      * Define Parameters of the SEAIR model *
    *******************************************/

/**
 * @brief Social distancing.
 */
template <typename FP = double>
struct AlphaA {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.2);
    }
    static std::string name()
    {
        return "AlphaA";
    }
};

/**
 * @brief Quarantining.
 */
template <typename FP = double>
struct AlphaI {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.2);
    }
    static std::string name()
    {
        return "AlphaI";
    }
};

/**
 * @brief Rate of testing.
 */
template <typename FP = double>
struct Kappa {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.2);
    }
    static std::string name()
    {
        return "Kappa";
    }
};

/**
 * @brief Recovery rate.
 */
template <typename FP = double>
struct Beta {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.0067);
    }
    static std::string name()
    {
        return "Beta";
    }
};

/**
 * @brief Death Rate.
 */
template <typename FP = double>
struct Mu {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.0041);
    }
    static std::string name()
    {
        return "Mu";
    }
};

/**
 * @brief Inverse of the latent period of the virus.
 */
template <typename FP = double>
struct TLatentInverse {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.5);
    }
    static std::string name()
    {
        return "TLatentInverse";
    }
};

/**
 * @brief Infectious period for unconfirmed infected people.
 */
template <typename FP = double>
struct Rho {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.1);
    }
    static std::string name()
    {
        return "Rho";
    }
};

/**
 * @brief Rate recovered people become susceptible again.
 */
template <typename FP = double>
struct Gamma {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "Gamma";
    }
};



template <typename FP = double>
using ParametersBase =
    ParameterSet<AlphaA<FP>, AlphaI<FP>, Kappa<FP>, Beta<FP>, Mu<FP>, TLatentInverse<FP>, Rho<FP>, Gamma<FP>>;

/**
 * @brief Parameters of an age-resolved SECIR/SECIHURD model.
 */
template <typename FP = double>
class Parameters : public ParametersBase<FP>
{
public:
    Parameters()
        : ParametersBase<FP>()
    {
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and logs an error 
     * if constraints are not satisfied.
     * @return Returns 1 if one constraint is not satisfied, otherwise 0.   
     */
    int check_constraints() const
    {
        return 0;
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
        BOOST_OUTCOME_TRY(base, ParametersBase<FP>::deserialize(io));
        return success(Parameters(std::move(base)));
    }
};

} // namespace oseair
} // namespace mio

#endif // SEAIR_PARAMETERS_H
