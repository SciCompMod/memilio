/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef MIO_COMPARTMENTS_STOCHASTIC_MODEL_H
#define MIO_COMPARTMENTS_STOCHASTIC_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/compartments/flow_model.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/utils/random_number_generator.h"

#include <cstddef>

namespace mio
{

/// @brief A CompartmentalModel with additional get_noise member.
template <typename FP, class Comp, class Pop, class Params, class Flows = void>
class StochasticModel
    : public std::conditional_t<std::is_same_v<Flows, void>, CompartmentalModel<FP, Comp, Pop, Params>,
                                FlowModel<FP, Comp, Pop, Params, Flows>>
{
public:
    using Base = std::conditional_t<std::is_same_v<Flows, void>, CompartmentalModel<FP, Comp, Pop, Params>,
                                    FlowModel<FP, Comp, Pop, Params, Flows>>;
    using Base::Base;

    /**
     * random noise part of a SDE model. represents the result of the deterministic noise matrix multiplied by a white noise vector.
     * In case of infectious disease models, the noise matrix may map noise contributions from each flow to their respective
     * compartments. In that case, the white noise vector has size #flows.
     * The resulting noise vector must have the same size as the state vector y, i.e. it must have Pop::Count entries. 
     * This is due to the fact that the integration of SDE models must happen on populations, not flows, as the applied noise
     * can occassionally push compartments into negative values. this is mitigated by the SdeIntegratorCore%s, but this
     * mitigation can not (in general) be applied to flows.
     * The noise must still be applied per flow, so both inflow and outflow use the same random value. otherwise, transitions
     * between compartments would no longer preserve the total population count. 
     */
    virtual void get_noise(Eigen::Ref<const Eigen::VectorX<FP>> /*pop*/, Eigen::Ref<const Eigen::VectorX<FP>> /*y*/,
                           FP /*t*/, Eigen::Ref<Eigen::VectorX<FP>> /*noise*/) const {};

    auto white_noise(Eigen::Index size) const
    {
        return Eigen::VectorX<FP>::NullaryExpr(size, 1, [this]() {
            return sample_standart_normal_distribution();
        });
    }

    FP sample_standart_normal_distribution() const
    {
        return FP{DistributionAdapter<std::normal_distribution<double>>::get_instance()(m_rng, 0.0, 1.0)};
    }

    RandomNumberGenerator& get_rng() const
    {
        return m_rng;
    }

private:
    mutable RandomNumberGenerator m_rng;
};

namespace details
{
template <typename FP, class M>
using get_noise_expr_t =
    decltype(std::declval<const M&>().get_noise(std::declval<Eigen::Ref<const Eigen::VectorX<FP>>>(),
                                                std::declval<Eigen::Ref<const Eigen::VectorX<FP>>>(),
                                                std::declval<FP>(), std::declval<Eigen::Ref<Eigen::VectorX<FP>>>()));
}

/**
 * Template meta function to check if a type is a valid stochastic model. 
 * Defines a static constant of name `value`. 
 * The constant `value` will be equal to true if M is a valid flow model type.
 * Otherwise, `value` will be equal to false.
 * @tparam FP A floating point type, e.g. double.
 * @tparam Model A type that may or may not be a flow model.
 */
template <typename FP, class Model>
using is_stochastic_model =
    std::integral_constant<bool, (is_expression_valid<details::get_noise_expr_t, FP, Model>::value &&
                                  is_compartment_model<FP, Model>::value)>;

} // namespace mio

#endif // MIO_COMPARTMENTS_STOCHASTIC_MODEL_H
