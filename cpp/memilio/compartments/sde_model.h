#pragma once

#include "memilio/compartments/compartmentalmodel.h"

namespace mio
{

template <typename FP, class BaseModel>
class StochasticModel : public BaseModel
{
public:
    static_assert(is_compartment_model<FP, BaseModel>::value, "StochasticModel needs to be a CompartmentalModel.");
    using Base = BaseModel;
    using Base::Base;

    virtual void get_noise(Eigen::Ref<const Eigen::VectorX<FP>> /*pop*/, Eigen::Ref<const Eigen::VectorX<FP>> /*y*/,
                           FP /*t*/, Eigen::Ref<Eigen::VectorX<FP>> /*noise*/) const {};
};

namespace details
{
template <typename FP, class M>
using get_noise_expr_t =
    decltype(std::declval<const M&>().get_flows(std::declval<Eigen::Ref<const Eigen::VectorX<FP>>>(),
                                                std::declval<Eigen::Ref<const Eigen::VectorX<FP>>>(),
                                                std::declval<FP>(), std::declval<Eigen::Ref<Eigen::VectorX<FP>>>()));
}

/**
 * Template meta function to check if a type is a valid flow model. 
 * Defines a static constant of name `value`. 
 * The constant `value` will be equal to true if M is a valid flow model type.
 * Otherwise, `value` will be equal to false.
 * @tparam FP floating point type, e.g. double.
 * @tparam Model A type that may or may not be a flow model.
 */
template <typename FP, class Model>
using is_stochsatic_model =
    std::integral_constant<bool, (is_expression_valid<details::get_noise_expr_t, FP, Model>::value &&
                                  is_compartment_model<FP, Model>::value)>;

} // namespace mio
