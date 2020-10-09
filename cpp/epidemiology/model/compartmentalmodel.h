#ifndef COMPARTMENTALMODEL_H
#define COMPARTMENTALMODEL_H

#include "epidemiology/model/ScalarType.h"
#include <vector>
#include <functional>
#include <Eigen/Core>

namespace
{

//some metaprogramming to transform a tuple into a parameter pack and use it as
//an argument in a function.
//Taken from https://stackoverflow.com/questions/7858817/unpacking-a-tuple-to-call-a-matching-function-pointer/9288547#9288547

template <typename Function, typename Tuple, size_t... I>
auto call(Function f, Tuple t, std::index_sequence<I...>)
{
    return f(std::get<I>(t)...);
}

template <typename Function, typename Tuple>
auto call(Function f, Tuple t)
{
    static constexpr auto size = std::tuple_size<Tuple>::value;
    return call(f, t, std::make_index_sequence<size>{});
}

} // namespace

namespace epi
{

template <class Populations, class ParameterSet>
struct CompartmentalModel {
public:
    using FlowFunction =
        std::function<ScalarType(ParameterSet const& p, Eigen::Ref<const Eigen::VectorXd> y, double t)>;
    using Flow = std::tuple<typename Populations::Index, typename Populations::Index, FlowFunction>;

    CompartmentalModel()
    {
    }

    void add_flow(typename Populations::Index from, typename Populations::Index to, FlowFunction f)
    {
        flows.push_back(Flow(from, to, f));
    }

    /**
     * @brief get_right_hand_side returns the right-hand-side f of the ODE dydt = f(y, t)
     * @param y the current state of the model as a flat array
     * @param t the current time
     * @param dydt a reference to the calculated output
     */
    void get_right_hand_side(Eigen::Ref<const Eigen::VectorXd> y, double t, Eigen::Ref<Eigen::VectorXd> dydt) const
    {
        for (size_t i = 0; i < y.size(); ++i) {
            dydt[i] = 0;
        }
        for (auto flow : flows) {
            ScalarType f = std::get<2>(flow)(parameters, y, t);
            dydt[call(Populations::get_flat_index, std::get<0>(flow))] += f;
            dydt[call(Populations::get_flat_index, std::get<1>(flow))] -= f;
        }
    }

    Eigen::VectorXd get_initial_values() const
    {
        return populations.get_compartments();
    }

    Populations populations{};
    ParameterSet parameters{};

private:
    std::vector<Flow> flows{};
};

} // namespace epi

#endif // COMPARTMENTALMODEL_H
