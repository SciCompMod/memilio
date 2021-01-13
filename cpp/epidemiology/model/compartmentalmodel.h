#ifndef COMPARTMENTALMODEL_H
#define COMPARTMENTALMODEL_H

#include "epidemiology/utils/ScalarType.h"
#include "epidemiology/utils/eigen.h"
#include <vector>
#include <functional>

#define USE_DERIV_FUNC 1

namespace
{

//some metaprogramming to transform a tuple into a parameter pack and use it as
//an argument in a function.
//Taken from https://stackoverflow.com/questions/7858817/unpacking-a-tuple-to-call-a-matching-function-pointer/9288547#9288547

template <typename Function, typename Tuple, size_t... I>
decltype(auto) call(Function f, Tuple t, std::index_sequence<I...>)
{
    return f(std::get<I>(t)...);
}

template <typename Function, typename Tuple>
decltype(auto) call(Function f, Tuple t)
{
    static constexpr auto size = std::tuple_size<Tuple>::value;
    return call(f, t, std::make_index_sequence<size>{});
}

} // namespace

namespace epi
{

/**
 * @brief ComppartmentalModel is a template for a compartmental model for an
 * array of initial populations and a parameter set
 *
 * The Populations must be a concrete class derived from the Populations template,
 * i.e. a multi-dimensional array of compartment populations where each dimension
 * corresponds to a category
 *
 * The ParameterSet must be a concrete class derived form the ParameterSet template,
 * i.e. a compile-time map of parameters used by the model. These can be referenced
 * when defining the flows between compartments and they can be used for parameter
 * studies
 *
 */
template <class Populations, class ParameterSet>
struct CompartmentalModel {
public:
    // The flow function takes a set of parameters, the current time t and the
    // snapshot y of all population sizes at time t, represented as a flat array and returns a scalar value
    // that represents a flow going from one compartment to another.
    using FlowFunction =
        std::function<ScalarType(ParameterSet const& p, Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t)>;

    // A flow is a tuple of a from-index corresponding to the departing compartment, a to-index
    // corresponding to the receiving compartment and a FlowFunction. The value returned by the flow
    // function will be subtracted from the time derivative of the populations at the flat index corresponding
    // to the from-compartment, and added to the time derivative of the populations at the flat index
    // corresponding to the to-compartment.
    using Flow = std::tuple<typename Populations::Index, typename Populations::Index, FlowFunction>;

    /**
     * @brief CompartmentalModel default constructor
     */
    CompartmentalModel()
    {
    }

    /**
     * @brief add_flow defines a flow from compartment A to another compartment B
     * @param from is the index of the departure compartment A
     * @param to is the index of the receiving compartment B
     * @param f is a function defining the flow given a set of parameters, the current time t and the
     * snapshot y of all population sizes at time t, represented as a flat array
     */
    void add_flow(typename Populations::Index from, typename Populations::Index to, FlowFunction f)
    {
        flows.push_back(Flow(from, to, f));
    }

#if USE_DERIV_FUNC
    //REMARK: Not pure virtual for easier java/python bindings
    virtual void get_derivatives(Eigen::Ref<const Eigen::VectorXd>,
                                 Eigen::Ref<const Eigen::VectorXd> /*y*/, double /*t*/,
                                 Eigen::Ref<Eigen::VectorXd> /*dydt*/) const {};
#endif  // USE_DERIV_FUNC

    /**
     * @brief eval_right_hand_side evaulates the right-hand-side f of the ODE dydt = f(y, t)
     *
     * The heart of the compartmental model is a first order ODE dydt = f(y,t), where y is a flat
     * representation of all the compartmental populations at time t. This function evaluates the
     * right-hand-side f of the ODE from the intercompartmental flows. It can be used in an ODE
     * solver
     *
     * @param y the current state of the model as a flat array
     * @param t the current time
     * @param dydt a reference to the calculated output
     */
    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t, Eigen::Ref<Eigen::VectorXd> dydt) const
    {
        dydt.setZero();

#if USE_DERIV_FUNC
        this->get_derivatives(pop, y, t, dydt);
#else // USE_DERIV_FUNC
        for (auto& flow : flows) {
            ScalarType f = std::get<2>(flow)(parameters, pop, y, t);
            dydt[call(Populations::get_flat_index, std::get<0>(flow))] -= f;
            dydt[call(Populations::get_flat_index, std::get<1>(flow))] += f;
        }
#endif // USE_DERIV_FUNC
    }

    /**
     * @brief get_initial_values returns the initial values for the compartmental populations.
     * This can be used as initial conditions in an ODE solver
     * @return the initial populatoins
     */
    Eigen::VectorXd get_initial_values() const
    {
        return populations.get_compartments();
    }

    void apply_constraints()
    {
        populations.apply_constraints();
        parameters.apply_constraints();
    }

    void check_constraints() const
    {
        populations.check_constraints();
        parameters.check_constraints();
    }

    Populations populations{};
    ParameterSet parameters{};

private:
    std::vector<Flow> flows{};
};

} // namespace epi

#endif // COMPARTMENTALMODEL_H
