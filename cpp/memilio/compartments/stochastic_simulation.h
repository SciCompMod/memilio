#ifndef MIO_STOCHASTIC_SIMULATION_H
#define MIO_STOCHASTIC_SIMULATION_H

// #include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/sde_model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler_maruyama.h"
#include "memilio/utils/time_series.h"

namespace mio
{

template <typename FP, class M>
class StochasticSimulation : public SimulationBase<FP, M, 2>
{
    static_assert(is_stochsatic_model<FP, M>::value, "Can only be used with StochasticModel.");

public:
    using Base  = SimulationBase<FP, M, 2>;
    using Model = M;

    /**
     * @brief Setup the simulation with an ODE solver.
     * @param[in] model An instance of a compartmental model
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration
     */
    StochasticSimulation(Model const& model, FP t0 = 0., FP dt = 0.1)
        : Base(model, std::make_shared<EulerMaruyamaIntegratorCore<FP>>(), t0, dt)
    {
    }

    /**
     * @brief advance simulation to tmax
     * tmax must be greater than get_result().get_last_time_point()
     * @param tmax next stopping point of simulation
     */
    Eigen::Ref<Eigen::VectorX<FP>> advance(FP tmax)
    {
        return Base::advance({[this](auto&& y, auto&& t, auto&& dydt) {
                                  Base::get_model().eval_right_hand_side(y, y, t, dydt);
                              },
                              [this](auto&& y, auto&& t, auto&& dydt) {
                                  dydt.setZero();
                                  Base::get_model().get_noise(y, y, t, dydt);
                              }},
                             tmax, Base::get_result());
    }
};

// FIXME: a flow simulation does not make sense, since the rescaler can only reasonably act on population.
// applying the result from rescaling to flows is - as always - generally impossible.
// We can, however, allow flow models - though their benefit is limited.

// template <typename FP, class M>
// class StochasticFlowSimulation : public FlowSimulationBase<FP, M, 2>
// {
//     // TODO: assert M is flow and has get_noise
// public:
//     using Base  = FlowSimulationBase<FP, M, 2>;
//     using Model = M;

//     /**
//      * @brief Setup the simulation with an ODE solver.
//      * @param[in] model An instance of a compartmental model
//      * @param[in] t0 Start time.
//      * @param[in] dt Initial step size of integration
//      */
//     StochasticFlowSimulation(Model const& model, FP t0 = 0., FP dt = 0.1)
//         : Base(model, std::make_shared<EulerMaruyamaIntegratorCore<FP>>(), t0, dt)
//         , m_t_pop_eval(t0 - FP{1.})
//         , m_pop(model.get_initial_values().size())
//     {
//     }

//     /**
//      * @brief advance simulation to tmax
//      * tmax must be greater than get_result().get_last_time_point()
//      * @param tmax next stopping point of simulation
//      */
//     Eigen::Ref<Eigen::VectorX<FP>> advance(FP tmax)
//     {
//         const auto result = Base::advance({[this](auto&& flows, auto&& t, auto&& dflows_dt) {
//                                                get_pop(flows, t);
//                                                dflows_dt.setZero();
//                                                this->get_model().get_flows(m_pop, m_pop, t, dflows_dt);
//                                            },
//                                            [this](auto&& flows, auto&& t, auto&& dnoise_dt) {
//                                                get_pop(flows, t);
//                                                dnoise_dt.setZero();
//                                                this->get_model().get_noise(m_pop, m_pop, t, dnoise_dt);
//                                            }},
//                                           tmax, Base::get_flows());
//         Base::compute_population_results();
//         return result;
//     }

// private:
//     void get_pop(Eigen::Ref<const Eigen::VectorX<FP>> flows, FP t)
//     {
//         const auto& pop_result = this->get_result();

//         if (m_t_pop_eval != t) {
//             this->get_model().get_derivatives(flows - Base::get_flows().get_value(pop_result.get_num_time_points() - 1),
//                                               m_pop);
//             m_pop += pop_result.get_last_value();
//             m_t_pop_eval = t;
//         }
//     }

//     FP m_t_pop_eval; ///< time at which m_pop was evaluated at, used to avoid recalculating the same population
//     Eigen::VectorX<FP> m_pop; ///< pre-allocated temporary, used in right_hand_side()
// };

template <typename FP, class Model, class Sim = StochasticSimulation<FP, Model>>
TimeSeries<FP> simulate_stochastic(FP t0, FP tmax, FP dt, Model const& model,
                                   std::shared_ptr<IntegratorCore<FP, 2>> integrator = nullptr)
{
    model.check_constraints();
    Sim sim(model, t0, dt);
    if (integrator) {
        sim.set_integrator(integrator);
    }
    sim.advance(tmax);
    return sim.get_result();
}

// template <typename FP, class Model, class Sim = StochasticFlowSimulation<FP, Model>>
// std::vector<TimeSeries<FP>> simulate_flows_stochastic(FP t0, FP tmax, FP dt, Model const& model,
//                                                       std::shared_ptr<IntegratorCore<FP, 2>> integrator = nullptr)
// {
//     model.check_constraints();
//     Sim sim(model, t0, dt);
//     if (integrator) {
//         sim.set_integrator(integrator);
//     }
//     sim.advance(tmax);
//     return {sim.get_result(), sim.get_flows()};
// }

} // namespace mio

#endif // MIO_STOCHASTIC_SIMULATION_H
