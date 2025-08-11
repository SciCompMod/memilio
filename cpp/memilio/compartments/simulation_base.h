/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Jan Kleinert, Daniel Abele, Rene Schmieding
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
#ifndef MIO_COMPARTMENTS_SIMULATION_BASE_H
#define MIO_COMPARTMENTS_SIMULATION_BASE_H

#include "memilio/compartments/compartmental_model.h"
#include "memilio/math/integrator.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace details
{

/**
 * @brief Base class to define a Simulation.
 * Provides a protected advance method that accepts the specified integrands, that must be exposed by the Derived
 * class. Also defines all relevant members and accessors for a Simulation.
 * @tparam FP A floating point type, e.g. double.
 * @tparam M An implementation of CompartmentalModel.
 * @tparam Integrands One or more function types used for defining the right hand side of a system of equations.
 */
template <typename FP, class M, class... Integrands>
class SimulationBase
{
    static_assert(is_compartment_model<FP, M>::value, "Template parameter must be a compartment model.");

public:
    using Model = M;
    using Core  = IntegratorCore<FP, Integrands...>;

    /**
     * @brief Create a SimulationBase.
     * @param[in] model An instance of a compartmental model
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration
     */
    SimulationBase(Model const& model, std::shared_ptr<Core> integrator, FP t0, FP dt)
        : m_integratorCore(integrator)
        , m_model(std::make_unique<Model>(model))
        , m_integrator(m_integratorCore)
        , m_result(t0, m_model->get_initial_values())
        , m_dt(dt)
    {
    }

    /** 
     * @brief Constructs a copy of another SimulationBase object.
     *
     * Performs a deep copy of the model, while sharing the same integrator core.
     * The time series m_result and step size m_dt are also copied.
     *
     * @param[in] other The SimulationBase object to copy from.
     */
    SimulationBase(SimulationBase const& other)
        : m_integratorCore(other.m_integratorCore)
        , m_model(std::make_unique<Model>(*other.m_model))
        , m_integrator(m_integratorCore)
        , m_result(other.m_result)
        , m_dt(other.m_dt)
    {
    }

    /**
     * @brief Assigns another SimulationBase object to this one.
     *
     * Performs a deep copy of the model, while sharing the same integrator core.
     * The time series m_result and step size m_dt are also copied.
     *
     * @param[in] other The SimulationBase object to assign from.
     * @return Reference to this SimulationBase object.
     */
    SimulationBase& operator=(SimulationBase const& other)
    {
        if (this != &other) {
            m_integratorCore = other.m_integratorCore;
            m_model          = std::make_unique<Model>(*other.m_model);
            m_integrator.set_integrator(m_integratorCore);
            m_result = other.m_result;
            m_dt     = other.m_dt;
        }
        return *this;
    }

    /**
     * @brief Set the integrator core used in the simulation.
     * @param[in] integrator A shared pointer to an object derived from IntegratorCore.
     */
    void set_integrator(std::shared_ptr<Core> integrator)
    {
        m_integratorCore = std::move(integrator);
        m_integrator.set_integrator(m_integratorCore);
    }

    /**
     * @brief Access the integrator core used in the simulation.
     * @return A reference to the integrator core used in the simulation
     * @{
     */
    Core& get_integrator()
    {
        return *m_integratorCore;
    }
    const Core& get_integrator() const
    {
        return *m_integratorCore;
    }
    /** @} */

    /**
     * @brief Returns the simulation result describing the model population in each time step.
     *
     * Which compartments are used by the model is defined by the Comp template argument for the CompartmentalModel
     * (usually an enum named InfectionState).
     *
     * @return A TimeSeries to represent a numerical solution for the population of the model.
     * For each simulated time step, the TimeSeries contains the population size in each compartment.
     * @{
     */
    TimeSeries<FP>& get_result()
    {
        return m_result;
    }
    const TimeSeries<FP>& get_result() const
    {
        return m_result;
    }
    /** @} */

    /**
     * @brief Get a reference to the model owned and used by the simulation.
     * @return The simulation model.
     * @{
     */
    const Model& get_model() const
    {
        return *m_model;
    }
    Model& get_model()
    {
        return *m_model;
    }
    /** @} */

    /**
     * @brief Returns the step size used by the integrator.
     * When using a integration scheme with adaptive time stepping, the integrator will store its estimate for the
     * next step size in this value.
     * @{
     */
    FP& get_dt()
    {
        return m_dt;
    }
    const FP& get_dt() const
    {
        return m_dt;
    }
    /** @} */

protected:
    /**
     * @brief Run the simulation up to a given time.
     * @param[in] fs Integrands passed to the integrator, e.g. a wrapper for `get_derivatives`.
     * @param[in] tmax Next stopping point of the simulation.
     * @param[in,out] results The TimeSeries to use as initial value and for storing integration results. 
     * @return The simulation result at tmax.
     */
    Eigen::Ref<Eigen::VectorX<FP>> advance(const Integrands&... fs, FP tmax, TimeSeries<FP>& results)
    {
        return m_integrator.advance(fs..., tmax, m_dt, results);
    }

private:
    std::shared_ptr<Core> m_integratorCore; ///< Defines the integration scheme via its step function.
    std::unique_ptr<Model> m_model; ///< The model defining the ODE system and initial conditions.
    SystemIntegrator<FP, Integrands...>
        m_integrator; ///< Integrates the DerivFunction (see advance) and stores resutls in m_result.
    TimeSeries<FP> m_result; ///< The simulation results.
    FP m_dt; ///< The time step used (and possibly set) by m_integratorCore::step.
};

/// @brief Specialization of SimulationBase that takes a SystemIntegrator instead of it's Integrands.
template <typename FP, class M, class... Integrands>
class SimulationBase<FP, M, SystemIntegrator<FP, Integrands...>> : public SimulationBase<FP, M, Integrands...>
{
    using SimulationBase<FP, M, Integrands...>::SimulationBase;
};

} // namespace details
} // namespace mio

#endif // MIO_COMPARTMENTS_SIMULATION_BASE_H
