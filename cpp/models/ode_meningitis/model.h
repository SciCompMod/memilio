/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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
#ifndef ODEMENG_MODEL_H
#define ODEMENG_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/epidemiology/populations.h"
#include "ode_meningitis/analyze_result.h"
#include "ode_meningitis/infection_state.h"
#include "ode_meningitis/parameters.h"
#include "memilio/math/smoother.h"
#include "memilio/math/eigen_util.h"
#include "memilio/math/interpolation.h"
#include "parameters.h"

#include <numbers>
#include <complex>

namespace mio
{
namespace omeng
{

// Create template specializations for the age resolved
// SECIHURD model
// clang-format off
using Flows = TypeList<Flow<InfectionState::Incoming,               InfectionState::SusceptibleHigh>,
                       Flow<InfectionState::Incoming,               InfectionState::SusceptibleLow>,
                       Flow<InfectionState::SusceptibleHigh,        InfectionState::Carrier>,
                       Flow<InfectionState::SusceptibleLow,         InfectionState::Carrier>,
                       Flow<InfectionState::Carrier,                InfectionState::Infected>,
                       Flow<InfectionState::Infected,               InfectionState::Recovered>,
                       Flow<InfectionState::Recovered,              InfectionState::SusceptibleHigh>,
                       Flow<InfectionState::Recovered,              InfectionState::SusceptibleLow>,
                       Flow<InfectionState::Infected,               InfectionState::Dead>,
                       Flow<InfectionState::Infected,               InfectionState::DeadNatural>>;
// clang-format on

template <typename FP>
class Model : public FlowModel<FP, InfectionState, Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>
{
    using Base = FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;

    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    Model(int num_agegroups)
        : Model(Populations({mio::AgeGroup(num_agegroups), InfectionState::Count}),
                ParameterSet(mio::AgeGroup(num_agegroups)))
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        auto const& params   = this->parameters;
        AgeGroup n_agegroups = params.get_num_groups();

        ContactMatrixGroup<FP> const& contact_matrix = params.template get<ContactPatterns<FP>>();

        for (AgeGroup i = 0; i < n_agegroups; i++) {

            size_t SHi = this->populations.get_flat_index({i, InfectionState::SusceptibleHigh});
            size_t SLi = this->populations.get_flat_index({i, InfectionState::SusceptibleLow});
            size_t Ci  = this->populations.get_flat_index({i, InfectionState::Carrier});
            size_t Ii  = this->populations.get_flat_index({i, InfectionState::Infected});
            size_t Ri  = this->populations.get_flat_index({i, InfectionState::Recovered});
            size_t DDi = this->populations.get_flat_index({i, InfectionState::DeadDisease});
            size_t DNi = this->populations.get_flat_index({i, InfectionState::DeadNatural});

            FP dummy_flow_S = 0;
            for (AgeGroup j = 0; j < n_agegroups; j++) {
                size_t SHj = this->populations.get_flat_index({j, InfectionState::SusceptibleHigh});
                size_t SLj = this->populations.get_flat_index({j, InfectionState::SusceptibleLow});
                size_t Cj  = this->populations.get_flat_index({j, InfectionState::Carrier});
                size_t Ij  = this->populations.get_flat_index({j, InfectionState::Infected});
                size_t Rj  = this->populations.get_flat_index({j, InfectionState::Recovered});
                size_t DDj = this->populations.get_flat_index({j, InfectionState::DeadDisease});
                size_t DNj = this->populations.get_flat_index({j, InfectionState::DeadNatural});

                FP cont_freq_eff = contact_matrix.get_matrix_at(SimulationTime<FP>(t))(
                    static_cast<Eigen::Index>((size_t)i), static_cast<Eigen::Index>((size_t)j));

                FP Nj          = pop[SHj] + pop[SLj] + pop[Cj] + pop[Ij] + pop[Rj]; // without deceased individuals
                const FP divNj = (Nj < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nj);
                dummy_flow_S += cont_freq_eff * divNj * params.template get<TransmissionProbabilityOnContact<FP>>()[i] *
                                (params.template get<RiskOfInfectionFromFromCarrier<FP>>()[j] * pop[Cj] +
                                 params.template get<RiskOfInfectionFromFromInfected<FP>>()[j] * pop[Ij]);
            }

            // Incoming from outside to SH
            flows[this->template get_flat_flow_index<InfectionState::Incoming, InfectionState::SusceptibleHigh>({i})] +=
                (1.0 - params.template get<IncomeFractionSusLow<typename FP><FP>>()[i]) *
                params.template get<IncomeRate<FP>>()[i];

            // Incoming from outside to SL
            flows[this->template get_flat_flow_index<InfectionState::Incoming, InfectionState::SusceptibleLow>({i})] +=
                params.template get<IncomeFractionSusLow<typename FP><FP>>()[i] *
                params.template get<IncomeRate<FP>>()[i];

            // SusceptibleHigh -> Carrier
            flows[this->template get_flat_flow_index<InfectionState::SusceptibleHigh, InfectionState::Carrier>({i})] +=
                dummy_flow_S * y[SHi];

            // SusceptibleLow -> Carrier
            flows[this->template get_flat_flow_index<InfectionState::SusceptibleLow, InfectionState::Carrier>({i})] +=
                (1.0 - params.template get<ModificationRate<FP>>()[i]) * dummy_flow_S * y[SLi];

            // Carrier -> Infected
            flows[this->template get_flat_flow_index<InfectionState::Carrier, InfectionState::Infected>({i})] =
                params.template get<RateCarrierToInfected<FP>>()[i] * y[Ci];

            // Carrier -> Recovered
            flows[this->template get_flat_flow_index<InfectionState::Carrier, InfectionState::Recovered>({i})] =
                params.template get<RateCarrierToRecovered<FP>>()[i] * y[Ci];

            // Infected -> Dead
            flows[this->template get_flat_flow_index<InfectionState::Infected, InfectionState::Dead>({i})] =
                params.template get<RateInfectedToDead<FP>>()[i] * y[Ii];

            // Infected -> Recovered
            flows[this->template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>({i})] =
                params.template get<RateInfectedToRecovered<FP>>()[i] * y[Ii];

            // Recovered -> SusceptibleLow
            flows[this->template get_flat_flow_index<InfectionState::Recovered, InfectionState::SusceptibleLow>({i})] =
                params.template get<RateImmunityLoss<FP>>()[i] *
                params.template get<ProbabilityImmunityLossSusLow<FP>>()[i] * y[Ri];

            // Recovered -> SusceptibleHigh
            flows[this->template get_flat_flow_index<InfectionState::Recovered, InfectionState::SusceptibleHigh>({i})] =
                params.template get<RateImmunityLoss<FP>>()[i] *
                (1.0 - params.template get<ProbabilityImmunityLossSusLow<FP>>()[i]) * y[Ri];

            // SusceptibleHigh -> DeadNatural
            flows[this->template get_flat_flow_index<InfectionState::SusceptibleHigh, InfectionState::DeadNatural>(
                {i})] = params.template get<RateNaturalDeath<FP>>()[i] * y[SHi];

            // SusceptibleLow -> DeadNatural
            flows[this->template get_flat_flow_index<InfectionState::SusceptibleLow, InfectionState::DeadNatural>(
                {i})] = params.template get<RateNaturalDeath<FP>>()[i] * y[SLi];

            // Carrier -> DeadNatural
            flows[this->template get_flat_flow_index<InfectionState::Carrier, InfectionState::DeadNatural>({i})] =
                params.template get<RateNaturalDeath<FP>>()[i] * y[Ci];

            // Infected -> DeadNatural
            flows[this->template get_flat_flow_index<InfectionState::Infected, InfectionState::DeadNatural>({i})] =
                params.template get<RateNaturalDeath<FP>>()[i] * y[Ii];

            // Recovered -> DeadNatural
            flows[this->template get_flat_flow_index<InfectionState::Recovered, InfectionState::DeadNatural>({i})] =
                params.template get<RateNaturalDeath<FP>>()[i] * y[Ri];
        }
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Model");
        obj.add_element("Parameters", this->parameters);
        obj.add_element("Populations", this->populations);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Model> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Model");
        auto par = obj.expect_element("Parameters", Tag<ParameterSet>{});
        auto pop = obj.expect_element("Populations", Tag<Populations>{});
        return apply(
            io,
            [](auto&& par_, auto&& pop_) {
                return Model{pop_, par_};
            },
            par, pop);
    }
};

//forward declaration, see below.
template <typename FP, class BaseT = mio::Simulation<FP, Model<FP>>>
class Simulation;

/**
 * get percentage of infections per total population.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @tparam FP floating point type, e.g., double.
 * @tparam Base simulation type that uses a secir compartment model. see Simulation.
 */
template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
FP get_infections_relative(const Simulation<FP, Base>& model, FP t, const Eigen::Ref<const Eigen::VectorX<FP>>& y);

/**
 * specialization of compartment model simulation for secir models.
 * @tparam FP floating point type, e.g., double.
 * @tparam BaseT simulation type that uses a secir compartment model. default mio::Simulation. For testing purposes only!
 */
template <typename FP, class BaseT>
class Simulation : public BaseT
{
public:
    /**
     * construct a simulation.
     * @param model the model to simulate.
     * @param t0 start time
     * @param dt time steps
     */
    Simulation(mio::omeng::Model<FP> const& model, FP t0 = 0., FP dt = 0.1)
        : BaseT(model, t0, dt)
        , m_t_last_npi_check(t0)
    {
    }

    /**
     * @brief advance simulation to tmax.
     * Overwrites Simulation::advance and includes a check for dynamic NPIs in regular intervals.
     * @see Simulation::advance
     * @param tmax next stopping point of simulation
     * @return value at tmax
     */
    Eigen::Ref<Eigen::VectorX<FP>> advance(FP tmax)
    {
        using std::min;
        auto& t_end_dyn_npis   = this->get_model().parameters.get_end_dynamic_npis();
        auto& dyn_npis         = this->get_model().parameters.template get<DynamicNPIsInfectedSymptoms<FP>>();
        auto& contact_patterns = this->get_model().parameters.template get<ContactPatterns<FP>>();

        FP delay_npi_implementation; // delay which is needed to implement a NPI that is criterion-dependent
        FP t        = BaseT::get_result().get_last_time();
        const FP dt = dyn_npis.get_thresholds().size() > 0 ? FP(dyn_npis.get_interval().get()) : FP(tmax);

        while (t < tmax) {
            FP dt_eff = min<FP>(dt, tmax - t);
            dt_eff    = min<FP>(dt_eff, m_t_last_npi_check + dt - t);

            BaseT::advance(t + dt_eff);
            if (t > 0) {
                delay_npi_implementation =
                    this->get_model().parameters.template get<DynamicNPIsImplementationDelay<FP>>();
            }
            else { // DynamicNPIs for t=0 are 'misused' to be from-start NPIs. I.e., do not enforce delay.
                delay_npi_implementation = 0;
            }
            t = t + dt_eff;

            if (dyn_npis.get_thresholds().size() > 0) {
                if (floating_point_greater_equal<FP>(t, m_t_last_npi_check + dt)) {
                    if (t < t_end_dyn_npis) {
                        auto inf_rel = get_infections_relative<FP>(*this, t, this->get_result().get_last_value()) *
                                       dyn_npis.get_base_value();
                        auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
                        if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
                            (exceeded_threshold->first > m_dynamic_npi.first ||
                             t > FP(m_dynamic_npi.second))) { //old npi was weaker or is expired

                            auto t_start  = SimulationTime<FP>(t + delay_npi_implementation);
                            auto t_end    = t_start + SimulationTime<FP>(FP(dyn_npis.get_duration()));
                            m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
                            implement_dynamic_npis<FP>(contact_patterns.get_cont_freq_mat(), exceeded_threshold->second,
                                                       t_start, t_end, [](auto& g) {
                                                           return make_contact_damping_matrix(g);
                                                       });
                        }
                    }
                    m_t_last_npi_check = t;
                }
            }
            else {
                m_t_last_npi_check = t;
            }
        }

        return this->get_result().get_last_value();
    }

private:
    FP m_t_last_npi_check;
    std::pair<FP, SimulationTime<FP>> m_dynamic_npi = {-std::numeric_limits<FP>::max(), mio::SimulationTime<FP>(0)};
};

/**
 * @brief Specialization of simulate for SECIR models using Simulation.
 *
 * @tparam FP floating point type, e.g., double.
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIR model to simulate.
 * @param[in] integrator_core optional IntegratorCore, uses rk45 if nullptr.
 * 
 * @return Returns the result of the simulation.
 */
template <typename FP>
inline auto simulate(FP t0, FP tmax, FP dt, const Model<FP>& model,
                     std::unique_ptr<OdeIntegratorCore<FP>>&& integrator_core = nullptr)
{
    return mio::simulate<FP, Model<FP>, Simulation<FP>>(t0, tmax, dt, model, std::move(integrator_core));
}

/**
 * @brief Specialization of simulate for SECIR models using the FlowSimulation.
 * 
 * @tparam FP floating point type, e.g., double.
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIR model to simulate.
 * @param[in] integrator_core optional IntegratorCore, uses rk45 if nullptr.
 * 
 * @return Returns the result of the Flowsimulation.
 */
template <typename FP>
inline auto simulate_flows(FP t0, FP tmax, FP dt, const Model<FP>& model,
                           std::unique_ptr<OdeIntegratorCore<FP>>&& integrator_core = nullptr)
{
    return mio::simulate_flows<FP, Model<FP>, Simulation<FP, mio::FlowSimulation<FP, Model<FP>>>>(
        t0, tmax, dt, model, std::move(integrator_core));
}

//see declaration above.
template <typename FP, class Base>
FP get_infections_relative(const Simulation<FP, Base>& sim, FP /* t*/, const Eigen::Ref<const Eigen::VectorX<FP>>& y)
{
    FP sum_inf = 0;
    for (auto i = AgeGroup(0); i < sim.get_model().parameters.get_num_groups(); ++i) {
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::Infected});
    }
    FP inf_rel = sum_inf / sim.get_model().populations.get_total();

    return inf_rel;
}

/**
 * Get mobility factors.
 * Used by mobility graph simulation.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @return vector expression, same size as y, with mobility factors per compartment.
 * @tparam FP floating point type, e.g., double.
 * @tparam Base simulation type that uses the meningitis compartment model; see Simulation.
 */
template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
auto get_mobility_factors(const Simulation<FP, Base>& /**/, FP /*t*/, const Eigen::Ref<const Eigen::VectorX<FP>>& y)
{
    auto factors = Eigen::VectorX<FP>::Ones(y.rows()).eval();

    return factors;
}

} // namespace omeng
} // namespace mio

#endif // ODEMENG_MODEL_H
