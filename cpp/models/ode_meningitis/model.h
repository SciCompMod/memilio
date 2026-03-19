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
                       Flow<InfectionState::Recovered,              InfectionState::InfectedHigh>,
                       Flow<InfectionState::Recovered,              InfectionState::InfectedLow>,
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

                FP Nj           = pop[SHj] + pop[SLj] + pop[Cj] + pop[Ij] + pop[Rj]; // without deceased individuals
                const FP divNj  = (Nj < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nj);
                FP dummy_flow_S = cont_freq_eff * divNj *
                                  params.template get<TransmissionProbabilityOnContact<FP>>()[i] *
                                  (params.template get<RiskOfInfectionFromFromCarrier<FP>>()[j] * pop[INSj] +
                                   params.template get<RiskOfInfectionFromFromInfected<FP>>()[j] * pop[ISyj]);
            }

            // Incoming from outside to SL

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
                params.template get<RateInfectedToDead<FP>>()[i] * y[Ci];

            // Infected -> Recovered
            flows[this->template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>({i})] =
                params.template get<RateInfectedToRecovered<FP>>()[i] * y[Ci];

            // Recovered -> SusceptibleLow
            flows[this->template get_flat_flow_index<InfectionState::Recovered, InfectionState::SusceptibleLow>({i})] =
                params.template get<RateImmunityLoss<FP>>()[i] *
                params.template get<ProbabilityImmunityLossSusLow<FP>>()[i] * y[Ci];

            // Recovered -> SusceptibleHigh
            flows[this->template get_flat_flow_index<InfectionState::Recovered, InfectionState::SusceptibleHigh>({i})] =
                params.template get<RateImmunityLoss<FP>>()[i] *
                (1.0 - params.template get<ProbabilityImmunityLossSusLow<FP>>()[i]) * y[Ci];
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
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedSymptoms});
    }
    FP inf_rel = sum_inf / sim.get_model().populations.get_total();

    return inf_rel;
}

/**
*@brief Computes the reproduction number at a given index time of the Model output obtained by the Simulation.
*@param t_idx The index time at which the reproduction number is computed.
*@param sim The simulation holding the SECIR model
*@tparam Base simulation type that uses a SECIR compartment model. see Simulation.
*@returns The computed reproduction number at the provided index time.
*/

template <typename FP, class Base>
IOResult<FP> get_reproduction_number(size_t t_idx, const Simulation<FP, Base>& sim)
{
    using std::abs;
    using std::sin;

    if (!(t_idx < static_cast<size_t>(sim.get_result().get_num_time_points()))) {
        return mio::failure(mio::StatusCode::OutOfRange, "t_idx is not a valid index for the TimeSeries");
    }

    auto const& params      = sim.get_model().parameters;
    const size_t num_groups = (size_t)sim.get_model().parameters.get_num_groups();
    //The infected compartments in the SECIR Model are: Exposed, Carrier, Infected, Hospitalized, ICU in respective agegroups
    const size_t num_infected_compartments   = 5;
    const size_t total_infected_compartments = num_infected_compartments * num_groups;
    const FP pi                              = std::numbers::pi_v<ScalarType>;

    //F encodes new Infections and V encodes transition times in the next-generation matrix calculation of R_t
    Eigen::MatrixX<FP> F(total_infected_compartments, total_infected_compartments);
    Eigen::MatrixX<FP> V(total_infected_compartments, total_infected_compartments);
    F = Eigen::MatrixX<FP>::Zero(total_infected_compartments,
                                 total_infected_compartments); //Initialize matrices F and V with zeroes
    V = Eigen::MatrixX<FP>::Zero(total_infected_compartments, total_infected_compartments);

    FP test_and_trace_required = 0.0;
    FP icu_occupancy           = 0.0;
    for (auto i = AgeGroup(0); i < (mio::AgeGroup)num_groups; ++i) {
        test_and_trace_required +=
            (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
            params.template get<TimeInfectedNoSymptoms<FP>>()[i] *
            sim.get_result().get_value(
                t_idx)[sim.get_model().populations.get_flat_index({i, InfectionState::InfectedNoSymptoms})];
        icu_occupancy += sim.get_result().get_value(
            t_idx)[sim.get_model().populations.get_flat_index({i, InfectionState::InfectedCritical})];
    }

    FP season_val =
        (1 +
         params.template get<Seasonality<FP>>() *
             sin(pi *
                 ((sim.get_model().parameters.template get<StartDay<FP>>() + sim.get_result().get_time(t_idx)) / 182.5 +
                  0.5)));
    ContactMatrixGroup<FP> const& contact_matrix = sim.get_model().parameters.template get<ContactPatterns<FP>>();

    Eigen::MatrixX<FP> cont_freq_eff(num_groups, num_groups);
    Eigen::MatrixX<FP> riskFromInfectedSymptomatic_derivatives(num_groups, num_groups);
    Eigen::VectorX<FP> divN(num_groups);
    Eigen::VectorX<FP> riskFromInfectedSymptomatic(num_groups);

    for (mio::AgeGroup k = 0; k < (mio::AgeGroup)num_groups; k++) {
        FP temp = sim.get_result().get_value(
                      t_idx)[sim.get_model().populations.get_flat_index({k, InfectionState::Susceptible})] +
                  sim.get_result().get_value(
                      t_idx)[sim.get_model().populations.get_flat_index({k, InfectionState::Exposed})] +
                  sim.get_result().get_value(
                      t_idx)[sim.get_model().populations.get_flat_index({k, InfectionState::InfectedNoSymptoms})] +
                  sim.get_result().get_value(
                      t_idx)[sim.get_model().populations.get_flat_index({k, InfectionState::InfectedSymptoms})] +
                  sim.get_result().get_value(
                      t_idx)[sim.get_model().populations.get_flat_index({k, InfectionState::InfectedSevere})] +
                  sim.get_result().get_value(
                      t_idx)[sim.get_model().populations.get_flat_index({k, InfectionState::InfectedCritical})] +
                  sim.get_result().get_value(
                      t_idx)[sim.get_model().populations.get_flat_index({k, InfectionState::Recovered})];
        if (temp < Limits<FP>::zero_tolerance()) {
            temp = 1.0;
        }
        divN[(size_t)k] = 1 / temp;

        riskFromInfectedSymptomatic[(size_t)k] = smoother_cosine<FP>(
            test_and_trace_required, params.template get<TestAndTraceCapacity<FP>>(),
            (params.template get<TestAndTraceCapacity<FP>>()) * params.template get<TestAndTraceCapacityMaxRisk<FP>>(),
            params.template get<RiskOfInfectionFromSymptomatic<FP>>()[k],
            params.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[k]);

        for (mio::AgeGroup l = 0; l < (mio::AgeGroup)num_groups; l++) {
            if (test_and_trace_required < params.template get<TestAndTraceCapacity<FP>>() ||
                test_and_trace_required > params.template get<TestAndTraceCapacityMaxRisk<FP>>() *
                                              params.template get<TestAndTraceCapacity<FP>>()) {
                riskFromInfectedSymptomatic_derivatives((size_t)k, (size_t)l) = 0;
            }
            else {
                riskFromInfectedSymptomatic_derivatives((size_t)k, (size_t)l) =
                    -0.5 * pi *
                    (params.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[k] -
                     params.template get<RiskOfInfectionFromSymptomatic<FP>>()[k]) /
                    (4 * params.template get<TestAndTraceCapacity<FP>>()) *
                    (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[l]) /
                    params.template get<TimeInfectedNoSymptoms<FP>>()[l] *
                    sin(pi / (4 * params.template get<TestAndTraceCapacity<FP>>()) *
                        (test_and_trace_required - params.template get<TestAndTraceCapacity<FP>>()));
            }
        }

        for (Eigen::Index l = 0; l < (Eigen::Index)num_groups; l++) {
            cont_freq_eff(l, (size_t)k) =
                season_val * contact_matrix.get_matrix_at(SimulationTime<FP>(static_cast<FP>(t_idx)))(
                                 static_cast<Eigen::Index>((size_t)l), static_cast<Eigen::Index>((size_t)k));
        }
    }

    //Check criterion if matrix V will be invertible by checking if subblock J is invertible
    Eigen::MatrixX<FP> J(num_groups, num_groups);
    J = Eigen::MatrixX<FP>::Zero(num_groups, num_groups);
    for (size_t i = 0; i < num_groups; i++) {
        J(i, i) = 1 / (params.template get<TimeInfectedCritical<FP>>()[(mio::AgeGroup)i]);

        if (!(icu_occupancy < 0.9 * params.template get<ICUCapacity<FP>>() ||
              icu_occupancy > params.template get<ICUCapacity<FP>>())) {
            for (size_t j = 0; j < num_groups; j++) {
                J(i, j) -= sim.get_result().get_value(t_idx)[sim.get_model().populations.get_flat_index(
                               {(mio::AgeGroup)i, InfectionState::InfectedSevere})] /
                           params.template get<TimeInfectedSevere<FP>>()[(mio::AgeGroup)i] * 5 *
                           params.template get<CriticalPerSevere<FP>>()[(mio::AgeGroup)i] * pi /
                           (params.template get<ICUCapacity<FP>>()) *
                           sin(pi / (0.1 * params.template get<ICUCapacity<FP>>()) *
                               (icu_occupancy - 0.9 * params.template get<ICUCapacity<FP>>()));
            }
        }
    }

    //Check, if J is invertible
    if (J.determinant() == 0) {
        return mio::failure(mio::StatusCode::UnknownError, "Matrix V is not invertible");
    }

    //Initialize the matrix F
    for (size_t i = 0; i < num_groups; i++) {
        for (size_t j = 0; j < num_groups; j++) {

            FP temp = 0;
            for (Eigen::Index k = 0; k < (Eigen::Index)num_groups; k++) {
                temp += cont_freq_eff(i, k) *
                        sim.get_result().get_value(t_idx)[sim.get_model().populations.get_flat_index(
                            {(mio::AgeGroup)k, InfectionState::InfectedSymptoms})] *
                        riskFromInfectedSymptomatic_derivatives(k, j) * divN[k];
            }

            F(i, j + num_groups) =
                sim.get_result().get_value(t_idx)[sim.get_model().populations.get_flat_index(
                    {(mio::AgeGroup)i, InfectionState::Susceptible})] *
                params.template get<TransmissionProbabilityOnContact<FP>>()[(mio::AgeGroup)i] *
                (cont_freq_eff(i, j) * params.template get<RelativeTransmissionNoSymptoms<FP>>()[(mio::AgeGroup)j] *
                     divN[(size_t)j] +
                 temp);
        }

        for (size_t j = 0; j < num_groups; j++) {
            F(i, j + 2 * num_groups) = sim.get_result().get_value(t_idx)[sim.get_model().populations.get_flat_index(
                                           {(mio::AgeGroup)i, InfectionState::Susceptible})] *
                                       params.template get<TransmissionProbabilityOnContact<FP>>()[(mio::AgeGroup)i] *
                                       cont_freq_eff(i, j) * riskFromInfectedSymptomatic[(size_t)j] * divN[(size_t)j];
        }
    }

    //Initialize the matrix V
    for (Eigen::Index i = 0; i < (Eigen::Index)num_groups; i++) {
        FP criticalPerSevereAdjusted = smoother_cosine<FP>(
            icu_occupancy, 0.90 * params.template get<ICUCapacity<FP>>(), params.template get<ICUCapacity<FP>>(),
            params.template get<CriticalPerSevere<FP>>()[(mio::AgeGroup)i], 0);

        V(i, i)                           = 1 / params.template get<TimeExposed<FP>>()[(mio::AgeGroup)i];
        V(i + num_groups, i)              = -1 / params.template get<TimeExposed<FP>>()[(mio::AgeGroup)i];
        V(i + num_groups, i + num_groups) = 1 / params.template get<TimeInfectedNoSymptoms<FP>>()[(mio::AgeGroup)i];
        V(i + 2 * num_groups, i + num_groups) =
            -(1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[(mio::AgeGroup)i]) /
            params.template get<TimeInfectedNoSymptoms<FP>>()[(mio::AgeGroup)i];
        V(i + 2 * num_groups, i + 2 * num_groups) =
            (1 / params.template get<TimeInfectedSymptoms<FP>>()[(mio::AgeGroup)i]);
        V(i + 3 * num_groups, i + 2 * num_groups) =
            -params.template get<SeverePerInfectedSymptoms<FP>>()[(mio::AgeGroup)i] /
            params.template get<TimeInfectedSymptoms<FP>>()[(mio::AgeGroup)i];
        V(i + 3 * num_groups, i + 3 * num_groups) =
            1 / (params.template get<TimeInfectedSevere<FP>>()[(mio::AgeGroup)i]);
        V(i + 4 * num_groups, i + 3 * num_groups) =
            -criticalPerSevereAdjusted / (params.template get<TimeInfectedSevere<FP>>()[(mio::AgeGroup)i]);

        for (size_t j = 0; j < num_groups; j++) {
            V(i + 4 * num_groups, j + 4 * num_groups) = J(i, j);
        }
    }

    V = V.inverse();

    Eigen::MatrixX<ScalarType> NextGenMatrix(num_infected_compartments * num_groups, 5 * num_groups);
    NextGenMatrix.noalias() = F * V;

    // Compute the largest eigenvalue in absolute value
    Eigen::ComplexEigenSolver<Eigen::MatrixX<ScalarType>> ces;
    ces.compute(NextGenMatrix);
    FP rho = ces.eigenvalues().cwiseAbs().maxCoeff();

    return mio::success(rho);
}

/**
*@brief Computes the reproduction number for all time points of the Model output obtained by the Simulation.
*@param sim The Model Simulation.
*@tparam Base simulation type that uses a SECIR compartment model. see Simulation.
*@returns Eigen::Vector containing all reproduction numbers
*/

template <typename FP, class Base>
Eigen::VectorX<FP> get_reproduction_numbers(const Simulation<FP, Base>& sim)
{
    Eigen::VectorX<FP> temp(sim.get_result().get_num_time_points());
    for (int i = 0; i < sim.get_result().get_num_time_points(); i++) {
        temp[i] = get_reproduction_number<ScalarType>((size_t)i, sim).value();
    }
    return temp;
}

/**
*@brief @brief Computes the reproduction number at a given time point of the Simulation. If the particular time point is not part of the output, a linearly interpolated value is returned.
*@param t_value The time point at which the reproduction number should be computed.
*@param sim The Model Simulation.
*@tparam Base simulation type that uses a SECIR compartment model. see Simulation.
*@returns The computed reproduction number at the provided time point, potentially using linear interpolation.
*/
template <typename FP, class Base>
IOResult<FP> get_reproduction_number(FP t_value, const Simulation<FP, Base>& sim)
{
    if (t_value < sim.get_result().get_time(0) || t_value > sim.get_result().get_last_time()) {
        return mio::failure(mio::StatusCode::OutOfRange,
                            "Cannot interpolate reproduction number outside computed horizon of the TimeSeries");
    }

    if (t_value == sim.get_result().get_time(0)) {
        return mio::success(get_reproduction_number((size_t)0, sim).value());
    }

    const auto& times = sim.get_result().get_times();

    auto time_late = std::distance(times.begin(), std::lower_bound(times.begin(), times.end(), t_value));

    FP y1 = get_reproduction_number(static_cast<size_t>(time_late - 1), sim).value();
    FP y2 = get_reproduction_number(static_cast<size_t>(time_late), sim).value();

    auto result = linear_interpolation(t_value, sim.get_result().get_time(time_late - 1),
                                       sim.get_result().get_time(time_late), y1, y2);
    return mio::success(static_cast<FP>(result));
}

/**
 * Get mobility factors.
 * Used by mobility graph simulation.
 * Like infection risk, mobility of infected individuals is reduced if they are well isolated.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @return vector expression, same size as y, with mobility factors per compartment.
 * @tparam FP floating point type, e.g., double.
 * @tparam Base simulation type that uses a secir compartment model; see Simulation.
 */
template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
auto get_mobility_factors(const Simulation<FP, Base>& sim, FP /*t*/, const Eigen::Ref<const Eigen::VectorX<FP>>& y)
{
    auto& params = sim.get_model().parameters;
    //parameters as arrays
    auto&& p_asymp   = params.template get<RecoveredPerInfectedNoSymptoms<FP>>().array().template cast<FP>();
    auto&& p_inf     = params.template get<RiskOfInfectionFromSymptomatic<FP>>().array().template cast<FP>();
    auto&& p_inf_max = params.template get<MaxRiskOfInfectionFromSymptomatic<FP>>().array().template cast<FP>();
    //slice of InfectedNoSymptoms
    auto y_INS = slice(y, {Eigen::Index(InfectionState::InfectedNoSymptoms),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)});

    //compute isolation, same as infection risk from main model
    auto test_and_trace_required =
        ((1 - p_asymp) / params.template get<TimeInfectedNoSymptoms<FP>>().array().template cast<FP>() * y_INS.array())
            .sum();
    auto test_and_trace_capacity          = FP(params.template get<TestAndTraceCapacity<FP>>());
    auto test_and_trace_capacity_max_risk = FP(params.template get<TestAndTraceCapacityMaxRisk<FP>>());
    auto riskFromInfectedSymptomatic      = smoother_cosine<FP>(test_and_trace_required, test_and_trace_capacity,
                                                           test_and_trace_capacity * test_and_trace_capacity_max_risk,
                                                           p_inf.matrix(), p_inf_max.matrix());

    //set factor for infected
    auto factors = Eigen::VectorX<FP>::Ones(y.rows()).eval();
    slice(factors, {Eigen::Index(InfectionState::InfectedSymptoms), Eigen::Index(size_t(params.get_num_groups())),
                    Eigen::Index(InfectionState::Count)})
        .array() = riskFromInfectedSymptomatic;
    return factors;
}

template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
auto test_commuters(Simulation<FP, Base>& sim, Eigen::Ref<Eigen::VectorX<FP>> mobile_population, FP time)
{
    auto& model     = sim.get_model();
    FP nondetection = 1.0;
    if (time >= model.parameters.get_start_commuter_detection() &&
        time < model.parameters.get_end_commuter_detection()) {
        nondetection = model.parameters.get_commuter_nondetection();
    }
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); ++i) {
        auto INSi  = model.populations.get_flat_index({i, InfectionState::InfectedNoSymptoms});
        auto INSCi = model.populations.get_flat_index({i, InfectionState::InfectedNoSymptomsConfirmed});
        auto ISyi  = model.populations.get_flat_index({i, InfectionState::InfectedSymptoms});
        auto ISyCi = model.populations.get_flat_index({i, InfectionState::InfectedSymptomsConfirmed});

        //put detected commuters in their own compartment so they don't contribute to infections in their home node
        sim.get_result().get_last_value()[INSi] -= mobile_population[INSi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSCi] += mobile_population[INSi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISyi] -= mobile_population[ISyi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISyCi] += mobile_population[ISyi] * (1 - nondetection);

        //reduce the number of commuters
        mobile_population[ISyi] *= nondetection;
        mobile_population[INSi] *= nondetection;
    }
}

} // namespace omeng
} // namespace mio

#endif // ODEMENG_MODEL_H
