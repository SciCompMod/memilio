/* 
* Copyright (C) 2020-2024 MEmilio
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
#ifndef ODESECIR_MODEL_H
#define ODESECIR_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/epidemiology/populations.h"
#include "ode_secir/analyze_result.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"
#include "memilio/math/smoother.h"
#include "memilio/math/eigen_util.h"
#include "memilio/math/interpolation.h"

namespace mio
{
namespace osecir
{

// Create template specializations for the age resolved
// SECIHURD model
// clang-format off
using Flows = TypeList<Flow<InfectionState::Susceptible,                 InfectionState::Exposed>,
                       Flow<InfectionState::Exposed,                     InfectionState::InfectedNoSymptoms>,
                       Flow<InfectionState::InfectedNoSymptoms,          InfectionState::InfectedSymptoms>,
                       Flow<InfectionState::InfectedNoSymptoms,          InfectionState::Recovered>,
                       Flow<InfectionState::InfectedNoSymptomsConfirmed, InfectionState::InfectedSymptomsConfirmed>,
                       Flow<InfectionState::InfectedNoSymptomsConfirmed, InfectionState::Recovered>,
                       Flow<InfectionState::InfectedSymptoms,            InfectionState::InfectedSevere>,
                       Flow<InfectionState::InfectedSymptoms,            InfectionState::Recovered>,
                       Flow<InfectionState::InfectedSymptomsConfirmed,   InfectionState::InfectedSevere>,
                       Flow<InfectionState::InfectedSymptomsConfirmed,   InfectionState::Recovered>,
                       Flow<InfectionState::InfectedSevere,              InfectionState::InfectedCritical>,
                       Flow<InfectionState::InfectedSevere,              InfectionState::Recovered>,
                       Flow<InfectionState::InfectedSevere,              InfectionState::Dead>,
                       Flow<InfectionState::InfectedCritical,            InfectionState::Dead>,
                       Flow<InfectionState::InfectedCritical,            InfectionState::Recovered>>;
// clang-format on

template <typename FP = ScalarType>
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

    void get_flows(Eigen::Ref<const Vector<FP>> pop, Eigen::Ref<const Vector<FP>> y, FP t,
                   Eigen::Ref<Vector<FP>> flows) const override
    {
        auto const& params   = this->parameters;
        AgeGroup n_agegroups = params.get_num_groups();

        ContactMatrixGroup const& contact_matrix = params.template get<ContactPatterns<FP>>();

        auto icu_occupancy           = 0.0;
        auto test_and_trace_required = 0.0;
        for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
            test_and_trace_required += (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                                       params.template get<TimeInfectedNoSymptoms<FP>>()[i] *
                                       this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptoms});
            icu_occupancy += this->populations.get_from(pop, {i, InfectionState::InfectedCritical});
        }

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t Si    = this->populations.get_flat_index({i, InfectionState::Susceptible});
            size_t Ei    = this->populations.get_flat_index({i, InfectionState::Exposed});
            size_t INSi  = this->populations.get_flat_index({i, InfectionState::InfectedNoSymptoms});
            size_t INSCi = this->populations.get_flat_index({i, InfectionState::InfectedNoSymptomsConfirmed});
            size_t ISyi  = this->populations.get_flat_index({i, InfectionState::InfectedSymptoms});
            size_t ISyCi = this->populations.get_flat_index({i, InfectionState::InfectedSymptomsConfirmed});
            size_t ISevi = this->populations.get_flat_index({i, InfectionState::InfectedSevere});
            size_t ICri  = this->populations.get_flat_index({i, InfectionState::InfectedCritical});

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {
                size_t Sj    = this->populations.get_flat_index({j, InfectionState::Susceptible});
                size_t Ej    = this->populations.get_flat_index({j, InfectionState::Exposed});
                size_t INSj  = this->populations.get_flat_index({j, InfectionState::InfectedNoSymptoms});
                size_t ISyj  = this->populations.get_flat_index({j, InfectionState::InfectedSymptoms});
                size_t ISevj = this->populations.get_flat_index({j, InfectionState::InfectedSevere});
                size_t ICrj  = this->populations.get_flat_index({j, InfectionState::InfectedCritical});
                size_t Rj    = this->populations.get_flat_index({j, InfectionState::Recovered});

                //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
                auto riskFromInfectedSymptomatic =
                    smoother_cosine(test_and_trace_required, params.template get<TestAndTraceCapacity<FP>>(),
                                    params.template get<TestAndTraceCapacity<FP>>() *
                                        params.template get<TestAndTraceCapacityMaxRisk<FP>>(),
                                    params.template get<RiskOfInfectionFromSymptomatic<FP>>()[j],
                                    params.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[j]);

                // effective contact rate by contact rate between groups i and j and damping j
                double season_val =
                    (1 + params.template get<Seasonality<FP>>() *
                             sin(3.141592653589793 * ((params.template get<StartDay>() + t) / 182.5 + 0.5)));
                double cont_freq_eff =
                    season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                 static_cast<Eigen::Index>((size_t)j));
                double Nj =
                    pop[Sj] + pop[Ej] + pop[INSj] + pop[ISyj] + pop[ISevj] + pop[ICrj] + pop[Rj]; // without died people
                double divNj   = 1.0 / Nj; // precompute 1.0/Nj
                double dummy_S = y[Si] * cont_freq_eff * divNj *
                                 params.template get<TransmissionProbabilityOnContact<FP>>()[i] *
                                 (params.template get<RelativeTransmissionNoSymptoms<FP>>()[j] * pop[INSj] +
                                  riskFromInfectedSymptomatic * pop[ISyj]);

                // Susceptible -> Exposed
                flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>({i})] +=
                    dummy_S;
            }

            // ICU capacity shortage is close
            double criticalPerSevereAdjusted = smoother_cosine(
                icu_occupancy, 0.90 * params.template get<ICUCapacity<FP>>(), params.template get<ICUCapacity<FP>>(),
                params.template get<CriticalPerSevere<FP>>()[i], 0);

            double deathsPerSevereAdjusted =
                params.template get<CriticalPerSevere<FP>>()[i] - criticalPerSevereAdjusted;

            // Exposed -> InfectedNoSymptoms
            flows[this->template get_flat_flow_index<InfectionState::Exposed, InfectionState::InfectedNoSymptoms>(
                {i})] = (1 / params.template get<TimeExposed<FP>>()[i]) * y[Ei];

            // InfectedNoSymptoms -> InfectedSymptoms / Recovered
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptoms,
                                                     InfectionState::InfectedSymptoms>({i})] =
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) *
                (1 / params.template get<TimeInfectedNoSymptoms<FP>>()[i]) * y[INSi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptoms, InfectionState::Recovered>(
                {i})] = params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] *
                        (1 / params.template get<TimeInfectedNoSymptoms<FP>>()[i]) * y[INSi];

            // InfectedNoSymptomsConfirmed -> InfectedSymptomsConfirmed / Recovered
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsConfirmed,
                                                     InfectionState::InfectedSymptomsConfirmed>({i})] =
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) *
                (1 / params.template get<TimeInfectedNoSymptoms<FP>>()[i]) * y[INSCi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsConfirmed,
                                                     InfectionState::Recovered>({i})] =
                params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] *
                (1 / params.template get<TimeInfectedNoSymptoms<FP>>()[i]) * y[INSCi];

            // InfectedSymptoms -> InfectedSevere / Recovered
            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptoms, InfectionState::InfectedSevere>(
                {i})] = params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                        params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptoms, InfectionState::Recovered>(
                {i})] = (1 - params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                        params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyi];

            // InfectedSymptomsConfirmed -> InfectedSevere / Recovered
            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsConfirmed,
                                                     InfectionState::InfectedSevere>({i})] =
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyCi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsConfirmed,
                                                     InfectionState::Recovered>({i})] =
                (1 - params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyCi];

            // InfectedSevere -> InfectedCritical / Recovered / Dead
            flows[this->template get_flat_flow_index<InfectionState::InfectedSevere, InfectionState::InfectedCritical>(
                {i})] = criticalPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedSevere, InfectionState::Recovered>({i})] =
                (1 - params.template get<CriticalPerSevere<FP>>()[i]) /
                params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedSevere, InfectionState::Dead>({i})] =
                deathsPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevi];

            // InfectedCritical -> Dead / Recovered
            flows[this->template get_flat_flow_index<InfectionState::InfectedCritical, InfectionState::Dead>({i})] =
                params.template get<DeathsPerCritical<FP>>()[i] / params.template get<TimeInfectedCritical<FP>>()[i] *
                y[ICri];
            flows[this->template get_flat_flow_index<InfectionState::InfectedCritical, InfectionState::Recovered>(
                {i})] = (1 - params.template get<DeathsPerCritical<FP>>()[i]) /
                        params.template get<TimeInfectedCritical<FP>>()[i] * y[ICri];
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
template <typename FP = ScalarType, class BaseT = mio::Simulation<FP, Model<FP>>>
class Simulation;

/**
 * get percentage of infections per total population.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @tparam FP floating point type, e.g., double.
 * @tparam Base simulation type that uses a secir compartment model. see Simulation.
 */
template <typename FP = ScalarType, class Base = mio::Simulation<FP, Model<FP>>>
double get_infections_relative(const Simulation<FP, Base>& model, FP t, const Eigen::Ref<const Vector<FP>>& y);

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
    Simulation(mio::osecir::Model<FP> const& model, FP t0 = 0., FP dt = 0.1)
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
    Eigen::Ref<Vector<FP>> advance(FP tmax)
    {
        auto& t_end_dyn_npis   = this->get_model().parameters.get_end_dynamic_npis();
        auto& dyn_npis         = this->get_model().parameters.template get<DynamicNPIsInfectedSymptoms<FP>>();
        auto& contact_patterns = this->get_model().parameters.template get<ContactPatterns<FP>>();

        FP delay_npi_implementation; // delay which is needed to implement a NPI that is criterion-dependent
        auto t        = BaseT::get_result().get_last_time();
        const auto dt = dyn_npis.get_thresholds().size() > 0 ? dyn_npis.get_interval().get() : tmax;

        while (t < tmax) {
            auto dt_eff = std::min({dt, tmax - t, m_t_last_npi_check + dt - t});

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
                if (floating_point_greater_equal(t, m_t_last_npi_check + dt)) {
                    if (t < t_end_dyn_npis) {
                        auto inf_rel = get_infections_relative<FP>(*this, t, this->get_result().get_last_value()) *
                                       dyn_npis.get_base_value();
                        auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
                        if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
                            (exceeded_threshold->first > m_dynamic_npi.first ||
                             t > ScalarType(m_dynamic_npi.second))) { //old npi was weaker or is expired

                            auto t_start  = SimulationTime(t + delay_npi_implementation);
                            auto t_end    = t_start + SimulationTime(ScalarType(dyn_npis.get_duration()));
                            m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
                            implement_dynamic_npis(contact_patterns.get_cont_freq_mat(), exceeded_threshold->second,
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
    double m_t_last_npi_check;
    std::pair<double, SimulationTime> m_dynamic_npi = {-std::numeric_limits<double>::max(), mio::SimulationTime(0)};
};

/**
 * @brief Specialization of simulate for SECIR models using Simulation.
 *
 * @tparam FP floating point type, e.g., double.
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIR model to simulate.
 * @param[in] integrator optional integrator, uses rk45 if nullptr.
 * 
 * @return Returns the result of the simulation.
 */
template <typename FP = ScalarType>
inline auto simulate(FP t0, FP tmax, FP dt, const Model<FP>& model,
                     std::shared_ptr<IntegratorCore<FP>> integrator = nullptr)
{
    return mio::simulate<FP, Model<FP>, Simulation<>>(t0, tmax, dt, model, integrator);
}

/**
 * @brief Specialization of simulate for SECIR models using the FlowSimulation.
 * 
 * @tparam FP floating point type, e.g., double.
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIR model to simulate.
 * @param[in] integrator optional integrator, uses rk45 if nullptr.
 * 
 * @return Returns the result of the Flowsimulation.
 */
template <typename FP = ScalarType>
inline auto simulate_flows(FP t0, FP tmax, FP dt, const Model<FP>& model,
                           std::shared_ptr<IntegratorCore<FP>> integrator = nullptr)
{
    return mio::simulate_flows<FP, Model<FP>, Simulation<FP, mio::FlowSimulation<FP, Model<FP>>>>(t0, tmax, dt, model,
                                                                                                  integrator);
}

//see declaration above.
template <typename FP, class Base>
double get_infections_relative(const Simulation<FP, Base>& sim, FP /* t*/, const Eigen::Ref<const Vector<FP>>& y)
{
    double sum_inf = 0;
    for (auto i = AgeGroup(0); i < sim.get_model().parameters.get_num_groups(); ++i) {
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedSymptoms});
    }
    auto inf_rel = sum_inf / sim.get_model().populations.get_total();

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

    if (!(t_idx < static_cast<size_t>(sim.get_result().get_num_time_points()))) {
        return mio::failure(mio::StatusCode::OutOfRange, "t_idx is not a valid index for the TimeSeries");
    }

    auto const& params      = sim.get_model().parameters;
    const size_t num_groups = (size_t)sim.get_model().parameters.get_num_groups();
    //The infected compartments in the SECIR Model are: Exposed, Carrier, Infected, Hospitalized, ICU in respective agegroups
    const size_t num_infected_compartments   = 5;
    const size_t total_infected_compartments = num_infected_compartments * num_groups;
    const double pi                          = std::acos(-1);

    //F encodes new Infections and V encodes transition times in the next-generation matrix calculation of R_t
    Eigen::MatrixXd F(total_infected_compartments, total_infected_compartments);
    Eigen::MatrixXd V(total_infected_compartments, total_infected_compartments);
    F = Eigen::MatrixXd::Zero(total_infected_compartments,
                              total_infected_compartments); //Initialize matrices F and V with zeroes
    V = Eigen::MatrixXd::Zero(total_infected_compartments, total_infected_compartments);

    auto test_and_trace_required = 0.0;
    auto icu_occupancy           = 0.0;
    for (auto i = AgeGroup(0); i < (mio::AgeGroup)num_groups; ++i) {
        test_and_trace_required +=
            (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
            params.template get<TimeInfectedNoSymptoms<FP>>()[i] *
            sim.get_result().get_value(
                t_idx)[sim.get_model().populations.get_flat_index({i, InfectionState::InfectedNoSymptoms})];
        icu_occupancy += sim.get_result().get_value(
            t_idx)[sim.get_model().populations.get_flat_index({i, InfectionState::InfectedCritical})];
    }

    double season_val =
        (1 + params.template get<Seasonality<FP>>() *
                 sin(pi *
                     ((sim.get_model().parameters.template get<StartDay>() + sim.get_result().get_time(t_idx)) / 182.5 +
                      0.5)));
    ContactMatrixGroup const& contact_matrix = sim.get_model().parameters.template get<ContactPatterns<FP>>();

    Eigen::MatrixXd cont_freq_eff(num_groups, num_groups);
    Eigen::MatrixXd riskFromInfectedSymptomatic_derivatives(num_groups, num_groups);
    Eigen::VectorXd divN(num_groups);
    Eigen::VectorXd riskFromInfectedSymptomatic(num_groups);

    for (mio::AgeGroup k = 0; k < (mio::AgeGroup)num_groups; k++) {
        double temp = sim.get_result().get_value(
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
        if (temp == 0) {
            temp = 1;
        }
        divN[(size_t)k] = 1 / temp;

        riskFromInfectedSymptomatic[(size_t)k] = smoother_cosine(
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
                    std::sin(pi / (4 * params.template get<TestAndTraceCapacity<FP>>()) *
                             (test_and_trace_required - params.template get<TestAndTraceCapacity<FP>>()));
            }
        }

        for (Eigen::Index l = 0; l < (Eigen::Index)num_groups; l++) {
            cont_freq_eff(l, (size_t)k) =
                season_val * contact_matrix.get_matrix_at(static_cast<double>(t_idx))(
                                 static_cast<Eigen::Index>((size_t)l), static_cast<Eigen::Index>((size_t)k));
        }
    }

    //Check criterion if matrix V will be invertible by checking if subblock J is invertible
    Eigen::MatrixXd J(num_groups, num_groups);
    J = Eigen::MatrixXd::Zero(num_groups, num_groups);
    for (size_t i = 0; i < num_groups; i++) {
        J(i, i) = 1 / (params.template get<TimeInfectedCritical<FP>>()[(mio::AgeGroup)i]);

        if (!(icu_occupancy < 0.9 * params.template get<ICUCapacity<FP>>() ||
              icu_occupancy > (double)(params.template get<ICUCapacity<FP>>()))) {
            for (size_t j = 0; j < num_groups; j++) {
                J(i, j) -= sim.get_result().get_value(t_idx)[sim.get_model().populations.get_flat_index(
                               {(mio::AgeGroup)i, InfectionState::InfectedSevere})] /
                           params.template get<TimeInfectedSevere<FP>>()[(mio::AgeGroup)i] * 5 *
                           params.template get<CriticalPerSevere<FP>>()[(mio::AgeGroup)i] * pi /
                           (params.template get<ICUCapacity<FP>>()) *
                           std::sin(pi / (0.1 * params.template get<ICUCapacity<FP>>()) *
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

            double temp = 0;
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
        double criticalPerSevereAdjusted = smoother_cosine(
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

    //Compute F*V
    Eigen::MatrixXd NextGenMatrix(num_infected_compartments * num_groups, 5 * num_groups);
    NextGenMatrix = F * V;

    //Compute the largest eigenvalue in absolute value
    Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces;

    ces.compute(NextGenMatrix);
    const Eigen::VectorXcd eigen_vals = ces.eigenvalues();

    Eigen::VectorXd eigen_vals_abs;
    eigen_vals_abs.resize(eigen_vals.size());

    for (int i = 0; i < eigen_vals.size(); i++) {
        eigen_vals_abs[i] = std::abs(eigen_vals[i]);
    }
    return mio::success(eigen_vals_abs.maxCoeff());
}
/**
*@brief Computes the reproduction number for all time points of the Model output obtained by the Simulation.
*@param sim The Model Simulation.
*@tparam Base simulation type that uses a SECIR compartment model. see Simulation.
*@returns Eigen::Vector containing all reproduction numbers
*/

template <typename FP, class Base>
Vector<FP> get_reproduction_numbers(const Simulation<FP, Base>& sim)
{
    Vector<FP> temp(sim.get_result().get_num_time_points());
    for (int i = 0; i < sim.get_result().get_num_time_points(); i++) {
        temp[i] = get_reproduction_number((size_t)i, sim).value();
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
template <class Base>
IOResult<ScalarType> get_reproduction_number(ScalarType t_value, const Simulation<Base>& sim)
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

    ScalarType y1 = get_reproduction_number(static_cast<size_t>(time_late - 1), sim).value();
    ScalarType y2 = get_reproduction_number(static_cast<size_t>(time_late), sim).value();

    auto result = linear_interpolation(t_value, sim.get_result().get_time(time_late - 1),
                                       sim.get_result().get_time(time_late), y1, y2);
    return mio::success(static_cast<ScalarType>(result));
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
template <typename FP = ScalarType, class Base = mio::Simulation<Model<FP>, FP>>
auto get_mobility_factors(const Simulation<Base>& sim, FP /*t*/, const Eigen::Ref<const Vector<FP>>& y)
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
    auto test_and_trace_capacity          = double(params.template get<TestAndTraceCapacity<FP>>());
    auto test_and_trace_capacity_max_risk = double(params.template get<TestAndTraceCapacityMaxRisk<FP>>());
    auto riskFromInfectedSymptomatic =
        smoother_cosine(test_and_trace_required, test_and_trace_capacity,
                        test_and_trace_capacity * test_and_trace_capacity_max_risk, p_inf.matrix(), p_inf_max.matrix());

    //set factor for infected
    auto factors = Eigen::VectorXd::Ones(y.rows()).eval();
    slice(factors, {Eigen::Index(InfectionState::InfectedSymptoms), Eigen::Index(size_t(params.get_num_groups())),
                    Eigen::Index(InfectionState::Count)})
        .array() = riskFromInfectedSymptomatic;
    return factors;
}

template <typename FP = ScalarType, class Base = mio::Simulation<Model<FP>, FP>>
auto test_commuters(Simulation<FP, Base>& sim, Eigen::Ref<Vector<FP>> mobile_population, FP time)
{
    auto& model       = sim.get_model();
    auto nondetection = 1.0;
    if (time >= model.parameters.get_start_commuter_detection() &&
        time < model.parameters.get_end_commuter_detection()) {
        nondetection = (double)model.parameters.get_commuter_nondetection();
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

} // namespace osecir
} // namespace mio

#endif // ODESECIR_MODEL_H
