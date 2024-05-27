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
#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/populations.h"
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

class Model : public FlowModel<InfectionState, Populations<AgeGroup, InfectionState>, Parameters, Flows>
{
    using Base = FlowModel<InfectionState, mio::Populations<AgeGroup, InfectionState>, Parameters, Flows>;

public:
    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    Model(int num_agegroups)
        : Model(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                   Eigen::Ref<Eigen::VectorXd> flows) const override
    {
        auto const& params   = this->parameters;
        AgeGroup n_agegroups = params.get_num_groups();

        ContactMatrixGroup const& contact_matrix = params.get<ContactPatterns>();

        auto icu_occupancy           = 0.0;
        auto test_and_trace_required = 0.0;
        for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
            test_and_trace_required += (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) /
                                       params.get<TimeInfectedNoSymptoms>()[i] *
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
                auto riskFromInfectedSymptomatic = smoother_cosine(
                    test_and_trace_required, params.get<TestAndTraceCapacity>(), params.get<TestAndTraceCapacity>() * 5,
                    params.get<RiskOfInfectionFromSymptomatic>()[j],
                    params.get<MaxRiskOfInfectionFromSymptomatic>()[j]);

                // effective contact rate by contact rate between groups i and j and damping j
                double season_val =
                    (1 + params.get<Seasonality>() *
                             sin(3.141592653589793 * (std::fmod((params.get<StartDay>() + t), 365.0) / 182.5 + 0.5)));
                double cont_freq_eff =
                    season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                 static_cast<Eigen::Index>((size_t)j));
                double Nj =
                    pop[Sj] + pop[Ej] + pop[INSj] + pop[ISyj] + pop[ISevj] + pop[ICrj] + pop[Rj]; // without died people
                double divNj   = 1.0 / Nj; // precompute 1.0/Nj
                double dummy_S = y[Si] * cont_freq_eff * divNj * params.get<TransmissionProbabilityOnContact>()[i] *
                                 (params.get<RelativeTransmissionNoSymptoms>()[j] * pop[INSj] +
                                  riskFromInfectedSymptomatic * pop[ISyj]);

                // Susceptible -> Exposed
                flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>({i})] += dummy_S;
            }

            // ICU capacity shortage is close
            double criticalPerSevereAdjusted =
                smoother_cosine(icu_occupancy, 0.90 * params.get<ICUCapacity>(), params.get<ICUCapacity>(),
                                params.get<CriticalPerSevere>()[i], 0);

            double deathsPerSevereAdjusted = params.get<CriticalPerSevere>()[i] - criticalPerSevereAdjusted;

            // Exposed -> InfectedNoSymptoms
            flows[get_flat_flow_index<InfectionState::Exposed, InfectionState::InfectedNoSymptoms>({i})] =
                (1 / params.get<TimeExposed>()[i]) * y[Ei];

            // InfectedNoSymptoms -> InfectedSymptoms / Recovered
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptoms, InfectionState::InfectedSymptoms>({i})] =
                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * (1 / params.get<TimeInfectedNoSymptoms>()[i]) *
                y[INSi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptoms, InfectionState::Recovered>({i})] =
                params.get<RecoveredPerInfectedNoSymptoms>()[i] * (1 / params.get<TimeInfectedNoSymptoms>()[i]) *
                y[INSi];

            // InfectedNoSymptomsConfirmed -> InfectedSymptomsConfirmed / Recovered
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsConfirmed,
                                      InfectionState::InfectedSymptomsConfirmed>({i})] =
                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * (1 / params.get<TimeInfectedNoSymptoms>()[i]) *
                y[INSCi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsConfirmed, InfectionState::Recovered>({i})] =
                params.get<RecoveredPerInfectedNoSymptoms>()[i] * (1 / params.get<TimeInfectedNoSymptoms>()[i]) *
                y[INSCi];

            // InfectedSymptoms -> InfectedSevere / Recovered
            flows[get_flat_flow_index<InfectionState::InfectedSymptoms, InfectionState::InfectedSevere>({i})] =
                params.get<SeverePerInfectedSymptoms>()[i] / params.get<TimeInfectedSymptoms>()[i] * y[ISyi];
            flows[get_flat_flow_index<InfectionState::InfectedSymptoms, InfectionState::Recovered>({i})] =
                (1 - params.get<SeverePerInfectedSymptoms>()[i]) / params.get<TimeInfectedSymptoms>()[i] * y[ISyi];

            // InfectedSymptomsConfirmed -> InfectedSevere / Recovered
            flows[get_flat_flow_index<InfectionState::InfectedSymptomsConfirmed, InfectionState::InfectedSevere>({i})] =
                params.get<SeverePerInfectedSymptoms>()[i] / params.get<TimeInfectedSymptoms>()[i] * y[ISyCi];
            flows[get_flat_flow_index<InfectionState::InfectedSymptomsConfirmed, InfectionState::Recovered>({i})] =
                (1 - params.get<SeverePerInfectedSymptoms>()[i]) / params.get<TimeInfectedSymptoms>()[i] * y[ISyCi];

            // InfectedSevere -> InfectedCritical / Recovered / Dead
            flows[get_flat_flow_index<InfectionState::InfectedSevere, InfectionState::InfectedCritical>({i})] =
                criticalPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevi];
            flows[get_flat_flow_index<InfectionState::InfectedSevere, InfectionState::Recovered>({i})] =
                (1 - params.get<CriticalPerSevere>()[i]) / params.get<TimeInfectedSevere>()[i] * y[ISevi];
            flows[get_flat_flow_index<InfectionState::InfectedSevere, InfectionState::Dead>({i})] =
                deathsPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevi];

            // InfectedCritical -> Dead / Recovered
            flows[get_flat_flow_index<InfectionState::InfectedCritical, InfectionState::Dead>({i})] =
                params.get<DeathsPerCritical>()[i] / params.get<TimeInfectedCritical>()[i] * y[ICri];
            flows[get_flat_flow_index<InfectionState::InfectedCritical, InfectionState::Recovered>({i})] =
                (1 - params.get<DeathsPerCritical>()[i]) / params.get<TimeInfectedCritical>()[i] * y[ICri];
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
        obj.add_element("Parameters", parameters);
        obj.add_element("Populations", populations);
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
template <class BaseT = mio::Simulation<Model>>
class Simulation;

template <class BaseT = mio::Simulation<Model>>
class FeedbackSimulation;

/**
 * get percentage of infections per total population.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @tparam Base simulation type that uses a secir compartment model. see Simulation.
 */
template <class Base = mio::Simulation<Model>>
double get_infections_relative(const Simulation<Base>& model, double t, const Eigen::Ref<const Eigen::VectorXd>& y);

/**
 * specialization of compartment model simulation for secir models.
 * @tparam BaseT simulation type that uses a secir compartment model. default mio::Simulation. For testing purposes only!
 */
template <class BaseT>
class Simulation : public BaseT
{
public:
    /**
     * construct a simulation.
     * @param model the model to simulate.
     * @param t0 start time
     * @param dt time steps
     */
    Simulation(Model const& model, double t0 = 0., double dt = 0.1)
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
    Eigen::Ref<Eigen::VectorXd> advance(double tmax)
    {
        auto& dyn_npis         = this->get_model().parameters.template get<DynamicNPIsInfectedSymptoms>();
        auto& contact_patterns = this->get_model().parameters.template get<ContactPatterns>();
        if (dyn_npis.get_thresholds().size() > 0) {
            auto t        = BaseT::get_result().get_last_time();
            const auto dt = dyn_npis.get_interval().get();

            while (t < tmax) {
                auto dt_eff = std::min({dt, tmax - t, m_t_last_npi_check + dt - t});

                BaseT::advance(t + dt_eff);
                t = t + dt_eff;

                if (floating_point_greater_equal(t, m_t_last_npi_check + dt)) {
                    auto inf_rel = get_infections_relative(*this, t, this->get_result().get_last_value()) *
                                   dyn_npis.get_base_value();
                    auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
                    if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
                        (exceeded_threshold->first > m_dynamic_npi.first ||
                         t > double(m_dynamic_npi.second))) { //old npi was weaker or is expired
                        auto t_end    = mio::SimulationTime(t + double(dyn_npis.get_duration()));
                        m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
                        mio::implement_dynamic_npis(contact_patterns.get_cont_freq_mat(), exceeded_threshold->second,
                                                    SimulationTime(t), t_end, [](auto& g) {
                                                        return mio::make_contact_damping_matrix(g);
                                                    });
                    }

                    m_t_last_npi_check = t;
                }
            }

            return this->get_result().get_last_value();
        }
        else {
            return BaseT::advance(tmax);
        }
    }

private:
    double m_t_last_npi_check;
    std::pair<double, SimulationTime> m_dynamic_npi = {-std::numeric_limits<double>::max(), mio::SimulationTime(0)};
};

template <class BaseT>
class FeedbackSimulation : public BaseT
{
public:
    /**
     * construct a simulation.
     * @param model the model to simulate.
     * @param t0 start time
     * @param dt time steps
     */
    FeedbackSimulation(Model const& model, double t0 = 0., double dt = 0.1)
        : BaseT(model, t0, dt)
        , m_t_last_npi_check(t0)
        , m_perceived_risk((size_t)model.parameters.get_num_groups())
    {
        m_model_ptr = std::make_shared<Model>(model);
        m_models.push_back(m_model_ptr);
    }

    /**
     * @brief Destructor for FeedbackSimulation.
     */
    ~FeedbackSimulation()
    {
        auto it = std::find(m_models.begin(), m_models.end(), m_model_ptr);
        if (it != m_models.end()) {
            m_models.erase(it);
        }
    }

    FeedbackSimulation(FeedbackSimulation&& other) noexcept
        : BaseT(std::move(other))
        , m_t_last_npi_check(other.m_t_last_npi_check)
        , m_dynamic_npi(std::move(other.m_dynamic_npi))
        , m_model_ptr(std::move(other.m_model_ptr))
        , m_perceived_risk(std::move(other.m_perceived_risk))
    {
    }

    /**
     * @brief Adjusts contact patterns based on ICU occupancy feedback.
     * 
     * @param t Current simulation time.
     * @param icu_regional Regional ICU occupancy data.
     * @param icu_national National ICU occupancy data.
     */
    void feedback_contacts(double t, Eigen::Ref<Eigen::MatrixXd> icu_regional, Eigen::Ref<Eigen::MatrixXd> icu_national)
    {
        auto& params                            = m_model_ptr->parameters;
        const size_t locations                  = params.template get<ContactReductionMax>().size();
        constexpr ScalarType threshold_softplus = 0.1;
        double perceived_risk_contacts =
            calc_risk_perceived(params.template get<alphaGammaContacts>(), params.template get<betaGammaContacts>(),
                                icu_regional, icu_national);

        auto group_weights_all = Eigen::VectorXd::Constant(size_t(params.get_num_groups()), 1.0);

        auto contact_reduction = [=](auto val) {
            auto v = mio::UncertainValue(val);
            return mio::DampingSampling(v, mio::DampingLevel(int(0)), mio::DampingType(int(0)), mio::SimulationTime(t),
                                        {size_t(0)}, group_weights_all);
        };

        for (size_t loc = 0; loc < locations; ++loc) {
            ScalarType reduc_fac_location = params.template get<ContactReductionMin>()[loc];
            if (perceived_risk_contacts > threshold_softplus) {
                reduc_fac_location =
                    (params.template get<ContactReductionMax>()[loc] -
                     params.template get<ContactReductionMin>()[loc]) *
                        params.template get<EpsilonContacts>() *
                        std::log(std::exp(perceived_risk_contacts / params.template get<EpsilonContacts>()) + 1.0) +
                    params.template get<ContactReductionMin>()[loc];
            }

            reduc_fac_location = std::min(reduc_fac_location, params.template get<ContactReductionMax>()[loc]);

            params.template get<ContactPatterns>().get_dampings().push_back(contact_reduction(reduc_fac_location));
            this->get_model().parameters.template get<ContactPatterns>().get_dampings().push_back(
                contact_reduction(reduc_fac_location));

            this->get_model().parameters.template get<ContactPatterns>().make_matrix();
            params.template get<ContactPatterns>().make_matrix();
        }
    }

    /**
     * @brief Calculates the perceived risk based on ICU occupancy data.
     * 
     * @param a Parameter alpha for the gamma distribution.
     * @param b Parameter beta for the gamma distribution.
     * @param icu_regional Regional ICU occupancy data.
     * @param icu_national National ICU occupancy data.
     * @return The perceived risk.
     */
    virtual ScalarType calc_risk_perceived(const ScalarType a, const ScalarType b,
                                           Eigen::Ref<Eigen::MatrixXd> icu_regional,
                                           Eigen::Ref<Eigen::MatrixXd> icu_national)
    {
        const auto& icu_occupancy    = m_model_ptr->parameters.template get<ICUOccupancyLocal>();
        const size_t num_time_points = icu_occupancy.get_num_time_points();
        const size_t n =
            std::min(static_cast<size_t>(num_time_points), m_model_ptr->parameters.template get<CutOffGamma>());
        ScalarType perceived_risk = 0.0;

        for (size_t i = num_time_points - n; i < num_time_points; ++i) {
            const size_t day       = i - (num_time_points - n);
            const ScalarType gamma = std::pow(b, a) * std::pow(day, a - 1) * std::exp(-b * day) / std::tgamma(a);
            auto mixed_regional_national_icu_occupancy =
                m_model_ptr->parameters.template get<BlendingFactorRegional>() * icu_regional.row(i).sum() +
                (1.0 - m_model_ptr->parameters.template get<BlendingFactorRegional>()) * icu_national.row(i).sum();

            auto icu_occupancy_adjusted_rel =
                (m_model_ptr->parameters.template get<BlendingFactorLocal>() * icu_occupancy.get_value(i).sum() +
                 (1.0 - m_model_ptr->parameters.template get<BlendingFactorRegional>()) *
                     mixed_regional_national_icu_occupancy) /
                m_model_ptr->parameters.template get<ICUCapacity>();

            // TODO: Should we allow values > 1.0? Maximal Reduction later is limited, but it may outweight other factors
            // icu_occupancy_adjusted_rel = std::min(icu_occupancy_adjusted_rel, 1.0); Now limited later
            perceived_risk += icu_occupancy_adjusted_rel * gamma;
        }
        perceived_risk  = std::min(perceived_risk, 1.0);
        auto num_groups = (size_t)m_model_ptr->parameters.get_num_groups();
        mio::unused(num_groups);
        m_perceived_risk.add_time_point(
            m_model_ptr->parameters.template get<ICUOccupancyLocal>().get_last_time(),
            Eigen::VectorXd::Constant((size_t)m_model_ptr->parameters.get_num_groups(), perceived_risk));
        return perceived_risk;
    }

    /**
     * @brief Adds ICU occupancy data at a given time point. The icu occupancy is calculated per 100,000 inhabitants.
     * 
     * @param t The current time point.
     */
    void add_icu_occupancy(double t)
    {
        auto& params           = m_model_ptr->parameters;
        size_t num_groups      = (size_t)params.get_num_groups();
        const auto& last_value = this->get_result().get_last_value();

        Eigen::VectorXd icu_occupancy(num_groups);
        for (size_t age = 0; age < num_groups; ++age) {
            auto indx_icu_naive =
                m_model_ptr->populations.get_flat_index({(AgeGroup)age, InfectionState::InfectedCritical});
            icu_occupancy[age] = last_value[indx_icu_naive];
        }

        icu_occupancy = icu_occupancy / m_model_ptr->populations.get_total() * 100000;
        params.template get<ICUOccupancyLocal>().add_time_point(t, icu_occupancy);
    }

    /**
     * @brief Calculates the global ICU occupancy from all models.
     * 
     * @return The global ICU occupancy data.
     */
    Eigen::MatrixXd calculate_global_icu_occupancy()
    {
        const size_t num_groups      = static_cast<size_t>(m_model_ptr->parameters.get_num_groups());
        const size_t num_time_points = m_model_ptr->parameters.template get<ICUOccupancyLocal>().get_num_time_points();
        Eigen::MatrixXd global_icu_occupancy = Eigen::MatrixXd::Zero(num_time_points, num_groups);
        // We need the total population to calculate the global ICU occupancy as weights for the local ICU occupancy
        ScalarType total_population = 0.0;

        for (const auto& model : m_models) {
            const auto& icu_occupancy_local = model->parameters.template get<ICUOccupancyLocal>();
            const auto model_population     = model->populations.get_total();
            total_population += model_population;
            for (size_t t = 0; t < num_time_points; ++t) {
                for (size_t age = 0; age < num_groups; ++age) {
                    global_icu_occupancy(t, age) += model_population * icu_occupancy_local.get_value(t)(age);
                }
            }
        }

        // Normalize the global ICU occupancy by the total population
        global_icu_occupancy /= total_population;

        return global_icu_occupancy;
    }

    /**
     * @brief Calculates the regional ICU occupancy for the same state ID as the current model.
     * 
     * @return The regional ICU occupancy data.
     */
    Eigen::MatrixXd calculate_regional_icu_occupancy()
    {
        const auto state_id          = m_model_ptr->parameters.template get<StateID>();
        const size_t num_groups      = static_cast<size_t>(m_model_ptr->parameters.get_num_groups());
        const size_t num_time_points = m_model_ptr->parameters.template get<ICUOccupancyLocal>().get_num_time_points();
        Eigen::MatrixXd regional_icu_occupancy = Eigen::MatrixXd::Zero(num_time_points, num_groups);
        // We need the total population per region to calculate the regional ICU occupancy as weights for the local ICU occupancy
        ScalarType regional_population = 0.0;

        for (const auto& model : m_models) {
            if (model->parameters.get<StateID>() == state_id) {
                const auto& icu_occupancy_local = model->parameters.template get<ICUOccupancyLocal>();
                const auto total_population     = model->populations.get_total();
                regional_population += total_population;
                for (size_t t = 0; t < num_time_points; ++t) {
                    for (size_t age = 0; age < num_groups; ++age) {
                        regional_icu_occupancy(t, age) += total_population * icu_occupancy_local.get_value(t)(age);
                    }
                }
            }
        }

        // Normalize the regional ICU occupancy by the total regional population
        regional_icu_occupancy /= regional_population;

        return regional_icu_occupancy;
    }

    /**
     * @brief Get the perceived risk time series.
     * 
     * @return The perceived risk time series.
     */
    auto& get_perceived_risk() const
    {
        return m_perceived_risk;
    }

    /**
     * @brief Advances the simulation to a specified time.
     * 
     * @param tmax The maximum time to advance to.
     * @return The last value of the simulation result.
     */
    Eigen::Ref<Eigen::VectorXd> advance(double tmax)
    {
        auto& dyn_npis         = m_model_ptr->parameters.template get<DynamicNPIsInfectedSymptoms>();
        auto& contact_patterns = m_model_ptr->parameters.template get<ContactPatterns>();
        auto t                 = BaseT::get_result().get_last_time();
        const auto dt          = dyn_npis.get_interval().get();

        while (t < tmax) {
            auto dt_eff = std::min({dt, tmax - t, m_t_last_npi_check + dt - t});
            if (dt_eff >= 1.0) {
                dt_eff = 1.0;
            }

            if (t + 0.5 + dt_eff - std::floor(t + 0.5) >= 1) {
                auto icu_regional = calculate_regional_icu_occupancy();
                auto icu_global   = calculate_global_icu_occupancy();
                this->feedback_contacts(t, icu_regional, icu_global);
            }
            BaseT::advance(t + dt_eff);
            if (t + 0.5 + dt_eff - std::floor(t + 0.5) >= 1) {
                this->add_icu_occupancy(t);
            }
            t = t + dt_eff;
            if (dyn_npis.get_thresholds().size() > 0) {
                if (floating_point_greater_equal(t, m_t_last_npi_check + dt)) {
                    auto inf_rel = get_infections_relative(*this, t, this->get_result().get_last_value()) *
                                   dyn_npis.get_base_value();
                    auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
                    if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
                        (exceeded_threshold->first > m_dynamic_npi.first ||
                         t > double(m_dynamic_npi.second))) { // old npi was weaker or is expired
                        auto t_end    = mio::SimulationTime(t + double(dyn_npis.get_duration()));
                        m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
                        mio::implement_dynamic_npis(contact_patterns.get_cont_freq_mat(), exceeded_threshold->second,
                                                    SimulationTime(t), t_end, [](auto& g) {
                                                        return mio::make_contact_damping_matrix(g);
                                                    });
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

    /**
     * @brief Get the model pointer.
     * 
     * @return The model pointer.
     */
    std::shared_ptr<Model> get_model_ptr() const
    {
        return m_model_ptr;
    }

private:
    double m_t_last_npi_check; /// Time of the last NPI check.
    std::pair<double, SimulationTime> m_dynamic_npi = {-std::numeric_limits<double>::max(),
                                                       mio::SimulationTime(0)}; /// Dynamic NPI data.
    std::shared_ptr<Model> m_model_ptr; /// Pointer to the model.
    inline static std::vector<std::shared_ptr<Model>> m_models; /// List of model pointers.
    mio::TimeSeries<ScalarType> m_perceived_risk; // = mio::TimeSeries<ScalarType>(6); /// Perceived risk time series.
};

/**
 * @brief Specialization of simulate for SECIR models using Simulation.
 * 
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIR model to simulate.
 * @param[in] integrator optional integrator, uses rk45 if nullptr.
 * 
 * @return Returns the result of the simulation.
 */
inline auto simulate(double t0, double tmax, double dt, const Model& model,
                     std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    return mio::simulate<Model, Simulation<>>(t0, tmax, dt, model, integrator);
}

/**
 * @brief Specialization of simulate for SECIR models using the FlowSimulation.
 * 
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIR model to simulate.
 * @param[in] integrator optional integrator, uses rk45 if nullptr.
 * 
 * @return Returns the result of the Flowsimulation.
 */
inline auto simulate_flows(double t0, double tmax, double dt, const Model& model,
                           std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    return mio::simulate_flows<Model, Simulation<mio::FlowSimulation<Model>>>(t0, tmax, dt, model, integrator);
}

inline auto simulate_feedback(double t0, double tmax, double dt, const Model& model,
                              std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    return mio::simulate<Model, FeedbackSimulation<>>(t0, tmax, dt, model, integrator);
}

//see declaration above.
template <class Base>
double get_infections_relative(const Simulation<Base>& sim, double /*t*/, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    double sum_inf = 0;
    for (auto i = AgeGroup(0); i < sim.get_model().parameters.get_num_groups(); ++i) {
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedSymptoms});
    }
    auto inf_rel = sum_inf / sim.get_model().populations.get_total();

    return inf_rel;
}

template <class Base>
double get_infections_relative(const FeedbackSimulation<Base>& sim, double /*t*/,
                               const Eigen::Ref<const Eigen::VectorXd>& y)
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

template <class Base>
IOResult<ScalarType> get_reproduction_number(size_t t_idx, const Simulation<Base>& sim)
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
            (1 - params.template get<RecoveredPerInfectedNoSymptoms>()[i]) /
            params.template get<TimeInfectedNoSymptoms>()[i] *
            sim.get_result().get_value(
                t_idx)[sim.get_model().populations.get_flat_index({i, InfectionState::InfectedNoSymptoms})];
        icu_occupancy += sim.get_result().get_value(
            t_idx)[sim.get_model().populations.get_flat_index({i, InfectionState::InfectedCritical})];
    }

    double season_val                        = (1 + params.template get<Seasonality>() *
                                 sin(pi * (std::fmod((sim.get_model().parameters.template get<StartDay>() +
                                                      sim.get_result().get_time(t_idx)),
                                                     365.0) /
                                               182.5 +
                                           0.5)));
    ContactMatrixGroup const& contact_matrix = sim.get_model().parameters.template get<ContactPatterns>();

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
            test_and_trace_required, params.template get<TestAndTraceCapacity>(),
            (params.template get<TestAndTraceCapacity>()) * 5, params.template get<RiskOfInfectionFromSymptomatic>()[k],
            params.template get<MaxRiskOfInfectionFromSymptomatic>()[k]);

        for (mio::AgeGroup l = 0; l < (mio::AgeGroup)num_groups; l++) {
            if (test_and_trace_required < params.template get<TestAndTraceCapacity>() ||
                test_and_trace_required > 5 * params.template get<TestAndTraceCapacity>()) {
                riskFromInfectedSymptomatic_derivatives((size_t)k, (size_t)l) = 0;
            }
            else {
                riskFromInfectedSymptomatic_derivatives((size_t)k, (size_t)l) =
                    -0.5 * pi *
                    (params.template get<MaxRiskOfInfectionFromSymptomatic>()[k] -
                     params.template get<RiskOfInfectionFromSymptomatic>()[k]) /
                    (4 * params.template get<TestAndTraceCapacity>()) *
                    (1 - params.template get<RecoveredPerInfectedNoSymptoms>()[l]) /
                    params.template get<TimeInfectedNoSymptoms>()[l] *
                    std::sin(pi / (4 * params.template get<TestAndTraceCapacity>()) *
                             (test_and_trace_required - params.template get<TestAndTraceCapacity>()));
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
        J(i, i) = 1 / (params.template get<TimeInfectedCritical>()[(mio::AgeGroup)i]);

        if (!(icu_occupancy < 0.9 * params.template get<ICUCapacity>() ||
              icu_occupancy > (double)(params.template get<ICUCapacity>()))) {
            for (size_t j = 0; j < num_groups; j++) {
                J(i, j) -= sim.get_result().get_value(t_idx)[sim.get_model().populations.get_flat_index(
                               {(mio::AgeGroup)i, InfectionState::InfectedSevere})] /
                           params.template get<TimeInfectedSevere>()[(mio::AgeGroup)i] * 5 *
                           params.template get<CriticalPerSevere>()[(mio::AgeGroup)i] * pi /
                           (params.template get<ICUCapacity>()) *
                           std::sin(pi / (0.1 * params.template get<ICUCapacity>()) *
                                    (icu_occupancy - 0.9 * params.template get<ICUCapacity>()));
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
                params.template get<TransmissionProbabilityOnContact>()[(mio::AgeGroup)i] *
                (cont_freq_eff(i, j) * params.template get<RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)j] *
                     divN[(size_t)j] +
                 temp);
        }

        for (size_t j = 0; j < num_groups; j++) {
            F(i, j + 2 * num_groups) = sim.get_result().get_value(t_idx)[sim.get_model().populations.get_flat_index(
                                           {(mio::AgeGroup)i, InfectionState::Susceptible})] *
                                       params.template get<TransmissionProbabilityOnContact>()[(mio::AgeGroup)i] *
                                       cont_freq_eff(i, j) * riskFromInfectedSymptomatic[(size_t)j] * divN[(size_t)j];
        }
    }

    //Initialize the matrix V
    for (Eigen::Index i = 0; i < (Eigen::Index)num_groups; i++) {

        double criticalPerSevereAdjusted = smoother_cosine(
            icu_occupancy, 0.90 * params.template get<ICUCapacity>(), params.template get<ICUCapacity>(),
            params.template get<CriticalPerSevere>()[(mio::AgeGroup)i], 0);

        V(i, i)                           = 1 / params.template get<TimeExposed>()[(mio::AgeGroup)i];
        V(i + num_groups, i)              = -1 / params.template get<TimeExposed>()[(mio::AgeGroup)i];
        V(i + num_groups, i + num_groups) = 1 / params.template get<TimeInfectedNoSymptoms>()[(mio::AgeGroup)i];
        V(i + 2 * num_groups, i + num_groups) =
            -(1 - params.template get<RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)i]) /
            params.template get<TimeInfectedNoSymptoms>()[(mio::AgeGroup)i];
        V(i + 2 * num_groups, i + 2 * num_groups) = (1 / params.template get<TimeInfectedSymptoms>()[(mio::AgeGroup)i]);
        V(i + 3 * num_groups, i + 2 * num_groups) =
            -params.template get<SeverePerInfectedSymptoms>()[(mio::AgeGroup)i] /
            params.template get<TimeInfectedSymptoms>()[(mio::AgeGroup)i];
        V(i + 3 * num_groups, i + 3 * num_groups) = 1 / (params.template get<TimeInfectedSevere>()[(mio::AgeGroup)i]);
        V(i + 4 * num_groups, i + 3 * num_groups) =
            -criticalPerSevereAdjusted / (params.template get<TimeInfectedSevere>()[(mio::AgeGroup)i]);

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

template <class Base>
Eigen::VectorXd get_reproduction_numbers(const Simulation<Base>& sim)
{
    Eigen::VectorXd temp(sim.get_result().get_num_time_points());
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
 * Get migration factors.
 * Used by migration graph simulation.
 * Like infection risk, migration of infected individuals is reduced if they are well isolated.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @return vector expression, same size as y, with migration factors per compartment.
 * @tparam Base simulation type that uses a secir compartment model; see Simulation.
 */
template <class Base = mio::Simulation<Model>>
auto get_migration_factors(const Simulation<Base>& sim, double /*t*/, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    auto& params = sim.get_model().parameters;
    //parameters as arrays
    auto&& p_asymp   = params.template get<RecoveredPerInfectedNoSymptoms>().array().template cast<double>();
    auto&& p_inf     = params.template get<RiskOfInfectionFromSymptomatic>().array().template cast<double>();
    auto&& p_inf_max = params.template get<MaxRiskOfInfectionFromSymptomatic>().array().template cast<double>();
    //slice of InfectedNoSymptoms
    auto y_INS = slice(y, {Eigen::Index(InfectionState::InfectedNoSymptoms),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)});

    //compute isolation, same as infection risk from main model
    auto test_and_trace_required =
        ((1 - p_asymp) / params.template get<TimeInfectedNoSymptoms>().array().template cast<double>() * y_INS.array())
            .sum();
    auto test_and_trace_capacity     = double(params.template get<TestAndTraceCapacity>());
    auto riskFromInfectedSymptomatic = smoother_cosine(test_and_trace_required, test_and_trace_capacity,
                                                       test_and_trace_capacity * 5, p_inf.matrix(), p_inf_max.matrix());

    //set factor for infected
    auto factors = Eigen::VectorXd::Ones(y.rows()).eval();
    slice(factors, {Eigen::Index(InfectionState::InfectedSymptoms), Eigen::Index(size_t(params.get_num_groups())),
                    Eigen::Index(InfectionState::Count)})
        .array() = riskFromInfectedSymptomatic;
    return factors;
}

template <class Base = mio::Simulation<Model>>
auto test_commuters(Simulation<Base>& sim, Eigen::Ref<Eigen::VectorXd> migrated, double time)
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
        sim.get_result().get_last_value()[INSi] -= migrated[INSi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSCi] += migrated[INSi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISyi] -= migrated[ISyi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISyCi] += migrated[ISyi] * (1 - nondetection);

        //reduce the number of commuters
        migrated[ISyi] *= nondetection;
        migrated[INSi] *= nondetection;
    }
}

} // namespace osecir
} // namespace mio

#endif // ODESECIR_MODEL_H