/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#ifndef SECIR_H
#define SECIR_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/populations.h"
#include "secir/infection_state.h"
#include "secir/secir_params.h"
#include "memilio/math/smoother.h"
#include "memilio/math/eigen_util.h"

namespace mio
{

// Create template specializations for the age resolved
// SECIHURD model

class SecirModel : public CompartmentalModel<InfectionState, Populations<AgeGroup, InfectionState>, SecirParams>
{
    using Base = CompartmentalModel<InfectionState, mio::Populations<AgeGroup, InfectionState>, SecirParams>;
    using Pa   = Base::ParameterSet;
    using Po   = Base::Populations;

public:
    SecirModel(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
#if !USE_DERIV_FUNC
        size_t n_agegroups = (size_t)AgeGroup::Count;
        for (size_t i = 0; i < n_agegroups; i++) {
            for (size_t j = 0; j < n_agegroups; j++) {

                // Si to Ei individually for each age group j
                this->add_flow(
                    std::make_tuple((AgeGroup)i, InfectionState::S), std::make_tuple((AgeGroup)i, InfectionState::E),
                    [i, j](Pa const& p, Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y,
                           double t) {
                        //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
                        auto test_and_trace_required = (1 - params.probabilities[i].get_asymp_per_infectious()) *
                                                       dummy_R3 * Po::get_from(pop, (AgeGroup)i, InfectionState::C);
                        auto risk_from_symptomatic =
                            smoother_cosine(test_and_trace_required, params.get_test_and_trace_capacity(),
                                            params.get_test_and_trace_capacity() * 5,
                                            params.probabilities[i].get_risk_from_symptomatic(),
                                            params.probabilities[i].get_test_and_trace_max_risk_from_symptomatic());

                        // effective contact rate by contact rate between groups i and j and damping j
                        ScalarType season_val =
                            (1 + p.get_seasonality() * sin(3.141592653589793 *
                                                           (std::fmod((p.get_start_day() + t), 365.0) / 182.5 + 0.5)));
                        ScalarType cont_freq_eff =
                            season_val *
                            contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                        ScalarType Nj = Po::get_from(pop, (AgeGroup)j, InfectionState::S) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::E) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::C) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::I) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::H) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::U) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::R); // without died people
                        ScalarType divNj = 1.0 / Nj; // precompute 1.0/Nj
                        ScalarType Si    = Po::get_from(y, (AgeGroup)i, InfectionState::S);
                        ScalarType Cj    = Po::get_from(pop, (AgeGroup)j, InfectionState::C);
                        ScalarType Ij    = Po::get_from(pop, (AgeGroup)j, InfectionState::I);
                        return Si * cont_freq_eff * divNj * p.probabilities[i].get_infection_from_contact() *
                               (p.probabilities[j].get_carrier_infectability() * Cj + risk_from_symptomatic * Ij);
                    });
            }

            // Ei to Ci
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::E),
                           std::make_tuple((AgeGroup)i, InfectionState::C),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::E) /
                                      (2 * p.times[i].get_serialinterval() - p.times[i].get_incubation());
                           });

            // Ci to Ii
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::C),
                           std::make_tuple((AgeGroup)i, InfectionState::I),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               double dummy_R3 = 0.5 / (p.times[i].get_incubation() - p.times[i].get_serialinterval());
                               double alpha    = p.probabilities[i].get_asymp_per_infectious();
                               return ((1 - alpha) * dummy_R3) * Po::get_from(y, (AgeGroup)i, InfectionState::C);
                           });

            // Ci to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::C),
                           std::make_tuple((AgeGroup)i, InfectionState::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               double alpha = p.probabilities[i].get_asymp_per_infectious();
                               return (alpha / p.times[i].get_infectious_asymp()) *
                                      Po::get_from(y, (AgeGroup)i, InfectionState::C);
                           });

            // Ii to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::I),
                           std::make_tuple((AgeGroup)i, InfectionState::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::I) *
                                      (1 - p.probabilities[i].get_hospitalized_per_infectious()) /
                                      p.times[i].get_infectious_mild();
                           });

            // Ii to Hi
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::I), std::make_tuple((AgeGroup)i, InfectionState::H),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y,
                    double /*t*/) {
                    return Po::get_from(y, (AgeGroup)i, InfectionState::I) *
                           p.probabilities[i].get_hospitalized_per_infectious() / p.times[i].get_home_to_hospitalized();
                });

            // Hi to Ui
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::H),
                           std::make_tuple((AgeGroup)i, InfectionState::U),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               ScalarType icu_occupancy = 0;
                               for (size_t j = 0; j < (size_t)AgeGroup::Count; ++j) {
                                   icu_occupancy += Po::get_from(y, (AgeGroup)j, InfectionState::U);
                               }

                               ScalarType prob_hosp2icu =
                                   smoother_cosine(icu_occupancy, 0.90 * p.get_icu_capacity(), p.get_icu_capacity(),
                                                   p.probabilities[i].get_icu_per_hospitalized(), 0);

                               return Po::get_from(y, (AgeGroup)i, InfectionState::H) * prob_hosp2icu /
                                      p.times[i].get_hospitalized_to_icu();
                           });

            // Hi to Di
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::H), std::make_tuple((AgeGroup)i, InfectionState::D),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y,
                    double /*t*/) {
                    ScalarType icu_occupancy = 0;
                    for (size_t j = 0; j < (size_t)AgeGroup::Count; ++j) {
                        icu_occupancy += Po::get_from(y, (AgeGroup)j, InfectionState::U);
                    }

                    ScalarType prob_hosp2icu =
                        smoother_cosine(icu_occupancy, 0.90 * p.get_icu_capacity(), p.get_icu_capacity(),
                                        p.probabilities[i].get_icu_per_hospitalized(), 0);
                    ScalarType prob_hosp2dead = p.probabilities[i].get_icu_per_hospitalized() - prob_hosp2icu;

                    return Po::get_from(y, (AgeGroup)i, InfectionState::H) * prob_hosp2dead /
                           p.times[i].get_hospitalized_to_icu();
                });

            // Hi to Ri
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::H), std::make_tuple((AgeGroup)i, InfectionState::R),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y,
                    double /*t*/) {
                    return Po::get_from(y, (AgeGroup)i, InfectionState::H) *
                           (1 - p.probabilities[i].get_icu_per_hospitalized()) / p.times[i].get_hospitalized_to_home();
                });

            // Ui to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::U),
                           std::make_tuple((AgeGroup)i, InfectionState::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::U) *
                                      (1 - p.probabilities[i].get_dead_per_icu()) / p.times[i].get_icu_to_home();
                           });

            // Ui to Di
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::U),
                           std::make_tuple((AgeGroup)i, InfectionState::D),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::U) *
                                      p.probabilities[i].get_dead_per_icu() / p.times[i].get_icu_to_dead();
                           });
        }
#endif
    }

    SecirModel(int num_agegroups)
        : SecirModel(Po({AgeGroup(num_agegroups), InfectionState::Count}), Pa(AgeGroup(num_agegroups)))
    {
    }

#if USE_DERIV_FUNC

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        // alpha  // percentage of asymptomatic cases
        // beta // risk of infection from the infected symptomatic patients
        // rho   // hospitalized per infectious
        // theta // icu per hospitalized
        // delta  // deaths per ICUs
        // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D
        auto const& params   = this->parameters;
        AgeGroup n_agegroups = params.get_num_groups();

        ContactMatrixGroup const& contact_matrix = params.get<mio::ContactPatterns>();

        auto icu_occupancy           = 0.0;
        auto test_and_trace_required = 0.0;
        for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
            auto dummy_R3 = 0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]);
            test_and_trace_required += (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 *
                                       this->populations.get_from(pop, {i, InfectionState::Carrier});
            icu_occupancy += this->populations.get_from(pop, {i, InfectionState::ICU});
        }

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            size_t Ei = this->populations.get_flat_index({i, InfectionState::Exposed});
            size_t Ci = this->populations.get_flat_index({i, InfectionState::Carrier});
            size_t Ii = this->populations.get_flat_index({i, InfectionState::Infected});
            size_t Hi = this->populations.get_flat_index({i, InfectionState::Hospitalized});
            size_t Ui = this->populations.get_flat_index({i, InfectionState::ICU});
            size_t Ri = this->populations.get_flat_index({i, InfectionState::Recovered});
            size_t Di = this->populations.get_flat_index({i, InfectionState::Dead});

            dydt[Si] = 0;
            dydt[Ei] = 0;

            double dummy_R2 =
                1.0 / (2 * params.get<SerialInterval>()[i] - params.get<IncubationTime>()[i]); // R2 = 1/(2SI-TINC)
            double dummy_R3 =
                0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]); // R3 = 1/(2(TINC-SI))

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {
                size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
                size_t Ej = this->populations.get_flat_index({j, InfectionState::Exposed});
                size_t Cj = this->populations.get_flat_index({j, InfectionState::Carrier});
                size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
                size_t Hj = this->populations.get_flat_index({j, InfectionState::Hospitalized});
                size_t Uj = this->populations.get_flat_index({j, InfectionState::ICU});
                size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

                //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
                auto risk_from_symptomatic = smoother_cosine(
                    test_and_trace_required, params.get<TestAndTraceCapacity>(), params.get<TestAndTraceCapacity>() * 5,
                    params.get<RiskOfInfectionFromSympomatic>()[j], params.get<MaxRiskOfInfectionFromSympomatic>()[j]);

                // effective contact rate by contact rate between groups i and j and damping j
                double season_val = (1 + params.get<mio::Seasonality>() *
                                             sin(3.141592653589793 *
                                                 (std::fmod((params.get<mio::StartDay>() + t), 365.0) / 182.5 + 0.5)));
                double cont_freq_eff =
                    season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                 static_cast<Eigen::Index>((size_t)j));
                double Nj = pop[Sj] + pop[Ej] + pop[Cj] + pop[Ij] + pop[Hj] + pop[Uj] + pop[Rj]; // without died people
                double divNj = 1.0 / Nj; // precompute 1.0/Nj
                double dummy_S =
                    y[Si] * cont_freq_eff * divNj * params.get<InfectionProbabilityFromContact>()[i] *
                    (params.get<RelativeCarrierInfectability>()[j] * pop[Cj] + risk_from_symptomatic * pop[Ij]);

                dydt[Si] -= dummy_S; // -R1*(C+beta*I)*S/N0
                dydt[Ei] += dummy_S; // R1*(C+beta*I)*S/N0-R2*E
            }

            // ICU capacity shortage is close
            double prob_hosp2icu =
                smoother_cosine(icu_occupancy, 0.90 * params.get<mio::ICUCapacity>(), params.get<mio::ICUCapacity>(),
                                params.get<ICUCasesPerHospitalized>()[i], 0);

            double prob_hosp2dead = params.get<ICUCasesPerHospitalized>()[i] - prob_hosp2icu;

            dydt[Ei] -= dummy_R2 * y[Ei]; // only exchange of E and C done here
            dydt[Ci] = dummy_R2 * y[Ei] -
                       ((1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                        params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i]) *
                           y[Ci];
            dydt[Ii] = (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 * y[Ci] -
                       ((1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] +
                        params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i]) *
                           y[Ii];
            dydt[Hi] =
                params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i] * y[Ii] -
                ((1 - params.get<ICUCasesPerHospitalized>()[i]) / params.get<HospitalizedToHomeTime>()[i] +
                 params.get<ICUCasesPerHospitalized>()[i] / params.get<HospitalizedToICUTime>()[i]) *
                    y[Hi];
            dydt[Ui] = -((1 - params.get<DeathsPerICU>()[i]) / params.get<ICUToHomeTime>()[i] +
                         params.get<DeathsPerICU>()[i] / params.get<ICUToDeathTime>()[i]) *
                       y[Ui];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[Ui] += prob_hosp2icu / params.get<HospitalizedToICUTime>()[i] * y[Hi];

            dydt[Ri] =
                params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i] * y[Ci] +
                (1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] * y[Ii] +
                (1 - params.get<ICUCasesPerHospitalized>()[i]) / params.get<HospitalizedToHomeTime>()[i] * y[Hi] +
                (1 - params.get<DeathsPerICU>()[i]) / params.get<ICUToHomeTime>()[i] * y[Ui];

            dydt[Di] = params.get<DeathsPerICU>()[i] / params.get<ICUToDeathTime>()[i] * y[Ui];
            // add potential, additional deaths due to ICU overflow
            dydt[Di] += prob_hosp2dead / params.get<HospitalizedToICUTime>()[i] * y[Hi];
        }
    }

#endif // USE_DERIV_FUNC

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Secir");
        obj.add_element("Parameters", parameters);
        obj.add_element("Populations", populations);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<SecirModel> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Secir");
        auto par = obj.expect_element("Parameters", Tag<ParameterSet>{});
        auto pop = obj.expect_element("Populations", Tag<Populations>{});
        return apply(
            io,
            [](auto&& par_, auto&& pop_) {
                return SecirModel{pop_, par_};
            },
            par, pop);
    }
};

//forward declaration, see below.
template <class Base = Simulation<SecirModel>>
class SecirSimulation;

/**
 * get percentage of infections per total population.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @tparam Base simulation type that uses a secir compartment model. see SecirSimulation.
 */
template <class Base = Simulation<SecirModel>>
double get_infections_relative(const SecirSimulation<Base>& model, double t,
                               const Eigen::Ref<const Eigen::VectorXd>& y);

/**
 * specialization of compartment model simulation for secir models.
 * @tparam Base simulation type that uses a secir compartment model. default mio::SecirSimulation. For testing purposes only!
 */
template <class Base>
class SecirSimulation : public Base
{
public:
    /**
     * construct a simulation.
     * @param model the model to simulate.
     * @param t0 start time
     * @param dt time steps
     */
    SecirSimulation(SecirModel const& model, double t0 = 0., double dt = 0.1)
        : Base(model, t0, dt)
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
        auto& dyn_npis         = this->get_model().parameters.template get<DynamicNPIsInfected>();
        auto& contact_patterns = this->get_model().parameters.template get<ContactPatterns>();
        if (dyn_npis.get_thresholds().size() > 0) {
            auto t        = Base::get_result().get_last_time();
            const auto dt = dyn_npis.get_interval().get();

            while (t < tmax) {
                auto dt_eff = std::min({dt, tmax - t, m_t_last_npi_check + dt - t});

                Base::advance(t + dt_eff);
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
            return Base::advance(tmax);
        }
    }

private:
    double m_t_last_npi_check;
    std::pair<double, SimulationTime> m_dynamic_npi = {-std::numeric_limits<double>::max(), mio::SimulationTime(0)};
};

/**
 * specialization of simulate for secir models using SecirSimulation.
 * @param t0 start time.
 * @param tmax end time.
 * @param dt time step.
 * @param model secir model to simulate.
 * @param integrator optional integrator, uses rk45 if nullptr.
 */
inline auto simulate(double t0, double tmax, double dt, const SecirModel& model,
                     std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    return simulate<SecirModel, SecirSimulation<>>(t0, tmax, dt, model, integrator);
}

//see declaration above.
template <class Base>
double get_infections_relative(const SecirSimulation<Base>& sim, double /*t*/,
                               const Eigen::Ref<const Eigen::VectorXd>& y)
{
    double sum_inf = 0;
    for (auto i = AgeGroup(0); i < sim.get_model().parameters.get_num_groups(); ++i) {
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::Infected});
    }
    auto inf_rel = sum_inf / sim.get_model().populations.get_total();

    return inf_rel;
}

/**
 * Get migration factors.
 * Used by migration graph simulation.
 * Like infection risk, migration of infected individuals is reduced if they are well isolated.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @return vector expression, same size as y, with migration factors per compartment.
 * @tparam Base simulation type that uses a secir compartment model. see SecirSimulation.
 */
template <class Base = Simulation<SecirModel>>
auto get_migration_factors(const SecirSimulation<Base>& sim, double /*t*/, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    auto& params = sim.get_model().parameters;
    //parameters as arrays
    auto&& t_inc     = params.template get<IncubationTime>().array().template cast<double>();
    auto&& t_ser     = params.template get<SerialInterval>().array().template cast<double>();
    auto&& p_asymp   = params.template get<AsymptoticCasesPerInfectious>().array().template cast<double>();
    auto&& p_inf     = params.template get<RiskOfInfectionFromSympomatic>().array().template cast<double>();
    auto&& p_inf_max = params.template get<MaxRiskOfInfectionFromSympomatic>().array().template cast<double>();
    //slice of carriers
    auto y_car = slice(y, {Eigen::Index(InfectionState::Carrier), Eigen::Index(size_t(params.get_num_groups())),
                           Eigen::Index(InfectionState::Count)});

    //compute isolation, same as infection risk from main model
    auto R3                      = 0.5 / (t_inc - t_ser);
    auto test_and_trace_required = ((1 - p_asymp) * R3 * y_car.array()).sum();
    auto test_and_trace_capacity = double(params.template get<TestAndTraceCapacity>());
    auto risk_from_symptomatic   = smoother_cosine(test_and_trace_required, test_and_trace_capacity,
                                                 test_and_trace_capacity * 5, p_inf.matrix(), p_inf_max.matrix());

    //set factor for infected
    auto factors = Eigen::VectorXd::Ones(y.rows()).eval();
    slice(factors, {Eigen::Index(InfectionState::Infected), Eigen::Index(size_t(params.get_num_groups())),
                    Eigen::Index(InfectionState::Count)})
        .array() = risk_from_symptomatic;
    return factors;
}

} // namespace mio

#endif
