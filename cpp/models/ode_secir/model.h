/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/populations.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"
#include "memilio/math/smoother.h"
#include "memilio/math/eigen_util.h"

namespace mio
{
namespace osecir
{

// Create template specializations for the age resolved
// SECIHURD model

class Model : public CompartmentalModel<InfectionState, Populations<AgeGroup, InfectionState>, Parameters>
{
    using Base = CompartmentalModel<InfectionState, mio::Populations<AgeGroup, InfectionState>, Parameters>;

public:
    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    Model(int num_agegroups)
        : Model(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
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

        ContactMatrixGroup const& contact_matrix = params.get<ContactPatterns>();

        auto icu_occupancy           = 0.0;
        auto test_and_trace_required = 0.0;
        for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
            auto rateINS = 0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]);
            test_and_trace_required += (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * rateINS *
                                       this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptoms});
            icu_occupancy += this->populations.get_from(pop, {i, InfectionState::InfectedCritical});
        }

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t Si    = this->populations.get_flat_index({i, InfectionState::Susceptible});
            size_t Ei    = this->populations.get_flat_index({i, InfectionState::Exposed});
            size_t INSi  = this->populations.get_flat_index({i, InfectionState::InfectedNoSymptoms});
            size_t ISyi  = this->populations.get_flat_index({i, InfectionState::InfectedSymptoms});
            size_t ISevi = this->populations.get_flat_index({i, InfectionState::InfectedSevere});
            size_t ICri  = this->populations.get_flat_index({i, InfectionState::InfectedCritical});
            size_t Ri    = this->populations.get_flat_index({i, InfectionState::Recovered});
            size_t Di    = this->populations.get_flat_index({i, InfectionState::Dead});

            dydt[Si] = 0;
            dydt[Ei] = 0;

            double rateE =
                1.0 / (2 * params.get<SerialInterval>()[i] - params.get<IncubationTime>()[i]); // R2 = 1/(2SI-TINC)
            double rateINS =
                0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]); // R3 = 1/(2(TINC-SI))

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

                dydt[Si] -= dummy_S;
                dydt[Ei] += dummy_S;
            }

            // ICU capacity shortage is close
            double criticalPerSevereAdjusted =
                smoother_cosine(icu_occupancy, 0.90 * params.get<ICUCapacity>(), params.get<ICUCapacity>(),
                                params.get<CriticalPerSevere>()[i], 0);

            double deathsPerSevereAdjusted = params.get<CriticalPerSevere>()[i] - criticalPerSevereAdjusted;

            dydt[Ei] -= rateE * y[Ei]; // only exchange of E and INS done here
            dydt[INSi] = rateE * y[Ei] - rateINS * y[INSi];
            dydt[ISyi] = (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * rateINS * y[INSi] -
                         (1 / params.get<TimeInfectedSymptoms>()[i]) * y[ISyi];
            dydt[ISevi] = params.get<SeverePerInfectedSymptoms>()[i] / params.get<TimeInfectedSymptoms>()[i] * y[ISyi] -
                          (1 / params.get<TimeInfectedSevere>()[i]) * y[ISevi];
            dydt[ICri] = -(1 / params.get<TimeInfectedCritical>()[i]) * y[ICri];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[ICri] += criticalPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevi];

            dydt[Ri] =
                params.get<RecoveredPerInfectedNoSymptoms>()[i] * rateINS * y[INSi] +
                (1 - params.get<SeverePerInfectedSymptoms>()[i]) / params.get<TimeInfectedSymptoms>()[i] * y[ISyi] +
                (1 - params.get<CriticalPerSevere>()[i]) / params.get<TimeInfectedSevere>()[i] * y[ISevi] +
                (1 - params.get<DeathsPerCritical>()[i]) / params.get<TimeInfectedCritical>()[i] * y[ICri];

            dydt[Di] = params.get<DeathsPerCritical>()[i] / params.get<TimeInfectedCritical>()[i] * y[ICri];
            // add potential, additional deaths due to ICU overflow
            dydt[Di] += deathsPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevi];
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
template <class Base = mio::Simulation<Model>>
class Simulation;

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
 * @tparam Base simulation type that uses a secir compartment model. default mio::Simulation. For testing purposes only!
 */
template <class Base>
class Simulation : public Base
{
public:
    /**
     * construct a simulation.
     * @param model the model to simulate.
     * @param t0 start time
     * @param dt time steps
     */
    Simulation(Model const& model, double t0 = 0., double dt = 0.1)
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
        auto& dyn_npis         = this->get_model().parameters.template get<DynamicNPIsInfectedSymptoms>();
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
 * specialization of simulate for secir models using Simulation.
 * @param t0 start time.
 * @param tmax end time.
 * @param dt time step.
 * @param model secir model to simulate.
 * @param integrator optional integrator, uses rk45 if nullptr.
 */
inline auto simulate(double t0, double tmax, double dt, const Model& model,
                     std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    return mio::simulate<Model, Simulation<>>(t0, tmax, dt, model, integrator);
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

//see declaration above
template <class Base>
double get_reproduction_number(const Simulation<Base>& sim, Eigen::Index timept)
{
    auto const& params = sim.get_model().parameters;
    size_t num_groups = (size_t)sim.get_model().parameters.get_num_groups();
    
    Eigen::MatrixXd F(5*num_groups, 5*num_groups);
    Eigen::MatrixXd V(5*num_groups, 5*num_groups);
    F = Eigen::MatrixXd::Zero(5*num_groups, 5*num_groups);//Initialize matrices F and V with zeroes
    V = Eigen::MatrixXd::Zero(5*num_groups, 5*num_groups);

        auto test_and_trace_required = 0.0;
        auto icu_occupancy = 0.0;
        for (auto i = AgeGroup(0); i < (mio::AgeGroup) num_groups; ++i) {
             auto rateINS = 0.5 / (params.template get<IncubationTime>()[i] - params.template get<SerialInterval>()[i]);
            test_and_trace_required += (1 - params.template get<RecoveredPerInfectedNoSymptoms>()[i]) * rateINS *
                                       sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::InfectedNoSymptoms*(num_groups)+(size_t)i];
            icu_occupancy += sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::InfectedCritical*(num_groups)+(size_t)i];
        }

    double season_val = (1 + params.template get<Seasonality>() *
                             sin(3.141592653589793 * (std::fmod((sim.get_model().parameters.template get<StartDay>() + timept), 365.0) / 182.5 + 0.5)));
    ContactMatrixGroup const& contact_matrix = sim.get_model().parameters.template get<ContactPatterns>();

    Eigen::MatrixXd cont_freq_eff(num_groups, num_groups);
    Eigen::VectorXd divN(num_groups);
    Eigen::VectorXd riskFromInfectedSymptomatic_derivatives(num_groups);
    Eigen::VectorXd riskFromInfectedSymptomatic(num_groups);

    for(mio::AgeGroup k = 0; k < (mio::AgeGroup)num_groups; k++){
        double temp = sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::Susceptible*(num_groups)+(size_t)k] + 
                    sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::Exposed*(num_groups)+(size_t)k] +
                    sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::InfectedNoSymptoms*(num_groups)+(size_t)k] +
                    sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::InfectedSymptoms*(num_groups)+(size_t)k] +
                    sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::InfectedSevere*(num_groups)+(size_t)k]+
                    sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::InfectedCritical*(num_groups)+(size_t)k]+
                    sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::Recovered*(num_groups)+(size_t)k];
        divN[(size_t)k] = 1/temp;

        riskFromInfectedSymptomatic[(size_t)k] = smoother_cosine(
                    test_and_trace_required, (size_t)params.template get<TestAndTraceCapacity>(), (size_t)(params.template get<TestAndTraceCapacity>()) * 5,
                    params.template get<RiskOfInfectionFromSymptomatic>()[k],
                    params.template get<MaxRiskOfInfectionFromSymptomatic>()[k]);

        auto rateINS = 0.5 / (params.template get<IncubationTime>()[k] - params.template get<SerialInterval>()[(mio::AgeGroup)k]);

        if(test_and_trace_required < params.template get<TestAndTraceCapacity>() || test_and_trace_required > 5*params.template get<TestAndTraceCapacity>()){
            riskFromInfectedSymptomatic_derivatives[(size_t)k] = 0;
        }
        else{
            riskFromInfectedSymptomatic_derivatives[(size_t)k] = -0.5*3.14159265358979323846*(params.template get<RiskOfInfectionFromSymptomatic>()[k]-params.template get<MaxRiskOfInfectionFromSymptomatic>()[k])/
                                                        (4*params.template get<TestAndTraceCapacity>())*(1-params.template get<RecoveredPerInfectedNoSymptoms>()[k])*
                                                        rateINS*std::sin(3.14159265358979323846/(4*params.template get<TestAndTraceCapacity>())*(test_and_trace_required - params.template get<TestAndTraceCapacity>()));
        }
        for(Eigen::Index l = 0; l < (Eigen::Index)num_groups; l++){
            cont_freq_eff(l,(size_t)k) = season_val * contact_matrix.get_matrix_at(timept)(static_cast<Eigen::Index>((size_t)l),
                                                                 static_cast<Eigen::Index>((size_t)k));
        }
    }

    //First initialize matrix F

    for(size_t i = 0; i < (size_t)params.get_num_groups(); i++){

        double temp = 0;
        for(Eigen::Index k = 0; k < (Eigen::Index)num_groups; k++){
            temp += cont_freq_eff(i,k)*
            sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::InfectedSymptoms*(num_groups)+(size_t)k]*
            riskFromInfectedSymptomatic_derivatives[k]*
            divN[k];
        }
    
        
        for(size_t j = 0; j < num_groups; j++){//Second block in first row of F (Denoted by A in the latex)
            F(i,j+num_groups) = sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::Susceptible*(num_groups)+(size_t)i]*
                    params.template get<TransmissionProbabilityOnContact>()[(mio::AgeGroup)i];
                    (cont_freq_eff(i,j)*params.template get<RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)j]*divN[(size_t)j]+temp);
        }
    
    

        for(size_t j = 0; j < num_groups; j++){
            F(i,j+2*num_groups) = sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::Susceptible*(num_groups)+(size_t)i]*
                    params.template get<TransmissionProbabilityOnContact>()[(mio::AgeGroup)i]*cont_freq_eff(i,j)*riskFromInfectedSymptomatic[(size_t)j]*divN[(size_t)j];    
        }
    }
    //Now initialize matrix V

    for(Eigen::Index i = 0; i < (Eigen::Index)num_groups; i++){

        double rateE =
                1.0 / (2 * params.template get<SerialInterval>()[(mio::AgeGroup)i] - params.template get<IncubationTime>()[(mio::AgeGroup)i]); // R2 = 1/(2SI-TINC)
        double rateINS =
                0.5 / (params.template get<IncubationTime>()[(mio::AgeGroup)i] - params.template get<SerialInterval>()[(mio::AgeGroup)i]); // R3 = 1/(2(TINC-SI))

        double criticalPerSevereAdjusted =
                smoother_cosine(icu_occupancy, 0.90 * params.template get<ICUCapacity>(), params.template get<ICUCapacity>(),
                                params.template get<CriticalPerSevere>()[(mio::AgeGroup)i], 0);
        
        V(i,i) = rateE;
        V(i+num_groups, i) = -rateE;
        V(i+num_groups,i+num_groups) = rateINS;
        V(i+2*num_groups, i+num_groups) = -(1-(1 - params.template get<RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)i]))*rateINS;
        V(i+2*num_groups, i+2*num_groups) = (1 / params.template get<TimeInfectedSymptoms>()[(mio::AgeGroup)i]);
        V(i+3*num_groups, i+2*num_groups) = -params.template get<SeverePerInfectedSymptoms>()[(mio::AgeGroup)i] / params.template get<TimeInfectedSymptoms>()[(mio::AgeGroup)i];
        V(i+3*num_groups, i+3*num_groups) = 1 /(params.template get<TimeInfectedSevere>()[(mio::AgeGroup)i]);
        V(i+4*num_groups, i+3*num_groups) = - criticalPerSevereAdjusted/(params.template get<TimeInfectedSevere>()[(mio::AgeGroup)i]);
        V(i+4*num_groups, i+4*num_groups) = 1/(params.template get<TimeInfectedCritical>()[(mio::AgeGroup)i]);

        if(!(icu_occupancy < 0.9*params.template get<ICUCapacity>() || icu_occupancy > (double)(params.template get<ICUCapacity>()))){
            for(size_t j = 0; j < num_groups; j++){
                V(i+4*num_groups, j+4*num_groups) -= sim.get_result().get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::InfectedSevere*num_groups +(size_t)i]/
                params.template get<TimeInfectedSevere>()[(mio::AgeGroup)i]*
                0.5*params.template get<CriticalPerSevere>()[(mio::AgeGroup)i]*3.141592653589793/(0.1*params.template get<ICUCapacity>())*
                std::sin(3.141592653589793/(0.1*params.template get<ICUCapacity>())*(icu_occupancy-0.9*params.template get<ICUCapacity>()));
            }
        }
    }
    //Now try to invert the matrix V
        V = V.inverse();

        //Compute F*V
        Eigen::MatrixXd NextGenMatrix(5*num_groups, 5*num_groups);
        NextGenMatrix = F*V;

        //Compute the biggest eigenvalue in absolute value
         Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces;

        ces.compute(NextGenMatrix);
        const Eigen::VectorXcd tempvector = ces.eigenvalues();

        Eigen::VectorXd tempvector1;
        tempvector1.resize(tempvector.size());
        //Do this later with some iterator
        for (int i = 0; i < tempvector.size(); i++) {
             tempvector1[i] = std::abs(tempvector[i]);
        }
        return tempvector1.maxCoeff();
}

//see declaration above
template <class Base>
Eigen::VectorXd get_reproduction_numbers(const Simulation<Base>& sim)
{
    Eigen::VectorXd temp(sim.get_result().get_num_time_points()); 
    for (int i = 0; i < sim.get_result().get_num_time_points(); i++) {
        temp[i] = get_reproduction_number(sim, (Eigen::Index)i);
    }
    return temp;
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
    auto&& t_inc     = params.template get<IncubationTime>().array().template cast<double>();
    auto&& t_ser     = params.template get<SerialInterval>().array().template cast<double>();
    auto&& p_asymp   = params.template get<RecoveredPerInfectedNoSymptoms>().array().template cast<double>();
    auto&& p_inf     = params.template get<RiskOfInfectionFromSymptomatic>().array().template cast<double>();
    auto&& p_inf_max = params.template get<MaxRiskOfInfectionFromSymptomatic>().array().template cast<double>();
    //slice of InfectedNoSymptoms
    auto y_car = slice(y, {Eigen::Index(InfectionState::InfectedNoSymptoms),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)});

    //compute isolation, same as infection risk from main model
    auto R3                          = 0.5 / (t_inc - t_ser);
    auto test_and_trace_required     = ((1 - p_asymp) * R3 * y_car.array()).sum();
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

} // namespace osecir
} // namespace mio

#endif // ODESECIR_MODEL_H
