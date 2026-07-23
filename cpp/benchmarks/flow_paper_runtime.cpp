/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License")
*/

#include "benchmark/benchmark.h"
#include "memilio/compartments/compartmental_model.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/math/smoother.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/utils/logging.h"
#include "ode_secir/model.h"
#include "ode_secirts/model.h"
#include "ode_seir/model.h"

#include "boost/numeric/odeint/stepper/runge_kutta4.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

using FP = ScalarType;

using SeirState            = mio::oseir::InfectionState;
using SeirPopulations      = mio::Populations<FP, mio::AgeGroup, SeirState>;
using SeirParameters       = mio::oseir::Parameters<FP>;
using SeirCompartmentModel = mio::CompartmentalModel<FP, SeirState, SeirPopulations, SeirParameters>;

using SecirState       = mio::osecir::InfectionState;
using SecirPopulations = mio::Populations<FP, mio::AgeGroup, SecirState>;
using SecirParameters  = mio::osecir::Parameters<FP>;
using SecirFlowModel   = mio::FlowModel<FP, SecirState, SecirPopulations, SecirParameters, mio::osecir::Flows>;

constexpr FP t0   = 0.0;
constexpr FP tmax = 600.0;
constexpr FP dt   = 0.2;

std::unique_ptr<mio::OdeIntegratorCore<FP>> make_fixed_rk4_integrator()
{
    return std::make_unique<mio::ExplicitStepperWrapper<FP, boost::numeric::odeint::runge_kutta4>>();
}

template <class Model>
auto make_compartment_simulation(const Model& model, Eigen::Index output_points)
{
    mio::Simulation<FP, Model> simulation(model, t0, dt);
    simulation.set_integrator_core(make_fixed_rk4_integrator());
    simulation.get_result().reserve(output_points);
    return simulation;
}

template <class Model>
auto make_flow_simulation(const Model& model, Eigen::Index output_points)
{
    mio::FlowSimulation<FP, Model> simulation(model, t0, dt);
    simulation.set_integrator_core(make_fixed_rk4_integrator());
    simulation.get_result().reserve(output_points);
    simulation.get_flows().reserve(output_points);
    return simulation;
}

void setup_runtime_model(mio::oseir::Model<FP>& model)
{
    const FP total_population = 10000.0;
    const auto num_groups     = model.parameters.get_num_groups();
    for (mio::AgeGroup group = 0; group < num_groups; ++group) {
        const auto denominator = static_cast<FP>(static_cast<size_t>(num_groups));
        model.populations[{group, mio::oseir::InfectionState::Exposed}]   = 100.0 / denominator;
        model.populations[{group, mio::oseir::InfectionState::Infected}]  = 100.0 / denominator;
        model.populations[{group, mio::oseir::InfectionState::Recovered}] = 100.0 / denominator;
        model.populations[{group, mio::oseir::InfectionState::Susceptible}] =
            total_population / denominator - model.populations[{group, mio::oseir::InfectionState::Exposed}] -
            model.populations[{group, mio::oseir::InfectionState::Infected}] -
            model.populations[{group, mio::oseir::InfectionState::Recovered}];
    }
    model.parameters.template set<mio::oseir::TimeExposed<FP>>(5.2);
    model.parameters.template set<mio::oseir::TimeInfected<FP>>(6.0);
    model.parameters.template set<mio::oseir::TransmissionProbabilityOnContact<FP>>(0.04);
    model.parameters.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat()[0].get_baseline().setConstant(
        10.0);
}

FP evaluate_seir_infection_rate(const SeirPopulations& populations, const SeirParameters& parameters,
                                Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                                mio::Index<mio::AgeGroup> target_group)
{
    const auto age_groups    = mio::reduce_index<mio::Index<mio::AgeGroup>>(populations.size());
    const size_t susceptible = populations.get_flat_index({target_group, SeirState::Susceptible});
    FP rate                  = 0.0;

    for (auto source_group : mio::make_index_range(age_groups)) {
        const size_t source_susceptible = populations.get_flat_index({source_group, SeirState::Susceptible});
        const size_t source_exposed     = populations.get_flat_index({source_group, SeirState::Exposed});
        const size_t source_infected    = populations.get_flat_index({source_group, SeirState::Infected});
        const size_t source_recovered   = populations.get_flat_index({source_group, SeirState::Recovered});
        const FP source_population =
            pop[source_susceptible] + pop[source_exposed] + pop[source_infected] + pop[source_recovered];
        const FP inverse_population =
            source_population < mio::Limits<FP>::zero_tolerance() ? FP(0.0) : FP(1.0 / source_population);
        const FP coefficient =
            parameters.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(
                mio::SimulationTime<FP>(t))(target_group.get(), source_group.get()) *
            parameters.template get<mio::oseir::TransmissionProbabilityOnContact<FP>>()[target_group] *
            inverse_population;
        rate += coefficient * y[susceptible] * pop[source_infected];
    }
    return rate;
}

enum class SeirRateHandling
{
    Reuse,
    EquationWise
};

template <SeirRateHandling RateHandling>
class DirectSeirModel : public SeirCompartmentModel
{
public:
    DirectSeirModel(const SeirPopulations& populations, const SeirParameters& parameters)
        : SeirCompartmentModel(populations, parameters)
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                         Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
    {
        const auto age_groups  = mio::reduce_index<mio::Index<mio::AgeGroup>>(this->populations.size());
        const auto& parameters = this->parameters;

        for (auto group : mio::make_index_range(age_groups)) {
            const size_t susceptible = this->populations.get_flat_index({group, SeirState::Susceptible});
            const size_t exposed     = this->populations.get_flat_index({group, SeirState::Exposed});
            const size_t infected    = this->populations.get_flat_index({group, SeirState::Infected});
            const size_t recovered   = this->populations.get_flat_index({group, SeirState::Recovered});

            if constexpr (RateHandling == SeirRateHandling::Reuse) {
                const FP infection = evaluate_seir_infection_rate(this->populations, parameters, pop, y, t, group);
                const FP progression =
                    (1.0 / parameters.template get<mio::oseir::TimeExposed<FP>>()[group]) * y[exposed];
                const FP recovery =
                    (1.0 / parameters.template get<mio::oseir::TimeInfected<FP>>()[group]) * y[infected];

                dydt[susceptible] -= infection;
                dydt[exposed] += infection;
                dydt[exposed] -= progression;
                dydt[infected] += progression;
                dydt[infected] -= recovery;
                dydt[recovered] += recovery;
            }
            else {
                dydt[susceptible] -= evaluate_seir_infection_rate(this->populations, parameters, pop, y, t, group);
                dydt[exposed] += evaluate_seir_infection_rate(this->populations, parameters, pop, y, t, group);
                dydt[exposed] -= (1.0 / parameters.template get<mio::oseir::TimeExposed<FP>>()[group]) * y[exposed];
                dydt[infected] += (1.0 / parameters.template get<mio::oseir::TimeExposed<FP>>()[group]) * y[exposed];
                dydt[infected] -= (1.0 / parameters.template get<mio::oseir::TimeInfected<FP>>()[group]) * y[infected];
                dydt[recovered] += (1.0 / parameters.template get<mio::oseir::TimeInfected<FP>>()[group]) * y[infected];
            }
        }
    }
};

using ReusedRateSeirModel       = DirectSeirModel<SeirRateHandling::Reuse>;
using EquationWiseRateSeirModel = DirectSeirModel<SeirRateHandling::EquationWise>;

void setup_runtime_model(mio::osecir::Model<FP>& model)
{
    using IS = mio::osecir::InfectionState;

    const FP total_population = 10000.0;
    const auto num_groups     = model.parameters.get_num_groups();
    model.populations.set_total(total_population);
    model.parameters.template get<mio::osecir::ContactPatterns<FP>>().get_cont_freq_mat()[0].get_baseline().setConstant(
        8.0);

    for (mio::AgeGroup group = 0; group < num_groups; ++group) {
        const auto denominator    = static_cast<FP>(static_cast<size_t>(num_groups));
        const FP group_population = total_population / denominator;

        model.parameters.template get<mio::osecir::TimeExposed<FP>>()[group]                       = 3.2;
        model.parameters.template get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[group]            = 2.0;
        model.parameters.template get<mio::osecir::TimeInfectedSymptoms<FP>>()[group]              = 6.0;
        model.parameters.template get<mio::osecir::TimeInfectedSevere<FP>>()[group]                = 8.0;
        model.parameters.template get<mio::osecir::TimeInfectedCritical<FP>>()[group]              = 7.0;
        model.parameters.template get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[group]  = 0.04;
        model.parameters.template get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[group]    = 0.7;
        model.parameters.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group]    = 0.25;
        model.parameters.template get<mio::osecir::RiskOfInfectionFromSymptomatic<FP>>()[group]    = 0.2;
        model.parameters.template get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>>()[group] = 0.4;
        model.parameters.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group]         = 0.08;
        model.parameters.template get<mio::osecir::CriticalPerSevere<FP>>()[group]                 = 0.2;
        model.parameters.template get<mio::osecir::DeathsPerSevere<FP>>()[group]                   = 0.01;
        model.parameters.template get<mio::osecir::DeathsPerCritical<FP>>()[group]                 = 0.25;

        model.populations[{group, IS::Exposed}]                     = 40.0 / denominator;
        model.populations[{group, IS::InfectedNoSymptoms}]          = 25.0 / denominator;
        model.populations[{group, IS::InfectedNoSymptomsConfirmed}] = 0.0;
        model.populations[{group, IS::InfectedSymptoms}]            = 20.0 / denominator;
        model.populations[{group, IS::InfectedSymptomsConfirmed}]   = 0.0;
        model.populations[{group, IS::InfectedSevere}]              = 4.0 / denominator;
        model.populations[{group, IS::InfectedCritical}]            = 1.0 / denominator;
        model.populations[{group, IS::Recovered}]                   = 250.0 / denominator;
        model.populations[{group, IS::Dead}]                        = 0.0;

        const FP assigned =
            model.populations[{group, IS::Exposed}] + model.populations[{group, IS::InfectedNoSymptoms}] +
            model.populations[{group, IS::InfectedNoSymptomsConfirmed}] +
            model.populations[{group, IS::InfectedSymptoms}] +
            model.populations[{group, IS::InfectedSymptomsConfirmed}] + model.populations[{group, IS::InfectedSevere}] +
            model.populations[{group, IS::InfectedCritical}] + model.populations[{group, IS::Recovered}] +
            model.populations[{group, IS::Dead}];
        model.populations[{group, IS::Susceptible}] = group_population - assigned;
    }
    model.parameters.template get<mio::osecir::ICUCapacity<FP>>()          = 1000.0;
    model.parameters.template get<mio::osecir::TestAndTraceCapacity<FP>>() = 10000.0;
    model.apply_constraints();
}

void setup_runtime_model(mio::osecirts::Model<FP>& model)
{
    using IS = mio::osecirts::InfectionState;

    const FP total_population = 10000.0;
    const auto num_groups     = model.parameters.get_num_groups();
    model.populations.set_total(total_population);
    model.parameters.template get<mio::osecirts::ContactPatterns<FP>>()
        .get_cont_freq_mat()[0]
        .get_baseline()
        .setConstant(8.0);

    constexpr size_t num_days = 700;
    model.parameters.template get<mio::osecirts::DailyPartialVaccinations<FP>>().resize(mio::SimulationDay(num_days));
    model.parameters.template get<mio::osecirts::DailyFullVaccinations<FP>>().resize(mio::SimulationDay(num_days));
    model.parameters.template get<mio::osecirts::DailyBoosterVaccinations<FP>>().resize(mio::SimulationDay(num_days));

    for (mio::AgeGroup group = 0; group < num_groups; ++group) {
        const auto denominator    = static_cast<FP>(static_cast<size_t>(num_groups));
        const FP group_population = total_population / denominator;

        model.parameters.template get<mio::osecirts::TimeExposed<FP>>()[group]                           = 3.2;
        model.parameters.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[group]                = 2.0;
        model.parameters.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[group]                  = 6.0;
        model.parameters.template get<mio::osecirts::TimeInfectedSevere<FP>>()[group]                    = 8.0;
        model.parameters.template get<mio::osecirts::TimeInfectedCritical<FP>>()[group]                  = 7.0;
        model.parameters.template get<mio::osecirts::TimeWaningPartialImmunity<FP>>()[group]             = 180.0;
        model.parameters.template get<mio::osecirts::TimeWaningImprovedImmunity<FP>>()[group]            = 240.0;
        model.parameters.template get<mio::osecirts::TimeTemporaryImmunityPI<FP>>()[group]               = 120.0;
        model.parameters.template get<mio::osecirts::TimeTemporaryImmunityII<FP>>()[group]               = 180.0;
        model.parameters.template get<mio::osecirts::TransmissionProbabilityOnContact<FP>>()[group]      = 0.04;
        model.parameters.template get<mio::osecirts::RelativeTransmissionNoSymptoms<FP>>()[group]        = 0.7;
        model.parameters.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[group]        = 0.25;
        model.parameters.template get<mio::osecirts::RiskOfInfectionFromSymptomatic<FP>>()[group]        = 0.2;
        model.parameters.template get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<FP>>()[group]     = 0.4;
        model.parameters.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[group]             = 0.08;
        model.parameters.template get<mio::osecirts::CriticalPerSevere<FP>>()[group]                     = 0.2;
        model.parameters.template get<mio::osecirts::DeathsPerSevere<FP>>()[group]                       = 0.01;
        model.parameters.template get<mio::osecirts::DeathsPerCritical<FP>>()[group]                     = 0.25;
        model.parameters.template get<mio::osecirts::ReducExposedPartialImmunity<FP>>()[group]           = 0.65;
        model.parameters.template get<mio::osecirts::ReducExposedImprovedImmunity<FP>>()[group]          = 0.35;
        model.parameters.template get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<FP>>()[group]  = 0.75;
        model.parameters.template get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<FP>>()[group] = 0.45;
        model.parameters.template get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[group] =
            0.55;
        model.parameters.template get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[group] =
            0.25;
        model.parameters.template get<mio::osecirts::ReducTimeInfectedMild<FP>>()[group]    = 0.8;
        model.parameters.template get<mio::osecirts::InfectiousnessNewVariant<FP>>()[group] = 1.0;

        for (auto day = mio::SimulationDay(0); day < mio::SimulationDay(num_days); ++day) {
            model.parameters.template get<mio::osecirts::DailyPartialVaccinations<FP>>()[{group, day}] = 0.0;
            model.parameters.template get<mio::osecirts::DailyFullVaccinations<FP>>()[{group, day}]    = 0.0;
            model.parameters.template get<mio::osecirts::DailyBoosterVaccinations<FP>>()[{group, day}] = 0.0;
        }

        model.populations[{group, IS::ExposedNaive}]                    = 40.0 / denominator;
        model.populations[{group, IS::InfectedNoSymptomsNaive}]         = 25.0 / denominator;
        model.populations[{group, IS::InfectedSymptomsNaive}]           = 20.0 / denominator;
        model.populations[{group, IS::InfectedSevereNaive}]             = 4.0 / denominator;
        model.populations[{group, IS::InfectedCriticalNaive}]           = 1.0 / denominator;
        model.populations[{group, IS::TemporaryImmunePartialImmunity}]  = 100.0 / denominator;
        model.populations[{group, IS::TemporaryImmuneImprovedImmunity}] = 150.0 / denominator;

        const FP assigned = model.populations[{group, IS::ExposedNaive}] +
                            model.populations[{group, IS::InfectedNoSymptomsNaive}] +
                            model.populations[{group, IS::InfectedSymptomsNaive}] +
                            model.populations[{group, IS::InfectedSevereNaive}] +
                            model.populations[{group, IS::InfectedCriticalNaive}] +
                            model.populations[{group, IS::TemporaryImmunePartialImmunity}] +
                            model.populations[{group, IS::TemporaryImmuneImprovedImmunity}];
        model.populations[{group, IS::SusceptibleNaive}] = group_population - assigned;
    }
    model.parameters.template get<mio::osecirts::ICUCapacity<FP>>()          = 1000.0;
    model.parameters.template get<mio::osecirts::TestAndTraceCapacity<FP>>() = 10000.0;
    model.parameters.template get<mio::osecirts::StartDay<FP>>()             = 0.0;
    model.parameters.template get<mio::osecirts::StartDayNewVariant<FP>>()   = std::numeric_limits<FP>::max();
    model.parameters.apply_constraints();
}

template <class Model>
Model make_runtime_model(int num_age_groups)
{
    Model model(num_age_groups);
    setup_runtime_model(model);
    return model;
}

template <class Model>
class AuxiliaryTrackingModel
{
public:
    AuxiliaryTrackingModel(Model model, size_t num_auxiliary_states)
        : m_model(std::move(model))
        , m_num_auxiliary_states(num_auxiliary_states)
        , m_compartment_dim(static_cast<size_t>(m_model.get_initial_values().size()))
        , m_initial(Eigen::VectorX<FP>::Zero(static_cast<Eigen::Index>(m_compartment_dim + m_num_auxiliary_states)))
        , m_flow_buffer(Eigen::VectorX<FP>::Zero(m_model.get_initial_flows().size()))
    {
        m_initial.head(static_cast<Eigen::Index>(m_compartment_dim)) = m_model.get_initial_values();
    }

    Eigen::VectorX<FP> get_initial_values() const
    {
        return m_initial;
    }

    bool check_constraints() const
    {
        return m_model.check_constraints();
    }

    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                              Eigen::Ref<Eigen::VectorX<FP>> dydt) const
    {
        dydt.setZero();
        const auto compartment_dim = static_cast<Eigen::Index>(m_compartment_dim);

        m_flow_buffer.setZero();
        m_model.get_flows(pop.head(compartment_dim), y.head(compartment_dim), t, m_flow_buffer);
        m_model.get_derivatives(m_flow_buffer, dydt.head(compartment_dim));

        for (size_t index = 0; index < m_num_auxiliary_states; ++index) {
            dydt[static_cast<Eigen::Index>(m_compartment_dim + index)] =
                m_flow_buffer[static_cast<Eigen::Index>(index)];
        }
    }

private:
    Model m_model;
    size_t m_num_auxiliary_states;
    size_t m_compartment_dim;
    Eigen::VectorX<FP> m_initial;
    mutable Eigen::VectorX<FP> m_flow_buffer;
};

template <class Model>
class PrimaryFlowTrackingModel
{
public:
    PrimaryFlowTrackingModel(Model model)
        : m_model(std::move(model))
        , m_initial_populations(m_model.get_initial_values())
        , m_initial_flows(m_model.get_initial_flows())
        , m_populations(Eigen::VectorX<FP>::Zero(m_initial_populations.size()))
        , m_flow_delta(Eigen::VectorX<FP>::Zero(m_initial_flows.size()))
    {
    }

    Eigen::VectorX<FP> get_initial_values() const
    {
        return m_initial_flows;
    }

    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorX<FP>>, Eigen::Ref<const Eigen::VectorX<FP>> flows, FP t,
                              Eigen::Ref<Eigen::VectorX<FP>> dflows_dt) const
    {
        m_flow_delta = flows - m_initial_flows;
        m_model.get_derivatives(m_flow_delta, m_populations);
        m_populations += m_initial_populations;
        dflows_dt.setZero();
        m_model.get_flows(m_populations, m_populations, t, dflows_dt);
    }

private:
    Model m_model;
    Eigen::VectorX<FP> m_initial_populations;
    Eigen::VectorX<FP> m_initial_flows;
    mutable Eigen::VectorX<FP> m_populations;
    mutable Eigen::VectorX<FP> m_flow_delta;
};

enum class SecirTransitionKind
{
    Infection,
    ExposedToNoSymptoms,
    NoSymptomsToSymptoms,
    NoSymptomsToRecovered,
    ConfirmedNoSymptomsToConfirmedSymptoms,
    ConfirmedNoSymptomsToRecovered,
    SymptomsToSevere,
    SymptomsToRecovered,
    ConfirmedSymptomsToSevere,
    ConfirmedSymptomsToRecovered,
    SevereToCritical,
    SevereToRecovered,
    SevereToDead,
    CriticalToDead,
    CriticalToRecovered
};

struct SecirTrackedTransition {
    mio::AgeGroup group;
    SecirTransitionKind kind;
};

constexpr std::array<SecirTransitionKind, 13> secir_auxiliary_transition_kinds{
    SecirTransitionKind::Infection,
    SecirTransitionKind::ExposedToNoSymptoms,
    SecirTransitionKind::NoSymptomsToSymptoms,
    SecirTransitionKind::NoSymptomsToRecovered,
    SecirTransitionKind::ConfirmedNoSymptomsToConfirmedSymptoms,
    SecirTransitionKind::ConfirmedNoSymptomsToRecovered,
    SecirTransitionKind::SymptomsToSevere,
    SecirTransitionKind::SymptomsToRecovered,
    SecirTransitionKind::ConfirmedSymptomsToSevere,
    SecirTransitionKind::ConfirmedSymptomsToRecovered,
    SecirTransitionKind::SevereToRecovered,
    SecirTransitionKind::CriticalToDead,
    SecirTransitionKind::CriticalToRecovered,
};

enum class SecirParameterization
{
    Times,
    Rates
};

struct SecirInverseTimes {
    explicit SecirInverseTimes(const SecirParameters& parameters)
        : exposed(parameters.get_num_groups().get())
        , no_symptoms(parameters.get_num_groups().get())
        , symptoms(parameters.get_num_groups().get())
        , severe(parameters.get_num_groups().get())
        , critical(parameters.get_num_groups().get())
    {
        for (mio::AgeGroup group = 0; group < parameters.get_num_groups(); ++group) {
            const auto index   = static_cast<Eigen::Index>((size_t)group);
            exposed[index]     = 1.0 / parameters.template get<mio::osecir::TimeExposed<FP>>()[group];
            no_symptoms[index] = 1.0 / parameters.template get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[group];
            symptoms[index]    = 1.0 / parameters.template get<mio::osecir::TimeInfectedSymptoms<FP>>()[group];
            severe[index]      = 1.0 / parameters.template get<mio::osecir::TimeInfectedSevere<FP>>()[group];
            critical[index]    = 1.0 / parameters.template get<mio::osecir::TimeInfectedCritical<FP>>()[group];
        }
    }

    Eigen::VectorX<FP> exposed;
    Eigen::VectorX<FP> no_symptoms;
    Eigen::VectorX<FP> symptoms;
    Eigen::VectorX<FP> severe;
    Eigen::VectorX<FP> critical;
};

template <SecirParameterization Parameterization>
class BenchmarkSecirFlowModel : public SecirFlowModel
{
public:
    using Base = SecirFlowModel;

    BenchmarkSecirFlowModel(const SecirPopulations& populations, const SecirParameters& parameters)
        : Base(populations, parameters)
        , m_inverse_times(parameters)
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const auto& params                                = this->parameters;
        const auto num_groups                             = params.get_num_groups();
        const mio::ContactMatrixGroup<FP>& contact_matrix = params.template get<mio::osecir::ContactPatterns<FP>>();

        FP icu_occupancy           = 0.0;
        FP test_and_trace_required = 0.0;
        for (mio::AgeGroup group = 0; group < num_groups; ++group) {
            const auto no_symptoms = this->populations.get_flat_index({group, SecirState::InfectedNoSymptoms});
            const auto critical    = this->populations.get_flat_index({group, SecirState::InfectedCritical});
            test_and_trace_required +=
                (1.0 - params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group]) *
                no_symptoms_rate(group) * pop[no_symptoms];
            icu_occupancy += pop[critical];
        }

        for (mio::AgeGroup group = 0; group < num_groups; ++group) {
            const auto susceptible = this->populations.get_flat_index({group, SecirState::Susceptible});
            const auto exposed     = this->populations.get_flat_index({group, SecirState::Exposed});
            const auto no_symptoms = this->populations.get_flat_index({group, SecirState::InfectedNoSymptoms});
            const auto confirmed_no_symptoms =
                this->populations.get_flat_index({group, SecirState::InfectedNoSymptomsConfirmed});
            const auto symptoms = this->populations.get_flat_index({group, SecirState::InfectedSymptoms});
            const auto confirmed_symptoms =
                this->populations.get_flat_index({group, SecirState::InfectedSymptomsConfirmed});
            const auto severe   = this->populations.get_flat_index({group, SecirState::InfectedSevere});
            const auto critical = this->populations.get_flat_index({group, SecirState::InfectedCritical});

            for (mio::AgeGroup source_group = 0; source_group < num_groups; ++source_group) {
                const auto source_susceptible =
                    this->populations.get_flat_index({source_group, SecirState::Susceptible});
                const auto source_exposed = this->populations.get_flat_index({source_group, SecirState::Exposed});
                const auto source_no_symptoms =
                    this->populations.get_flat_index({source_group, SecirState::InfectedNoSymptoms});
                const auto source_symptoms =
                    this->populations.get_flat_index({source_group, SecirState::InfectedSymptoms});
                const auto source_severe = this->populations.get_flat_index({source_group, SecirState::InfectedSevere});
                const auto source_critical =
                    this->populations.get_flat_index({source_group, SecirState::InfectedCritical});
                const auto source_recovered = this->populations.get_flat_index({source_group, SecirState::Recovered});
                const FP source_population  = pop[source_susceptible] + pop[source_exposed] + pop[source_no_symptoms] +
                                             pop[source_symptoms] + pop[source_severe] + pop[source_critical] +
                                             pop[source_recovered];
                const FP inverse_population =
                    source_population < mio::Limits<FP>::zero_tolerance() ? FP(0.0) : FP(1.0 / source_population);
                const FP symptomatic_risk = mio::smoother_cosine<FP>(
                    test_and_trace_required, params.template get<mio::osecir::TestAndTraceCapacity<FP>>(),
                    params.template get<mio::osecir::TestAndTraceCapacity<FP>>() *
                        params.template get<mio::osecir::TestAndTraceCapacityMaxRisk<FP>>(),
                    params.template get<mio::osecir::RiskOfInfectionFromSymptomatic<FP>>()[source_group],
                    params.template get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>>()[source_group]);
                const FP seasonality =
                    1.0 + params.template get<mio::osecir::Seasonality<FP>>() *
                              std::sin(std::numbers::pi_v<FP> *
                                       ((params.template get<mio::osecir::StartDay<FP>>() + t) / 182.5 + 0.5));
                const FP contact_rate = seasonality * contact_matrix.get_matrix_at(mio::SimulationTime<FP>(t))(
                                                          static_cast<Eigen::Index>((size_t)group),
                                                          static_cast<Eigen::Index>((size_t)source_group));
                flows[this->template get_flat_flow_index<SecirState::Susceptible, SecirState::Exposed>({group})] +=
                    y[susceptible] * contact_rate * inverse_population *
                    params.template get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[group] *
                    (params.template get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[source_group] *
                         pop[source_no_symptoms] +
                     symptomatic_risk * pop[source_symptoms]);
            }

            const FP adjusted_critical =
                mio::smoother_cosine<FP>(icu_occupancy, 0.9 * params.template get<mio::osecir::ICUCapacity<FP>>(),
                                         params.template get<mio::osecir::ICUCapacity<FP>>(),
                                         params.template get<mio::osecir::CriticalPerSevere<FP>>()[group], 0.0);
            const FP adjusted_deaths =
                params.template get<mio::osecir::CriticalPerSevere<FP>>()[group] - adjusted_critical;
            const FP recovered_no_symptoms =
                params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group];
            const FP severe_symptoms = params.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group];
            const FP critical_severe = params.template get<mio::osecir::CriticalPerSevere<FP>>()[group];
            const FP deaths_severe   = params.template get<mio::osecir::DeathsPerSevere<FP>>()[group];
            const FP deaths_critical = params.template get<mio::osecir::DeathsPerCritical<FP>>()[group];

            flows[this->template get_flat_flow_index<SecirState::Exposed, SecirState::InfectedNoSymptoms>({group})] =
                exposed_rate(group) * y[exposed];
            flows[this->template get_flat_flow_index<SecirState::InfectedNoSymptoms, SecirState::InfectedSymptoms>(
                {group})] = (1.0 - recovered_no_symptoms) * no_symptoms_rate(group) * y[no_symptoms];
            flows[this->template get_flat_flow_index<SecirState::InfectedNoSymptoms, SecirState::Recovered>({group})] =
                recovered_no_symptoms * no_symptoms_rate(group) * y[no_symptoms];
            flows[this->template get_flat_flow_index<SecirState::InfectedNoSymptomsConfirmed,
                                                     SecirState::InfectedSymptomsConfirmed>({group})] =
                (1.0 - recovered_no_symptoms) * no_symptoms_rate(group) * y[confirmed_no_symptoms];
            flows[this->template get_flat_flow_index<SecirState::InfectedNoSymptomsConfirmed, SecirState::Recovered>(
                {group})] = recovered_no_symptoms * no_symptoms_rate(group) * y[confirmed_no_symptoms];
            flows[this->template get_flat_flow_index<SecirState::InfectedSymptoms, SecirState::InfectedSevere>(
                {group})] = severe_symptoms * symptoms_rate(group) * y[symptoms];
            flows[this->template get_flat_flow_index<SecirState::InfectedSymptoms, SecirState::Recovered>({group})] =
                (1.0 - severe_symptoms) * symptoms_rate(group) * y[symptoms];
            flows[this->template get_flat_flow_index<SecirState::InfectedSymptomsConfirmed, SecirState::InfectedSevere>(
                {group})] = severe_symptoms * symptoms_rate(group) * y[confirmed_symptoms];
            flows[this->template get_flat_flow_index<SecirState::InfectedSymptomsConfirmed, SecirState::Recovered>(
                {group})] = (1.0 - severe_symptoms) * symptoms_rate(group) * y[confirmed_symptoms];
            flows[this->template get_flat_flow_index<SecirState::InfectedSevere, SecirState::InfectedCritical>(
                {group})] = adjusted_critical * severe_rate(group) * y[severe];
            flows[this->template get_flat_flow_index<SecirState::InfectedSevere, SecirState::Recovered>({group})] =
                (1.0 - critical_severe - deaths_severe) * severe_rate(group) * y[severe];
            flows[this->template get_flat_flow_index<SecirState::InfectedSevere, SecirState::Dead>({group})] =
                (deaths_severe + adjusted_deaths) * severe_rate(group) * y[severe];
            flows[this->template get_flat_flow_index<SecirState::InfectedCritical, SecirState::Dead>({group})] =
                deaths_critical * critical_rate(group) * y[critical];
            flows[this->template get_flat_flow_index<SecirState::InfectedCritical, SecirState::Recovered>({group})] =
                (1.0 - deaths_critical) * critical_rate(group) * y[critical];
        }
    }

    FP get_transition_rate(SecirTransitionKind kind, mio::AgeGroup group, Eigen::Ref<const Eigen::VectorX<FP>> pop,
                           Eigen::Ref<const Eigen::VectorX<FP>> y, FP t) const
    {
        const auto& params     = this->parameters;
        const auto state_value = [&](SecirState state) {
            return y[this->populations.get_flat_index({group, state})];
        };

        switch (kind) {
        case SecirTransitionKind::Infection:
            return infection_rate(group, pop, y, t);
        case SecirTransitionKind::ExposedToNoSymptoms:
            return exposed_rate(group) * state_value(SecirState::Exposed);
        case SecirTransitionKind::NoSymptomsToSymptoms:
            return (1.0 - params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group]) *
                   no_symptoms_rate(group) * state_value(SecirState::InfectedNoSymptoms);
        case SecirTransitionKind::NoSymptomsToRecovered:
            return params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group] *
                   no_symptoms_rate(group) * state_value(SecirState::InfectedNoSymptoms);
        case SecirTransitionKind::ConfirmedNoSymptomsToConfirmedSymptoms:
            return (1.0 - params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group]) *
                   no_symptoms_rate(group) * state_value(SecirState::InfectedNoSymptomsConfirmed);
        case SecirTransitionKind::ConfirmedNoSymptomsToRecovered:
            return params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group] *
                   no_symptoms_rate(group) * state_value(SecirState::InfectedNoSymptomsConfirmed);
        case SecirTransitionKind::SymptomsToSevere:
            return params.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group] * symptoms_rate(group) *
                   state_value(SecirState::InfectedSymptoms);
        case SecirTransitionKind::SymptomsToRecovered:
            return (1.0 - params.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group]) *
                   symptoms_rate(group) * state_value(SecirState::InfectedSymptoms);
        case SecirTransitionKind::ConfirmedSymptomsToSevere:
            return params.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group] * symptoms_rate(group) *
                   state_value(SecirState::InfectedSymptomsConfirmed);
        case SecirTransitionKind::ConfirmedSymptomsToRecovered:
            return (1.0 - params.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group]) *
                   symptoms_rate(group) * state_value(SecirState::InfectedSymptomsConfirmed);
        case SecirTransitionKind::SevereToCritical:
            return adjusted_critical_per_severe(group, pop) * severe_rate(group) *
                   state_value(SecirState::InfectedSevere);
        case SecirTransitionKind::SevereToRecovered:
            return (1.0 - params.template get<mio::osecir::CriticalPerSevere<FP>>()[group] -
                    params.template get<mio::osecir::DeathsPerSevere<FP>>()[group]) *
                   severe_rate(group) * state_value(SecirState::InfectedSevere);
        case SecirTransitionKind::SevereToDead: {
            const FP adjusted_critical = adjusted_critical_per_severe(group, pop);
            const FP adjusted_deaths =
                params.template get<mio::osecir::CriticalPerSevere<FP>>()[group] - adjusted_critical;
            return (params.template get<mio::osecir::DeathsPerSevere<FP>>()[group] + adjusted_deaths) *
                   severe_rate(group) * state_value(SecirState::InfectedSevere);
        }
        case SecirTransitionKind::CriticalToDead:
            return params.template get<mio::osecir::DeathsPerCritical<FP>>()[group] * critical_rate(group) *
                   state_value(SecirState::InfectedCritical);
        case SecirTransitionKind::CriticalToRecovered:
            return (1.0 - params.template get<mio::osecir::DeathsPerCritical<FP>>()[group]) * critical_rate(group) *
                   state_value(SecirState::InfectedCritical);
        }
        throw std::runtime_error("Unhandled SECIR transition kind.");
    }

private:
    FP exposed_rate(mio::AgeGroup group) const
    {
        if constexpr (Parameterization == SecirParameterization::Rates) {
            return m_inverse_times.exposed[static_cast<Eigen::Index>((size_t)group)];
        }
        return 1.0 / this->parameters.template get<mio::osecir::TimeExposed<FP>>()[group];
    }

    FP no_symptoms_rate(mio::AgeGroup group) const
    {
        if constexpr (Parameterization == SecirParameterization::Rates) {
            return m_inverse_times.no_symptoms[static_cast<Eigen::Index>((size_t)group)];
        }
        return 1.0 / this->parameters.template get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[group];
    }

    FP symptoms_rate(mio::AgeGroup group) const
    {
        if constexpr (Parameterization == SecirParameterization::Rates) {
            return m_inverse_times.symptoms[static_cast<Eigen::Index>((size_t)group)];
        }
        return 1.0 / this->parameters.template get<mio::osecir::TimeInfectedSymptoms<FP>>()[group];
    }

    FP severe_rate(mio::AgeGroup group) const
    {
        if constexpr (Parameterization == SecirParameterization::Rates) {
            return m_inverse_times.severe[static_cast<Eigen::Index>((size_t)group)];
        }
        return 1.0 / this->parameters.template get<mio::osecir::TimeInfectedSevere<FP>>()[group];
    }

    FP critical_rate(mio::AgeGroup group) const
    {
        if constexpr (Parameterization == SecirParameterization::Rates) {
            return m_inverse_times.critical[static_cast<Eigen::Index>((size_t)group)];
        }
        return 1.0 / this->parameters.template get<mio::osecir::TimeInfectedCritical<FP>>()[group];
    }

    FP infection_rate(mio::AgeGroup group, Eigen::Ref<const Eigen::VectorX<FP>> pop,
                      Eigen::Ref<const Eigen::VectorX<FP>> y, FP t) const
    {
        const auto& params         = this->parameters;
        const auto num_groups      = params.get_num_groups();
        FP test_and_trace_required = 0.0;
        for (mio::AgeGroup source_group = 0; source_group < num_groups; ++source_group) {
            const auto no_symptoms = this->populations.get_flat_index({source_group, SecirState::InfectedNoSymptoms});
            test_and_trace_required +=
                (1.0 - params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[source_group]) *
                no_symptoms_rate(source_group) * pop[no_symptoms];
        }

        const auto susceptible = this->populations.get_flat_index({group, SecirState::Susceptible});
        const FP seasonality =
            1.0 + params.template get<mio::osecir::Seasonality<FP>>() *
                      std::sin(std::numbers::pi_v<FP> *
                               ((params.template get<mio::osecir::StartDay<FP>>() + t) / 182.5 + 0.5));
        const mio::ContactMatrixGroup<FP>& contact_matrix = params.template get<mio::osecir::ContactPatterns<FP>>();
        FP rate                                           = 0.0;
        for (mio::AgeGroup source_group = 0; source_group < num_groups; ++source_group) {
            const auto source_susceptible = this->populations.get_flat_index({source_group, SecirState::Susceptible});
            const auto source_exposed     = this->populations.get_flat_index({source_group, SecirState::Exposed});
            const auto source_no_symptoms =
                this->populations.get_flat_index({source_group, SecirState::InfectedNoSymptoms});
            const auto source_symptoms = this->populations.get_flat_index({source_group, SecirState::InfectedSymptoms});
            const auto source_severe   = this->populations.get_flat_index({source_group, SecirState::InfectedSevere});
            const auto source_critical = this->populations.get_flat_index({source_group, SecirState::InfectedCritical});
            const auto source_recovered = this->populations.get_flat_index({source_group, SecirState::Recovered});
            const FP source_population  = pop[source_susceptible] + pop[source_exposed] + pop[source_no_symptoms] +
                                         pop[source_symptoms] + pop[source_severe] + pop[source_critical] +
                                         pop[source_recovered];
            const FP inverse_population =
                source_population < mio::Limits<FP>::zero_tolerance() ? FP(0.0) : FP(1.0 / source_population);
            const FP symptomatic_risk = mio::smoother_cosine<FP>(
                test_and_trace_required, params.template get<mio::osecir::TestAndTraceCapacity<FP>>(),
                params.template get<mio::osecir::TestAndTraceCapacity<FP>>() *
                    params.template get<mio::osecir::TestAndTraceCapacityMaxRisk<FP>>(),
                params.template get<mio::osecir::RiskOfInfectionFromSymptomatic<FP>>()[source_group],
                params.template get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>>()[source_group]);
            const FP contact_rate = seasonality * contact_matrix.get_matrix_at(mio::SimulationTime<FP>(t))(
                                                      static_cast<Eigen::Index>((size_t)group),
                                                      static_cast<Eigen::Index>((size_t)source_group));
            rate += y[susceptible] * contact_rate * inverse_population *
                    params.template get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[group] *
                    (params.template get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[source_group] *
                         pop[source_no_symptoms] +
                     symptomatic_risk * pop[source_symptoms]);
        }
        return rate;
    }

    FP adjusted_critical_per_severe(mio::AgeGroup group, Eigen::Ref<const Eigen::VectorX<FP>> pop) const
    {
        FP icu_occupancy = 0.0;
        for (mio::AgeGroup source_group = 0; source_group < this->parameters.get_num_groups(); ++source_group) {
            icu_occupancy += pop[this->populations.get_flat_index({source_group, SecirState::InfectedCritical})];
        }
        return mio::smoother_cosine<FP>(
            icu_occupancy, 0.9 * this->parameters.template get<mio::osecir::ICUCapacity<FP>>(),
            this->parameters.template get<mio::osecir::ICUCapacity<FP>>(),
            this->parameters.template get<mio::osecir::CriticalPerSevere<FP>>()[group], 0.0);
    }

    SecirInverseTimes m_inverse_times;
};

template <SecirParameterization Parameterization>
class EquationWiseSecirModel
{
public:
    explicit EquationWiseSecirModel(BenchmarkSecirFlowModel<Parameterization> model)
        : m_model(std::move(model))
    {
    }

    Eigen::VectorX<FP> get_initial_values() const
    {
        return m_model.get_initial_values();
    }

    bool check_constraints() const
    {
        return m_model.check_constraints();
    }

    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                              Eigen::Ref<Eigen::VectorX<FP>> dydt) const
    {
        dydt.setZero();
        for (mio::AgeGroup group = 0; group < m_model.parameters.get_num_groups(); ++group) {
            add_transition<SecirState::Susceptible, SecirState::Exposed>(SecirTransitionKind::Infection, group, pop, y,
                                                                         t, dydt);
            add_transition<SecirState::Exposed, SecirState::InfectedNoSymptoms>(
                SecirTransitionKind::ExposedToNoSymptoms, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedNoSymptoms, SecirState::InfectedSymptoms>(
                SecirTransitionKind::NoSymptomsToSymptoms, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedNoSymptoms, SecirState::Recovered>(
                SecirTransitionKind::NoSymptomsToRecovered, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedNoSymptomsConfirmed, SecirState::InfectedSymptomsConfirmed>(
                SecirTransitionKind::ConfirmedNoSymptomsToConfirmedSymptoms, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedNoSymptomsConfirmed, SecirState::Recovered>(
                SecirTransitionKind::ConfirmedNoSymptomsToRecovered, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedSymptoms, SecirState::InfectedSevere>(
                SecirTransitionKind::SymptomsToSevere, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedSymptoms, SecirState::Recovered>(
                SecirTransitionKind::SymptomsToRecovered, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedSymptomsConfirmed, SecirState::InfectedSevere>(
                SecirTransitionKind::ConfirmedSymptomsToSevere, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedSymptomsConfirmed, SecirState::Recovered>(
                SecirTransitionKind::ConfirmedSymptomsToRecovered, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedSevere, SecirState::InfectedCritical>(
                SecirTransitionKind::SevereToCritical, group, pop, y, t, dydt);
            add_transition<SecirState::InfectedSevere, SecirState::Recovered>(SecirTransitionKind::SevereToRecovered,
                                                                              group, pop, y, t, dydt);
            add_transition<SecirState::InfectedSevere, SecirState::Dead>(SecirTransitionKind::SevereToDead, group, pop,
                                                                         y, t, dydt);
            add_transition<SecirState::InfectedCritical, SecirState::Dead>(SecirTransitionKind::CriticalToDead, group,
                                                                           pop, y, t, dydt);
            add_transition<SecirState::InfectedCritical, SecirState::Recovered>(
                SecirTransitionKind::CriticalToRecovered, group, pop, y, t, dydt);
        }
    }

private:
    template <SecirState Source, SecirState Target>
    void add_transition(SecirTransitionKind kind, mio::AgeGroup group, Eigen::Ref<const Eigen::VectorX<FP>> pop,
                        Eigen::Ref<const Eigen::VectorX<FP>> y, FP t, Eigen::Ref<Eigen::VectorX<FP>> dydt) const
    {
        const auto source = m_model.populations.get_flat_index({group, Source});
        const auto target = m_model.populations.get_flat_index({group, Target});
        dydt[source] -= m_model.get_transition_rate(kind, group, pop, y, t);
        dydt[target] += m_model.get_transition_rate(kind, group, pop, y, t);
    }

    BenchmarkSecirFlowModel<Parameterization> m_model;
};

std::vector<SecirTrackedTransition> make_secir_transition_order(size_t num_age_groups,
                                                                const std::string& selection_order)
{
    std::vector<SecirTrackedTransition> transitions;
    transitions.reserve(secir_auxiliary_transition_kinds.size() * num_age_groups);

    const auto append_infections = [&]() {
        for (size_t group = 0; group < num_age_groups; ++group) {
            transitions.push_back({mio::AgeGroup(group), SecirTransitionKind::Infection});
        }
    };
    const auto append_remaining_group_wise = [&]() {
        for (size_t group = 0; group < num_age_groups; ++group) {
            for (const auto kind : secir_auxiliary_transition_kinds) {
                if (kind != SecirTransitionKind::Infection) {
                    transitions.push_back({mio::AgeGroup(group), kind});
                }
            }
        }
    };

    if (selection_order == "group_wise") {
        for (size_t group = 0; group < num_age_groups; ++group) {
            for (const auto kind : secir_auxiliary_transition_kinds) {
                transitions.push_back({mio::AgeGroup(group), kind});
            }
        }
    }
    else if (selection_order == "infection_first") {
        append_infections();
        append_remaining_group_wise();
    }
    else if (selection_order == "other_first") {
        append_remaining_group_wise();
        append_infections();
    }
    else {
        throw std::runtime_error("Unknown SECIR transition selection order: " + selection_order);
    }
    return transitions;
}

size_t get_secir_auxiliary_transition_rank(SecirTransitionKind kind)
{
    const auto position =
        std::find(secir_auxiliary_transition_kinds.begin(), secir_auxiliary_transition_kinds.end(), kind);
    if (position == secir_auxiliary_transition_kinds.end()) {
        throw std::runtime_error("SECIR transition is not part of the auxiliary benchmark.");
    }
    return static_cast<size_t>(std::distance(secir_auxiliary_transition_kinds.begin(), position));
}

void sort_secir_transitions_canonically(std::vector<SecirTrackedTransition>& transitions)
{
    std::sort(transitions.begin(), transitions.end(), [](const auto& left, const auto& right) {
        const auto left_group  = static_cast<size_t>(left.group);
        const auto right_group = static_cast<size_t>(right.group);
        if (left_group != right_group) {
            return left_group < right_group;
        }
        return get_secir_auxiliary_transition_rank(left.kind) < get_secir_auxiliary_transition_rank(right.kind);
    });
}

std::vector<SecirTrackedTransition> select_secir_transitions(size_t num_age_groups, const std::string& selection_order,
                                                             size_t tracked_transitions)
{
    auto transitions = make_secir_transition_order(num_age_groups, selection_order);
    if (tracked_transitions > transitions.size()) {
        throw std::runtime_error("Invalid SECIR auxiliary transition selection.");
    }
    transitions.erase(transitions.begin() + static_cast<std::ptrdiff_t>(tracked_transitions), transitions.end());

    // The preference order determines the selected set only. A fixed evaluation
    // order prevents branch-prediction effects from being mistaken for rate cost.
    sort_secir_transitions_canonically(transitions);
    return transitions;
}

FP recompute_secir_infection_rate(const mio::osecir::Model<FP>& model, mio::AgeGroup group,
                                  Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y,
                                  FP t)
{
    using IS                  = mio::osecir::InfectionState;
    const auto& params        = model.parameters;
    const auto num_age_groups = params.get_num_groups();

    FP test_and_trace_required = 0.0;
    for (mio::AgeGroup source_group = 0; source_group < num_age_groups; ++source_group) {
        const auto infected_no_symptoms = model.populations.get_flat_index({source_group, IS::InfectedNoSymptoms});
        test_and_trace_required +=
            (1.0 - params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[source_group]) /
            params.template get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[source_group] * pop[infected_no_symptoms];
    }

    const auto susceptible = model.populations.get_flat_index({group, IS::Susceptible});
    const FP seasonality   = 1.0 + params.template get<mio::osecir::Seasonality<FP>>() *
                                     std::sin(std::numbers::pi_v<FP> *
                                              ((params.template get<mio::osecir::StartDay<FP>>() + t) / 182.5 + 0.5));
    const mio::ContactMatrixGroup<FP>& contact_matrix = params.template get<mio::osecir::ContactPatterns<FP>>();

    FP infection_rate = 0.0;
    for (mio::AgeGroup source_group = 0; source_group < num_age_groups; ++source_group) {
        const auto source_susceptible = model.populations.get_flat_index({source_group, IS::Susceptible});
        const auto source_exposed     = model.populations.get_flat_index({source_group, IS::Exposed});
        const auto source_no_symptoms = model.populations.get_flat_index({source_group, IS::InfectedNoSymptoms});
        const auto source_symptoms    = model.populations.get_flat_index({source_group, IS::InfectedSymptoms});
        const auto source_severe      = model.populations.get_flat_index({source_group, IS::InfectedSevere});
        const auto source_critical    = model.populations.get_flat_index({source_group, IS::InfectedCritical});
        const auto source_recovered   = model.populations.get_flat_index({source_group, IS::Recovered});
        const FP source_population    = pop[source_susceptible] + pop[source_exposed] + pop[source_no_symptoms] +
                                     pop[source_symptoms] + pop[source_severe] + pop[source_critical] +
                                     pop[source_recovered];
        const FP inverse_population =
            source_population < mio::Limits<FP>::zero_tolerance() ? FP(0.0) : FP(1.0 / source_population);
        const FP symptomatic_risk = mio::smoother_cosine<FP>(
            test_and_trace_required, params.template get<mio::osecir::TestAndTraceCapacity<FP>>(),
            params.template get<mio::osecir::TestAndTraceCapacity<FP>>() *
                params.template get<mio::osecir::TestAndTraceCapacityMaxRisk<FP>>(),
            params.template get<mio::osecir::RiskOfInfectionFromSymptomatic<FP>>()[source_group],
            params.template get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>>()[source_group]);
        const FP contact_rate = seasonality * contact_matrix.get_matrix_at(mio::SimulationTime<FP>(t))(
                                                  static_cast<Eigen::Index>((size_t)group),
                                                  static_cast<Eigen::Index>((size_t)source_group));
        infection_rate += y[susceptible] * contact_rate * inverse_population *
                          params.template get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[group] *
                          (params.template get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[source_group] *
                               pop[source_no_symptoms] +
                           symptomatic_risk * pop[source_symptoms]);
    }
    return infection_rate;
}

FP recompute_secir_critical_per_severe(const mio::osecir::Model<FP>& model, mio::AgeGroup group,
                                       Eigen::Ref<const Eigen::VectorX<FP>> pop)
{
    using IS           = mio::osecir::InfectionState;
    const auto& params = model.parameters;
    FP icu_occupancy   = 0.0;
    for (mio::AgeGroup source_group = 0; source_group < params.get_num_groups(); ++source_group) {
        icu_occupancy += pop[model.populations.get_flat_index({source_group, IS::InfectedCritical})];
    }
    return mio::smoother_cosine<FP>(icu_occupancy, 0.9 * params.template get<mio::osecir::ICUCapacity<FP>>(),
                                    params.template get<mio::osecir::ICUCapacity<FP>>(),
                                    params.template get<mio::osecir::CriticalPerSevere<FP>>()[group], 0.0);
}

FP recompute_secir_transition_rate(const mio::osecir::Model<FP>& model, const SecirTrackedTransition& transition,
                                   Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y,
                                   FP t)
{
    using IS           = mio::osecir::InfectionState;
    const auto group   = transition.group;
    const auto& params = model.parameters;

    const auto state_value = [&](IS state) {
        return y[model.populations.get_flat_index({group, state})];
    };

    switch (transition.kind) {
    case SecirTransitionKind::Infection:
        return recompute_secir_infection_rate(model, group, pop, y, t);
    case SecirTransitionKind::ExposedToNoSymptoms:
        return state_value(IS::Exposed) / params.template get<mio::osecir::TimeExposed<FP>>()[group];
    case SecirTransitionKind::NoSymptomsToSymptoms:
        return (1.0 - params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group]) /
               params.template get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[group] *
               state_value(IS::InfectedNoSymptoms);
    case SecirTransitionKind::NoSymptomsToRecovered:
        return params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group] /
               params.template get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[group] *
               state_value(IS::InfectedNoSymptoms);
    case SecirTransitionKind::ConfirmedNoSymptomsToConfirmedSymptoms:
        return (1.0 - params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group]) /
               params.template get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[group] *
               state_value(IS::InfectedNoSymptomsConfirmed);
    case SecirTransitionKind::ConfirmedNoSymptomsToRecovered:
        return params.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group] /
               params.template get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[group] *
               state_value(IS::InfectedNoSymptomsConfirmed);
    case SecirTransitionKind::SymptomsToSevere:
        return params.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group] /
               params.template get<mio::osecir::TimeInfectedSymptoms<FP>>()[group] * state_value(IS::InfectedSymptoms);
    case SecirTransitionKind::SymptomsToRecovered:
        return (1.0 - params.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group]) /
               params.template get<mio::osecir::TimeInfectedSymptoms<FP>>()[group] * state_value(IS::InfectedSymptoms);
    case SecirTransitionKind::ConfirmedSymptomsToSevere:
        return params.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group] /
               params.template get<mio::osecir::TimeInfectedSymptoms<FP>>()[group] *
               state_value(IS::InfectedSymptomsConfirmed);
    case SecirTransitionKind::ConfirmedSymptomsToRecovered:
        return (1.0 - params.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group]) /
               params.template get<mio::osecir::TimeInfectedSymptoms<FP>>()[group] *
               state_value(IS::InfectedSymptomsConfirmed);
    case SecirTransitionKind::SevereToCritical:
        return recompute_secir_critical_per_severe(model, group, pop) /
               params.template get<mio::osecir::TimeInfectedSevere<FP>>()[group] * state_value(IS::InfectedSevere);
    case SecirTransitionKind::SevereToRecovered:
        return (1.0 - params.template get<mio::osecir::CriticalPerSevere<FP>>()[group] -
                params.template get<mio::osecir::DeathsPerSevere<FP>>()[group]) /
               params.template get<mio::osecir::TimeInfectedSevere<FP>>()[group] * state_value(IS::InfectedSevere);
    case SecirTransitionKind::SevereToDead: {
        const FP adjusted_critical = recompute_secir_critical_per_severe(model, group, pop);
        const FP adjusted_deaths = params.template get<mio::osecir::CriticalPerSevere<FP>>()[group] - adjusted_critical;
        return (params.template get<mio::osecir::DeathsPerSevere<FP>>()[group] + adjusted_deaths) /
               params.template get<mio::osecir::TimeInfectedSevere<FP>>()[group] * state_value(IS::InfectedSevere);
    }
    case SecirTransitionKind::CriticalToDead:
        return params.template get<mio::osecir::DeathsPerCritical<FP>>()[group] /
               params.template get<mio::osecir::TimeInfectedCritical<FP>>()[group] * state_value(IS::InfectedCritical);
    case SecirTransitionKind::CriticalToRecovered:
        return (1.0 - params.template get<mio::osecir::DeathsPerCritical<FP>>()[group]) /
               params.template get<mio::osecir::TimeInfectedCritical<FP>>()[group] * state_value(IS::InfectedCritical);
    }
    throw std::runtime_error("Unhandled SECIR transition kind.");
}

size_t get_secir_transition_index(const mio::osecir::Model<FP>& model, const SecirTrackedTransition& transition)
{
    using IS = mio::osecir::InfectionState;
    switch (transition.kind) {
    case SecirTransitionKind::Infection:
        return model.template get_flat_flow_index<IS::Susceptible, IS::Exposed>({transition.group});
    case SecirTransitionKind::ExposedToNoSymptoms:
        return model.template get_flat_flow_index<IS::Exposed, IS::InfectedNoSymptoms>({transition.group});
    case SecirTransitionKind::NoSymptomsToSymptoms:
        return model.template get_flat_flow_index<IS::InfectedNoSymptoms, IS::InfectedSymptoms>({transition.group});
    case SecirTransitionKind::NoSymptomsToRecovered:
        return model.template get_flat_flow_index<IS::InfectedNoSymptoms, IS::Recovered>({transition.group});
    case SecirTransitionKind::ConfirmedNoSymptomsToConfirmedSymptoms:
        return model.template get_flat_flow_index<IS::InfectedNoSymptomsConfirmed, IS::InfectedSymptomsConfirmed>(
            {transition.group});
    case SecirTransitionKind::ConfirmedNoSymptomsToRecovered:
        return model.template get_flat_flow_index<IS::InfectedNoSymptomsConfirmed, IS::Recovered>({transition.group});
    case SecirTransitionKind::SymptomsToSevere:
        return model.template get_flat_flow_index<IS::InfectedSymptoms, IS::InfectedSevere>({transition.group});
    case SecirTransitionKind::SymptomsToRecovered:
        return model.template get_flat_flow_index<IS::InfectedSymptoms, IS::Recovered>({transition.group});
    case SecirTransitionKind::ConfirmedSymptomsToSevere:
        return model.template get_flat_flow_index<IS::InfectedSymptomsConfirmed, IS::InfectedSevere>(
            {transition.group});
    case SecirTransitionKind::ConfirmedSymptomsToRecovered:
        return model.template get_flat_flow_index<IS::InfectedSymptomsConfirmed, IS::Recovered>({transition.group});
    case SecirTransitionKind::SevereToCritical:
        return model.template get_flat_flow_index<IS::InfectedSevere, IS::InfectedCritical>({transition.group});
    case SecirTransitionKind::SevereToRecovered:
        return model.template get_flat_flow_index<IS::InfectedSevere, IS::Recovered>({transition.group});
    case SecirTransitionKind::SevereToDead:
        return model.template get_flat_flow_index<IS::InfectedSevere, IS::Dead>({transition.group});
    case SecirTransitionKind::CriticalToDead:
        return model.template get_flat_flow_index<IS::InfectedCritical, IS::Dead>({transition.group});
    case SecirTransitionKind::CriticalToRecovered:
        return model.template get_flat_flow_index<IS::InfectedCritical, IS::Recovered>({transition.group});
    }
    throw std::runtime_error("Unhandled SECIR transition kind.");
}

void validate_recomputed_secir_rates(const mio::osecir::Model<FP>& model)
{
    const auto state = model.get_initial_values();
    const auto transitions =
        make_secir_transition_order(static_cast<size_t>(model.parameters.get_num_groups()), "group_wise");
    for (const FP time : {FP(0.0), FP(30.0), FP(90.0)}) {
        auto flow_rates = model.get_initial_flows();
        model.get_flows(state, state, time, flow_rates);
        for (const auto& transition : transitions) {
            const FP expected   = flow_rates[static_cast<Eigen::Index>(get_secir_transition_index(model, transition))];
            const FP recomputed = recompute_secir_transition_rate(model, transition, state, state, time);
            const FP tolerance  = 1e-12 * std::max(FP(1.0), std::abs(expected));
            if (std::abs(expected - recomputed) > tolerance) {
                throw std::runtime_error("Recomputed SECIR transition rate does not match the model rate.");
            }
        }
    }
}

class RecomputedAuxiliarySecirModel
{
public:
    RecomputedAuxiliarySecirModel(mio::osecir::Model<FP> model, std::vector<SecirTrackedTransition> tracked_transitions)
        : m_model(std::move(model))
        , m_tracked_transitions(std::move(tracked_transitions))
        , m_compartment_dim(static_cast<size_t>(m_model.get_initial_values().size()))
        , m_initial(
              Eigen::VectorX<FP>::Zero(static_cast<Eigen::Index>(m_compartment_dim + m_tracked_transitions.size())))
        , m_flow_buffer(Eigen::VectorX<FP>::Zero(m_model.get_initial_flows().size()))
    {
        m_initial.head(static_cast<Eigen::Index>(m_compartment_dim)) = m_model.get_initial_values();
    }

    Eigen::VectorX<FP> get_initial_values() const
    {
        return m_initial;
    }

    bool check_constraints() const
    {
        return m_model.check_constraints();
    }

    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                              Eigen::Ref<Eigen::VectorX<FP>> dydt) const
    {
        dydt.setZero();
        const auto compartment_dim  = static_cast<Eigen::Index>(m_compartment_dim);
        const auto pop_compartments = pop.head(compartment_dim);
        const auto y_compartments   = y.head(compartment_dim);

        m_flow_buffer.setZero();
        m_model.get_flows(pop_compartments, y_compartments, t, m_flow_buffer);
        m_model.get_derivatives(m_flow_buffer, dydt.head(compartment_dim));

        for (size_t index = 0; index < m_tracked_transitions.size(); ++index) {
            dydt[static_cast<Eigen::Index>(m_compartment_dim + index)] = recompute_secir_transition_rate(
                m_model, m_tracked_transitions[index], pop_compartments, y_compartments, t);
        }
    }

private:
    mio::osecir::Model<FP> m_model;
    std::vector<SecirTrackedTransition> m_tracked_transitions;
    size_t m_compartment_dim;
    Eigen::VectorX<FP> m_initial;
    mutable Eigen::VectorX<FP> m_flow_buffer;
};

struct DirectComparisonDiagnostics {
    size_t time_points;
    FP max_abs_difference;
    FP max_population_fraction;
};

template <class Model>
DirectComparisonDiagnostics validate_direct_comparison(const Model& model)
{
    const auto output_points    = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;
    auto compartment_simulation = make_compartment_simulation(model, output_points);
    auto flow_simulation        = make_flow_simulation(model, output_points);
    compartment_simulation.advance(tmax);
    flow_simulation.advance(tmax);

    const auto& compartment_result = compartment_simulation.get_result();
    const auto& flow_result        = flow_simulation.get_result();
    if (compartment_result.get_num_time_points() != flow_result.get_num_time_points()) {
        throw std::runtime_error("Direct comparison produced different output grids.");
    }

    FP max_difference = 0.0;
    for (Eigen::Index index = 0; index < compartment_result.get_num_time_points(); ++index) {
        max_difference = std::max(
            max_difference, (compartment_result.get_value(index) - flow_result.get_value(index)).cwiseAbs().maxCoeff());
    }
    const FP population = model.get_initial_values().sum();
    if (max_difference > 1e-7) {
        throw std::runtime_error("Direct compartment and flow results differ beyond the validation tolerance.");
    }
    return {static_cast<size_t>(compartment_result.get_num_time_points()), max_difference, max_difference / population};
}

DirectComparisonDiagnostics validate_seir_rate_reuse_comparison(const mio::oseir::Model<FP>& flow_model,
                                                                const ReusedRateSeirModel& reused_model,
                                                                const EquationWiseRateSeirModel& equation_wise_model)
{
    const auto output_points      = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;
    auto flow_simulation          = make_flow_simulation(flow_model, output_points);
    auto reused_simulation        = make_compartment_simulation(reused_model, output_points);
    auto equation_wise_simulation = make_compartment_simulation(equation_wise_model, output_points);
    flow_simulation.advance(tmax);
    reused_simulation.advance(tmax);
    equation_wise_simulation.advance(tmax);

    const auto& reference_result     = flow_simulation.get_result();
    const auto& reused_result        = reused_simulation.get_result();
    const auto& equation_wise_result = equation_wise_simulation.get_result();
    if (reference_result.get_num_time_points() != reused_result.get_num_time_points() ||
        reference_result.get_num_time_points() != equation_wise_result.get_num_time_points()) {
        throw std::runtime_error("SEIR rate-reuse comparison produced different output grids.");
    }

    FP max_difference = 0.0;
    for (Eigen::Index index = 0; index < reference_result.get_num_time_points(); ++index) {
        max_difference = std::max(
            max_difference, (reused_result.get_value(index) - reference_result.get_value(index)).cwiseAbs().maxCoeff());
        max_difference =
            std::max(max_difference,
                     (equation_wise_result.get_value(index) - reference_result.get_value(index)).cwiseAbs().maxCoeff());
    }
    const FP population = flow_model.get_initial_values().sum();
    if (max_difference > 1e-7) {
        throw std::runtime_error("SEIR rate-reuse implementations differ beyond the validation tolerance.");
    }
    return {static_cast<size_t>(reference_result.get_num_time_points()), max_difference, max_difference / population};
}

DirectComparisonDiagnostics validate_secir_implementation_comparison(
    const mio::osecir::Model<FP>& reference_model,
    const BenchmarkSecirFlowModel<SecirParameterization::Times>& time_model,
    const EquationWiseSecirModel<SecirParameterization::Times>& time_equation_wise_model,
    const BenchmarkSecirFlowModel<SecirParameterization::Rates>& rate_model,
    const EquationWiseSecirModel<SecirParameterization::Rates>& rate_equation_wise_model)
{
    const auto output_points           = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;
    auto reference_simulation          = make_compartment_simulation(reference_model, output_points);
    auto time_flow_simulation          = make_flow_simulation(time_model, output_points);
    auto time_cached_simulation        = make_compartment_simulation(time_model, output_points);
    auto time_equation_wise_simulation = make_compartment_simulation(time_equation_wise_model, output_points);
    auto rate_flow_simulation          = make_flow_simulation(rate_model, output_points);
    auto rate_cached_simulation        = make_compartment_simulation(rate_model, output_points);
    auto rate_equation_wise_simulation = make_compartment_simulation(rate_equation_wise_model, output_points);

    reference_simulation.advance(tmax);
    time_flow_simulation.advance(tmax);
    time_cached_simulation.advance(tmax);
    time_equation_wise_simulation.advance(tmax);
    rate_flow_simulation.advance(tmax);
    rate_cached_simulation.advance(tmax);
    rate_equation_wise_simulation.advance(tmax);

    const auto& reference = reference_simulation.get_result();
    const std::array<const mio::TimeSeries<FP>*, 6> comparisons{
        &time_flow_simulation.get_result(),          &time_cached_simulation.get_result(),
        &time_equation_wise_simulation.get_result(), &rate_flow_simulation.get_result(),
        &rate_cached_simulation.get_result(),        &rate_equation_wise_simulation.get_result(),
    };
    for (const auto* result : comparisons) {
        if (result->get_num_time_points() != reference.get_num_time_points()) {
            throw std::runtime_error("SECIR implementation comparison produced different output grids.");
        }
    }

    FP max_difference = 0.0;
    for (Eigen::Index index = 0; index < reference.get_num_time_points(); ++index) {
        for (const auto* result : comparisons) {
            max_difference =
                std::max(max_difference, (result->get_value(index) - reference.get_value(index)).cwiseAbs().maxCoeff());
        }
    }
    const FP population = time_model.get_initial_values().sum();
    if (max_difference > 1e-7) {
        throw std::runtime_error("SECIR implementation variants differ beyond the validation tolerance.");
    }
    return {static_cast<size_t>(reference.get_num_time_points()), max_difference, max_difference / population};
}

template <class Model>
void validate_reused_auxiliary_comparison(const Model& base_model)
{
    const auto output_points   = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;
    const auto compartment_dim = static_cast<Eigen::Index>(base_model.get_initial_values().size());
    const auto flow_dim        = static_cast<Eigen::Index>(base_model.get_initial_flows().size());
    const AuxiliaryTrackingModel<Model> auxiliary_model(base_model, static_cast<size_t>(flow_dim));

    auto auxiliary_simulation = make_compartment_simulation(auxiliary_model, output_points);
    auto flow_simulation      = make_flow_simulation(base_model, output_points);
    auxiliary_simulation.advance(tmax);
    flow_simulation.advance(tmax);

    const auto& auxiliary_result   = auxiliary_simulation.get_result();
    const auto& compartment_result = flow_simulation.get_result();
    const auto& flow_result        = flow_simulation.get_flows();
    if (auxiliary_result.get_num_time_points() != compartment_result.get_num_time_points() ||
        auxiliary_result.get_num_time_points() != flow_result.get_num_time_points()) {
        throw std::runtime_error("Auxiliary and flow simulations produced different output grids.");
    }

    FP max_compartment_difference = 0.0;
    FP max_flow_difference        = 0.0;
    for (Eigen::Index index = 0; index < auxiliary_result.get_num_time_points(); ++index) {
        const auto& auxiliary_value = auxiliary_result.get_value(index);
        max_compartment_difference  = std::max(
            max_compartment_difference,
            (auxiliary_value.head(compartment_dim) - compartment_result.get_value(index)).cwiseAbs().maxCoeff());
        max_flow_difference = std::max(
            max_flow_difference, (auxiliary_value.tail(flow_dim) - flow_result.get_value(index)).cwiseAbs().maxCoeff());
    }
    if (max_compartment_difference > 1e-7 || max_flow_difference > 1e-7) {
        throw std::runtime_error("Auxiliary and flow results differ beyond the validation tolerance.");
    }
}

template <class Model>
void validate_primary_flow_comparison(const Model& base_model)
{
    const auto output_points = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;
    const PrimaryFlowTrackingModel<Model> primary_flow_model(base_model);
    auto primary_simulation = make_compartment_simulation(primary_flow_model, output_points);
    auto flow_simulation    = make_flow_simulation(base_model, output_points);
    primary_simulation.advance(tmax);
    flow_simulation.advance(tmax);

    const auto& primary_result = primary_simulation.get_result();
    const auto& flow_result    = flow_simulation.get_flows();
    if (primary_result.get_num_time_points() != flow_result.get_num_time_points()) {
        throw std::runtime_error("Primary-state and complete flow simulations produced different output grids.");
    }
    FP max_difference = 0.0;
    for (Eigen::Index index = 0; index < primary_result.get_num_time_points(); ++index) {
        max_difference = std::max(
            max_difference, (primary_result.get_value(index) - flow_result.get_value(index)).cwiseAbs().maxCoeff());
    }
    if (max_difference > 1e-7) {
        throw std::runtime_error("Primary-state and complete flow results differ beyond the validation tolerance.");
    }
}

template <class Model>
void benchmark_compartment_trajectory(benchmark::State& state, const Model& model, size_t compartment_dim,
                                      size_t flow_dim, size_t tracked_transitions)
{
    const auto output_points = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;
    FP checksum              = 0.0;

    for (auto _ : state) {
        state.PauseTiming();
        {
            auto simulation = make_compartment_simulation(model, output_points);
            state.ResumeTiming();
            simulation.advance(tmax);
            state.PauseTiming();
            checksum += simulation.get_result().get_last_value().sum();
            benchmark::DoNotOptimize(checksum);
            benchmark::ClobberMemory();
        }
        state.ResumeTiming();
    }

    state.counters["N_C"]                 = static_cast<double>(compartment_dim);
    state.counters["N_T"]                 = static_cast<double>(flow_dim);
    state.counters["integrated_dim"]      = static_cast<double>(model.get_initial_values().size());
    state.counters["tracked_transitions"] = static_cast<double>(tracked_transitions);
    state.counters["time_points"]         = static_cast<double>(output_points);
    state.SetItemsProcessed(state.iterations());
}

template <class Model>
void benchmark_final_state_only(benchmark::State& state, const Model& model, size_t compartment_dim, size_t flow_dim)
{
    const auto num_steps = static_cast<size_t>(std::llround((tmax - t0) / dt));
    FP checksum          = 0.0;

    for (auto _ : state) {
        state.PauseTiming();
        auto integrator                              = make_fixed_rk4_integrator();
        auto current                                 = model.get_initial_values();
        Eigen::VectorX<FP> next                      = Eigen::VectorX<FP>::Zero(current.size());
        const mio::DerivFunction<FP> right_hand_side = [&model](auto&& y, auto&& t, auto&& dydt) {
            model.eval_right_hand_side(y, y, t, dydt);
        };
        FP time      = t0;
        FP step_size = dt;
        state.ResumeTiming();
        for (size_t step = 0; step < num_steps; ++step) {
            integrator->step(right_hand_side, current, time, step_size, next);
            current.swap(next);
        }
        state.PauseTiming();
        checksum += current.sum();
        benchmark::DoNotOptimize(checksum);
        benchmark::ClobberMemory();
        state.ResumeTiming();
    }

    state.counters["N_C"]            = static_cast<double>(compartment_dim);
    state.counters["N_T"]            = static_cast<double>(flow_dim);
    state.counters["integrated_dim"] = static_cast<double>(model.get_initial_values().size());
    state.counters["steps"]          = static_cast<double>(num_steps);
    state.SetItemsProcessed(state.iterations());
}

template <class Model>
void benchmark_compartment_simulation(benchmark::State& state, int num_age_groups,
                                      const DirectComparisonDiagnostics& diagnostics)
{
    const auto model = make_runtime_model<Model>(num_age_groups);
    benchmark_compartment_trajectory(state, model, static_cast<size_t>(model.get_initial_values().size()),
                                     static_cast<size_t>(model.get_initial_flows().size()), 0);
    state.counters["max_abs_difference"]      = static_cast<double>(diagnostics.max_abs_difference);
    state.counters["max_population_fraction"] = static_cast<double>(diagnostics.max_population_fraction);
}

template <class Model>
void benchmark_flow_trajectory(benchmark::State& state, const Model& model,
                               const DirectComparisonDiagnostics& diagnostics)
{
    const auto output_points = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;
    FP checksum              = 0.0;

    for (auto _ : state) {
        state.PauseTiming();
        {
            auto simulation = make_flow_simulation(model, output_points);
            state.ResumeTiming();
            simulation.advance(tmax);
            state.PauseTiming();
            checksum += simulation.get_result().get_last_value().sum();
            checksum += simulation.get_flows().get_last_value().sum();
            benchmark::DoNotOptimize(checksum);
            benchmark::ClobberMemory();
        }
        state.ResumeTiming();
    }

    state.counters["N_C"]                     = static_cast<double>(model.get_initial_values().size());
    state.counters["N_T"]                     = static_cast<double>(model.get_initial_flows().size());
    state.counters["integrated_dim"]          = static_cast<double>(model.get_initial_flows().size());
    state.counters["max_abs_difference"]      = static_cast<double>(diagnostics.max_abs_difference);
    state.counters["max_population_fraction"] = static_cast<double>(diagnostics.max_population_fraction);
    state.counters["tracked_transitions"]     = static_cast<double>(model.get_initial_flows().size());
    state.counters["time_points"]             = static_cast<double>(output_points);
    state.SetItemsProcessed(state.iterations());
}

template <class Model>
void benchmark_flow_simulation(benchmark::State& state, int num_age_groups,
                               const DirectComparisonDiagnostics& diagnostics)
{
    const auto model = make_runtime_model<Model>(num_age_groups);
    benchmark_flow_trajectory(state, model, diagnostics);
}

template <class Model>
void benchmark_reused_auxiliary_simulation(benchmark::State& state, int num_age_groups, size_t tracked_transitions)
{
    const auto base_model        = make_runtime_model<Model>(num_age_groups);
    const size_t compartment_dim = static_cast<size_t>(base_model.get_initial_values().size());
    const size_t flow_dim        = static_cast<size_t>(base_model.get_initial_flows().size());
    if (tracked_transitions > flow_dim) {
        throw std::runtime_error("Requested more auxiliary counters than model transitions.");
    }
    const AuxiliaryTrackingModel<Model> model(base_model, tracked_transitions);
    benchmark_compartment_trajectory(state, model, compartment_dim, flow_dim, tracked_transitions);
}

void benchmark_recomputed_auxiliary_secir_simulation(benchmark::State& state, int num_age_groups,
                                                     const std::string& selection_order, size_t tracked_transitions)
{
    auto base_model              = make_runtime_model<mio::osecir::Model<FP>>(num_age_groups);
    const size_t compartment_dim = static_cast<size_t>(base_model.get_initial_values().size());
    const size_t flow_dim        = static_cast<size_t>(base_model.get_initial_flows().size());
    const size_t tracked_transition_pool =
        make_secir_transition_order(static_cast<size_t>(num_age_groups), selection_order).size();
    auto selected_transitions =
        select_secir_transitions(static_cast<size_t>(num_age_groups), selection_order, tracked_transitions);
    const RecomputedAuxiliarySecirModel model(std::move(base_model), std::move(selected_transitions));
    benchmark_compartment_trajectory(state, model, compartment_dim, flow_dim, tracked_transitions);
    state.counters["tracked_transition_pool"] = static_cast<double>(tracked_transition_pool);
}

template <class Model>
void register_direct_benchmarks(const std::string& model_name)
{
    for (int num_age_groups = 1; num_age_groups <= 6; ++num_age_groups) {
        const auto model       = make_runtime_model<Model>(num_age_groups);
        const auto diagnostics = validate_direct_comparison(model);
        const auto prefix      = "flow_paper/" + model_name + "/" + std::to_string(num_age_groups);
        benchmark::RegisterBenchmark((prefix + "/compartment").c_str(), [num_age_groups,
                                                                         diagnostics](benchmark::State& state) {
            benchmark_compartment_simulation<Model>(state, num_age_groups, diagnostics);
        })->Unit(benchmark::kMillisecond);
        benchmark::RegisterBenchmark((prefix + "/flow").c_str(), [num_age_groups,
                                                                  diagnostics](benchmark::State& state) {
            benchmark_flow_simulation<Model>(state, num_age_groups, diagnostics);
        })->Unit(benchmark::kMillisecond);
    }
}

void register_seir_rate_reuse_benchmarks()
{
    for (int num_age_groups = 1; num_age_groups <= 6; ++num_age_groups) {
        const auto flow_model = make_runtime_model<mio::oseir::Model<FP>>(num_age_groups);
        const ReusedRateSeirModel reused_model(flow_model.populations, flow_model.parameters);
        const EquationWiseRateSeirModel equation_wise_model(flow_model.populations, flow_model.parameters);
        const auto diagnostics = validate_seir_rate_reuse_comparison(flow_model, reused_model, equation_wise_model);
        const size_t compartment_dim = static_cast<size_t>(flow_model.get_initial_values().size());
        const size_t flow_dim        = static_cast<size_t>(flow_model.get_initial_flows().size());
        const auto prefix            = "flow_paper_rate_reuse/SEIR/" + std::to_string(num_age_groups);

        benchmark::RegisterBenchmark((prefix + "/compartment_cached").c_str(), [reused_model, compartment_dim, flow_dim,
                                                                                diagnostics](benchmark::State& state) {
            benchmark_compartment_trajectory(state, reused_model, compartment_dim, flow_dim, 0);
            state.counters["max_abs_difference"]      = static_cast<double>(diagnostics.max_abs_difference);
            state.counters["max_population_fraction"] = static_cast<double>(diagnostics.max_population_fraction);
        })->Unit(benchmark::kMillisecond);
        benchmark::RegisterBenchmark((prefix + "/compartment_equation_wise").c_str(), [equation_wise_model,
                                                                                       compartment_dim, flow_dim,
                                                                                       diagnostics](
                                                                                          benchmark::State& state) {
            benchmark_compartment_trajectory(state, equation_wise_model, compartment_dim, flow_dim, 0);
            state.counters["max_abs_difference"]      = static_cast<double>(diagnostics.max_abs_difference);
            state.counters["max_population_fraction"] = static_cast<double>(diagnostics.max_population_fraction);
        })->Unit(benchmark::kMillisecond);
        benchmark::RegisterBenchmark((prefix + "/flow").c_str(), [num_age_groups,
                                                                  diagnostics](benchmark::State& state) {
            benchmark_flow_simulation<mio::oseir::Model<FP>>(state, num_age_groups, diagnostics);
        })->Unit(benchmark::kMillisecond);
    }
}

template <SecirParameterization Parameterization>
void register_secir_implementation_parameterization(const std::string& parameterization_name, int num_age_groups,
                                                    const BenchmarkSecirFlowModel<Parameterization>& flow_model,
                                                    const EquationWiseSecirModel<Parameterization>& equation_wise_model,
                                                    const DirectComparisonDiagnostics& diagnostics)
{
    const size_t compartment_dim = static_cast<size_t>(flow_model.get_initial_values().size());
    const size_t flow_dim        = static_cast<size_t>(flow_model.get_initial_flows().size());
    const auto prefix =
        "flow_paper_implementation/SECIR/" + std::to_string(num_age_groups) + "/" + parameterization_name;

    benchmark::RegisterBenchmark((prefix + "/compartment_cached").c_str(), [flow_model, compartment_dim, flow_dim,
                                                                            diagnostics](benchmark::State& state) {
        benchmark_compartment_trajectory(state, flow_model, compartment_dim, flow_dim, 0);
        state.counters["max_abs_difference"]      = static_cast<double>(diagnostics.max_abs_difference);
        state.counters["max_population_fraction"] = static_cast<double>(diagnostics.max_population_fraction);
    })->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark((prefix + "/compartment_equation_wise").c_str(), [equation_wise_model, compartment_dim,
                                                                                   flow_dim, diagnostics](
                                                                                      benchmark::State& state) {
        benchmark_compartment_trajectory(state, equation_wise_model, compartment_dim, flow_dim, 0);
        state.counters["max_abs_difference"]      = static_cast<double>(diagnostics.max_abs_difference);
        state.counters["max_population_fraction"] = static_cast<double>(diagnostics.max_population_fraction);
    })->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark((prefix + "/flow").c_str(), [flow_model, diagnostics](benchmark::State& state) {
        benchmark_flow_trajectory(state, flow_model, diagnostics);
    })->Unit(benchmark::kMillisecond);
}

void register_secir_implementation_benchmarks()
{
    for (int num_age_groups = 1; num_age_groups <= 6; ++num_age_groups) {
        const auto source_model = make_runtime_model<mio::osecir::Model<FP>>(num_age_groups);
        const BenchmarkSecirFlowModel<SecirParameterization::Times> time_model(source_model.populations,
                                                                               source_model.parameters);
        const BenchmarkSecirFlowModel<SecirParameterization::Rates> rate_model(source_model.populations,
                                                                               source_model.parameters);
        const EquationWiseSecirModel<SecirParameterization::Times> time_equation_wise_model(time_model);
        const EquationWiseSecirModel<SecirParameterization::Rates> rate_equation_wise_model(rate_model);
        const auto diagnostics = validate_secir_implementation_comparison(
            source_model, time_model, time_equation_wise_model, rate_model, rate_equation_wise_model);

        register_secir_implementation_parameterization("times", num_age_groups, time_model, time_equation_wise_model,
                                                       diagnostics);
        register_secir_implementation_parameterization("rates", num_age_groups, rate_model, rate_equation_wise_model,
                                                       diagnostics);
    }
}

template <class Model>
void register_reused_auxiliary_benchmarks(const std::string& model_name)
{
    for (const int num_age_groups : {1, 3, 6}) {
        const auto model = make_runtime_model<Model>(num_age_groups);
        validate_reused_auxiliary_comparison(model);
        const size_t flow_dim = static_cast<size_t>(model.get_initial_flows().size());
        const auto prefix     = "flow_paper/" + model_name + "/" + std::to_string(num_age_groups) + "/auxiliary_reuse";
        for (size_t tracked_transitions = 0; tracked_transitions <= flow_dim; ++tracked_transitions) {
            benchmark::RegisterBenchmark((prefix + "/" + std::to_string(tracked_transitions)).c_str(),
                                         [num_age_groups, tracked_transitions](benchmark::State& state) {
                                             benchmark_reused_auxiliary_simulation<Model>(state, num_age_groups,
                                                                                          tracked_transitions);
                                         })
                ->Unit(benchmark::kMillisecond);
        }
    }
}

template <class Model>
void register_primary_state_diagnostics(const std::string& model_name)
{
    for (const int num_age_groups : {1, 3, 6}) {
        const auto base_model = make_runtime_model<Model>(num_age_groups);
        validate_primary_flow_comparison(base_model);
        const auto diagnostics       = validate_direct_comparison(base_model);
        const size_t compartment_dim = static_cast<size_t>(base_model.get_initial_values().size());
        const size_t flow_dim        = static_cast<size_t>(base_model.get_initial_flows().size());
        const PrimaryFlowTrackingModel<Model> flow_model(base_model);
        const AuxiliaryTrackingModel<Model> auxiliary_model(base_model, flow_dim);
        const auto prefix = "flow_paper_core/" + model_name + "/" + std::to_string(num_age_groups);

        benchmark::RegisterBenchmark((prefix + "/flow_primary").c_str(), [flow_model, compartment_dim,
                                                                          flow_dim](benchmark::State& state) {
            benchmark_compartment_trajectory(state, flow_model, compartment_dim, flow_dim, flow_dim);
        })->Unit(benchmark::kMillisecond);
        benchmark::RegisterBenchmark((prefix + "/auxiliary_full").c_str(), [auxiliary_model, compartment_dim,
                                                                            flow_dim](benchmark::State& state) {
            benchmark_compartment_trajectory(state, auxiliary_model, compartment_dim, flow_dim, flow_dim);
        })->Unit(benchmark::kMillisecond);
        benchmark::RegisterBenchmark((prefix + "/flow_complete").c_str(), [num_age_groups,
                                                                           diagnostics](benchmark::State& state) {
            benchmark_flow_simulation<Model>(state, num_age_groups, diagnostics);
        })->Unit(benchmark::kMillisecond);
        benchmark::RegisterBenchmark((prefix + "/flow_final").c_str(), [flow_model, compartment_dim,
                                                                        flow_dim](benchmark::State& state) {
            benchmark_final_state_only(state, flow_model, compartment_dim, flow_dim);
        })->Unit(benchmark::kMillisecond);
        benchmark::RegisterBenchmark((prefix + "/auxiliary_final").c_str(), [auxiliary_model, compartment_dim,
                                                                             flow_dim](benchmark::State& state) {
            benchmark_final_state_only(state, auxiliary_model, compartment_dim, flow_dim);
        })->Unit(benchmark::kMillisecond);
    }
}

void register_recomputed_auxiliary_secir_benchmarks()
{
    const std::array<std::string, 3> selection_orders{"group_wise", "other_first", "infection_first"};
    for (const int num_age_groups : {1, 3, 6}) {
        const auto model = make_runtime_model<mio::osecir::Model<FP>>(num_age_groups);
        validate_recomputed_secir_rates(model);
        const size_t tracked_transition_pool =
            make_secir_transition_order(static_cast<size_t>(num_age_groups), "group_wise").size();
        const auto full_reference =
            select_secir_transitions(static_cast<size_t>(num_age_groups), "group_wise", tracked_transition_pool);
        for (const auto& selection_order : selection_orders) {
            const auto full_selection =
                select_secir_transitions(static_cast<size_t>(num_age_groups), selection_order, tracked_transition_pool);
            const bool identical_full_selection =
                std::equal(full_reference.begin(), full_reference.end(), full_selection.begin(),
                           [](const auto& left, const auto& right) {
                               return static_cast<size_t>(left.group) == static_cast<size_t>(right.group) &&
                                      left.kind == right.kind;
                           });
            if (!identical_full_selection) {
                throw std::runtime_error("Full SECIR auxiliary selections differ between preference orders.");
            }
            const auto prefix =
                "flow_paper/SECIR/" + std::to_string(num_age_groups) + "/auxiliary_recomputed/" + selection_order;
            for (size_t tracked_transitions = 0; tracked_transitions <= tracked_transition_pool;
                 ++tracked_transitions) {
                benchmark::RegisterBenchmark(
                    (prefix + "/" + std::to_string(tracked_transitions)).c_str(),
                    [num_age_groups, selection_order, tracked_transitions](benchmark::State& state) {
                        benchmark_recomputed_auxiliary_secir_simulation(state, num_age_groups, selection_order,
                                                                        tracked_transitions);
                    })
                    ->Unit(benchmark::kMillisecond);
            }
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::critical);
    register_direct_benchmarks<mio::oseir::Model<FP>>("SEIR");
    register_direct_benchmarks<mio::osecir::Model<FP>>("SECIR");
    register_direct_benchmarks<mio::osecirts::Model<FP>>("SECIRTS");
    register_seir_rate_reuse_benchmarks();
    register_secir_implementation_benchmarks();
    register_reused_auxiliary_benchmarks<mio::oseir::Model<FP>>("SEIR");
    register_reused_auxiliary_benchmarks<mio::osecir::Model<FP>>("SECIR");
    register_reused_auxiliary_benchmarks<mio::osecirts::Model<FP>>("SECIRTS");
    register_recomputed_auxiliary_secir_benchmarks();
    if (std::getenv("FLOW_PAPER_CORE_DIAGNOSTIC") != nullptr) {
        register_primary_state_diagnostics<mio::oseir::Model<FP>>("SEIR");
        register_primary_state_diagnostics<mio::osecir::Model<FP>>("SECIR");
        register_primary_state_diagnostics<mio::osecirts::Model<FP>>("SECIRTS");
    }

    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }
    benchmark::AddCustomContext("integrator", "classical RK4");
    benchmark::AddCustomContext("t0", std::to_string(t0));
    benchmark::AddCustomContext("tmax", std::to_string(tmax));
    benchmark::AddCustomContext("dt", std::to_string(dt));
    benchmark::AddCustomContext("timed_region", "integration with benchmark-specific output retention");
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
