/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
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

#include "memilio/compartments/compartmental_model.h"
#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/math/euler.h"
#include "memilio/math/interpolation.h"
#include "memilio/utils/logging.h"
#include "ode_secir/model_reduced.h"
#include "ode_secirts/model.h"
#include "ode_seir/model.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/model.h"
#include "ode_sir/parameters.h"

#include "boost/numeric/odeint/stepper/runge_kutta4.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

using FP = ScalarType;

using PaperSecirState      = mio::osecir::reduced::InfectionState;
using PaperSecirModel      = mio::osecir::reduced::Model<FP>;
using SecirTransitionKind  = mio::osecir::reduced::Transition;

using SirState  = mio::osir::InfectionState;
using SirPop    = mio::Populations<FP, mio::AgeGroup, SirState>;
using SirParams = mio::osir::Parameters<FP>;
using SirFlows  = mio::TypeList<mio::Flow<SirState::Susceptible, SirState::Infected>,
                               mio::Flow<SirState::Infected, SirState::Recovered>>;

enum class AuxiliarySirState
{
    Susceptible,
    Infected,
    Recovered,
    CumulativeInfections,
    CumulativeRecoveries,
    Count
};

using AuxiliarySirPop = mio::Populations<FP, mio::AgeGroup, AuxiliarySirState>;

struct FlowValues {
    FP infections = 0.0;
    FP recoveries = 0.0;
};

struct SecirSelectedFlowValues {
    FP infections          = 0.0;
    FP hospital_admissions = 0.0;
    FP deaths              = 0.0;
};

struct CalendarDate {
    int year  = 2020;
    int month = 1;
    int day   = 1;
};

constexpr int policy_reduction_day = 30;

// All runtime experiments use the same trajectory and accumulate several trajectories per timing block.
constexpr FP runtime_t0                         = 0.0;
constexpr FP runtime_tmax                       = 600.0;
constexpr FP runtime_dt                         = 0.2;
constexpr size_t runtime_repetitions            = 10;
constexpr size_t runtime_trajectories_per_block = 32;

std::unique_ptr<mio::OdeIntegratorCore<FP>> make_euler_integrator()
{
    return std::make_unique<mio::EulerIntegratorCore<FP>>();
}

std::unique_ptr<mio::OdeIntegratorCore<FP>> make_fixed_rk4_integrator()
{
    return std::make_unique<mio::ExplicitStepperWrapper<FP, boost::numeric::odeint::runge_kutta4>>();
}

template <class Model>
auto make_compartment_runtime_simulation(FP t0, FP dt, const Model& model)
{
    mio::Simulation<FP, Model> simulation(model, t0, dt);
    simulation.set_integrator_core(make_fixed_rk4_integrator());
    return simulation;
}

template <class Model>
auto make_flow_runtime_simulation(FP t0, FP dt, const Model& model)
{
    mio::FlowSimulation<FP, Model> simulation(model, t0, dt);
    simulation.set_integrator_core(make_fixed_rk4_integrator());
    return simulation;
}

size_t get_environment_size(const char* name, size_t fallback)
{
    const char* value = std::getenv(name);
    if (value == nullptr) {
        return fallback;
    }
    return static_cast<size_t>(std::stoull(value));
}

std::vector<int> get_environment_ints(const char* name, std::vector<int> fallback)
{
    const char* value = std::getenv(name);
    if (value == nullptr) {
        return fallback;
    }
    std::string text(value);
    std::replace(text.begin(), text.end(), ',', ' ');
    std::istringstream stream(text);
    std::vector<int> values;
    for (int parsed; stream >> parsed;) {
        values.push_back(parsed);
    }
    return values.empty() ? fallback : values;
}

bool environment_contains(const char* name, const std::string& candidate)
{
    const char* value = std::getenv(name);
    if (value == nullptr) {
        return true;
    }
    std::string text(value);
    std::replace(text.begin(), text.end(), ',', ' ');
    std::istringstream stream(text);
    for (std::string parsed; stream >> parsed;) {
        if (parsed == candidate) {
            return true;
        }
    }
    return false;
}

void write_header(std::ofstream& file, const std::vector<std::string>& columns)
{
    for (size_t i = 0; i < columns.size(); ++i) {
        if (i > 0) {
            file << ",";
        }
        file << columns[i];
    }
    file << "\n";
}

std::ofstream open_csv(const std::filesystem::path& path)
{
    std::ofstream file(path);
    file << std::setprecision(14);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open " + path.string());
    }
    return file;
}

bool is_leap_year(int year)
{
    return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}

int days_in_month(int year, int month)
{
    static const std::array<int, 12> days{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (month == 2 && is_leap_year(year)) {
        return 29;
    }
    return days[static_cast<size_t>(month - 1)];
}

CalendarDate parse_date(const std::string& text)
{
    if (text.size() < 10) {
        throw std::runtime_error("Malformed date: " + text);
    }
    return {std::stoi(text.substr(0, 4)), std::stoi(text.substr(5, 2)), std::stoi(text.substr(8, 2))};
}

CalendarDate add_days(CalendarDate date, int offset)
{
    date.day += offset;
    while (date.day > days_in_month(date.year, date.month)) {
        date.day -= days_in_month(date.year, date.month);
        ++date.month;
        if (date.month > 12) {
            date.month = 1;
            ++date.year;
        }
    }
    return date;
}

int day_of_year(CalendarDate date)
{
    int doy = date.day;
    for (int month = 1; month < date.month; ++month) {
        doy += days_in_month(date.year, month);
    }
    return doy;
}

std::string format_date(CalendarDate date)
{
    std::ostringstream out;
    out << date.year << "-" << std::setw(2) << std::setfill('0') << date.month << "-" << std::setw(2) << date.day
        << std::setfill(' ');
    return out.str();
}

std::string date_plus(CalendarDate start_date, int offset)
{
    return format_date(add_days(start_date, offset));
}

std::string npi_phase_for_day(int day)
{
    if (day < policy_reduction_day) {
        return "baseline_contacts";
    }
    if (day < policy_reduction_day + 1) {
        return "contact_reduction_start";
    }
    return "contact_reduction";
}

FP contact_multiplier_for_day(int day)
{
    if (day < policy_reduction_day + 1) {
        return 1.0;
    }
    return 0.4;
}

bool assign_secir_population_value(PaperSecirModel& model, mio::AgeGroup group, const std::string& state_name,
                                   FP value)
{
    using IS = PaperSecirState;
    if (state_name == "Susceptible") {
        model.populations[{group, IS::Susceptible}] += value;
    }
    else if (state_name == "Exposed") {
        model.populations[{group, IS::Exposed}] += value;
    }
    else if (state_name == "InfectedNoSymptoms" || state_name == "InfectedNoSymptomsConfirmed") {
        model.populations[{group, IS::InfectedNoSymptoms}] += value;
    }
    else if (state_name == "InfectedSymptoms" || state_name == "InfectedSymptomsConfirmed") {
        model.populations[{group, IS::InfectedSymptoms}] += value;
    }
    else if (state_name == "InfectedSevere") {
        model.populations[{group, IS::InfectedSevere}] += value;
    }
    else if (state_name == "InfectedCritical") {
        model.populations[{group, IS::InfectedCritical}] += value;
    }
    else if (state_name == "Recovered") {
        model.populations[{group, IS::Recovered}] += value;
    }
    else if (state_name == "Dead") {
        model.populations[{group, IS::Dead}] += value;
    }
    else {
        return false;
    }
    return true;
}

bool read_germany_initial_conditions(PaperSecirModel& model, const std::filesystem::path& path)
{
    std::ifstream file(path);
    if (!file.is_open()) {
        return false;
    }

    model.populations.array().setConstant(0.0);
    std::string line;
    std::getline(file, line);

    size_t num_rows = 0;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        std::stringstream row(line);
        std::string age_text;
        std::string state_name;
        std::string value_text;
        std::getline(row, age_text, ',');
        std::getline(row, state_name, ',');
        std::getline(row, value_text, ',');
        if (age_text.empty() || state_name.empty() || value_text.empty()) {
            throw std::runtime_error("Malformed Germany initial-condition row: " + line);
        }
        if (!value_text.empty() && value_text.back() == '\r') {
            value_text.pop_back();
        }

        const auto group = mio::AgeGroup(static_cast<size_t>(std::stoul(age_text)));
        const FP value   = static_cast<FP>(std::stod(value_text));
        if (!assign_secir_population_value(model, group, state_name, value)) {
            throw std::runtime_error("Unknown SECIR state in Germany initial-condition row: " + state_name);
        }
        ++num_rows;
    }

    return num_rows > 0;
}

CalendarDate read_germany_start_date(const std::filesystem::path& path)
{
    std::ifstream file(path);
    if (!file.is_open()) {
        return {2020, 10, 15};
    }

    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::stringstream row(line);
        std::string key;
        std::string value;
        std::getline(row, key, ',');
        std::getline(row, value, ',');
        if (key == "start_date") {
            if (!value.empty() && value.back() == '\r') {
                value.pop_back();
            }
            return parse_date(value);
        }
    }
    return {2020, 10, 15};
}

FP total_population(Eigen::Ref<const Eigen::VectorX<FP>> y, size_t s, size_t i, size_t r)
{
    return y[s] + y[i] + y[r];
}

FP infection_coefficient(const SirParams& params, FP t, mio::AgeGroup group)
{
    const auto& contact_matrix = params.template get<mio::osir::ContactPatterns<FP>>().get_cont_freq_mat();
    const auto contact_rate    = contact_matrix.get_matrix_at(mio::SimulationTime<FP>(t))(
        static_cast<Eigen::Index>((size_t)group), static_cast<Eigen::Index>((size_t)group));
    return contact_rate * params.template get<mio::osir::TransmissionProbabilityOnContact<FP>>()[group];
}

FlowValues compute_sir_rates(const SirParams& params, Eigen::Ref<const Eigen::VectorX<FP>> pop,
                             Eigen::Ref<const Eigen::VectorX<FP>> y, FP t, size_t s, size_t i, size_t r)
{
    const FP n = total_population(pop, s, i, r);
    if (n < mio::Limits<FP>::zero_tolerance()) {
        return {};
    }

    const auto group        = mio::AgeGroup(0);
    const FP infection_rate = infection_coefficient(params, t, group) * y[s] * pop[i] / n;
    const FP recovery_rate  = y[i] / params.template get<mio::osir::TimeInfected<FP>>()[group];
    return {infection_rate, recovery_rate};
}

class FlowSirModel : public mio::FlowModel<FP, SirState, SirPop, SirParams, SirFlows>
{
public:
    using Base = mio::FlowModel<FP, SirState, SirPop, SirParams, SirFlows>;

    FlowSirModel(const SirPop& initial_populations, const SirParams& model_parameters)
        : Base(initial_populations, model_parameters)
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const auto group = mio::AgeGroup(0);
        const size_t s   = this->populations.get_flat_index({group, SirState::Susceptible});
        const size_t i   = this->populations.get_flat_index({group, SirState::Infected});
        const size_t r   = this->populations.get_flat_index({group, SirState::Recovered});

        const auto rates = compute_sir_rates(this->parameters, pop, y, t, s, i, r);
        flows[this->template get_flat_flow_index<SirState::Susceptible, SirState::Infected>(group)] = rates.infections;
        flows[this->template get_flat_flow_index<SirState::Infected, SirState::Recovered>(group)]   = rates.recoveries;
    }
};

class AuxiliarySirModel : public mio::CompartmentalModel<FP, AuxiliarySirState, AuxiliarySirPop, SirParams>
{
public:
    using Base = mio::CompartmentalModel<FP, AuxiliarySirState, AuxiliarySirPop, SirParams>;

    AuxiliarySirModel(const AuxiliarySirPop& initial_populations, const SirParams& model_parameters)
        : Base(initial_populations, model_parameters)
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                         Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
    {
        const auto group = mio::AgeGroup(0);
        const size_t s   = this->populations.get_flat_index({group, AuxiliarySirState::Susceptible});
        const size_t i   = this->populations.get_flat_index({group, AuxiliarySirState::Infected});
        const size_t r   = this->populations.get_flat_index({group, AuxiliarySirState::Recovered});
        const size_t cumulative_infections =
            this->populations.get_flat_index({group, AuxiliarySirState::CumulativeInfections});
        const size_t cumulative_recoveries =
            this->populations.get_flat_index({group, AuxiliarySirState::CumulativeRecoveries});

        const auto rates = compute_sir_rates(this->parameters, pop, y, t, s, i, r);
        dydt[s] -= rates.infections;
        dydt[i] += rates.infections - rates.recoveries;
        dydt[r] += rates.recoveries;
        dydt[cumulative_infections] += rates.infections;
        dydt[cumulative_recoveries] += rates.recoveries;
    }
};

SirParams make_sir_parameters()
{
    SirParams params(mio::AgeGroup(1));
    params.set<mio::osir::TimeInfected<FP>>(4.0);
    params.set<mio::osir::TransmissionProbabilityOnContact<FP>>(0.08);
    params.get<mio::osir::ContactPatterns<FP>>().get_cont_freq_mat()[0].get_baseline().setConstant(8.0);
    return params;
}

SirPop make_sir_populations()
{
    SirPop populations({mio::AgeGroup(1), SirState::Count});
    const auto group                            = mio::AgeGroup(0);
    populations[{group, SirState::Susceptible}] = 9990.0;
    populations[{group, SirState::Infected}]    = 10.0;
    populations[{group, SirState::Recovered}]   = 0.0;
    return populations;
}

AuxiliarySirPop make_auxiliary_sir_populations()
{
    AuxiliarySirPop populations({mio::AgeGroup(1), AuxiliarySirState::Count}, 0.0);
    const auto group                                     = mio::AgeGroup(0);
    populations[{group, AuxiliarySirState::Susceptible}] = 9990.0;
    populations[{group, AuxiliarySirState::Infected}]    = 10.0;
    populations[{group, AuxiliarySirState::Recovered}]   = 0.0;
    return populations;
}

mio::osir::Model<FP> make_classical_sir_model()
{
    mio::osir::Model<FP> model(make_sir_populations(), make_sir_parameters());
    model.check_constraints();
    return model;
}

FlowSirModel make_flow_sir_model()
{
    FlowSirModel model(make_sir_populations(), make_sir_parameters());
    model.check_constraints();
    return model;
}

AuxiliarySirModel make_auxiliary_sir_model()
{
    AuxiliarySirModel model(make_auxiliary_sir_populations(), make_sir_parameters());
    model.check_constraints();
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
        const auto comp_dim = static_cast<Eigen::Index>(m_compartment_dim);
        auto dydt_comp      = dydt.head(comp_dim);

        m_flow_buffer.setZero();
        m_model.get_flows(pop.head(comp_dim), y.head(comp_dim), t, m_flow_buffer);
        m_model.get_derivatives(m_flow_buffer, dydt_comp);

        if (m_num_auxiliary_states == 0) {
            return;
        }

        for (size_t i = 0; i < m_num_auxiliary_states; ++i) {
            dydt[static_cast<Eigen::Index>(m_compartment_dim + i)] =
                m_flow_buffer[static_cast<Eigen::Index>(i)];
        }
    }

private:
    Model m_model;
    size_t m_num_auxiliary_states;
    size_t m_compartment_dim;
    Eigen::VectorX<FP> m_initial;
    mutable Eigen::VectorX<FP> m_flow_buffer;
};

struct SecirTrackedTransition {
    mio::AgeGroup group;
    SecirTransitionKind kind;
};

inline constexpr auto secir_transition_kinds = mio::osecir::reduced::transitions;

std::vector<SecirTrackedTransition> make_secir_transition_order(size_t num_age_groups,
                                                                const std::string& selection_order)
{
    std::vector<SecirTrackedTransition> transitions;
    transitions.reserve(secir_transition_kinds.size() * num_age_groups);

    const auto append_infections = [&]() {
        for (size_t group = 0; group < num_age_groups; ++group) {
            transitions.push_back({mio::AgeGroup(group), SecirTransitionKind::Infection});
        }
    };
    const auto append_remaining_group_wise = [&]() {
        for (size_t group = 0; group < num_age_groups; ++group) {
            for (const auto kind : secir_transition_kinds) {
                if (kind != SecirTransitionKind::Infection) {
                    transitions.push_back({mio::AgeGroup(group), kind});
                }
            }
        }
    };

    if (selection_order == "group_wise") {
        for (size_t group = 0; group < num_age_groups; ++group) {
            for (const auto kind : secir_transition_kinds) {
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

FP recompute_secir_transition_rate(const PaperSecirModel& model, const SecirTrackedTransition& transition,
                                   Eigen::Ref<const Eigen::VectorX<FP>> pop,
                                   Eigen::Ref<const Eigen::VectorX<FP>> y, FP t)
{
    return model.get_transition_rate(transition.kind, transition.group, pop, y, t);
}

size_t get_secir_transition_index(const PaperSecirModel& model, const SecirTrackedTransition& transition)
{
    return model.get_transition_index(transition.kind, transition.group);
}

void validate_recomputed_secir_rates(const PaperSecirModel& model)
{
    const auto state = model.get_initial_values();
    const auto transitions =
        make_secir_transition_order(static_cast<size_t>(model.parameters.get_num_groups()), "group_wise");
    for (const FP time : {FP(0.0), FP(30.0), FP(90.0)}) {
        auto flow_rates = model.get_initial_flows();
        model.get_flows(state, state, time, flow_rates);
        for (const auto& transition : transitions) {
            const FP expected = flow_rates[static_cast<Eigen::Index>(get_secir_transition_index(model, transition))];
            const FP recomputed = recompute_secir_transition_rate(model, transition, state, state, time);
            const FP tolerance  = 1e-12 * std::max(FP(1.0), std::abs(expected));
            if (std::abs(expected - recomputed) > tolerance) {
                throw std::runtime_error("Naive SECIR transition recomputation does not match the model rate.");
            }
        }
    }
}

class NaiveAuxiliarySecirModel
{
public:
    NaiveAuxiliarySecirModel(PaperSecirModel model, std::vector<SecirTrackedTransition> tracked_transitions)
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
    PaperSecirModel m_model;
    std::vector<SecirTrackedTransition> m_tracked_transitions;
    size_t m_compartment_dim;
    Eigen::VectorX<FP> m_initial;
    mutable Eigen::VectorX<FP> m_flow_buffer;
};

FlowValues last_flow_values(const mio::TimeSeries<FP>& flows)
{
    const auto& last = flows.get_last_value();
    return {last[0], last[1]};
}

FlowValues auxiliary_flow_values(const mio::TimeSeries<FP>& result)
{
    const auto& last = result.get_last_value();
    return {last[3], last[4]};
}

FlowValues sir_balance_flows(const mio::TimeSeries<FP>& result, Eigen::Index idx)
{
    const auto& initial = result.get_value(0);
    const auto& value   = result.get_value(idx);
    return {initial[0] - value[0], value[2] - initial[2]};
}

Eigen::Vector2<FP> compute_rates_from_classical_sir_state(const mio::osir::Model<FP>& model,
                                                          Eigen::Ref<const Eigen::VectorX<FP>> y, FP t)
{
    const auto group = mio::AgeGroup(0);
    const size_t s   = model.populations.get_flat_index({group, SirState::Susceptible});
    const size_t i   = model.populations.get_flat_index({group, SirState::Infected});
    const size_t r   = model.populations.get_flat_index({group, SirState::Recovered});

    const auto rates = compute_sir_rates(model.parameters, y, y, t, s, i, r);
    return Eigen::Vector2<FP>(rates.infections, rates.recoveries);
}

mio::TimeSeries<FP> reconstruct_sir_flows_a_posteriori(const mio::osir::Model<FP>& model,
                                                       const mio::TimeSeries<FP>& result)
{
    mio::TimeSeries<FP> flows(result.get_time(0), Eigen::Vector2<FP>::Zero());
    for (Eigen::Index k = 1; k < result.get_num_time_points(); ++k) {
        const FP t_prev = result.get_time(k - 1);
        const FP t      = result.get_time(k);
        const FP dt     = t - t_prev;

        const auto rates_prev = compute_rates_from_classical_sir_state(model, result.get_value(k - 1), t_prev);
        const auto rates      = compute_rates_from_classical_sir_state(model, result.get_value(k), t);
        const auto next_value = (flows.get_last_value() + 0.5 * dt * (rates_prev + rates)).eval();
        flows.add_time_point(t, next_value);
    }
    return flows;
}

mio::TimeSeries<FP> subsample_time_series(const mio::TimeSeries<FP>& result, Eigen::Index stride)
{
    mio::TimeSeries<FP> subsampled(result.get_time(0), result.get_value(0));
    for (Eigen::Index k = stride; k < result.get_num_time_points(); k += stride) {
        subsampled.add_time_point(result.get_time(k), result.get_value(k));
    }
    if (subsampled.get_last_time() < result.get_last_time()) {
        subsampled.add_time_point(result.get_last_time(), result.get_last_value());
    }
    return subsampled;
}

void generate_sir_method_comparison(const std::filesystem::path& data_dir)
{
    const FP t0   = 0.0;
    const FP tmax = 20.0;
    const FP dt   = 0.1;

    auto classical_model = make_classical_sir_model();
    auto flow_model      = make_flow_sir_model();
    auto auxiliary_model = make_auxiliary_sir_model();

    const auto classical_result = mio::simulate<FP>(t0, tmax, dt, classical_model, make_euler_integrator());
    const auto flow_result      = mio::simulate_flows<FP>(t0, tmax, dt, flow_model, make_euler_integrator());
    const auto auxiliary_result = mio::simulate<FP>(t0, tmax, dt, auxiliary_model, make_euler_integrator());
    const auto aposteriori      = reconstruct_sir_flows_a_posteriori(classical_model, classical_result);

    auto summary = open_csv(data_dir / "sir_method_comparison.csv");
    write_header(summary, {"method", "cumulative_infections", "cumulative_recoveries"});
    const auto direct_flows      = last_flow_values(flow_result[1]);
    const auto auxiliary         = auxiliary_flow_values(auxiliary_result);
    const auto aposteriori_final = last_flow_values(aposteriori);
    const auto balances          = sir_balance_flows(classical_result, classical_result.get_num_time_points() - 1);
    summary << "Flow," << direct_flows.infections << "," << direct_flows.recoveries << "\n";
    summary << "Auxiliary," << auxiliary.infections << "," << auxiliary.recoveries << "\n";
    summary << "A-posteriori," << aposteriori_final.infections << "," << aposteriori_final.recoveries << "\n";
    summary << "Balance," << balances.infections << "," << balances.recoveries << "\n";

    auto series = open_csv(data_dir / "sir_method_timeseries.csv");
    write_header(series,
                 {"time", "flow_infections", "flow_recoveries", "auxiliary_infections", "auxiliary_recoveries",
                  "aposteriori_infections", "aposteriori_recoveries", "balance_infections", "balance_recoveries"});
    for (Eigen::Index k = 0; k < flow_result[1].get_num_time_points(); ++k) {
        const auto balance = sir_balance_flows(classical_result, k);
        series << flow_result[1].get_time(k) << "," << flow_result[1].get_value(k)[0] << ","
               << flow_result[1].get_value(k)[1] << "," << auxiliary_result.get_value(k)[3] << ","
               << auxiliary_result.get_value(k)[4] << "," << aposteriori.get_value(k)[0] << ","
               << aposteriori.get_value(k)[1] << "," << balance.infections << "," << balance.recoveries << "\n";
    }
}

void generate_grid_sensitivity(const std::filesystem::path& data_dir)
{
    const FP t0           = 0.0;
    const FP tmax         = 20.0;
    const FP solver_dt    = 0.01;
    const auto model      = make_classical_sir_model();
    const auto flow_model = make_flow_sir_model();

    const auto fine_compartments = mio::simulate<FP>(t0, tmax, solver_dt, model, make_euler_integrator());
    const auto fine_flows        = mio::simulate_flows<FP>(t0, tmax, solver_dt, flow_model, make_euler_integrator())[1];
    const auto reference         = last_flow_values(fine_flows);

    auto csv = open_csv(data_dir / "sir_grid_sensitivity.csv");
    write_header(csv, {"output_dt", "num_output_points", "infection_error_abs", "infection_error_rel",
                       "recovery_error_abs", "recovery_error_rel"});
    for (Eigen::Index stride : {1, 2, 5, 10, 25, 50, 100, 200, 500, 1000}) {
        const auto sampled         = subsample_time_series(fine_compartments, stride);
        const auto reconstructed   = reconstruct_sir_flows_a_posteriori(model, sampled);
        const auto estimated       = last_flow_values(reconstructed);
        const auto infection_error = std::abs(estimated.infections - reference.infections);
        const auto recovery_error  = std::abs(estimated.recoveries - reference.recoveries);
        csv << stride * solver_dt << "," << sampled.get_num_time_points() << "," << infection_error << ","
            << infection_error / reference.infections << "," << recovery_error << ","
            << recovery_error / reference.recoveries << "\n";
    }

    auto paths = open_csv(data_dir / "sir_grid_reconstruction_paths.csv");
    write_header(paths, {"output_dt", "time", "flow_infections", "flow_recoveries", "aposteriori_infections",
                         "aposteriori_recoveries"});
    for (Eigen::Index stride : {10, 100, 500}) {
        const auto sampled       = subsample_time_series(fine_compartments, stride);
        const auto reconstructed = reconstruct_sir_flows_a_posteriori(model, sampled);
        for (Eigen::Index k = 0; k < sampled.get_num_time_points(); ++k) {
            const auto fine_idx = static_cast<Eigen::Index>(std::llround(sampled.get_time(k) / solver_dt));
            paths << stride * solver_dt << "," << sampled.get_time(k) << "," << fine_flows.get_value(fine_idx)[0] << ","
                  << fine_flows.get_value(fine_idx)[1] << "," << reconstructed.get_value(k)[0] << ","
                  << reconstructed.get_value(k)[1] << "\n";
        }
    }
}

PaperSecirModel make_policy_secir_model(const std::filesystem::path& data_dir)
{
    using IS = PaperSecirState;

    const size_t num_groups = 6;
    PaperSecirModel model(static_cast<int>(num_groups));
    const auto start_date = read_germany_start_date(data_dir / "germany_secir_initial_summary.csv");

    const FP total_population = 1000000.0;
    const std::array<FP, 6> age_share{0.055, 0.095, 0.245, 0.335, 0.205, 0.065};
    const std::array<FP, 6> exposed_share{0.05, 0.08, 0.29, 0.34, 0.18, 0.06};
    const std::array<FP, 6> severe_prob{0.006, 0.008, 0.02, 0.055, 0.14, 0.24};
    const std::array<FP, 6> critical_prob{0.05, 0.05, 0.08, 0.13, 0.27, 0.38};
    const std::array<FP, 6> death_critical{0.02, 0.03, 0.08, 0.16, 0.34, 0.58};

    const FP exposed_total     = 2000.0;
    const FP carriers_total    = 1300.0;
    const FP symptomatic_total = 950.0;
    const FP severe_total      = 240.0;
    const FP critical_total    = 80.0;
    const FP recovered_total   = 70000.0;

    model.parameters.template set<mio::osecir::StartDay<FP>>(static_cast<FP>(day_of_year(start_date)));
    model.parameters.set<mio::osecir::Seasonality<FP>>(0.20);
    model.parameters.get<mio::osecir::TestAndTraceCapacity<FP>>() = 2500.0;

    mio::ContactMatrixGroup<FP>& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<FP>>();
    Eigen::MatrixX<FP> baseline(num_groups, num_groups);
    baseline.setConstant(1.45);
    for (Eigen::Index i = 0; i < baseline.rows(); ++i) {
        baseline(i, i) = 9.2;
    }
    contact_matrix[0] = mio::ContactMatrix<FP>(baseline);
    contact_matrix[0].add_damping(0.6, mio::SimulationTime<FP>(policy_reduction_day + 1.0));

    model.populations.set_total(total_population);
    for (size_t group_idx = 0; group_idx < num_groups; ++group_idx) {
        const auto group = mio::AgeGroup(group_idx);
        const FP scale   = exposed_share[group_idx];

        model.parameters.get<mio::osecir::TimeExposed<FP>>()[group]            = 3.2;
        model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[group] = 2.0;
        model.parameters.get<mio::osecir::TimeInfectedSymptoms<FP>>()[group]   = 6.0;
        model.parameters.get<mio::osecir::TimeInfectedSevere<FP>>()[group]     = 8.5 + 0.8 * group_idx;
        model.parameters.get<mio::osecir::TimeInfectedCritical<FP>>()[group]   = 7.0 + 0.5 * group_idx;

        model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[group]  = 0.024;
        model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[group]    = 0.7;
        model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[group]    = 0.22;
        model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<FP>>()[group]    = 0.18;
        model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>>()[group] = 0.38;
        model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[group]         = severe_prob[group_idx];
        model.parameters.get<mio::osecir::CriticalPerSevere<FP>>()[group]                 = critical_prob[group_idx];
        model.parameters.get<mio::osecir::DeathsPerSevere<FP>>()[group]                   = 0.0;
        model.parameters.get<mio::osecir::DeathsPerCritical<FP>>()[group]                 = death_critical[group_idx];

        model.populations[{group, IS::Exposed}]            = exposed_total * scale;
        model.populations[{group, IS::InfectedNoSymptoms}] = carriers_total * scale;
        model.populations[{group, IS::InfectedSymptoms}]   = symptomatic_total * scale;
        model.populations[{group, IS::InfectedSevere}] =
            severe_total * age_share[group_idx] * (0.4 + 4.0 * severe_prob[group_idx]);
        model.populations[{group, IS::InfectedCritical}] =
            critical_total * age_share[group_idx] * (0.3 + 3.5 * critical_prob[group_idx]);
        model.populations[{group, IS::Recovered}] = recovered_total * age_share[group_idx];
        model.populations[{group, IS::Dead}]      = 0.0;
        const FP assigned = model.populations[{group, IS::Exposed}] +
                            model.populations[{group, IS::InfectedNoSymptoms}] +
                            model.populations[{group, IS::InfectedSymptoms}] +
                            model.populations[{group, IS::InfectedSevere}] +
                            model.populations[{group, IS::InfectedCritical}] +
                            model.populations[{group, IS::Recovered}] + model.populations[{group, IS::Dead}];
        model.populations[{group, IS::Susceptible}] =
            total_population * age_share[group_idx] - assigned;
    }

    read_germany_initial_conditions(model, data_dir / "germany_secir_initial_conditions.csv");
    model.apply_constraints();
    return model;
}

template <PaperSecirState Source, PaperSecirState Target>
size_t secir_flow_index(const PaperSecirModel& model, mio::AgeGroup group)
{
    return model.template get_flat_flow_index<Source, Target>({group});
}

size_t secir_state_index(const PaperSecirModel& model, mio::AgeGroup group, PaperSecirState state)
{
    return model.populations.get_flat_index({group, state});
}

template <PaperSecirState Source, PaperSecirState Target>
FP aggregate_secir_flow_value(const PaperSecirModel& model, Eigen::Ref<const Eigen::VectorX<FP>> flow_values)
{
    FP sum = 0.0;
    for (mio::AgeGroup group = 0; group < model.parameters.get_num_groups(); ++group) {
        sum += flow_values[static_cast<Eigen::Index>(secir_flow_index<Source, Target>(model, group))];
    }
    return sum;
}

template <PaperSecirState Source, PaperSecirState Target>
FP aggregate_secir_flow_difference(const PaperSecirModel& model,
                                   Eigen::Ref<const Eigen::VectorX<FP>> flow_values,
                                   Eigen::Ref<const Eigen::VectorX<FP>> previous_flow_values)
{
    return aggregate_secir_flow_value<Source, Target>(model, flow_values) -
           aggregate_secir_flow_value<Source, Target>(model, previous_flow_values);
}

FP aggregate_secir_state_value(const PaperSecirModel& model, Eigen::Ref<const Eigen::VectorX<FP>> compartments,
                               PaperSecirState state)
{
    FP sum = 0.0;
    for (mio::AgeGroup group = 0; group < model.parameters.get_num_groups(); ++group) {
        sum += compartments[static_cast<Eigen::Index>(secir_state_index(model, group, state))];
    }
    return sum;
}

std::vector<mio::TimeSeries<FP>> simulate_policy_flows_adaptively(const PaperSecirModel& model, FP t0, FP tmax,
                                                                  FP dt)
{
    return mio::simulate_flows<FP>(t0, tmax, dt, model);
}

SecirSelectedFlowValues aggregate_selected_secir_flow_values(const PaperSecirModel& model,
                                                             Eigen::Ref<const Eigen::VectorX<FP>> flow_values)
{
    using IS = PaperSecirState;

    SecirSelectedFlowValues values;
    for (mio::AgeGroup group = 0; group < model.parameters.get_num_groups(); ++group) {
        values.infections +=
            flow_values[static_cast<Eigen::Index>(secir_flow_index<IS::Susceptible, IS::Exposed>(model, group))];
        values.hospital_admissions += flow_values[static_cast<Eigen::Index>(
            secir_flow_index<IS::InfectedSymptoms, IS::InfectedSevere>(model, group))];
        values.deaths += flow_values[static_cast<Eigen::Index>(
            secir_flow_index<IS::InfectedCritical, IS::Dead>(model, group))];
    }
    return values;
}

Eigen::Vector3<FP> selected_secir_flow_vector(const SecirSelectedFlowValues& values)
{
    return Eigen::Vector3<FP>(values.infections, values.hospital_admissions, values.deaths);
}

Eigen::Vector3<FP> compute_selected_secir_rates_from_state(const PaperSecirModel& model,
                                                           Eigen::Ref<const Eigen::VectorX<FP>> y, FP t)
{
    Eigen::VectorX<FP> flow_rates = Eigen::VectorX<FP>::Zero(model.get_initial_flows().size());
    model.get_flows(y, y, t, flow_rates);
    return selected_secir_flow_vector(aggregate_selected_secir_flow_values(model, flow_rates));
}

mio::TimeSeries<FP> reconstruct_secir_flows_a_posteriori(const PaperSecirModel& model,
                                                         const mio::TimeSeries<FP>& compartments)
{
    mio::TimeSeries<FP> reconstructed(compartments.get_time(0), Eigen::Vector3<FP>::Zero());
    for (Eigen::Index k = 1; k < compartments.get_num_time_points(); ++k) {
        const FP t_prev = compartments.get_time(k - 1);
        const FP t      = compartments.get_time(k);
        const FP dt     = t - t_prev;

        const auto rates_prev = compute_selected_secir_rates_from_state(model, compartments.get_value(k - 1), t_prev);
        reconstructed.add_time_point(t, (reconstructed.get_last_value() + dt * rates_prev).eval());
    }
    return reconstructed;
}

Eigen::Vector3<FP> aggregate_selected_secir_cumulative_flows(const PaperSecirModel& model,
                                                             const mio::TimeSeries<FP>& flows, Eigen::Index idx)
{
    return selected_secir_flow_vector(aggregate_selected_secir_flow_values(model, flows.get_value(idx)));
}

void generate_secir_grid_sensitivity(const std::filesystem::path& data_dir)
{
    const FP t0           = 0.0;
    const FP tmax         = 60.0;
    const FP solver_dt    = 0.05;
    const auto model      = make_policy_secir_model(data_dir);
    const auto flow_model = make_policy_secir_model(data_dir);

    const auto fine_compartments = mio::simulate<FP>(t0, tmax, solver_dt, model, make_euler_integrator());
    const auto fine_flows        = mio::simulate_flows<FP>(t0, tmax, solver_dt, flow_model, make_euler_integrator())[1];
    const auto reference =
        aggregate_selected_secir_cumulative_flows(flow_model, fine_flows, fine_flows.get_num_time_points() - 1);

    auto csv = open_csv(data_dir / "secir_grid_sensitivity.csv");
    write_header(csv, {"output_dt", "num_output_points", "infection_error_abs", "infection_error_rel",
                       "hospital_admission_error_abs", "hospital_admission_error_rel", "death_error_abs",
                       "death_error_rel"});
    for (Eigen::Index stride : {1, 2, 5, 10, 20, 40, 100, 200}) {
        const auto sampled       = subsample_time_series(fine_compartments, stride);
        const auto reconstructed = reconstruct_secir_flows_a_posteriori(model, sampled);
        const auto estimated     = reconstructed.get_last_value();
        const auto error         = (estimated - reference).cwiseAbs();

        csv << stride * solver_dt << "," << sampled.get_num_time_points() << "," << error[0] << ","
            << error[0] / reference[0] << "," << error[1] << "," << error[1] / reference[1] << "," << error[2] << ","
            << error[2] / reference[2] << "\n";
    }

    auto paths = open_csv(data_dir / "secir_grid_reconstruction_paths.csv");
    write_header(paths, {"output_dt", "time", "flow_infections", "flow_hospital_admissions", "flow_deaths",
                         "aposteriori_infections", "aposteriori_hospital_admissions", "aposteriori_deaths"});
    for (Eigen::Index stride : {5, 40, 200}) {
        const auto sampled       = subsample_time_series(fine_compartments, stride);
        const auto reconstructed = reconstruct_secir_flows_a_posteriori(model, sampled);
        for (Eigen::Index k = 0; k < sampled.get_num_time_points(); ++k) {
            const auto fine_idx = static_cast<Eigen::Index>(std::llround(sampled.get_time(k) / solver_dt));
            const auto direct   = aggregate_selected_secir_cumulative_flows(flow_model, fine_flows, fine_idx);
            const auto estimate = reconstructed.get_value(k);
            paths << stride * solver_dt << "," << sampled.get_time(k) << "," << direct[0] << "," << direct[1] << ","
                  << direct[2] << "," << estimate[0] << "," << estimate[1] << "," << estimate[2] << "\n";
        }
    }
}

void generate_policy_indicators(const std::filesystem::path& data_dir)
{
    using IS = PaperSecirState;

    auto model            = make_policy_secir_model(data_dir);
    const FP t0           = 0.0;
    const FP tmax         = 60.0;
    const FP dt           = 0.25;
    const FP popsize      = model.populations.get_total();
    const auto start_date = read_germany_start_date(data_dir / "germany_secir_initial_summary.csv");

    const auto result        = simulate_policy_flows_adaptively(model, t0, tmax, dt);
    const auto& compartments = result[0];
    const auto& flows        = result[1];

    auto csv = open_csv(data_dir / "secir_policy_indicators.csv");
    write_header(csv, {"day", "date", "npi_phase", "contact_multiplier", "hospital_occupancy", "icu_occupancy",
                       "daily_infections", "daily_hospital_admissions", "daily_icu_admissions", "daily_deaths",
                       "incidence_7d_per_100k", "hospital_admissions_7d_per_100k", "deaths_7d_per_100k"});

    for (Eigen::Index day = 0; day <= static_cast<Eigen::Index>(tmax); ++day) {
        const FP current_day  = static_cast<FP>(day);
        const FP previous_day = static_cast<FP>(std::max<Eigen::Index>(0, day - 1));
        const FP previous_7d  = static_cast<FP>(std::max<Eigen::Index>(0, day - 7));

        const auto compartments_day = mio::linear_interpolation(current_day, compartments);
        const auto flows_day        = mio::linear_interpolation(current_day, flows);
        const auto flows_prev       = mio::linear_interpolation(previous_day, flows);
        const auto flows_prev_7d    = mio::linear_interpolation(previous_7d, flows);

        const FP daily_infections =
            aggregate_secir_flow_difference<IS::Susceptible, IS::Exposed>(model, flows_day, flows_prev);
        const FP daily_hosp =
            aggregate_secir_flow_difference<IS::InfectedSymptoms, IS::InfectedSevere>(model, flows_day, flows_prev);
        const FP daily_icu =
            aggregate_secir_flow_difference<IS::InfectedSevere, IS::InfectedCritical>(model, flows_day, flows_prev);
        const FP daily_deaths =
            aggregate_secir_flow_difference<IS::InfectedCritical, IS::Dead>(model, flows_day, flows_prev);

        const FP infections_7d =
            aggregate_secir_flow_difference<IS::Susceptible, IS::Exposed>(model, flows_day, flows_prev_7d);
        const FP hosp_7d =
            aggregate_secir_flow_difference<IS::InfectedSymptoms, IS::InfectedSevere>(model, flows_day, flows_prev_7d);
        const FP deaths_7d =
            aggregate_secir_flow_difference<IS::InfectedCritical, IS::Dead>(model, flows_day, flows_prev_7d);

        csv << day << "," << date_plus(start_date, static_cast<int>(day)) << ","
            << npi_phase_for_day(static_cast<int>(day)) << "," << contact_multiplier_for_day(static_cast<int>(day))
            << "," << aggregate_secir_state_value(model, compartments_day, IS::InfectedSevere) << ","
            << aggregate_secir_state_value(model, compartments_day, IS::InfectedCritical) << "," << daily_infections
            << "," << daily_hosp << "," << daily_icu << "," << daily_deaths << "," << 100000.0 * infections_7d / popsize
            << "," << 100000.0 * hosp_7d / popsize << "," << 100000.0 * deaths_7d / popsize << "\n";
    }
}

void generate_flow_decomposition(const std::filesystem::path& data_dir)
{
    using IS = PaperSecirState;

    auto model            = make_policy_secir_model(data_dir);
    const FP t0           = 0.0;
    const FP tmax         = 60.0;
    const FP dt           = 0.25;
    const FP popsize      = model.populations.get_total();
    const auto start_date = read_germany_start_date(data_dir / "germany_secir_initial_summary.csv");

    const auto result        = simulate_policy_flows_adaptively(model, t0, tmax, dt);
    const auto& compartments = result[0];
    const auto& flows        = result[1];

    auto csv = open_csv(data_dir / "secir_flow_decomposition.csv");
    write_header(csv, {"day", "date", "symptomatic_occupancy", "symptomatic_inflow", "symptomatic_to_severe",
                       "symptomatic_to_recovered", "symptomatic_net_change", "severe_occupancy", "severe_inflow",
                       "severe_to_icu", "severe_to_recovered", "daily_infections_per_100k"});

    for (Eigen::Index day = 0; day <= static_cast<Eigen::Index>(tmax); ++day) {
        const FP current_day  = static_cast<FP>(day);
        const FP previous_day = static_cast<FP>(std::max<Eigen::Index>(0, day - 1));

        const auto compartments_day = mio::linear_interpolation(current_day, compartments);
        const auto flows_day        = mio::linear_interpolation(current_day, flows);
        const auto flows_prev       = mio::linear_interpolation(previous_day, flows);

        const FP symptomatic_occupancy =
            aggregate_secir_state_value(model, compartments_day, IS::InfectedSymptoms);
        const FP symptomatic_inflow =
            aggregate_secir_flow_difference<IS::InfectedNoSymptoms, IS::InfectedSymptoms>(model, flows_day, flows_prev);
        const FP symptomatic_to_severe =
            aggregate_secir_flow_difference<IS::InfectedSymptoms, IS::InfectedSevere>(model, flows_day, flows_prev);
        const FP symptomatic_to_recovered =
            aggregate_secir_flow_difference<IS::InfectedSymptoms, IS::Recovered>(model, flows_day, flows_prev);

        const FP severe_occupancy = aggregate_secir_state_value(model, compartments_day, IS::InfectedSevere);
        const FP severe_inflow    = symptomatic_to_severe;
        const FP severe_to_icu =
            aggregate_secir_flow_difference<IS::InfectedSevere, IS::InfectedCritical>(model, flows_day, flows_prev);
        const FP severe_to_recovered =
            aggregate_secir_flow_difference<IS::InfectedSevere, IS::Recovered>(model, flows_day, flows_prev);
        const FP daily_infections =
            aggregate_secir_flow_difference<IS::Susceptible, IS::Exposed>(model, flows_day, flows_prev);

        csv << day << "," << date_plus(start_date, static_cast<int>(day)) << "," << symptomatic_occupancy << ","
            << symptomatic_inflow << "," << symptomatic_to_severe << "," << symptomatic_to_recovered << ","
            << symptomatic_inflow - symptomatic_to_severe - symptomatic_to_recovered << "," << severe_occupancy << ","
            << severe_inflow << "," << severe_to_icu << "," << severe_to_recovered << ","
            << 100000.0 * daily_infections / popsize << "\n";
    }
}

class FlowlessSeirModel
    : public mio::CompartmentalModel<FP, mio::oseir::InfectionState,
                                     mio::Populations<FP, mio::AgeGroup, mio::oseir::InfectionState>,
                                     mio::oseir::Parameters<FP>>
{
    using Base = mio::CompartmentalModel<FP, mio::oseir::InfectionState,
                                         mio::Populations<FP, mio::AgeGroup, mio::oseir::InfectionState>,
                                         mio::oseir::Parameters<FP>>;

public:
    FlowlessSeirModel(int num_agegroups)
        : Base(typename Base::Populations({mio::AgeGroup(num_agegroups), mio::oseir::InfectionState::Count}),
               typename Base::ParameterSet(mio::AgeGroup(num_agegroups)))
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                         Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
    {
        auto const& params    = this->parameters;
        const auto age_groups = mio::reduce_index<mio::Index<mio::AgeGroup>>(this->populations.size());

        for (auto i : mio::make_index_range(age_groups)) {
            const size_t Si = this->populations.get_flat_index({i, mio::oseir::InfectionState::Susceptible});
            const size_t Ei = this->populations.get_flat_index({i, mio::oseir::InfectionState::Exposed});
            const size_t Ii = this->populations.get_flat_index({i, mio::oseir::InfectionState::Infected});
            const size_t Ri = this->populations.get_flat_index({i, mio::oseir::InfectionState::Recovered});

            for (auto j : mio::make_index_range(age_groups)) {
                const size_t Sj = this->populations.get_flat_index({j, mio::oseir::InfectionState::Susceptible});
                const size_t Ej = this->populations.get_flat_index({j, mio::oseir::InfectionState::Exposed});
                const size_t Ij = this->populations.get_flat_index({j, mio::oseir::InfectionState::Infected});
                const size_t Rj = this->populations.get_flat_index({j, mio::oseir::InfectionState::Recovered});

                const FP Nj    = pop[Sj] + pop[Ej] + pop[Ij] + pop[Rj];
                const FP divNj = (Nj < mio::Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nj);
                const FP coeffStoE =
                    params.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(
                        mio::SimulationTime<FP>(t))(i.get(), j.get()) *
                    params.template get<mio::oseir::TransmissionProbabilityOnContact<FP>>()[i] * divNj;

                dydt[Si] -= y[Si] * pop[Ij] * coeffStoE;
                dydt[Ei] += y[Si] * pop[Ij] * coeffStoE;
            }

            dydt[Ei] -= (1.0 / params.template get<mio::oseir::TimeExposed<FP>>()[i]) * y[Ei];
            dydt[Ii] += (1.0 / params.template get<mio::oseir::TimeExposed<FP>>()[i]) * y[Ei];
            dydt[Ii] -= (1.0 / params.template get<mio::oseir::TimeInfected<FP>>()[i]) * y[Ii];
            dydt[Ri] += (1.0 / params.template get<mio::oseir::TimeInfected<FP>>()[i]) * y[Ii];
        }
    }
};

template <class Model>
void setup_seir_runtime_model(Model& model)
{
    const FP total_population = 10000.0;
    const auto num_groups     = model.parameters.get_num_groups();
    for (mio::AgeGroup i = 0; i < num_groups; i++) {
        model.populations[{i, mio::oseir::InfectionState::Exposed}]   = 100.0 / static_cast<size_t>(num_groups);
        model.populations[{i, mio::oseir::InfectionState::Infected}]  = 100.0 / static_cast<size_t>(num_groups);
        model.populations[{i, mio::oseir::InfectionState::Recovered}] = 100.0 / static_cast<size_t>(num_groups);
        model.populations[{i, mio::oseir::InfectionState::Susceptible}] =
            total_population / static_cast<size_t>(num_groups) -
            model.populations[{i, mio::oseir::InfectionState::Exposed}] -
            model.populations[{i, mio::oseir::InfectionState::Infected}] -
            model.populations[{i, mio::oseir::InfectionState::Recovered}];
    }
    model.parameters.template set<mio::oseir::TimeExposed<FP>>(5.2);
    model.parameters.template set<mio::oseir::TimeInfected<FP>>(6.0);
    model.parameters.template set<mio::oseir::TransmissionProbabilityOnContact<FP>>(0.04);
    auto& contact_matrix = model.parameters.template get<mio::oseir::ContactPatterns<FP>>();
    contact_matrix.get_cont_freq_mat()[0].get_baseline().setConstant(10.0);
}

void setup_secir_runtime_model(PaperSecirModel& model)
{
    const FP total_population = 10000.0;
    const auto num_groups     = model.parameters.get_num_groups();
    model.populations.set_total(total_population);
    auto& contact_matrix = model.parameters.template get<mio::osecir::ContactPatterns<FP>>();
    contact_matrix.get_cont_freq_mat()[0].get_baseline().setConstant(8.0);

    for (mio::AgeGroup i = 0; i < num_groups; i++) {
        const FP group_population = total_population / static_cast<size_t>(num_groups);
        model.parameters.template get<mio::osecir::TimeExposed<FP>>()[i]                       = 3.2;
        model.parameters.template get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[i]            = 2.0;
        model.parameters.template get<mio::osecir::TimeInfectedSymptoms<FP>>()[i]              = 6.0;
        model.parameters.template get<mio::osecir::TimeInfectedSevere<FP>>()[i]                = 8.0;
        model.parameters.template get<mio::osecir::TimeInfectedCritical<FP>>()[i]              = 7.0;
        model.parameters.template get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[i]  = 0.04;
        model.parameters.template get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[i]    = 0.7;
        model.parameters.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[i]    = 0.25;
        model.parameters.template get<mio::osecir::RiskOfInfectionFromSymptomatic<FP>>()[i]    = 0.2;
        model.parameters.template get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>>()[i] = 0.4;
        model.parameters.template get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[i]         = 0.08;
        model.parameters.template get<mio::osecir::CriticalPerSevere<FP>>()[i]                 = 0.2;
        model.parameters.template get<mio::osecir::DeathsPerSevere<FP>>()[i]                   = 0.0;
        model.parameters.template get<mio::osecir::DeathsPerCritical<FP>>()[i]                 = 0.25;

        model.populations[{i, PaperSecirState::Exposed}] = 40.0 / static_cast<size_t>(num_groups);
        model.populations[{i, PaperSecirState::InfectedNoSymptoms}] =
            25.0 / static_cast<size_t>(num_groups);
        model.populations[{i, PaperSecirState::InfectedSymptoms}] = 20.0 / static_cast<size_t>(num_groups);
        model.populations[{i, PaperSecirState::InfectedSevere}]   = 4.0 / static_cast<size_t>(num_groups);
        model.populations[{i, PaperSecirState::InfectedCritical}] = 1.0 / static_cast<size_t>(num_groups);
        model.populations[{i, PaperSecirState::Recovered}]        = 250.0 / static_cast<size_t>(num_groups);
        model.populations[{i, PaperSecirState::Dead}]             = 0.0;
        const FP assigned = model.populations[{i, PaperSecirState::Exposed}] +
                            model.populations[{i, PaperSecirState::InfectedNoSymptoms}] +
                            model.populations[{i, PaperSecirState::InfectedSymptoms}] +
                            model.populations[{i, PaperSecirState::InfectedSevere}] +
                            model.populations[{i, PaperSecirState::InfectedCritical}] +
                            model.populations[{i, PaperSecirState::Recovered}] +
                            model.populations[{i, PaperSecirState::Dead}];
        model.populations[{i, PaperSecirState::Susceptible}] = group_population - assigned;
    }
    model.parameters.template get<mio::osecir::TestAndTraceCapacity<FP>>() = 10000.0;
    model.apply_constraints();
}

void setup_secirts_runtime_model(mio::osecirts::Model<FP>& model)
{
    const FP total_population = 10000.0;
    const auto num_groups     = model.parameters.get_num_groups();
    model.populations.set_total(total_population);
    auto& contact_matrix = model.parameters.template get<mio::osecirts::ContactPatterns<FP>>();
    contact_matrix.get_cont_freq_mat()[0].get_baseline().setConstant(8.0);

    const size_t num_days = 700;
    model.parameters.template get<mio::osecirts::DailyPartialVaccinations<FP>>().resize(mio::SimulationDay(num_days));
    model.parameters.template get<mio::osecirts::DailyFullVaccinations<FP>>().resize(mio::SimulationDay(num_days));
    model.parameters.template get<mio::osecirts::DailyBoosterVaccinations<FP>>().resize(mio::SimulationDay(num_days));

    for (mio::AgeGroup i = 0; i < num_groups; i++) {
        const FP group_population = total_population / static_cast<size_t>(num_groups);
        model.parameters.template get<mio::osecirts::TimeExposed<FP>>()[i]                                     = 3.2;
        model.parameters.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i]                          = 2.0;
        model.parameters.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i]                            = 6.0;
        model.parameters.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i]                              = 8.0;
        model.parameters.template get<mio::osecirts::TimeInfectedCritical<FP>>()[i]                            = 7.0;
        model.parameters.template get<mio::osecirts::TimeWaningPartialImmunity<FP>>()[i]                       = 180.0;
        model.parameters.template get<mio::osecirts::TimeWaningImprovedImmunity<FP>>()[i]                      = 240.0;
        model.parameters.template get<mio::osecirts::TimeTemporaryImmunityPI<FP>>()[i]                         = 120.0;
        model.parameters.template get<mio::osecirts::TimeTemporaryImmunityII<FP>>()[i]                         = 180.0;
        model.parameters.template get<mio::osecirts::TransmissionProbabilityOnContact<FP>>()[i]                = 0.04;
        model.parameters.template get<mio::osecirts::RelativeTransmissionNoSymptoms<FP>>()[i]                  = 0.7;
        model.parameters.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]                  = 0.25;
        model.parameters.template get<mio::osecirts::RiskOfInfectionFromSymptomatic<FP>>()[i]                  = 0.2;
        model.parameters.template get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<FP>>()[i]               = 0.4;
        model.parameters.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i]                       = 0.08;
        model.parameters.template get<mio::osecirts::CriticalPerSevere<FP>>()[i]                               = 0.2;
        model.parameters.template get<mio::osecirts::DeathsPerSevere<FP>>()[i]                                 = 0.01;
        model.parameters.template get<mio::osecirts::DeathsPerCritical<FP>>()[i]                               = 0.25;
        model.parameters.template get<mio::osecirts::ReducExposedPartialImmunity<FP>>()[i]                     = 0.65;
        model.parameters.template get<mio::osecirts::ReducExposedImprovedImmunity<FP>>()[i]                    = 0.35;
        model.parameters.template get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<FP>>()[i]            = 0.75;
        model.parameters.template get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<FP>>()[i]           = 0.45;
        model.parameters.template get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i]  = 0.55;
        model.parameters.template get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i] = 0.25;
        model.parameters.template get<mio::osecirts::ReducTimeInfectedMild<FP>>()[i]                           = 0.8;
        model.parameters.template get<mio::osecirts::InfectiousnessNewVariant<FP>>()[i]                        = 1.0;

        for (auto d = mio::SimulationDay(0); d < mio::SimulationDay(num_days); ++d) {
            model.parameters.template get<mio::osecirts::DailyPartialVaccinations<FP>>()[{i, d}] = 0.0;
            model.parameters.template get<mio::osecirts::DailyFullVaccinations<FP>>()[{i, d}]    = 0.0;
            model.parameters.template get<mio::osecirts::DailyBoosterVaccinations<FP>>()[{i, d}] = 0.0;
        }

        using IS                                                    = mio::osecirts::InfectionState;
        model.populations[{i, IS::ExposedNaive}]                    = 40.0 / static_cast<size_t>(num_groups);
        model.populations[{i, IS::InfectedNoSymptomsNaive}]         = 25.0 / static_cast<size_t>(num_groups);
        model.populations[{i, IS::InfectedSymptomsNaive}]           = 20.0 / static_cast<size_t>(num_groups);
        model.populations[{i, IS::InfectedSevereNaive}]             = 4.0 / static_cast<size_t>(num_groups);
        model.populations[{i, IS::InfectedCriticalNaive}]           = 1.0 / static_cast<size_t>(num_groups);
        model.populations[{i, IS::TemporaryImmunePartialImmunity}]  = 100.0 / static_cast<size_t>(num_groups);
        model.populations[{i, IS::TemporaryImmuneImprovedImmunity}] = 150.0 / static_cast<size_t>(num_groups);
        const FP assigned =
            model.populations[{i, IS::ExposedNaive}] + model.populations[{i, IS::InfectedNoSymptomsNaive}] +
            model.populations[{i, IS::InfectedSymptomsNaive}] + model.populations[{i, IS::InfectedSevereNaive}] +
            model.populations[{i, IS::InfectedCriticalNaive}] +
            model.populations[{i, IS::TemporaryImmunePartialImmunity}] +
            model.populations[{i, IS::TemporaryImmuneImprovedImmunity}];
        model.populations[{i, IS::SusceptibleNaive}] = group_population - assigned;
    }
    model.parameters.template get<mio::osecirts::ICUCapacity<FP>>()          = 1000.0;
    model.parameters.template get<mio::osecirts::TestAndTraceCapacity<FP>>() = 10000.0;
    model.parameters.template get<mio::osecirts::StartDay<FP>>()             = 0.0;
    model.parameters.template get<mio::osecirts::StartDayNewVariant<FP>>()   = std::numeric_limits<FP>::max();
    model.parameters.apply_constraints();
}

struct RuntimeStats {
    FP mean_ms   = 0.0;
    FP median_ms = 0.0;
    FP min_ms    = 0.0;
    FP max_ms    = 0.0;
    FP checksum  = 0.0;
};

struct PairedRuntimeStats {
    RuntimeStats compartment;
    RuntimeStats flow;
    FP ratio_mean        = 0.0;
    FP ratio_median      = 0.0;
    FP ratio_q10         = 0.0;
    FP ratio_q90         = 0.0;
    FP ratio_min         = 0.0;
    FP ratio_max         = 0.0;
    FP difference_mean   = 0.0;
    FP difference_median = 0.0;
    FP difference_q10    = 0.0;
    FP difference_q90    = 0.0;
    FP difference_min    = 0.0;
    FP difference_max    = 0.0;
};

struct PairedRuntimeMeasurement {
    PairedRuntimeStats stats;
    std::vector<FP> compartment_ms;
    std::vector<FP> flow_ms;
    std::vector<FP> runtime_ratios;
    std::vector<FP> runtime_differences_ms;
};

struct DirectComparisonDiagnostics {
    size_t time_points         = 0;
    FP max_abs_difference      = 0.0;
    FP max_population_fraction = 0.0;
};

template <class Model>
DirectComparisonDiagnostics get_direct_comparison_diagnostics(const Model& model, FP t0, FP tmax, FP dt)
{
    auto compartment_simulation = make_compartment_runtime_simulation(t0, dt, model);
    auto flow_simulation        = make_flow_runtime_simulation(t0, dt, model);
    compartment_simulation.advance(tmax);
    flow_simulation.advance(tmax);
    const auto& compartment_result = compartment_simulation.get_result();
    const auto& flow_result        = flow_simulation.get_result();
    if (compartment_result.get_num_time_points() != flow_result.get_num_time_points()) {
        throw std::runtime_error("Fixed-step comparison produced different output grids.");
    }

    FP max_difference = 0.0;
    for (Eigen::Index i = 0; i < compartment_result.get_num_time_points(); ++i) {
        max_difference = std::max(max_difference,
                                  (compartment_result.get_value(i) - flow_result.get_value(i)).cwiseAbs().maxCoeff());
    }
    const FP population = model.get_initial_values().sum();

    return {static_cast<size_t>(compartment_result.get_num_time_points()), max_difference, max_difference / population};
}

template <class SetupFn, class Fn>
RuntimeStats measure_runtime_ms(size_t repetitions, size_t trajectories_per_block, SetupFn&& setup, Fn&& fn)
{
    if (repetitions == 0 || trajectories_per_block == 0) {
        throw std::invalid_argument("Runtime benchmark repetitions and trajectories per block must be positive.");
    }

    std::vector<FP> run_times;
    run_times.reserve(repetitions);

    FP checksum = 0.0;
    {
        auto warmup = setup();
        checksum    = fn(warmup);
    }

    for (size_t i = 0; i < repetitions; ++i) {
        FP block_time = 0.0;
        for (size_t trajectory = 0; trajectory < trajectories_per_block; ++trajectory) {
            auto item        = setup();
            const auto start = std::chrono::steady_clock::now();
            checksum += fn(item);
            const auto stop = std::chrono::steady_clock::now();
            block_time += std::chrono::duration<FP, std::milli>(stop - start).count();
        }
        run_times.push_back(block_time);
    }

    auto sorted_times = run_times;
    std::sort(sorted_times.begin(), sorted_times.end());
    const auto middle = sorted_times.size() / 2;
    const FP median =
        sorted_times.size() % 2 == 0 ? 0.5 * (sorted_times[middle - 1] + sorted_times[middle]) : sorted_times[middle];
    const FP mean = std::accumulate(run_times.begin(), run_times.end(), FP(0.0)) / static_cast<FP>(run_times.size());

    return {mean, median, sorted_times.front(), sorted_times.back(), checksum};
}

FP quantile(const std::vector<FP>& sorted_values, FP probability)
{
    const FP position = probability * static_cast<FP>(sorted_values.size() - 1);
    const auto lower  = static_cast<size_t>(std::floor(position));
    const auto upper  = static_cast<size_t>(std::ceil(position));
    const FP weight   = position - static_cast<FP>(lower);
    return (1.0 - weight) * sorted_values[lower] + weight * sorted_values[upper];
}

RuntimeStats summarize_runtime(const std::vector<FP>& run_times, FP checksum)
{
    auto sorted_times = run_times;
    std::sort(sorted_times.begin(), sorted_times.end());
    const auto middle = sorted_times.size() / 2;
    const FP median =
        sorted_times.size() % 2 == 0 ? 0.5 * (sorted_times[middle - 1] + sorted_times[middle]) : sorted_times[middle];
    const FP mean = std::accumulate(run_times.begin(), run_times.end(), FP(0.0)) / static_cast<FP>(run_times.size());
    return {mean, median, sorted_times.front(), sorted_times.back(), checksum};
}

template <class CompartmentSetupFn, class CompartmentFn, class FlowSetupFn, class FlowFn>
PairedRuntimeMeasurement measure_paired_runtime_ms(size_t repetitions, size_t trajectories_per_block,
                                                   CompartmentSetupFn&& compartment_setup,
                                                   CompartmentFn&& compartment_fn, FlowSetupFn&& flow_setup,
                                                   FlowFn&& flow_fn)
{
    if (repetitions == 0 || trajectories_per_block == 0) {
        throw std::invalid_argument("Runtime benchmark repetitions and trajectories per block must be positive.");
    }

    PairedRuntimeMeasurement measurement;
    measurement.compartment_ms.reserve(repetitions);
    measurement.flow_ms.reserve(repetitions);
    measurement.runtime_ratios.reserve(repetitions);
    measurement.runtime_differences_ms.reserve(repetitions);

    FP compartment_checksum = 0.0;
    FP flow_checksum        = 0.0;

    const auto measure_one = [](auto& item, auto& fn, FP& checksum) {
        const auto start = std::chrono::steady_clock::now();
        checksum += fn(item);
        const auto stop = std::chrono::steady_clock::now();
        return std::chrono::duration<FP, std::milli>(stop - start).count();
    };

    {
        auto compartment_warmup = compartment_setup();
        auto flow_warmup        = flow_setup();
        compartment_checksum += compartment_fn(compartment_warmup);
        flow_checksum += flow_fn(flow_warmup);
    }

    for (size_t run = 0; run < repetitions; ++run) {
        FP compartment_time   = 0.0;
        FP flow_time          = 0.0;

        for (size_t trajectory = 0; trajectory < trajectories_per_block; ++trajectory) {
            auto compartment_item = compartment_setup();
            auto flow_item        = flow_setup();
            if ((run * trajectories_per_block + trajectory) % 2 == 0) {
                compartment_time += measure_one(compartment_item, compartment_fn, compartment_checksum);
                flow_time += measure_one(flow_item, flow_fn, flow_checksum);
            }
            else {
                flow_time += measure_one(flow_item, flow_fn, flow_checksum);
                compartment_time += measure_one(compartment_item, compartment_fn, compartment_checksum);
            }
        }

        measurement.compartment_ms.push_back(compartment_time);
        measurement.flow_ms.push_back(flow_time);
        measurement.runtime_ratios.push_back(flow_time / compartment_time);
        measurement.runtime_differences_ms.push_back(flow_time - compartment_time);
    }

    auto sorted_ratios = measurement.runtime_ratios;
    std::sort(sorted_ratios.begin(), sorted_ratios.end());
    auto sorted_differences = measurement.runtime_differences_ms;
    std::sort(sorted_differences.begin(), sorted_differences.end());
    const FP ratio_mean =
        std::accumulate(measurement.runtime_ratios.begin(), measurement.runtime_ratios.end(), FP(0.0)) /
        static_cast<FP>(measurement.runtime_ratios.size());
    const FP difference_mean =
        std::accumulate(measurement.runtime_differences_ms.begin(), measurement.runtime_differences_ms.end(), FP(0.0)) /
        static_cast<FP>(measurement.runtime_differences_ms.size());

    measurement.stats = {summarize_runtime(measurement.compartment_ms, compartment_checksum),
                         summarize_runtime(measurement.flow_ms, flow_checksum),
                         ratio_mean,
                         quantile(sorted_ratios, 0.5),
                         quantile(sorted_ratios, 0.1),
                         quantile(sorted_ratios, 0.9),
                         sorted_ratios.front(),
                         sorted_ratios.back(),
                         difference_mean,
                         quantile(sorted_differences, 0.5),
                         quantile(sorted_differences, 0.1),
                         quantile(sorted_differences, 0.9),
                         sorted_differences.front(),
                         sorted_differences.back()};
    return measurement;
}

void write_direct_runtime_runs(std::ofstream& csv, const std::string& model_name, size_t num_age_groups,
                               const PairedRuntimeMeasurement& measurement)
{
    for (size_t run = 0; run < measurement.runtime_ratios.size(); ++run) {
        csv << model_name << "," << num_age_groups << "," << run << "," << measurement.compartment_ms[run] << ","
            << measurement.flow_ms[run] << "," << measurement.runtime_ratios[run] << ","
            << measurement.runtime_differences_ms[run] << "\n";
    }
}

void write_auxiliary_runtime_runs(std::ofstream& csv, const std::string& model_name, size_t num_age_groups,
                                  const std::string& method, size_t tracked_transitions,
                                  const PairedRuntimeMeasurement& measurement)
{
    for (size_t run = 0; run < measurement.runtime_ratios.size(); ++run) {
        csv << model_name << "," << num_age_groups << "," << method << "," << tracked_transitions << "," << run << ","
            << measurement.compartment_ms[run] << "," << measurement.flow_ms[run] << ","
            << measurement.runtime_ratios[run] << "," << measurement.runtime_differences_ms[run] << "\n";
    }
}

void write_naive_auxiliary_runtime_runs(std::ofstream& csv, const std::string& model_name, size_t num_age_groups,
                                        const std::string& method, const std::string& selection_order,
                                        size_t tracked_transitions, const PairedRuntimeMeasurement& measurement)
{
    for (size_t run = 0; run < measurement.runtime_ratios.size(); ++run) {
        csv << model_name << "," << num_age_groups << "," << method << "," << selection_order << ","
            << tracked_transitions << "," << run << "," << measurement.compartment_ms[run] << ","
            << measurement.flow_ms[run] << "," << measurement.runtime_ratios[run] << ","
            << measurement.runtime_differences_ms[run] << "\n";
    }
}

void write_naive_auxiliary_summary_row(std::ofstream& csv, const std::string& model_name, size_t num_age_groups,
                                       const std::string& method, const std::string& selection_order,
                                       size_t tracked_transitions, size_t integrated_dim, size_t compartment_dim,
                                       size_t flow_dim, size_t repetitions, const PairedRuntimeMeasurement& measurement)
{
    const auto& paired   = measurement.stats;
    const auto& baseline = paired.compartment;
    const auto& runtime  = paired.flow;
    csv << model_name << "," << num_age_groups << "," << method << "," << selection_order << ","
        << tracked_transitions << "," << integrated_dim << "," << compartment_dim << "," << flow_dim << ","
        << runtime.mean_ms << ","
        << runtime.median_ms << "," << runtime.min_ms << "," << runtime.max_ms << "," << repetitions << ",1,"
        << runtime.checksum << "," << baseline.mean_ms << "," << baseline.median_ms << "," << paired.ratio_mean << ","
        << paired.ratio_median << "," << paired.ratio_q10 << "," << paired.ratio_q90 << "," << paired.ratio_min << ","
        << paired.ratio_max << "," << baseline.checksum << "\n";
}

template <class Model>
PairedRuntimeMeasurement measure_direct_runtime_pair(const Model& model, FP t0, FP tmax, FP dt, size_t repetitions,
                                                     size_t trajectories_per_block)
{
    const auto output_points = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;
    return measure_paired_runtime_ms(
        repetitions, trajectories_per_block,
        [&]() {
            auto simulation = make_compartment_runtime_simulation(t0, dt, model);
            simulation.get_result().reserve(output_points);
            return simulation;
        },
        [&](auto& simulation) {
            simulation.advance(tmax);
            return simulation.get_result().get_last_value().sum();
        },
        [&]() {
            auto simulation = make_flow_runtime_simulation(t0, dt, model);
            simulation.get_result().reserve(output_points);
            simulation.get_flows().reserve(output_points);
            return simulation;
        },
        [&](auto& simulation) {
            simulation.advance(tmax);
            return simulation.get_result().get_last_value().sum() + simulation.get_flows().get_last_value().sum();
        });
}

template <class Model>
void write_auxiliary_scaling_rows(std::ofstream& csv, std::ofstream& run_csv, const std::string& model_name,
                                  const Model& base_model, size_t num_age_groups, FP t0, FP tmax, FP dt,
                                  size_t repetitions, size_t trajectories_per_block)
{
    const size_t compartment_dim = static_cast<size_t>(base_model.get_initial_values().size());
    const size_t flow_dim        = static_cast<size_t>(base_model.get_initial_flows().size());
    const auto output_points     = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;

    const auto make_compartment_simulation = [&](const auto& runtime_model) {
        auto simulation = make_compartment_runtime_simulation(t0, dt, runtime_model);
        simulation.get_result().reserve(output_points);
        return simulation;
    };
    const auto advance_compartment_simulation = [&](auto& simulation) {
        simulation.advance(tmax);
        return simulation.get_result().get_last_value().sum();
    };
    const auto make_flow_simulation = [&]() {
        auto simulation = make_flow_runtime_simulation(t0, dt, base_model);
        simulation.get_result().reserve(output_points);
        simulation.get_flows().reserve(output_points);
        return simulation;
    };
    const auto advance_flow_simulation = [&](auto& simulation) {
        simulation.advance(tmax);
        return simulation.get_result().get_last_value().sum() + simulation.get_flows().get_last_value().sum();
    };

    const auto flow_measurement = measure_paired_runtime_ms(
        repetitions, trajectories_per_block, make_flow_simulation, advance_flow_simulation, make_flow_simulation,
        advance_flow_simulation);
    const auto& flow_pair       = flow_measurement.stats;
    const auto& flow_runtime    = flow_pair.compartment;
    write_auxiliary_runtime_runs(run_csv, model_name, num_age_groups, "Flow based", flow_dim, flow_measurement);
    csv << model_name << "," << num_age_groups << ",Flow based," << flow_dim << "," << flow_dim << ","
        << compartment_dim << "," << flow_dim << "," << flow_runtime.mean_ms << "," << flow_runtime.median_ms << ","
        << flow_runtime.min_ms << "," << flow_runtime.max_ms << "," << repetitions << ",1,Flow based,"
        << flow_runtime.checksum << "," << flow_runtime.mean_ms << "," << flow_runtime.median_ms << ","
        << flow_pair.ratio_mean << "," << flow_pair.ratio_median << "," << flow_pair.ratio_q10 << ","
        << flow_pair.ratio_q90 << "," << flow_pair.ratio_min << "," << flow_pair.ratio_max << ","
        << flow_runtime.checksum << "\n";

    for (size_t num_auxiliary_states = 0; num_auxiliary_states <= flow_dim; ++num_auxiliary_states) {
        AuxiliaryTrackingModel<Model> auxiliary_model(base_model, num_auxiliary_states);
        const auto paired_measurement = measure_paired_runtime_ms(
            repetitions, trajectories_per_block, make_flow_simulation, advance_flow_simulation,
            [&]() {
                return make_compartment_simulation(auxiliary_model);
            },
            advance_compartment_simulation);
        const auto& paired   = paired_measurement.stats;
        const auto& baseline = paired.compartment;
        const auto& runtime  = paired.flow;
        write_auxiliary_runtime_runs(run_csv, model_name, num_age_groups, "Auxiliary ODE", num_auxiliary_states,
                                     paired_measurement);
        csv << model_name << "," << num_age_groups << ",Auxiliary ODE," << num_auxiliary_states << ","
            << auxiliary_model.get_initial_values().size() << "," << compartment_dim << "," << flow_dim << ","
            << runtime.mean_ms << "," << runtime.median_ms << "," << runtime.min_ms << "," << runtime.max_ms << ","
            << repetitions << ",1,Flow based," << runtime.checksum << "," << baseline.mean_ms << ","
            << baseline.median_ms << "," << paired.ratio_mean << "," << paired.ratio_median << "," << paired.ratio_q10
            << "," << paired.ratio_q90 << "," << paired.ratio_min << "," << paired.ratio_max << "," << baseline.checksum
            << "\n";
    }
}

void generate_runtime_overhead(const std::filesystem::path& data_dir)
{
    mio::set_log_level(mio::LogLevel::critical);
    const FP t0           = runtime_t0;
    const FP tmax         = runtime_tmax;
    const FP dt           = runtime_dt;
    const size_t repeats  = get_environment_size("FLOW_PAPER_DIRECT_REPETITIONS", runtime_repetitions);
    const size_t trajectories_per_block =
        get_environment_size("FLOW_PAPER_TRAJECTORIES_PER_BLOCK", runtime_trajectories_per_block);
    const auto age_groups = get_environment_ints("FLOW_PAPER_DIRECT_AGE_GROUPS", {1, 2, 3, 4, 5, 6});

    auto config = open_csv(data_dir / "runtime_direct_config.csv");
    write_header(config, {"integrator", "t0", "tmax", "dt", "stages", "repetitions", "warmup_runs", "timed_region",
                          "pairing", "output_allocation", "trajectories_per_timing_block"});
    config << "classical RK4," << t0 << "," << tmax << "," << dt << ",4," << repeats
           << ",1,advance only,alternating paired trajectories,preallocated," << trajectories_per_block << "\n";

    auto csv = open_csv(data_dir / "runtime_overhead.csv");
    write_header(csv, {"model", "age_groups", "mode", "compartment_dim", "flow_dim", "mean_ms", "median_ms", "min_ms",
                       "max_ms", "repetitions", "warmup_runs", "checksum"});
    auto model_summary = open_csv(data_dir / "runtime_model_overhead.csv");
    write_header(model_summary, {"model",
                                 "age_groups",
                                 "compartment_dim",
                                 "flow_dim",
                                 "compartment_mean_ms",
                                 "compartment_median_ms",
                                 "flow_mean_ms",
                                 "flow_median_ms",
                                 "runtime_ratio_mean",
                                 "runtime_ratio_median",
                                 "runtime_ratio_q10",
                                 "runtime_ratio_q90",
                                 "runtime_ratio_min",
                                 "runtime_ratio_max",
                                 "runtime_difference_mean_ms",
                                 "runtime_difference_median_ms",
                                 "runtime_difference_q10_ms",
                                 "runtime_difference_q90_ms",
                                 "runtime_difference_min_ms",
                                 "runtime_difference_max_ms",
                                 "overhead_pct",
                                 "flow_to_compartment_dim_ratio",
                                 "time_points",
                                 "max_abs_difference",
                                 "max_population_fraction",
                                 "repetitions",
                                 "warmup_runs"});
    auto run_csv = open_csv(data_dir / "runtime_direct_runs.csv");
    write_header(run_csv,
                 {"model", "age_groups", "run", "compartment_ms", "flow_ms", "runtime_ratio", "runtime_difference_ms"});

    for (int num_groups : age_groups) {
        if (!environment_contains("FLOW_PAPER_DIRECT_MODELS", "SEIR")) {
            continue;
        }
        FlowlessSeirModel flowless(num_groups);
        setup_seir_runtime_model(flowless);

        mio::oseir::Model<FP> flow_aware(num_groups);
        setup_seir_runtime_model(flow_aware);
        const auto diagnostics = get_direct_comparison_diagnostics(flow_aware, t0, tmax, dt);

        const auto handwritten = measure_runtime_ms(
            repeats, trajectories_per_block,
            [&]() {
                auto simulation          = make_compartment_runtime_simulation(t0, dt, flowless);
                const auto output_points = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;
                simulation.get_result().reserve(output_points);
                return simulation;
            },
            [&](auto& simulation) {
                simulation.advance(tmax);
                return simulation.get_result().get_last_value().sum();
            });
        const auto paired_measurement =
            measure_direct_runtime_pair(flow_aware, t0, tmax, dt, repeats, trajectories_per_block);
        const auto& paired            = paired_measurement.stats;
        const auto& compartment       = paired.compartment;
        const auto& flow              = paired.flow;
        write_direct_runtime_runs(run_csv, "SEIR", static_cast<size_t>(num_groups), paired_measurement);

        csv << "SEIR," << num_groups << ",handwritten compartment RHS," << flowless.get_initial_values().size() << ",0,"
            << handwritten.mean_ms << "," << handwritten.median_ms << "," << handwritten.min_ms << ","
            << handwritten.max_ms << "," << repeats << ",1," << handwritten.checksum << "\n";
        csv << "SEIR," << num_groups << ",FlowModel as compartment RHS," << flow_aware.get_initial_values().size()
            << ",0," << compartment.mean_ms << "," << compartment.median_ms << "," << compartment.min_ms << ","
            << compartment.max_ms << "," << repeats << ",1," << compartment.checksum << "\n";
        csv << "SEIR," << num_groups << ",Flow based," << flow_aware.get_initial_values().size() << ","
            << flow_aware.get_initial_flows().size() << "," << flow.mean_ms << "," << flow.median_ms << ","
            << flow.min_ms << "," << flow.max_ms << "," << repeats << ",1," << flow.checksum << "\n";
        model_summary << "SEIR," << num_groups << "," << flow_aware.get_initial_values().size() << ","
                      << flow_aware.get_initial_flows().size() << "," << compartment.mean_ms << ","
                      << compartment.median_ms << "," << flow.mean_ms << "," << flow.median_ms << ","
                      << paired.ratio_mean << "," << paired.ratio_median << "," << paired.ratio_q10 << ","
                      << paired.ratio_q90 << "," << paired.ratio_min << "," << paired.ratio_max << ","
                      << paired.difference_mean << "," << paired.difference_median << "," << paired.difference_q10
                      << "," << paired.difference_q90 << "," << paired.difference_min << "," << paired.difference_max
                      << "," << 100.0 * (paired.ratio_median - 1.0) << ","
                      << static_cast<FP>(flow_aware.get_initial_flows().size()) /
                             static_cast<FP>(flow_aware.get_initial_values().size())
                      << "," << diagnostics.time_points << "," << diagnostics.max_abs_difference << ","
                      << diagnostics.max_population_fraction << "," << repeats << ",1\n";
    }

    for (int num_groups : age_groups) {
        if (!environment_contains("FLOW_PAPER_DIRECT_MODELS", "SECIR")) {
            continue;
        }
        PaperSecirModel model(num_groups);
        setup_secir_runtime_model(model);
        const auto diagnostics = get_direct_comparison_diagnostics(model, t0, tmax, dt);

        const auto paired_measurement =
            measure_direct_runtime_pair(model, t0, tmax, dt, repeats, trajectories_per_block);
        const auto& paired            = paired_measurement.stats;
        const auto& compartment       = paired.compartment;
        const auto& flow              = paired.flow;
        write_direct_runtime_runs(run_csv, "SECIR", static_cast<size_t>(num_groups), paired_measurement);
        csv << "SECIR," << num_groups << ",FlowModel as compartment RHS," << model.get_initial_values().size() << ",0,"
            << compartment.mean_ms << "," << compartment.median_ms << "," << compartment.min_ms << ","
            << compartment.max_ms << "," << repeats << ",1," << compartment.checksum << "\n";
        csv << "SECIR," << num_groups << ",Flow based," << model.get_initial_values().size() << ","
            << model.get_initial_flows().size() << "," << flow.mean_ms << "," << flow.median_ms << "," << flow.min_ms
            << "," << flow.max_ms << "," << repeats << ",1," << flow.checksum << "\n";
        model_summary << "SECIR," << num_groups << "," << model.get_initial_values().size() << ","
                      << model.get_initial_flows().size() << "," << compartment.mean_ms << "," << compartment.median_ms
                      << "," << flow.mean_ms << "," << flow.median_ms << "," << paired.ratio_mean << ","
                      << paired.ratio_median << "," << paired.ratio_q10 << "," << paired.ratio_q90 << ","
                      << paired.ratio_min << "," << paired.ratio_max << "," << paired.difference_mean << ","
                      << paired.difference_median << "," << paired.difference_q10 << "," << paired.difference_q90 << ","
                      << paired.difference_min << "," << paired.difference_max << ","
                      << 100.0 * (paired.ratio_median - 1.0) << ","
                      << static_cast<FP>(model.get_initial_flows().size()) /
                             static_cast<FP>(model.get_initial_values().size())
                      << "," << diagnostics.time_points << "," << diagnostics.max_abs_difference << ","
                      << diagnostics.max_population_fraction << "," << repeats << ",1\n";
    }

    for (int num_groups : age_groups) {
        if (!environment_contains("FLOW_PAPER_DIRECT_MODELS", "SECIRTS")) {
            continue;
        }
        mio::osecirts::Model<FP> model(num_groups);
        setup_secirts_runtime_model(model);
        const auto diagnostics = get_direct_comparison_diagnostics(model, t0, tmax, dt);

        const auto paired_measurement =
            measure_direct_runtime_pair(model, t0, tmax, dt, repeats, trajectories_per_block);
        const auto& paired            = paired_measurement.stats;
        const auto& compartment       = paired.compartment;
        const auto& flow              = paired.flow;
        write_direct_runtime_runs(run_csv, "SECIRTS", static_cast<size_t>(num_groups), paired_measurement);
        csv << "SECIRTS," << num_groups << ",FlowModel as compartment RHS," << model.get_initial_values().size()
            << ",0," << compartment.mean_ms << "," << compartment.median_ms << "," << compartment.min_ms << ","
            << compartment.max_ms << "," << repeats << ",1," << compartment.checksum << "\n";
        csv << "SECIRTS," << num_groups << ",Flow based," << model.get_initial_values().size() << ","
            << model.get_initial_flows().size() << "," << flow.mean_ms << "," << flow.median_ms << "," << flow.min_ms
            << "," << flow.max_ms << "," << repeats << ",1," << flow.checksum << "\n";
        model_summary << "SECIRTS," << num_groups << "," << model.get_initial_values().size() << ","
                      << model.get_initial_flows().size() << "," << compartment.mean_ms << "," << compartment.median_ms
                      << "," << flow.mean_ms << "," << flow.median_ms << "," << paired.ratio_mean << ","
                      << paired.ratio_median << "," << paired.ratio_q10 << "," << paired.ratio_q90 << ","
                      << paired.ratio_min << "," << paired.ratio_max << "," << paired.difference_mean << ","
                      << paired.difference_median << "," << paired.difference_q10 << "," << paired.difference_q90 << ","
                      << paired.difference_min << "," << paired.difference_max << ","
                      << 100.0 * (paired.ratio_median - 1.0) << ","
                      << static_cast<FP>(model.get_initial_flows().size()) /
                             static_cast<FP>(model.get_initial_values().size())
                      << "," << diagnostics.time_points << "," << diagnostics.max_abs_difference << ","
                      << diagnostics.max_population_fraction << "," << repeats << ",1\n";
    }
}

void generate_auxiliary_scaling(const std::filesystem::path& data_dir)
{
    const FP t0           = runtime_t0;
    const FP tmax         = runtime_tmax;
    const FP dt           = runtime_dt;
    const size_t repeats  = get_environment_size("FLOW_PAPER_AUXILIARY_REPETITIONS", runtime_repetitions);
    const size_t trajectories_per_block =
        get_environment_size("FLOW_PAPER_TRAJECTORIES_PER_BLOCK", runtime_trajectories_per_block);
    const auto age_groups = get_environment_ints("FLOW_PAPER_AUXILIARY_AGE_GROUPS", {1, 3, 6});

    auto config = open_csv(data_dir / "runtime_auxiliary_config.csv");
    write_header(config, {"integrator", "t0", "tmax", "dt", "stages", "repetitions", "warmup_runs", "timed_region",
                          "pairing", "output_allocation", "trajectories_per_timing_block", "ratio_reference"});
    config << "classical RK4," << t0 << "," << tmax << "," << dt << ",4," << repeats
           << ",1,advance only,alternating paired trajectories,preallocated," << trajectories_per_block
           << ",Flow based\n";

    auto csv = open_csv(data_dir / "auxiliary_scaling.csv");
    write_header(csv, {"model",
                       "age_groups",
                       "method",
                       "tracked_transitions",
                       "integrated_dim",
                       "compartment_dim",
                       "flow_dim",
                       "mean_ms",
                       "median_ms",
                       "min_ms",
                       "max_ms",
                       "repetitions",
                       "warmup_runs",
                       "ratio_reference",
                       "checksum",
                       "baseline_mean_ms",
                       "baseline_median_ms",
                       "runtime_ratio_mean",
                       "runtime_ratio_median",
                       "runtime_ratio_q10",
                       "runtime_ratio_q90",
                       "runtime_ratio_min",
                       "runtime_ratio_max",
                       "baseline_checksum"});
    auto run_csv = open_csv(data_dir / "runtime_auxiliary_runs.csv");
    write_header(run_csv, {"model", "age_groups", "method", "tracked_transitions", "run", "baseline_ms",
                           "comparison_ms", "runtime_ratio", "runtime_difference_ms"});

    for (int num_groups : age_groups) {
        if (environment_contains("FLOW_PAPER_AUXILIARY_MODELS", "SEIR")) {
            mio::oseir::Model<FP> seir_model(num_groups);
            setup_seir_runtime_model(seir_model);
            write_auxiliary_scaling_rows(csv, run_csv, "SEIR", seir_model, static_cast<size_t>(num_groups), t0, tmax,
                                         dt, repeats, trajectories_per_block);
        }

        if (environment_contains("FLOW_PAPER_AUXILIARY_MODELS", "SECIR")) {
            PaperSecirModel secir_model(num_groups);
            setup_secir_runtime_model(secir_model);
            write_auxiliary_scaling_rows(csv, run_csv, "SECIR", secir_model, static_cast<size_t>(num_groups), t0, tmax,
                                         dt, repeats, trajectories_per_block);
        }

        if (environment_contains("FLOW_PAPER_AUXILIARY_MODELS", "SECIRTS")) {
            mio::osecirts::Model<FP> secirts_model(num_groups);
            setup_secirts_runtime_model(secirts_model);
            write_auxiliary_scaling_rows(csv, run_csv, "SECIRTS", secirts_model, static_cast<size_t>(num_groups), t0,
                                         tmax, dt, repeats, trajectories_per_block);
        }
    }
}

void generate_naive_auxiliary_scaling(const std::filesystem::path& data_dir)
{
    mio::set_log_level(mio::LogLevel::critical);
    const FP t0              = runtime_t0;
    const FP tmax            = runtime_tmax;
    const FP dt              = runtime_dt;
    const size_t repetitions = get_environment_size("FLOW_PAPER_NAIVE_REPETITIONS", runtime_repetitions);
    const size_t trajectories_per_block =
        get_environment_size("FLOW_PAPER_TRAJECTORIES_PER_BLOCK", runtime_trajectories_per_block);
    const auto age_groups    = get_environment_ints("FLOW_PAPER_NAIVE_AGE_GROUPS", {1, 3, 6});

    auto config = open_csv(data_dir / "runtime_naive_auxiliary_config.csv");
    write_header(config, {"integrator", "t0", "tmax", "dt", "stages", "repetitions", "warmup_runs", "timed_region",
                          "pairing", "output_allocation", "trajectories_per_timing_block", "model",
                          "selection_orders", "ratio_reference"});
    config << "classical RK4," << t0 << "," << tmax << "," << dt << ",4," << repetitions
           << ",1,advance only,alternating paired trajectories,preallocated," << trajectories_per_block
           << ",SECIR,group_wise|other_first|infection_first,flow-only\n";

    auto csv = open_csv(data_dir / "naive_auxiliary_scaling.csv");
    write_header(csv, {"model",
                       "age_groups",
                       "method",
                       "selection_order",
                       "tracked_transitions",
                       "integrated_dim",
                       "compartment_dim",
                       "flow_dim",
                       "mean_ms",
                       "median_ms",
                       "min_ms",
                       "max_ms",
                       "repetitions",
                       "warmup_runs",
                       "checksum",
                       "baseline_mean_ms",
                       "baseline_median_ms",
                       "runtime_ratio_mean",
                       "runtime_ratio_median",
                       "runtime_ratio_q10",
                       "runtime_ratio_q90",
                       "runtime_ratio_min",
                       "runtime_ratio_max",
                       "baseline_checksum"});
    auto run_csv = open_csv(data_dir / "runtime_naive_auxiliary_runs.csv");
    write_header(run_csv, {"model", "age_groups", "method", "selection_order", "tracked_transitions", "run",
                           "baseline_ms", "comparison_ms", "runtime_ratio", "runtime_difference_ms"});

    for (int num_groups : age_groups) {
        PaperSecirModel model(num_groups);
        setup_secir_runtime_model(model);
        validate_recomputed_secir_rates(model);

        const size_t compartment_dim = static_cast<size_t>(model.get_initial_values().size());
        const size_t flow_dim        = static_cast<size_t>(model.get_initial_flows().size());
        const auto output_points     = static_cast<Eigen::Index>(std::ceil((tmax - t0) / dt)) + 1;

        const auto make_compartment_simulation = [&](const auto& runtime_model) {
            auto simulation = make_compartment_runtime_simulation(t0, dt, runtime_model);
            simulation.get_result().reserve(output_points);
            return simulation;
        };
        const auto advance_compartment_simulation = [&](auto& simulation) {
            simulation.advance(tmax);
            return simulation.get_result().get_last_value().sum();
        };
        const auto make_flow_simulation = [&]() {
            auto simulation = make_flow_runtime_simulation(t0, dt, model);
            simulation.get_result().reserve(output_points);
            simulation.get_flows().reserve(output_points);
            return simulation;
        };
        const auto advance_flow_simulation = [&](auto& simulation) {
            simulation.advance(tmax);
            return simulation.get_result().get_last_value().sum() + simulation.get_flows().get_last_value().sum();
        };
        const auto measure_against_flow = [&](const auto& comparison_model) {
            return measure_paired_runtime_ms(
                repetitions, trajectories_per_block, make_flow_simulation, advance_flow_simulation,
                [&]() {
                    return make_compartment_simulation(comparison_model);
                },
                advance_compartment_simulation);
        };

        std::vector<std::string> selection_orders;
        for (const std::string selection_order : {"group_wise", "other_first", "infection_first"}) {
            if (!environment_contains("FLOW_PAPER_NAIVE_SELECTION_ORDERS", selection_order)) {
                continue;
            }
            selection_orders.push_back(selection_order);
        }

        const NaiveAuxiliarySecirModel no_counter_model(model, {});
        const auto no_counter_measurement = measure_against_flow(no_counter_model);
        for (const auto& selection_order : selection_orders) {
            write_naive_auxiliary_runtime_runs(run_csv, "SECIR", static_cast<size_t>(num_groups),
                                               "Recomputed auxiliary", selection_order, 0, no_counter_measurement);
            write_naive_auxiliary_summary_row(csv, "SECIR", static_cast<size_t>(num_groups),
                                              "Recomputed auxiliary", selection_order, 0, compartment_dim,
                                              compartment_dim, flow_dim, repetitions, no_counter_measurement);
        }

        for (const auto& selection_order : selection_orders) {
            const auto ordered_transitions =
                make_secir_transition_order(static_cast<size_t>(num_groups), selection_order);
            for (size_t tracked_transitions = 1; tracked_transitions < flow_dim; ++tracked_transitions) {
                const std::vector<SecirTrackedTransition> selected_transitions(
                    ordered_transitions.begin(),
                    ordered_transitions.begin() + static_cast<std::ptrdiff_t>(tracked_transitions));
                const NaiveAuxiliarySecirModel naive_model(model, selected_transitions);
                const auto naive_measurement = measure_against_flow(naive_model);
                write_naive_auxiliary_runtime_runs(run_csv, "SECIR", static_cast<size_t>(num_groups),
                                                   "Recomputed auxiliary", selection_order, tracked_transitions,
                                                   naive_measurement);
                write_naive_auxiliary_summary_row(
                    csv, "SECIR", static_cast<size_t>(num_groups), "Recomputed auxiliary", selection_order,
                    tracked_transitions, compartment_dim + tracked_transitions, compartment_dim, flow_dim,
                    repetitions, naive_measurement);
            }
        }

        const auto all_transitions =
            make_secir_transition_order(static_cast<size_t>(num_groups), "group_wise");
        const NaiveAuxiliarySecirModel all_counter_model(model, all_transitions);
        const auto all_counter_measurement = measure_against_flow(all_counter_model);
        for (const auto& selection_order : selection_orders) {
            write_naive_auxiliary_runtime_runs(run_csv, "SECIR", static_cast<size_t>(num_groups),
                                               "Recomputed auxiliary", selection_order, flow_dim,
                                               all_counter_measurement);
            write_naive_auxiliary_summary_row(
                csv, "SECIR", static_cast<size_t>(num_groups), "Recomputed auxiliary", selection_order, flow_dim,
                compartment_dim + flow_dim, compartment_dim, flow_dim, repetitions, all_counter_measurement);
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::warn);

    std::filesystem::path data_dir = "submitted/data/numerics";
    if (argc > 1) {
        data_dir = argv[1];
    }
    const std::string section = argc > 2 ? argv[2] : "all";
    std::filesystem::create_directories(data_dir);

    if (section == "all" || section == "scenarios") {
        generate_sir_method_comparison(data_dir);
        generate_grid_sensitivity(data_dir);
        generate_secir_grid_sensitivity(data_dir);
        generate_flow_decomposition(data_dir);
        generate_policy_indicators(data_dir);
    }
    if (section == "all" || section == "runtime" || section == "runtime-direct") {
        generate_runtime_overhead(data_dir);
    }
    if (section == "all" || section == "runtime" || section == "runtime-auxiliary") {
        generate_auxiliary_scaling(data_dir);
    }
    if (section == "all" || section == "runtime-naive-auxiliary") {
        generate_naive_auxiliary_scaling(data_dir);
    }
    if (section != "all" && section != "scenarios" && section != "runtime" && section != "runtime-direct" &&
        section != "runtime-auxiliary" && section != "runtime-naive-auxiliary") {
        std::cerr << "Unknown section " << section << "\n";
        return 1;
    }

    std::cout << "Wrote numerical data to " << data_dir << "\n";
    return 0;
}
