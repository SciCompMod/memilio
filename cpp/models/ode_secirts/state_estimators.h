#include "ode_secirts/model.h"
#include "ode_secirts/infection_state.h"
#include "ode_secirts/parameters.h"
#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/math/eigen.h"
#include <Eigen/src/Core/ReturnByValue.h>

namespace mio
{
namespace osecirts
{

using ScalarType = double;
enum class CommuterType
{
    NonCommuter = 0,
    CommuterBase
};

using InfectionStateSecirtsExplicit = mio::osecirts::InfectionState;
using ParametersSecirtsExplicit     = mio::osecirts::Parameters<ScalarType>;
using FlowsSecirtsExplicit          = mio::osecirts::Flows;

template <typename FP = ScalarType>
class ModelSecirtsExplicit
    : public mio::FlowModel<FP, InfectionStateSecirtsExplicit,
                            mio::Populations<FP, mio::AgeGroup, CommuterType, InfectionStateSecirtsExplicit>,
                            ParametersSecirtsExplicit, FlowsSecirtsExplicit>
{
    using Base = mio::FlowModel<FP, InfectionStateSecirtsExplicit,
                                mio::Populations<FP, mio::AgeGroup, CommuterType, InfectionStateSecirtsExplicit>,
                                ParametersSecirtsExplicit, FlowsSecirtsExplicit>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;
    using CommuterIndex = CommuterType;

protected:
    int m_num_commuter_groups;

public:
    ModelSecirtsExplicit(int num_agegroups, int num_commuter_groups)
        : Base(Populations({mio::AgeGroup(num_agegroups), CommuterIndex(num_commuter_groups + 1),
                            InfectionStateSecirtsExplicit::Count}),
               ParametersSecirtsExplicit(num_agegroups))
        , m_num_commuter_groups(num_commuter_groups)
    {
    }

    int get_num_commuter_groups() const
    {
        return m_num_commuter_groups;
    }

    /**
    * @brief Calculates smoothed vaccinations for a given time point.
    */
    Eigen::VectorXd
    vaccinations_at(const FP t,
                    const mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::SimulationDay>& daily_vaccinations,
                    const FP eps = 0.15) const
    {
        auto const& params = this->parameters;
        const FP ub        = (size_t)t + 1.0;
        const FP lb        = ub - eps;

        const auto max_time = static_cast<size_t>(daily_vaccinations.template size<mio::SimulationDay>()) - 1;

        Eigen::VectorXd smoothed_vaccinations((size_t)params.get_num_groups());
        smoothed_vaccinations.setZero();

        // if daily_vaccinations is not available for the current time point, we return zero vaccinations.
        if (max_time <= (size_t)t) {
            mio::log_warning("Vaccination data not available for time point ", t, ". Returning zero vaccinations.");
            return smoothed_vaccinations;
        }
        if (t >= lb) {
            for (mio::AgeGroup age = mio::AgeGroup(0); age < params.get_num_groups(); age++) {
                // if ub + 1 is out of bounds, we use the value at ub
                auto ubp1 = static_cast<size_t>(ub + 1);
                if (max_time < ubp1) {
                    ubp1 = static_cast<size_t>(ub);
                }
                const auto num_vaccinations_ub =
                    daily_vaccinations[{age, mio::SimulationDay(static_cast<size_t>(ubp1))}] -
                    daily_vaccinations[{age, mio::SimulationDay(static_cast<size_t>(ub))}];
                const auto num_vaccinations_lb =
                    daily_vaccinations[{age, mio::SimulationDay(static_cast<size_t>(lb + 1))}] -
                    daily_vaccinations[{age, mio::SimulationDay(static_cast<size_t>(lb))}];
                smoothed_vaccinations[(size_t)age] =
                    mio::smoother_cosine(t, lb, ub, num_vaccinations_lb, num_vaccinations_ub);
            }
        }
        else {
            for (auto age = mio::AgeGroup(0); age < params.get_num_groups(); age++) {
                smoothed_vaccinations[(size_t)age] = daily_vaccinations[{age, mio::SimulationDay((size_t)t + 1)}] -
                                                     daily_vaccinations[{age, mio::SimulationDay((size_t)t)}];
            }
        }
        return smoothed_vaccinations;
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        mio::unused(pop);
        flows.setZero();
        auto const& params                            = this->parameters;
        mio::AgeGroup n_agegroups                     = params.get_num_groups();
        mio::ContactMatrixGroup const& contact_matrix = params.template get<mio::osecirts::ContactPatterns<FP>>();

        FP test_and_trace_required = 0.0;
        FP icu_occupancy           = 0.0;

        for (auto i = mio::AgeGroup(0); i < n_agegroups; ++i) {
            for (int c = 0; c <= m_num_commuter_groups; ++c) {
                auto c_type = static_cast<CommuterType>(c);
                using IS    = mio::osecirts::InfectionState;

                test_and_trace_required +=
                    (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                    params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] *
                    (y[this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsNaive})] +
                     y[this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsNaiveConfirmed})]);

                test_and_trace_required +=
                    (params.template get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<FP>>()[i] /
                     params.template get<mio::osecirts::ReducExposedPartialImmunity<FP>>()[i]) *
                    (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] *
                     params.template get<mio::osecirts::ReducTimeInfectedMild<FP>>()[i]) *
                    (y[this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsPartialImmunity})] +
                     y[this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsPartialImmunityConfirmed})]);

                test_and_trace_required +=
                    (params.template get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<FP>>()[i] /
                     params.template get<mio::osecirts::ReducExposedImprovedImmunity<FP>>()[i]) *
                    (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] *
                     params.template get<mio::osecirts::ReducTimeInfectedMild<FP>>()[i]) *
                    (y[this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsImprovedImmunity})] +
                     y[this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsImprovedImmunityConfirmed})]);

                icu_occupancy += y[this->populations.get_flat_index({i, c_type, IS::InfectedCriticalNaive})] +
                                 y[this->populations.get_flat_index({i, c_type, IS::InfectedCriticalPartialImmunity})] +
                                 y[this->populations.get_flat_index({i, c_type, IS::InfectedCriticalImprovedImmunity})];
            }
        }

        auto const partial_vaccination =
            vaccinations_at(t, params.template get<mio::osecirts::DailyPartialVaccinations<FP>>());
        auto const full_vaccination =
            vaccinations_at(t, params.template get<mio::osecirts::DailyFullVaccinations<FP>>());
        auto const booster_vaccination =
            vaccinations_at(t, params.template get<mio::osecirts::DailyBoosterVaccinations<FP>>());

        // calculate force of infection and flows for each age group
        for (auto i = mio::AgeGroup(0); i < n_agegroups; i++) {

            FP criticalPerSevereAdjusted =
                mio::smoother_cosine(icu_occupancy, 0.90 * params.template get<mio::osecirts::ICUCapacity<FP>>(),
                                     params.template get<mio::osecirts::ICUCapacity<FP>>(),
                                     params.template get<mio::osecirts::CriticalPerSevere<FP>>()[i], 0.0);
            FP deathsPerSevereAdjusted =
                params.template get<mio::osecirts::CriticalPerSevere<FP>>()[i] - criticalPerSevereAdjusted;

            auto riskFromInfectedSymptomatic = mio::smoother_cosine(
                test_and_trace_required, params.template get<mio::osecirts::TestAndTraceCapacity<FP>>(),
                params.template get<mio::osecirts::TestAndTraceCapacity<FP>>() *
                    params.template get<mio::osecirts::TestAndTraceCapacityMaxRiskSymptoms<FP>>(),
                params.template get<mio::osecirts::RiskOfInfectionFromSymptomatic<FP>>()[i],
                params.template get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<FP>>()[i]);

            auto riskFromInfectedNoSymptoms = mio::smoother_cosine(
                test_and_trace_required, params.template get<mio::osecirts::TestAndTraceCapacity<FP>>(),
                params.template get<mio::osecirts::TestAndTraceCapacity<FP>>() *
                    params.template get<mio::osecirts::TestAndTraceCapacityMaxRiskNoSymptoms<FP>>(),
                params.template get<mio::osecirts::RelativeTransmissionNoSymptoms<FP>>()[i], 1.0);

            FP ext_inf_force_dummy = 0.0;
            for (auto j = mio::AgeGroup(0); j < n_agegroups; j++) {
                FP Nj_total    = 0.0;
                FP Ij_no_sympt = 0.0;
                FP Ij_sympt    = 0.0;

                using IS = mio::osecirts::InfectionState;
                for (int c2 = 0; c2 <= m_num_commuter_groups; ++c2) {
                    auto c_type_j = static_cast<CommuterType>(c2);
                    for (size_t state = 0; state < static_cast<size_t>(IS::Count); ++state) {
                        if (state != static_cast<size_t>(IS::DeadNaive) &&
                            state != static_cast<size_t>(IS::DeadPartialImmunity) &&
                            state != static_cast<size_t>(IS::DeadImprovedImmunity) &&
                            state != static_cast<size_t>(IS::TemporaryImmunePartialImmunity) &&
                            state != static_cast<size_t>(IS::TemporaryImmuneImprovedImmunity)) {
                            // WICHTIG: Nutze 'y'
                            Nj_total += y[this->populations.get_flat_index({j, c_type_j, static_cast<IS>(state)})];
                        }
                    }
                    Ij_no_sympt +=
                        y[this->populations.get_flat_index({j, c_type_j, IS::InfectedNoSymptomsNaive})] +
                        y[this->populations.get_flat_index({j, c_type_j, IS::InfectedNoSymptomsPartialImmunity})] +
                        y[this->populations.get_flat_index({j, c_type_j, IS::InfectedNoSymptomsImprovedImmunity})];

                    Ij_sympt +=
                        y[this->populations.get_flat_index({j, c_type_j, IS::InfectedSymptomsNaive})] +
                        y[this->populations.get_flat_index({j, c_type_j, IS::InfectedSymptomsPartialImmunity})] +
                        y[this->populations.get_flat_index({j, c_type_j, IS::InfectedSymptomsImprovedImmunity})];
                }

                FP season_val =
                    (1.0 +
                     params.template get<mio::osecirts::Seasonality<FP>>() *
                         sin(3.141592653589793 *
                             (std::fmod((params.template get<mio::osecirts::StartDay>() + t), 365.0) / 182.5 + 0.5)));
                FP cont_freq_eff = season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                                static_cast<Eigen::Index>((size_t)j));
                const FP divNj   = (Nj_total < 1e-12) ? 0.0 : 1.0 / Nj_total;

                ext_inf_force_dummy +=
                    cont_freq_eff * divNj *
                    params.template get<mio::osecirts::TransmissionProbabilityOnContact<FP>>()[i] *
                    (riskFromInfectedNoSymptoms * Ij_no_sympt + riskFromInfectedSymptomatic * Ij_sympt);
            }

            // flows for each commuter group
            for (int c = 0; c <= m_num_commuter_groups; ++c) {
                auto c_type = static_cast<CommuterType>(c);
                using IS    = mio::osecirts::InfectionState;

                size_t SNi    = this->populations.get_flat_index({i, c_type, IS::SusceptibleNaive});
                size_t ENi    = this->populations.get_flat_index({i, c_type, IS::ExposedNaive});
                size_t INSNi  = this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsNaive});
                size_t ISyNi  = this->populations.get_flat_index({i, c_type, IS::InfectedSymptomsNaive});
                size_t ISevNi = this->populations.get_flat_index({i, c_type, IS::InfectedSevereNaive});
                size_t ICrNi  = this->populations.get_flat_index({i, c_type, IS::InfectedCriticalNaive});
                size_t INSNCi = this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsNaiveConfirmed});
                size_t ISyNCi = this->populations.get_flat_index({i, c_type, IS::InfectedSymptomsNaiveConfirmed});

                size_t SPIi    = this->populations.get_flat_index({i, c_type, IS::SusceptiblePartialImmunity});
                size_t EPIi    = this->populations.get_flat_index({i, c_type, IS::ExposedPartialImmunity});
                size_t INSPIi  = this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsPartialImmunity});
                size_t ISyPIi  = this->populations.get_flat_index({i, c_type, IS::InfectedSymptomsPartialImmunity});
                size_t ISevPIi = this->populations.get_flat_index({i, c_type, IS::InfectedSeverePartialImmunity});
                size_t ICrPIi  = this->populations.get_flat_index({i, c_type, IS::InfectedCriticalPartialImmunity});
                size_t INSPICi =
                    this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsPartialImmunityConfirmed});
                size_t ISyPICi =
                    this->populations.get_flat_index({i, c_type, IS::InfectedSymptomsPartialImmunityConfirmed});

                size_t EIIi    = this->populations.get_flat_index({i, c_type, IS::ExposedImprovedImmunity});
                size_t INSIIi  = this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsImprovedImmunity});
                size_t ISyIIi  = this->populations.get_flat_index({i, c_type, IS::InfectedSymptomsImprovedImmunity});
                size_t ISevIIi = this->populations.get_flat_index({i, c_type, IS::InfectedSevereImprovedImmunity});
                size_t ICrIIi  = this->populations.get_flat_index({i, c_type, IS::InfectedCriticalImprovedImmunity});
                size_t INSIICi =
                    this->populations.get_flat_index({i, c_type, IS::InfectedNoSymptomsImprovedImmunityConfirmed});
                size_t ISyIICi =
                    this->populations.get_flat_index({i, c_type, IS::InfectedSymptomsImprovedImmunityConfirmed});

                size_t TImm1 = this->populations.get_flat_index({i, c_type, IS::TemporaryImmunePartialImmunity});
                size_t TImm2 = this->populations.get_flat_index({i, c_type, IS::TemporaryImmuneImprovedImmunity});
                size_t SIIi  = this->populations.get_flat_index({i, c_type, IS::SusceptibleImprovedImmunity});

                FP reducEPI   = params.template get<mio::osecirts::ReducExposedPartialImmunity<FP>>()[i];
                FP reducEII   = params.template get<mio::osecirts::ReducExposedImprovedImmunity<FP>>()[i];
                FP reducISyPI = params.template get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<FP>>()[i];
                FP reducISyII = params.template get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<FP>>()[i];
                FP reducSevPI =
                    params.template get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i];
                FP reducSevII =
                    params.template get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i];
                FP reducTimeM = params.template get<mio::osecirts::ReducTimeInfectedMild<FP>>()[i];

                FP dummy_SN  = y[SNi] * ext_inf_force_dummy;
                FP dummy_SPI = y[SPIi] * reducEPI * ext_inf_force_dummy;
                FP dummy_SII = y[SIIi] * reducEII * ext_inf_force_dummy;

                flows[this->template get_flat_flow_index<IS::SusceptibleNaive, IS::ExposedNaive>({i, c_type})] =
                    dummy_SN;
                flows[this->template get_flat_flow_index<IS::SusceptiblePartialImmunity, IS::ExposedPartialImmunity>(
                    {i, c_type})] = dummy_SPI;
                flows[this->template get_flat_flow_index<IS::SusceptibleImprovedImmunity, IS::ExposedImprovedImmunity>(
                    {i, c_type})] = dummy_SII;

                FP total_SN = 0.0, total_SPI = 0.0, total_SII = 0.0;
                for (int c_all = 0; c_all <= m_num_commuter_groups; ++c_all) {
                    auto ct_all = static_cast<CommuterType>(c_all);
                    total_SN += y[this->populations.get_flat_index({i, ct_all, IS::SusceptibleNaive})];
                    total_SPI += y[this->populations.get_flat_index({i, ct_all, IS::SusceptiblePartialImmunity})];
                    total_SII += y[this->populations.get_flat_index({i, ct_all, IS::SusceptibleImprovedImmunity})];
                }

                FP frac_SN  = (total_SN > 1e-12) ? (y[SNi] / total_SN) : 0.0;
                FP frac_SPI = (total_SPI > 1e-12) ? (y[SPIi] / total_SPI) : 0.0;
                FP frac_SII = (total_SII > 1e-12) ? (y[SIIi] / total_SII) : 0.0;

                flows[this->template get_flat_flow_index<IS::SusceptibleNaive, IS::TemporaryImmunePartialImmunity>(
                    {i, c_type})] =
                    std::min(std::max(0.0, y[SNi] - dummy_SN), partial_vaccination[static_cast<size_t>(i)] * frac_SN);

                flows[this->template get_flat_flow_index<IS::SusceptiblePartialImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    std::min(std::max(0.0, y[SPIi] - dummy_SPI), full_vaccination[static_cast<size_t>(i)] * frac_SPI);

                flows[this->template get_flat_flow_index<IS::SusceptibleImprovedImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    std::min(std::max(0.0, y[SIIi] - dummy_SII),
                             booster_vaccination[static_cast<size_t>(i)] * frac_SII);

                flows[this->template get_flat_flow_index<IS::ExposedNaive, IS::InfectedNoSymptomsNaive>({i, c_type})] =
                    y[ENi] / params.template get<mio::osecirts::TimeExposed<FP>>()[i];

                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsNaive,
                                                         IS::TemporaryImmunePartialImmunity>({i, c_type})] =
                    params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i] *
                    (1.0 / params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i]) * y[INSNi];
                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsNaive, IS::InfectedSymptomsNaive>(
                    {i, c_type})] =
                    (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                    params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * y[INSNi];

                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsNaiveConfirmed,
                                                         IS::InfectedSymptomsNaiveConfirmed>({i, c_type})] =
                    (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                    params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * y[INSNCi];
                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsNaiveConfirmed,
                                                         IS::TemporaryImmunePartialImmunity>({i, c_type})] =
                    params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i] *
                    (1.0 / params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i]) * y[INSNCi];

                flows[this->template get_flat_flow_index<IS::InfectedSymptomsNaive, IS::InfectedSevereNaive>(
                    {i, c_type})] = params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i] /
                                    params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * y[ISyNi];
                flows[this->template get_flat_flow_index<IS::InfectedSymptomsNaive, IS::TemporaryImmunePartialImmunity>(
                    {i, c_type})] = (1.0 - params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i]) /
                                    params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * y[ISyNi];

                flows[this->template get_flat_flow_index<IS::InfectedSymptomsNaiveConfirmed, IS::InfectedSevereNaive>(
                    {i, c_type})] = params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i] /
                                    params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * y[ISyNCi];
                flows[this->template get_flat_flow_index<IS::InfectedSymptomsNaiveConfirmed,
                                                         IS::TemporaryImmunePartialImmunity>({i, c_type})] =
                    (1.0 - params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i]) /
                    params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * y[ISyNCi];

                flows[this->template get_flat_flow_index<IS::InfectedSevereNaive, IS::InfectedCriticalNaive>(
                    {i, c_type})] = criticalPerSevereAdjusted /
                                    params.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i] * y[ISevNi];
                flows[this->template get_flat_flow_index<IS::InfectedSevereNaive, IS::TemporaryImmunePartialImmunity>(
                    {i, c_type})] = (1.0 - params.template get<mio::osecirts::CriticalPerSevere<FP>>()[i]) /
                                    params.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i] * y[ISevNi];
                flows[this->template get_flat_flow_index<IS::InfectedSevereNaive, IS::DeadNaive>({i, c_type})] =
                    deathsPerSevereAdjusted / params.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i] *
                    y[ISevNi];

                flows[this->template get_flat_flow_index<IS::InfectedCriticalNaive, IS::DeadNaive>({i, c_type})] =
                    params.template get<mio::osecirts::DeathsPerCritical<FP>>()[i] /
                    params.template get<mio::osecirts::TimeInfectedCritical<FP>>()[i] * y[ICrNi];
                flows[this->template get_flat_flow_index<IS::InfectedCriticalNaive, IS::TemporaryImmunePartialImmunity>(
                    {i, c_type})] = (1.0 - params.template get<mio::osecirts::DeathsPerCritical<FP>>()[i]) /
                                    params.template get<mio::osecirts::TimeInfectedCritical<FP>>()[i] * y[ICrNi];

                flows[this->template get_flat_flow_index<IS::SusceptiblePartialImmunity, IS::SusceptibleNaive>(
                    {i, c_type})] =
                    1.0 / params.template get<mio::osecirts::TimeWaningPartialImmunity<FP>>()[i] * y[SPIi];
                flows[this->template get_flat_flow_index<IS::TemporaryImmunePartialImmunity,
                                                         IS::SusceptiblePartialImmunity>({i, c_type})] =
                    1.0 / params.template get<mio::osecirts::TimeTemporaryImmunityPI<FP>>()[i] * y[TImm1];

                flows[this->template get_flat_flow_index<IS::ExposedPartialImmunity,
                                                         IS::InfectedNoSymptomsPartialImmunity>({i, c_type})] =
                    y[EPIi] / params.template get<mio::osecirts::TimeExposed<FP>>()[i];

                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsPartialImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducISyPI / reducEPI) *
                               (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i])) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * reducTimeM) * y[INSPIi];
                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsPartialImmunity,
                                                         IS::InfectedSymptomsPartialImmunity>({i, c_type})] =
                    (reducISyPI / reducEPI) *
                    (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * reducTimeM) * y[INSPIi];

                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsPartialImmunityConfirmed,
                                                         IS::InfectedSymptomsPartialImmunityConfirmed>({i, c_type})] =
                    (reducISyPI / reducEPI) *
                    (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * reducTimeM) * y[INSPICi];
                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsPartialImmunityConfirmed,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducISyPI / reducEPI) *
                               (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i])) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * reducTimeM) * y[INSPICi];

                flows[this->template get_flat_flow_index<IS::InfectedSymptomsPartialImmunity,
                                                         IS::InfectedSeverePartialImmunity>({i, c_type})] =
                    reducSevPI / reducISyPI * params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i] /
                    (params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * reducTimeM) * y[ISyPIi];
                flows[this->template get_flat_flow_index<IS::InfectedSymptomsPartialImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducSevPI / reducISyPI) *
                               params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * reducTimeM) * y[ISyPIi];

                flows[this->template get_flat_flow_index<IS::InfectedSymptomsPartialImmunityConfirmed,
                                                         IS::InfectedSeverePartialImmunity>({i, c_type})] =
                    reducSevPI / reducISyPI * params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i] /
                    (params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * reducTimeM) * y[ISyPICi];
                flows[this->template get_flat_flow_index<IS::InfectedSymptomsPartialImmunityConfirmed,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducSevPI / reducISyPI) *
                               params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * reducTimeM) * y[ISyPICi];

                flows[this->template get_flat_flow_index<IS::InfectedSeverePartialImmunity,
                                                         IS::InfectedCriticalPartialImmunity>({i, c_type})] =
                    (reducSevPI / reducSevPI) * criticalPerSevereAdjusted /
                    params.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i] * y[ISevPIi];
                flows[this->template get_flat_flow_index<IS::InfectedSeverePartialImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducSevPI / reducSevPI) * params.template get<mio::osecirts::CriticalPerSevere<FP>>()[i]) /
                    params.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i] * y[ISevPIi];
                flows[this->template get_flat_flow_index<IS::InfectedSeverePartialImmunity, IS::DeadPartialImmunity>(
                    {i, c_type})] = (reducSevPI / reducSevPI) * deathsPerSevereAdjusted /
                                    params.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i] * y[ISevPIi];

                flows[this->template get_flat_flow_index<IS::InfectedCriticalPartialImmunity, IS::DeadPartialImmunity>(
                    {i, c_type})] = (reducSevPI / reducSevPI) *
                                    params.template get<mio::osecirts::DeathsPerCritical<FP>>()[i] /
                                    params.template get<mio::osecirts::TimeInfectedCritical<FP>>()[i] * y[ICrPIi];
                flows[this->template get_flat_flow_index<IS::InfectedCriticalPartialImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducSevPI / reducSevPI) * params.template get<mio::osecirts::DeathsPerCritical<FP>>()[i]) /
                    params.template get<mio::osecirts::TimeInfectedCritical<FP>>()[i] * y[ICrPIi];

                flows[this->template get_flat_flow_index<IS::ExposedImprovedImmunity,
                                                         IS::InfectedNoSymptomsImprovedImmunity>({i, c_type})] =
                    y[EIIi] / params.template get<mio::osecirts::TimeExposed<FP>>()[i];

                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsImprovedImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducISyII / reducEII) *
                               (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i])) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * reducTimeM) * y[INSIIi];
                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsImprovedImmunity,
                                                         IS::InfectedSymptomsImprovedImmunity>({i, c_type})] =
                    (reducISyII / reducEII) *
                    (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * reducTimeM) * y[INSIIi];

                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsImprovedImmunityConfirmed,
                                                         IS::InfectedSymptomsImprovedImmunityConfirmed>({i, c_type})] =
                    (reducISyII / reducEII) *
                    (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * reducTimeM) * y[INSIICi];
                flows[this->template get_flat_flow_index<IS::InfectedNoSymptomsImprovedImmunityConfirmed,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducISyII / reducEII) *
                               (1.0 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<FP>>()[i])) /
                    (params.template get<mio::osecirts::TimeInfectedNoSymptoms<FP>>()[i] * reducTimeM) * y[INSIICi];

                flows[this->template get_flat_flow_index<IS::InfectedSymptomsImprovedImmunity,
                                                         IS::InfectedSevereImprovedImmunity>({i, c_type})] =
                    reducSevII / reducISyII * params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i] /
                    (params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * reducTimeM) * y[ISyIIi];
                flows[this->template get_flat_flow_index<IS::InfectedSymptomsImprovedImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducSevII / reducISyII) *
                               params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * reducTimeM) * y[ISyIIi];

                flows[this->template get_flat_flow_index<IS::InfectedSymptomsImprovedImmunityConfirmed,
                                                         IS::InfectedSevereImprovedImmunity>({i, c_type})] =
                    reducSevII / reducISyII * params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i] /
                    (params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * reducTimeM) * y[ISyIICi];
                flows[this->template get_flat_flow_index<IS::InfectedSymptomsImprovedImmunityConfirmed,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducSevII / reducISyII) *
                               params.template get<mio::osecirts::SeverePerInfectedSymptoms<FP>>()[i]) /
                    (params.template get<mio::osecirts::TimeInfectedSymptoms<FP>>()[i] * reducTimeM) * y[ISyIICi];

                flows[this->template get_flat_flow_index<IS::InfectedSevereImprovedImmunity,
                                                         IS::InfectedCriticalImprovedImmunity>({i, c_type})] =
                    (reducSevII / reducSevII) * criticalPerSevereAdjusted /
                    params.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i] * y[ISevIIi];
                flows[this->template get_flat_flow_index<IS::InfectedSevereImprovedImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducSevII / reducSevII) * params.template get<mio::osecirts::CriticalPerSevere<FP>>()[i]) /
                    params.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i] * y[ISevIIi];
                flows[this->template get_flat_flow_index<IS::InfectedSevereImprovedImmunity, IS::DeadImprovedImmunity>(
                    {i, c_type})] = (reducSevII / reducSevII) * deathsPerSevereAdjusted /
                                    params.template get<mio::osecirts::TimeInfectedSevere<FP>>()[i] * y[ISevIIi];

                flows[this->template get_flat_flow_index<IS::InfectedCriticalImprovedImmunity,
                                                         IS::DeadImprovedImmunity>({i, c_type})] =
                    (reducSevII / reducSevII) * params.template get<mio::osecirts::DeathsPerCritical<FP>>()[i] /
                    params.template get<mio::osecirts::TimeInfectedCritical<FP>>()[i] * y[ICrIIi];
                flows[this->template get_flat_flow_index<IS::InfectedCriticalImprovedImmunity,
                                                         IS::TemporaryImmuneImprovedImmunity>({i, c_type})] =
                    (1.0 - (reducSevII / reducSevII) * params.template get<mio::osecirts::DeathsPerCritical<FP>>()[i]) /
                    params.template get<mio::osecirts::TimeInfectedCritical<FP>>()[i] * y[ICrIIi];

                flows[this->template get_flat_flow_index<IS::TemporaryImmuneImprovedImmunity,
                                                         IS::SusceptibleImprovedImmunity>({i, c_type})] =
                    1.0 / params.template get<mio::osecirts::TimeTemporaryImmunityII<FP>>()[i] * y[TImm2];
                flows[this->template get_flat_flow_index<IS::SusceptibleImprovedImmunity,
                                                         IS::SusceptiblePartialImmunity>({i, c_type})] =
                    1.0 / params.template get<mio::osecirts::TimeWaningImprovedImmunity<FP>>()[i] * y[SIIi];
            }
        }
    }
};

using StandardModelLagrangianSecirts = ModelSecirtsExplicit<ScalarType>;
using StandardLagrangianSimSecirts   = mio::FlowSimulation<ScalarType, StandardModelLagrangianSecirts>;

using StandardModelSecirts = mio::osecirts::Model<ScalarType>;
using StandardSimSecirts   = mio::FlowSimulation<ScalarType, StandardModelSecirts>;

struct SecirtsIntensities {
    std::vector<double> foi;
    std::vector<double> crit_per_sev;
    std::vector<double> dead_per_sev;
    std::vector<double> vac_rate_naive;
    std::vector<double> vac_rate_pi;
    std::vector<double> vac_rate_ii;

    void resize(size_t NG)
    {
        foi.assign(NG, 0.0);
        crit_per_sev.assign(NG, 0.0);
        dead_per_sev.assign(NG, 0.0);
        vac_rate_naive.assign(NG, 0.0);
        vac_rate_pi.assign(NG, 0.0);
        vac_rate_ii.assign(NG, 0.0);
    }
};

namespace
{

struct ModelEvaluatorSecirts {
    const mio::osecirts::Model<double>& model;
    std::size_t NG;
    std::size_t NC;

    ModelEvaluatorSecirts(const mio::osecirts::Model<double>& m)
        : model(m)
    {
        using IS = mio::osecirts::InfectionState;
        NG       = static_cast<std::size_t>(model.parameters.get_num_groups());
        NC       = static_cast<std::size_t>(IS::Count) * NG; // 29 * NG
    }

    Eigen::VectorXd
    vaccinations_at(const double t,
                    const mio::CustomIndexArray<double, mio::AgeGroup, mio::SimulationDay>& daily_vaccinations,
                    const double eps = 0.15) const
    {
        const double ub     = (size_t)t + 1.0;
        const double lb     = ub - eps;
        const auto max_time = static_cast<size_t>(daily_vaccinations.template size<mio::SimulationDay>()) - 1;

        Eigen::VectorXd smoothed_vaccinations((size_t)NG);
        smoothed_vaccinations.setZero();

        if (max_time <= (size_t)t)
            return smoothed_vaccinations;

        if (t >= lb) {
            for (mio::AgeGroup age = mio::AgeGroup(0); age < mio::AgeGroup(NG); age++) {
                auto ubp1 = static_cast<size_t>(ub + 1);
                if (max_time < ubp1)
                    ubp1 = static_cast<size_t>(ub);
                const auto num_vaccinations_ub = daily_vaccinations[{age, mio::SimulationDay(ubp1)}] -
                                                 daily_vaccinations[{age, mio::SimulationDay(static_cast<size_t>(ub))}];
                const auto num_vaccinations_lb =
                    daily_vaccinations[{age, mio::SimulationDay(static_cast<size_t>(lb + 1))}] -
                    daily_vaccinations[{age, mio::SimulationDay(static_cast<size_t>(lb))}];
                smoothed_vaccinations[(size_t)age] =
                    mio::smoother_cosine(t, lb, ub, num_vaccinations_lb, num_vaccinations_ub);
            }
        }
        else {
            for (auto age = mio::AgeGroup(0); age < mio::AgeGroup(NG); age++) {
                smoothed_vaccinations[(size_t)age] = daily_vaccinations[{age, mio::SimulationDay((size_t)t + 1)}] -
                                                     daily_vaccinations[{age, mio::SimulationDay((size_t)t)}];
            }
        }
        return smoothed_vaccinations;
    }

    // Total rhs
    void get_totals_rhs(const Eigen::VectorXd& y, double t, Eigen::VectorXd& dy) const
    {
        dy.resize(NC);
        model.get_derivatives(y, y, t, dy);
    }

    // intensities
    void get_intensities(const Eigen::VectorXd& y_tot, double t, SecirtsIntensities& intensities) const
    {
        intensities.resize(NG);
        auto const& params                            = model.parameters;
        mio::ContactMatrixGroup const& contact_matrix = params.template get<mio::osecirts::ContactPatterns<double>>();
        using IS                                      = mio::osecirts::InfectionState;

        double test_and_trace_required = 0.0;
        double icu_occupancy           = 0.0;

        for (size_t i = 0; i < NG; ++i) {
            auto ag        = mio::AgeGroup(i);
            size_t INSNi   = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsNaive});
            size_t INSNCi  = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsNaiveConfirmed});
            size_t INSPIi  = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsPartialImmunity});
            size_t INSPICi = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsPartialImmunityConfirmed});
            size_t INSIIi  = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsImprovedImmunity});
            size_t INSIICi = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsImprovedImmunityConfirmed});

            test_and_trace_required +=
                (1 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>()[ag]) /
                params.template get<mio::osecirts::TimeInfectedNoSymptoms<double>>()[ag] *
                (y_tot[INSNi] + y_tot[INSNCi]);

            test_and_trace_required +=
                (params.template get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>()[ag] /
                 params.template get<mio::osecirts::ReducExposedPartialImmunity<double>>()[ag]) *
                (1 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>()[ag]) /
                (params.template get<mio::osecirts::TimeInfectedNoSymptoms<double>>()[ag] *
                 params.template get<mio::osecirts::ReducTimeInfectedMild<double>>()[ag]) *
                (y_tot[INSPIi] + y_tot[INSPICi]);

            test_and_trace_required +=
                (params.template get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>()[ag] /
                 params.template get<mio::osecirts::ReducExposedImprovedImmunity<double>>()[ag]) *
                (1 - params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>()[ag]) /
                (params.template get<mio::osecirts::TimeInfectedNoSymptoms<double>>()[ag] *
                 params.template get<mio::osecirts::ReducTimeInfectedMild<double>>()[ag]) *
                (y_tot[INSIIi] + y_tot[INSIICi]);

            icu_occupancy += y_tot[model.populations.get_flat_index({ag, IS::InfectedCriticalNaive})] +
                             y_tot[model.populations.get_flat_index({ag, IS::InfectedCriticalPartialImmunity})] +
                             y_tot[model.populations.get_flat_index({ag, IS::InfectedCriticalImprovedImmunity})];
        }

        auto partial_vac = vaccinations_at(t, params.template get<mio::osecirts::DailyPartialVaccinations<double>>());
        auto full_vac    = vaccinations_at(t, params.template get<mio::osecirts::DailyFullVaccinations<double>>());
        auto booster_vac = vaccinations_at(t, params.template get<mio::osecirts::DailyBoosterVaccinations<double>>());

        for (size_t i = 0; i < NG; ++i) {
            auto ag = mio::AgeGroup(i);

            intensities.crit_per_sev[i] =
                mio::smoother_cosine(icu_occupancy, 0.90 * params.template get<mio::osecirts::ICUCapacity<double>>(),
                                     params.template get<mio::osecirts::ICUCapacity<double>>(),
                                     params.template get<mio::osecirts::CriticalPerSevere<double>>()[ag], 0.0);
            intensities.dead_per_sev[i] =
                params.template get<mio::osecirts::CriticalPerSevere<double>>()[ag] - intensities.crit_per_sev[i];

            auto riskFromInfectedSymptomatic = mio::smoother_cosine(
                test_and_trace_required, params.template get<mio::osecirts::TestAndTraceCapacity<double>>(),
                params.template get<mio::osecirts::TestAndTraceCapacity<double>>() *
                    params.template get<mio::osecirts::TestAndTraceCapacityMaxRiskSymptoms<double>>(),
                params.template get<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>()[ag],
                params.template get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<double>>()[ag]);

            auto riskFromInfectedNoSymptoms = mio::smoother_cosine(
                test_and_trace_required, params.template get<mio::osecirts::TestAndTraceCapacity<double>>(),
                params.template get<mio::osecirts::TestAndTraceCapacity<double>>() *
                    params.template get<mio::osecirts::TestAndTraceCapacityMaxRiskNoSymptoms<double>>(),
                params.template get<mio::osecirts::RelativeTransmissionNoSymptoms<double>>()[ag], 1.0);

            double lambda = 0.0;
            for (size_t j = 0; j < NG; ++j) {
                auto ag_j = mio::AgeGroup(j);
                double Nj = 0.0;
                for (size_t state = 0; state < static_cast<size_t>(IS::Count); ++state) {
                    if (state != static_cast<size_t>(IS::DeadNaive) &&
                        state != static_cast<size_t>(IS::DeadPartialImmunity) &&
                        state != static_cast<size_t>(IS::DeadImprovedImmunity) &&
                        state != static_cast<size_t>(IS::TemporaryImmunePartialImmunity) &&
                        state != static_cast<size_t>(IS::TemporaryImmuneImprovedImmunity)) {
                        Nj += y_tot[model.populations.get_flat_index({ag_j, static_cast<IS>(state)})];
                    }
                }

                double Ij_no_sympt =
                    y_tot[model.populations.get_flat_index({ag_j, IS::InfectedNoSymptomsNaive})] +
                    y_tot[model.populations.get_flat_index({ag_j, IS::InfectedNoSymptomsPartialImmunity})] +
                    y_tot[model.populations.get_flat_index({ag_j, IS::InfectedNoSymptomsImprovedImmunity})];

                double Ij_sympt = y_tot[model.populations.get_flat_index({ag_j, IS::InfectedSymptomsNaive})] +
                                  y_tot[model.populations.get_flat_index({ag_j, IS::InfectedSymptomsPartialImmunity})] +
                                  y_tot[model.populations.get_flat_index({ag_j, IS::InfectedSymptomsImprovedImmunity})];

                double season_val =
                    (1 +
                     params.template get<mio::osecirts::Seasonality<double>>() *
                         sin(3.141592653589793 *
                             (std::fmod((params.template get<mio::osecirts::StartDay>() + t), 365.0) / 182.5 + 0.5)));
                double cont_freq = season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>(i),
                                                                                static_cast<Eigen::Index>(j));
                double divNj     = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;

                lambda += cont_freq * divNj *
                          params.template get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[ag] *
                          (riskFromInfectedNoSymptoms * Ij_no_sympt + riskFromInfectedSymptomatic * Ij_sympt);
            }
            intensities.foi[i] = lambda;

            size_t SNi  = model.populations.get_flat_index({ag, IS::SusceptibleNaive});
            size_t SPIi = model.populations.get_flat_index({ag, IS::SusceptiblePartialImmunity});
            size_t SIIi = model.populations.get_flat_index({ag, IS::SusceptibleImprovedImmunity});

            double dummy_SN = y_tot[SNi] * lambda;
            double dummy_SPI =
                y_tot[SPIi] * params.template get<mio::osecirts::ReducExposedPartialImmunity<double>>()[ag] * lambda;
            double dummy_SII =
                y_tot[SIIi] * params.template get<mio::osecirts::ReducExposedImprovedImmunity<double>>()[ag] * lambda;

            double vac_N  = std::min(std::max(0.0, y_tot[SNi] - dummy_SN), partial_vac[i]);
            double vac_PI = std::min(std::max(0.0, y_tot[SPIi] - dummy_SPI), full_vac[i]);
            double vac_II = std::min(std::max(0.0, y_tot[SIIi] - dummy_SII), booster_vac[i]);

            intensities.vac_rate_naive[i] = (y_tot[SNi] > 1e-12) ? (vac_N / y_tot[SNi]) : 0.0;
            intensities.vac_rate_pi[i]    = (y_tot[SPIi] > 1e-12) ? (vac_PI / y_tot[SPIi]) : 0.0;
            intensities.vac_rate_ii[i]    = (y_tot[SIIi] > 1e-12) ? (vac_II / y_tot[SIIi]) : 0.0;
        }
    }

    // commuter reconstruction
    void get_commuter_rhs(const SecirtsIntensities& intensities, const Eigen::VectorXd& Xc, Eigen::VectorXd& dxc) const
    {
        dxc.setZero(NC);
        using IS           = mio::osecirts::InfectionState;
        auto const& params = model.parameters;

        auto apply_flow = [&](size_t src, size_t tgt, double val) {
            dxc[src] -= val;
            dxc[tgt] += val;
        };

        for (size_t i = 0; i < NG; ++i) {
            auto ag = mio::AgeGroup(i);

            size_t SNi    = model.populations.get_flat_index({ag, IS::SusceptibleNaive});
            size_t ENi    = model.populations.get_flat_index({ag, IS::ExposedNaive});
            size_t INSNi  = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsNaive});
            size_t ISyNi  = model.populations.get_flat_index({ag, IS::InfectedSymptomsNaive});
            size_t ISevNi = model.populations.get_flat_index({ag, IS::InfectedSevereNaive});
            size_t ICrNi  = model.populations.get_flat_index({ag, IS::InfectedCriticalNaive});
            size_t INSNCi = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsNaiveConfirmed});
            size_t ISyNCi = model.populations.get_flat_index({ag, IS::InfectedSymptomsNaiveConfirmed});
            size_t DNi    = model.populations.get_flat_index({ag, IS::DeadNaive});

            size_t SPIi    = model.populations.get_flat_index({ag, IS::SusceptiblePartialImmunity});
            size_t EPIi    = model.populations.get_flat_index({ag, IS::ExposedPartialImmunity});
            size_t INSPIi  = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsPartialImmunity});
            size_t ISyPIi  = model.populations.get_flat_index({ag, IS::InfectedSymptomsPartialImmunity});
            size_t ISevPIi = model.populations.get_flat_index({ag, IS::InfectedSeverePartialImmunity});
            size_t ICrPIi  = model.populations.get_flat_index({ag, IS::InfectedCriticalPartialImmunity});
            size_t INSPICi = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsPartialImmunityConfirmed});
            size_t ISyPICi = model.populations.get_flat_index({ag, IS::InfectedSymptomsPartialImmunityConfirmed});
            size_t DPIi    = model.populations.get_flat_index({ag, IS::DeadPartialImmunity});

            size_t EIIi    = model.populations.get_flat_index({ag, IS::ExposedImprovedImmunity});
            size_t INSIIi  = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsImprovedImmunity});
            size_t ISyIIi  = model.populations.get_flat_index({ag, IS::InfectedSymptomsImprovedImmunity});
            size_t ISevIIi = model.populations.get_flat_index({ag, IS::InfectedSevereImprovedImmunity});
            size_t ICrIIi  = model.populations.get_flat_index({ag, IS::InfectedCriticalImprovedImmunity});
            size_t INSIICi = model.populations.get_flat_index({ag, IS::InfectedNoSymptomsImprovedImmunityConfirmed});
            size_t ISyIICi = model.populations.get_flat_index({ag, IS::InfectedSymptomsImprovedImmunityConfirmed});
            size_t DIIi    = model.populations.get_flat_index({ag, IS::DeadImprovedImmunity});

            size_t TImm1 = model.populations.get_flat_index({ag, IS::TemporaryImmunePartialImmunity});
            size_t TImm2 = model.populations.get_flat_index({ag, IS::TemporaryImmuneImprovedImmunity});
            size_t SIIi  = model.populations.get_flat_index({ag, IS::SusceptibleImprovedImmunity});

            double rEPI  = params.template get<mio::osecirts::ReducExposedPartialImmunity<double>>()[ag];
            double rEII  = params.template get<mio::osecirts::ReducExposedImprovedImmunity<double>>()[ag];
            double rISPI = params.template get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>()[ag];
            double rISII = params.template get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>()[ag];
            double rSevPI =
                params.template get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[ag];
            double rSevII =
                params.template get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[ag];
            double rTimeM = params.template get<mio::osecirts::ReducTimeInfectedMild<double>>()[ag];

            double timeE    = params.template get<mio::osecirts::TimeExposed<double>>()[ag];
            double timeINS  = params.template get<mio::osecirts::TimeInfectedNoSymptoms<double>>()[ag];
            double timeISy  = params.template get<mio::osecirts::TimeInfectedSymptoms<double>>()[ag];
            double timeISev = params.template get<mio::osecirts::TimeInfectedSevere<double>>()[ag];
            double timeICr  = params.template get<mio::osecirts::TimeInfectedCritical<double>>()[ag];
            double timeWPI  = params.template get<mio::osecirts::TimeWaningPartialImmunity<double>>()[ag];
            double timeWII  = params.template get<mio::osecirts::TimeWaningImprovedImmunity<double>>()[ag];
            double timeTPI  = params.template get<mio::osecirts::TimeTemporaryImmunityPI<double>>()[ag];
            double timeTII  = params.template get<mio::osecirts::TimeTemporaryImmunityII<double>>()[ag];

            double recINS  = params.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>()[ag];
            double sevISy  = params.template get<mio::osecirts::SeverePerInfectedSymptoms<double>>()[ag];
            double deadICr = params.template get<mio::osecirts::DeathsPerCritical<double>>()[ag];
            double critSev = params.template get<mio::osecirts::CriticalPerSevere<double>>()[ag];

            apply_flow(SNi, ENi, intensities.foi[i] * Xc[SNi]);
            apply_flow(SPIi, EPIi, intensities.foi[i] * rEPI * Xc[SPIi]);
            apply_flow(SIIi, EIIi, intensities.foi[i] * rEII * Xc[SIIi]);

            apply_flow(SNi, TImm1, intensities.vac_rate_naive[i] * Xc[SNi]);
            apply_flow(SPIi, TImm2, intensities.vac_rate_pi[i] * Xc[SPIi]);
            apply_flow(SIIi, TImm2, intensities.vac_rate_ii[i] * Xc[SIIi]);

            apply_flow(ENi, INSNi, Xc[ENi] / timeE);
            apply_flow(INSNi, TImm1, recINS / timeINS * Xc[INSNi]);
            apply_flow(INSNi, ISyNi, (1 - recINS) / timeINS * Xc[INSNi]);
            apply_flow(INSNCi, ISyNCi, (1 - recINS) / timeINS * Xc[INSNCi]);
            apply_flow(INSNCi, TImm1, recINS / timeINS * Xc[INSNCi]);

            apply_flow(ISyNi, ISevNi, sevISy / timeISy * Xc[ISyNi]);
            apply_flow(ISyNi, TImm1, (1 - sevISy) / timeISy * Xc[ISyNi]);
            apply_flow(ISyNCi, ISevNi, sevISy / timeISy * Xc[ISyNCi]);
            apply_flow(ISyNCi, TImm1, (1 - sevISy) / timeISy * Xc[ISyNCi]);

            apply_flow(ISevNi, ICrNi, intensities.crit_per_sev[i] / timeISev * Xc[ISevNi]);
            apply_flow(ISevNi, TImm1, (1 - critSev) / timeISev * Xc[ISevNi]);
            apply_flow(ISevNi, DNi, intensities.dead_per_sev[i] / timeISev * Xc[ISevNi]);

            apply_flow(ICrNi, DNi, deadICr / timeICr * Xc[ICrNi]);
            apply_flow(ICrNi, TImm1, (1 - deadICr) / timeICr * Xc[ICrNi]);

            apply_flow(SPIi, SNi, Xc[SPIi] / timeWPI);
            apply_flow(TImm1, SPIi, Xc[TImm1] / timeTPI);

            apply_flow(EPIi, INSPIi, Xc[EPIi] / timeE);
            double pi_TImm = (1 - (rISPI / rEPI) * (1 - recINS)) / (timeINS * rTimeM);
            double pi_ISy  = (rISPI / rEPI) * (1 - recINS) / (timeINS * rTimeM);
            apply_flow(INSPIi, TImm2, pi_TImm * Xc[INSPIi]);
            apply_flow(INSPIi, ISyPIi, pi_ISy * Xc[INSPIi]);
            apply_flow(INSPICi, ISyPICi, pi_ISy * Xc[INSPICi]);
            apply_flow(INSPICi, TImm2, pi_TImm * Xc[INSPICi]);

            double pi_ISev     = (rSevPI / rISPI) * sevISy / (timeISy * rTimeM);
            double pi_ISy_TImm = (1 - (rSevPI / rISPI) * sevISy) / (timeISy * rTimeM);
            apply_flow(ISyPIi, ISevPIi, pi_ISev * Xc[ISyPIi]);
            apply_flow(ISyPIi, TImm2, pi_ISy_TImm * Xc[ISyPIi]);
            apply_flow(ISyPICi, ISevPIi, pi_ISev * Xc[ISyPICi]);
            apply_flow(ISyPICi, TImm2, pi_ISy_TImm * Xc[ISyPICi]);

            apply_flow(ISevPIi, ICrPIi, (rSevPI / rSevPI) * intensities.crit_per_sev[i] / timeISev * Xc[ISevPIi]);
            apply_flow(ISevPIi, TImm2, (1 - (rSevPI / rSevPI) * critSev) / timeISev * Xc[ISevPIi]);
            apply_flow(ISevPIi, DPIi, (rSevPI / rSevPI) * intensities.dead_per_sev[i] / timeISev * Xc[ISevPIi]);

            apply_flow(ICrPIi, DPIi, (rSevPI / rSevPI) * deadICr / timeICr * Xc[ICrPIi]);
            apply_flow(ICrPIi, TImm2, (1 - (rSevPI / rSevPI) * deadICr) / timeICr * Xc[ICrPIi]);

            apply_flow(EIIi, INSIIi, Xc[EIIi] / timeE);
            double ii_TImm = (1 - (rISII / rEII) * (1 - recINS)) / (timeINS * rTimeM);
            double ii_ISy  = (rISII / rEII) * (1 - recINS) / (timeINS * rTimeM);
            apply_flow(INSIIi, TImm2, ii_TImm * Xc[INSIIi]);
            apply_flow(INSIIi, ISyIIi, ii_ISy * Xc[INSIIi]);
            apply_flow(INSIICi, ISyIICi, ii_ISy * Xc[INSIICi]);
            apply_flow(INSIICi, TImm2, ii_TImm * Xc[INSIICi]);

            double ii_ISev     = (rSevII / rISII) * sevISy / (timeISy * rTimeM);
            double ii_ISy_TImm = (1 - (rSevII / rISII) * sevISy) / (timeISy * rTimeM);
            apply_flow(ISyIIi, ISevIIi, ii_ISev * Xc[ISyIIi]);
            apply_flow(ISyIIi, TImm2, ii_ISy_TImm * Xc[ISyIIi]);
            apply_flow(ISyIICi, ISevIIi, ii_ISev * Xc[ISyIICi]);
            apply_flow(ISyIICi, TImm2, ii_ISy_TImm * Xc[ISyIICi]);

            apply_flow(ISevIIi, ICrIIi, (rSevII / rSevII) * intensities.crit_per_sev[i] / timeISev * Xc[ISevIIi]);
            apply_flow(ISevIIi, TImm2, (1 - (rSevII / rSevII) * critSev) / timeISev * Xc[ISevIIi]);
            apply_flow(ISevIIi, DIIi, (rSevII / rSevII) * intensities.dead_per_sev[i] / timeISev * Xc[ISevIIi]);

            apply_flow(ICrIIi, DIIi, (rSevII / rSevII) * deadICr / timeICr * Xc[ICrIIi]);
            apply_flow(ICrIIi, TImm2, (1 - (rSevII / rSevII) * deadICr) / timeICr * Xc[ICrIIi]);

            apply_flow(TImm2, SIIi, Xc[TImm2] / timeTII);
            apply_flow(SIIi, SPIi, Xc[SIIi] / timeWII);
        }
    }
};

} // namespace

inline void flow_based_mobility_returns_secirts_rk1(Eigen::Ref<Eigen::VectorXd> commuter,
                                                    Eigen::Ref<Eigen::VectorXd> totals,
                                                    const mio::osecirts::Model<double>& model, double t, double dt)
{
    mio::osecirts::ModelEvaluatorSecirts eval(model);
    const std::size_t NC = eval.NC;

    Eigen::VectorXd k1_tot(NC), k1_com(NC);
    mio::osecirts::SecirtsIntensities intensities;

    // --- Stage 1 ---
    eval.get_totals_rhs(totals, t, k1_tot);
    eval.get_intensities(totals, t, intensities);
    eval.get_commuter_rhs(intensities, commuter, k1_com);

    // --- Final update ---
    totals += dt * k1_tot;
    commuter += dt * k1_com;
}

inline void flow_based_mobility_returns_secirts_rk2(Eigen::Ref<Eigen::VectorXd> commuter,
                                                    Eigen::Ref<Eigen::VectorXd> totals,
                                                    const mio::osecirts::Model<double>& model, double t, double dt)
{
    mio::osecirts::ModelEvaluatorSecirts eval(model);
    const std::size_t NC = eval.NC;

    Eigen::VectorXd k1_tot(NC), k2_tot(NC);
    Eigen::VectorXd k1_com(NC), k2_com(NC);
    mio::osecirts::SecirtsIntensities intensities;

    // --- Stage 1 ---
    eval.get_totals_rhs(totals, t, k1_tot);
    eval.get_intensities(totals, t, intensities);
    eval.get_commuter_rhs(intensities, commuter, k1_com);

    // --- Stage 2 ---
    Eigen::VectorXd y2 = totals + (dt * 0.5) * k1_tot;
    Eigen::VectorXd c2 = commuter + (dt * 0.5) * k1_com;

    eval.get_totals_rhs(y2, t + 0.5 * dt, k2_tot);
    eval.get_intensities(y2, t + 0.5 * dt, intensities);
    eval.get_commuter_rhs(intensities, c2, k2_com);

    // --- Final update ---
    totals += dt * k2_tot;
    commuter += dt * k2_com;
}

inline void flow_based_mobility_returns_secirts_rk3(Eigen::Ref<Eigen::VectorXd> commuter,
                                                    Eigen::Ref<Eigen::VectorXd> totals,
                                                    const mio::osecirts::Model<double>& model, double t, double dt)
{
    mio::osecirts::ModelEvaluatorSecirts eval(model);
    const std::size_t NC = eval.NC;

    Eigen::VectorXd k1_tot(NC), k2_tot(NC), k3_tot(NC);
    Eigen::VectorXd k1_com(NC), k2_com(NC), k3_com(NC);
    mio::osecirts::SecirtsIntensities intensities;

    // --- Stage 1 ---
    eval.get_totals_rhs(totals, t, k1_tot);
    eval.get_intensities(totals, t, intensities);
    eval.get_commuter_rhs(intensities, commuter, k1_com);

    // --- Stage 2 ---
    Eigen::VectorXd y2 = totals + (dt * 0.5) * k1_tot;
    Eigen::VectorXd c2 = commuter + (dt * 0.5) * k1_com;

    eval.get_totals_rhs(y2, t + 0.5 * dt, k2_tot);
    eval.get_intensities(y2, t + 0.5 * dt, intensities);
    eval.get_commuter_rhs(intensities, c2, k2_com);

    // --- Stage 3 ---
    Eigen::VectorXd y3 = totals + dt * (-k1_tot + 2.0 * k2_tot);
    Eigen::VectorXd c3 = commuter + dt * (-k1_com + 2.0 * k2_com);

    eval.get_totals_rhs(y3, t + dt, k3_tot);
    eval.get_intensities(y3, t + dt, intensities);
    eval.get_commuter_rhs(intensities, c3, k3_com);

    // --- Final update ---
    totals += (dt / 6.0) * (k1_tot + 4.0 * k2_tot + k3_tot);
    commuter += (dt / 6.0) * (k1_com + 4.0 * k2_com + k3_com);
}

inline void flow_based_mobility_returns_secirts_rk4(Eigen::Ref<Eigen::VectorXd> commuter,
                                                    Eigen::Ref<Eigen::VectorXd> totals,
                                                    const mio::osecirts::Model<double>& model, double t, double dt)
{
    mio::osecirts::ModelEvaluatorSecirts eval(model);
    const std::size_t NC = eval.NC;

    Eigen::VectorXd k1_tot(NC), k2_tot(NC), k3_tot(NC), k4_tot(NC);
    Eigen::VectorXd k1_com(NC), k2_com(NC), k3_com(NC), k4_com(NC);
    mio::osecirts::SecirtsIntensities intensities;

    // --- Stage 1 ---
    eval.get_totals_rhs(totals, t, k1_tot);
    eval.get_intensities(totals, t, intensities);
    eval.get_commuter_rhs(intensities, commuter, k1_com);

    // --- Stage 2 ---
    Eigen::VectorXd y2 = totals + (dt * 0.5) * k1_tot;
    Eigen::VectorXd c2 = commuter + (dt * 0.5) * k1_com;

    eval.get_totals_rhs(y2, t + 0.5 * dt, k2_tot);
    eval.get_intensities(y2, t + 0.5 * dt, intensities);
    eval.get_commuter_rhs(intensities, c2, k2_com);

    // --- Stage 3 ---
    Eigen::VectorXd y3 = totals + (dt * 0.5) * k2_tot;
    Eigen::VectorXd c3 = commuter + (dt * 0.5) * k2_com;

    eval.get_totals_rhs(y3, t + 0.5 * dt, k3_tot);
    eval.get_intensities(y3, t + 0.5 * dt, intensities);
    eval.get_commuter_rhs(intensities, c3, k3_com);

    // --- Stage 4 ---
    Eigen::VectorXd y4 = totals + dt * k3_tot;
    Eigen::VectorXd c4 = commuter + dt * k3_com;

    eval.get_totals_rhs(y4, t + dt, k4_tot);
    eval.get_intensities(y4, t + dt, intensities);
    eval.get_commuter_rhs(intensities, c4, k4_com);

    // --- Final update ---
    totals += (dt / 6.0) * (k1_tot + 2.0 * k2_tot + 2.0 * k3_tot + k4_tot);
    commuter += (dt / 6.0) * (k1_com + 2.0 * k2_com + 2.0 * k3_com + k4_com);
}

} // namespace osecirts
} // namespace mio