/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef FEEDBACK_SIMULATION_H
#define FEEDBACK_SIMULATION_H

#include <cassert>
#include "memilio/compartments/simulation.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/damping_sampling.h"

namespace mio
{
/**
 * @brief Daily local ICU occupancy per age group.
 * 
 * This parameter stores all historic (local) ICU occupancy values normalized per 100,000 inhabitants.
 * The values are extracted from the simulation results based on the provided ICU compartment indices
 * and are updated at each feedback step.
 */
template <typename FP>
struct ICUOccupancyHistory {
    using Type = mio::TimeSeries<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(size.get());
    }
    static std::string name()
    {
        return "ICUOccupancyHistory";
    }
};

/**
 * @brief Shape parameter of the gamma distribution.
 */
template <typename FP>
struct GammaShapeParameter {
    using Type = FP;
    static Type get_default(AgeGroup)
    {
        return FP(6);
    }
    static std::string name()
    {
        return "GammaShapeParameter";
    }
};

/**
 * @brief Scale parameter of the gamma distribution.
 */
template <typename FP>
struct GammaScaleParameter {
    using Type = FP;
    static Type get_default(AgeGroup)
    {
        return FP(0.4);
    }
    static std::string name()
    {
        return "GammaScaleParameter";
    }
};

/**
 * @brief Number of days in the past considered for the gamma distribution.
 */
struct GammaCutOff {
    using Type = size_t;
    static Type get_default(AgeGroup)
    {
        return 45;
    }
    static std::string name()
    {
        return "GammaCutOff";
    }
};

/**
 * @brief Maximum allowed contact reduction factors per location.
 */
template <typename FP>
struct ContactReductionMax {
    using Type = std::vector<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(1, FP(1.0));
    }
    static std::string name()
    {
        return "ContactReductionMax";
    }
};

/**
 * @brief Minimum allowed contact reduction factors per location.
 */
template <typename FP>
struct ContactReductionMin {
    using Type = std::vector<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(1, FP(0.0));
    }
    static std::string name()
    {
        return "ContactReductionMin";
    }
};

/**
 * @brief Soft-plus curvature parameter for contact adjustment.
 */
template <typename FP>
struct SoftPlusCurvatureParameter {
    using Type = FP;
    static Type get_default(AgeGroup)
    {
        return FP(0.1);
    }
    static std::string name()
    {
        return "SoftPlusCurvatureParameter";
    }
};

/**
 * @brief Nominal ICU capacity.
 */
template <typename FP>
struct NominalICUCapacity {
    using Type = FP;
    static Type get_default(AgeGroup)
    {
        return FP(9.0);
    }
    static std::string name()
    {
        return "NominalICUCapacity";
    }
};

/**
 * @brief Blending factor for local ICU occupancy. 
 * If we only have a local simulation, this factor is set to 1.0 per default.
 */
template <typename FP>
struct BlendingFactorLocal {
    using Type = FP;
    static Type get_default(AgeGroup)
    {
        return FP(1.0);
    }
    static std::string name()
    {
        return "BlendingFactorLocal";
    }
};

/**
 * @brief Blending factor for regional ICU occupancy.
 */
template <typename FP>
struct BlendingFactorRegional {
    using Type = FP;
    static Type get_default(AgeGroup)
    {
        return FP(0.0);
    }
    static std::string name()
    {
        return "BlendingFactorRegional";
    }
};

template <typename FP>
using FeedbackSimulationParameters =
    ParameterSet<ICUOccupancyHistory<FP>, GammaShapeParameter<FP>, GammaScaleParameter<FP>, GammaCutOff,
                 ContactReductionMax<FP>, ContactReductionMin<FP>, SoftPlusCurvatureParameter<FP>,
                 NominalICUCapacity<FP>, BlendingFactorLocal<FP>, BlendingFactorRegional<FP>>;

/**
 * @brief A generic feedback simulation extending existing simulations with a feedback mechanism.
 *
 * This class wraps any simulation (e.g. Simulation or FlowSimulation) and applies additional
 * feedback logic in a model-independent way. Model-specific details—such as the ICU compartments—
 * are provided via arguments.
 *
 * @tparam FP The floating point type.
 * @tparam Sim The simulation type.
 * @tparam ContactPatterns The model-specific contact patterns type.
 */
template <typename FP, typename Sim, typename ContactPatterns>
class FeedbackSimulation
{
public:
    using Model = typename Sim::Model;

    /**
     * @brief Constructs the FeedbackSimulation by taking ownership of an existing simulation instance.
     *
     * @param sim The simulation instance to be extended with feedback mechanism.
     * @param icu_indices A vector of indices indicating ICU compartments for specific model.
     */
    explicit FeedbackSimulation(Sim&& sim, const std::vector<size_t>& icu_indices)
        : m_simulation(std::move(sim))
        , m_icu_indices(icu_indices)
        , m_feedback_parameters(m_simulation.get_model().parameters.get_num_groups())
        , m_perceived_risk(static_cast<Eigen::Index>(m_simulation.get_model().parameters.get_num_groups().get()))
        , m_regional_icu_occupancy(
              static_cast<Eigen::Index>(m_simulation.get_model().parameters.get_num_groups().get()))
        , m_global_icu_occupancy(static_cast<Eigen::Index>(m_simulation.get_model().parameters.get_num_groups().get()))
    {
    }

    /**
     * @brief Advances the simulation until tmax while applying feedback at fixed intervals.
     *
     * The simulation is advanced in steps of dt_feedback. At each step, feedback
     * is applied, then the simulation is advanced, and afterwards the current ICU occupancy is stored.
     * 
     * Note that the simulation may make additional substeps depending on its own
     * timestep dt. When using fixed-step integrators, dt_feedback should be an integer multiple of
     * the simulation timestep dt.
     *
     * @param tmax The maximum simulation time.
     * @param dt_feedback The feedback time step (default 1.0).
     * @return The result in the last time step of the simulation.
     */
    auto advance(const FP tmax, const FP dt_feedback = 1.0)
    {
        FP t = m_simulation.get_result().get_last_time();
        while (t < tmax) {
            FP dt_eff = std::min(dt_feedback, tmax - t);
            apply_feedback(t + dt_eff);
            m_simulation.advance(t + dt_eff);
            add_icu_occupancy(t + dt_eff);
            t += dt_eff;
        }
        return m_simulation.get_result().get_last_value();
    }

    /**
     * @brief Returns the simulation result.
     */
    auto& get_result()
    {
        return m_simulation.get_result();
    }

    /**
     * @brief Returns the model used in the simulation.
     */
    auto& get_model()
    {
        return m_simulation.get_model();
    }

    const auto& get_model() const
    {
        return m_simulation.get_model();
    }

    /**
     * @brief Returns the perceived risk time series.
     */
    auto& get_perceived_risk() const
    {
        return m_perceived_risk;
    }
    auto& get_perceived_risk()
    {
        return m_perceived_risk;
    }

    /**
     * @brief Returns the local feedback parameters.
     */
    auto& get_parameters()
    {
        return m_feedback_parameters;
    }

    /**
     * @brief Calculates the perceived risk based on ICU occupancy data.
     *
     * The risk is computed as a weighted sum of ICU occupancy over recent time points using a gamma-distribution
     * acting as memory kernel.
     *
     * @return The computed (local) perceived risk.
     */
    FP calc_risk_perceived()
    {
        const auto& local_icu_occ = m_feedback_parameters.template get<ICUOccupancyHistory<FP>>();
        size_t num_time_points    = local_icu_occ.get_num_time_points();
        size_t n = std::min(static_cast<size_t>(num_time_points), m_feedback_parameters.template get<GammaCutOff>());
        FP perceived_risk = 0.0;
        const auto& a     = m_feedback_parameters.template get<GammaShapeParameter<FP>>();
        const auto& b     = m_feedback_parameters.template get<GammaScaleParameter<FP>>();

        const auto& blending_local    = m_feedback_parameters.template get<BlendingFactorLocal<FP>>();
        const auto& blending_regional = m_feedback_parameters.template get<BlendingFactorRegional<FP>>();
        const auto& blending_global   = 1 - blending_local - blending_regional;

        // assert that that each blending factor is positive
        assert(blending_local >= 0 && blending_regional >= 0 && blending_global >= 0 &&
               "Blending factors must be non-negative.");

        for (size_t i = num_time_points - n; i < num_time_points; ++i) {
            size_t day   = i - (num_time_points - n);
            FP gamma     = std::pow(b, a) * std::pow(day, a - 1) * std::exp(-b * day) / std::tgamma(a);
            FP perc_risk = 0.0;

            if (blending_local > 0) {
                perc_risk += blending_local * local_icu_occ.get_value(i).sum() /
                             m_feedback_parameters.template get<NominalICUCapacity<FP>>();
            }
            if (blending_regional > 0 && m_regional_icu_occupancy.get_num_time_points() > 0) {
                perc_risk += blending_regional * m_regional_icu_occupancy.get_value(i).sum() /
                             m_feedback_parameters.template get<NominalICUCapacity<FP>>();
            }
            if (blending_global > 0 && m_global_icu_occupancy.get_num_time_points() > 0) {
                perc_risk += blending_global * m_global_icu_occupancy.get_value(i).sum() /
                             m_feedback_parameters.template get<NominalICUCapacity<FP>>();
            }

            perc_risk = std::min(perc_risk, 1.0);
            perceived_risk += perc_risk * gamma;
        }
        return perceived_risk;
    }

    /**
     * @brief Adds the current (local) ICU occupancy into the parameter containing all historical ICU occupancy values.
     *
     * This function use the latest simulation results and extracts the ICU occupancy
     * based on the indices given. The occupancy values are then normalized to  100,000 inhabitants,
     * and then stored (as a new time point) in the ICUOccupancyHistory parameter.
     *
     * @param t The current simulation time at which the ICU occupancy is recorded.
     */
    void add_icu_occupancy(FP t)
    {
        auto& model             = m_simulation.get_model();
        size_t num_groups       = static_cast<size_t>(model.parameters.get_num_groups());
        const auto& last_values = m_simulation.get_result().get_last_value();
        Eigen::VectorXd icu_occ(num_groups);
        for (size_t i = 0; i < m_icu_indices.size(); ++i) {
            icu_occ[i] = last_values[m_icu_indices[i]];
        }
        // normalize by total population and scale to 100,000 inhabitants
        icu_occ = icu_occ / model.populations.get_total() * 100000;
        m_feedback_parameters.template get<ICUOccupancyHistory<FP>>().add_time_point(t, icu_occ);
    }

    /**
     * @brief Sets the regional ICU occupancy time series. 
     * This is used in the graph simulation to initialize and update the regional ICU occupancy.
     * 
     * @param icu_regional The regional ICU occupancy time series.
     */
    void set_regional_icu_occupancy(const mio::TimeSeries<FP>& icu_regional)
    {
        m_regional_icu_occupancy = icu_regional;
    }

    /**
     * @brief Sets the global ICU occupancy time series.
     * This is used in the graph simulation to initialize and update the global ICU occupancy.
     * @param icu_global The global ICU occupancy time series.
     */
    void set_global_icu_occupancy(const mio::TimeSeries<FP>& icu_global)
    {
        m_global_icu_occupancy = icu_global;
    }

    /**
     * @brief Transforms the perceived risk into a contact reduction factor and applies it to the contact patterns.
     *
     * This function computes a contact reduction factor for each location based on the perceived risk. For smooth transitions,
     * we use a softplus function as described in Dönges et al. (doi.org/10.3389/fphy.2022.842180).
     *
     * @param t The simulation time at which the feedback adjustments are applied.
     */
    void apply_feedback(FP t)
    {
        // get model parameters
        auto& params          = m_simulation.get_model().parameters;
        const auto num_groups = (size_t)params.get_num_groups();

        // get feedback parameters
        const auto& contactReductionMax = m_feedback_parameters.template get<ContactReductionMax<FP>>();
        const auto& contactReductionMin = m_feedback_parameters.template get<ContactReductionMin<FP>>();
        const auto epsilon              = m_feedback_parameters.template get<SoftPlusCurvatureParameter<FP>>();

        const size_t num_locations = contactReductionMax.size();

        // calculate local perceived risk
        FP perceived_risk = calc_risk_perceived();

        // add the perceived risk to the time series
        m_perceived_risk.add_time_point(t, Eigen::VectorXd::Constant(num_groups, perceived_risk));

        auto group_weights_all = Eigen::VectorXd::Constant(num_groups, 1.0);

        // define a lambda function to create a damping sampling object given a reduction factor and location
        auto set_contact_reduction = [=](FP val, size_t location) {
            auto v = mio::UncertainValue(val);
            return mio::DampingSampling(v, mio::DampingLevel(0), mio::DampingType(0), mio::SimulationTime(t),
                                        std::vector<size_t>{location}, group_weights_all);
        };

        auto& dampings = params.template get<ContactPatterns>().get_dampings();

        // for each location, compute the effective contact reduction factor and add the corresponding damping sample
        for (size_t loc = 0; loc < num_locations; ++loc) {

            // compute the effective reduction factor using a softplus function.
            FP reduc_fac = (contactReductionMax[loc] - contactReductionMin[loc]) * epsilon *
                               std::log(std::exp(perceived_risk / epsilon) + 1.0) +
                           contactReductionMin[loc];

            // clamp the reduction factor to the maximum allowed value.
            reduc_fac = std::min(reduc_fac, contactReductionMax[loc]);

            // generate the damping for the current location.
            auto damping = set_contact_reduction(reduc_fac, loc);

            // add the damping to the contact pattern.
            dampings.push_back(damping);
        }
        // update the contact matrices after all dampings have been added.
        params.template get<ContactPatterns>().make_matrix();
    }

private:
    Sim m_simulation; ///< The simulation instance.
    std::vector<size_t> m_icu_indices; ///< The ICU compartment indices from the model.
    FeedbackSimulationParameters<FP> m_feedback_parameters; ///< The feedback parameters.
    mio::TimeSeries<FP> m_perceived_risk; ///< The perceived risk time series.
    mio::TimeSeries<FP> m_regional_icu_occupancy; ///< The regional ICU occupancy time series.
    mio::TimeSeries<FP> m_global_icu_occupancy; ///< The global ICU occupancy time series.
};

} // namespace mio

#endif // FEEDBACK_SIMULATION_H
