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
#ifndef SEIRV_PARAMETERS_H
#define SEIRV_PARAMETERS_H

#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/parameter_set.h"

namespace mio
{
namespace oseirv
{

/**
 * @brief Baseline transmissibility R_e (dimensionless).
 *
 * Multiplies the force of infection λ(t) after contact-matrix normalization. Main scalar for transmission
 * intensity. Default: 1.0.
 */
template <class FP>
struct BaselineTransmissibility {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "BaselineTransmissibility";
    }
};

/**
 * @brief Transition/recovery rate γ in 1/week.
 *
 * Controls E→I and I→Rand is included in λ(t) to decouple from R_e.
 * Use the model time unit (weeks). Default: 1/2.0
 */
template <class FP>
struct Gamma {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(1.0 / 2.0);
    }
    static std::string name()
    {
        return "Gamma";
    }
};

/**
 * @brief Seasonality amplitude δ (dimensionless).
 *
 * Amplitude of the seasonal modulation exp(δ·sin(2π(t/52 − t_z + t_s))) in λ(t).
 * Set to 0.0 to disable seasonality. Default: 0.0.
 */
template <class FP>
struct SeasonalityAmplitude {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "SeasonalityAmplitude";
    }
};

/**
 * @brief Subtype-specific seasonal shift t_z.
 *
 * Phase shift for the sinusoid in λ(t); may be adjusted by a seasonal-shift correction if implemented. Default: 0.0.
 */
template <class FP>
struct ShiftTZ {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "ShiftTZ";
    }
};

/**
 * @brief Season-specific fine shift t_s.
 *
 * Additional small seasonal timing correction per season. Default: 0.0.
 */
template <class FP>
struct ShiftTS {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "ShiftTS";
    }
};

/**
 * @brief Outside force of infection λ₀ in 1/week.
 *
 * Additive external hazard that can seed infections when no infectives are present.
 * Default: 0.0.
 */
template <class FP>
struct OutsideFoI {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "OutsideFoI";
    }
};

/**
 * @brief Clustering/concavity parameter ρ (dimensionless).
 *
 * Exponent on the infectious fraction in λ(t). ρ=1 recovers mass-action; 0<ρ<1 yields concave,
 * clustered transmission. Must be > 0. Default: 1.0.
 */
template <class FP>
struct ClusteringRho {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "ClusteringRho";
    }
};

/**
 * @brief Mixing parameter m for “sick” contacts (dimensionless).
 *
 * Weight for the symptomatically sick contact matrix in β_eff = ν(m)^{-1}(β_healthy + m·β_sick). 
 * Combines share of symptomatic cases and their higher infectiousness. Default: 1.0.
 */
template <class FP>
struct SickMixingM {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "SickMixingM";
    }
};

/**
 * @brief Contact patterns of healthy people (age-structured contact frequencies).
 *
 * Used in β_eff, can be time-varying via damping. Unit: contacts per model time step (week).
 */
template <class FP>
struct ContactPatternsHealthy {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "ContactPatternsHealthy";
    }
};

/**
 * @brief Contact patterns of symptomatically sick people (age-structured contact frequencies).
 *
 * Used in β_eff together with m, can be time-varying via damping.
 * Unit: contacts per model time step (week).
 */
template <class FP>
struct ContactPatternsSick {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "ContactPatternsSick";
    }
};

/**
 * @brief Age-specific baseline susceptibility σ_i.
 *
 * Represents pre-existing subtype/age-specific immunity, used with Φ at t0 to form φ·σ_i.
 * Default: 1.0 (fully susceptible).
 */
template <class FP>
struct SigmaByAge {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "SigmaByAge";
    }
};

/**
 * @brief Vaccination coverage VC_i at t0 (dimensionless, in [0,1]).
 *
 * Age-specific share vaccinated at season start; used in initial conditions.
 * Default: 0.0.
 */
template <class FP>
struct VaccineCoverage {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "VaccineCoverage";
    }
};

/**
 * @brief Vaccine effectiveness VE_i (dimensionless, in [0,1]).
 *
 * Reduces effective susceptibility of vaccinated individuals in the initial state via (1 − VE_i). Default: 0.0.
 */
template <class FP>
struct VaccineEffectiveness {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "VaccineEffectiveness";
    }
};

/**
 * @brief Φ (phi): fraction of the population that remains susceptible at t0 (dimensionless, typically in [0,1]).
 *
 * Scales the susceptible pool at season start: φ·σ_i). Default: 1.0.
 */
template <class FP>
struct Phi {
    using Type = UncertainValue<FP>;
    static Type get_default(AgeGroup)
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "Phi";
    }
};

template <class FP>
using ParametersBase =
    ParameterSet<BaselineTransmissibility<FP>, Gamma<FP>, SeasonalityAmplitude<FP>, ShiftTZ<FP>, ShiftTS<FP>,
                 OutsideFoI<FP>, ClusteringRho<FP>, SickMixingM<FP>, ContactPatternsHealthy<FP>,
                 ContactPatternsSick<FP>, SigmaByAge<FP>, VaccineCoverage<FP>, VaccineEffectiveness<FP>, Phi<FP>>;

/**
 * @brief Parameter set for the age-resolved SEIRV model (S,E,I,R plus vaccinated states) as per the appendix.
 *
 * Contains seasonality, contact-pattern and immunity/vaccination parameters.
 * Model time unit is week; contact matrices may be time-dependent (damping).
 */
template <class FP>
class Parameters : public ParametersBase<FP>
{
public:
    /**
     * @brief Construct with the number of age groups.
     * @param ng Number of age groups.
     */
    Parameters(AgeGroup ng)
        : ParametersBase<FP>(ng)
        , m_num_groups{ng}
    {
    }

    /**
     * @brief Returns the number of age groups.
     */
    AgeGroup get_num_groups() const
    {
        return m_num_groups;
    }

    /**
     * @brief Apply minimal constraints to keep parameters valid.
     *
     * @return true if one or more parameters were corrected; false otherwise.
     */
    bool apply_constraints()
    {
        bool corrected = false;
        if (this->template get<ClusteringRho<FP>>() <= 0.0) {
            this->template get<ClusteringRho<FP>>() = 1.0;
            corrected                               = true;
        }
        if (this->template get<BaselineTransmissibility<FP>>() < 0.0) {
            this->template get<BaselineTransmissibility<FP>>() = 0.0;
            corrected                                          = true;
        }
        if (this->template get<OutsideFoI<FP>>() < 0.0) {
            this->template get<OutsideFoI<FP>>() = 0.0;
            corrected                            = true;
        }
        return corrected;
    }

private:
    AgeGroup m_num_groups;
};

} // namespace oseirv
} // namespace mio
#endif
