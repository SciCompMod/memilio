/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef EPI_ABM_PARAMETERS_H
#define EPI_ABM_PARAMETERS_H

#include "epidemiology/abm/age.h"
#include "state.h"
#include "epidemiology/abm/time.h"
#include "epidemiology/utils/custom_index_array.h"
#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/parameter_set.h"
#include "epidemiology/secir/damping.h"
#include "epidemiology/secir/contact_matrix.h"

#include <limits>

namespace epi
{

struct IncubationPeriod
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct SusceptibleToExposedByCarrier
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct SusceptibleToExposedByInfected
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct CarrierToInfected
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct CarrierToRecovered
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct InfectedToRecovered
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct InfectedToSevere
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct SevereToCritical
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct SevereToRecovered
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct CriticalToRecovered
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};


struct CriticalToDead
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct RecoveredToSusceptible
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 1.);
    }
};

struct DetectInfection
{
    using Type = CustomIndexArray<double, AbmAgeGroup, VacinationState>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count, VacinationState::Count}, 0.5);
    }
};


/**
 * parameters of the infection that are the same everywhere within the world.
 */
using GlobalInfectionParameters = ParameterSet<IncubationPeriod,
                                               SusceptibleToExposedByCarrier,
                                               SusceptibleToExposedByInfected,
                                               CarrierToInfected,
                                               CarrierToRecovered,
                                               InfectedToRecovered,
                                               InfectedToSevere,
                                               SevereToCritical,
                                               SevereToRecovered,
                                               CriticalToDead,
                                               CriticalToRecovered,
                                               RecoveredToSusceptible,
                                               DetectInfection>;


struct EffectiveContacts
{
    using Type = double;
    static constexpr Type get_default()
    {
        return std::numeric_limits<double>::max();
    }
};

/**
 * parameters of the infection that depend on the location.
 */
using LocalInfectionParameters = ParameterSet<EffectiveContacts>;

/**
 * parameters that govern the migration between locations
 */
struct LockdownDate {
    using Type = TimePoint;
    static auto get_default()
    {
        return TimePoint(std::numeric_limits<int>::max());
    }
};
struct BasicShoppingRate {
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static auto get_default()
    {
        return CustomIndexArray<double, AbmAgeGroup>(AbmAgeGroup::Count, 1.0);
    }
};
struct WorkRatio {
    using Type = DampingMatrixExpression<Dampings<Damping<ColumnVectorShape>>>;
    static auto get_default()
    {
        return Type(Eigen::VectorXd::Constant(1,1.0));
    }
};
struct SchoolRatio {
    using Type = DampingMatrixExpression<Dampings<Damping<ColumnVectorShape>>>;
    static auto get_default()
    {
        return Type(Eigen::VectorXd::Constant(1,1.0));
    }
};
struct SocialEventRate {
    using Type = DampingMatrixExpression<Dampings<Damping<ColumnVectorShape>>>;
    static auto get_default()
    {
        return Type(Eigen::VectorXd::Constant((size_t)AbmAgeGroup::Count,1.0));
    }
};

/**
 * parameters that control the migration between locations.
 */
using AbmMigrationParameters =
    ParameterSet<LockdownDate, SocialEventRate, BasicShoppingRate, WorkRatio, SchoolRatio>;

} // namespace epi
#endif
