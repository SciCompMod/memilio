/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth, Khoa Nguyen
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

#include "abm/age.h"
#include "abm/mask_type.h"
#include "abm/time.h"
#include "abm/virus_variant.h"
#include "abm/vaccine.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/contact_matrix.h"
#include <limits>

namespace mio
{
namespace abm
{

/**
 * @brief Time that a Person is infected but not yet infectious.
 */
struct IncubationPeriod {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "IncubationPeriod";
    }
};

struct InfectedNoSymptomsToSymptoms {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "InfectedNoSymptomsToSymptoms";
    }
};

struct InfectedNoSymptomsToRecovered {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "InfectedNoSymptomsToRecovered";
    }
};

struct InfectedSymptomsToRecovered {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "InfectedSymptomsToRecovered";
    }
};

struct InfectedSymptomsToSevere {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "InfectedSymptomsToSevere";
    }
};

struct SevereToCritical {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "SevereToCritical";
    }
};

struct SevereToRecovered {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "SevereToRecovered";
    }
};

struct CriticalToRecovered {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "CriticalToRecovered";
    }
};

struct CriticalToDead {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "CriticalToDead";
    }
};

struct RecoveredToSusceptible {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 0.);
    }
    static std::string name()
    {
        return "RecoveredToSusceptible";
    }
};
/**
 * @brief Parameters for the ViralLoad course. Default values taken as constant values from the average from
 * https://github.com/VirologyCharite/SARS-CoV-2-VL-paper/tree/main
 * Section 3.3.1 or see also supplementary materials Fig. S5.
*/
struct ViralLoadDistributionsParameters {
    UniformDistribution<double>::ParamType viral_load_peak;
    UniformDistribution<double>::ParamType viral_load_incline;
    UniformDistribution<double>::ParamType viral_load_decline;
};

struct ViralLoadDistributions {
    using Type = CustomIndexArray<ViralLoadDistributionsParameters, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        Type default_val({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count},
                         ViralLoadDistributionsParameters{{8.1, 8.1}, {2., 2.}, {-0.17, -0.17}});
        return default_val;
    }
    static std::string name()
    {
        return "ViralLoadDistributions";
    }
};

/**
 * @brief Parameters for the Infectivity. Default values taken as constant values that match the graph 2C from
 * https://github.com/VirologyCharite/SARS-CoV-2-VL-paper/tree/main
*/
struct InfectivityDistributionsParameters {
    UniformDistribution<double>::ParamType infectivity_alpha;
    UniformDistribution<double>::ParamType infectivity_beta;
};

struct InfectivityDistributions {
    using Type = CustomIndexArray<InfectivityDistributionsParameters, VirusVariant, AgeGroup>;
    static Type get_default()
    {
        Type default_val({VirusVariant::Count, AgeGroup::Count},
                         InfectivityDistributionsParameters{{-7., -7.}, {1., 1.}});
        return default_val;
    }
    static std::string name()
    {
        return "InfectivityDistributions";
    }
};

/**
 * @brief Probability that an Infection is detected.
 */
struct DetectInfection {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup, VaccinationState>;
    static Type get_default()
    {
        return Type({VirusVariant::Count, AgeGroup::Count, VaccinationState::Count}, 0.5);
    }
    static std::string name()
    {
        return "DetectInfection";
    }
};

/**
 * @brief Effectiveness of a Mask of a certain MaskType against an Infection.
 */
struct MaskProtection {
    using Type = CustomIndexArray<UncertainValue, MaskType>;
    static auto get_default()
    {
        return Type({MaskType::Count}, 1.);
    }
    static std::string name()
    {
        return "MaskProtection";
    }
};

/**
 * @brief Parameters of the Infection that are the same everywhere within the World.
 */
using GlobalInfectionParameters =
    ParameterSet<IncubationPeriod, InfectedNoSymptomsToSymptoms, InfectedNoSymptomsToRecovered,
                 InfectedSymptomsToRecovered, InfectedSymptomsToSevere, SevereToCritical, SevereToRecovered,
                 CriticalToDead, CriticalToRecovered, RecoveredToSusceptible, ViralLoadDistributions,
                 InfectivityDistributions, DetectInfection, MaskProtection>;

/**
 * @brief Maximum number of Person%s an infectious Person can infect at the respective Location.
 */
struct MaximumContacts {
    using Type = ScalarType;
    static constexpr Type get_default()
    {
        return std::numeric_limits<ScalarType>::max();
    }
    static std::string name()
    {
        return "MaximumContacts";
    }
};

/**
 * contact rates
*/
struct ContactRates {
    using Type = CustomIndexArray<ScalarType, AgeGroup, AgeGroup>;
    static Type get_default()
    {
        return Type({AgeGroup::Count, AgeGroup::Count},
                    1.0); // amount of contacts from AgeGroup a to AgeGroup b per day
    }
    static std::string name()
    {
        return "ContactRates";
    }
};

/**
 * aerosol transmission rates
*/
struct AerosolTransmissionRates {
    using Type = CustomIndexArray<ScalarType, VirusVariant>;
    static Type get_default()
    {
        return Type({VirusVariant::Count}, 1.0); // amount of infections per m^3 per day
    }
    static std::string name()
    {
        return "AerosolTransmissionRates";
    }
};

/**
 * @brief Parameters of the Infection that depend on the Location.
 */
using LocalInfectionParameters = ParameterSet<MaximumContacts, ContactRates, AerosolTransmissionRates>;

/**
 * @brief Parameters that describe the reliability of a test.
 */
struct TestParameters {
    UncertainValue sensitivity;
    UncertainValue specificity;
};

struct GenericTest {
    using Type = TestParameters;
    static Type get_default()
    {
        return Type{0.9, 0.99};
    }
    static std::string name()
    {
        return "GenericTest";
    }
};

/**
 * @brief Reliability of an AntigenTest.
 */
struct AntigenTest : public GenericTest {
    using Type = TestParameters;
    static Type get_default()
    {
        return Type{0.8, 0.88};
    }
    static std::string name()
    {
        return "AntigenTest";
    }
};

/**
 * @brief Reliability of a PCRTest.
 */
struct PCRTest : public GenericTest {
    using Type = TestParameters;
    static Type get_default()
    {
        return Type{0.9, 0.99};
    }
    static std::string name()
    {
        return "PCRTest";
    }
};

/**
 * @brief Starting date of interventions.
 */
struct LockdownDate {
    using Type = TimePoint;
    static auto get_default()
    {
        return TimePoint(std::numeric_limits<int>::max());
    }
    static std::string name()
    {
        return "LockdownDate";
    }
};

/**
 * @brief Parameter for the exponential distribution to decide if a Person goes shopping.
 */
struct BasicShoppingRate {
    using Type = CustomIndexArray<UncertainValue, AgeGroup>;
    static auto get_default()
    {
        return Type({AgeGroup::Count}, 1.0);
    }
    static std::string name()
    {
        return "BasicShoppingRate";
    }
};

/**
 * @brief Percentage of Person%s of the respective age going to work.
 */
struct WorkRatio {
    using Type = DampingMatrixExpression<Dampings<Damping<ColumnVectorShape>>>;
    static auto get_default()
    {
        return Type(Eigen::VectorXd::Constant(1, 1.0));
    }
    static std::string name()
    {
        return "WorkRatio";
    }
};

/**
 * @brief Percentage of Person%s of the respective age going to school.
 */
struct SchoolRatio {
    using Type = DampingMatrixExpression<Dampings<Damping<ColumnVectorShape>>>;
    static auto get_default()
    {
        return Type(Eigen::VectorXd::Constant(1, 1.0));
    }
    static std::string name()
    {
        return "SchoolRatio";
    }
};

/**
 * @brief Parameter for the exponential distribution to decide if a Person goes to a social event.
 */
struct SocialEventRate {
    using Type = DampingMatrixExpression<Dampings<Damping<ColumnVectorShape>>>;
    static auto get_default()
    {
        return Type(Eigen::VectorXd::Constant((size_t)AgeGroup::Count, 1.0));
    }
    static std::string name()
    {
        return "SocialEventRate";
    }
};

/**
 * @brief Earliest time that a Person can go to work.
 */
struct GotoWorkTimeMinimum {
    using Type = CustomIndexArray<TimeSpan, AgeGroup>;
    static auto get_default()
    {
        return CustomIndexArray<TimeSpan, AgeGroup>(AgeGroup::Count, hours(6));
    }
    static std::string name()
    {
        return "GotoWorkTimeMinimum";
    }
};

/**
 * @brief Latest time that a Person can go to work.
 */
struct GotoWorkTimeMaximum {
    using Type = CustomIndexArray<TimeSpan, AgeGroup>;
    static auto get_default()
    {
        return CustomIndexArray<TimeSpan, AgeGroup>(AgeGroup::Count, hours(9));
    }
    static std::string name()
    {
        return "GotoWorkTimeMaximum";
    }
};

/**
 * @brief Earliest time that a Person can go to school.
 */
struct GotoSchoolTimeMinimum {
    using Type = CustomIndexArray<TimeSpan, AgeGroup>;
    static auto get_default()
    {
        return CustomIndexArray<TimeSpan, AgeGroup>(AgeGroup::Count, hours(6));
    }
    static std::string name()
    {
        return "GotoSchoolTimeMinimum";
    }
};

/**
 * @brief Latest time that a Person can go to school.
 */
struct GotoSchoolTimeMaximum {
    using Type = CustomIndexArray<TimeSpan, AgeGroup>;
    static auto get_default()
    {
        return CustomIndexArray<TimeSpan, AgeGroup>(AgeGroup::Count, hours(9));
    }
    static std::string name()
    {
        return "GotoSchoolTimeMaximum";
    }
};

/**
 * @brief Parameters that control the migration between Location%s.
 */
using MigrationParameters =
    ParameterSet<LockdownDate, SocialEventRate, BasicShoppingRate, WorkRatio, SchoolRatio, GotoWorkTimeMinimum,
                 GotoWorkTimeMaximum, GotoSchoolTimeMinimum, GotoSchoolTimeMaximum>;

} // namespace abm
} // namespace mio
#endif
