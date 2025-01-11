/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, Khoa Nguyen, David Kerkmann, Julia Bicker
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
#ifndef MIO_ABM_PARAMETERS_H
#define MIO_ABM_PARAMETERS_H

#include "abm/mask_type.h"
#include "abm/time.h"
#include "abm/virus_variant.h"
#include "abm/protection_event.h"
#include "abm/protection_event.h"
#include "abm/test_type.h"
#include "memilio/config.h"
#include "memilio/io/default_serialize.h"
#include "memilio/io/io.h"
#include "memilio/math/time_series_functor.h"
#include "memilio/utils/parameter_distribution_wrapper.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/index_range.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/contact_matrix.h"

#include <algorithm>
#include <limits>
#include <string>

namespace mio
{
namespace abm
{

// Distribution that can be used for the time spend in InfectionStates
using InfectionStateTimesDistributionsParameters = LogNormalDistribution<double>::ParamType;

/**
 * @brief Time that a Person is infected but not yet infectious in day unit
 */
struct IncubationPeriod {
    using Type = CustomIndexArray<ParameterDistributionWrapper, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ParameterDistributionLogNormal log_norm(1., 1.);
        return Type({VirusVariant::Count, size}, ParameterDistributionWrapper(log_norm));
    }
    static std::string name()
    {
        return "IncubationPeriod";
    }
};

/**
* @brief Time that a Person is infected but presymptomatic in day unit
*/
struct TimeInfectedNoSymptomsToSymptoms {
    using Type = CustomIndexArray<ParameterDistributionWrapper, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ParameterDistributionLogNormal log_norm(1., 1.);
        return Type({VirusVariant::Count, size}, ParameterDistributionWrapper(log_norm));
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptomsToSymptoms";
    }
};

/**
* @brief Time that a Person is infected when staying asymptomatic in day unit
*/
struct TimeInfectedNoSymptomsToRecovered {
    using Type = CustomIndexArray<ParameterDistributionWrapper, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ParameterDistributionLogNormal log_norm(1., 1.);
        return Type({VirusVariant::Count, size}, ParameterDistributionWrapper(log_norm));
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptomsToRecovered";
    }
};

/**
* @brief Time that a Person is infected and symptomatic but
*        who do not need to be hospitalized yet in day unit
*/
struct TimeInfectedSymptomsToSevere {
    using Type = CustomIndexArray<ParameterDistributionWrapper, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ParameterDistributionLogNormal log_norm(1., 1.);

        return Type({VirusVariant::Count, size}, ParameterDistributionWrapper(log_norm));
    }
    static std::string name()
    {
        return "TimeInfectedSymptomsToSevere";
    }
};

/**
* @brief Time that a Person is infected and symptomatic who will recover in day unit
*/
struct TimeInfectedSymptomsToRecovered {
    using Type = CustomIndexArray<ParameterDistributionWrapper, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ParameterDistributionLogNormal log_norm(1., 1.);
        return Type({VirusVariant::Count, size}, ParameterDistributionWrapper(log_norm));
    }
    static std::string name()
    {
        return "TimeInfectedSymptomsToRecovered";
    }
};

/**
 * @brief Time that a Person is infected and 'simply' hospitalized before becoming critical in day unit
 */
struct TimeInfectedSevereToCritical {
    using Type = CustomIndexArray<ParameterDistributionWrapper, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ParameterDistributionLogNormal log_norm(1., 1.);
        return Type({VirusVariant::Count, size}, ParameterDistributionWrapper(log_norm));
    }
    static std::string name()
    {
        return "TimeInfectedSevereToCritical";
    }
};

/**
 * @brief Time that a Person is infected and 'simply' hospitalized before recovering in day unit
 */
struct TimeInfectedSevereToRecovered {
    using Type = CustomIndexArray<ParameterDistributionWrapper, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ParameterDistributionLogNormal log_norm(1., 1.);
        return Type({VirusVariant::Count, size}, ParameterDistributionWrapper(log_norm));
    }
    static std::string name()
    {
        return "TimeInfectedSevereToRecovered";
    }
};

/**
 * @brief Time that a Person is treated by ICU before dying in day unit
 */
struct TimeInfectedCriticalToDead {
    using Type = CustomIndexArray<ParameterDistributionWrapper, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ParameterDistributionLogNormal log_norm(1., 1.);
        return Type({VirusVariant::Count, size}, ParameterDistributionWrapper(log_norm));
    }
    static std::string name()
    {
        return "TimeInfectedCriticalToDead";
    }
};

/**
 * @brief Time that a Person is treated by ICU before recovering in day unit
 */
struct TimeInfectedCriticalToRecovered {
    using Type = CustomIndexArray<ParameterDistributionWrapper, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ParameterDistributionLogNormal log_norm(1., 1.);
        return Type({VirusVariant::Count, size}, ParameterDistributionWrapper(log_norm));
    }
    static std::string name()
    {
        return "TimeInfectedCriticalToRecovered";
    }
};

/**
* @brief the percentage of symptomatic cases
*/
struct SymptomsPerInfectedNoSymptoms {
    using Type = CustomIndexArray<UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, .5);
    }
    static std::string name()
    {
        return "SymptomaticPerInfectedNoSymptoms";
    }
};

/**
* @brief the percentage of hospitalized cases per infected cases
*/
struct SeverePerInfectedSymptoms {
    using Type = CustomIndexArray<UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, .5);
    }
    static std::string name()
    {
        return "SeverePerInfectedSymptoms";
    }
};

/**
* @brief the percentage of ICU cases per hospitalized cases
*/
struct CriticalPerInfectedSevere {
    using Type = CustomIndexArray<UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, .5);
    }
    static std::string name()
    {
        return "CriticalPerInfectedSevere";
    }
};

/**
* @brief the percentage of dead cases per ICU cases
*/
struct DeathsPerInfectedCritical {
    using Type = CustomIndexArray<UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, .5);
    }
    static std::string name()
    {
        return "DeathsPerInfectedCritical";
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

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("ViralLoadDistributionsParameters")
            .add("viral_load_peak", viral_load_peak)
            .add("viral_load_incline", viral_load_incline)
            .add("viral_load_decline", viral_load_decline);
    }
};

struct ViralLoadDistributions {
    using Type = CustomIndexArray<ViralLoadDistributionsParameters, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        Type default_val({VirusVariant::Count, size},
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

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("InfectivityDistributionsParameters")
            .add("infectivity_alpha", infectivity_alpha)
            .add("infectivity_beta", infectivity_beta);
    }
};

struct InfectivityDistributions {
    using Type = CustomIndexArray<InfectivityDistributionsParameters, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        Type default_val({VirusVariant::Count, size}, InfectivityDistributionsParameters{{-7., -7.}, {1., 1.}});
        return default_val;
    }
    static std::string name()
    {
        return "InfectivityDistributions";
    }
};

/**
 * @brief Individual virus shed factor to account for variability in infectious viral load spread.
*/
struct VirusShedFactor {
    using Type = CustomIndexArray<UniformDistribution<double>::ParamType, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        Type default_val({VirusVariant::Count, size}, UniformDistribution<double>::ParamType{0., 0.28});
        return default_val;
    }
    static std::string name()
    {
        return "VirusShedFactor";
    }
};

/**
 * @brief Probability that an Infection is detected.
 */
struct DetectInfection {
    using Type = CustomIndexArray<UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "DetectInfection";
    }
};

/**
 * @brief Effectiveness of a Mask of a certain MaskType% against an Infection%.
 */
struct MaskProtection {
    using Type = CustomIndexArray<UncertainValue<>, MaskType>;
    static Type get_default(AgeGroup /*size*/)
    {
        Type defaut_value = Type(MaskType::Count, 0.0);
        // Initial values according to http://dx.doi.org/10.15585/mmwr.mm7106e1
        defaut_value[MaskType::FFP2]      = 0.83;
        defaut_value[MaskType::Surgical]  = 0.66;
        defaut_value[MaskType::Community] = 0.56;
        return defaut_value;
    }
    static std::string name()
    {
        return "MaskProtection";
    }
};

/**
 * @brief Aerosol transmission rates. 
*/
struct AerosolTransmissionRates {
    using Type = CustomIndexArray<ScalarType, VirusVariant>;
    static Type get_default(AgeGroup /*size*/)
    {
        return Type({VirusVariant::Count}, 1.0);
    }
    static std::string name()
    {
        return "AerosolTransmissionRates";
    }
};

/**
 * @brief Personal protection factor against #Infection% after #Infection and vaccination, which depends on #ProtectionType,
 * #AgeGroup and #VirusVariant. Its value is between 0 and 1.
 */
struct InfectionProtectionFactor {
    using Type = CustomIndexArray<TimeSeriesFunctor<ScalarType>, ProtectionType, AgeGroup, VirusVariant>;
    static auto get_default(AgeGroup size)
    {
        return Type({ProtectionType::Count, size, VirusVariant::Count}, TimeSeriesFunctor<ScalarType>());
    }
    static std::string name()
    {
        return "InfectionProtectionFactor";
    }
};

/**
 * @brief Personal protective factor against severe symptoms after #Infection and vaccination, which depends on #ProtectionType,
 * #AgeGroup and #VirusVariant. Its value is between 0 and 1.
 */
struct SeverityProtectionFactor {
    using Type = CustomIndexArray<TimeSeriesFunctor<ScalarType>, ProtectionType, AgeGroup, VirusVariant>;
    static auto get_default(AgeGroup size)
    {
        return Type({ProtectionType::Count, size, VirusVariant::Count}, TimeSeriesFunctor<ScalarType>());
    }
    static std::string name()
    {
        return "SeverityProtectionFactor";
    }
};

/**
 * @brief Personal protective factor against high viral load, which depends on #ProtectionType,
 * #AgeGroup and #VirusVariant. Its value is between 0 and 1.
 */
struct HighViralLoadProtectionFactor {
    using Type = CustomIndexArray<TimeSeriesFunctor<ScalarType>, ProtectionType, AgeGroup, VirusVariant>;
    static auto get_default(AgeGroup size)
    {
        return Type({ProtectionType::Count, size, VirusVariant::Count}, TimeSeriesFunctor<ScalarType>());
    }
    static std::string name()
    {
        return "HighViralLoadProtectionFactor";
    }
};

/**
 * @brief Parameters that describe the reliability of a test.
 */
struct TestParameters {
    UncertainValue<> sensitivity;
    UncertainValue<> specificity;
    TimeSpan required_time;
    TestType type;

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("TestParameters")
            .add("sensitivity", sensitivity)
            .add("specificity", specificity)
            .add("required_time", required_time)
            .add("test_type", type);
    }
};

/**
 * @brief Store a map from the TestTypes to their TestParameters.
 */
struct TestData {
    using Type = CustomIndexArray<TestParameters, TestType>;
    static auto get_default(AgeGroup /*size*/)
    {
        Type default_val                 = Type({TestType::Count});
        default_val[{TestType::Generic}] = TestParameters{0.9, 0.99, hours(48), TestType::Generic};
        default_val[{TestType::Antigen}] = TestParameters{0.8, 0.88, minutes(30), TestType::Antigen};
        default_val[{TestType::PCR}]     = TestParameters{0.9, 0.99, hours(48), TestType::PCR};
        return default_val;
    }
    static std::string name()
    {
        return "TestData";
    }
};

/**
 * @brief Starting date of interventions.
 */
struct LockdownDate {
    using Type = TimePoint;
    static auto get_default(AgeGroup /*size*/)
    {
        return TimePoint(std::numeric_limits<int>::max());
    }
    static std::string name()
    {
        return "LockdownDate";
    }
};

/**
 * @brief Duration of quarantine.
 */
struct QuarantineDuration {
    using Type = TimeSpan;
    static auto get_default(AgeGroup /*size*/)
    {
        return days(10);
    }
    static std::string name()
    {
        return "QuarantineDuration";
    }
};

/**
 * @brief Parameter for the exponential distribution to decide if a Person goes shopping.
 */
struct BasicShoppingRate {
    using Type = CustomIndexArray<UncertainValue<>, AgeGroup>;
    static auto get_default(AgeGroup size)
    {
        return Type({size}, 1.0);
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
    static auto get_default(AgeGroup /*size*/)
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
    static auto get_default(AgeGroup /*size*/)
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
    static auto get_default(AgeGroup size)
    {
        return Type(Eigen::VectorXd::Constant((size_t)size, 1.0));
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
    static auto get_default(AgeGroup size)
    {
        return CustomIndexArray<TimeSpan, AgeGroup>(size, hours(6));
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
    static auto get_default(AgeGroup size)
    {
        return CustomIndexArray<TimeSpan, AgeGroup>(size, hours(9));
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
    static auto get_default(AgeGroup size)
    {
        return CustomIndexArray<TimeSpan, AgeGroup>(size, hours(6));
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
    static auto get_default(AgeGroup size)
    {
        return CustomIndexArray<TimeSpan, AgeGroup>(size, hours(9));
    }
    static std::string name()
    {
        return "GotoSchoolTimeMaximum";
    }
};

/**
 * @brief The set of AgeGroups that can go to school.
 */
struct AgeGroupGotoSchool {
    using Type = CustomIndexArray<bool, AgeGroup>;
    static Type get_default(AgeGroup num_agegroups)
    {
        return Type(num_agegroups, false);
    }
    static std::string name()
    {
        return "AgeGroupGotoSchool";
    }
};

/**
 * @brief The set of AgeGroups that can go to work.
 */
struct AgeGroupGotoWork {
    using Type = CustomIndexArray<bool, AgeGroup>;
    static Type get_default(AgeGroup num_agegroups)
    {
        return Type(num_agegroups, false);
    }
    static std::string name()
    {
        return "AgeGroupGotoWork";
    }
};

using ParametersBase =
    ParameterSet<IncubationPeriod, TimeInfectedNoSymptomsToSymptoms, TimeInfectedNoSymptomsToRecovered,
                 TimeInfectedSymptomsToSevere, TimeInfectedSymptomsToRecovered, TimeInfectedSevereToCritical,
                 TimeInfectedSevereToRecovered, TimeInfectedCriticalToDead, TimeInfectedCriticalToRecovered,
                 SymptomsPerInfectedNoSymptoms, SeverePerInfectedSymptoms, CriticalPerInfectedSevere,
                 DeathsPerInfectedCritical, ViralLoadDistributions, InfectivityDistributions, VirusShedFactor,
                 DetectInfection, MaskProtection, AerosolTransmissionRates, LockdownDate, QuarantineDuration,
                 SocialEventRate, BasicShoppingRate, WorkRatio, SchoolRatio, GotoWorkTimeMinimum, GotoWorkTimeMaximum,
                 GotoSchoolTimeMinimum, GotoSchoolTimeMaximum, AgeGroupGotoSchool, AgeGroupGotoWork,
                 InfectionProtectionFactor, SeverityProtectionFactor, HighViralLoadProtectionFactor, TestData>;

/**
 * @brief Maximum number of Person%s an infectious Person can infect at the respective Location.
 */
struct MaximumContacts {
    using Type = ScalarType;
    static Type get_default(AgeGroup /*size*/)
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
    static Type get_default(AgeGroup size)
    {
        return Type({size, size},
                    1.0); // amount of contacts from AgeGroup a to AgeGroup b per day
    }
    static std::string name()
    {
        return "ContactRates";
    }
};

// If true, consider the capacity of the Cell%s of this Location for the computation of relative transmission risk.
struct UseLocationCapacityForTransmissions {
    using Type = bool;
    static Type get_default(AgeGroup /*size*/)
    {
        return false;
    }
    static std::string name()
    {
        return "UseLocationCapacityForTransmissions";
    }
};

/**
 * @brief Parameters of the Infection that depend on the Location.
 */
using LocalInfectionParameters = ParameterSet<MaximumContacts, ContactRates, UseLocationCapacityForTransmissions>;

/**
 * @brief Parameters of the simulation that are the same everywhere within the Model.
 */
class Parameters : public ParametersBase
{
public:
    Parameters(size_t num_agegroups)
        : ParametersBase(AgeGroup(num_agegroups))
        , m_num_groups(num_agegroups)
    {
    }

private:
    Parameters(ParametersBase&& base)
        : ParametersBase(std::move(base))
        , m_num_groups(this->get<AgeGroupGotoWork>().size<AgeGroup>().get())
    {
    }

public:
    /**
    * @brief Get the number of the age groups.
    */
    size_t get_num_groups() const
    {
        return m_num_groups;
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and logs an error 
     * if constraints are not satisfied.
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false.
     */
    bool check_constraints() const
    {
        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {
            for (auto&& v : enum_members<VirusVariant>()) {

                if (this->get<IncubationPeriod>()[{v, i}].params()[0] < 0) {
                    log_error("Constraint check: Mean of parameter IncubationPeriod of virus variant {} and "
                              "age group {:.0f} smaller "
                              "than {:.4f}",
                              (uint32_t)v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedNoSymptomsToSymptoms>()[{v, i}].params()[0] < 0.0) {
                    log_error("Constraint check: Mean of parameter TimeInfectedNoSymptomsToSymptoms "
                              "of virus variant "
                              "{} and age group {:.0f} smaller "
                              "than {:d}",
                              (uint32_t)v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedNoSymptomsToRecovered>()[{v, i}].params()[0] < 0.0) {
                    log_error("Constraint check: Mean of parameter TimeInfectedNoSymptomsToRecovered of "
                              "virus variant "
                              "{} and age group {:.0f} smaller "
                              "than {:d}",
                              (uint32_t)v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedSymptomsToSevere>()[{v, i}].params()[0] < 0.0) {
                    log_error("Constraint check: Mean of parameter TimeInfectedSymptomsToSevere of virus "
                              "variant {} "
                              "and age group {:.0f} smaller "
                              "than {:d}",
                              (uint32_t)v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedSymptomsToRecovered>()[{v, i}].params()[0] < 0.0) {
                    log_error("Constraint check: Mean of parameter TimeInfectedSymptomsToRecovered of virus "
                              "variant {} "
                              "and age group {:.0f} smaller "
                              "than {:d}",
                              (uint32_t)v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedSevereToCritical>()[{v, i}].params()[0] < 0.0) {
                    log_error("Constraint check: Mean of parameter TimeInfectedSevereToCritical of virus "
                              "variant {} "
                              "and age group {:.0f} smaller "
                              "than {:d}",
                              (uint32_t)v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedSevereToRecovered>()[{v, i}].params()[0] < 0.0) {
                    log_error("Constraint check: Mean of parameter TimeInfectedSevereToRecovered of virus "
                              "variant {} "
                              "and age group {:.0f} smaller "
                              "than {:d}",
                              (uint32_t)v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedCriticalToDead>()[{v, i}].params()[0] < 0.0) {
                    log_error("Constraint check: Mean of parameter TimeInfectedCriticalToDead of virus variant {} "
                              "and age group {:.0f} smaller "
                              "than {:d}",
                              (uint32_t)v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedCriticalToRecovered>()[{v, i}].params()[0] < 0.0) {
                    log_error("Constraint check: Mean of parameter TimeInfectedCriticalToRecovered of virus "
                              "variant {} "
                              "and age group {:.0f} smaller "
                              "than {:d}",
                              (uint32_t)v, (size_t)i, 0);
                    return true;
                }

                if (this->get<SymptomsPerInfectedNoSymptoms>()[{v, i}] < 0.0 ||
                    this->get<SymptomsPerInfectedNoSymptoms>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter SymptomsPerInfectedNoSymptoms of virus variant {} and age "
                              "group {:.0f} smaller than {:d} or larger than {:d}",
                              (uint32_t)v, (size_t)i, 0, 1);
                    return true;
                }

                if (this->get<SeverePerInfectedSymptoms>()[{v, i}] < 0.0 ||
                    this->get<SeverePerInfectedSymptoms>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter SeverePerInfectedSymptoms of virus variant {} and age group "
                              "{:.0f} smaller than {:d} or larger than {:d}",
                              (uint32_t)v, (size_t)i, 0, 1);
                    return true;
                }

                if (this->get<CriticalPerInfectedSevere>()[{v, i}] < 0.0 ||
                    this->get<CriticalPerInfectedSevere>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter CriticalPerInfectedSevere of virus variant {} and age group "
                              "{:.0f} smaller than {:d} or larger than {:d}",
                              (uint32_t)v, (size_t)i, 0, 1);
                    return true;
                }

                if (this->get<DeathsPerInfectedCritical>()[{v, i}] < 0.0 ||
                    this->get<DeathsPerInfectedCritical>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter DeathsPerInfectedCritical of age group {:.0f} smaller than "
                              "{:d} or larger than {:d}",
                              (uint32_t)v, (size_t)i, 0, 1);
                    return true;
                }

                if (this->get<DetectInfection>()[{v, i}] < 0.0 || this->get<DetectInfection>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter DetectInfection of virus variant {} and age group {:.0f} "
                              "smaller than {:d} or "
                              "larger than {:d}",
                              (uint32_t)v, (size_t)i, 0, 1);
                    return true;
                }
            }

            if (this->get<GotoWorkTimeMinimum>()[i].seconds() < 0.0 ||
                this->get<GotoWorkTimeMinimum>()[i].seconds() > this->get<GotoWorkTimeMaximum>()[i].seconds()) {
                log_error("Constraint check: Parameter GotoWorkTimeMinimum of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, this->get<GotoWorkTimeMaximum>()[i].seconds());
                return true;
            }

            if (this->get<GotoWorkTimeMaximum>()[i].seconds() < this->get<GotoWorkTimeMinimum>()[i].seconds() ||
                this->get<GotoWorkTimeMaximum>()[i] > days(1)) {
                log_error("Constraint check: Parameter GotoWorkTimeMaximum of age group {:.0f} smaller {:d} or larger "
                          "than one day time span",
                          (size_t)i, this->get<GotoWorkTimeMinimum>()[i].seconds());
                return true;
            }

            if (this->get<GotoSchoolTimeMinimum>()[i].seconds() < 0.0 ||
                this->get<GotoSchoolTimeMinimum>()[i].seconds() > this->get<GotoSchoolTimeMaximum>()[i].seconds()) {
                log_error("Constraint check: Parameter GotoSchoolTimeMinimum of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, this->get<GotoWorkTimeMaximum>()[i].seconds());
                return true;
            }

            if (this->get<GotoSchoolTimeMaximum>()[i].seconds() < this->get<GotoSchoolTimeMinimum>()[i].seconds() ||
                this->get<GotoSchoolTimeMaximum>()[i] > days(1)) {
                log_error("Constraint check: Parameter GotoWorkTimeMaximum of age group {:.0f} smaller {:d} or larger "
                          "than one day time span",
                          (size_t)i, this->get<GotoSchoolTimeMinimum>()[i].seconds());
                return true;
            }
        }

        if (this->get<MaskProtection>()[MaskType::Community] < 0.0 ||
            this->get<MaskProtection>()[MaskType::Community] > 1.0) {
            log_error(
                "Constraint check: Parameter MaskProtection for MaskType Community is smaller {:d} or larger {:d}", 0,
                1);
            return true;
        }

        if (this->get<MaskProtection>()[MaskType::FFP2] < 0.0 || this->get<MaskProtection>()[MaskType::FFP2] > 1.0) {
            log_error("Constraint check: Parameter MaskProtection for MaskType FFP2 is smaller {:d} or larger {:d}", 0,
                      1);
            return true;
        }

        if (this->get<MaskProtection>()[MaskType::Surgical] < 0.0 ||
            this->get<MaskProtection>()[MaskType::Surgical] > 1.0) {
            log_error("Constraint check: Parameter MaskProtection for MaskType Surgical smaller {:d} or larger {:d}", 0,
                      1);
            return true;
        }

        if (this->get<LockdownDate>().seconds() < 0.0) {
            log_error("Constraint check: Parameter LockdownDate smaller {:d}", 0);
            return true;
        }

        return false;
    }

    /**
     * deserialize an object of this class.
     * @see epi::deserialize
     */
    template <class IOContext>
    static IOResult<Parameters> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& base, ParametersBase::deserialize(io));
        return success(Parameters(std::move(base)));
    }

private:
    size_t m_num_groups;
};

} // namespace abm
} // namespace mio
#endif
