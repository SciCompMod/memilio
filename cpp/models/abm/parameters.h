/* 
* Copyright (C) 2020-2024 MEmilio
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

#include "abm/mask_type.h"
#include "abm/time.h"
#include "abm/virus_variant.h"
#include "abm/vaccine.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/contact_matrix.h"
#include <limits>
#include <set>

namespace mio
{
namespace abm
{

// Distribution that can be used for the time spend in InfectionStates
using InfectionStateTimesDistributionsParameters = UniformDistribution<double>::ParamType;

/**
 * @brief Time that a Person is infected but not yet infectious in day unit
 */
struct IncubationPeriod {
    using Type = CustomIndexArray<InfectionStateTimesDistributionsParameters, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, InfectionStateTimesDistributionsParameters{1., 1.});
    }
    static std::string name()
    {
        return "IncubationPeriod";
    }
};

/**
* @brief Time that a Person is infected but nonsymptomatic in day unit
*/
struct TimeInfectedNoSymptoms {
    using Type = CustomIndexArray<InfectionStateTimesDistributionsParameters, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, InfectionStateTimesDistributionsParameters{1., 1.});
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptoms";
    }
};

/**
* @brief Time that a Person is infected and symptomatic but
*        who do not need to be hospitalized (yet) in day unit
*/
struct TimeInfectedSymptoms {
    using Type = CustomIndexArray<InfectionStateTimesDistributionsParameters, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, InfectionStateTimesDistributionsParameters{1., 1.});
    }
    static std::string name()
    {
        return "TimeInfectedSymptoms";
    }
};

/**
 * @brief Time that a Person is infected and 'simply' hospitalized in day unit
 */
struct TimeInfectedSevere {
    using Type = CustomIndexArray<InfectionStateTimesDistributionsParameters, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, InfectionStateTimesDistributionsParameters{1., 1.});
    }
    static std::string name()
    {
        return "TimeInfectedSevere";
    }
};

/**
 * @brief Time that a Person is treated by ICU in day unit
 */
struct TimeInfectedCritical {
    using Type = CustomIndexArray<InfectionStateTimesDistributionsParameters, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, InfectionStateTimesDistributionsParameters{1., 1.});
    }
    static std::string name()
    {
        return "TimeInfectedCritical";
    }
};

/**
* @brief the percentage of symptomatic cases
*/
struct SymptomsPerInfectedNoSymptoms {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup>;
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
* @brief the percentage of hospitalized patients per infected patients
*/
struct SeverePerInfectedSymptoms {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup>;
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
* @brief the percentage of ICU patients per hospitalized patients
*/
struct CriticalPerInfectedSevere {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup>;
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
* @brief the percentage of dead patients per ICU patients
*/
struct DeathsPerInfectedCritical {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup>;
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
 * @brief Probability that an Infection is detected.
 */
struct DetectInfection {
    using Type = CustomIndexArray<UncertainValue, VirusVariant, AgeGroup>;
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
    using Type = CustomIndexArray<UncertainValue, MaskType>;
    static Type get_default(AgeGroup /*size*/)
    {
        return Type({MaskType::Count}, 1.);
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
        return Type({VirusVariant::Count}, 1.0); // amount of infections per m^3 per day
    }
    static std::string name()
    {
        return "AerosolTransmissionRates";
    }
};

using InputFunctionForProtectionLevel = std::function<ScalarType(ScalarType)>;

/**
 * @brief Personal protection factor against #Infection% after #Infection and #Vaccination, which depends on #ExposureType,
 * #AgeGroup and #VirusVariant. Its value is between 0 and 1.
 */
struct InfectionProtectionFactor {
    using Type = CustomIndexArray<InputFunctionForProtectionLevel, ExposureType, AgeGroup, VirusVariant>;
    static auto get_default(AgeGroup size)
    {
        return Type({ExposureType::Count, size, VirusVariant::Count}, [](ScalarType /*days*/) -> ScalarType {
            return 0;
        });
    }
    static std::string name()
    {
        return "InfectionProtectionFactor";
    }
};

/**
 * @brief Personal protective factor against severe symptoms after #Infection and #Vaccination, which depends on #ExposureType,
 * #AgeGroup and #VirusVariant. Its value is between 0 and 1.
 */
struct SeverityProtectionFactor {
    using Type = CustomIndexArray<InputFunctionForProtectionLevel, ExposureType, AgeGroup, VirusVariant>;
    static auto get_default(AgeGroup size)
    {
        return Type({ExposureType::Count, size, VirusVariant::Count}, [](ScalarType /*days*/) -> ScalarType {
            return 0;
        });
    }
    static std::string name()
    {
        return "SeverityProtectionFactor";
    }
};

/**
 * @brief Personal protective factor against high viral load. Its value is between 0 and 1.
 */
struct HighViralLoadProtectionFactor {
    using Type = InputFunctionForProtectionLevel;
    static auto get_default()
    {
        return Type([](ScalarType /*days*/) -> ScalarType {
            return 0;
        });
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
    using Type = std::set<AgeGroup>;
    static Type get_default(AgeGroup /*size*/)
    {
        return std::set<AgeGroup>{AgeGroup(1)};
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
    using Type = std::set<AgeGroup>;
    static Type get_default(AgeGroup /*size*/)
    {
        return std::set<AgeGroup>{AgeGroup(2), AgeGroup(3)};
    }
    static std::string name()
    {
        return "AgeGroupGotoWork";
    }
};

using ParametersBase =
    ParameterSet<IncubationPeriod, TimeInfectedNoSymptoms, TimeInfectedSymptoms, TimeInfectedSevere,
                 TimeInfectedCritical, SymptomsPerInfectedNoSymptoms, SeverePerInfectedSymptoms,
                 CriticalPerInfectedSevere, DeathsPerInfectedCritical, ViralLoadDistributions, InfectivityDistributions,
                 DetectInfection, MaskProtection, AerosolTransmissionRates, LockdownDate, SocialEventRate,
                 BasicShoppingRate, WorkRatio, SchoolRatio, GotoWorkTimeMinimum, GotoWorkTimeMaximum,
                 GotoSchoolTimeMinimum, GotoSchoolTimeMaximum, AgeGroupGotoSchool, AgeGroupGotoWork,
                 InfectionProtectionFactor, SeverityProtectionFactor, HighViralLoadProtectionFactor>;

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

/**
 * @brief Parameters of the Infection that depend on the Location.
 */
using LocalInfectionParameters = ParameterSet<MaximumContacts, ContactRates>;

/**
 * @brief Parameters of the simulation that are the same everywhere within the World.
 */
class Parameters : public ParametersBase
{
public:
    Parameters(size_t num_agegroups)
        : ParametersBase(AgeGroup(num_agegroups))
        , m_num_groups(num_agegroups)
    {
    }

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

                if (this->get<IncubationPeriod>()[{v, i}].params.a() < 0) {
                    log_error("Constraint check: Lower end of parameter range IncubationPeriod of virus variant {} and "
                              "age group {:.0f} smaller "
                              "than {:.4f}",
                              v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedNoSymptoms>()[{v, i}].params.a() < 0.0) {
                    log_error("Constraint check: Lower end of parameter range TimeInfectedNoSymptoms of virus variant "
                              "{} and age group {:.0f} smaller "
                              "than {:d}",
                              v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedSymptoms>()[{v, i}].params.a() < 0.0) {
                    log_error("Constraint check: Lower end of parameter range TimeInfectedSymptoms of virus variant {} "
                              "and age group {:.0f} smaller "
                              "than {:d}",
                              v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedSevere>()[{v, i}].params.a() < 0.0) {
                    log_error("Constraint check: Lower end of parameter range TimeInfectedSevere of virus variant {} "
                              "and age group {:.0f} smaller "
                              "than {:d}",
                              v, (size_t)i, 0);
                    return true;
                }

                if (this->get<TimeInfectedCritical>()[{v, i}].params.a() < 0.0) {
                    log_error("Constraint check: Lower end of parameter range TimeInfectedCritical of virus variant {} "
                              "and age group {:.0f} smaller "
                              "than {:d}",
                              v, (size_t)i, 0);
                    return true;
                }

                if (this->get<SymptomsPerInfectedNoSymptoms>()[{v, i}] < 0.0 ||
                    this->get<SymptomsPerInfectedNoSymptoms>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter SymptomsPerInfectedNoSymptoms of virus variant {} and age "
                              "group {:.0f} smaller than {:d} or larger than {:d}",
                              v, (size_t)i, 0, 1);
                    return true;
                }

                if (this->get<SeverePerInfectedSymptoms>()[{v, i}] < 0.0 ||
                    this->get<SeverePerInfectedSymptoms>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter SeverePerInfectedSymptoms of virus variant {} and age group "
                              "{:.0f} smaller than {:d} or larger than {:d}",
                              v, (size_t)i, 0, 1);
                    return true;
                }

                if (this->get<CriticalPerInfectedSevere>()[{v, i}] < 0.0 ||
                    this->get<CriticalPerInfectedSevere>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter CriticalPerInfectedSevere of virus variant {} and age group "
                              "{:.0f} smaller than {:d} or larger than {:d}",
                              v, (size_t)i, 0, 1);
                    return true;
                }

                if (this->get<DeathsPerInfectedCritical>()[{v, i}] < 0.0 ||
                    this->get<DeathsPerInfectedCritical>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter DeathsPerInfectedCritical of age group {:.0f} smaller than "
                              "{:d} or larger than {:d}",
                              v, (size_t)i, 0, 1);
                    return true;
                }

                if (this->get<DetectInfection>()[{v, i}] < 0.0 || this->get<DetectInfection>()[{v, i}] > 1.0) {
                    log_error("Constraint check: Parameter DetectInfection of virus variant {} and age group {:.0f} "
                              "smaller than {:d} or "
                              "larger than {:d}",
                              v, (size_t)i, 0, 1);
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

private:
    size_t m_num_groups;
};

} // namespace abm
} // namespace mio
#endif
