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
#include "abm/test_type.h"
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

/**
 * @brief Time that a Person is infected but not yet infectious.
 */
struct IncubationPeriod {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "IncubationPeriod";
    }
};

struct InfectedNoSymptomsToSymptoms {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "InfectedNoSymptomsToSymptoms";
    }
};

struct InfectedNoSymptomsToRecovered {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "InfectedNoSymptomsToRecovered";
    }
};

struct InfectedSymptomsToRecovered {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "InfectedSymptomsToRecovered";
    }
};

struct InfectedSymptomsToSevere {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "InfectedSymptomsToSevere";
    }
};

struct SevereToCritical {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "SevereToCritical";
    }
};

struct SevereToRecovered {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "SevereToRecovered";
    }
};

struct CriticalToRecovered {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "CriticalToRecovered";
    }
};

struct CriticalToDead {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
    }
    static std::string name()
    {
        return "CriticalToDead";
    }
};

struct RecoveredToSusceptible {
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type({VirusVariant::Count, size}, 1.);
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
    using Type = CustomIndexArray< UncertainValue<>, VirusVariant, AgeGroup>;
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
    using Type = CustomIndexArray< UncertainValue<>, MaskType>;
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
        return Type({VirusVariant::Count}, 1.0);
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
     UncertainValue<> sensitivity;
     UncertainValue<> specificity;

     /**
     * serialize this. 
     * @see mio::serialize
     */
     template <class IOContext>
     void serialize(IOContext& io) const
     {
         auto obj = io.create_object("TestParameters");
         obj.add_element("Sensitivity", sensitivity);
         obj.add_element("Specificity", specificity);
     }

     /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
     template <class IOContext>
     static IOResult<TestParameters> deserialize(IOContext& io)
     {
         auto obj  = io.expect_object("TestParameters");
         auto sens = obj.expect_element("Sensitivity", mio::Tag<UncertainValue<>>{});
         auto spec = obj.expect_element("Specificity", mio::Tag<UncertainValue<>>{});
         return apply(
             io,
             [](auto&& sens_, auto&& spec_) {
                 return TestParameters{sens_, spec_};
             },
             sens, spec);
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
        default_val[{TestType::Generic}] = TestParameters{0.9, 0.99};
        default_val[{TestType::Antigen}] = TestParameters{0.8, 0.88};
        default_val[{TestType::PCR}]     = TestParameters{0.9, 0.99};
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
    using Type = CustomIndexArray< UncertainValue<>, AgeGroup>;
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
    ParameterSet<IncubationPeriod, InfectedNoSymptomsToSymptoms, InfectedNoSymptomsToRecovered,
                 InfectedSymptomsToRecovered, InfectedSymptomsToSevere, SevereToCritical, SevereToRecovered,
                 CriticalToDead, CriticalToRecovered, RecoveredToSusceptible, ViralLoadDistributions,
                 InfectivityDistributions, DetectInfection, MaskProtection, AerosolTransmissionRates, LockdownDate,
                 QuarantineDuration, SocialEventRate, BasicShoppingRate, WorkRatio, SchoolRatio, GotoWorkTimeMinimum,
                 GotoWorkTimeMaximum, GotoSchoolTimeMinimum, GotoSchoolTimeMaximum, AgeGroupGotoSchool,
                 AgeGroupGotoWork, InfectionProtectionFactor, SeverityProtectionFactor, HighViralLoadProtectionFactor,
                 TestData>;

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

            if (this->get<IncubationPeriod>()[{VirusVariant::Wildtype, i}] < 0) {
                log_error("Constraint check: Parameter IncubationPeriod of age group {:.0f} smaller than {:.4f}",
                          (size_t)i, 0);
                return true;
            }

            if (this->get<InfectedNoSymptomsToSymptoms>()[{VirusVariant::Wildtype, i}] < 0.0) {
                log_error("Constraint check: Parameter InfectedNoSymptomsToSymptoms of age group {:.0f} smaller "
                          "than {:d}",
                          (size_t)i, 0);
                return true;
            }

            if (this->get<InfectedNoSymptomsToRecovered>()[{VirusVariant::Wildtype, i}] < 0.0) {
                log_error("Constraint check: Parameter InfectedNoSymptomsToRecovered of age group {:.0f} smaller "
                          "than {:d}",
                          (size_t)i, 0);
                return true;
            }

            if (this->get<InfectedSymptomsToRecovered>()[{VirusVariant::Wildtype, i}] < 0.0) {
                log_error(
                    "Constraint check: Parameter InfectedSymptomsToRecovered of age group {:.0f} smaller than {:d}",
                    (size_t)i, 0);
                return true;
            }

            if (this->get<InfectedSymptomsToSevere>()[{VirusVariant::Wildtype, i}] < 0.0) {
                log_error("Constraint check: Parameter InfectedSymptomsToSevere of age group {:.0f} smaller than {:d}",
                          (size_t)i, 0);
                return true;
            }

            if (this->get<SevereToCritical>()[{VirusVariant::Wildtype, i}] < 0.0) {
                log_error("Constraint check: Parameter SevereToCritical of age group {:.0f} smaller than {:d}",
                          (size_t)i, 0);
                return true;
            }

            if (this->get<SevereToRecovered>()[{VirusVariant::Wildtype, i}] < 0.0) {
                log_error("Constraint check: Parameter SevereToRecovered of age group {:.0f} smaller than {:d}",
                          (size_t)i, 0);
                return true;
            }

            if (this->get<CriticalToDead>()[{VirusVariant::Wildtype, i}] < 0.0) {
                log_error("Constraint check: Parameter CriticalToDead of age group {:.0f} smaller than {:d}", (size_t)i,
                          0);
                return true;
            }

            if (this->get<CriticalToRecovered>()[{VirusVariant::Wildtype, i}] < 0.0) {
                log_error("Constraint check: Parameter CriticalToRecovered of age group {:.0f} smaller than {:d}",
                          (size_t)i, 0);
                return true;
            }

            if (this->get<RecoveredToSusceptible>()[{VirusVariant::Wildtype, i}] < 0.0) {
                log_error("Constraint check: Parameter RecoveredToSusceptible of age group {:.0f} smaller than {:d}",
                          (size_t)i, 0);
                return true;
            }

            if (this->get<DetectInfection>()[{VirusVariant::Wildtype, i}] < 0.0 ||
                this->get<DetectInfection>()[{VirusVariant::Wildtype, i}] > 1.0) {
                log_error("Constraint check: Parameter DetectInfection of age group {:.0f} smaller than {:d} or "
                          "larger than {:d}",
                          (size_t)i, 0, 1);
                return true;
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
