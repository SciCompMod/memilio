/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
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
#include "abm/state.h"
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
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "IncubationPeriod";
    }
};

struct SusceptibleToExposedByCarrier {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "SusceptibleToExposedByCarrier";
    }
};

struct SusceptibleToExposedByInfected {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "SusceptibleToExposedByInfected";
    }
};

struct CarrierToInfected {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "CarrierToInfected";
    }
};

struct CarrierToRecovered {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "CarrierToRecovered";
    }
};

struct InfectedToRecovered {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "InfectedToRecovered";
    }
};

struct InfectedToSevere {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "InfectedToSevere";
    }
};

struct SevereToCritical {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "SevereToCritical";
    }
};

struct SevereToRecovered {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "SevereToRecovered";
    }
};

struct CriticalToRecovered {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "CriticalToRecovered";
    }
};

struct CriticalToDead {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "CriticalToDead";
    }
};

struct RecoveredToSusceptible {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
    }
    static std::string name()
    {
        return "RecoveredToSusceptible";
    }
};

/**
 * @brief Probability that an Infection is detected.
 */
struct DetectInfection {
    using Type = CustomIndexArray<UncertainValue, AgeGroup, VaccinationState>;
    static Type get_default(AgeGroup size)
    {
        return Type({size, VaccinationState::Count}, 1.);
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
 * @brief Parameters of the Infection that depend on the Location.
 */
using LocalInfectionParameters = ParameterSet<MaximumContacts>;

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
 * @brief The vector of AgeGroups that can go to school.
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
 * @brief The vector of AgeGroups that can go to Work.
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
    ParameterSet<IncubationPeriod, SusceptibleToExposedByCarrier, SusceptibleToExposedByInfected, CarrierToInfected,
                 CarrierToRecovered, InfectedToRecovered, InfectedToSevere, SevereToCritical, SevereToRecovered,
                 CriticalToDead, CriticalToRecovered, RecoveredToSusceptible, DetectInfection, MaskProtection,
                 LockdownDate, SocialEventRate, BasicShoppingRate, WorkRatio, SchoolRatio, GotoWorkTimeMinimum,
                 GotoWorkTimeMaximum, GotoSchoolTimeMinimum, GotoSchoolTimeMaximum, AgeGroupGotoSchool,
                 AgeGroupGotoWork>;

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
     * @return Returns 1 if one constraint is not satisfied, otherwise 0.   
     */
    int check_constraints() const
    {
        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {

            if (this->get<IncubationPeriod>()[{i, VaccinationState::Unvaccinated}] < 0) {
                log_error("Constraint check: Parameter IncubationPeriod of age group {:.0f} smaller {:.4f}", (size_t)i,
                          0);
                return 1;
            }

            if (this->get<SusceptibleToExposedByCarrier>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<SusceptibleToExposedByCarrier>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error(
                    "Constraint check: Parameter SusceptibleToExposedByCarrier of age group {:.0f} smaller {:d} or "
                    "larger {:d}",
                    (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<SusceptibleToExposedByInfected>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<SusceptibleToExposedByInfected>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error(
                    "Constraint check: Parameter SusceptibleToExposedByInfected of age group {:.0f} smaller {:d} or "
                    "larger {:d}",
                    (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<CarrierToInfected>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<CarrierToInfected>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter CarrierToInfected of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<CarrierToRecovered>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<CarrierToRecovered>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter CarrierToInfected of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<InfectedToRecovered>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<InfectedToRecovered>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter InfectedToRecovered of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<InfectedToSevere>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<InfectedToSevere>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter InfectedToSevere of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<SevereToCritical>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<SevereToCritical>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter SevereToCritical of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<SevereToRecovered>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<SevereToRecovered>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter SevereToRecovered of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<CriticalToDead>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<CriticalToDead>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter CriticalToDead of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<CriticalToRecovered>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<CriticalToRecovered>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter CriticalToRecovered of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<RecoveredToSusceptible>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<RecoveredToSusceptible>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter RecoveredToSusceptible of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<DetectInfection>()[{i, VaccinationState::Unvaccinated}] < 0.0 ||
                this->get<DetectInfection>()[{i, VaccinationState::Unvaccinated}] > 1.0) {
                log_error("Constraint check: Parameter DetectInfection of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, 1);
                return 1;
            }

            if (this->get<GotoWorkTimeMinimum>()[i].seconds() < 0.0 ||
                this->get<GotoWorkTimeMinimum>()[i].seconds() > this->get<GotoWorkTimeMaximum>()[i].seconds()) {
                log_error("Constraint check: Parameter GotoWorkTimeMinimum of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, this->get<GotoWorkTimeMaximum>()[i].seconds());
                return 1;
            }

            if (this->get<GotoWorkTimeMaximum>()[i].seconds() < this->get<GotoWorkTimeMinimum>()[i].seconds()) {
                log_error("Constraint check: Parameter GotoWorkTimeMaximum of age group {:.0f} smaller {:d}", (size_t)i,
                          this->get<GotoWorkTimeMinimum>()[i].seconds());
                return 1;
            }

            if (this->get<GotoSchoolTimeMinimum>()[i].seconds() < 0.0 ||
                this->get<GotoSchoolTimeMinimum>()[i].seconds() > this->get<GotoSchoolTimeMaximum>()[i].seconds()) {
                log_error("Constraint check: Parameter GotoSchoolTimeMinimum of age group {:.0f} smaller {:d} or "
                          "larger {:d}",
                          (size_t)i, 0, this->get<GotoWorkTimeMaximum>()[i].seconds());
                return 1;
            }

            if (this->get<GotoSchoolTimeMaximum>()[i].seconds() < this->get<GotoSchoolTimeMinimum>()[i].seconds()) {
                log_error("Constraint check: Parameter GotoWorkTimeMaximum of age group {:.0f} smaller {:d}", (size_t)i,
                          this->get<GotoSchoolTimeMinimum>()[i].seconds());
                return 1;
            }
        }

        if (this->get<MaskProtection>()[MaskType::Count] < 0.0 || this->get<MaskProtection>()[MaskType::Count] > 1.0) {
            log_error("Constraint check: Parameter MaskProtection smaller {:d} or larger {:d}", 0, 1);
            return 1;
        }

        if (this->get<LockdownDate>().seconds() < 0.0) {
            log_error("Constraint check: Parameter LockdownDate smaller {:d}", 0);
            return 1;
        }

        return 0;
    }

private:
    size_t m_num_groups;
};

} // namespace abm
} // namespace mio
#endif
