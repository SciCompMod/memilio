/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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
#ifndef SECIR_PARAMS_H
#define SECIR_PARAMS_H

#include "memilio/math/eigen.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/math/adapt_rk.h"
#include "secir/age_group.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/dynamic_npis.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/custom_index_array.h"

#include <vector>

namespace mio
{
namespace vaccinated
{

    /*******************************************
 * Define Parameters of the SECIHURD model *
 *******************************************/

    /**
 * @brief the start day in the SECIR model
 * The start day defines in which season the simulation can be started
 * If the start day is 180 and simulation takes place from t0=0 to
 * tmax=100 the days 180 to 280 of the year are simulated
 */
    struct StartDay {
        using Type = double;
        static Type get_default(AgeGroup)
        {
            return 0.;
        }
        static std::string name()
        {
            return "StartDay";
        }
    };

    /**
 * @brief the seasonality in the SECIR model
 * the seasonality is given as (1+k*sin()) where the sine
 * curve is below one in summer and above one in winter
 */
    struct Seasonality {
        using Type = UncertainValue;
        static Type get_default(AgeGroup)
        {
            return Type(0.);
        }
        static std::string name()
        {
            return "Seasonality";
        }
    };

    /**
 * @brief the icu capacity in the SECIR model
 */
    struct ICUCapacity {
        using Type = UncertainValue;
        static Type get_default(AgeGroup)
        {
            return Type(std::numeric_limits<double>::max());
        }
        static std::string name()
        {
            return "ICUCapacity";
        }
    };

    /**
 * @brief the incubation time in the SECIR model
 * @param tinc incubation time in day unit
 */
    struct IncubationTime {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "IncubationTime";
        }
    };

    /**
 * @brief the infectious time for symptomatic cases that are infected but
 *        who do not need to be hsopitalized in the SECIR model in day unit
 */
    struct InfectiousTimeMild {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "InfectiousTimeMild";
        }
    };

    /**
 * @brief the serial interval in the SECIR model in day unit
 */
    struct SerialInterval {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "SerialInterval";
        }
    };

    /**
 * @brief the time people are 'simply' hospitalized before returning home in the SECIR model
 *        in day unit
 */
    struct HospitalizedToHomeTime {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "HospitalizedToHomeTime";
        }
    };

    /**
 * @brief the time people are infectious at home before 'simply' hospitalized in the SECIR model
 *        in day unit
 */
    struct HomeToHospitalizedTime {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "HomeToHospitalizedTime";
        }
    };

    /**
 * @brief the time people are 'simply' hospitalized before being treated by ICU in the SECIR model
 *        in day unit
 */
    struct HospitalizedToICUTime {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "HospitalizedToICUTime";
        }
    };

    /**
 * @brief the time people are treated by ICU before returning home in the SECIR model
 *        in day unit
 */
    struct ICUToHomeTime {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "ICUToHomeTime";
        }
    };

    /**
 * @brief the infectious time for asymptomatic cases in the SECIR model
 *        in day unit
 */
    struct InfectiousTimeAsymptomatic {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "InfectiousTimeAsymptomatic";
        }
    };

    /**
 * @brief the time people are treated by ICU before dying in the SECIR model
 *        in day unit
 */
    struct ICUToDeathTime {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "ICUToDeathTime";
        }
    };

    /**
* @brief probability of getting infected from a contact
*/
    struct InfectionProbabilityFromContact {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "InfectionProbabilityFromContact";
        }
    };

    /**
* @brief the relative carrier infectability
*/
    struct RelativeCarrierInfectability {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 1.);
        }
        static std::string name()
        {
            return "RelativeCarrierInfectability";
        }
    };

    /**
* @brief the percentage of asymptomatic cases in the SECIR model
*/
    struct AsymptoticCasesPerInfectious {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.);
        }
        static std::string name()
        {
            return "AsymptoticCasesPerInfectious";
        }
    };

    /**
* @brief the risk of infection from symptomatic cases in the SECIR model
*/
    struct RiskOfInfectionFromSympomatic {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.);
        }
        static std::string name()
        {
            return "RiskOfInfectionFromSympomatic";
        }
    };

    /**
* @brief risk of infection from symptomatic cases increases as test and trace capacity is exceeded.
*/
    struct MaxRiskOfInfectionFromSympomatic {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.);
        }
        static std::string name()
        {
            return "MaxRiskOfInfectionFromSympomatic";
        }
    };

    /**
* @brief the percentage of hospitalized patients per infected patients in the SECIR model
*/
    struct HospitalizedCasesPerInfectious {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.);
        }
        static std::string name()
        {
            return "HospitalizedCasesPerInfectious";
        }
    };

    /**
* @brief the percentage of ICU patients per hospitalized patients in the SECIR model
*/
    struct ICUCasesPerHospitalized {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.);
        }
        static std::string name()
        {
            return "ICUCasesPerHospitalized";
        }
    };

    /**
* @brief the percentage of dead patients per ICU patients in the SECIR model
*/
    struct DeathsPerHospitalized {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.);
        }
        static std::string name()
        {
            return "DeathsPerHospitalized";
        }
    };

    /**
 * @brief the contact patterns within the society are modelled using an UncertainContactMatrix
 */
    struct ContactPatterns {
        using Type = UncertainContactMatrix;
        static Type get_default(AgeGroup size)
        {
            return Type(1, static_cast<Eigen::Index>((size_t)size));
        }
        static std::string name()
        {
            return "ContactPatterns";
        }
    };

    /**
 * @brief the NPIs that are enacted if certain infection thresholds are exceeded.
 */
    struct DynamicNPIsInfected {
        using Type = DynamicNPIs;
        static Type get_default(AgeGroup /*size*/)
        {
            return {};
        }
        static std::string name()
        {
            return "DynamicNPIsInfected";
        }
    };

    /**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
    struct TestAndTraceCapacity {
        using Type = UncertainValue;
        static Type get_default(AgeGroup)
        {
            return Type(std::numeric_limits<double>::max());
        }
        static std::string name()
        {
            return "TestAndTraceCapacity";
        }
    };

    /**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
    struct VaccineGrowthFirst {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "VaccineGrowthFirst";
        }
    };

    /**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
    struct VaccineGrowthFull {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "VaccineGrowthFull";
        }
    };

    /**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
    struct BaseInfB117 {
        using Type = CustomIndexArray<double, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "BaseInfB117";
        }
    };

    struct BaseInfB161 {
        using Type = CustomIndexArray<double, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "BaseInfB161";
        }
    };

    struct ReducVaccExp {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "ReducVaccExp";
        }
    };

    struct ReducImmuneExp {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "ReducImmuneExp";
        }
    };

    struct ReducExpInf {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "ReducExpInf";
        }
    };

    struct ReducImmuneExpInf {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "ReducImmuneExpInf";
        }
    };

    struct ReducInfHosp {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "ReducInfHosp";
        }
    };

    struct ReducImmuneInfHosp {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.0);
        }
        static std::string name()
        {
            return "ReducImmuneInfHosp";
        }
    };

    struct ReducTime {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 0.5);
        }
        static std::string name()
        {
            return "ReducTime";
        }
    };

    /**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
    struct VaccinationGap {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 49.0);
        }
        static std::string name()
        {
            return "VaccinationGap";
        }
    };

    /**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
    struct DaysUntilEffective {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 14.0);
        }
        static std::string name()
        {
            return "DaysUntilEffective";
        }
    };

    /**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
    struct DaysUntilEffectiveFull {
        using Type = CustomIndexArray<UncertainValue, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size, 7.0);
        }
        static std::string name()
        {
            return "DaysUntilEffectiveFull";
        }
    };

    /**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
    struct DailyFirstVaccination {
        using Type = CustomIndexArray<std::vector<double>, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size);
        }
        static std::string name()
        {
            return "DailyFirstVaccination";
        }
    };

    /**
 * @brief capacity to test and trace contacts of infected for quarantine per day.
 */
    struct DailyFullVaccination {
        using Type = CustomIndexArray<std::vector<double>, AgeGroup>;
        static Type get_default(AgeGroup size)
        {
            return Type(size);
        }
        static std::string name()
        {
            return "DailyFullVaccination";
        }
    };

    using SecirParamsBase =
        ParameterSet<StartDay, Seasonality, ICUCapacity, TestAndTraceCapacity, ContactPatterns, DynamicNPIsInfected,
                     IncubationTime, InfectiousTimeMild, InfectiousTimeAsymptomatic, SerialInterval,
                     HospitalizedToHomeTime, HomeToHospitalizedTime, HospitalizedToICUTime, ICUToHomeTime,
                     ICUToDeathTime, InfectionProbabilityFromContact, RelativeCarrierInfectability,
                     AsymptoticCasesPerInfectious, RiskOfInfectionFromSympomatic, MaxRiskOfInfectionFromSympomatic,
                     HospitalizedCasesPerInfectious, ICUCasesPerHospitalized, DeathsPerHospitalized, VaccineGrowthFirst,
                     VaccineGrowthFull, VaccinationGap, DaysUntilEffective, DaysUntilEffectiveFull, BaseInfB117,
                     BaseInfB161, DailyFullVaccination, DailyFirstVaccination, ReducVaccExp, ReducImmuneExp,
                     ReducExpInf, ReducImmuneExpInf, ReducInfHosp, ReducImmuneInfHosp, ReducTime>;

    /**
 * @brief Parameters of an age-resolved SECIR/SECIHURD model.
 */
    class SecirParams : public SecirParamsBase
    {
    public:
        SecirParams(AgeGroup num_agegroups)
            : SecirParamsBase(num_agegroups)
            , m_num_groups{num_agegroups}
        {
        }

        AgeGroup get_num_groups() const
        {
            return m_num_groups;
        }

        /**
     * @brief checks whether all Parameters satisfy their corresponding constraints and applies them, if they do not
     */
        void apply_constraints()
        {
            if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
                log_warning("Constraint check: Parameter Seasonality changed from {:0.4f} to {:d}",
                            this->get<Seasonality>(), 0);
                this->set<Seasonality>(0);
            }

            if (this->get<ICUCapacity>() < 0.0) {
                log_warning("Constraint check: Parameter ICUCapacity changed from {:0.4f} to {:d}",
                            this->get<ICUCapacity>(), 0);
                this->set<ICUCapacity>(0);
            }

            for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {

                if (this->get<IncubationTime>()[i] < 2.0) {
                    log_warning("Constraint check: Parameter IncubationTime changed from {:.4f} to {:.4f}",
                                this->get<IncubationTime>()[i], 2.0);
                    this->get<IncubationTime>()[i] = 2.0;
                }

                if (2 * this->get<SerialInterval>()[i] < this->get<IncubationTime>()[i] + 1.0) {
                    log_warning("Constraint check: Parameter SerialInterval changed from {:.4f} to {:.4f}",
                                this->get<SerialInterval>()[i], 0.5 * this->get<IncubationTime>()[i] + 0.5);
                    this->get<SerialInterval>()[i] = 0.5 * this->get<IncubationTime>()[i] + 0.5;
                }
                else if (this->get<SerialInterval>()[i] > this->get<IncubationTime>()[i] - 0.5) {
                    log_warning("Constraint check: Parameter SerialInterval changed from {:.4f} to {:.4f}",
                                this->get<SerialInterval>()[i], this->get<IncubationTime>()[i] - 0.5);
                    this->get<SerialInterval>()[i] = this->get<IncubationTime>()[i] - 0.5;
                }

                if (this->get<InfectiousTimeMild>()[i] < 1.0) {
                    log_warning("Constraint check: Parameter InfectiousTimeMild changed from {:.4f} to {:.4f}",
                                this->get<InfectiousTimeMild>()[i], 1.0);
                    this->get<InfectiousTimeMild>()[i] = 1.0;
                }

                if (this->get<HospitalizedToHomeTime>()[i] < 1.0) {
                    log_warning("Constraint check: Parameter HospitalizedToHomeTime changed from {:.4f} to {:.4f}",
                                this->get<HospitalizedToHomeTime>()[i], 1.0);
                    this->get<HospitalizedToHomeTime>()[i] = 1.0;
                }

                if (this->get<HomeToHospitalizedTime>()[i] < 1.0) {
                    log_warning("Constraint check: Parameter HomeToHospitalizedTime changed from {:.4f} to {:.4f}",
                                this->get<HomeToHospitalizedTime>()[i], 1.0);
                    this->get<HomeToHospitalizedTime>()[i] = 1.0;
                }

                if (this->get<HospitalizedToICUTime>()[i] < 1.0) {
                    log_warning("Constraint check: Parameter HospitalizedToICUTime changed from {:.4f} to {:.4f}",
                                this->get<HospitalizedToICUTime>()[i], 1.0);
                    this->get<HospitalizedToICUTime>()[i] = 1.0;
                }

                if (this->get<ICUToHomeTime>()[i] < 1.0) {
                    log_warning("Constraint check: Parameter ICUToHomeTime changed from {:.4f} to {:.4f}",
                                this->get<ICUToHomeTime>()[i], 1.0);
                    this->get<ICUToHomeTime>()[i] = 1.0;
                }

                if (this->get<ICUToDeathTime>()[i] < 1.0) {
                    log_warning("Constraint check: Parameter ICUToDeathTime changed from {:.4f} to {:.4f}",
                                this->get<ICUToDeathTime>()[i], 1.0);
                    this->get<ICUToDeathTime>()[i] = 1.0;
                }

                auto t_inf_asymp = 1.0 / (0.5 / (this->get<IncubationTime>()[i] - this->get<SerialInterval>()[i])) +
                                   0.5 * this->get<InfectiousTimeMild>()[i];
                if (abs(this->get<InfectiousTimeAsymptomatic>()[i] - t_inf_asymp) > 1e-12) {
                    log_info("Constraint check: Parameter InfectiousTimeAsymptomatic set as fully dependent on tinc, "
                             "tserint and tinfmild, as proposed by "
                             "https://www.medrxiv.org/content/10.1101/2020.04.04.20053637v1.");
                    this->get<InfectiousTimeAsymptomatic>()[i] = t_inf_asymp;
                }

                if (this->get<InfectionProbabilityFromContact>()[i] < 0.0) {
                    log_warning(
                        "Constraint check: Parameter InfectionProbabilityFromContact changed from {:0.4f} to {:d} ",
                        this->get<InfectionProbabilityFromContact>()[i], 0);
                    this->get<InfectionProbabilityFromContact>()[i] = 0;
                }

                if (this->get<RelativeCarrierInfectability>()[i] < 0.0) {
                    log_warning(
                        "Constraint check: Parameter RelativeCarrierInfectability changed from {:0.4f} to {:d} ",
                        this->get<RelativeCarrierInfectability>()[i], 0);
                    this->get<RelativeCarrierInfectability>()[i] = 0;
                }

                if (this->get<AsymptoticCasesPerInfectious>()[i] < 0.0 ||
                    this->get<AsymptoticCasesPerInfectious>()[i] > 1.0) {
                    log_warning(
                        "Constraint check: Parameter AsymptoticCasesPerInfectious changed from {:0.4f} to {:d} ",
                        this->get<AsymptoticCasesPerInfectious>()[i], 0);
                    this->get<AsymptoticCasesPerInfectious>()[i] = 0;
                }

                if (this->get<RiskOfInfectionFromSympomatic>()[i] < 0.0 ||
                    this->get<RiskOfInfectionFromSympomatic>()[i] > 1.0) {
                    log_warning(
                        "Constraint check: Parameter RiskOfInfectionFromSympomatic changed from {:0.4f} to {:d}",
                        this->get<RiskOfInfectionFromSympomatic>()[i], 0);
                    this->get<RiskOfInfectionFromSympomatic>()[i] = 0;
                }

                if (this->get<HospitalizedCasesPerInfectious>()[i] < 0.0 ||
                    this->get<HospitalizedCasesPerInfectious>()[i] > 1.0) {
                    log_warning(
                        "Constraint check: Parameter HospitalizedCasesPerInfectious changed from {:0.4f} to {:d}",
                        this->get<HospitalizedCasesPerInfectious>()[i], 0);
                    this->get<HospitalizedCasesPerInfectious>()[i] = 0;
                }

                if (this->get<ICUCasesPerHospitalized>()[i] < 0.0 || this->get<ICUCasesPerHospitalized>()[i] > 1.0) {
                    log_warning("Constraint check: Parameter ICUCasesPerHospitalized changed from {:0.4f} to {:d}",
                                this->get<ICUCasesPerHospitalized>()[i], 0);
                    this->get<ICUCasesPerHospitalized>()[i] = 0;
                }

                if (this->get<DeathsPerHospitalized>()[i] < 0.0 || this->get<DeathsPerHospitalized>()[i] > 1.0) {
                    log_warning("Constraint check: Parameter DeathsPerHospitalized changed from {:0.4f} to {:d}",
                                this->get<DeathsPerHospitalized>()[i], 0);
                    this->get<DeathsPerHospitalized>()[i] = 0;
                }
            }
        }

        /**
     * @brief checks whether all Parameters satisfy their corresponding constraints and throws errors, if they do not
     */
        void check_constraints() const
        {
            if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
                log_warning("Constraint check: Parameter m_seasonality smaller {:d} or larger {:d}", 0, 0.5);
            }

            if (this->get<ICUCapacity>() < 0.0) {
                log_warning("Constraint check: Parameter m_icu_capacity smaller {:d}", 0);
            }

            for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {

                if (this->get<IncubationTime>()[i] < 2.0) {
                    log_error("Constraint check: Parameter IncubationTime {:.4f} smaller {:.4f}",
                              this->get<IncubationTime>()[i], 2.0);
                }

                if (2 * this->get<SerialInterval>()[i] < this->get<IncubationTime>()[i] + 1.0) {
                    log_error("Constraint check: Parameter SerialInterval {:.4f} smaller {:.4f}",
                              this->get<SerialInterval>()[i], 0.5 * this->get<IncubationTime>()[i] + 0.5);
                }
                else if (this->get<SerialInterval>()[i] > this->get<IncubationTime>()[i] - 0.5) {
                    log_error("Constraint check: Parameter SerialInterval {:.4f} smaller {:.4f}",
                              this->get<SerialInterval>()[i], this->get<IncubationTime>()[i] - 0.5);
                }

                if (this->get<InfectiousTimeMild>()[i] < 1.0) {
                    log_error("Constraint check: Parameter InfectiousTimeMild {:.4f} smaller {:.4f}",
                              this->get<InfectiousTimeMild>()[i], 1.0);
                }

                if (this->get<HospitalizedToHomeTime>()[i] < 1.0) {
                    log_error("Constraint check: Parameter HospitalizedToHomeTime {:.4f} smaller {:.4f}",
                              this->get<HospitalizedToHomeTime>()[i], 1.0);
                }

                if (this->get<HomeToHospitalizedTime>()[i] < 1.0) {
                    log_error("Constraint check: Parameter HomeToHospitalizedTime {:.4f} smaller {:.4f}",
                              this->get<HomeToHospitalizedTime>()[i], 1.0);
                }

                if (this->get<HospitalizedToICUTime>()[i] < 1.0) {
                    log_error("Constraint check: Parameter HospitalizedToICUTime {:.4f} smaller {:.4f}",
                              this->get<HospitalizedToICUTime>()[i], 1.0);
                }

                if (this->get<ICUToHomeTime>()[i] < 1.0) {
                    log_error("Constraint check: Parameter ICUToHomeTime {:.4f} smaller {:.4f}",
                              this->get<ICUToHomeTime>()[i], 1.0);
                }

                if (abs(this->get<InfectiousTimeAsymptomatic>()[i] -
                        1.0 / (0.5 / (this->get<IncubationTime>()[i] - this->get<SerialInterval>()[i])) +
                        0.5 * this->get<InfectiousTimeMild>()[i]) > 1e-12) {
                    log_warning("Constraint check: Parameter InfectiousTimeAsymptomatic not set as fully dependent on "
                                "tinc, tserint and tinfmild, as proposed by "
                                "https://www.medrxiv.org/content/10.1101/2020.04.04.20053637v1.");
                }

                if (this->get<ICUToDeathTime>()[i] < 1.0) {
                    log_error("Constraint check: Parameter ICUToDeathTime {:.4f} smaller {:.4f}",
                              this->get<ICUToDeathTime>()[i], 1.0);
                }

                if (this->get<InfectionProbabilityFromContact>()[i] < 0.0) {
                    log_warning("Constraint check: Parameter InfectionProbabilityFromContact smaller {:d}", 0);
                }

                if (this->get<RelativeCarrierInfectability>()[i] < 0.0) {
                    log_warning("Constraint check: Parameter RelativeCarrierInfectability smaller {:d}", 0);
                }

                if (this->get<AsymptoticCasesPerInfectious>()[i] < 0.0 ||
                    this->get<AsymptoticCasesPerInfectious>()[i] > 1.0) {
                    log_warning("Constraint check: Parameter AsymptoticCasesPerInfectious smaller {:d} or larger {:d}",
                                0, 1);
                }

                if (this->get<RiskOfInfectionFromSympomatic>()[i] < 0.0 ||
                    this->get<RiskOfInfectionFromSympomatic>()[i] > 1.0) {
                    log_warning("Constraint check: Parameter RiskOfInfectionFromSympomatic smaller {:d} or larger {:d}",
                                0, 1);
                }

                if (this->get<HospitalizedCasesPerInfectious>()[i] < 0.0 ||
                    this->get<HospitalizedCasesPerInfectious>()[i] > 1.0) {
                    log_warning(
                        "Constraint check: Parameter HospitalizedCasesPerInfectious smaller {:d} or larger {:d}", 0, 1);
                }

                if (this->get<ICUCasesPerHospitalized>()[i] < 0.0 || this->get<ICUCasesPerHospitalized>()[i] > 1.0) {
                    log_warning("Constraint check: Parameter ICUCasesPerHospitalized smaller {:d} or larger {:d}", 0,
                                1);
                }

                if (this->get<DeathsPerHospitalized>()[i] < 0.0 || this->get<DeathsPerHospitalized>()[i] > 1.0) {
                    log_warning("Constraint check: Parameter DeathsPerHospitalized smaller {:d} or larger {:d}", 0, 1);
                }
            }
        }

    private:
        SecirParams(SecirParamsBase&& base)
            : SecirParamsBase(std::move(base))
            , m_num_groups(get<ContactPatterns>().get_cont_freq_mat().get_num_groups())
        {
        }

    public:
        /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
        template <class IOContext>
        static IOResult<SecirParams> deserialize(IOContext& io)
        {
            BOOST_OUTCOME_TRY(base, SecirParamsBase::deserialize(io));
            return success(SecirParams(std::move(base)));
        }

    private:
        AgeGroup m_num_groups;
        double m_commuter_nondetection    = 0.0;
        double m_start_commuter_detection = 0.0;
        double m_end_commuter_detection   = 0.0;
    };

    /**
 * @brief WIP !! TO DO: returns the actual, approximated reproduction rate 
 */
    //double get_reprod_rate(SecirParams const& params, double t, std::vector<double> const& yt);

} // namespace vaccinated
} // namespace mio

#endif // SECIR_H
