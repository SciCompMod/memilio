/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "secir_vaccine/parameter_space.h"
#include "memilio/utils/parameter_distributions.h"
#include "secir_vaccine/secir.h"

namespace mio
{
namespace vaccinated
{

    void set_params_distributions_normal(SecirModel& model, double t0, double tmax, double dev_rel)
    {
        auto set_distribution = [dev_rel](UncertainValue& v, double min_val = 0.001) {
            v.set_distribution(ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * v),
                                                           std::max(min_val, (1 + dev_rel * 2.6) * v),
                                                           std::max(min_val, double(v)), dev_rel * v));
        };

        set_distribution(model.parameters.get<Seasonality>(), 0.0);
        set_distribution(model.parameters.get<ICUCapacity>());
        set_distribution(model.parameters.get<TestAndTraceCapacity>());

        // populations
        for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
            for (auto j = Index<InfectionState>(0); j < Index<InfectionState>(InfectionState::Count); j++) {

                // don't touch S and D
                if (j == InfectionState::Susceptible || j == InfectionState::Dead) {
                    continue;
                }

                // variably sized groups
                set_distribution(model.populations[{i, j}], 0.0);
            }
        }

        // times
        for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {

            set_distribution(model.parameters.get<IncubationTime>()[i]);
            set_distribution(model.parameters.get<SerialInterval>()[i]);
            set_distribution(model.parameters.get<InfectiousTimeMild>()[i]);
            set_distribution(model.parameters.get<HospitalizedToHomeTime>()[i]);
            set_distribution(model.parameters.get<HomeToHospitalizedTime>()[i]);
            set_distribution(model.parameters.get<InfectiousTimeAsymptomatic>()[i]);
            set_distribution(model.parameters.get<HospitalizedToICUTime>()[i]);
            set_distribution(model.parameters.get<ICUToHomeTime>()[i]);
            set_distribution(model.parameters.get<ICUToDeathTime>()[i]);

            set_distribution(model.parameters.get<InfectionProbabilityFromContact>()[i]);
            set_distribution(model.parameters.get<RelativeCarrierInfectability>()[i]);
            set_distribution(model.parameters.get<AsymptoticCasesPerInfectious>()[i]);
            set_distribution(model.parameters.get<RiskOfInfectionFromSympomatic>()[i]);
            set_distribution(model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[i]);
            set_distribution(model.parameters.get<DeathsPerHospitalized>()[i]);
            set_distribution(model.parameters.get<HospitalizedCasesPerInfectious>()[i]);
            set_distribution(model.parameters.get<ICUCasesPerHospitalized>()[i]);
        }

        // dampings
        auto matrices = std::vector<size_t>();
        for (size_t i = 0; i < model.parameters.get<ContactPatterns>().get_cont_freq_mat().get_num_matrices(); ++i) {
            matrices.push_back(i);
        }
        auto groups = Eigen::VectorXd::Constant(Eigen::Index(model.parameters.get_num_groups().get()), 1.0);
        model.parameters.get<ContactPatterns>().get_dampings().emplace_back(
            mio::UncertainValue(0.5), mio::DampingLevel(0), mio::DampingType(0),
            mio::SimulationTime(t0 + (tmax - t0) / 2), matrices, groups);
        set_distribution(model.parameters.get<ContactPatterns>().get_dampings()[0].get_value(), 0.0);
    }

    void draw_sample_demographics(SecirModel& model)
    {
        model.parameters.get<ICUCapacity>().draw_sample();
        model.parameters.get<TestAndTraceCapacity>().draw_sample();

        for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
            double group_total = model.populations.get_group_total(i);

            model.populations[{i, InfectionState::Exposed}].draw_sample();
            model.populations[{i, InfectionState::Carrier}].draw_sample();
            model.populations[{i, InfectionState::Infected}].draw_sample();
            model.populations[{i, InfectionState::Hospitalized}].draw_sample();
            model.populations[{i, InfectionState::ICU}].draw_sample();
            model.populations[{i, InfectionState::Recovered}].draw_sample();

            // no sampling for dead and total numbers
            // [...]

            model.populations[{i, InfectionStateV::Susceptible}] = 0;
            double group_total_dummy                             = model.populations.get_group_total(i);
            if (group_total_dummy < group_total) {
                model.populations.set_difference_from_group_total<AgeGroup>({i, epi::InfectionStateV::Susceptible},
                                                                            group_total);
            }
            else {
                double diff = group_total_dummy - group_total;
                model.populations[{i, InfectionStateV::Recovered}] =
                    model.populations[{i, InfectionStateV::Recovered}] - diff;
                assert(std::abs(group_total - model.populations.get_group_total(i)) < 1e-10);
            }

            model.populations.set_difference_from_group_total<AgeGroup>({i, InfectionStateV::Susceptible},
                                                                        model.populations.get_group_total(i));
        }
    }

    void draw_sample_infection(SecirModel& model)
    {
        model.parameters.get<Seasonality>().draw_sample();

        //not age dependent
        model.parameters.get<IncubationTime>()[AgeGroup(0)].draw_sample();
        model.parameters.get<SerialInterval>()[AgeGroup(0)].draw_sample();
        model.parameters.get<InfectiousTimeMild>()[AgeGroup(0)].draw_sample();
        model.parameters.get<HospitalizedToICUTime>()[AgeGroup(0)].draw_sample();
        model.parameters.get<RelativeCarrierInfectability>()[AgeGroup(0)].draw_sample();
        model.parameters.get<RiskOfInfectionFromSympomatic>()[AgeGroup(0)].draw_sample();
        model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[AgeGroup(0)].draw_sample();

        model.parameters.get<ReducVaccExp>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducImmuneExp>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducExpInf>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducImmuneExpInf>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducInfHosp>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducImmuneInfHosp>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducTime>()[AgeGroup(0)].draw_sample();

        for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
            //not age dependent
            model.parameters.get<IncubationTime>()[i]     = model.parameters.get<IncubationTime>()[AgeGroup(0)];
            model.parameters.get<SerialInterval>()[i]     = model.parameters.get<SerialInterval>()[AgeGroup(0)];
            model.parameters.get<InfectiousTimeMild>()[i] = model.parameters.get<InfectiousTimeMild>()[AgeGroup(0)];
            model.parameters.get<HospitalizedToICUTime>()[i] =
                model.parameters.get<HospitalizedToICUTime>()[AgeGroup(0)];
            model.parameters.get<RelativeCarrierInfectability>()[i] =
                model.parameters.get<RelativeCarrierInfectability>()[AgeGroup(0)];
            model.parameters.get<RiskOfInfectionFromSympomatic>()[i] =
                model.parameters.get<RiskOfInfectionFromSympomatic>()[AgeGroup(0)];
            model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[i] =
                model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[AgeGroup(0)];

            model.parameters.get<ReducVaccExp>()[i]       = model.parameters.get<ReducVaccExp>()[AgeGroup(0)];
            model.parameters.get<ReducImmuneExp>()[i]     = model.parameters.get<ReducImmuneExp>()[AgeGroup(0)];
            model.parameters.get<ReducExpInf>()[i]        = model.parameters.get<ReducExpInf>()[AgeGroup(0)];
            model.parameters.get<ReducImmuneExpInf>()[i]  = model.parameters.get<ReducImmuneExpInf>()[AgeGroup(0)];
            model.parameters.get<ReducInfHosp>()[i]       = model.parameters.get<ReducInfHosp>()[AgeGroup(0)];
            model.parameters.get<ReducImmuneInfHosp>()[i] = model.parameters.get<ReducImmuneInfHosp>()[AgeGroup(0)];
            model.parameters.get<ReducTime>()[i]          = model.parameters.get<ReducTime>()[AgeGroup(0)];

            //age dependent
            model.parameters.get<HospitalizedToHomeTime>()[i].draw_sample(); // here: home=recovered
            model.parameters.get<HomeToHospitalizedTime>()[i].draw_sample(); // here: home=infectious
            model.parameters.get<InfectiousTimeAsymptomatic>()[i].draw_sample();
            model.parameters.get<ICUToDeathTime>()[i].draw_sample();
            model.parameters.get<ICUToHomeTime>()[i].draw_sample();

            model.parameters.get<InfectionProbabilityFromContact>()[i].draw_sample();
            model.parameters.get<AsymptoticCasesPerInfectious>()[i].draw_sample();
            model.parameters.get<DeathsPerHospitalized>()[i].draw_sample();
            model.parameters.get<HospitalizedCasesPerInfectious>()[i].draw_sample();
            model.parameters.get<ICUCasesPerHospitalized>()[i].draw_sample();
        }
    }

    void draw_sample(SecirModel& model)
    {
        draw_sample_infection(model);
        draw_sample_demographics(model);
        model.parameters.get<ContactPatterns>().draw_sample();
        model.apply_constraints();
    }

} // namespace vaccinated
} // namespace mio
