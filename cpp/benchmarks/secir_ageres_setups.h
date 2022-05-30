/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Rene Schmieding
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
#ifndef SECIR_AGERES_SETUPS_H_
#define SECIR_AGERES_SETUPS_H_

#include "memilio/compartments/simulation.h"
#include "memilio/mobility/mobility.h"
#include "models/secir/secir.h"
#include "models/secir/parameter_space.h"
#include "models/secir/parameter_studies.h"

namespace mio
{
namespace benchmark
{
    namespace detail
    {
        /**
         * @brief Helper function to create a secir model with consistent setup for use in benchmarking.
         */
        mio::SecirModel make_model(int num)
        {
            double tinc    = 5.2, // R_2^(-1)+R_3^(-1)
                tinfmild   = 6, // 4-14  (=R4^(-1))
                tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
                thosp2home = 12, // 7-16 (=R5^(-1))
                thome2hosp = 5, // 2.5-7 (=R6^(-1))
                thosp2icu  = 2, // 1-3.5 (=R7^(-1))
                ticu2home  = 8, // 5-16 (=R8^(-1))
                // tinfasy    = 6.2, // (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
                ticu2death = 5; // 3.5-7 (=R5^(-1))

            double cont_freq = 10, // see Polymod study
                inf_prob = 0.05, carr_infec = 0.67,
                   alpha = 0.09, // 0.01-0.16
                beta     = 0.25, // 0.05-0.5
                delta    = 0.3, // 0.15-0.77
                rho      = 0.2, // 0.1-0.35
                theta    = 0.25; // 0.15-0.4

            double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20,
                   nb_icu_t0 = 10, nb_rec_t0 = 10, nb_dead_t0 = 0;

            mio::SecirModel model(num);
            auto nb_groups = model.parameters.get_num_groups();
            double fact    = 1.0 / (double)(size_t)nb_groups;

            auto& params = model.parameters;

            params.set<mio::ICUCapacity>(std::numeric_limits<double>::max());
            params.set<mio::StartDay>(0);
            params.set<mio::Seasonality>(0);

            for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
                params.get<mio::IncubationTime>()[i]         = tinc;
                params.get<mio::InfectiousTimeMild>()[i]     = tinfmild;
                params.get<mio::SerialInterval>()[i]         = tserint;
                params.get<mio::HospitalizedToHomeTime>()[i] = thosp2home;
                params.get<mio::HomeToHospitalizedTime>()[i] = thome2hosp;
                params.get<mio::HospitalizedToICUTime>()[i]  = thosp2icu;
                params.get<mio::ICUToHomeTime>()[i]          = ticu2home;
                params.get<mio::ICUToDeathTime>()[i]         = ticu2death;

                model.populations[{i, mio::InfectionState::Exposed}]      = fact * nb_exp_t0;
                model.populations[{i, mio::InfectionState::Carrier}]      = fact * nb_car_t0;
                model.populations[{i, mio::InfectionState::Infected}]     = fact * nb_inf_t0;
                model.populations[{i, mio::InfectionState::Hospitalized}] = fact * nb_hosp_t0;
                model.populations[{i, mio::InfectionState::ICU}]          = fact * nb_icu_t0;
                model.populations[{i, mio::InfectionState::Recovered}]    = fact * nb_rec_t0;
                model.populations[{i, mio::InfectionState::Dead}]         = fact * nb_dead_t0;
                model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::InfectionState::Susceptible},
                                                                                 fact * nb_total_t0);

                params.get<mio::InfectionProbabilityFromContact>()[i] = inf_prob;
                params.get<mio::RelativeCarrierInfectability>()[i]    = carr_infec;
                params.get<mio::AsymptomaticCasesPerInfectious>()[i]    = alpha;
                params.get<mio::RiskOfInfectionFromSymptomatic>()[i]   = beta;
                params.get<mio::HospitalizedCasesPerInfectious>()[i]  = rho;
                params.get<mio::ICUCasesPerHospitalized>()[i]         = theta;
                params.get<mio::DeathsPerICU>()[i]                    = delta;
            }

            mio::ContactMatrixGroup& contact_matrix = params.get<mio::ContactPatterns>();
            contact_matrix[0] =
                mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                                       mio::SimulationTime(30.));

            model.apply_constraints();

            return model;
        }
    } // namespace detail

    namespace model
    {
        /**
         * @brief Secir model with consistent setup for use in benchmarking.
         */
        mio::SecirModel SecirAgeres(size_t num_agegroups)
        {
            mio::SecirModel model = mio::benchmark::detail::make_model(num_agegroups);

            auto nb_groups   = model.parameters.get_num_groups();
            double cont_freq = 10, fact = 1.0 / (double)(size_t)nb_groups;
            mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::ContactPatterns>();
            contact_matrix[0] =
                mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                                       mio::SimulationTime(30.));

            return model;
        }
        /**
         * @brief Secir model with consistent setup for use in benchmarking with added dampings.
         */
        mio::SecirModel SecirAgeresDampings(size_t num_agegroups)
        {
            mio::SecirModel model = mio::benchmark::detail::make_model(num_agegroups);

            auto nb_groups   = model.parameters.get_num_groups();
            double cont_freq = 10, fact = 1.0 / (double)(size_t)nb_groups;
            mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::ContactPatterns>();
            contact_matrix[0] =
                mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                                       mio::SimulationTime(25.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.3),
                                       mio::SimulationTime(40.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.8),
                                       mio::SimulationTime(60.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.5),
                                       mio::SimulationTime(75.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 1.0),
                                       mio::SimulationTime(95.));

            return model;
        }
        /**
         * @brief Secir model with consistent setup for use in benchmarking with added dampings.
         * Dampings are set up to challenge the integrator, not to be realistic.
         */
        mio::SecirModel SecirAgeresAbsurdDampings(size_t num_agegroups)
        {
            mio::SecirModel model = mio::benchmark::detail::make_model(num_agegroups);

            auto nb_groups   = model.parameters.get_num_groups();
            double cont_freq = 10, fact = 1.0 / (double)(size_t)nb_groups;
            mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::ContactPatterns>();
            contact_matrix[0] =
                mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.8),
                                       mio::SimulationTime(10.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.5),
                                       mio::SimulationTime(11.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.2),
                                       mio::SimulationTime(12.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.1),
                                       mio::SimulationTime(13.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.9),
                                       mio::SimulationTime(30.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.2),
                                       mio::SimulationTime(30.5));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                                       mio::SimulationTime(31.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.2),
                                       mio::SimulationTime(31.5));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.8),
                                       mio::SimulationTime(32.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.1),
                                       mio::SimulationTime(40.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.001),
                                       mio::SimulationTime(44.));
            contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.9),
                                       mio::SimulationTime(46.));

            return model;
        }
    } // namespace model
} // namespace benchmark

} // namespace mio

#endif