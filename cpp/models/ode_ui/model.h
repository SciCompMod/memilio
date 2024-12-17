/* 
* Copyright (C) 2020-2024 MEmilio
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

#ifndef ODEUI_MODEL_H
#define ODEUI_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_ui/infection_state.h"
#include "ode_ui/parameters.h"



namespace mio
{
namespace oui
{

/********************
 * define the model *
 ********************/

template <typename FP = ScalarType>
class Model
    : public mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>>
{
    using Base =
        mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;

    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    Model(int num_agegroups)
        : Base(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_derivatives(Eigen::Ref<const Vector<FP>> pop, Eigen::Ref<const Vector<FP>> y, FP t,
                         Eigen::Ref<Vector<FP>> dydt) const override
    {
        auto params                              = this->parameters;
        FP step_size                           = 0.1;
        AgeGroup n_agegroups                     = params.get_num_groups();
        ContactMatrixGroup const& contact_matrix = params.template get<ContactPatterns<FP>>();
        Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic> const& recovery_matrix_v1 = params.template get<RecoveryPatternsV1<FP>>();
        Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic> const& recovery_matrix_v2 = params.template get<RecoveryPatternsV2<FP>>();
        Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic> const& recovery_matrix_v3 = params.template get<RecoveryPatternsV3<FP>>();
        const FP divN = 1.0 / populations.get_total();


        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            size_t I1i = this->populations.get_flat_index({i, InfectionState::InfectedV1});
            size_t I2i = this->populations.get_flat_index({i, InfectionState::InfectedV2});
            size_t I3i = this->populations.get_flat_index({i, InfectionState::InfectedV3});

            FP coeffStoI1 = contact_matrix.get_matrix_at(t)(0, 0) * 
                                        params.template get<TransmissionProbabilityOnContactV1<FP>>() * divN;

            FP coeffStoI2 = contact_matrix.get_matrix_at(t)(0, 0) * 
                                        params.template get<TransmissionProbabilityOnContactV2<FP>>() * divN;

            FP coeffStoI3 = contact_matrix.get_matrix_at(t)(0, 0) * 
                                        params.template get<TransmissionProbabilityOnContactV3<FP>>() * divN;

            FP I1_sum = 0;
            FP I2_sum = 0;
            FP I3_sum = 0;

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {

                size_t I1j = this->populations.get_flat_index({j, InfectionState::InfectedV1});
                size_t I2j = this->populations.get_flat_index({j, InfectionState::InfectedV2});
                size_t I3j = this->populations.get_flat_index({j, InfectionState::InfectedV3});

                I1_sum += pop[I1j];
                I2_sum += pop[I2j];
                I3_sum += pop[I3j];

                FP coeffStoR1 = recovery_matrix_v1(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV1<FP>>());
                FP coeffStoR2 = recovery_matrix_v2(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV2<FP>>());
                FP coeffStoR3 = recovery_matrix_v3(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV3<FP>>());

                FP is_I1 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);
                FP is_I2 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);
                FP is_I3 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);

                FP flow_is_i1 = coeffStoR1 * pop[I1j] + sqrt(coeffStoR1 * pop[I1j] / step_size) * is_I1;
                FP flow_is_i2 = coeffStoR2 * pop[I2j] + sqrt(coeffStoR2 * pop[I2j] / step_size) * is_I2;
                FP flow_is_i3 = coeffStoR3 * pop[I3j] + sqrt(coeffStoR3 * pop[I3j] / step_size) * is_I3;

                
                dydt[I1j] += -flow_is_i1; 
                dydt[I2j] += -flow_is_i2; 
                dydt[I3j] += -flow_is_i3; 
                dydt[Si] += flow_is_i1 + flow_is_i2 + flow_is_i3; 

            }
            FP si_V1 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);
            FP si_V2 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);
            FP si_V3 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);

            FP flow_si_v1 = coeffStoI1 * y[Si] * I1_sum + sqrt(coeffStoI1 * y[Si] * I1_sum / step_size) * si_V1;
            FP flow_si_v2 = coeffStoI2 * y[Si] * I2_sum + sqrt(coeffStoI2 * y[Si] * I2_sum / step_size) * si_V2;
            FP flow_si_v3 = coeffStoI3 * y[Si] * I3_sum + sqrt(coeffStoI3 * y[Si] * I3_sum / step_size) * si_V3;

            dydt[Si] += -flow_si_v1 - flow_si_v2 - flow_si_v3;
            dydt[I1i] += flow_si_v1;
            dydt[I2i] += flow_si_v2;
            dydt[I3i] += flow_si_v3;

            FP ii_V1V2 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);
            FP ii_V1V3 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);
            FP ii_V2V1 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);
            FP ii_V2V3 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);
            FP ii_V3V1 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);
            FP ii_V3V2 = mio::DistributionAdapter<std::normal_distribution<FP>>::get_instance()(rng, 0.0, 1.0);

            FP flow_ii_v1v2 = 0.01 * y[I1i] + sqrt(0.01 * y[I1i] / step_size) * ii_V1V2;
            FP flow_ii_v1v3 = 0.0001 * y[I1i] + sqrt(0.0001 * y[I1i] / step_size) * ii_V1V3;
            FP flow_ii_v2v1 = 0.01 * y[I2i] + sqrt(0.01 * y[I2i] / step_size) * ii_V2V1;
            FP flow_ii_v2v3 = 0.01 * y[I2i] + sqrt(0.01 * y[I2i] / step_size) * ii_V2V3;
            FP flow_ii_v3v1 = 0.0001 * y[I3i] + sqrt(0.0001 * y[I3i] / step_size) * ii_V3V1;
            FP flow_ii_v3v2 = 0.01 * y[I3i] + sqrt(0.01 * y[I3i] / step_size) * ii_V3V2;

            dydt[I1i] += -flow_ii_v1v2 - flow_ii_v1v3 + flow_ii_v2v1 + flow_ii_v3v1;
            dydt[I2i] += -flow_ii_v2v1 - flow_ii_v2v3 + flow_ii_v1v2 + flow_ii_v3v2;
            dydt[I3i] += -flow_ii_v3v1 - flow_ii_v3v1 + flow_ii_v1v3 + flow_ii_v2v3;

            /*std::cout << "Si: "  <<  y[Si] << "  " << dydt[Si] << " I1i: " << y[I1i] << "  " << dydt[I1i] << " I2i: " << y[I2i] 
                        << "  "  <<  dydt[I2i] << " I3i: " << y[I3i] << "  "  << dydt[I3i] << "\n";*/
            
        }
        /*std::cout << "\n";
        getchar();*/
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Model");
        obj.add_element("Parameters", this->parameters);
        obj.add_element("Populations", this->populations);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Model> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Model");
        auto par = obj.expect_element("Parameters", Tag<ParameterSet>{});
        auto pop = obj.expect_element("Populations", Tag<Populations>{});
        return apply(
            io,
            [](auto&& par_, auto&& pop_) {
                return Model{pop_, par_};
            },
            par, pop);
    }
    mutable RandomNumberGenerator rng;
};

} // namespace oui
} // namespace mio

#endif // ODEUI_MODEL_H
