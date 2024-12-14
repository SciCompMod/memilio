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

#ifndef ODEUI10_MODEL_H
#define ODEUI10_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "ode_ui10/infection_state.h"
#include "ode_ui10/parameters.h"



namespace mio
{
namespace oui10
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
        AgeGroup n_agegroups                     = params.get_num_groups();
        /*ContactMatrixGroup const& contact_matrix = params.template get<ContactPatterns<FP>>();*/
        ContactMatrixGroup const& recovery_matrix_v1 = params.template get<RecoveryPatternsV1<FP>>();
        ContactMatrixGroup const& recovery_matrix_v2 = params.template get<RecoveryPatternsV2<FP>>();
        ContactMatrixGroup const& recovery_matrix_v3 = params.template get<RecoveryPatternsV3<FP>>();
        ContactMatrixGroup const& recovery_matrix_v4 = params.template get<RecoveryPatternsV4<FP>>();
        ContactMatrixGroup const& recovery_matrix_v5 = params.template get<RecoveryPatternsV5<FP>>();
        ContactMatrixGroup const& recovery_matrix_v6 = params.template get<RecoveryPatternsV6<FP>>();
        ContactMatrixGroup const& recovery_matrix_v7 = params.template get<RecoveryPatternsV7<FP>>();
        ContactMatrixGroup const& recovery_matrix_v8 = params.template get<RecoveryPatternsV8<FP>>();
        ContactMatrixGroup const& recovery_matrix_v9 = params.template get<RecoveryPatternsV9<FP>>();
        ContactMatrixGroup const& recovery_matrix_v10 = params.template get<RecoveryPatternsV10<FP>>();
        const ScalarType divN = 1.0 / populations.get_total();

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            size_t I1i = this->populations.get_flat_index({i, InfectionState::InfectedV1});
            size_t I2i = this->populations.get_flat_index({i, InfectionState::InfectedV2});
            size_t I3i = this->populations.get_flat_index({i, InfectionState::InfectedV3});
            size_t I4i = this->populations.get_flat_index({i, InfectionState::InfectedV4});
            size_t I5i = this->populations.get_flat_index({i, InfectionState::InfectedV5});
            size_t I6i = this->populations.get_flat_index({i, InfectionState::InfectedV6});
            size_t I7i = this->populations.get_flat_index({i, InfectionState::InfectedV7});
            size_t I8i = this->populations.get_flat_index({i, InfectionState::InfectedV8});
            size_t I9i = this->populations.get_flat_index({i, InfectionState::InfectedV9});
            size_t I10i = this->populations.get_flat_index({i, InfectionState::InfectedV10});

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {

                size_t I1j = this->populations.get_flat_index({j, InfectionState::InfectedV1});
                size_t I2j = this->populations.get_flat_index({j, InfectionState::InfectedV2});
                size_t I3j = this->populations.get_flat_index({j, InfectionState::InfectedV3});
                size_t I4j = this->populations.get_flat_index({j, InfectionState::InfectedV4});
                size_t I5j = this->populations.get_flat_index({j, InfectionState::InfectedV5});
                size_t I6j = this->populations.get_flat_index({j, InfectionState::InfectedV6});
                size_t I7j = this->populations.get_flat_index({j, InfectionState::InfectedV7});
                size_t I8j = this->populations.get_flat_index({j, InfectionState::InfectedV8});
                size_t I9j = this->populations.get_flat_index({j, InfectionState::InfectedV9});
                size_t I10j = this->populations.get_flat_index({j, InfectionState::InfectedV10});

                ScalarType coeffStoI1 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV1<FP>>()[i] * divN;
                ScalarType coeffStoI2 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV2<FP>>()[i] * divN;
                ScalarType coeffStoI3 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV3<FP>>()[i] * divN;                
                ScalarType coeffStoI4 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV4<FP>>()[i] * divN;
                ScalarType coeffStoI5 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV5<FP>>()[i] * divN;                
                ScalarType coeffStoI6 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV6<FP>>()[i] * divN;
                ScalarType coeffStoI7 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV7<FP>>()[i] * divN;
                ScalarType coeffStoI8 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV8<FP>>()[i] * divN;
                ScalarType coeffStoI9 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV9<FP>>()[i] * divN;
                ScalarType coeffStoI10 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV10<FP>>()[i] * divN;


                dydt[Si] += y[Si] * (-coeffStoI1 * pop[I1j] - coeffStoI2 * pop[I2j] - coeffStoI3 * pop[I3j] 
                                     -coeffStoI4 * pop[I4j] - coeffStoI5 * pop[I5j] - coeffStoI6 * pop[I6j]
                                     -coeffStoI7 * pop[I7j] - coeffStoI8 * pop[I8j] - coeffStoI9 * pop[I9j]
                                     -coeffStoI10 * pop[I10j]);
                dydt[I1i] += coeffStoI1 * y[Si] * pop[I1j];
                dydt[I2i] += coeffStoI2 * y[Si] * pop[I2j];
                dydt[I3i] += coeffStoI3 * y[Si] * pop[I3j];
                dydt[I4i] += coeffStoI4 * y[Si] * pop[I4j];
                dydt[I5i] += coeffStoI5 * y[Si] * pop[I5j];
                dydt[I6i] += coeffStoI6 * y[Si] * pop[I6j];
                dydt[I7i] += coeffStoI7 * y[Si] * pop[I7j];
                dydt[I8i] += coeffStoI8 * y[Si] * pop[I8j];
                dydt[I9i] += coeffStoI9 * y[Si] * pop[I9j];
                dydt[I10i] += coeffStoI10 * y[Si] * pop[I10j];

                ScalarType coeffStoR1 = recovery_matrix_v1.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV1<FP>>()[i]);
                ScalarType coeffStoR2 = recovery_matrix_v2.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV2<FP>>()[i]);
                ScalarType coeffStoR3 = recovery_matrix_v3.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV3<FP>>()[i]);
                ScalarType coeffStoR4 = recovery_matrix_v4.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV4<FP>>()[i]);
                ScalarType coeffStoR5 = recovery_matrix_v5.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV5<FP>>()[i]);
                ScalarType coeffStoR6 = recovery_matrix_v6.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV6<FP>>()[i]);
                ScalarType coeffStoR7 = recovery_matrix_v7.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV7<FP>>()[i]);
                ScalarType coeffStoR8 = recovery_matrix_v8.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV8<FP>>()[i]);
                ScalarType coeffStoR9 = recovery_matrix_v9.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV9<FP>>()[i]);
                ScalarType coeffStoR10 = recovery_matrix_v10.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV10<FP>>()[i]);

                /*dydt[I1j] += -coeffStoR1 * pop[I1j]; 
                dydt[I2j] += -coeffStoR2 * pop[I2j]; 
                dydt[I3j] += -coeffStoR3 * pop[I3j]; 
                dydt[I4j] += -coeffStoR4 * pop[I4j];
                dydt[I5j] += -coeffStoR5 * pop[I5j];
                dydt[I6j] += -coeffStoR6 * pop[I6j];
                dydt[I7j] += -coeffStoR7 * pop[I7j];
                dydt[I8j] += -coeffStoR8 * pop[I8j];
                dydt[I9j] += -coeffStoR9 * pop[I9j];
                dydt[I10j] += -coeffStoR10 * pop[I10j];
                dydt[Si] += coeffStoR1 * pop[I1j] + coeffStoR2 * pop[I2j] + coeffStoR3 * pop[I3j] 
                            + coeffStoR4 * pop[I4j] + coeffStoR5 * pop[I5j] + coeffStoR6 * pop[I6j]
                            + coeffStoR7 * pop[I7j] + coeffStoR8 * pop[I8j] + coeffStoR9 * pop[I9j]
                            + coeffStoR10 * pop[I10j];*/

                dydt[I1j] += -coeffStoR1 * y[I1j]; 
                dydt[I2j] += -coeffStoR2 * y[I2j]; 
                dydt[I3j] += -coeffStoR3 * y[I3j]; 
                dydt[I4j] += -coeffStoR4 * y[I4j];
                dydt[I5j] += -coeffStoR5 * y[I5j];
                dydt[I6j] += -coeffStoR6 * y[I6j];
                dydt[I7j] += -coeffStoR7 * y[I7j];
                dydt[I8j] += -coeffStoR8 * y[I8j];
                dydt[I9j] += -coeffStoR9 * y[I9j];
                dydt[I10j] += -coeffStoR10 * y[I10j];
                dydt[Si] += coeffStoR1 * y[I1j] + coeffStoR2 * y[I2j] + coeffStoR3 * y[I3j] 
                            + coeffStoR4 * y[I4j] + coeffStoR5 * y[I5j] + coeffStoR6 * y[I6j]
                            + coeffStoR7 * y[I7j] + coeffStoR8 * y[I8j] + coeffStoR9 * y[I9j]
                            + coeffStoR10 * y[I10j];
            }
            ScalarType mut_p = 0.001;
            ScalarType mut_p2 = mut_p * mut_p;

            /*dydt[I1i] += (-mut_p - mut_p2) * y[I1i] + mut_p * y[I2i] + mut_p2 * y[I3i];
            dydt[I2i] += (-2.0 * mut_p - mut_p2) * y[I2i] + mut_p * (y[I1i] + y[I3i]) + mut_p2 * y[I4i];
            dydt[I3i] += (-mut_p - mut_p2) * 2.0 * y[I3i] + mut_p * (y[I2i] + y[I4i]) + mut_p2 * (y[I1i] + y[I5i]);
            dydt[I4i] += (-mut_p - mut_p2) * 2.0 * y[I4i] + mut_p * (y[I3i] + y[I5i]) + mut_p2 * (y[I2i] + y[I6i]);
            dydt[I5i] += (-mut_p - mut_p2) * 2.0 * y[I5i] + mut_p * (y[I4i] + y[I6i]) + mut_p2 * (y[I3i] + y[I7i]);
            dydt[I6i] += (-mut_p - mut_p2) * 2.0 * y[I6i] + mut_p * (y[I5i] + y[I7i]) + mut_p2 * (y[I4i] + y[I8i]);
            dydt[I7i] += (-mut_p - mut_p2) * 2.0 * y[I7i] + mut_p * (y[I6i] + y[I8i]) + mut_p2 * (y[I5i] + y[I9i]);
            dydt[I8i] += (-mut_p - mut_p2) * 2.0 * y[I8i] + mut_p * (y[I7i] + y[I9i]) + mut_p2 * (y[I6i] + y[I10i]);
            dydt[I9i] += (-2.0 * mut_p - mut_p2) * y[I9i] + mut_p * (y[I8i] + y[I10i]) + mut_p2 * y[I7i];
            dydt[I10i] += (-mut_p - mut_p2) * y[I10i] + mut_p * y[I9i] + mut_p2 * y[I8i];*/

            dydt[I1i] += (-mut_p - mut_p2) * pop[I1i] + mut_p * pop[I2i] + mut_p2 * pop[I3i];
            dydt[I2i] += (-2.0 * mut_p - mut_p2) * pop[I2i] + mut_p * (pop[I1i] + pop[I3i]) + mut_p2 * pop[I4i];
            dydt[I3i] += (-mut_p - mut_p2) * 2.0 * pop[I3i] + mut_p * (pop[I2i] + pop[I4i]) + mut_p2 * (pop[I1i] + pop[I5i]);
            dydt[I4i] += (-mut_p - mut_p2) * 2.0 * pop[I4i] + mut_p * (pop[I3i] + pop[I5i]) + mut_p2 * (pop[I2i] + pop[I6i]);
            dydt[I5i] += (-mut_p - mut_p2) * 2.0 * pop[I5i] + mut_p * (pop[I4i] + pop[I6i]) + mut_p2 * (pop[I3i] + pop[I7i]);
            dydt[I6i] += (-mut_p - mut_p2) * 2.0 * pop[I6i] + mut_p * (pop[I5i] + pop[I7i]) + mut_p2 * (pop[I4i] + pop[I8i]);
            dydt[I7i] += (-mut_p - mut_p2) * 2.0 * pop[I7i] + mut_p * (pop[I6i] + pop[I8i]) + mut_p2 * (pop[I5i] + pop[I9i]);
            dydt[I8i] += (-mut_p - mut_p2) * 2.0 * pop[I8i] + mut_p * (pop[I7i] + pop[I9i]) + mut_p2 * (pop[I6i] + pop[I10i]);
            dydt[I9i] += (-2.0 * mut_p - mut_p2) * pop[I9i] + mut_p * (pop[I8i] + pop[I10i]) + mut_p2 * pop[I7i];
            dydt[I10i] += (-mut_p - mut_p2) * pop[I10i] + mut_p * pop[I9i] + mut_p2 * pop[I8i];

            /*dydt[I2i] += (-mut_p) * y[I2i] + mut_p * y[I1i];
            dydt[I3i] += (-mut_p - mut_p2) * y[I3i] + mut_p * y[I2i] + mut_p2 * y[I1i];
            dydt[I4i] += (-mut_p - mut_p2) * y[I4i] + mut_p * y[I3i] + mut_p2 * y[I2i];
            dydt[I5i] += (-mut_p - mut_p2) * y[I5i] + mut_p * y[I4i] + mut_p2 * y[I3i];
            dydt[I6i] += (-mut_p - mut_p2) * y[I6i] + mut_p * y[I5i] + mut_p2 * y[I4i];
            dydt[I7i] += (-mut_p - mut_p2) * y[I7i] + mut_p * y[I6i] + mut_p2 * y[I5i];
            dydt[I8i] += (-mut_p - mut_p2) * y[I8i] + mut_p * y[I7i] + mut_p2 * y[I6i];
            dydt[I9i] += (-mut_p - mut_p2) * y[I9i] + mut_p * y[I8i] + mut_p2 * y[I7i];
            dydt[I10i] += (-mut_p - mut_p2) * y[I10i] + mut_p * y[I9i] + mut_p2 * y[I8i];*/

            /*dydt[I2i] += (-mut_p) * pop[I2i] + mut_p * pop[I1i];
            dydt[I3i] += (-mut_p - mut_p2) * pop[I3i] + mut_p * pop[I2i] + mut_p2 * pop[I1i];
            dydt[I4i] += (-mut_p - mut_p2) * pop[I4i] + mut_p * pop[I3i] + mut_p2 * pop[I2i];
            dydt[I5i] += (-mut_p - mut_p2) * pop[I5i] + mut_p * pop[I4i] + mut_p2 * pop[I3i];
            dydt[I6i] += (-mut_p - mut_p2) * pop[I6i] + mut_p * pop[I5i] + mut_p2 * pop[I4i];
            dydt[I7i] += (-mut_p - mut_p2) * pop[I7i] + mut_p * pop[I6i] + mut_p2 * pop[I5i];
            dydt[I8i] += (-mut_p - mut_p2) * pop[I8i] + mut_p * pop[I7i] + mut_p2 * pop[I6i];
            dydt[I9i] += (-mut_p - mut_p2) * pop[I9i] + mut_p * pop[I8i] + mut_p2 * pop[I7i];
            dydt[I10i] += (-mut_p - mut_p2) * pop[I10i] + mut_p * pop[I9i] + mut_p2 * pop[I8i];*/
            
        }
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
};

} // namespace oui
} // namespace mio

#endif // ODEUI_MODEL_H
