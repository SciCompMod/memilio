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
        AgeGroup n_agegroups                     = params.get_num_groups();
        /*ContactMatrixGroup const& contact_matrix = params.template get<ContactPatterns<FP>>();*/
        ContactMatrixGroup const& recovery_matrix_v1 = params.template get<RecoveryPatternsV1<FP>>();
        ContactMatrixGroup const& recovery_matrix_v2 = params.template get<RecoveryPatternsV2<FP>>();
        ContactMatrixGroup const& recovery_matrix_v3 = params.template get<RecoveryPatternsV3<FP>>();
        const ScalarType divN = 1.0 / populations.get_total();

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            size_t I1i = this->populations.get_flat_index({i, InfectionState::InfectedV1});
            size_t I2i = this->populations.get_flat_index({i, InfectionState::InfectedV2});
            size_t I3i = this->populations.get_flat_index({i, InfectionState::InfectedV3});

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {

                size_t I1j = this->populations.get_flat_index({j, InfectionState::InfectedV1});
                size_t I2j = this->populations.get_flat_index({j, InfectionState::InfectedV2});
                size_t I3j = this->populations.get_flat_index({j, InfectionState::InfectedV3});

                ScalarType coeffStoI1 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV1<FP>>()[i] * divN;

                

                ScalarType coeffStoI2 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV2<FP>>()[i] * divN;

                ScalarType coeffStoI3 = /*contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * */
                                       params.template get<TransmissionProbabilityOnContactV3<FP>>()[i] * divN;

                
                dydt[Si] += -coeffStoI1 * y[Si] * pop[I1j] - coeffStoI2 * y[Si] * pop[I2j] - coeffStoI3 * y[Si] * pop[I3j];
                dydt[I1i] += coeffStoI1 * y[Si] * pop[I1j];
                dydt[I2i] += coeffStoI2 * y[Si] * pop[I2j];
                dydt[I3i] += coeffStoI3 * y[Si] * pop[I3j];


                ScalarType coeffStoR1 = recovery_matrix_v1.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV1<FP>>()[i]);
                ScalarType coeffStoR2 = recovery_matrix_v2.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV2<FP>>()[i]);
                ScalarType coeffStoR3 = recovery_matrix_v3.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                       static_cast<Eigen::Index>((size_t)j)) * 
                                        (1.0 / params.template get<TimeInfectedV3<FP>>()[i]);

                dydt[I1j] += -coeffStoR1 * pop[I1j]; 
                dydt[I2j] += -coeffStoR2 * pop[I2j]; 
                dydt[I3j] += -coeffStoR3 * pop[I3j]; 
                dydt[Si] += coeffStoR1 * pop[I1j] + coeffStoR2 * pop[I2j] + coeffStoR3 * pop[I3j]; 

            }
            dydt[I1i] += -0.01 * y[I1i] - 0.0001 * y[I1i] + 0.01 * y[I2i] + 0.0001 * y[I3i];
            dydt[I2i] += -0.01 * y[I2i] - 0.01 * y[I2i] + 0.01 * y[I1i] + 0.01 * y[I3i];
            dydt[I3i] += -0.01 * y[I3i] - 0.0001 * y[I3i] + 0.0001 * y[I1i] + 0.01 * y[I2i];
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
