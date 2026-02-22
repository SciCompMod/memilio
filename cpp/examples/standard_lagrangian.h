#ifndef MIO_Standard_Lagrangian
#define MIO_Standard_Lagrangian

#include "state_estimators.h"

namespace mio
{
namespace examples
{

template <typename FP = ScalarType>
class ModelStandardLagrangian : public ModelExplicit<FP>
{
    using Base = ModelExplicit<FP>;

public:
    using Base::Base; // Inherit constructors

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const auto& params = this->parameters;

        // E->I and I->R flows (same as optimized version)
        for (mio::AgeGroup g(0); g < mio::AgeGroup(params.get_num_groups()); ++g) {
            {
                CommuterType commuter_type = CommuterType::NonCommuter;
                const size_t E_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::E});
                const size_t I_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::I});

                flows[Base::template get_flat_flow_index<InfectionStateExplicit::E, InfectionStateExplicit::I>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeExposed<FP>>()[g]) * y[E_idx];

                flows[Base::template get_flat_flow_index<InfectionStateExplicit::I, InfectionStateExplicit::R>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeInfected<FP>>()[g]) * y[I_idx];
            }

            for (int c = 0; c < this->m_num_commuter_groups; ++c) {
                CommuterType commuter_type =
                    static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c);

                const size_t E_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::E});
                const size_t I_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::I});

                flows[Base::template get_flat_flow_index<InfectionStateExplicit::E, InfectionStateExplicit::I>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeExposed<FP>>()[g]) * y[E_idx];

                flows[Base::template get_flat_flow_index<InfectionStateExplicit::I, InfectionStateExplicit::R>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeInfected<FP>>()[g]) * y[I_idx];
            }
        }

        for (mio::AgeGroup i(0); i < mio::AgeGroup(params.get_num_groups()); ++i) {

            std::vector<FP> N_per_age_local(static_cast<size_t>(params.get_num_groups()));
            std::vector<FP> I_per_age_local(static_cast<size_t>(params.get_num_groups()));

            for (mio::AgeGroup j(0); j < mio::AgeGroup(params.get_num_groups()); ++j) {
                FP N_j = 0;
                FP I_j = 0;

                // NonCommuter
                for (size_t state = 0; state < static_cast<size_t>(InfectionStateExplicit::Count); ++state) {
                    N_j += pop[this->populations.get_flat_index(
                        {j, CommuterType::NonCommuter, static_cast<InfectionStateExplicit>(state)})];
                }
                I_j = pop[this->populations.get_flat_index({j, CommuterType::NonCommuter, InfectionStateExplicit::I})];

                // Commuters
                for (int c2 = 0; c2 < this->m_num_commuter_groups; ++c2) {
                    CommuterType ct_j = static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c2);
                    for (size_t state = 0; state < static_cast<size_t>(InfectionStateExplicit::Count); ++state) {
                        N_j += pop[this->populations.get_flat_index(
                            {j, ct_j, static_cast<InfectionStateExplicit>(state)})];
                    }
                    I_j += pop[this->populations.get_flat_index({j, ct_j, InfectionStateExplicit::I})];
                }

                N_per_age_local[j.get()] = N_j;
                I_per_age_local[j.get()] = I_j;
            }

            // Non-commuter flow
            {
                CommuterType commuter_type_i = CommuterType::NonCommuter;
                const size_t Si = this->populations.get_flat_index({i, commuter_type_i, InfectionStateExplicit::S});
                auto S_to_E_flow_idx =
                    Base::template get_flat_flow_index<InfectionStateExplicit::S, InfectionStateExplicit::E>(
                        {i, commuter_type_i});

                flows[S_to_E_flow_idx] = 0;

                for (mio::AgeGroup j(0); j < mio::AgeGroup(params.get_num_groups()); ++j) {
                    std::vector<FP> temp_div(1);
                    temp_div[0] = (N_per_age_local[j.get()] < mio::Limits<FP>::zero_tolerance())
                                      ? 0.0
                                      : FP(1.0) / N_per_age_local[j.get()];
                    const FP divN_j = temp_div[0];

                    const FP coeffStoE =
                        params.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(t)(
                            i.get(), j.get()) *
                        params.template get<mio::oseir::TransmissionProbabilityOnContact<FP>>()[i] * divN_j;

                    flows[S_to_E_flow_idx] += coeffStoE * y[Si] * I_per_age_local[j.get()];
                }
            }

            // Commuter flows
            for (int c = 0; c < this->m_num_commuter_groups; ++c) {
                CommuterType commuter_type_i =
                    static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c);
                const size_t Si = this->populations.get_flat_index({i, commuter_type_i, InfectionStateExplicit::S});
                auto S_to_E_flow_idx =
                    Base::template get_flat_flow_index<InfectionStateExplicit::S, InfectionStateExplicit::E>(
                        {i, commuter_type_i});

                flows[S_to_E_flow_idx] = 0;

                for (mio::AgeGroup j(0); j < mio::AgeGroup(params.get_num_groups()); ++j) {
                    std::vector<FP> temp_div(1);
                    temp_div[0] = (N_per_age_local[j.get()] < mio::Limits<FP>::zero_tolerance())
                                      ? 0.0
                                      : FP(1.0) / N_per_age_local[j.get()];
                    const FP divN_j = temp_div[0];

                    const FP coeffStoE =
                        params.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(t)(
                            i.get(), j.get()) *
                        params.template get<mio::oseir::TransmissionProbabilityOnContact<FP>>()[i] * divN_j;

                    flows[S_to_E_flow_idx] += coeffStoE * y[Si] * I_per_age_local[j.get()];
                }
            }
        }
    }
};

// Type aliases for use in benchmarks
using StandardModelLagrangian = ModelStandardLagrangian<ScalarType>;
using StandardLagrangianSim   = mio::FlowSimulation<ScalarType, StandardModelLagrangian>;

} // namespace examples
} // namespace mio

#endif // MIO_Standard_Lagrangian
