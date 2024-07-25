
#ifndef ODESIRMOBILITY_MODEL_H
#define ODESIRMOBILITY_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/epidemiology/populations.h"
#include "ode_sir_mobility/infection_state.h"
#include "ode_sir_mobility/parameters.h"
#include "ode_sir_mobility/regions.h"
#include "memilio/epidemiology/age_group.h"

namespace mio
{
namespace osirmobility
{

/********************
    * define the model *
    ********************/

using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Infected>,
                       Flow<InfectionState::Infected, InfectionState::Recovered>>;

template <typename FP = ScalarType>
class Model : public FlowModel<FP, InfectionState, mio::Populations<FP, Region, AgeGroup, InfectionState>,
                               Parameters<FP>, Flows>
{

    using Base =
        FlowModel<FP, InfectionState, mio::Populations<FP, Region, AgeGroup, InfectionState>, Parameters<FP>, Flows>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;

    Model(int num_regions, int num_agegroups)
        : Base(Populations({Region(num_regions), AgeGroup(num_agegroups), InfectionState::Count}),
               ParameterSet(Region(num_regions), AgeGroup(num_agegroups)))
    {
    }
    // Einmal über den Vektor und später nochmal über die Regions

    void get_flows(Eigen::Ref<const Vector<FP>> pop, Eigen::Ref<const Vector<FP>> y, FP t,
                   Eigen::Ref<Vector<FP>> flows) const override
    {
        const auto& params     = this->parameters;
        const auto& population = this->populations;

        const Index<AgeGroup> n_age_groups = reduce_index<Index<AgeGroup>>(params.get_num_agegroups());
        const Index<Region> n_regions      = reduce_index<Index<Region>>(params.get_num_regions());

        for (auto age_i : make_index_range(n_age_groups)) {
            for (auto age_j : make_index_range(n_age_groups)) {
                double coeffStoI = params.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(t)(
                                       age_i.get(), age_j.get()) *
                                   params.template get<TransmissionProbabilityOnContact<FP>>()[age_i] /
                                   population.get_group_total(age_j);
                for (auto edge : params.template get<CommutingRatio>()) {
                    auto start_region = get<0>(edge);
                    auto end_region   = get<1>(edge);
                    auto strength     = get<double>(edge);
                    if (start_region == end_region) {
                        continue;
                    }
                    // s_n += h_mn/P_m * i_m
                    flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(
                        {start_region, age_i})] +=
                        strength * pop[population.get_flat_index({end_region, age_j, InfectionState::Infected})];
                    // s_m += h_mn/P_m * i_n
                    flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(
                        {end_region, age_i})] +=
                        strength * pop[population.get_flat_index({start_region, age_j, InfectionState::Infected})];

                    // s_n += gamma * h_nm/P_n * sum(h_km/P_k * p_nm,k * i_k)
                    for (auto edge_commuter : params.template get<CommutingRatio>()) {
                        auto start_region_commuter = get<0>(edge_commuter);
                        auto end_region_commuter   = get<1>(edge_commuter);
                        auto strength_commuter     = get<double>(edge_commuter);
                        if (end_region_commuter != end_region || start_region_commuter == start_region ||
                            ((std::find(params.template get<PathIntersections>()[{start_region, end_region}].begin(),
                                        params.template get<PathIntersections>()[{start_region, end_region}].end(),
                                        start_region_commuter)) ==
                             params.template get<PathIntersections>()[{start_region, end_region}].end())) {
                            continue;
                        }
                        flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(
                            {start_region, age_i})] +=
                            params.template get<ImpactTransmissionDuringCommuting<FP>>() * strength *
                            strength_commuter *
                            pop[population.get_flat_index({start_region_commuter, age_j, InfectionState::Infected})];
                    }
                }
                for (auto region : make_index_range(n_regions)) {
                    flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(
                        {region, age_i})] += pop[population.get_flat_index({region, age_j, InfectionState::Infected})];
                    flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(
                        {region, age_i})] *=
                        coeffStoI * y[population.get_flat_index({region, age_j, InfectionState::Susceptible})];
                }
            }

            for (auto region : make_index_range(n_regions)) {
                flows[Base::template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>(
                    {region, age_i})] = (1.0 / params.template get<TimeInfected<FP>>()[age_i]) *
                                        y[population.get_flat_index({region, age_i, InfectionState::Infected})];
            }
        }
    }
};

} // namespace osirmobility
} // namespace mio

#endif // ODESIRMOBILITY_MODEL_H
