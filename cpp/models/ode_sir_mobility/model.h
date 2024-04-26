
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

class Model : public FlowModel<InfectionState, Populations<Region, AgeGroup, InfectionState>, Parameters, Flows>
{

    using Base = FlowModel<InfectionState, mio::Populations<Region, AgeGroup, InfectionState>, Parameters, Flows>;

public:
    Model(int num_regions, int num_agegroups)
        : Base(Populations({Region(num_regions), AgeGroup(num_agegroups), InfectionState::Count}, 0.),
               ParameterSet(num_regions, num_agegroups))
    {
    }
    // Einmal über den Vektor und später nochmal über die Regions

    void get_flows(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                   Eigen::Ref<Eigen::VectorXd> flows) const override
    {
        auto& params     = this->parameters;
        auto& population = this->populations;

        double coeffStoI = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                           params.get<TransmissionProbabilityOnContact>() / population.get_total();

        Region n_regions     = params.get_num_regions();
        AgeGroup n_agegroups = params.get_num_agegroups();
        for (auto age = AgeGroup(0); age < n_agegroups; age++) {
            for (auto edge : params.get<CommutingRatio>()) {
                auto start_region = get<0>(edge);
                auto end_region   = get<1>(edge);
                auto strength     = get<double>(edge);
                if (start_region == end_region) {
                    continue;
                }
                // s_n += h_mn/P_m * i_m
                flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(
                    {start_region, age})] +=
                    strength * pop[population.get_flat_index({end_region, age, InfectionState::Infected})];
                // s_m += h_mn/P_m * i_n
                flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>({end_region, age})] +=
                    strength * pop[population.get_flat_index({start_region, age, InfectionState::Infected})];

                // s_n += gamma * h_nm/P_n * sum(h_km/P_k * p_nm,k * i_k)
                for (auto edge_commuter : params.get<CommutingRatio>()) {
                    auto start_region_commuter = get<0>(edge_commuter);
                    auto end_region_commuter   = get<1>(edge_commuter);
                    auto strength_commuter     = get<double>(edge_commuter);
                    if (end_region_commuter != end_region || start_region_commuter == start_region ||
                        ((std::find(params.get<PathIntersections>()[{start_region, end_region}].begin(),
                                    params.get<PathIntersections>()[{start_region, end_region}].end(),
                                    start_region_commuter)) ==
                         params.get<PathIntersections>()[{start_region, end_region}].end())) {
                        continue;
                    }
                    auto test = params.get<PathIntersections>()[{start_region, end_region}];
                    flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(
                        {start_region, age})] +=
                        params.get<ImpactCommuters>() * strength * strength_commuter *
                        pop[population.get_flat_index({start_region_commuter, age, InfectionState::Infected})];
                }
            }

            for (auto i = Region(0); i < n_regions; i++) {
                flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>({i, age})] +=
                    pop[population.get_flat_index({i, age, InfectionState::Infected})];
                flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>({i, age})] *=
                    coeffStoI * y[population.get_flat_index({i, age, InfectionState::Susceptible})];
                flows[get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>({i, age})] =
                    (1.0 / params.get<TimeInfected>()) *
                    y[population.get_flat_index({i, age, InfectionState::Infected})];
            }
        }
    }
};

} // namespace osirmobility
} // namespace mio

#endif // ODESIRMOBILITY_MODEL_H
