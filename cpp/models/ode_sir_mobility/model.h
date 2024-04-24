
#ifndef ODESIRMOBILITY_MODEL_H
#define ODESIRMOBILITY_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/epidemiology/populations.h"
#include "ode_sir_mobility/infection_state.h"
#include "ode_sir_mobility/parameters.h"
#include "ode_sir_mobility/regions.h"

namespace mio
{
namespace osirmobility
{

/********************
    * define the model *
    ********************/

using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Infected>,
                       Flow<InfectionState::Infected, InfectionState::Recovered>>;

class Model : public FlowModel<InfectionState, Populations<Region, InfectionState>, Parameters, Flows>
{

    using Base = FlowModel<InfectionState, mio::Populations<Region, InfectionState>, Parameters, Flows>;

public:
    Model(int num_regions)
        : Base(Populations({Region(num_regions), InfectionState::Count}, 0.), ParameterSet(num_regions))
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

        Region n_regions = params.get_num_regions();

        for (auto edge : params.get<CommutingRatio>()) {
            auto start_region = get<0>(edge);
            auto end_region   = get<1>(edge);
            auto strength     = get<double>(edge);
            // s_n += h_mn/P_m * i_m
            flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(start_region)] +=
                strength * pop[population.get_flat_index({end_region, InfectionState::Infected})];
            // s_m += h_mn/P_m * i_n
            flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>(end_region)] +=
                strength * pop[population.get_flat_index({start_region, InfectionState::Infected})];
        }

        for (auto i = Region(0); i < n_regions; i++) {
            flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>({i})] +=
                pop[population.get_flat_index({i, InfectionState::Infected})];
            flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infected>({i})] *=
                coeffStoI * y[population.get_flat_index({i, InfectionState::Susceptible})];
            flows[get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>({i})] =
                (1.0 / params.get<TimeInfected>()) * y[population.get_flat_index({i, InfectionState::Infected})];
        }
    }
};

} // namespace osirmobility
} // namespace mio

#endif // ODESIRMOBILITY_MODEL_H
