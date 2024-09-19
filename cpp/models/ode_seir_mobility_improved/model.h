
#ifndef ODESEIRMOBILITYIMPROVED_MODEL_H
#define ODESEIRMOBILITYIMPROVED_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/epidemiology/populations.h"
#include "models/ode_seir_mobility_improved/infection_state.h"
#include "models/ode_seir_mobility_improved/parameters.h"
#include "models/ode_seir_mobility_improved/regions.h"
#include "memilio/epidemiology/age_group.h"

namespace mio
{
namespace oseirmobilityimproved
{

/********************
    * define the model *
    ********************/

using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Exposed>,
                       Flow<InfectionState::Exposed, InfectionState::Infected>,
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
        const auto& params              = this->parameters;
        const auto& population          = this->populations;
        const auto& commuting_strengths = params.template get<CommutingStrengths<>>();

        const Index<AgeGroup> n_age_groups = reduce_index<Index<AgeGroup>>(params.get_num_agegroups());
        const Index<Region> n_regions      = reduce_index<Index<Region>>(params.get_num_regions());

        CustomIndexArray<double, Region> infectives_per_region(n_regions);
        for (auto age_i : make_index_range(n_age_groups)) {
            for (auto age_j : make_index_range(n_age_groups)) {
                for (int n = 0; n < commuting_strengths.outerSize(); ++n) {
                    infectives_per_region[Region(n)] =
                        0.5 * pop[population.get_flat_index({Region(n), age_j, InfectionState::Infected})] *
                        commuting_strengths.coeff(n, n);
                    for (Eigen::SparseMatrix<double>::InnerIterator it(commuting_strengths, n); it; ++it) {
                        infectives_per_region[Region(n)] +=
                            0.5 * pop[population.get_flat_index({Region(it.row()), age_j, InfectionState::Infected})] *
                            it.value();
                    }
                }
                double coeffStoI = params.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(t)(
                                       age_i.get(), age_j.get()) *
                                   params.template get<TransmissionProbabilityOnContact<FP>>()[age_i];
                for (int n = 0; n < commuting_strengths.outerSize(); ++n) {
                    for (Eigen::SparseMatrix<double>::InnerIterator it(commuting_strengths, n); it; ++it) {
                        flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>(
                            {Region(it.row()), age_i})] += it.value() * infectives_per_region[Region(it.col())] /
                                                           params.template get<PopulationSizes<FP>>()[Region(it.col())];
                    }
                }
                for (auto region : make_index_range(n_regions)) {
                    flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>(
                        {region, age_i})] *=
                        coeffStoI * y[population.get_flat_index({region, age_j, InfectionState::Susceptible})];
                    flows[Base::template get_flat_flow_index<InfectionState::Exposed, InfectionState::Infected>(
                        {region, age_i})] = (1.0 / params.template get<TimeExposed<FP>>()[age_i]) *
                                            y[population.get_flat_index({region, age_i, InfectionState::Exposed})];
                    flows[Base::template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>(
                        {region, age_i})] = (1.0 / params.template get<TimeInfected<FP>>()[age_i]) *
                                            y[population.get_flat_index({region, age_i, InfectionState::Infected})];
                }
            }
        }
    }
}; // namespace oseirmobilityimproved

} // namespace oseirmobilityimproved
} // namespace mio

#endif // ODESEIRMOBILITY_MODEL_H
