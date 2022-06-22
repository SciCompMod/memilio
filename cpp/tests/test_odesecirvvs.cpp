#include "memilio/epidemiology/age_group.h"
#include "memilio/mobility/mobility.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/uncertain_value.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameter_space.h"
#include "ode_secirvvs/parameters.h"

#include "gtest/gtest.h"
#include "gmock/gmock-matchers.h"
#include <gmock/gmock-generated-matchers.h>

void set_distribution(mio::UncertainValue& uv, double v, double min, double max)
{
    uv        = v;
    auto dist = mio::ParameterDistributionUniform(min, max);
    uv.set_distribution(dist);
}

void set_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array, double v, double min, double max)
{
    for (auto& uv : array) {
        set_distribution(uv, v, min, max);
    }
}

template<class... T>
void set_distribution(T&&...)
{
    //overload for parameters that aren't sampled
}

bool check_sample(const mio::UncertainValue& uv, double min, double max)
{
    return uv.value() >= min && uv.value() <= max;
}

bool check_sample(const mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array, double min, double max)
{
    return std::all_of(array.begin(), array.end(), [=](auto&& uv) { return check_sample(uv, min, max); });
}

template<class... T>
bool check_sample(T&&...)
{
    //overload for parameters that aren't sampled
    return true;
}

MATCHER_P2(IsSampled, min, max, std::string(negation ? "isn't" : "is") + " sampled correctly.")
{
    mio::unused(result_listener);
    return check_sample(arg, min, max);
}

TEST(TestOdeSECIRVVS, draw_sample)
{
    mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> graph;

    auto num_age_groups = 2;
    graph.add_node(0, num_age_groups);
    graph.add_node(1, num_age_groups);
    graph.add_edge(0, 1, Eigen::VectorXd::Zero(Eigen::Index(mio::osecirvvs::InfectionState::Count) * num_age_groups));
    for (auto& n : graph.nodes()) {
        mio::foreach(n.property.parameters, [](auto&& p, auto&& /*t*/) {
            set_distribution(p, 0.0, 0.9, 1.1);
        });
    }

    auto sampled_graph = mio::osecirvvs::draw_sample(graph, true, false);

    //check sampling
    auto& parameters0 = sampled_graph.nodes()[0].property.parameters;
    mio::foreach(parameters0, [](auto&& p, auto&& /*t*/) {
        ASSERT_THAT(p, IsSampled(0.9, 1.1));
    });

    //check parameters that should be equal between nodes
    auto& parameters1 = sampled_graph.nodes()[1].property.parameters;
    ASSERT_THAT(parameters1.get<mio::osecirvvs::InfectiousTimeMild>(),
                testing::ElementsAreArray(parameters0.get<mio::osecirvvs::InfectiousTimeMild>().array().data(), num_age_groups));
}