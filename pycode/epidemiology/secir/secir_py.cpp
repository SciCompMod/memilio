#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

#include <epidemiology/secir/secir.h>
#include <epidemiology/secir/damping.h>
#include <epidemiology/utils/time_series.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology/secir/analyze_result.h>

#include <Eigen/Core>
#include <vector>

namespace py = pybind11;

namespace
{

template <class Model>
std::vector<epi::TimeSeries<double>> filter_graph_results(
    const std::vector<epi::Graph<epi::ModelNode<epi::Simulation<Model>>, epi::MigrationEdge>>& graph_results)
{
    std::vector<epi::TimeSeries<double>> results;
    results.reserve(graph_results.size());
    std::transform(graph_results.begin(), graph_results.end(), std::back_inserter(results), [](auto&& graph) {
        return graph.nodes()[0].property.get_result();
    });
    return results;
}

// the following functions help bind class template realizations
//https://stackoverflow.com/questions/64552878/how-can-i-automatically-bind-templated-member-functions-of-variadic-class-templa
template <typename T>
std::string pretty_name()
{
    std::ostringstream o;
    o << typeid(T).name();
    return o.str();
}
template <> std::string pretty_name<epi::AgeGroup1>(){ return "AgeGroup"; }
template <> std::string pretty_name<epi::AgeGroup2>(){ return "AgeGroup"; }
template <> std::string pretty_name<epi::AgeGroup3>(){ return "AgeGroup"; }
template <> std::string pretty_name<epi::AgeGroup8>(){ return "AgeGroup"; }
template <> std::string pretty_name<epi::InfectionType>(){ return "InfectionType"; }

template <class C>
void bind_populations_members_for_all_cats(py::class_<C>&){}

template <class C, class T, class... Ts>
void bind_populations_members_for_all_cats(py::class_<C>& c)
{
    std::string tname = pretty_name<T>();
    c.def(("set_difference_from_group_total_" + tname).c_str(), &C::template set_difference_from_group_total<T>)
     .def(("set_group_total_" + tname).c_str(), &C::template set_group_total<T>)
     .def(("get_group_total_" + tname).c_str(), &C::template get_group_total<T>);

    // recursively bind the member for each type
    bind_populations_members_for_all_cats<C, Ts...>(c);
}

/*
 * @brief bind Populations class template for any choice of categories
 */
template<class... Cats>
void bind_Populations(py::module& m, std::string const& name)
{
    py::class_<epi::Populations<Cats...>> c(m, name.c_str());
    c.def(py::init<>())
        .def_static("get_num_compartments", &epi::Populations<Cats...>::get_num_compartments)
        .def("get_compartments", &epi::Populations<Cats...>::get_compartments)
        .def("get", py::overload_cast<Cats...>(&epi::Populations<Cats...>::get),
            py::return_value_policy::reference_internal)
        .def("get", py::overload_cast<Cats...>(&epi::Populations<Cats...>::get, py::const_),
            py::return_value_policy::reference_internal)
        .def("get_total", &epi::Populations<Cats...>::get_total)
        .def("set", py::overload_cast<typename epi::Populations<Cats...>::Type const &, Cats...>(&epi::Populations<Cats...>::set))
        .def("set_total", &epi::Populations<Cats...>::set_total)
        .def("set_difference_from_total", &epi::Populations<Cats...>::set_difference_from_total)
        .def("get_flat_index", &epi::Populations<Cats...>::get_flat_index);

        //get_group_total, set_group_total and set_difference_from_group_total
        bind_populations_members_for_all_cats<epi::Populations<Cats...>, Cats...>(c);
}

/*
 * @brief bind StageTimes for any number of AgeGroups N
 */
template<int N>
void bind_StageTimes(py::module& m, std::string const& name)
{
    py::class_<typename epi::SecirParams<N>::StageTimes>(m, name.c_str())
            .def(py::init<>())
            .def("set_incubation", py::overload_cast<double>(&epi::SecirParams<N>::StageTimes::set_incubation))
            .def("set_incubation",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::StageTimes::set_incubation))
            .def("set_infectious_mild", py::overload_cast<double>(&epi::SecirParams<N>::StageTimes::set_infectious_mild))
            .def("set_infectious_mild",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::StageTimes::set_infectious_mild))
            .def("set_serialinterval", py::overload_cast<double>(&epi::SecirParams<N>::StageTimes::set_serialinterval))
            .def("set_serialinterval",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::StageTimes::set_serialinterval))
            .def("set_hospitalized_to_home",
                 py::overload_cast<double>(&epi::SecirParams<N>::StageTimes::set_hospitalized_to_home))
            .def("set_hospitalized_to_home", py::overload_cast<const epi::ParameterDistribution&>(
                                                 &epi::SecirParams<N>::StageTimes::set_hospitalized_to_home))
            .def("set_home_to_hospitalized",
                 py::overload_cast<double>(&epi::SecirParams<N>::StageTimes::set_home_to_hospitalized))
            .def("set_home_to_hospitalized", py::overload_cast<const epi::ParameterDistribution&>(
                                                 &epi::SecirParams<N>::StageTimes::set_home_to_hospitalized))
            .def("set_hospitalized_to_icu",
                 py::overload_cast<double>(&epi::SecirParams<N>::StageTimes::set_hospitalized_to_icu))
            .def("set_hospitalized_to_icu", py::overload_cast<const epi::ParameterDistribution&>(
                                                &epi::SecirParams<N>::StageTimes::set_hospitalized_to_icu))
            .def("set_icu_to_home", py::overload_cast<double>(&epi::SecirParams<N>::StageTimes::set_icu_to_home))
            .def("set_icu_to_home",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::StageTimes::set_icu_to_home))
            .def("set_infectious_asymp", py::overload_cast<double>(&epi::SecirParams<N>::StageTimes::set_infectious_asymp))
            .def("set_infectious_asymp",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::StageTimes::set_infectious_asymp))
            .def("set_icu_to_death", py::overload_cast<double>(&epi::SecirParams<N>::StageTimes::set_icu_to_death))
            .def("set_icu_to_death",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::StageTimes::set_icu_to_death))

            .def("get_incubation", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_incubation),
                 py::return_value_policy::reference_internal)
            .def("get_incubation", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_incubation, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_infectious_mild", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_infectious_mild),
                 py::return_value_policy::reference_internal)
            .def("get_infectious_mild", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_infectious_mild, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_serialinterval", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_serialinterval),
                 py::return_value_policy::reference_internal)
            .def("get_serialinterval", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_serialinterval, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_hospitalized_to_home", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_hospitalized_to_home),
                 py::return_value_policy::reference_internal)
            .def("get_hospitalized_to_home",
                 py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_hospitalized_to_home, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_home_to_hospitalized", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_home_to_hospitalized),
                 py::return_value_policy::reference_internal)
            .def("get_home_to_hospitalized",
                 py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_home_to_hospitalized, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_hospitalized_to_icu", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_hospitalized_to_icu),
                 py::return_value_policy::reference_internal)
            .def("get_hospitalized_to_icu",
                 py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_hospitalized_to_icu, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_icu_to_home", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_icu_to_home),
                 py::return_value_policy::reference_internal)
            .def("get_icu_to_home", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_icu_to_home, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_infectious_asymp", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_infectious_asymp),
                 py::return_value_policy::reference_internal)
            .def("get_infectious_asymp",
                 py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_infectious_asymp, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_icu_to_dead", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_icu_to_dead),
                 py::return_value_policy::reference_internal)
            .def("get_icu_to_dead", py::overload_cast<>(&epi::SecirParams<N>::StageTimes::get_icu_to_dead, py::const_),
                 py::return_value_policy::reference_internal);
}

/*
 * @brief bind Probabilities for any number of AgeGroups N
 */
template<int N>
void bind_Probabilities(py::module& m, std::string const& name)
{
    py::class_<typename epi::SecirParams<N>::Probabilities>(m, name.c_str())
            .def(py::init<>())
            .def("set_infection_from_contact",
                 py::overload_cast<double>(&epi::SecirParams<N>::Probabilities::set_infection_from_contact))
            .def("set_infection_from_contact", py::overload_cast<const epi::ParameterDistribution&>(
                                                   &epi::SecirParams<N>::Probabilities::set_infection_from_contact))
            .def("set_carrier_infectability",
                 py::overload_cast<double>(&epi::SecirParams<N>::Probabilities::set_carrier_infectability))
            .def("set_carrier_infectability", py::overload_cast<const epi::ParameterDistribution&>(
                                                  &epi::SecirParams<N>::Probabilities::set_carrier_infectability))
            .def("set_asymp_per_infectious",
                 py::overload_cast<double>(&epi::SecirParams<N>::Probabilities::set_asymp_per_infectious))
            .def("set_asymp_per_infectious", py::overload_cast<const epi::ParameterDistribution&>(
                                                 &epi::SecirParams<N>::Probabilities::set_asymp_per_infectious))
            .def("set_risk_from_symptomatic",
                 py::overload_cast<double>(&epi::SecirParams<N>::Probabilities::set_risk_from_symptomatic))
            .def("set_risk_from_symptomatic", py::overload_cast<const epi::ParameterDistribution&>(
                                                  &epi::SecirParams<N>::Probabilities::set_risk_from_symptomatic))
            .def("set_hospitalized_per_infectious",
                 py::overload_cast<double>(&epi::SecirParams<N>::Probabilities::set_hospitalized_per_infectious))
            .def("set_hospitalized_per_infectious", py::overload_cast<const epi::ParameterDistribution&>(
                                                        &epi::SecirParams<N>::Probabilities::set_hospitalized_per_infectious))
            .def("set_icu_per_hospitalized",
                 py::overload_cast<double>(&epi::SecirParams<N>::Probabilities::set_icu_per_hospitalized))
            .def("set_icu_per_hospitalized", py::overload_cast<const epi::ParameterDistribution&>(
                                                 &epi::SecirParams<N>::Probabilities::set_icu_per_hospitalized))
            .def("set_dead_per_icu", py::overload_cast<double>(&epi::SecirParams<N>::Probabilities::set_dead_per_icu))
            .def("set_dead_per_icu",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::Probabilities::set_dead_per_icu))

            .def("get_infection_from_contact",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_infection_from_contact),
                 py::return_value_policy::reference_internal)
            .def("get_infection_from_contact",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_infection_from_contact, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_carrier_infectability",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_carrier_infectability),
                 py::return_value_policy::reference_internal)
            .def("get_carrier_infectability",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_carrier_infectability, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_asymp_per_infectious",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_asymp_per_infectious),
                 py::return_value_policy::reference_internal)
            .def("get_asymp_per_infectious",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_asymp_per_infectious, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_risk_from_symptomatic",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_risk_from_symptomatic),
                 py::return_value_policy::reference_internal)
            .def("get_test_and_trace_max_risk_from_symptomatic",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_test_and_trace_max_risk_from_symptomatic),
                 py::return_value_policy::reference_internal)
            .def("get_risk_from_symptomatic",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_risk_from_symptomatic, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_hospitalized_per_infectious",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_hospitalized_per_infectious),
                 py::return_value_policy::reference_internal)
            .def("get_hospitalized_per_infectious",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_hospitalized_per_infectious, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_icu_per_hospitalized",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_icu_per_hospitalized),
                 py::return_value_policy::reference_internal)
            .def("get_icu_per_hospitalized",
                 py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_icu_per_hospitalized, py::const_),
                 py::return_value_policy::reference_internal)
            .def("get_dead_per_icu", py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_dead_per_icu),
                 py::return_value_policy::reference_internal)
            .def("get_dead_per_icu", py::overload_cast<>(&epi::SecirParams<N>::Probabilities::get_dead_per_icu, py::const_),
                 py::return_value_policy::reference_internal);
}

/*
 * @brief bind SecirParams for any number of AgeGroups N
 */
template<int N>
void bind_SecirParams(py::module& m, std::string const& name)
{
    py::class_<typename epi::SecirParams<N>>(m, name.c_str())
            .def(py::init<>())
            .def_readwrite("times", &epi::SecirParams<N>::times)
            .def_readwrite("probabilities", &epi::SecirParams<N>::probabilities)
            .def("get_icu_capacity", py::overload_cast<>(&epi::SecirParams<N>::get_icu_capacity),
                 py::return_value_policy::reference_internal)
            .def("get_icu_capacity", py::overload_cast<>(&epi::SecirParams<N>::get_icu_capacity, py::const_),
                 py::return_value_policy::reference_internal)
            .def("set_icu_capacity", py::overload_cast<double>(&epi::SecirParams<N>::set_icu_capacity))
            .def("set_icu_capacity",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::set_icu_capacity))
            .def("get_start_day", &epi::SecirParams<N>::get_start_day)
            .def("set_start_day", &epi::SecirParams<N>::set_start_day)
            .def("get_seasonality", py::overload_cast<>(&epi::SecirParams<N>::get_seasonality),
                 py::return_value_policy::reference_internal)
            .def("get_seasonality", py::overload_cast<>(&epi::SecirParams<N>::get_seasonality, py::const_),
                 py::return_value_policy::reference_internal)
            .def("set_seasonality", py::overload_cast<double>(&epi::SecirParams<N>::set_seasonality))
            .def("set_seasonality",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::set_seasonality))
            .def("get_test_and_trace_capacity", py::overload_cast<>(&epi::SecirParams<N>::get_test_and_trace_capacity),
                 py::return_value_policy::reference_internal)
            .def("get_test_and_trace_capacity", py::overload_cast<>(&epi::SecirParams<N>::get_test_and_trace_capacity, py::const_),
                 py::return_value_policy::reference_internal)
            .def("set_test_and_trace_capacity", py::overload_cast<double>(&epi::SecirParams<N>::set_test_and_trace_capacity))
            .def("set_test_and_trace_capacity",
                 py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams<N>::set_test_and_trace_capacity))
            .def("check_constraints", &epi::SecirParams<N>::check_constraints)
            .def("apply_constraints", &epi::SecirParams<N>::apply_constraints)
            .def("get_contact_patterns", py::overload_cast<>(&epi::SecirParams<N>::get_contact_patterns),
                 py::return_value_policy::reference_internal)
            .def("get_contact_patterns", py::overload_cast<>(&epi::SecirParams<N>::get_contact_patterns, py::const_),
                 py::return_value_policy::reference_internal);
}

/*
 * @brief bind a compartmental model for any Populations and Parameters class
 */
template<class Populations, class Parameters>
void bind_CompartmentalModel(py::module& m, std::string const& name)
{
    using Model = epi::CompartmentalModel<Populations, Parameters>;
    py::class_<Model>(m, name.c_str())
            .def(py::init<>())
            .def("apply_constraints", &Model::apply_constraints)
            .def("check_constraints", &Model::check_constraints)
            .def("get_initial_values", &Model::get_initial_values)
            .def_property("populations",
                [](const Model& self) -> auto& { return self.populations; },
                [](Model& self, Populations& p) { self.populations = p; })
            .def_property("parameters",
                [](const Model& self) -> auto& { return self.parameters; },
                [](Model& self, Parameters& p) { self.parameters = p; });
}

template<class AgeGroup>
void bind_SecirModel(py::module& m, std::string const& name)
{
    using Pa = epi::SecirParams<(size_t)AgeGroup::Count>;
    using Po = epi::Populations<AgeGroup, epi::InfectionType>;
    using Model = epi::CompartmentalModel<Po, Pa>;
    py::class_<epi::SecirModel<AgeGroup>, Model>(m, name.c_str())
            .def(py::init<>());
}


/*
 * @brief bind Simulation for any number model
 */
template<class Model>
void bind_Simulation(py::module& m, std::string const& name)
{
    py::class_<epi::Simulation<Model>>(m, name.c_str())
            .def(py::init<const Model&, double, double>(), py::arg("model"), py::arg("t0") = 0,
                 py::arg("dt") = 0.1)
            .def_property_readonly("result", py::overload_cast<>(&epi::Simulation<Model>::get_result, py::const_),
                                   py::return_value_policy::reference_internal)
            .def("advance", &epi::Simulation<Model>::advance, py::arg("tmax"));
}

template <typename Model>
void bind_SecirModelNode(py::module& m, std::string const& name)
{
    py::class_<epi::Node<Model>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const epi::Node<Model>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property", [](const epi::Node<Model>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);
}

template <typename Model>
void bind_SecirSimulationNode(py::module& m, std::string const& name)
{

    using Simulation = epi::Simulation<Model>;
    py::class_<epi::Node<epi::ModelNode<Simulation>>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const epi::Node<Simulation>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property", [](const epi::Node<epi::ModelNode<Simulation>>& self) -> auto& { return self.property.model; },
            py::return_value_policy::reference_internal);
}

/*
 * @brief bind Graph for any node and edge type
 */
template<class Model>
void bind_SecirModelGraph(py::module& m, std::string const& name)
{
    using G = epi::Graph<Model, epi::MigrationEdge>;
    py::class_<G>(m, name.c_str())
            .def(py::init<>())
            .def("add_node", &G::template add_node<const Model&>,
                 py::return_value_policy::reference_internal)
            .def("add_edge", &G::template add_edge<const epi::MigrationEdge&>,
                 py::return_value_policy::reference_internal)
            .def("add_edge", &G::template add_edge<const Eigen::VectorXd&>,
                 py::return_value_policy::reference_internal)
            .def_property_readonly("num_nodes",
                                   [](const G& self) {
                                       return self.nodes().size();
                                   })
            .def(
                "get_node", [](const G& self, size_t node_idx) -> auto& { return self.nodes()[node_idx]; },
                py::return_value_policy::reference_internal)
            .def_property_readonly("num_edges",
                                   [](const G& self) {
                                       return self.edges().size();
                                   })
            .def(
                "get_edge", [](const G& self, size_t edge_idx) -> auto& { return self.edges()[edge_idx]; },
                py::return_value_policy::reference_internal)
            .def("get_num_out_edges",
                 [](const G& self, size_t node_idx) {
                     return self.out_edges(node_idx).size();
                 })
            .def(
                "get_out_edge", [](const G& self, size_t node_idx, size_t edge_idx) -> auto& {
                    return self.out_edges(node_idx)[edge_idx];
                },
                py::return_value_policy::reference_internal);

}

template<class Simulation, class Model>
void bind_MigrationGraph(py::module& m, std::string const& name)
{
    using G = epi::Graph<epi::ModelNode<Simulation>, epi::MigrationEdge>;
    py::class_<G>(m, name.c_str())
            .def(py::init<>())
            .def(
                "add_node", [](G & self, int id, const Model& p, double t0, double dt) -> auto& {
                    return self.add_node(id, p, t0, dt);
                },
                py::arg("id"), py::arg("model"), py::arg("t0") = 0.0, py::arg("dt") = 0.1, py::return_value_policy::reference_internal)
            .def("add_edge", &G::template add_edge<const epi::MigrationEdge&>,
                 py::return_value_policy::reference_internal)
            .def("add_edge", &G::template add_edge<const Eigen::VectorXd&>, py::return_value_policy::reference_internal)
            .def_property_readonly("num_nodes",
                                   [](const G& self) {
                                       return self.nodes().size();
                                   })
            .def(
                "get_node",
                [](const G& self, size_t node_idx) -> auto& { return self.nodes()[node_idx]; },
                py::return_value_policy::reference_internal)
            .def_property_readonly("num_edges",
                                   [](const G& self) {
                                       return self.edges().size();
                                   })
            .def(
                "get_edge", [](const G& self, size_t edge_idx) -> auto& { return self.edges()[edge_idx]; },
                py::return_value_policy::reference_internal)
            .def("get_num_out_edges",
                 [](const G& self, size_t node_idx) {
                     return self.out_edges(node_idx).size();
                 })
            .def(
                "get_out_edge", [](const G& self, size_t node_idx, size_t edge_idx) -> auto& {
                    return self.out_edges(node_idx)[edge_idx];
                },
                py::return_value_policy::reference_internal);

}

/*
 * @brief bind GraphSimulation for any node and edge type
 */
template<class Graph>
void bind_GraphSimulation(py::module& m, std::string const& name)
{
    using GS = epi::GraphSimulation<Graph>;
    py::class_<GS>(m, name.c_str())
        .def(py::init([](const Graph& graph, double t0, double dt) {
                 return std::make_unique<GS>(epi::make_migration_sim(t0, dt, graph));
             }),
             py::arg("graph"), py::arg("t0") = 0.0, py::arg("dt") = 1.0)
        .def_property_readonly("graph",
                               py::overload_cast<>(&GS::get_graph, py::const_),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("t", &GS::get_t)
        .def("advance", &GS::advance, py::arg("tmax"));
}

/*
 * @brief bind ParameterStudy for any model
 */
template<class Model>
void bind_ParameterStudy(py::module& m, std::string const& name)
{
    py::class_<epi::ParameterStudy<Model>>(m, name.c_str())
        .def(py::init<const Model&, double, double, size_t>(), py::arg("model"), py::arg("t0"),
             py::arg("tmax"), py::arg("num_runs"))
        .def(py::init<const Model&, double, double, double, size_t>(), py::arg("model"), py::arg("t0"),
             py::arg("tmax"), py::arg("dev_rel"), py::arg("num_runs"))
        .def(py::init<const epi::Graph<Model, epi::MigrationEdge>&, double, double, double, size_t>(), py::arg("model_graph"),
             py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("num_runs"))
        .def_property("num_runs", &epi::ParameterStudy<Model>::get_num_runs, &epi::ParameterStudy<Model>::set_num_runs)
        .def_property("tmax", &epi::ParameterStudy<Model>::get_tmax, &epi::ParameterStudy<Model>::set_tmax)
        .def_property("t0", &epi::ParameterStudy<Model>::get_t0, &epi::ParameterStudy<Model>::set_t0)
        .def_property_readonly("model", py::overload_cast<>(&epi::ParameterStudy<Model>::get_model),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("model", py::overload_cast<>(&epi::ParameterStudy<Model>::get_model, py::const_),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph", py::overload_cast<>(&epi::ParameterStudy<Model>::get_secir_model_graph),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               py::overload_cast<>(&epi::ParameterStudy<Model>::get_secir_model_graph, py::const_),
                               py::return_value_policy::reference_internal)
        .def("run", &epi::ParameterStudy<Model>::run)
        .def("run",
             [](epi::ParameterStudy<Model>& self) { //default argument doesn't seem to work with functions
                 return self.run();
             })
        .def(
            "run_single",
            [](epi::ParameterStudy<Model>& self, typename epi::ParameterStudy<Model>::HandleSimulationResultFunction handle_result_func) {
                return filter_graph_results(self.run(handle_result_func));
            },
            py::arg("handle_result_func"))
        .def("run_single", [](epi::ParameterStudy<Model>& self) {
            return filter_graph_results(self.run());
        });
}

/*
 *@brief bind the relevant classes/functions of an age resolved secir model for any agegroup enum
 */
template <class AgeGroup>
void bind_secir_ageres(py::module& m)
{
    size_t constexpr N = (size_t)AgeGroup::Count;

    py::enum_<AgeGroup> agegroup_enum(m, ("AgeGroup" + std::to_string(N)).c_str());
    for (size_t i=0; i < (size_t)AgeGroup::Count; ++i) {
        agegroup_enum.value(("Group" + std::to_string(i)).c_str(), (AgeGroup)i);
    }
    agegroup_enum.value("Count", AgeGroup::Count).export_values();

    bind_Populations<AgeGroup, epi::InfectionType>(m, "Populations" + std::to_string(N));

    bind_SecirParams<N>(m, "SecirParams" + std::to_string(N));

    //TODO: Currently, StageTimes and Probabilies are subclasses of the class template SecirParams.
    // If this were not the case, we would not really need a StageTimes class per
    // number of AgeGroups, but since we are planning to remove the SecirParams class
    // this is a valid workaround for now.
    bind_StageTimes<N>(m, "StageTimes" + std::to_string(N));
    bind_Probabilities<N>(m, "Probabilities" + std::to_string(N));

    using Populations = epi::Populations<AgeGroup, epi::InfectionType>;
    using SecirParams = epi::SecirParams<N>;
    bind_CompartmentalModel<Populations, SecirParams>(m, "SecirModelBase" + std::to_string(N));  //<- no flows
    bind_SecirModel<AgeGroup>(m, "SecirModel" + std::to_string(N));  //<-- flows defined in derived constructor at rt

    using SecirModel = epi::SecirModel<AgeGroup>;
    bind_Simulation<SecirModel>(m, "SecirSimulation" + std::to_string(N));

    m.def("simulate", [](double t0, double tmax, double dt, const SecirModel& model) { return epi::simulate<SecirModel>(t0, tmax, dt, model); },
          "Simulates a SecirModel1 from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
          py::arg("model"));

    bind_SecirModelNode<SecirModel>(m, "SecirModelNode" + std::to_string(N));
    bind_SecirSimulationNode<SecirModel>(m, "SecirSimulationNode" + std::to_string(N));
    bind_SecirModelGraph<SecirModel>(m, "SecirModelGraph" + std::to_string(N));


    using Simulation = epi::Simulation<SecirModel>;
    bind_MigrationGraph<Simulation, SecirModel>(m, "MigrationGraph" + std::to_string(N));

    using MigrationGraph = epi::Graph<epi::ModelNode<Simulation>, epi::MigrationEdge>;
    bind_GraphSimulation<MigrationGraph>(m, "MigrationSimulation" + std::to_string(N));

    bind_ParameterStudy<SecirModel>(m, "ParameterStudy" + std::to_string(N));

    m.def("set_params_distributions_normal", &epi::set_params_distributions_normal<AgeGroup>, py::arg("model"), py::arg("t0"),
          py::arg("tmax"), py::arg("dev_rel"));

    m.def("draw_sample", &epi::draw_sample<AgeGroup>, py::arg("model"));

    m.def("interpolate_simulation_result",
          py::overload_cast<const MigrationGraph&>(&epi::interpolate_simulation_result<Simulation>));

    m.def("interpolate_ensemble_results", &epi::interpolate_ensemble_results<MigrationGraph>);
}

} // namespace

PYBIND11_MODULE(_secir, m)
{

    py::class_<epi::Damping>(m, "Damping")
        .def(py::init([](const Eigen::Ref<const Eigen::MatrixXd>& c, double t, int level, int type) {
                 return epi::Damping(c, epi::DampingLevel(level), epi::DampingType(type), epi::SimulationTime(t));
             }),
             py::arg("coeffs"), py::arg("t"), py::arg("level") = 0, py::arg("type") = 0)
        .def_property(
            "coeffs", [](const epi::Damping& self) -> const auto& { return self.get_coeffs(); },
            [](epi::Damping& self, const Eigen::Ref<const Eigen::MatrixXd>& v) {
                self.get_coeffs() = v;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "time",
            [](const epi::Damping& self) {
                return self.get_time();
            },
            [](epi::Damping& self, double v) {
                self.get_time() = epi::SimulationTime(v);
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "type",
            [](const epi::Damping& self) {
                return self.get_type();
            },
            [](epi::Damping& self, int v) {
                self.get_type() = epi::DampingType(v);
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "level",
            [](const epi::Damping& self) {
                return self.get_level();
            },
            [](epi::Damping& self, int v) {
                self.get_level() = epi::DampingLevel(v);
            },
            py::return_value_policy::reference_internal);

    py::class_<epi::Dampings>(m, "Dampings")
        .def(py::init<Eigen::Index>())
        .def("add",
             [](epi::Dampings& self, const epi::Damping& d) {
                 self.add(d);
             })
        .def("get_matrix_at", [](const epi::Dampings& self, double t) {
            return self.get_matrix_at(t);
        });

    py::class_<epi::TimeSeries<double>>(m, "TimeSeries")
        .def(py::init<Eigen::Index>(), py::arg("num_elements"))
        .def("get_num_time_points", &epi::TimeSeries<double>::get_num_time_points)
        .def("get_num_elements", &epi::TimeSeries<double>::get_num_elements)
        .def("get_time", py::overload_cast<Eigen::Index>(&epi::TimeSeries<double>::get_time), py::arg("index"))
        .def("get_last_time", py::overload_cast<>(&epi::TimeSeries<double>::get_last_time))
        .def("get_value", py::overload_cast<Eigen::Index>(&epi::TimeSeries<double>::get_value), py::arg("index"))
        .def("get_last_value", py::overload_cast<>(&epi::TimeSeries<double>::get_last_value))
        .def(
            "__getitem__",
            [](epi::TimeSeries<double>& self, Eigen::Index i) {
                if (i >= 0 && i < self.get_num_time_points()) {
                    return self[i];
                }
                else {
                    throw pybind11::index_error("Index out of range."); //needs to throw exception for iterable
                }
            },
            py::is_operator(), py::arg("index"))
        .def(
            "__setitem__",
            [](epi::TimeSeries<double>& self, Eigen::Index i, Eigen::Ref<const epi::TimeSeries<double>::Vector> expr) {
                if (i >= 0 && i < self.get_num_time_points()) {
                    self[i] = expr;
                }
                else {
                    throw pybind11::index_error("Index out of range."); //needs to throw exception for iterable
                }
            },
            py::is_operator(), py::arg("index"), py::arg("v"))
        .def("add_time_point",
             [](epi::TimeSeries<double>& self) {
                 return self.add_time_point();
             })
        .def("add_time_point",
             [](epi::TimeSeries<double>& self, double t) {
                 return self.add_time_point(t);
             })
        .def("add_time_point",
             [](epi::TimeSeries<double>& self, double t, Eigen::Ref<const epi::TimeSeries<double>::Vector> expr) {
                 return self.add_time_point(t, expr);
             })
        .def("as_ndarray", [](epi::TimeSeries<double>& self) {
            auto m = Eigen::Map<epi::TimeSeries<double>::Matrix>(self.data(), self.get_num_rows(),
                                                                 self.get_num_time_points());
            return Eigen::Ref<epi::TimeSeries<double>::Matrix>(m);
        });

    py::class_<epi::ParameterDistribution>(m, "ParameterDistribution")
        .def_property("lower_bound", &epi::ParameterDistribution::get_lower_bound,
                      &epi::ParameterDistribution::set_lower_bound)
        .def_property("upper_bound", &epi::ParameterDistribution::get_upper_bound,
                      &epi::ParameterDistribution::set_upper_bound)
        .def("add_predefined_sample", &epi::ParameterDistribution::add_predefined_sample)
        .def("remove_predefined_samples", &epi::ParameterDistribution::remove_predefined_samples)
        .def("get_sample", &epi::ParameterDistribution::get_sample);

    py::class_<epi::ParameterDistributionNormal, epi::ParameterDistribution>(m, "ParameterDistributionNormal")
        .def(py::init<double, double, double, double>(), py::arg("lb"), py::arg("ub"), py::arg("mean"),
             py::arg("std_dev"))
        .def(py::init<double, double, double>(), py::arg("lb"), py::arg("ub"), py::arg("mean"))
        .def_property("mean", &epi::ParameterDistributionNormal::get_mean, &epi::ParameterDistributionNormal::set_mean)
        .def_property("standard_dev", &epi::ParameterDistributionNormal::get_standard_dev,
                      &epi::ParameterDistributionNormal::set_standard_dev);

    py::class_<epi::ParameterDistributionUniform, epi::ParameterDistribution>(m, "ParameterDistributionUniform")
        .def(py::init<>())
        .def(py::init<double, double>(), py::arg("lb"), py::arg("ub"));

    py::class_<epi::UncertainValue>(m, "UncertainValue")
        .def(py::init<ScalarType>(), py::arg("value") = 0.0)
        .def_property(
            "value",
            [](epi::UncertainValue& self) {
                return ScalarType(self);
            },
            [](epi::UncertainValue& self, ScalarType v) {
                self = v;
            })
        .def("set_distribution", //a property would be nicer but getter and setter use a different type
             &epi::UncertainValue::set_distribution)
        .def(
            "get_distribution",
            [](const epi::UncertainValue& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def(
            "get_distribution",
            [](epi::UncertainValue& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def("draw_sample", &epi::UncertainValue::draw_sample);


    py::class_<epi::ContactMatrix>(m, "ContactMatrix")
        .def(py::init<const Eigen::Ref<const Eigen::MatrixXd>&, const Eigen::Ref<const Eigen::MatrixXd>&>(),
             py::arg("baseline"), py::arg("minimum"))
        .def(py::init<const Eigen::Ref<const Eigen::MatrixXd>&>(), py::arg("baseline"))
        .def(py::init<Eigen::Index>())
        .def("add_damping",
             [](epi::ContactMatrix& self, const epi::Damping& d) {
                 self.add_damping(d);
             })
        .def_property(
            "baseline", [](const epi::ContactMatrix& self) -> auto& { return self.get_baseline(); },
            [](epi::ContactMatrix& self, const Eigen::Ref<const Eigen::MatrixXd>& v) {
                self.get_baseline() = v;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "minimum", [](const epi::ContactMatrix& self) -> auto& { return self.get_minimum(); },
            [](epi::ContactMatrix& self, const Eigen::Ref<const Eigen::MatrixXd>& v) {
                self.get_minimum() = v;
            },
            py::return_value_policy::reference_internal)
        .def_property_readonly("num_groups",
                               [](const epi::ContactMatrix& self) {
                                   return self.get_num_groups();
                               })
        .def("get_dampings",
             [](const epi::ContactMatrix& self) {
                 return std::vector<epi::Damping>(self.get_dampings().begin(), self.get_dampings().end());
             })
        .def("get_matrix_at", [](const epi::ContactMatrix& self, double t) {
            return self.get_matrix_at(t);
        });

    py::class_<epi::ContactMatrixGroup>(m, "ContactMatrixGroup")
        .def(py::init<Eigen::Index, size_t>(), py::arg("num_groups") = 1, py::arg("num_matrices") = 1)
        .def("add_damping",
             [](epi::ContactMatrixGroup& self, const epi::Damping& d) {
                 self.add_damping(d);
             })
        .def_property_readonly("num_groups",
                               [](const epi::ContactMatrixGroup& self) {
                                   return self.get_num_groups();
                               })
        .def_property_readonly("num_matrices",
                               [](const epi::ContactMatrixGroup& self) {
                                   return self.get_num_matrices();
                               })
        .def(
            "__getitem__", [](epi::ContactMatrixGroup & self, size_t i) -> auto& {
                if (i < 0 || i >= self.get_num_matrices()) {
                    throw py::index_error("index out of range");
                }
                return self[i];
            },
            py::return_value_policy::reference_internal)
        .def("__setitem__",
             [](epi::ContactMatrixGroup& self, size_t i, const epi::ContactMatrix& m) {
                 if (i < 0 && i >= self.get_num_matrices()) {
                     throw py::index_error("index out of range");
                 }
                 self[i] = m;
             })
        .def("get_matrix_at", [](const epi::ContactMatrixGroup& self, double t) {
            return self.get_matrix_at(t);
        });

    py::class_<epi::UncertainContactMatrix>(m, "UncertainContactMatrix")
        .def(py::init<>())
        .def(py::init<const epi::ContactMatrixGroup&>())
        .def_property(
            "cont_freq_mat",
            [](const epi::UncertainContactMatrix& self) {
                return self.get_cont_freq_mat();
            },
            [](epi::UncertainContactMatrix& self, const epi::ContactMatrixGroup& c) {
                self.get_cont_freq_mat() = c;
            },
            py::return_value_policy::reference_internal);

    py::class_<epi::MigrationEdge>(m, "MigrationParams")
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def_property(
            "coefficients", [](const epi::MigrationEdge& self) -> auto { return self.get_coefficients(); },
            [](epi::MigrationEdge& self, const Eigen::VectorXd& v) {
                self.get_coefficients() = v;
            },
            py::return_value_policy::reference_internal);

    py::class_<epi::Edge<epi::MigrationEdge>>(m, "MigrationEdge")
        .def_property_readonly("start_node_idx",
                               [](const epi::Edge<epi::MigrationEdge>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const epi::Edge<epi::MigrationEdge>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property", [](const epi::Edge<epi::MigrationEdge>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);

    // https://github.com/pybind/pybind11/issues/1153
    m.def("interpolate_simulation_result",
          static_cast<epi::TimeSeries<double> (*)(const epi::TimeSeries<double>&)>(&epi::interpolate_simulation_result));


    m.def("interpolate_ensemble_results", &epi::interpolate_ensemble_results<epi::TimeSeries<double>>);

    m.def("ensemble_mean", &epi::ensemble_mean);
    m.def("ensemble_percentile", &epi::ensemble_percentile);

    py::enum_<epi::InfectionType>(m, "InfectionType")
        .value("S", epi::InfectionType::S)
        .value("E", epi::InfectionType::E)
        .value("C", epi::InfectionType::C)
        .value("I", epi::InfectionType::I)
        .value("H", epi::InfectionType::H)
        .value("U", epi::InfectionType::U)
        .value("R", epi::InfectionType::R)
        .value("D", epi::InfectionType::D)
        .value("Count", epi::InfectionType::Count)
        .export_values();

    bind_secir_ageres<epi::AgeGroup1>(m);
    bind_secir_ageres<epi::AgeGroup2>(m);
    bind_secir_ageres<epi::AgeGroup3>(m);
    bind_secir_ageres<epi::AgeGroup8>(m);

    m.attr("__version__") = "dev";
}
