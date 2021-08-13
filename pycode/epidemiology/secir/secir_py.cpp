/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin Siggel, Daniel Abele, Martin J. Kuehn, Jan Kleinert
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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

#include <epidemiology/utils/custom_index_array.h>
#include <epidemiology/secir/secir.h>
#include <epidemiology/secir/damping.h>
#include <epidemiology/utils/regions.h>
#include <epidemiology/utils/time_series.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology/secir/analyze_result.h>

#include <Eigen/Core>
#include <vector>

namespace py = pybind11;

namespace
{

//select only the first node of the graph of each run, used for parameterstudy with single nodes
template<class Sim>
std::vector<Sim> filter_graph_results(
    const std::vector<epi::Graph<epi::SimulationNode<Sim>, epi::MigrationEdge>>& graph_results)
{
    std::vector<Sim> results;
    results.reserve(graph_results.size());
    std::transform(graph_results.begin(), graph_results.end(), std::back_inserter(results), [](auto&& graph) {
        return graph.nodes()[0].property.get_simulation();
    });
    return results;
}

// bind an index for a single tag
template <class Tag>
void bind_Index(py::module& m, std::string const& name)
{
    py::class_<epi::Index<Tag>> c(m, name.c_str());
    c.def(py::init<size_t>(), py::arg("value"));
}

// helper function for implicitly casting from py::tuple to Index in Python.
// This extracts an Index from a py::tuple of Indices from the correct position,
// given the corresponding Index type
template <typename Tag, class Tuple>
epi::Index<Tag> extract_index(py::tuple& t)
{
    return t[epi::details::IndexPosition<Tag, Tuple>::value].template cast<epi::Index<Tag>>();
}

// bind an index for more than one tag
template <class... Tags>
void bind_MultiIndex(py::module&m, std::string const& name)
{
    using C = epi::Index<Tags...>;
    py::class_<C> c(m, name.c_str());
    c.def(py::init<epi::Index<Tags> const&...>())
            .def(py::init([](py::tuple t){
                     return C(extract_index<Tags, C>(t)...);
                 }
                 ));

    py::implicitly_convertible<py::tuple, C>();
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
template <> std::string pretty_name<epi::InfectionState>(){ return "InfectionState"; }
template <> std::string pretty_name<epi::AgeGroup>(){ return "AgeGroup"; }


template <class C>
void bind_templated_members_CustomIndexArray(py::class_<C>&){}

template <class C, class T, class... Ts>
void bind_templated_members_CustomIndexArray(py::class_<C>& c)
{
    std::string tname = pretty_name<T>();
    c.def(("size_" + tname).c_str(), &C::template size<T>);

    // recursively bind the member for each type
    bind_templated_members_CustomIndexArray<C, Ts...>(c);
}

template <class Type, class... Tags>
void bind_CustomIndexArray(py::module& m, std::string const& name)
{
    using C = typename epi::CustomIndexArray<Type, Tags...>;
    using Index = typename epi::CustomIndexArray<Type, Tags...>::Index;
    py::class_<C> c(m, name.c_str());
    c.def(py::init([](Index const& sizes, Type const& val){ return C(sizes, val); }))
     .def(py::init([](Index const& sizes){ return C(sizes); }))
     .def("numel", &C::numel)
     .def("__getitem__",
          [](const C& self, Index const& idx) -> auto& { return self[idx]; },
         py::return_value_policy::reference_internal)
     .def("__setitem__",
          [](C& self, Index const& idx, double value) { self[idx] = value; })
     .def("__iter__", [](const C &s) { return py::make_iterator(s.begin(), s.end()); },
          py::keep_alive<0, 1>())
     .def("get_flat_index", &C::get_flat_index);

    // Not supported in Python yet: Slicing

    // bind all templated members for types in Tags...
    bind_templated_members_CustomIndexArray<C, Tags...>(c);
}


template <class C, class Base>
void bind_templated_members_Population(py::class_<C, Base>&){}

template <class C, class Base, class T, class... Ts>
void bind_templated_members_Population(py::class_<C, Base>& c)
{
    std::string tname = pretty_name<T>();
    c.def(("set_difference_from_group_total_" + tname).c_str(), &C::template set_difference_from_group_total<T>)
     .def(("set_group_total_" + tname).c_str(), &C::template set_group_total<T>)
     .def(("get_group_total_" + tname).c_str(), &C::template get_group_total<T>);

    // recursively bind the member for each type
    bind_templated_members_Population<C, Base, Ts...>(c);
}

/*
 * @brief bind Populations class template for any choice of categories
 */
template<class... Cats>
void bind_Population(py::module& m, std::string const& name)
{
    using C = epi::Populations<Cats...>;
    using Base = epi::CustomIndexArray<epi::UncertainValue, Cats...>;
    py::class_<C, Base> c(m, name.c_str());
    c.def(py::init([](epi::Index<Cats...> const& sizes, double val){ return C(sizes, val); }))
        .def(py::init([](epi::Index<Cats...> const& sizes){ return C(sizes); }))
        .def("get_num_compartments", &C::get_num_compartments)
        .def("get_compartments", &C::get_compartments)
        .def("get_total", &C::get_total)
        .def("set_total", &C::set_total)
        .def("set_difference_from_total", &C::set_difference_from_total);

        //get_group_total, set_group_total and set_difference_from_group_total
        bind_templated_members_Population<C, Base, Cats...>(c);
}

template <class ParameterSet>
void bind_ParameterSet(py::module& m, std::string const& name)
{
    py::class_<ParameterSet> c(m, name.c_str());
    epi::foreach_tag<ParameterSet>([&c](auto t) {
        using Tag = decltype(t);

        //CAUTION: This requires ParameterTag::name() to be unique within the ParameterSet
        c.def_property(
            Tag::name().c_str(), [](const ParameterSet& self) -> auto& { return self.template get<Tag>(); },
            [](ParameterSet& self, typename Tag::Type const& v) {
                self.template get<Tag>() = v;
            },
            py::return_value_policy::reference_internal);
    });
}

/*
 * @brief bind a compartmental model for any Populations and Parameters class
 */
template<class Populations, class Parameters>
void bind_CompartmentalModel(py::module& m, std::string const& name)
{
    using Model = epi::CompartmentalModel<Populations, Parameters>;
    py::class_<Model>(m, name.c_str())
            .def(py::init<Populations const&, Parameters const&>())
            .def("apply_constraints", &Model::template apply_constraints<Parameters>)
            .def("check_constraints", &Model::template check_constraints<Parameters>)
            .def("get_initial_values", &Model::get_initial_values)
            .def_property("populations",
                [](const Model& self) -> auto& { return self.populations; },
                [](Model& self, Populations& p) { self.populations = p; })
            .def_property("parameters",
                [](const Model& self) -> auto& { return self.parameters; },
                [](Model& self, Parameters& p) { self.parameters = p; });
}


/*
 * @brief bind Simulation for any number model
 */
template <class Simulation>
void bind_Simulation(py::module& m, std::string const& name)
{
    py::class_<Simulation>(m, name.c_str())
        .def(py::init<const typename Simulation::Model&, double, double>(), py::arg("model"), py::arg("t0") = 0,
             py::arg("dt") = 0.1)
        .def_property_readonly("result", py::overload_cast<>(&Simulation::get_result, py::const_),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("model", py::overload_cast<>(&Simulation::get_model, py::const_),
                               py::return_value_policy::reference_internal)
        .def("advance", &Simulation::advance, py::arg("tmax"));
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

template <typename Simulation>
void bind_SecirSimulationNode(py::module& m, std::string const& name)
{
    py::class_<epi::Node<epi::SimulationNode<Simulation>>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const epi::Node<Simulation>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property",
            [](const epi::Node<epi::SimulationNode<Simulation>>& self) -> auto& { return self.property.get_simulation(); },
            py::return_value_policy::reference_internal);
}

/*
 * @brief bind Graph for any node and edge type
 */
template<class Model>
void bind_SecirModelGraph(py::module& m, std::string const& name)
{
    using G = epi::Graph<Model, epi::MigrationParameters>;
    py::class_<G>(m, name.c_str())
            .def(py::init<>())
            .def("add_node", &G::template add_node<const Model&>,
                 py::return_value_policy::reference_internal)
            .def("add_edge", &G::template add_edge<const epi::MigrationParameters&>,
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

template <class Simulation>
void bind_MigrationGraph(py::module& m, std::string const& name)
{
    using G = epi::Graph<epi::SimulationNode<Simulation>, epi::MigrationEdge>;
    py::class_<G>(m, name.c_str())
        .def(py::init<>())
        .def(
            "add_node", [](G & self, int id, const typename Simulation::Model& p, double t0, double dt) -> auto& {
                return self.add_node(id, p, t0, dt);
            },
            py::arg("id"), py::arg("model"), py::arg("t0") = 0.0, py::arg("dt") = 0.1,
            py::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const epi::MigrationEdge&>, py::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const Eigen::VectorXd&>, py::return_value_policy::reference_internal)
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
            "get_out_edge",
            [](const G& self, size_t node_idx, size_t edge_idx) -> auto& { return self.out_edges(node_idx)[edge_idx]; },
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
        .def_property_readonly(
            "graph", [](const GS& self) -> const Graph& { return self.get_graph(); },
            py::return_value_policy::reference_internal)
        .def_property_readonly("t", &GS::get_t)
        .def("advance", &GS::advance, py::arg("tmax"));
}

/*
 * @brief bind ParameterStudy for any model
 */
template <class Simulation>
void bind_ParameterStudy(py::module& m, std::string const& name)
{
    py::class_<epi::ParameterStudy<Simulation>>(m, name.c_str())
        .def(py::init<const typename Simulation::Model&, double, double, size_t>(), py::arg("model"), py::arg("t0"),
             py::arg("tmax"), py::arg("num_runs"))
        .def(py::init<const typename Simulation::Model&, double, double, double, size_t>(), py::arg("model"),
             py::arg("t0"), py::arg("tmax"), py::arg("dev_rel"), py::arg("num_runs"))
        .def(py::init<const epi::Graph<typename Simulation::Model, epi::MigrationParameters>&, double, double, double,
                      size_t>(),
             py::arg("model_graph"), py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("num_runs"))
        .def_property("num_runs", &epi::ParameterStudy<Simulation>::get_num_runs,
                      &epi::ParameterStudy<Simulation>::set_num_runs)
        .def_property("tmax", &epi::ParameterStudy<Simulation>::get_tmax, &epi::ParameterStudy<Simulation>::set_tmax)
        .def_property("t0", &epi::ParameterStudy<Simulation>::get_t0, &epi::ParameterStudy<Simulation>::set_t0)
        .def_property_readonly("model", py::overload_cast<>(&epi::ParameterStudy<Simulation>::get_model),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("model", py::overload_cast<>(&epi::ParameterStudy<Simulation>::get_model, py::const_),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               py::overload_cast<>(&epi::ParameterStudy<Simulation>::get_secir_model_graph),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               py::overload_cast<>(&epi::ParameterStudy<Simulation>::get_secir_model_graph, py::const_),
                               py::return_value_policy::reference_internal)
        .def(
            "run",
            [](epi::ParameterStudy<Simulation>& self, std::function<void(const epi::Graph<epi::SimulationNode<Simulation>, epi::MigrationEdge>&)> handle_result) {
                self.run([&handle_result](auto&& g) { handle_result(g); });
            },
            py::arg("handle_result_func"))
        .def("run",
             [](epi::ParameterStudy<Simulation>& self) { //default argument doesn't seem to work with functions
                 return self.run();
             })
        .def(
            "run_single",
            [](epi::ParameterStudy<Simulation>& self, std::function<void(Simulation)> handle_result) {
                self.run([handle_result](auto r) {
                    handle_result(r.nodes()[0].property.get_simulation());
                });
            },
            py::arg("handle_result_func"))
        .def("run_single", [](epi::ParameterStudy<Simulation>& self) {
            return filter_graph_results(self.run());
        });
}

/**
 * bind a constructor that has variable number of matrix shape arguments.
 * same number of arguments as the constructor of Shape.
 * @tparam C instance of pybind class_
 * @tparam ArgTuples tuples of string and some other type.
 * @param cl class_ that the constructor is defined for.
 * @param arg_tuples tuples that define additional arguments before the shape arguments.
 *                   tuples (s, t) where s is a string, the name of the argument, and t
 *                   is a value of Type T, where T is the type of the argument.
 */
template <class C, class... ArgTuples,
          class = std::enable_if_t<(std::is_same<typename C::type::Shape, epi::SquareMatrixShape>::value ||
                                    std::is_same<typename C::type::Shape, epi::ColumnVectorShape>::value),
                                   void>>
void bind_shape_constructor(C& cl, ArgTuples... arg_tuples)
{
    cl.def(py::init<Eigen::Index, std::tuple_element_t<1, ArgTuples>...>(), py::arg(std::get<0>(arg_tuples))...,
           py::arg("size"));
}

/**
 * binds a property that returns the shape of matrix valued object.
 * @tparam C instance of pybind class_.
 * @param cl class that the property is added to.
 */
template <class C>
void bind_shape_property(C& cl)
{
    cl.def_property_readonly("shape", [](typename C::type& self) {
        auto tup = py::tuple(2);
        tup[0]   = self.get_shape().rows();
        tup[1]   = self.get_shape().cols();
        return tup;
    });
}

/**
 * binds all members for an instance of epi::Damping.
 * @tparam DampingClass instance of pybind class_.
 * @param damping_class class that the members are added to.
 */
template <class DampingClass>
void bind_damping_members(DampingClass& damping_class)
{
    using Damping = typename DampingClass::type;
    using Matrix  = typename Damping::Matrix;
    using Shape   = typename Damping::Shape;

    bind_shape_constructor(damping_class);
    bind_shape_property(damping_class);

    damping_class
        .def(py::init([](const Eigen::Ref<const Matrix>& c, double t, int level, int type) {
                 return Damping(c, epi::DampingLevel(level), epi::DampingType(type), epi::SimulationTime(t));
             }),
             py::arg("coeffs"), py::arg("t"), py::arg("level") = 0, py::arg("type") = 0)
        .def_property(
            "coeffs", [](const Damping& self) -> const auto& { return self.get_coeffs(); },
            [](Damping& self, const Eigen::Ref<const Matrix>& v) {
                self.get_coeffs() = v;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "time",
            [](const Damping& self) {
                return self.get_time();
            },
            [](Damping& self, double v) {
                self.get_time() = epi::SimulationTime(v);
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "type",
            [](const Damping& self) {
                return self.get_type();
            },
            [](Damping& self, int v) {
                self.get_type() = epi::DampingType(v);
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "level",
            [](const Damping& self) {
                return self.get_level();
            },
            [](Damping& self, int v) {
                self.get_level() = epi::DampingLevel(v);
            },
            py::return_value_policy::reference_internal);
}

/**
 * binds all members for an instance of epi::Dampings.
 * @tparam DampingsClass instance of pybind class_.
 * @param dampings_class class that the members are added to.
 */
template <class DampingsClass>
void bind_dampings_members(DampingsClass& dampings_class)
{
    using Dampings = typename DampingsClass::type;
    using Damping  = typename Dampings::value_type;
    using Matrix   = typename Dampings::Matrix;
    using Shape    = typename Dampings::Shape;

    bind_shape_constructor(dampings_class);
    bind_shape_property(dampings_class);

    dampings_class
        .def("add",
             [](Dampings& self, const Damping& d) {
                 self.add(d);
             })
        .def("get_matrix_at", [](const Dampings& self, double t) {
            return self.get_matrix_at(t);
        });
}

/**
 * binds all members for an instance of epi::DampingExpression.
 * @tparam DampingExpressionClass instance of pybind class_.
 * @param damping_expression_class class that the members are added to.
 */
template <class DampingExpressionClass>
void bind_damping_expression_members(DampingExpressionClass& damping_expression_class)
{
    using DampingExpression = typename DampingExpressionClass::type;
    using Dampings          = typename DampingExpression::DampingsType;
    using Damping           = typename Dampings::value_type;
    using Matrix            = typename DampingExpression::Matrix;
    using Shape             = typename DampingExpression::Shape;


    //matrix constructors have to be defined before shape constructors.
    //otherwise 1x1 numpy matrices/vectors are converted to scalars and used as shape arguments
    damping_expression_class
        .def(py::init<const Eigen::Ref<const Matrix>&, const Eigen::Ref<const Matrix>&>(), py::arg("baseline"),
             py::arg("minimum"))
        .def(py::init<const Eigen::Ref<const Matrix>&>(), py::arg("baseline"));
    bind_shape_constructor(damping_expression_class);

    damping_expression_class
        .def("add_damping",
             [](DampingExpression& self, const Damping& d) {
                 self.add_damping(d);
             })
        .def_property(
            "baseline", [](const DampingExpression& self) -> auto& { return self.get_baseline(); },
            [](DampingExpression& self, const Eigen::Ref<const Matrix>& v) {
                self.get_baseline() = v;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "minimum", [](const DampingExpression& self) -> auto& { return self.get_minimum(); },
            [](DampingExpression& self, const Eigen::Ref<const Matrix>& v) {
                self.get_minimum() = v;
            },
            py::return_value_policy::reference_internal)
        .def("get_dampings",
             [](const DampingExpression& self) {
                 return std::vector<Damping>(self.get_dampings().begin(), self.get_dampings().end());
             })
        .def("get_matrix_at", [](const DampingExpression& self, double t) {
            return self.get_matrix_at(t);
        });
    bind_shape_property(damping_expression_class);
}

/**
 * binds all members for an instance of epi::DampingExpressionGroup.
 * @tparam DampingExpressionGroupClass instance of pybind class_.
 * @param cl class that the members are added to.
 */
template <class DampingExpressionGroupClass>
void bind_damping_expression_group_members(DampingExpressionGroupClass& cl)
{
    using DampingExpressionGroup = typename DampingExpressionGroupClass::type;
    using DampingExpression      = typename DampingExpressionGroup::value_type;
    using Dampings               = typename DampingExpression::DampingsType;
    using Damping                = typename Dampings::value_type;
    using Matrix                 = typename Damping::Matrix;
    using Shape                  = typename Damping::Shape;

    bind_shape_constructor(cl, std::make_tuple("num_matrices", size_t(0)));
    bind_shape_property(cl);

    cl.def("add_damping",
           [](DampingExpressionGroup& self, const Damping& d) {
               self.add_damping(d);
           })
        .def_property_readonly("num_matrices",
                               [](const DampingExpressionGroup& self) {
                                   return self.get_num_matrices();
                               })
        .def(
            "__getitem__", [](DampingExpressionGroup & self, size_t i) -> auto& {
                if (i < 0 || i >= self.get_num_matrices()) {
                    throw py::index_error("index out of range");
                }
                return self[i];
            },
            py::return_value_policy::reference_internal)
        .def("__setitem__",
             [](DampingExpressionGroup& self, size_t i, const DampingExpression& m) {
                 if (i < 0 && i >= self.get_num_matrices()) {
                     throw py::index_error("index out of range");
                 }
                 self[i] = m;
             })
        .def("get_matrix_at", [](const DampingExpressionGroup& self, double t) {
            return self.get_matrix_at(t);
        });
}

} // namespace

PYBIND11_MODULE(_secir, m)
{
    py::class_<epi::Date>(m, "Date")
        .def(py::init<int, int, int>(), py::arg("year"), py::arg("month"), py::arg("day"))
        .def_readwrite("year", &epi::Date::year)
        .def_readwrite("month", &epi::Date::month)
        .def_readwrite("day", &epi::Date::day)
        .def(py::self == py::self)
        .def(py::self != py::self);

    auto damping_class = py::class_<epi::SquareDamping>(m, "Damping");
    bind_damping_members(damping_class);

    auto dampings_class = py::class_<epi::SquareDampings>(m, "Dampings");
    bind_dampings_members(dampings_class);

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

    auto contact_matrix_class = py::class_<epi::ContactMatrix>(m, "ContactMatrix");
    bind_damping_expression_members(contact_matrix_class);
    contact_matrix_class.def_property_readonly("num_groups", &epi::ContactMatrix::get_num_groups);

    auto contact_matrix_group_class = py::class_<epi::ContactMatrixGroup>(m, "ContactMatrixGroup");
    bind_damping_expression_group_members(contact_matrix_group_class);
    contact_matrix_group_class.def_property_readonly("num_groups", &epi::ContactMatrixGroup::get_num_groups);

    py::class_<epi::DampingSampling>(m, "DampingSampling")
        .def(py::init([](const epi::UncertainValue& value, int level, int type, double time,
                         const std::vector<size_t>& matrices, const Eigen::Ref<const Eigen::VectorXd>& groups) {
                 return epi::DampingSampling(value, epi::DampingLevel(level), epi::DampingType(type),
                                             epi::SimulationTime(time), matrices, groups);
             }),
             py::arg("value"), py::arg("level"), py::arg("type"), py::arg("time"), py::arg("matrix_indices"),
             py::arg("group_weights"))
        .def_property("value", py::overload_cast<>(&epi::DampingSampling::get_value), &epi::DampingSampling::set_value,
                      py::return_value_policy::reference_internal)
        .def_property(
            "level",
            [](const epi::DampingSampling& self) {
                return int(self.get_level());
            },
            [](epi::DampingSampling& self, int lvl) {
                self.set_level(epi::DampingLevel(lvl));
            })
        .def_property(
            "type",
            [](const epi::DampingSampling& self) {
                return int(self.get_type());
            },
            [](epi::DampingSampling& self, int typ) {
                self.set_type(epi::DampingType(typ));
            })
        .def_property(
            "time",
            [](const epi::DampingSampling& self) {
                return double(self.get_time());
            },
            [](epi::DampingSampling& self, double t) {
                self.set_time(epi::SimulationTime(t));
            })
        .def_property("matrix_indices", &epi::DampingSampling::get_matrix_indices,
                      &epi::DampingSampling::set_matrix_indices)
        .def_property(
            "group_weights", &epi::DampingSampling::get_group_weights,
            [](epi::DampingSampling& self, const Eigen::Ref<const Eigen::VectorXd>& v) {
                self.set_group_weights(v);
            },
            py::return_value_policy::reference_internal);

    py::class_<epi::UncertainContactMatrix>(m, "UncertainContactMatrix")
        .def(py::init<>())
        .def(py::init<const epi::ContactMatrixGroup&>())
        .def_property(
            "cont_freq_mat", py::overload_cast<>(&epi::UncertainContactMatrix::get_cont_freq_mat),
            [](epi::UncertainContactMatrix& self, const epi::ContactMatrixGroup& c) {
                self.get_cont_freq_mat() = c;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "dampings", py::overload_cast<>(&epi::UncertainContactMatrix::get_dampings),
            [](epi::UncertainContactMatrix& self, const std::vector<epi::DampingSampling>& v) {
                self.get_dampings() = v;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "school_holidays",
            [](const epi::UncertainContactMatrix& self) {
                std::vector<std::pair<double, double>> v(self.get_school_holidays().size());
                std::transform(self.get_school_holidays().begin(), self.get_school_holidays().end(), v.begin(),
                               [](auto& p) {
                                   return std::make_pair(double(p.first), double(p.second));
                               });
                return v;
            },
            [](epi::UncertainContactMatrix& self, const std::vector<std::pair<double, double>>& v) {
                self.get_school_holidays().resize(v.size());
                std::transform(v.begin(), v.end(), self.get_school_holidays().begin(), [](auto& p) {
                    return std::make_pair(epi::SimulationTime(p.first), epi::SimulationTime(p.second));
                });
            })
        .def_property("school_holiday_damping",
                      py::overload_cast<>(&epi::UncertainContactMatrix::get_school_holiday_damping),
                      [](epi::UncertainContactMatrix& self, const epi::DampingSampling& v) {
                          self.get_school_holiday_damping() = v;
                      });

    auto migration_damping_class = py::class_<epi::VectorDamping>(m, "MigrationDamping");
    bind_damping_members(migration_damping_class);

    auto migration_dampings_class = py::class_<epi::VectorDampings>(m, "MigrationDampings");
    bind_dampings_members(migration_dampings_class);

    auto migration_coeffs_class = py::class_<epi::MigrationCoefficients>(m, "MigrationCoefficients");
    bind_damping_expression_members(migration_coeffs_class);

    auto migration_coeff_group_class = py::class_<epi::MigrationCoefficientGroup>(m, "MigrationCoefficientGroup");
    bind_damping_expression_group_members(migration_coeff_group_class);

    py::class_<epi::MigrationParameters>(m, "MigrationParameters")
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def(py::init<const epi::MigrationCoefficientGroup&>(), py::arg("coeffs"))
        .def_property(
            "coefficients", py::overload_cast<>(&epi::MigrationParameters::get_coefficients),
            [](epi::MigrationParameters& self, const epi::MigrationCoefficientGroup& v) {
                self.get_coefficients() = v;
            },
            py::return_value_policy::reference_internal);

    py::class_<epi::Edge<epi::MigrationParameters>>(m, "MigrationParameterEdge")
        .def_property_readonly("start_node_idx",
                               [](const epi::Edge<epi::MigrationParameters>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const epi::Edge<epi::MigrationParameters>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property", [](const epi::Edge<epi::MigrationEdge>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);

    py::class_<epi::MigrationEdge>(m, "Migration")
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def(py::init<const epi::MigrationParameters&>(), py::arg("params"))
        .def_property_readonly(
            "parameters", [](const epi::MigrationEdge& self) -> auto& { return self.get_parameters(); },
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
    m.def("interpolate_simulation_result", static_cast<epi::TimeSeries<double> (*)(const epi::TimeSeries<double>&)>(
                                               &epi::interpolate_simulation_result));

    m.def("interpolate_ensemble_results", &epi::interpolate_ensemble_results<epi::TimeSeries<double>>);

    m.def("ensemble_mean", &epi::ensemble_mean);
    m.def("ensemble_percentile", &epi::ensemble_percentile);

    py::enum_<epi::InfectionState>(m, "InfectionState")
        .value("Susceptible", epi::InfectionState::Susceptible)
        .value("Exposed", epi::InfectionState::Exposed)
        .value("Carrier", epi::InfectionState::Carrier)
        .value("Infected", epi::InfectionState::Infected)
        .value("Hospitalized", epi::InfectionState::Hospitalized)
        .value("ICU", epi::InfectionState::ICU)
        .value("Recovered", epi::InfectionState::Recovered)
        .value("Dead", epi::InfectionState::Dead)
        .value("Count", epi::InfectionState::Count)
        .export_values();

    bind_Index<epi::InfectionState>(m, "Index_InfectionState");
    bind_Index<epi::AgeGroup>(m, "Index_AgeGroup");

    py::class_<epi::AgeGroup, epi::Index<epi::AgeGroup>>(m, "AgeGroup")
        .def(py::init<size_t>());

    bind_MultiIndex<epi::AgeGroup, epi::InfectionState>(m, "Index_Agegroup_InfectionState");
    bind_CustomIndexArray<epi::UncertainValue, epi::AgeGroup, epi::InfectionState>(m, "SecirPopulationArray");
    bind_CustomIndexArray<epi::UncertainValue, epi::AgeGroup>(m, "AgeGroupArray");

    bind_ParameterSet<epi::SecirParamsBase>(m, "SecirParamsBase");

    py::class_<epi::SecirParams, epi::SecirParamsBase>(m, "SecirParams")
        .def(py::init<epi::AgeGroup>())
        .def("check_constraints", &epi::SecirParams::check_constraints)
        .def("apply_constraints", &epi::SecirParams::apply_constraints);

    bind_Population<epi::AgeGroup, epi::InfectionState>(m, "SecirPopulation");

    using SecirPopulations = epi::Populations<epi::AgeGroup, epi::InfectionState>;
    bind_CompartmentalModel<SecirPopulations, epi::SecirParams>(m, "SecirModelBase");
    py::class_<epi::SecirModel, epi::CompartmentalModel<SecirPopulations, epi::SecirParams>>(m, "SecirModel")
            .def(py::init<size_t>(), py::arg("num_agegroups"));

    bind_Simulation<epi::SecirSimulation<>>(m, "SecirSimulation");

    m.def("simulate", [](double t0, double tmax, double dt, const epi::SecirModel& model) { return epi::simulate(t0, tmax, dt, model); },
          "Simulates a SecirModel1 from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
          py::arg("model"));

    bind_SecirModelNode<epi::SecirModel>(m, "SecirModelNode");
    bind_SecirSimulationNode<epi::SecirSimulation<>>(m, "SecirSimulationNode");
    bind_SecirModelGraph<epi::SecirModel>(m, "SecirModelGraph");
    using Simulation = epi::SecirSimulation<>;
    bind_MigrationGraph<Simulation>(m, "MigrationGraph");
    using MigrationGraph = epi::Graph<epi::SimulationNode<Simulation>, epi::MigrationEdge>;
    bind_GraphSimulation<MigrationGraph>(m, "MigrationSimulation");

    bind_ParameterStudy<epi::SecirSimulation<>>(m, "ParameterStudy");

    m.def("set_params_distributions_normal", &epi::set_params_distributions_normal, py::arg("model"), py::arg("t0"),
          py::arg("tmax"), py::arg("dev_rel"));

    m.def("draw_sample", &epi::draw_sample, py::arg("model"));

    m.def("interpolate_simulation_result",
          py::overload_cast<const MigrationGraph&>(&epi::interpolate_simulation_result<Simulation>));

    m.def("interpolate_ensemble_results", &epi::interpolate_ensemble_results<MigrationGraph>);

    m.def(
        "get_state_id_de",
        [](int county) {
            return int(epi::regions::de::get_state_id(epi::regions::de::CountyId(county)));
        },
        py::arg("county_id"));
    m.def(
        "get_holidays_de",
        [](int state, epi::Date start_date, epi::Date end_date) {
            auto h = epi::regions::de::get_holidays(epi::regions::de::StateId(state), start_date, end_date);
            return std::vector<std::pair<epi::Date, epi::Date>>(h.begin(), h.end());
        },
        py::arg("state_id"), py::arg("start_date") = epi::Date(std::numeric_limits<int>::min(), 1, 1),
        py::arg("end_date") = epi::Date(std::numeric_limits<int>::max(), 1, 1));

    m.attr("__version__") = "dev";
}
