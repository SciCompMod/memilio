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


#include <memilio/utils/custom_index_array.h>
#include <memilio/utils/time_series.h>
#include <memilio/epidemiology/damping.h>
#include <memilio/epidemiology/regions.h>
#include <secir/secir.h>
#include <secir/parameter_studies.h>
#include <secir/analyze_result.h>
#include <pickle_serializer.h>

#include <Eigen/Core>
#include <vector>

namespace py = pybind11;

namespace
{

//select only the first node of the graph of each run, used for parameterstudy with single nodes
template<class Sim>
std::vector<Sim> filter_graph_results(
    const std::vector<mio::Graph<mio::SimulationNode<Sim>, mio::MigrationEdge>>& graph_results)
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
    py::class_<mio::Index<Tag>> c(m, name.c_str());
    c.def(py::init<size_t>(), py::arg("value"));
}

// helper function for implicitly casting from py::tuple to Index in Python.
// This extracts an Index from a py::tuple of Indices from the correct position,
// given the corresponding Index type
template <typename Tag, class Tuple>
mio::Index<Tag> extract_index(py::tuple& t)
{
    return t[mio::details::IndexPosition<Tag, Tuple>::value].template cast<mio::Index<Tag>>();
}

// bind an index for more than one tag
template <class... Tags>
void bind_MultiIndex(py::module& m, std::string const& name)
{
    using C = mio::Index<Tags...>;
    py::class_<C> c(m, name.c_str());
    c.def(py::init<mio::Index<Tags> const&...>()).def(py::init([](py::tuple t) {
        return C(extract_index<Tags, C>(t)...);
    }));

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
template <>
std::string pretty_name<mio::InfectionState>()
{
    return "InfectionState";
}
template <>
std::string pretty_name<mio::AgeGroup>()
{
    return "AgeGroup";
}

template <class C>
void bind_templated_members_CustomIndexArray(py::class_<C>&)
{
}

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
    using C     = typename mio::CustomIndexArray<Type, Tags...>;
    using Index = typename mio::CustomIndexArray<Type, Tags...>::Index;
    py::class_<C> c(m, name.c_str());
    c.def(py::init([](Index const& sizes, Type const& val) {
         return C(sizes, val);
     }))
        .def(py::init([](Index const& sizes) {
            return C(sizes);
        }))
        .def("numel", &C::numel)
        .def(
            "__getitem__", [](const C& self, Index const& idx) -> auto& { return self[idx]; },
            py::return_value_policy::reference_internal)
        .def("__setitem__",
             [](C& self, Index const& idx, double value) {
                 self[idx] = value;
             })
        .def(
            "__iter__",
            [](const C& s) {
                return py::make_iterator(s.begin(), s.end());
            },
            py::keep_alive<0, 1>())
        .def("get_flat_index", &C::get_flat_index);

    // Not supported in Python yet: Slicing

    // bind all templated members for types in Tags...
    bind_templated_members_CustomIndexArray<C, Tags...>(c);
}

template <class C, class Base>
void bind_templated_members_Population(py::class_<C, Base>&)
{
}

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
template <class... Cats>
void bind_Population(py::module& m, std::string const& name)
{
    using C    = mio::Populations<Cats...>;
    using Base = mio::CustomIndexArray<mio::UncertainValue, Cats...>;
    py::class_<C, Base> c(m, name.c_str());
    c.def(py::init([](mio::Index<Cats...> const& sizes, double val) {
         return C(sizes, val);
     }))
        .def(py::init([](mio::Index<Cats...> const& sizes) {
            return C(sizes);
        }))
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
    mio::foreach_tag<ParameterSet>([&c](auto t) {
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
template <class Populations, class Parameters>
void bind_CompartmentalModel(py::module& m, std::string const& name)
{
    using Model = mio::CompartmentalModel<Populations, Parameters>;
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
    py::class_<mio::Node<Model>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const mio::Node<Model>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property", [](const mio::Node<Model>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);
}

template <typename Simulation>
void bind_SecirSimulationNode(py::module& m, std::string const& name)
{
    py::class_<mio::Node<mio::SimulationNode<Simulation>>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const mio::Node<Simulation>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property",
            [](const mio::Node<mio::SimulationNode<Simulation>>& self) -> auto& { return self.property.get_simulation(); },
            py::return_value_policy::reference_internal);
}

/*
 * @brief bind Graph for any node and edge type
 */
template <class Model>
void bind_SecirModelGraph(py::module& m, std::string const& name)
{
    using G = mio::Graph<Model, mio::MigrationParameters>;
    py::class_<G>(m, name.c_str())
        .def(py::init<>())
        .def("add_node", &G::template add_node<const Model&>, py::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const mio::MigrationParameters&>,
             py::return_value_policy::reference_internal)
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

template <class Simulation>
void bind_MigrationGraph(py::module& m, std::string const& name)
{
    using G = mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>;
    py::class_<G>(m, name.c_str())
        .def(py::init<>())
        .def(
            "add_node", [](G & self, int id, const typename Simulation::Model& p, double t0, double dt) -> auto& {
                return self.add_node(id, p, t0, dt);
            },
            py::arg("id"), py::arg("model"), py::arg("t0") = 0.0, py::arg("dt") = 0.1,
            py::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const mio::MigrationEdge&>, py::return_value_policy::reference_internal)
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
template <class Graph>
void bind_GraphSimulation(py::module& m, std::string const& name)
{
    using GS = mio::GraphSimulation<Graph>;
    py::class_<GS>(m, name.c_str())
        .def(py::init([](const Graph& graph, double t0, double dt) {
                 return std::make_unique<GS>(mio::make_migration_sim(t0, dt, graph));
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
    py::class_<mio::ParameterStudy<Simulation>>(m, name.c_str())
        .def(py::init<const typename Simulation::Model&, double, double, size_t>(), py::arg("model"), py::arg("t0"),
             py::arg("tmax"), py::arg("num_runs"))
        .def(py::init<const typename Simulation::Model&, double, double, double, size_t>(), py::arg("model"),
             py::arg("t0"), py::arg("tmax"), py::arg("dev_rel"), py::arg("num_runs"))
        .def(py::init<const mio::Graph<typename Simulation::Model, mio::MigrationParameters>&, double, double, double,
                      size_t>(),
             py::arg("model_graph"), py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("num_runs"))
        .def_property("num_runs", &mio::ParameterStudy<Simulation>::get_num_runs,
                      &mio::ParameterStudy<Simulation>::set_num_runs)
        .def_property("tmax", &mio::ParameterStudy<Simulation>::get_tmax, &mio::ParameterStudy<Simulation>::set_tmax)
        .def_property("t0", &mio::ParameterStudy<Simulation>::get_t0, &mio::ParameterStudy<Simulation>::set_t0)
        .def_property_readonly("model", py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_model),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("model", py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_model, py::const_),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_secir_model_graph),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_secir_model_graph, py::const_),
                               py::return_value_policy::reference_internal)
        .def(
            "run",
            [](mio::ParameterStudy<Simulation>& self, std::function<void(const mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>&)> handle_result) {
                self.run([&handle_result](auto&& g) { handle_result(g); });
            },
            py::arg("handle_result_func"))
        .def("run",
             [](mio::ParameterStudy<Simulation>& self) { //default argument doesn't seem to work with functions
                 return self.run();
             })
        .def(
            "run_single",
            [](mio::ParameterStudy<Simulation>& self, std::function<void(Simulation)> handle_result) {
                self.run([handle_result](auto r) {
                    handle_result(r.nodes()[0].property.get_simulation());
                });
            },
            py::arg("handle_result_func"))
        .def("run_single", [](mio::ParameterStudy<Simulation>& self) {
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
          class = std::enable_if_t<(std::is_same<typename C::type::Shape, mio::SquareMatrixShape>::value ||
                                    std::is_same<typename C::type::Shape, mio::ColumnVectorShape>::value),
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
 * binds all members for an instance of mio::Damping.
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
                 return Damping(c, mio::DampingLevel(level), mio::DampingType(type), mio::SimulationTime(t));
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
                self.get_time() = mio::SimulationTime(v);
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "type",
            [](const Damping& self) {
                return self.get_type();
            },
            [](Damping& self, int v) {
                self.get_type() = mio::DampingType(v);
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "level",
            [](const Damping& self) {
                return self.get_level();
            },
            [](Damping& self, int v) {
                self.get_level() = mio::DampingLevel(v);
            },
            py::return_value_policy::reference_internal);
}

/**
 * binds all members for an instance of mio::Dampings.
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
 * binds all members for an instance of mio::DampingExpression.
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
 * binds all members for an instance of mio::DampingExpressionGroup.
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

template <class T, class ... Args>
decltype(auto) pybind_pickle_class(py::module &m, const char* name)
{
    decltype(auto) pickle_class = py::class_<T, Args...>(m, name);
    pickle_class.def(py::pickle(
             [](const T &object) { // __getstate__
                auto tuple = mio::serialize_pickle(object);
                if (tuple)
                {
                    return std::move(tuple).value();
                }
                else
                {
                    throw std::runtime_error(tuple.error().formatted_message());
                }
            },
            [](const py::tuple t) { // __setstate__

                auto object = mio::deserialize_pickle(t,mio::Tag<T>{});
                if (object)
                {
                    return std::move(object).value();
                }
                else
                {
                    throw std::runtime_error(object.error().formatted_message());
                }
            }
    ));
    return pickle_class;
}

PYBIND11_MODULE(_secir, m)
{
    pybind_pickle_class<mio::Date>(m, "Date")
        .def(py::init<int, int, int>(), py::arg("year"), py::arg("month"), py::arg("day"))
        .def_readwrite("year", &mio::Date::year)
        .def_readwrite("month", &mio::Date::month)
        .def_readwrite("day", &mio::Date::day)
        .def(py::self == py::self)
        .def(py::self != py::self);

    auto damping_class = py::class_<mio::SquareDamping>(m, "Damping");
    bind_damping_members(damping_class);

    auto dampings_class = py::class_<mio::SquareDampings>(m, "Dampings");
    bind_dampings_members(dampings_class);

    py::class_<mio::TimeSeries<double>>(m, "TimeSeries")
        .def(py::init<Eigen::Index>(), py::arg("num_elements"))
        .def("get_num_time_points", &mio::TimeSeries<double>::get_num_time_points)
        .def("get_num_elements", &mio::TimeSeries<double>::get_num_elements)
        .def("get_time", py::overload_cast<Eigen::Index>(&mio::TimeSeries<double>::get_time), py::arg("index"))
        .def("get_last_time", py::overload_cast<>(&mio::TimeSeries<double>::get_last_time))
        .def("get_value", py::overload_cast<Eigen::Index>(&mio::TimeSeries<double>::get_value), py::arg("index"))
        .def("get_last_value", py::overload_cast<>(&mio::TimeSeries<double>::get_last_value))
        .def(
            "__getitem__",
            [](mio::TimeSeries<double>& self, Eigen::Index i) {
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
            [](mio::TimeSeries<double>& self, Eigen::Index i, Eigen::Ref<const mio::TimeSeries<double>::Vector> expr) {
                if (i >= 0 && i < self.get_num_time_points()) {
                    self[i] = expr;
                }
                else {
                    throw pybind11::index_error("Index out of range."); //needs to throw exception for iterable
                }
            },
            py::is_operator(), py::arg("index"), py::arg("v"))
        .def("add_time_point",
             [](mio::TimeSeries<double>& self) {
                 return self.add_time_point();
             })
        .def("add_time_point",
             [](mio::TimeSeries<double>& self, double t) {
                 return self.add_time_point(t);
             })
        .def("add_time_point",
             [](mio::TimeSeries<double>& self, double t, Eigen::Ref<const mio::TimeSeries<double>::Vector> expr) {
                 return self.add_time_point(t, expr);
             })
        .def("as_ndarray", [](mio::TimeSeries<double>& self) {
            auto m = Eigen::Map<mio::TimeSeries<double>::Matrix>(self.data(), self.get_num_rows(),
                                                                 self.get_num_time_points());
            return Eigen::Ref<mio::TimeSeries<double>::Matrix>(m);
        });

    py::class_<mio::ParameterDistribution>(m, "ParameterDistribution")
        .def_property("lower_bound", &mio::ParameterDistribution::get_lower_bound,
                      &mio::ParameterDistribution::set_lower_bound)
        .def_property("upper_bound", &mio::ParameterDistribution::get_upper_bound,
                      &mio::ParameterDistribution::set_upper_bound)
        .def("add_predefined_sample", &mio::ParameterDistribution::add_predefined_sample)
        .def("remove_predefined_samples", &mio::ParameterDistribution::remove_predefined_samples)
        .def("get_sample", &mio::ParameterDistribution::get_sample);

    pybind_pickle_class<mio::ParameterDistributionNormal, mio::ParameterDistribution>(m, "ParameterDistributionNormal")
        .def(py::init<double, double, double, double>(), py::arg("lb"), py::arg("ub"), py::arg("mean"),
             py::arg("std_dev"))
        .def(py::init<double, double, double>(), py::arg("lb"), py::arg("ub"), py::arg("mean"))
        .def_property("mean", &mio::ParameterDistributionNormal::get_mean, &mio::ParameterDistributionNormal::set_mean)
        .def_property("standard_dev", &mio::ParameterDistributionNormal::get_standard_dev,
                      &mio::ParameterDistributionNormal::set_standard_dev);

    pybind_pickle_class<mio::ParameterDistributionUniform, mio::ParameterDistribution>(m, "ParameterDistributionUniform")
        .def(py::init<>())
        .def(py::init<double, double>(), py::arg("lb"), py::arg("ub"));

    pybind_pickle_class<mio::UncertainValue>(m, "UncertainValue")
        .def(py::init<ScalarType>(), py::arg("value") = 0.0)
        .def_property(
            "value",
            [](mio::UncertainValue& self) {
                return ScalarType(self);
            },
            [](mio::UncertainValue& self, ScalarType v) {
                self = v;
            })
        .def("set_distribution", //a property would be nicer but getter and setter use a different type
             &mio::UncertainValue::set_distribution)
        .def(
            "get_distribution",
            [](const mio::UncertainValue& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def(
            "get_distribution",
            [](mio::UncertainValue& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def("draw_sample", &mio::UncertainValue::draw_sample);

    auto contact_matrix_class = py::class_<mio::ContactMatrix>(m, "ContactMatrix");
    bind_damping_expression_members(contact_matrix_class);
    contact_matrix_class.def_property_readonly("num_groups", &mio::ContactMatrix::get_num_groups);

    auto contact_matrix_group_class = py::class_<mio::ContactMatrixGroup>(m, "ContactMatrixGroup");
    bind_damping_expression_group_members(contact_matrix_group_class);
    contact_matrix_group_class.def_property_readonly("num_groups", &mio::ContactMatrixGroup::get_num_groups);

    pybind_pickle_class<mio::DampingSampling>(m, "DampingSampling")
        .def(py::init([](const mio::UncertainValue& value, int level, int type, double time,
                         const std::vector<size_t>& matrices, const Eigen::Ref<const Eigen::VectorXd>& groups) {
                 return mio::DampingSampling(value, mio::DampingLevel(level), mio::DampingType(type),
                                             mio::SimulationTime(time), matrices, groups);
             }),
             py::arg("value"), py::arg("level"), py::arg("type"), py::arg("time"), py::arg("matrix_indices"),
             py::arg("group_weights"))
        .def_property("value", py::overload_cast<>(&mio::DampingSampling::get_value), &mio::DampingSampling::set_value,
                      py::return_value_policy::reference_internal)
        .def_property(
            "level",
            [](const mio::DampingSampling& self) {
                return int(self.get_level());
            },
            [](mio::DampingSampling& self, int lvl) {
                self.set_level(mio::DampingLevel(lvl));
            })
        .def_property(
            "type",
            [](const mio::DampingSampling& self) {
                return int(self.get_type());
            },
            [](mio::DampingSampling& self, int typ) {
                self.set_type(mio::DampingType(typ));
            })
        .def_property(
            "time",
            [](const mio::DampingSampling& self) {
                return double(self.get_time());
            },
            [](mio::DampingSampling& self, double t) {
                self.set_time(mio::SimulationTime(t));
            })
        .def_property("matrix_indices", &mio::DampingSampling::get_matrix_indices,
                      &mio::DampingSampling::set_matrix_indices)
        .def_property(
            "group_weights", &mio::DampingSampling::get_group_weights,
            [](mio::DampingSampling& self, const Eigen::Ref<const Eigen::VectorXd>& v) {
                self.set_group_weights(v);
            },
            py::return_value_policy::reference_internal);

    py::class_<mio::UncertainContactMatrix>(m, "UncertainContactMatrix")
        .def(py::init<>())
        .def(py::init<const mio::ContactMatrixGroup&>())
        .def_property(
            "cont_freq_mat", py::overload_cast<>(&mio::UncertainContactMatrix::get_cont_freq_mat),
            [](mio::UncertainContactMatrix& self, const mio::ContactMatrixGroup& c) {
                self.get_cont_freq_mat() = c;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "dampings", py::overload_cast<>(&mio::UncertainContactMatrix::get_dampings),
            [](mio::UncertainContactMatrix& self, const std::vector<mio::DampingSampling>& v) {
                self.get_dampings() = v;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "school_holidays",
            [](const mio::UncertainContactMatrix& self) {
                std::vector<std::pair<double, double>> v(self.get_school_holidays().size());
                std::transform(self.get_school_holidays().begin(), self.get_school_holidays().end(), v.begin(),
                               [](auto& p) {
                                   return std::make_pair(double(p.first), double(p.second));
                               });
                return v;
            },
            [](mio::UncertainContactMatrix& self, const std::vector<std::pair<double, double>>& v) {
                self.get_school_holidays().resize(v.size());
                std::transform(v.begin(), v.end(), self.get_school_holidays().begin(), [](auto& p) {
                    return std::make_pair(mio::SimulationTime(p.first), mio::SimulationTime(p.second));
                });
            })
        .def_property("school_holiday_damping",
                      py::overload_cast<>(&mio::UncertainContactMatrix::get_school_holiday_damping),
                      [](mio::UncertainContactMatrix& self, const mio::DampingSampling& v) {
                          self.get_school_holiday_damping() = v;
                      });

    auto migration_damping_class = py::class_<mio::VectorDamping>(m, "MigrationDamping");
    bind_damping_members(migration_damping_class);

    auto migration_dampings_class = py::class_<mio::VectorDampings>(m, "MigrationDampings");
    bind_dampings_members(migration_dampings_class);

    auto migration_coeffs_class = py::class_<mio::MigrationCoefficients>(m, "MigrationCoefficients");
    bind_damping_expression_members(migration_coeffs_class);

    auto migration_coeff_group_class = py::class_<mio::MigrationCoefficientGroup>(m, "MigrationCoefficientGroup");
    bind_damping_expression_group_members(migration_coeff_group_class);

    py::class_<mio::MigrationParameters>(m, "MigrationParameters")
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def(py::init<const mio::MigrationCoefficientGroup&>(), py::arg("coeffs"))
        .def_property(
            "coefficients", py::overload_cast<>(&mio::MigrationParameters::get_coefficients),
            [](mio::MigrationParameters& self, const mio::MigrationCoefficientGroup& v) {
                self.get_coefficients() = v;
            },
            py::return_value_policy::reference_internal);

    py::class_<mio::Edge<mio::MigrationParameters>>(m, "MigrationParameterEdge")
        .def_property_readonly("start_node_idx",
                               [](const mio::Edge<mio::MigrationParameters>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const mio::Edge<mio::MigrationParameters>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property", [](const mio::Edge<mio::MigrationEdge>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);

    py::class_<mio::MigrationEdge>(m, "Migration")
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def(py::init<const mio::MigrationParameters&>(), py::arg("params"))
        .def_property_readonly(
            "parameters", [](const mio::MigrationEdge& self) -> auto& { return self.get_parameters(); },
            py::return_value_policy::reference_internal);

    py::class_<mio::Edge<mio::MigrationEdge>>(m, "MigrationEdge")
        .def_property_readonly("start_node_idx",
                               [](const mio::Edge<mio::MigrationEdge>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const mio::Edge<mio::MigrationEdge>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property", [](const mio::Edge<mio::MigrationEdge>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);

    // https://github.com/pybind/pybind11/issues/1153
    m.def("interpolate_simulation_result", static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&)>(
                                               &mio::interpolate_simulation_result));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<mio::TimeSeries<double>>);

    m.def("ensemble_mean", &mio::ensemble_mean);
    m.def("ensemble_percentile", &mio::ensemble_percentile);

    py::enum_<mio::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::InfectionState::Susceptible)
        .value("Exposed", mio::InfectionState::Exposed)
        .value("Carrier", mio::InfectionState::Carrier)
        .value("Infected", mio::InfectionState::Infected)
        .value("Hospitalized", mio::InfectionState::Hospitalized)
        .value("ICU", mio::InfectionState::ICU)
        .value("Recovered", mio::InfectionState::Recovered)
        .value("Dead", mio::InfectionState::Dead)
        .value("Count", mio::InfectionState::Count)
        .export_values();

    bind_Index<mio::InfectionState>(m, "Index_InfectionState");
    bind_Index<mio::AgeGroup>(m, "Index_AgeGroup");

    py::class_<mio::AgeGroup, mio::Index<mio::AgeGroup>>(m, "AgeGroup").def(py::init<size_t>());

    bind_MultiIndex<mio::AgeGroup, mio::InfectionState>(m, "Index_Agegroup_InfectionState");
    bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup, mio::InfectionState>(m, "SecirPopulationArray");
    bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, "AgeGroupArray");

    bind_ParameterSet<mio::SecirParamsBase>(m, "SecirParamsBase");

    py::class_<mio::SecirParams, mio::SecirParamsBase>(m, "SecirParams")
        .def(py::init<mio::AgeGroup>())
        .def("check_constraints", &mio::SecirParams::check_constraints)
        .def("apply_constraints", &mio::SecirParams::apply_constraints);

    bind_Population<mio::AgeGroup, mio::InfectionState>(m, "SecirPopulation");

    using SecirPopulations = mio::Populations<mio::AgeGroup, mio::InfectionState>;
    bind_CompartmentalModel<SecirPopulations, mio::SecirParams>(m, "SecirModelBase");
    py::class_<mio::SecirModel, mio::CompartmentalModel<SecirPopulations, mio::SecirParams>>(m, "SecirModel")
        .def(py::init<size_t>(), py::arg("num_agegroups"));

    bind_Simulation<mio::SecirSimulation<>>(m, "SecirSimulation");

    m.def("simulate", [](double t0, double tmax, double dt, const mio::SecirModel& model) { return mio::simulate(t0, tmax, dt, model); },
          "Simulates a SecirModel1 from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
          py::arg("model"));

    bind_SecirModelNode<mio::SecirModel>(m, "SecirModelNode");
    bind_SecirSimulationNode<mio::SecirSimulation<>>(m, "SecirSimulationNode");
    bind_SecirModelGraph<mio::SecirModel>(m, "SecirModelGraph");
    using Simulation = mio::SecirSimulation<>;
    bind_MigrationGraph<Simulation>(m, "MigrationGraph");
    using MigrationGraph = mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>;
    bind_GraphSimulation<MigrationGraph>(m, "MigrationSimulation");

    bind_ParameterStudy<mio::SecirSimulation<>>(m, "ParameterStudy");

    m.def("set_params_distributions_normal", &mio::set_params_distributions_normal, py::arg("model"), py::arg("t0"),
          py::arg("tmax"), py::arg("dev_rel"));

    m.def("draw_sample", &mio::draw_sample, py::arg("model"));

    m.def("interpolate_simulation_result",
          py::overload_cast<const MigrationGraph&>(&mio::interpolate_simulation_result<Simulation>));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<MigrationGraph>);

    m.def(
        "get_state_id_de",
        [](int county) {
            return int(mio::regions::de::get_state_id(mio::regions::de::CountyId(county)));
        },
        py::arg("county_id"));
    m.def(
        "get_holidays_de",
        [](int state, mio::Date start_date, mio::Date end_date) {
            auto h = mio::regions::de::get_holidays(mio::regions::de::StateId(state), start_date, end_date);
            return std::vector<std::pair<mio::Date, mio::Date>>(h.begin(), h.end());
        },
        py::arg("state_id"), py::arg("start_date") = mio::Date(std::numeric_limits<int>::min(), 1, 1),
        py::arg("end_date") = mio::Date(std::numeric_limits<int>::max(), 1, 1));

    m.attr("__version__") = "dev";
}
