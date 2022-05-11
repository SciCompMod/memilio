#ifndef PYMIO_TEMPLATES_H
#define PYMIO_TEMPLATES_H

#include "memilio/utils/parameter_set.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/mobility.h"
#include "pickle_serializer.h"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#include <type_traits>

namespace pymio
{

//bind class and add pickling based on memilio serialization framework
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

// bind an index for a single tag
template <class Tag> 
void bind_Index(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::Index<Tag>> c(m, name.c_str());
    c.def(pybind11::init<size_t>(), pybind11::arg("value"));
    c.def(pybind11::self == pybind11::self);
    c.def(pybind11::self != pybind11:: self);
}

// helper function for implicitly casting from pybind11::tuple to Index in Python.
// This extracts an Index from a pybind11::tuple of Indices from the correct position,
// given the corresponding Index type
template <typename Tag, class Tuple>
mio::Index<Tag> extract_index(pybind11::tuple& t)
{
    return t[mio::details::IndexPosition<Tag, Tuple>::value].template cast<mio::Index<Tag>>();
}

// bind an index for more than one tag
template <class... Tags>
void bind_MultiIndex(pybind11::module& m, std::string const& name)
{
    using C = mio::Index<Tags...>;
    pybind11::class_<C> c(m, name.c_str());
    c.def(pybind11::init<mio::Index<Tags> const&...>()).def(pybind11::init([](pybind11::tuple t) {
        return C(extract_index<Tags, C>(t)...);
    }));

    pybind11::implicitly_convertible<pybind11::tuple, C>();
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

template <class C>
void bind_templated_members_CustomIndexArray(pybind11::class_<C>&)
{
}

template <class C, class T, class... Ts>
void bind_templated_members_CustomIndexArray(pybind11::class_<C>& c)
{
    std::string tname = pretty_name<T>();
    c.def(("size_" + tname).c_str(), &C::template size<T>);

    // recursively bind the member for each type
    bind_templated_members_CustomIndexArray<C, Ts...>(c);
}

template <class Type, class... Tags>
void bind_CustomIndexArray(pybind11::module& m, std::string const& name)
{
    using C     = typename mio::CustomIndexArray<Type, Tags...>;
    using Index = typename mio::CustomIndexArray<Type, Tags...>::Index;
    pybind11::class_<C> c(m, name.c_str());
    c.def(pybind11::init([](Index const& sizes, Type const& val) {
         return C(sizes, val);
     }))
        .def(pybind11::init([](Index const& sizes) {
            return C(sizes);
        }))
        .def("numel", &C::numel)
        .def(
            "__getitem__", [](const C& self, Index const& idx) -> auto& { return self[idx]; },
            pybind11::return_value_policy::reference_internal)
        .def(
            "__getitem__", [](const C& self, std::tuple<mio::Index<Tags>...> idx) -> auto& { //python natively handles multi-indices as tuples
                return self[{std::get<mio::Index<Tags>>(idx)...}];
            }, pybind11::return_value_policy::reference_internal)
        .def("__setitem__",
             [](C& self, Index const& idx, double value) {
                 self[idx] = value;
             })
        .def("__setitem__",
             [](C& self, std::tuple<mio::Index<Tags>...> idx, double value) {
                 self[{std::get<mio::Index<Tags>>(idx)...}] = value;
             })
        .def(
            "__iter__",
            [](const C& s) {
                return pybind11::make_iterator(s.begin(), s.end());
            },
            pybind11::keep_alive<0, 1>())
        .def("get_flat_index", &C::get_flat_index);

    // Not supported in Python yet: Slicing

    // bind all templated members for types in Tags...
    bind_templated_members_CustomIndexArray<C, Tags...>(c);
}

template <class C, class Base>
void bind_templated_members_Population(pybind11::class_<C, Base>&)
{
}

template <class C, class Base, class T, class... Ts>
void bind_templated_members_Population(pybind11::class_<C, Base>& c)
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
void bind_Population(pybind11::module& m, std::string const& name)
{
    using C    = mio::Populations<Cats...>;
    using Base = mio::CustomIndexArray<mio::UncertainValue, Cats...>;
    pybind11::class_<C, Base> c(m, name.c_str());
    c.def(pybind11::init([](mio::Index<Cats...> const& sizes, double val) {
         return C(sizes, val);
     }))
        .def(pybind11::init([](mio::Index<Cats...> const& sizes) {
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
auto bind_ParameterSet(pybind11::module& m, std::string const& name)
{
    pybind11::class_<ParameterSet> c(m, name.c_str());
    mio::foreach_tag<ParameterSet>([&c](auto t) {
        using Tag = decltype(t);

        //CAUTION: This requires ParameterTag::name() to be unique within the ParameterSet
        c.def_property(
            Tag::name().c_str(), [](const ParameterSet& self) -> auto& { return self.template get<Tag>(); },
            [](ParameterSet& self, typename Tag::Type const& v) {
                self.template get<Tag>() = v;
            },
            pybind11::return_value_policy::reference_internal);
    });
    return c;
}

/*
 * @brief bind a compartmental model for any Populations and Parameters class
 */
template <class Populations, class Parameters>
void bind_CompartmentalModel(pybind11::module& m, std::string const& name)
{
    using Model = mio::CompartmentalModel<Populations, Parameters>;
    pybind11::class_<Model>(m, name.c_str())
            .def(pybind11::init<Populations const&, Parameters const&>())
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
void bind_Simulation(pybind11::module& m, std::string const& name)
{
    pybind11::class_<Simulation>(m, name.c_str())
        .def(pybind11::init<const typename Simulation::Model&, double, double>(), pybind11::arg("model"), pybind11::arg("t0") = 0,
             pybind11::arg("dt") = 0.1)
        .def_property_readonly("result", pybind11::overload_cast<>(&Simulation::get_result, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def_property_readonly("model", pybind11::overload_cast<>(&Simulation::get_model, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def("advance", &Simulation::advance, pybind11::arg("tmax"));
}

template <typename Model>
void bind_SecirModelNode(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::Node<Model>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const mio::Node<Model>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property", [](const mio::Node<Model>& self) -> auto& { return self.property; },
            pybind11::return_value_policy::reference_internal);
}

template <typename Simulation>
void bind_SecirSimulationNode(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::Node<mio::SimulationNode<Simulation>>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const mio::Node<Simulation>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property",
            [](const mio::Node<mio::SimulationNode<Simulation>>& self) -> auto& { return self.property.get_simulation(); },
            pybind11::return_value_policy::reference_internal);
}

/*
 * @brief bind Graph for any node and edge type
 */
template <class Model>
void bind_SecirModelGraph(pybind11::module& m, std::string const& name)
{
    using G = mio::Graph<Model, mio::MigrationParameters>;
    pybind11::class_<G>(m, name.c_str())
        .def(pybind11::init<>())
        .def("add_node", &G::template add_node<const Model&>, py::arg("id"), py::arg("model"), pybind11::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const mio::MigrationParameters&>, py::arg("start_node_idx"), py::arg("end_node_idx"), py::arg("migration_parameters"),
             pybind11::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const Eigen::VectorXd&>, pybind11::return_value_policy::reference_internal)
        .def_property_readonly("num_nodes",
                               [](const G& self) {
                                   return self.nodes().size();
                               })
        .def(
            "get_node", [](const G& self, size_t node_idx) -> auto& { return self.nodes()[node_idx]; },
            pybind11::return_value_policy::reference_internal)
        .def_property_readonly("num_edges",
                               [](const G& self) {
                                   return self.edges().size();
                               })
        .def(
            "get_edge", [](const G& self, size_t edge_idx) -> auto& { return self.edges()[edge_idx]; },
            pybind11::return_value_policy::reference_internal)
        .def("get_num_out_edges",
             [](const G& self, size_t node_idx) {
                 return self.out_edges(node_idx).size();
             })
        .def(
            "get_out_edge",
            [](const G& self, size_t node_idx, size_t edge_idx) -> auto& { return self.out_edges(node_idx)[edge_idx]; },
            pybind11::return_value_policy::reference_internal);
}

template <class Simulation>
void bind_MigrationGraph(pybind11::module& m, std::string const& name)
{
    using G = mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>;
    pybind11::class_<G>(m, name.c_str())
        .def(pybind11::init<>())
        .def(
            "add_node", [](G & self, int id, const typename Simulation::Model& p, double t0, double dt) -> auto& {
                return self.add_node(id, p, t0, dt);
            },
            pybind11::arg("id"), pybind11::arg("model"), pybind11::arg("t0") = 0.0, pybind11::arg("dt") = 0.1,
            pybind11::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const mio::MigrationParameters&>, pybind11::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const Eigen::VectorXd&>, pybind11::return_value_policy::reference_internal)
        .def_property_readonly("num_nodes",
                               [](const G& self) {
                                   return self.nodes().size();
                               })
        .def(
            "get_node", [](const G& self, size_t node_idx) -> auto& { return self.nodes()[node_idx]; },
            pybind11::return_value_policy::reference_internal)
        .def_property_readonly("num_edges",
                               [](const G& self) {
                                   return self.edges().size();
                               })
        .def(
            "get_edge", [](const G& self, size_t edge_idx) -> auto& { return self.edges()[edge_idx]; },
            pybind11::return_value_policy::reference_internal)
        .def("get_num_out_edges",
             [](const G& self, size_t node_idx) {
                 return self.out_edges(node_idx).size();
             })
        .def(
            "get_out_edge",
            [](const G& self, size_t node_idx, size_t edge_idx) -> auto& { return self.out_edges(node_idx)[edge_idx]; },
            pybind11::return_value_policy::reference_internal);
}

/*
 * @brief bind GraphSimulation for any node and edge type
 */
template <class Graph>
void bind_GraphSimulation(pybind11::module& m, std::string const& name)
{
    using GS = mio::GraphSimulation<Graph>;
    pybind11::class_<GS>(m, name.c_str())
        .def(pybind11::init([](Graph& graph, double t0, double dt) {
                 return std::make_unique<GS>(mio::make_migration_sim(t0, dt, std::move(graph)));
             }),
             pybind11::arg("graph"), pybind11::arg("t0") = 0.0, pybind11::arg("dt") = 1.0)
        .def_property_readonly(
            "graph", [](GS& self) -> Graph& { return self.get_graph(); },
            pybind11::return_value_policy::reference_internal)
        .def_property_readonly("t", &GS::get_t)
        .def("advance", &GS::advance, pybind11::arg("tmax"));
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
    cl.def(pybind11::init<Eigen::Index, std::tuple_element_t<1, ArgTuples>...>(),
           pybind11::arg(std::get<0>(arg_tuples))..., pybind11::arg("size"));
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
        auto tup = pybind11::tuple(2);
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
        .def(pybind11::init([](const Eigen::Ref<const Matrix>& c, double t, int level, int type) {
                 return Damping(c, mio::DampingLevel(level), mio::DampingType(type), mio::SimulationTime(t));
             }),
             pybind11::arg("coeffs"), pybind11::arg("t"), pybind11::arg("level") = 0, pybind11::arg("type") = 0)
        .def_property(
            "coeffs", [](const Damping& self) -> const auto& { return self.get_coeffs(); },
            [](Damping& self, const Eigen::Ref<const Matrix>& v) {
                self.get_coeffs() = v;
            },
            pybind11::return_value_policy::reference_internal)
        .def_property(
            "time",
            [](const Damping& self) {
                return self.get_time();
            },
            [](Damping& self, double v) {
                self.get_time() = mio::SimulationTime(v);
            },
            pybind11::return_value_policy::reference_internal)
        .def_property(
            "type",
            [](const Damping& self) {
                return self.get_type();
            },
            [](Damping& self, int v) {
                self.get_type() = mio::DampingType(v);
            },
            pybind11::return_value_policy::reference_internal)
        .def_property(
            "level",
            [](const Damping& self) {
                return self.get_level();
            },
            [](Damping& self, int v) {
                self.get_level() = mio::DampingLevel(v);
            },
            pybind11::return_value_policy::reference_internal);
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
        .def(pybind11::init<const Eigen::Ref<const Matrix>&, const Eigen::Ref<const Matrix>&>(), pybind11::arg("baseline"),
             pybind11::arg("minimum"))
        .def(pybind11::init<const Eigen::Ref<const Matrix>&>(), pybind11::arg("baseline"));
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
            pybind11::return_value_policy::reference_internal)
        .def_property(
            "minimum", [](const DampingExpression& self) -> auto& { return self.get_minimum(); },
            [](DampingExpression& self, const Eigen::Ref<const Matrix>& v) {
                self.get_minimum() = v;
            },
            pybind11::return_value_policy::reference_internal)
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
                    throw pybind11::index_error("index out of range");
                }
                return self[i];
            },
            pybind11::return_value_policy::reference_internal)
        .def("__setitem__",
             [](DampingExpressionGroup& self, size_t i, const DampingExpression& m) {
                 if (i < 0 && i >= self.get_num_matrices()) {
                     throw pybind11::index_error("index out of range");
                 }
                 self[i] = m;
             })
        .def("get_matrix_at", [](const DampingExpressionGroup& self, double t) {
            return self.get_matrix_at(t);
        });
}

template <class Range>
auto bind_Range(py::module& m, const std::string& class_name)
{
    //bindings for iterator for the range
    struct Iterator {
        typename Range::Iterators iter_pair;
    };
    py::class_<Iterator>(m, (std::string("_Iter") + class_name).c_str())
        .def(
            "__next__", [](Iterator& self) -> auto&& {
                if (self.iter_pair.first != self.iter_pair.second) {
                    auto&& ref = *self.iter_pair.first;
                    ++self.iter_pair.first;
                    return ref;
                }
                throw py::stop_iteration();
            },
            py::return_value_policy::reference_internal);

    //bindings for the range itself
    py::class_<Range>(m, class_name.c_str())
        .def(
            "__iter__",
            [](Range& self) {
                return Iterator{{self.begin(), self.end()}};
            },
            py::keep_alive<1, 0>{}) //keep alive the Range as long as there is an iterator
        .def(
            "__getitem__", [](Range& self, size_t idx) -> auto&& { return self[idx]; },
            py::return_value_policy::reference_internal)
        .def("__len__", &Range::size);
}

//bind an enum class that can be iterated over
//requires the class to have a member `Count`
//adds a static `values` method to the enum class that returns an iterable list of the values
template<class E, class... Args>
auto iterable_enum(pybind11::module& m, const std::string& name, Args&&... args)
{
    using T = std::underlying_type_t<E>;

    struct Values
    {
    };
    pybind11::class_<Values>(m, (name + "Values").c_str(), std::forward<Args>(args)...)
        .def("__iter__", [](Values& /*self*/) {
            return E(0);
        })
        .def("__len__", [](Values& /*self*/) {
            return E::Count;
        });

    auto enum_class = pybind11::enum_<E>(m, name.c_str(), std::forward<Args>(args)...);
    enum_class.def_static("values", [](){
        return Values{};
    });
    enum_class.def("__next__", [](E& self) {
        if (self < E::Count) {
            auto current = self;
            self = E(T(self) + T(1));
            return current;
        } else {
            throw pybind11::stop_iteration();
        }
    });
    return enum_class;
}

} // namespace pymio

#endif //PYMIO_TEMPLATES_H