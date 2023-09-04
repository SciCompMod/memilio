/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin Siggel, Daniel Abele, Martin J. Kuehn, Jan Kleinert, Maximilian Betz
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
#ifndef PYMIO_PYBIND_UTIL_H
#define PYMIO_PYBIND_UTIL_H

#include "memilio/math/matrix_shape.h"
#include "pickle_serializer.h"

#include "pybind11/pybind11.h"

namespace pymio
{

//bind class and add pickling based on memilio serialization framework
template <class T, class... Args>
decltype(auto) pybind_pickle_class(pybind11::module_& m, const char* name)
{
    decltype(auto) pickle_class = pybind11::class_<T, Args...>(m, name);
    pickle_class.def(pybind11::pickle(
        [](const T& object) { // __getstate__
            auto tuple = mio::serialize_pickle(object);
            if (tuple) {
                return std::move(tuple).value();
            }
            else {
                throw std::runtime_error(tuple.error().formatted_message());
            }
        },
        [](const pybind11::tuple t) { // __setstate__
            auto object = mio::deserialize_pickle(t, mio::Tag<T>{});
            if (object) {
                return std::move(object).value();
            }
            else {
                throw std::runtime_error(object.error().formatted_message());
            }
        }));
    return pickle_class;
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

template <class Range>
auto bind_Range(pybind11::module_& m, const std::string& class_name)
{
    //bindings for iterator for the range
    struct Iterator {
        typename Range::Iterators iter_pair;
    };
    pybind11::class_<Iterator>(m, (std::string("_Iter") + class_name).c_str())
        .def(
            "__next__", [](Iterator & self) -> auto&& {
                if (self.iter_pair.first != self.iter_pair.second) {
                    auto&& ref = *self.iter_pair.first;
                    ++self.iter_pair.first;
                    return ref;
                }
                throw pybind11::stop_iteration();
            },
            pybind11::return_value_policy::reference_internal);

    //bindings for the range itself
    pybind11::class_<Range>(m, class_name.c_str())
        .def(
            "__iter__",
            [](Range& self) {
                return Iterator{{self.begin(), self.end()}};
            },
            pybind11::keep_alive<1, 0>{}) //keep alive the Range as long as there is an iterator
        .def(
            "__getitem__", [](Range & self, size_t idx) -> auto&& { return self[idx]; },
            pybind11::return_value_policy::reference_internal)
        .def("__len__", &Range::size);
}

//bind an enum class that can be iterated over
//requires the class to have a member `Count`
//adds a static `values` method to the enum class that returns an iterable list of the values
template <class E, class... Args>
auto iterable_enum(pybind11::module_& m, const std::string& name, Args&&... args)
{
    using T = std::underlying_type_t<E>;

    //dummy type that provides iteration of enums
    //not meant to be used directly by users, so name starts with _
    struct Values {
    };
    pybind11::class_<Values>(m, ("_" + name + "Values").c_str(), std::forward<Args>(args)...)
        .def("__iter__",
             [](Values& /*self*/) {
                 return E(0);
             })
        .def("__len__", [](Values& /*self*/) {
            return (size_t)E::Count; //len() expects integer
        });

    auto enum_class = pybind11::enum_<E>(m, name.c_str(), std::forward<Args>(args)...);
    enum_class.def_static("values", []() {
        return Values{};
    });
    enum_class.def("__next__", [](E& self) {
        if (self < E::Count) {
            auto current = self;
            self         = E(T(self) + T(1));
            return current;
        }
        else {
            throw pybind11::stop_iteration();
        }
    });
    return enum_class;
}

// If the python object is None: returns empty optional.
// Otherwise: casts the python object to type T.
// Throws exception if obj cannot be cast to T.
template <class T, class Obj>
boost::optional<T> cast_or_none(Obj&& obj)
{
    if (obj.is_none()) {
        return {};
    }
    else {
        return obj.template cast<T>();
    }
}

} // namespace pymio

#endif //PYMIO_PYBIND_UTIL_H
