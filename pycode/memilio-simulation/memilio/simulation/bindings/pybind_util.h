/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/io/json_serializer.h"
#include "memilio/io/io.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/math/eigen_util.h"

#include "pybind11/pybind11.h"

namespace pymio
{

// Enum for controlling pickling support behavior
enum class EnablePickling {
    Never,      // Pickling support is disabled.
    IfAvailable, // Pickling support is enabled if the class matches SFINEA-clause of `serialize_internal` 
    Required    // Pickling support is required, and an error is raised if not available.
};

// Tag dispatch to convert value into type
template <EnablePickling F>
struct PicklingTag {};

// Detect serialization by matching a `serialize_internal` function
template <class IOContext, class T>
using serialize_internal_t = decltype(mio::serialize_internal(std::declval<IOContext&>(), std::declval<T&>()));
template <class IOContext, class T>
using has_serialize_internal = mio::is_expression_valid<serialize_internal_t, IOContext, T>;

template <class T, class... Args>
void pybind_pickle_class(pybind11::class_<T, Args...>& cls)
{
    cls.def(pybind11::pickle(
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
}

/**
 * Call binding for a class with bind_class<...>(...).
 * Strategy for pickling is deduced depending on the EnablePickling value.
 * Compile-time errors originating from this function are likely due to issues with the serialization of the given class in the C++ library.
 * 
 * Here's a small guideline on when to use each strategy:
 * 
 * 1. Use EnablePickling::Never for classes:
 *    - Containing sensitive or non-serializable data.
 *    - With complex ownership semantics, especially when serializing the object might lead to unintended side effects.
 *    - Where pickling/unpickling is not feasible or doesn't make sense.
 * 
 * 2. Use EnablePickling::IfAvailable for classes:
 *    - That are generic and may become picklable in the future or in specific use cases.
 *    - When you want to allow pickling if the necessary serialization functions are provided in the future.
 * 
 *    Important Note:
 *    This strategy attempts to deduce if serialization is possible during compile-time by checking for a matching serialize_internal function.
 *    The C++ library does not check pickling for every class during compile-time, which can result in classes unexpectedly lacking serialization support.
 *    This approach does not account for subclasses contained within the class you intend to pickle. As a result, classes can be deduced as serializable even if they are not, 
 *    leading to a compile-time error. In such cases, it may be either an error in the C++ class implementation or the class should not be serializable and should use the EnablePickling::Never strategy during binding.
 * 
 * 3. Use EnablePickling::Required for classes:
 *    - That are intended to be used with pickling, and it's critical to have pickling support.
 *    - When you want to ensure serialization in the C++ library.
 * 
 * @{
 */
template<class T, EnablePickling F, class... Args>
struct BindClassHelper
{
private:
    template <class... Options>
    auto _bind_class(pybind11::module& m, std::string const& name, PicklingTag<EnablePickling::Never> /*tags*/, Options&&... options) const {
        auto cls = pybind11::class_<T, Args...>(m, name.c_str(), std::forward<Options>(options)...);
        return cls;
    }

    template <class... Options>
    auto _bind_class(pybind11::module& m, std::string const& name, PicklingTag<EnablePickling::IfAvailable> /*tags*/, Options&&... options) const {
        auto cls = pybind11::class_<T, Args...>(m, name.c_str(), std::forward<Options>(options)...);
        // Bind the class depending on its features
        if constexpr (has_serialize_internal<mio::PickleSerializer, T>::value) {
            pybind_pickle_class<T, Args...>(cls);
        }
        return cls;
    }

    template <class... Options>
    auto _bind_class(pybind11::module& m, std::string const& name, PicklingTag<EnablePickling::Required> /*tags*/, Options&&... options) const {
        auto cls = pybind11::class_<T, Args...>(m, name.c_str(), std::forward<Options>(options)...);
        pybind_pickle_class<T, Args...>(cls);
        return cls;
    }

public:
    template <class... Options>
    auto operator()(pybind11::module& m, std::string const& name, Options&&... options) const {
        return _bind_class(m, name, PicklingTag<F>{}, std::forward<Options>(options)...);
    }
};

/**
 * bind_class is defined as a function-like object using the overloaded call operator of BindClassHelper.
 * It is a nibloid, that means:
 *      - Explicit template argument lists may not be specified when calling it.
 *      - It is visible to argument-dependent lookup, because it is a function object.
 *      - When it is found by normal unqualified lookup for the name to the left of the function-call operator, it inhibits argument-dependent lookup.
 * In this case, nibloid is used to have a function-like structure with two template parameter packs, by hiding one inside the operator of BindClassHelper.
 * The first template parameter pack 'Args' needs to be defined explicitly, while the second template parameter pack 'Options' is deduced automatically.
 * @tparam T class for binding
 * @tparam F value of enum EnablePickling defining pickling behaviour
 * @tparam Args base class of T.
 * @param m the pybind11 module.
 * @param name of the class in python.
 * @param options optional arguments for pybind11::class_.
 * @return instance of pybind class_.
 */
template<class T, EnablePickling F, class... Args>
constexpr BindClassHelper<T, F, Args...> bind_class;
/**@}*/

template <typename T>
T check_and_throw(mio::IOResult<T>& result)
{
    if (result.has_error()) {
        auto status = result.error();
        if (status.code() == std::errc::no_such_file_or_directory) {
            throw pybind11::value_error(status.message());
        } else {
            throw std::runtime_error(status.message());
        }
    } else {
        return result.value();
    }
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
          class = std::enable_if_t<(std::is_same<typename C::type::Shape, mio::SquareMatrixShape<double>>::value ||
                                    std::is_same<typename C::type::Shape, mio::ColumnVectorShape<double>>::value),
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
* Bind a specialization of mio::Range class template.
* The python class will be a read-only container/iterable.
* You probably also want to use PYMIO_IGNORE_VALUE_TYPE to
* enable copying and moving in all cases.
*/
template <class Range>
auto bind_Range(pybind11::module_& m, const std::string& class_name)
{
    //bindings for iterator for the range
    struct Iterator {
        typename Range::Iterators iter_pair;
    };
    bind_class<Iterator, EnablePickling::Never>(m, (std::string("_Iter") + class_name).c_str())
        .def(
            "__next__",
            [](Iterator& self) -> auto&& {
                if (self.iter_pair.first != self.iter_pair.second) {
                    auto&& ref = *self.iter_pair.first;
                    ++self.iter_pair.first;
                    return ref;
                }
                throw pybind11::stop_iteration();
            },
            pybind11::return_value_policy::reference_internal);

    //bindings for the range itself
    bind_class<Range, EnablePickling::Never>(m, class_name.c_str())
        .def(
            "__iter__",
            [](Range& self) {
                return Iterator{{self.begin(), self.end()}};
            },
            pybind11::keep_alive<1, 0>{}) //keep alive the Range as long as there is an iterator
        .def(
            "__getitem__",
            [](Range& self, size_t idx) -> auto&& {
                return self[idx];
            },
            pybind11::return_value_policy::reference_internal)
        .def("__len__", &Range::size);
}

/**
* A Range looks like a container, so pybind11 checks the value_type alias to see if 
* the Range can be copied or moved. But since a Range is just a view and does not own its contents,
* it can always be copied and moved, even if the value_type is uncopieable/unmovable.
* This macro stops the value_type check.
* see {Pybind11_SRC_DIR}/include/pybind11/detail/type_caster_base.h and {Pybind11_SRC_DIR}/tests/test_stl_binders.cpp
*/
#define PYMIO_IGNORE_VALUE_TYPE(Range)                                                                                 \
    namespace pybind11                                                                                                 \
    {                                                                                                                  \
    namespace detail                                                                                                   \
    {                                                                                                                  \
    template <typename SFINAE>                                                                                         \
    struct recursive_container_traits<Range, SFINAE> {                                                                 \
        using type_to_check_recursively = recursive_bottom;                                                            \
    };                                                                                                                 \
    }                                                                                                                  \
}

//bind an enum class that can be iterated over
//requires the class to have a member `Count`
//adds a static `values` method to the enum class that returns an iterable list of the values
template <class E, class... Args>
auto iterable_enum(pybind11::module_& m, const std::string& name, Args&&... args)
{
    using T = std::underlying_type_t<E>;
    auto enum_class = pybind11::enum_<E>(m, name.c_str(), std::forward<Args>(args)...);
    
    //dummy type that provides iteration of enums
    //not meant to be used directly by users, so name starts with _
    struct Values {
    };
    bind_class<Values, EnablePickling::Never>(m, ("_" + name + "Values").c_str(), std::forward<Args>(args)...)
        .def("__iter__",
             [](Values& /*self*/) {
                 return E(0);
             })
        .def("__len__", [](Values& /*self*/) {
            return (size_t)E::Count; //len() expects integer
        });

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
