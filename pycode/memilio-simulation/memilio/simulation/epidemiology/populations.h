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
#ifndef PYMIO_POPULATIONS_H
#define PYMIO_POPULATIONS_H

#include "pybind_util.h"
#include "utils/custom_index_array.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/epidemiology/populations.h"

#include "pybind11/pybind11.h"

#include <stdexcept>

namespace pymio
{

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
void bind_Population(pybind11::module_& m, std::string const& name, mio::Tag<mio::Populations<Cats...>> /*tags*/)
{
    using C    = mio::Populations<Cats...>;
    using Base = mio::CustomIndexArray<mio::UncertainValue, Cats...>;

    // Catch warning ImportError: generic_type: type "" is already registered!
    try {
        bind_CustomIndexArray<mio::UncertainValue, Cats...>(m, (name + "Array").c_str());
    }
    catch (std::runtime_error& e) {
    }

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

} // namespace pymio

#endif //PYMIO_POPULATIONS_H
