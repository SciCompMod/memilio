#ifndef PYMIO_POPULATIONS_H
#define PYMIO_POPULATIONS_H

#include "pybind_util.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/epidemiology/populations.h"

#include "pybind11/pybind11.h"

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

} // namespace pymio

#endif //PYMIO_POPULATIONS_H