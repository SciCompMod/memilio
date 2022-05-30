#ifndef PYMIO_DAMPING_H
#define PYMIO_DAMPING_H

#include "memilio/epidemiology/damping.h"
#include "pybind_util.h"

#include "pybind11/pybind11.h"

namespace pymio
{

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

} // namespace pymio

#endif //PYMIO_DAMPING_H