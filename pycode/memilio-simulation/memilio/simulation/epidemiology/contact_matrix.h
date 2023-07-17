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
#ifndef PYMIO_CONTACT_MATRIX_H
#define PYMIO_CONTACT_MATRIX_H

#include "memilio/epidemiology/contact_matrix.h"
#include "pybind_util.h"

#include "pybind11/pybind11.h"

namespace pymio
{

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
        .def(pybind11::init<const Eigen::Ref<const Matrix>&, const Eigen::Ref<const Matrix>&>(),
             pybind11::arg("baseline"), pybind11::arg("minimum"))
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

} // namespace pymio

#endif //PYMIO_CONTACT_MATRIX_H
