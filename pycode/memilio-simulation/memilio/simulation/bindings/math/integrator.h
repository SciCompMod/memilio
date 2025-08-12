/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Maximilian Betz
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
#ifndef PYMIO_INTEGRATOR_H
#define PYMIO_INTEGRATOR_H

#include "memilio/math/integrator.h"
#include "memilio/math/euler.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/math/stepper_wrapper.h"
#include <Eigen/Dense>
#include "pybind_util.h"

#include "pybind11/pybind11.h"
#include <pybind11/eigen.h>

namespace pymio
{

void bind_Integrator_Core(pybind11::module_& m)
{
    pymio::bind_class<mio::OdeIntegratorCore<double>, pymio::EnablePickling::Never, 
                      pybind11::smart_holder>(m, "IntegratorCore")
        .def_property("dt_max",
                      pybind11::overload_cast<>(&mio::OdeIntegratorCore<double>::get_dt_max, pybind11::const_),
                      [](mio::OdeIntegratorCore<double>& self, double dt_max) {
                          self.get_dt_max() = dt_max;
                      })
        .def_property("dt_min",
                      pybind11::overload_cast<>(&mio::OdeIntegratorCore<double>::get_dt_min, pybind11::const_),
                      [](mio::OdeIntegratorCore<double>& self, double dt_min) {
                          self.get_dt_min() = dt_min;
                      });

    pymio::bind_class<mio::EulerIntegratorCore<double>, pymio::EnablePickling::Never, mio::OdeIntegratorCore<double>,
                      pybind11::smart_holder>(m, "EulerIntegratorCore")
        .def(pybind11::init<>())
        .def(
            "step",
            [](const mio::EulerIntegratorCore<double>& self, pybind11::function f, Eigen::Ref<const Eigen::VectorXd> yt,
               double t, double dt, Eigen::Ref<Eigen::VectorXd> ytp1) {
                bool result = self.step(
                    [f](Eigen::Ref<const Eigen::VectorXd> y, double t, Eigen::Ref<Eigen::VectorXd> dydt) {
                        f(y, t, dydt);
                    },
                    yt, t, dt, ytp1);
                return result;
            },
            pybind11::arg("f"), pybind11::arg("yt"), pybind11::arg("t"), pybind11::arg("dt"), pybind11::arg("ytp1"));

    using RungeKuttaCashKarp54Integrator =
        mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>;
    pymio::bind_class<RungeKuttaCashKarp54Integrator, pymio::EnablePickling::Never, mio::OdeIntegratorCore<double>,
                      pybind11::smart_holder>(m, "RungeKuttaCashKarp54IntegratorCore")
        .def(pybind11::init<>())
        .def(pybind11::init<const double, const double, const double, const double>(), pybind11::arg("abs_tol"),
             pybind11::arg("rel_tol"), pybind11::arg("dt_min"), pybind11::arg("dt_max"))
        .def("set_abs_tolerance", &RungeKuttaCashKarp54Integrator::set_abs_tolerance, pybind11::arg("tol"))
        .def("set_rel_tolerance", &RungeKuttaCashKarp54Integrator::set_rel_tolerance, pybind11::arg("tol"));

    pymio::bind_class<mio::RKIntegratorCore<double>, pymio::EnablePickling::Never, mio::OdeIntegratorCore<double>,
                      pybind11::smart_holder>(m, "RKIntegratorCore")
        .def(pybind11::init<>())
        .def(pybind11::init<double, double, double, double>(), pybind11::arg("abs_tol") = 1e-10,
             pybind11::arg("rel_tol") = 1e-5, pybind11::arg("dt_min") = std::numeric_limits<double>::min(),
             pybind11::arg("dt_max") = std::numeric_limits<double>::max())
        .def("set_abs_tolerance", &mio::RKIntegratorCore<double>::set_abs_tolerance, pybind11::arg("tol"))
        .def("set_rel_tolerance", &mio::RKIntegratorCore<double>::set_rel_tolerance, pybind11::arg("tol"));
}
} // namespace pymio

#endif //PYMIO_INTEGRATOR_H
