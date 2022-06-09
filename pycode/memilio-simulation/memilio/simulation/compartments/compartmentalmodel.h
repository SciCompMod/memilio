#ifndef PYMIO_COMPARTMENTALMODEL_H
#define PYMIO_COMPARTMENTALMODEL_H

#include "memilio/compartments/compartmentalmodel.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{

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

} // namespace pymio

#endif //PYMIO_COMPARTMENTALMODEL_H