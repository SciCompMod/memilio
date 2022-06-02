#ifndef PYMIO_DATE_H
#define PYMIO_DATE_H

#include "pybind_util.h"
#include "memilio/utils/date.h"

#include "pybind11/pybind11.h"

namespace pymio
{

void bind_date(pybind11::module& m, std::string const& name);

} // namespace pymio

#endif //PYMIO_DATE_H