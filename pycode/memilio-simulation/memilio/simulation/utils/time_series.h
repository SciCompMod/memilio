#ifndef PYMIO_TIME_SERIES_H
#define PYMIO_TIME_SERIES_H

#include "pybind11/pybind11.h"

namespace pymio
{

void bind_time_series(pybind11::module& m, std::string const& name);

} // namespace pymio

#endif //PYMIO_TIME_SERIES_H