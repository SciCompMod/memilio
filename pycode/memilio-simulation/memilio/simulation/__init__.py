#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Daniel Abele
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################

"""
The Python bindings to the MEmilio C++ library.
"""

from memilio.simulation._simulation import *


def __getattr__(attr):
    """ The __getattr__ function is used here to implement lazy loading for the submodules within the memilio.simulation package. 
    Submodules are only imported when they are first accessed, which can save memory and reduce startup time, if not all submodules 
    are needed for every execution.
    """
    if attr == "abm":
        import memilio.simulation.abm as abm
        return abm
    elif attr == "osir":
        import memilio.simulation.osir as osir
        return osir
    elif attr == "oseir":
        import memilio.simulation.oseir as oseir
        return oseir
    elif attr == "osecir":
        import memilio.simulation.osecir as osecir
        return osecir
    elif attr == "osecirvvs":
        import memilio.simulation.osecirvvs as osecirvvs
        return osecirvvs
    elif attr == "ssir":
        import memilio.simulation.ssir as ssir
        return ssir

    elif attr == "ssirs":
        import memilio.simulation.ssirs as ssirs
        return ssirs

    raise AttributeError("module {!r} has no attribute "
                         "{!r}".format(__name__, attr))
