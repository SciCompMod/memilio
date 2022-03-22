#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
MEmilio main namespace package.
"""

__path__ = __import__('pkgutil').extend_path(__path__, __name__)

import os
if os.name == 'nt': # can be nt, posix, or java
    # Windows uses nt, which does not support carriage returns by default
    # the following Windows specific module should fix this
    try:
        from ctypes import windll
        k = windll.kernel32
        # GetStdHandle(-11) is the standart output device, the flag 7 is equal to
        # ENABLE_PROCESSED_OUTPUT | ENABLE_WRAP_AT_EOL_OUTPUT | ENABLE_VIRTUAL_TERMINAL_PROCESSING
        # where the last flag enables several control sequences, like \r
        k.SetConsoleMode(k.GetStdHandle(-11), 7)
    except:
        print("Error: Failed to set console mode for 'nt' system (e.g. Windows).")
        print("       Some output may be displayed incorrectly.")