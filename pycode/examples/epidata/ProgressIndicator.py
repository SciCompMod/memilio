#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Rene Schmieding
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

from memilio.epidata import progress_indicator
import time

print("This is only a usage example, and does not actually do anything.")
# Also, the following values for delay, sleep etc. are chosen arbitrary,
# and have no further relevancy other than to demonstrate the indicator.

# using start/stop
p = progress_indicator.Dots(message="waiting", delay=0.5)
p.start()
time.sleep(1.6)
p.stop()

# using with-as block
with progress_indicator.Percentage(message="download 1", delay=0.4) as p:
    for i in range(13):
        time.sleep(0.1467)
        p.set_progress((i + 1) / 13)

with progress_indicator.Percentage(message="download 2", use_bar=False,
                                   delay=0, keep_output=False) as p:
    for i in range(97):
        time.sleep(0.0367)
        p.set_progress((i + 1) / 97)

# using with block
# the 'as' is only required for calling e.g. message() or set_progress()
with progress_indicator.Spinner(message="finish"):
    time.sleep(2)
