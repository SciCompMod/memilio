#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Maximilian Betz
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
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from memilio.simulation.abm import forward_pass 
import time

start = time.time()

for i in range(100):
    test = forward_pass(0.08, 5)
end = time.time()
print(end - start)
#df = pd.DataFrame(test)
#df.to_csv("test.csv")

#for i in range(1,11):
#    plt.plot(df[0],df[i])
#plt.show()

#print(test)