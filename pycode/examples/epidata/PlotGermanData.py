#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors:
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
"""@PlotGermanData.py
WARNING: This file is currently not tested and maintained.
"""
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_json('cases_all_age.json')

for i in range(len(df)):
    df.Date[i] = df.Date[i].date()


ages = df.Age_RKI.unique()
num_groups = len(ages)
time = range(126)
dates = pd.date_range(start=datetime(2020, 1, 28), end=datetime(2020, 6, 2))
group_data = np.zeros((len(time), num_groups, 3))

for age, i in zip(ages, range(num_groups)):
    for date, j in zip(dates, range(len(time))):
        if date.date() in [x.date() for x in df.Date[df.Age_RKI == age]]:
            group_data[j, i, 0] = df.Confirmed[(
                df.Age_RKI == age) & (df.Date == date)].values
            group_data[j, i, 1] = df.Deaths[(
                df.Age_RKI == age) & (df.Date == date)].values
            group_data[j, i, 2] = df.Recovered[(
                df.Age_RKI == age) & (df.Date == date)].values

datelist = np.array(pd.date_range(datetime(2020, 1, 28),
                    periods=len(time), freq='D').strftime('%m-%d').tolist())
tick_range = np.arange(int(len(time)/10)+1)*10
compartiments = ['Confirmed', 'Dead', 'Recovered']
fig, ax = plt.subplots(3, 1, figsize=(8, 10))
for comp, j in zip(compartiments, range(len(compartiments))):
    for age, i in zip(ages, range(len(ages))):
        ax[j].plot(time, group_data[:, i, j], label=age)
        ax[j].legend()
        ax[j].set_title(comp)
        ax[j].set_xticks([])
ax[2].set_xticks(tick_range)
ax[2].set_xticklabels(datelist[tick_range], rotation=45)
fig.tight_layout()
fig.savefig('Cases_Groups.pdf')


all_data = np.sum(group_data, axis=1)
new_inf = []
for i in range(1, len(all_data[:, 0])):
    new_inf.append(all_data[i, 0] - all_data[i-1, 0])

fig, ax = plt.subplots(4, 1, figsize=(8, 14))
for comp, j in zip(compartiments, range(len(compartiments))):
    ax[j].plot(time, np.sum(group_data[:, :, j], axis=1))
    ax[j].set_title(comp)
    ax[j].set_xticks([])
ax[3].plot(time[:-1], new_inf)
ax[3].set_title('New Daily Infections')
ax[3].set_xticks(tick_range)
ax[3].set_xticklabels(datelist[tick_range], rotation=45)
fig.tight_layout()
fig.savefig('Cases_All.pdf')
plt.show()
