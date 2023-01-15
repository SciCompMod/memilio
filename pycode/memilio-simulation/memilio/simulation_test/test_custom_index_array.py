#############################################################################
# Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
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

from unittest import TestCase, main

from memilio.simulation import UncertainValue
from memilio.simulation.secir import (AgeGroup, AgeGroupArray, InfectionState,
                                      SecirPopulationArray)


class TestCustomIndexArray(TestCase):
    def test_init(self):
        dims = (AgeGroup(5), InfectionState(len(InfectionState.values())))
        array = SecirPopulationArray(dims)
        self.assertEqual(array.size_AgeGroup(), dims[0])
        self.assertEqual(array.size_InfectionState(), dims[1])
        self.assertEqual(array.size(), dims)

    def test_init_single_index(self):
        dims = AgeGroup(3)
        array = AgeGroupArray(dims)
        self.assertEqual(array.size_AgeGroup(), dims)
        self.assertEqual(array.size(), dims)

    def test_assign(self):
        dims = (AgeGroup(5), InfectionState(len(InfectionState.values())))
        array = SecirPopulationArray(dims)
        array[:, :] = 1.0
        array[:, InfectionState.Exposed] = 2.0
        array[AgeGroup(0), InfectionState.InfectedNoSymptoms] = 3.0
        self.assertEqual(
            array[AgeGroup(0), InfectionState.Susceptible].value, 1.0)
        self.assertEqual(array[AgeGroup(1), InfectionState.Exposed].value, 2.0)
        self.assertEqual(
            array[AgeGroup(0), InfectionState.InfectedNoSymptoms].value, 3.0)

    def test_assign_steps(self):
        dims = (AgeGroup(5), InfectionState(len(InfectionState.values())))
        array = SecirPopulationArray(dims)
        array[::AgeGroup(1), InfectionState.InfectedNoSymptoms] = 1.0
        array[::AgeGroup(2), InfectionState.Exposed] = 2.0
        self.assertEqual(
            array[AgeGroup(0), InfectionState.InfectedNoSymptoms].value, 1.0)
        self.assertEqual(
            array[AgeGroup(1), InfectionState.InfectedNoSymptoms].value, 1.0)
        self.assertEqual(array[AgeGroup(0), InfectionState.Exposed].value, 2.0)
        self.assertEqual(array[AgeGroup(1), InfectionState.Exposed].value, 0.0)

    def test_iter(self):
        dims = AgeGroup(3)
        array = AgeGroupArray(dims)
        for v in array:
            v.value = 1
        self.assertEqual(array[AgeGroup(0)].value, 1)
        self.assertEqual(array[AgeGroup(2)].value, 1)


if __name__ == '__main__':
    main()
