#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Anna Wendler
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


def compute_optimal_num_subcompartments(mean, std):

    optimal_num_subcompartments = mean**2/std**2

    return optimal_num_subcompartments


def get_weighted_value(prob_1, value_1, value_2):

    weighted_value = prob_1*value_1 + (1-prob_1)*value_2

    return weighted_value


def main():

    # Mean and standard deviation of lognormal distributions per transiton from Covasim paper.
    mean_and_std_per_transition = {"ExposedToInfectedNoSymptoms": [4.5, 1.5],
                                   "InfectedNoSymptomsToInfectedSymptoms": [1.1, 0.9],
                                   "InfectedNoSymptomsToRecovered": [8.0, 2.0],
                                   "InfectedSymptomsToInfectedSevere": [6.6, 4.9],
                                   "InfectedSymptomsToRecovered": [8.0, 2.0],
                                   "InfectedSevereToInfectedCritical": [1.5, 2.0],
                                   "InfectedSevereToRecovered": [18.1, 6.3],
                                   "InfectedCriticalToDead": [10.7, 4.8],
                                   "InfectedCriticalToRecovered": [18.1, 6.3]}

    # Transition probabilities from Assessment paper, see LCT paper.
    # Store only the probabilities that are not considering the transition to the Recovered compartment as this is
    # determined by the given probabilities.
    transition_probabilities_non_recovered = {"ExposedToInfectedNoSymptoms": [1., 1., 1., 1., 1., 1.],
                                              "InfectedNoSymptomsToInfectedSymptoms": [0.75, 0.75, 0.8, 0.8, 0.8, 0.8],
                                              "InfectedSymptomsToInfectedSevere": [0.0075, 0.0075, 0.019, 0.0615, 0.0615, 0.225],
                                              "InfectedSevereToInfectedCritical": [0.075, 0.075, 0.075, 0.15, 0.3, 0.4],
                                              "InfectedCriticalToDead": [0.05, 0.05, 0.14, 0.14, 0.4, 0.6]}

    compartments = ["Exposed", "InfectedNoSymptoms",
                    "InfectedSymptoms", "InfectedSevere", "InfectedCritical"]

    # Get the number of subcompartment for the LCT model.
    # To get the numvber of subcompartments, we first calculate the optimal number of subcompartments that is not
    # necessarily an integer. In the LCT model we cannot consider different distributions for e.g. the transition
    # for InfectedNoSymptomsToInfectedSymptoms and InfectedNoSymptomsToRecovered but only have one transition for
    # leaving the InfectedNoSymptoms compartments. This is why we weigh the optimal number of subcompartments
    # of both transitions with the corresponding transition probability.
    #
    # For the same reason, we weigh the mean of the corresponding transitions from Covasim by the corresponding
    # transition probabilities.
    # Note that the mean is the same for the LCT and the ODE model, only the number of subcompartments differs and is
    # always one for the ODE model.

    print("Number of subcompartments for the LCT and mean for the LCT and ODE models:")
    print()

    # Compute the optimal number of subcompartments for each transition.
    optimal_num_subcompartments = []
    for i in mean_and_std_per_transition.keys():
        optimal_num_subcompartment = compute_optimal_num_subcompartments(
            mean_and_std_per_transition[i][0], mean_and_std_per_transition[i][1])
        optimal_num_subcompartments.append(optimal_num_subcompartment)

    for i, transition in enumerate(transition_probabilities_non_recovered.keys()):

        # For the transition ExposedToInfectedNoSymptoms we do not need to weigh anything as there is only one
        # possible transition from Exposed.
        if transition == "ExposedToInfectedNoSymptoms":
            weighted_num_subcompartments_per_age = []
            weighted_mean_per_age = []

            for age in range(len(transition_probabilities_non_recovered[transition])):
                weighted_num_subcompartments_per_age.append(round(
                    optimal_num_subcompartments[i]))
                weighted_mean_per_age.append(
                    mean_and_std_per_transition[transition][0])

        else:
            weight = transition_probabilities_non_recovered[transition]
            # Get the indices of the two transitions that we need to weigh.
            first_index = 2*i-1
            second_index = 2*i
            # print(
            #     f"Weighing transitions {keys_list[first_index]} and {keys_list[second_index]} with weight {weight}.")

            # Get the corresponding keys of the transitions.
            keys_list = list(mean_and_std_per_transition.keys())
            first_key = keys_list[first_index]
            second_key = keys_list[second_index]

            # Weigh the optimal number of subcompartments of both transitions with the corresponding transition probability.
            weighted_num_subcompartments_per_age = []
            weighted_mean_per_age = []
            for age in range(len(transition_probabilities_non_recovered[transition])):
                # Round to integer.
                weighted_num_subcompartments_per_age.append(round(get_weighted_value(
                    weight[age], optimal_num_subcompartments[first_index], optimal_num_subcompartments[second_index])))
                # Round to 8 digits.
                weighted_mean_per_age.append(round(get_weighted_value(
                    weight[age], mean_and_std_per_transition[first_key][0], mean_and_std_per_transition[second_key][0]), 8))

        print(compartments[i])
        print(
            f"Num subcompartments: {weighted_num_subcompartments_per_age}")
        print(f"Mean: {weighted_mean_per_age}")
        print()


if __name__ == '__main__':
    main()
