#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Daniel Richter
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
Default mapping dictionary used for generating or configuring model-related code components.

This dictionary provides default values for various model-related terms and file names, which can be used throughout the code generation process. 
It serves as a centralized reference for consistent naming conventions and can be easily modified if needed.
"""

default_dict = {
    "model": "Model",
    "agegroup": "AgeGroup",
    "emptystring": "",
    "simulation": "Simulation",
    "flowmodel": "FlowModel",
    "compartmentalmodel": "CompartmentaModel",
    "modelfile": "model.h",
    "infectionstatefile": "infection_state.h",
    "parameterspacefile": "parameter_space.h",
    "analyzeresultfile": "analyze_result.h",
    "namespace": "namespace",
    "mio": "mio",
    "drawsample": "draw_sample"
}

""" 
Configuration list defining which bindings should be generated.

Each entry specifies a class or function to be created, including its name,
type, and associated metadata.
For class entries, the "methods" field lists the member functions that
should be generated as part of the binding.
"""

general_bindings_dict = [
    {
        "type": "class",
        "name": "Parameters",
        "cursorkind": "CLASS_TEMPLATE",
        "methods": [
            "check_constraints",
            "apply_constraints",
            "get_start_commuter_detection",
            "get_end_commuter_detection",
            "get_commuter_nondetection"
        ]
    },

    {
        "type": "class",
        "name": "ModelBase",
        "cursorkind": "CLASS_TEMPLATE",
        "methods": []
    },
    {
        "type": "function",
        "name": "draw_sample",
        "cursorkind": "FUNCTION_TEMPLATE",
    },
    {
        "type": "function",
        "name": "simulate_flows",
        "cursorkind": "FUNCTION_TEMPLATE"
    },
    {
        "type": "function",
        "name": "simulate",
        "cursorkind": "FUNCTION_TEMPLATE"
    }


]
