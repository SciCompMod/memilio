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
"""
:strong:`intermediate_representation.py`
Dataclass to represent the needed information of a model. Interface between the scanner and generator.
"""
from dataclasses import dataclass, field
from typing import Any, Dict, Union

from typing_extensions import Self
from memilio.generation import Generator


@dataclass
class IntermediateRepresentation:
    """
    Dataclass storing the model features. Serves as interface between Scanner and Generator.
    """
    namespace: str = ""
    model_class: str = ""
    python_module_name: str = ""
    parameterset: str = ""
    parameterset_wrapper: str = ""
    simulation: bool = False
    is_compartmentalmodel: bool = False
    is_flowmodel: bool = False
    has_age_group: bool = False
    has_draw_sample: bool = False
    scalartype: str = "double"
    python_generation_module_path: str = ""
    target_folder: str = ""
    enum_populations: dict = field(default_factory=dict)
    model_init: list = field(default_factory=list)
    model_base: list = field(default_factory=list)
    model_base_templates: str = ""
    population_groups: list = field(default_factory=list)
    include_list: list = field(default_factory=list)
    age_group: dict = field(default_factory=dict)

    def set_attribute(self: Self, attribute_name: str, value: Any) -> None:
        """Setter for the attributes of this class.

        :param attribute_name: Name of the attribute in IntermediateRepresentation to be set.
        :param value: Value the attribute is set to. Needs to be of right type.
        :param self: Self: 
        :param attribute_name: str: 
        :param value: Any: 

        """
        self.__setattr__(attribute_name, value)

    def check_model_base(self: Self) -> None:
        """
        Check if the model_base is set. If not, set it to the model_class.
        """
        if len(self.model_base) > 0:
            self.model_base_templates = ", ".join(
                entry[0] for entry in self.model_base if len(entry) > 0
            )
        else:
            raise IndexError("model_base is empty. No base classes found.")

    def check_complete_data(self: Self, optional: Dict
                            [str, Union[str, bool]]) -> None:
        """Check for missing data in the IntermediateRepresentation.
        Called by the Scanner as last step of the data extraction.

        :param optional: Dictionary of the optional data from the config.json. (Default value = Dict[str, Union[str, bool]])
        :param self: Self: 


        """
        assert (self.model_class != None), "set a model name"
        assert (self.namespace != None), "set a model namespace"
        assert (self.parameterset != None), "set a parameterset"
        if optional.get("parameterset_wrapper"):
            assert (self.parameterset_wrapper !=
                    None), "No parameterset_wrapper found. If None is used in this model, then set parameterset_wrapper to false in config.json"
        if optional.get("age_group"):
            assert (self.age_group != {
            }), "No age_group found. If None is used in this model set age_group to false in config.json"
