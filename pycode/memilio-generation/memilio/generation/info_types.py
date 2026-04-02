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

from dataclasses import dataclass, field


@dataclass
class method_type_info:
    """ Dataclass to represent the needed information of a method. Used for methods of classes in binding_type_info."""
    type: str
    name: str
    cursorkind: str
    namespace: str
    return_type: str
    arg_types: list[str] = field(default_factory=list)
    arg_names: list[str] = field(default_factory=list)
    parent_name: str = ""
    is_const: bool = False
    is_member: bool = False

    def signature_key(self) -> tuple:
        """ Create a unique key for the method. This can be used to identify methods."""
        return (
            self.name,
            self.cursorkind,
            tuple(self.arg_types),
            tuple(self.arg_names),
            self.parent_name
        )


@dataclass
class binding_type_info:
    """ Dataclass to represent the needed information of a class or function for binding generation."""
    type: str
    name: str
    cursorkind: str
    namespace: str
    return_type: str
    arg_types: list[str] = field(default_factory=list)
    arg_names: list[str] = field(default_factory=list)
    parent_name: str = ""
    is_const: bool = False
    is_member: bool = False
    methods: list[method_type_info] = field(default_factory=list)
    base_classes: list[str] = field(default_factory=list)
    init: list[dict[str, str]] = field(default_factory=list)
    template_args: list[str] = field(default_factory=list)

    def signature_key(self) -> tuple:
        """ Create a unique key for the binding. This can be used to identify bindings."""
        return (
            self.name,
            self.cursorkind,
            tuple(self.arg_types),
            tuple(self.arg_names),
            self.parent_name
        )
