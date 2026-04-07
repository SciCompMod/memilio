#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Henrik Zunker
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
Dataclass definitions that represent a parsed model configuration.

These are the internal representations produced by parsing a YAML file.
The `Generator` consumes these objects and passes them to the
Jinja2 templates.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple


# Transition types

class TransitionType:
    """Symbolic constants for the supported flow types."""
    INFECTION = "infection"
    """Force-of-infection flow using contact matrix and S*I/N."""
    LINEAR = "linear"
    """Simple outflow: (1 / parameter) * source_compartment."""
    CUSTOM = "custom"
    """Placeholder. User must supply the expression manually."""

    ALL = (INFECTION, LINEAR, CUSTOM)


class ParameterType:
    """Symbolic constants for built-in parameter storage types."""
    PROBABILITY = "probability"
    """Scalar in [0, 1] per age group; stored as UncertainValue."""
    TIME = "time"
    """Positive duration in days per age group; stored as UncertainValue."""
    CUSTOM = "custom"
    """User-defined; no automatic constraint check is generated."""

    ALL = (PROBABILITY, TIME, CUSTOM)


@dataclass
class ModelMeta:

    name: str
    """Human-readable model name, e.g. ``"SEIR"``."""

    namespace: str
    """Inner C++ namespace, e.g. ``"oseir"`` → ``mio::oseir``."""

    prefix: str
    """Folder and CMake target prefix, e.g. ``"ode_seir"``."""

    @property
    def guard_prefix(self) -> str:
        """Upper-case version of ``prefix`` used in include guards."""
        return self.prefix.upper()


@dataclass
class ParameterConfig:
    """Configuration for a single model parameter."""

    name: str
    """C++ struct name, e.g. ``"TransmissionProbabilityOnContact"``."""

    description: str
    """Short description used in the Doxygen comment."""

    type: str
    """One of `ParameterType.ALL`."""

    default: float
    """Default value passed to ``get_default``."""

    per_age_group: bool = True
    """If ``True`` the storage type is ``CustomIndexArray<UncertainValue<FP>, AgeGroup>``."""

    bounds: tuple[float | None, float | None] = field(
        default_factory=lambda: (None, None))
    """(lower, upper) bounds used in the constraint checks. ``None`` means unchecked."""


@dataclass
class TransitionConfig:
    """Configuration for a single compartment flow."""

    from_state: str
    """Source compartment name."""

    to_state: str
    """Target compartment name."""

    type: str
    """One of `TransitionType.ALL`."""

    parameter: str | None = None
    """Name of the `ParameterConfig` that drives this flow."""

    infectious_state: str | None = None
    """
    For ``type == "infection"``: the compartment whose population drives
    infection (typically ``"Infected"``). Kept as a compatibility alias
    for the first entry of ``infectious_states``.
    """

    infectious_states: list[str] = field(default_factory=list)
    """
    For ``type == "infection"``: list of compartments whose populations
    are summed to drive infection.
    """

    custom_formula: str | None = None
    """
    For ``type == "custom"``: an optional hint that is placed in a
    ``TODO`` comment next to the placeholder.
    """


@dataclass
class ModelConfig:
    """Complete parsed model configuration."""

    meta: ModelMeta
    infection_states: list[str]
    parameters: list[ParameterConfig]
    transitions: list[TransitionConfig]

    @property
    def has_infection_transition(self) -> bool:
        """``True`` if at least one transition uses the force-of-infection."""
        return any(t.type == TransitionType.INFECTION for t in self.transitions)

    @property
    def all_parameters(self) -> list[ParameterConfig]:
        """
        Full parameter list including the implicitly added
        ``ContactPatterns`` when any infection transition is present.
        """
        if not self.has_infection_transition:
            return self.parameters
        # ContactPatterns is added at the end; the generator inserts it
        # directly into the template so we only expose the user-defined ones
        # here.  The template accesses has_infection_transition separately.
        return self.parameters

    def parameters_for_constraint_check(self) -> list[ParameterConfig]:
        """Return parameters that have explicit bound constraints."""
        return [
            p for p in self.parameters
            if p.type in (ParameterType.PROBABILITY, ParameterType.TIME)
        ]
