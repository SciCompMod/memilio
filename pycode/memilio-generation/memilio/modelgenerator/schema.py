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

import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple


# Transition types

class TransitionType:
    """Symbolic constants for the supported flow types."""
    INFECTION = "infection"
    """Force-of-infection flow using contact matrix and S*I/N."""
    LINEAR = "linear"
    """Simple outflow: (1 / parameter) * source_compartment."""
    RATE = "rate"
    """Simple outflow: parameter * source_compartment."""
    CUSTOM = "custom"
    """Placeholder. User must supply the expression manually."""

    ALL = (INFECTION, LINEAR, RATE, CUSTOM)


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
class DerivedQuantityConfig:
    """Configuration for a local formula-derived quantity."""

    name: str
    """C++ local variable name."""

    formula: str
    """Formula using state, parameter, and earlier derived quantity names."""


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

    rate: str | None = None
    """For ``type == "rate"``: formula for the source-proportional rate."""

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
    derived_quantities: list[DerivedQuantityConfig]
    transitions: list[TransitionConfig]

    @property
    def has_infection_transition(self) -> bool:
        """``True`` if at least one transition uses the force-of-infection."""
        return any(t.type == TransitionType.INFECTION for t in self.transitions)

    @property
    def all_parameters(self) -> list[ParameterConfig]:
        """
        User-defined parameter list.

        ``ContactPatterns`` is emitted separately by templates when
        ``has_infection_transition`` is true, so it is not represented as a
        ``ParameterConfig`` here.
        """
        return self.parameters

    def parameter_by_name(self, name: str) -> ParameterConfig:
        """Return the parameter config with the given name."""
        for parameter in self.parameters:
            if parameter.name == name:
                return parameter
        raise KeyError(name)

    @property
    def uses_safe_div(self) -> bool:
        """``True`` if generated formulas need the safe_div helper."""
        return any("safe_div" in d.formula for d in self.derived_quantities) or any(
            t.rate is not None and "safe_div" in t.rate for t in self.transitions)

    def formula_to_cpp(self, formula: str, index_name: str = "i") -> str:
        """Translate a formula by replacing known names with C++ expressions."""
        states = set(self.infection_states)
        parameters = {p.name: p for p in self.parameters}
        derived = {d.name for d in self.derived_quantities}

        def replace(match):
            name = match.group(0)
            if name in states:
                return f"y[idx_{name}_{index_name}]"
            if name in parameters:
                suffix = f"[{index_name}]" if parameters[name].per_age_group else ""
                return f"params.template get<{name}<FP>>(){suffix}"
            return name if name in derived or name in {"t", "safe_div"} else name

        return re.sub(r"\b[A-Za-z_][A-Za-z0-9_]*\b", replace, formula)

    def parameters_for_constraint_check(self) -> list[ParameterConfig]:
        """Return parameters that have explicit bound constraints."""
        return [
            p for p in self.parameters
            if p.type in (ParameterType.PROBABILITY, ParameterType.TIME)
        ]
