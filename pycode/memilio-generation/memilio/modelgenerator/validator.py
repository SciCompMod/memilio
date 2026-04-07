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
Validation of a raw YAML dictionary before it is converted to `ModelConfig`.

All errors are collected and raised together as a single
`ValidationError` so the user sees the full list at once.
"""

from __future__ import annotations

from typing import Any, Dict, List

from .schema import ParameterType, TransitionType


class ValidationError(Exception):
    """Raised when one or more validation errors are found."""

    def __init__(self, errors: list[str]):
        self.errors = errors
        bullet_list = "\n".join(f"  • {e}" for e in errors)
        super().__init__(f"Model configuration is invalid:\n{bullet_list}")


class Validator:
    """
    Validates a raw dictionary loaded from YAML.

    Usage::

        Validator.validate(raw_dict)   # raises ValidationError on failure
    """

    @staticmethod
    def validate(data: dict[str, Any]) -> None:
        """
        Validate *data* and raise `ValidationError` if any problem
        is found.

        Parameters
        ----------
        data:
            Dictionary as returned by ``yaml.safe_load``.
        """
        errors: list[str] = []

        # model
        model = data.get("model")
        if not isinstance(model, dict):
            errors.append("'model' section is missing or not a mapping.")
        else:
            for key in ("name", "namespace", "prefix"):
                if not isinstance(
                        model.get(key),
                        str) or not model[key].strip():
                    errors.append(f"'model.{key}' must be a non-empty string.")

        # infection_states
        states = data.get("infection_states")
        if not isinstance(states, list) or len(states) < 2:
            errors.append(
                "'infection_states' must be a list with at least 2 entries.")
            states = []
        else:
            for i, s in enumerate(states):
                if not isinstance(s, str) or not s.strip():
                    errors.append(
                        f"'infection_states[{i}]' must be a non-empty string.")
            if len(states) != len(set(states)):
                errors.append("'infection_states' contains duplicate entries.")

        state_set = set(states)

        # parameters
        params = data.get("parameters")
        if not isinstance(params, list) or len(params) == 0:
            errors.append("'parameters' must be a non-empty list.")
            params = []

        param_names: list[str] = []
        for i, p in enumerate(params):
            loc = f"parameters[{i}]"
            if not isinstance(p, dict):
                errors.append(f"'{loc}' must be a mapping.")
                continue

            name = p.get("name")
            if not isinstance(name, str) or not name.strip():
                errors.append(f"'{loc}.name' must be a non-empty string.")
            else:
                param_names.append(name)

            if not isinstance(p.get("description"), str):
                errors.append(f"'{loc}.description' must be a string.")

            ptype = p.get("type")
            if ptype not in ParameterType.ALL:
                errors.append(
                    f"'{loc}.type' must be one of {ParameterType.ALL}, got {ptype!r}."
                )

            default = p.get("default")
            if not isinstance(default, (int, float)):
                errors.append(f"'{loc}.default' must be a number.")

            bounds = p.get("bounds")
            if bounds is not None:
                if not (
                    isinstance(bounds, (list, tuple))
                    and len(bounds) == 2
                    and all(b is None or isinstance(b, (int, float)) for b in bounds)
                ):
                    errors.append(
                        f"'{loc}.bounds' must be a list of two numbers or null, e.g. [0.0, 1.0]."
                    )

        if len(param_names) != len(set(param_names)):
            errors.append("'parameters' contains duplicate 'name' entries.")

        param_name_set = set(param_names)

        # transitions
        transitions = data.get("transitions")
        if not isinstance(transitions, list) or len(transitions) == 0:
            errors.append("'transitions' must be a non-empty list.")
            transitions = []

        for i, t in enumerate(transitions):
            loc = f"transitions[{i}]"
            if not isinstance(t, dict):
                errors.append(f"'{loc}' must be a mapping.")
                continue

            from_state = t.get("from")
            to_state = t.get("to")
            ttype = t.get("type")

            if from_state not in state_set:
                errors.append(
                    f"'{loc}.from' references unknown state {from_state!r}."
                )
            if to_state not in state_set:
                errors.append(
                    f"'{loc}.to' references unknown state {to_state!r}."
                )
            if from_state == to_state and from_state is not None:
                errors.append(
                    f"'{loc}': 'from' and 'to' must differ (got {from_state!r})."
                )

            if ttype not in TransitionType.ALL:
                errors.append(
                    f"'{loc}.type' must be one of {TransitionType.ALL}, got {ttype!r}."
                )
                continue

            if ttype in (TransitionType.INFECTION, TransitionType.LINEAR):
                param = t.get("parameter")
                if param not in param_name_set:
                    errors.append(
                        f"'{loc}.parameter' references unknown parameter {param!r}."
                    )

            if ttype == TransitionType.INFECTION:
                has_singular = "infectious_state" in t
                has_plural = "infectious_states" in t
                if has_singular and has_plural:
                    errors.append(
                        f"'{loc}' must define only one of 'infectious_state' or 'infectious_states'."
                    )
                    continue

                key = "infectious_states" if has_plural else "infectious_state"
                inf_raw = t.get(key)

                if inf_raw is None:
                    errors.append(
                        f"'{loc}.{key}' must be provided for infection transitions."
                    )
                    continue

                if isinstance(inf_raw, str):
                    inf_states = [inf_raw]
                elif isinstance(inf_raw, list):
                    if len(inf_raw) == 0:
                        errors.append(f"'{loc}.{key}' must not be empty.")
                        continue
                    inf_states = []
                    for j, s in enumerate(inf_raw):
                        if not isinstance(s, str) or not s.strip():
                            errors.append(
                                f"'{loc}.{key}[{j}]' must be a non-empty string."
                            )
                        else:
                            inf_states.append(s)
                    if len(inf_states) != len(set(inf_states)):
                        errors.append(
                            f"'{loc}.{key}' contains duplicate entries."
                        )
                else:
                    errors.append(
                        f"'{loc}.{key}' must be a string or a non-empty list of strings."
                    )
                    continue

                for s in inf_states:
                    if s not in state_set:
                        errors.append(
                            f"'{loc}.{key}' references unknown state {s!r}."
                        )

        if errors:
            raise ValidationError(errors)
