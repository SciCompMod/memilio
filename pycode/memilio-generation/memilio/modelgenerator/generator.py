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
Core generator: parses a YAML config, builds the internal `ModelConfig`
representation, and renders all Jinja2 templates into strings.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Optional
import re

import yaml

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

if sys.version_info >= (3, 9):
    import importlib.resources as importlib_resources
else:
    import importlib_resources

from jinja2 import Environment, PackageLoader, StrictUndefined

from .schema import (
    ModelConfig,
    ModelMeta,
    ParameterConfig,
    ParameterType,
    TransitionConfig,
    TransitionType,
)
from .validator import Validator


class Generator:
    """
    Parses a YAML configuration and renders all model source files.

    Parameters
    ----------
    config:
        Fully-validated `ModelConfig` instance.
    """

    def __init__(self, config: ModelConfig):
        self._config = config
        self._env = Environment(
            loader=PackageLoader("memilio.modelgenerator", "templates"),
            undefined=StrictUndefined,
            keep_trailing_newline=True,
            trim_blocks=True,
            lstrip_blocks=True,
        )

    @classmethod
    def from_yaml(cls, yaml_path: str | Path) -> Generator:
        """
        Build a `Generator` from a YAML file.

        Parameters
        ----------
        yaml_path:
            Path to the ``.yaml`` configuration file.

        Returns
        -------
        Generator
        """
        with open(yaml_path, encoding="utf-8") as fh:
            raw = yaml.safe_load(fh)

        Validator.validate(raw)
        config = cls._parse(raw)
        return cls(config)

    @classmethod
    def from_toml(cls, toml_path: str | Path) -> Generator:
        """
        Build a `Generator` from a TOML file.

        Parameters
        ----------
        toml_path:
            Path to the ``.toml`` configuration file.

        Returns
        -------
        Generator
        """
        with open(toml_path, "rb") as fh:
            raw = tomllib.load(fh)

        Validator.validate(raw)
        config = cls._parse(raw)
        return cls(config)

    @classmethod
    def from_dict(cls, raw: dict) -> Generator:
        """
        Build a `Generator` from an already-loaded dictionary.

        Parameters
        ----------
        raw:
            Dictionary as returned by ``yaml.safe_load``.
        """
        Validator.validate(raw)
        config = cls._parse(raw)
        return cls(config)

    def render(self) -> dict[str, str]:
        """
        Render all new files and return a mapping of

        ``relative_output_path  to  file_content``

        The paths are relative to the MEmilio repository root.
        Use `render_patches` for the in-place edits to existing
        CMakeLists files.
        """
        cfg = self._config
        prefix = cfg.meta.prefix

        return {
            f"cpp/models/{prefix}/infection_state.h": self._render("infection_state_h.jinja2"),
            f"cpp/models/{prefix}/parameters.h": self._render("parameters_h.jinja2"),
            f"cpp/models/{prefix}/model.h": self._render("model_h.jinja2"),
            f"cpp/models/{prefix}/model.cpp": self._render("model_cpp.jinja2"),
            f"cpp/models/{prefix}/CMakeLists.txt": self._render("CMakeLists_model_txt.jinja2"),
            (
                f"pycode/memilio-simulation/memilio/simulation/bindings/models/{prefix}.cpp"
            ): self._render("pybindings_cpp.jinja2"),
            f"pycode/examples/simulation/{prefix}_simple.py": self._render("example_py.jinja2"),
            (
                f"pycode/memilio-simulation/memilio/simulation/{cfg.meta.namespace}.py"
            ): self._render("simulation_py.jinja2"),
        }

    _CPP_CMAKE = "cpp/CMakeLists.txt"
    _SIM_CMAKE = "pycode/memilio-simulation/CMakeLists.txt"
    _SIM_INIT = "pycode/memilio-simulation/memilio/simulation/__init__.py"

    def render_patches(self, output_dir: Path) -> dict[str, str | None]:
        """
        Compute the patched content of the two existing CMakeLists files.

        Returns a dict ``{relative_path: new_content | None}`` where
        ``None`` means the entry is already present (no change needed).
        """
        prefix = self._config.meta.prefix
        namespace = self._config.meta.namespace
        results: dict[str, str | None] = {}

        # cpp/CMakeLists.txt
        cpp_cmake = output_dir / self._CPP_CMAKE
        if cpp_cmake.exists():
            text = cpp_cmake.read_text(encoding="utf-8")
            entry = f"    add_subdirectory(models/{prefix})"
            if entry not in text:
                # Insert after the last add_subdirectory(models/…) line
                pattern = r"(    add_subdirectory\(models/[^)]+\))(?!.*add_subdirectory\(models/)"
                m = re.search(pattern, text, re.DOTALL)
                if m:
                    insert_at = m.end()
                    text = text[:insert_at] + "\n" + entry + text[insert_at:]
                results[self._CPP_CMAKE] = text
            else:
                results[self._CPP_CMAKE] = None  # already present

        # pycode/memilio-simulation/CMakeLists.txt
        sim_cmake = output_dir / self._SIM_CMAKE
        if sim_cmake.exists():
            text = sim_cmake.read_text(encoding="utf-8")
            module_name = f"_simulation_{namespace}"
            block = (
                f"add_pymio_module({module_name}\n"
                f"    LINKED_LIBRARIES memilio {prefix}\n"
                f"    SOURCES memilio/simulation/bindings/models/{prefix}.cpp\n"
                f")")
            if f"add_pymio_module({module_name}\n" not in text:
                # Insert before the "# install all shared" comment
                marker = "# install all shared memilio libraries"
                text = text.replace(marker, block + "\n\n" + marker)
                results[self._SIM_CMAKE] = text
            else:
                results[self._SIM_CMAKE] = None  # already present

        # pycode/memilio-simulation/memilio/simulation/__init__.py
        sim_init = output_dir / self._SIM_INIT
        if sim_init.exists():
            text = sim_init.read_text(encoding="utf-8")
            lazy_entry = (
                f"    elif attr == \"{namespace}\":\n"
                f"        import memilio.simulation.{namespace} as {namespace}\n"
                f"        return {namespace}\n"
            )
            if f'attr == "{namespace}"' not in text:
                text = text.replace(
                    "    raise AttributeError",
                    lazy_entry + "    raise AttributeError"
                )
                results[self._SIM_INIT] = text
            else:
                results[self._SIM_INIT] = None  # already present

        return results

    def write(self, output_dir: str | Path, overwrite: bool = False) -> None:
        """
        Write all rendered files under *output_dir* and patch the two
        existing CMakeLists files.

        Directories are created as needed.

        Parameters
        ----------
        output_dir:
            Root of the MEmilio repository (or any target directory).
        overwrite:
            If ``False`` (default) and the model directory
            ``cpp/models/<prefix>`` already exists, an error is raised to
            prevent accidentally overwriting an existing handwritten model.
            Set to ``True`` to allow overwriting.

        Raises
        ------
        FileExistsError
            When *overwrite* is ``False`` and the target model directory
            already exists.
        """
        output_dir = Path(output_dir)
        prefix = self._config.meta.prefix
        model_dir = output_dir / "cpp" / "models" / prefix
        if model_dir.exists() and not overwrite:
            raise FileExistsError(
                f"Model directory already exists: {model_dir}\n"
                f"Pass overwrite=True (or --force on the CLI) to overwrite it."
            )

        for rel_path, content in self.render().items():
            target = output_dir / rel_path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content, encoding="utf-8")
            print(f"  wrote  {rel_path}")

        for rel_path, content in self.render_patches(output_dir).items():
            if content is None:
                print(f"  skip   {rel_path}  (entry already present)")
            else:
                (output_dir / rel_path).write_text(content, encoding="utf-8")
                print(f"  patched {rel_path}")

    def _render(self, template_name: str) -> str:
        tmpl = self._env.get_template(template_name)
        return tmpl.render(cfg=self._config)

    @staticmethod
    def _parse(raw: dict) -> ModelConfig:
        meta = ModelMeta(
            name=raw["model"]["name"],
            namespace=raw["model"]["namespace"],
            prefix=raw["model"]["prefix"],
        )

        states: list[str] = raw["infection_states"]

        parameters = []
        for p in raw["parameters"]:
            bounds_raw = p.get("bounds")
            if bounds_raw is not None:
                bounds = (bounds_raw[0], bounds_raw[1])
            else:
                if p["type"] == ParameterType.PROBABILITY:
                    bounds = (0.0, 1.0)
                elif p["type"] == ParameterType.TIME:
                    bounds = (1e-1, None)
                else:
                    bounds = (None, None)

            parameters.append(
                ParameterConfig(
                    name=p["name"],
                    description=p.get("description", ""),
                    type=p["type"],
                    default=float(p["default"]),
                    per_age_group=bool(p.get("per_age_group", True)),
                    bounds=bounds,
                )
            )

        transitions = []
        for t in raw["transitions"]:
            raw_infectious_states = t.get("infectious_states")
            raw_infectious_state = t.get("infectious_state")
            if isinstance(raw_infectious_states, list):
                infectious_states = list(raw_infectious_states)
            elif raw_infectious_states is not None:
                infectious_states = [raw_infectious_states]
            elif isinstance(raw_infectious_state, list):
                infectious_states = list(raw_infectious_state)
            elif raw_infectious_state is not None:
                infectious_states = [raw_infectious_state]
            else:
                infectious_states = []

            transitions.append(
                TransitionConfig(
                    from_state=t["from"],
                    to_state=t["to"],
                    type=t["type"],
                    parameter=t.get("parameter"),
                    infectious_state=infectious_states[0] if infectious_states else None,
                    infectious_states=infectious_states,
                    custom_formula=t.get("custom_formula"),
                )
            )

        return ModelConfig(
            meta=meta,
            infection_states=states,
            parameters=parameters,
            transitions=transitions,
        )
