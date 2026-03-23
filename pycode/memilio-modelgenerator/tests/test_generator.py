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

import os
import tempfile
import textwrap
import unittest

from memilio.modelgenerator import Generator
from memilio.modelgenerator.validator import ValidationError

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES_DIR = os.path.join(HERE, "..", "..", "examples", "modelgenerator")

SEIR_YAML = os.path.join(EXAMPLES_DIR, "seir.yaml")
SEIRD_YAML = os.path.join(EXAMPLES_DIR, "seird.yaml")


def _render(yaml_path: str) -> dict:
    return Generator.from_yaml(yaml_path).render()


class TestParsing(unittest.TestCase):

    def test_seir_meta(self):
        gen = Generator.from_yaml(SEIR_YAML)
        cfg = gen._config
        self.assertEqual(cfg.meta.name, "SEIR")
        self.assertEqual(cfg.meta.namespace, "oseir")
        self.assertEqual(cfg.meta.prefix, "ode_seir")

    def test_seir_states(self):
        gen = Generator.from_yaml(SEIR_YAML)
        self.assertEqual(gen._config.infection_states,
                         ["Susceptible", "Exposed", "Infected", "Recovered"])

    def test_seir_parameters(self):
        gen = Generator.from_yaml(SEIR_YAML)
        names = [p.name for p in gen._config.parameters]
        self.assertIn("TransmissionProbabilityOnContact", names)
        self.assertIn("TimeExposed", names)
        self.assertIn("TimeInfected", names)

    def test_seir_transitions(self):
        gen = Generator.from_yaml(SEIR_YAML)
        types = [t.type for t in gen._config.transitions]
        self.assertIn("infection", types)
        self.assertIn("linear", types)

    def test_seird_has_custom_transition(self):
        gen = Generator.from_yaml(SEIRD_YAML)
        custom = [t for t in gen._config.transitions if t.type == "custom"]
        self.assertEqual(len(custom), 1)
        self.assertEqual(custom[0].from_state, "Infected")
        self.assertEqual(custom[0].to_state, "Dead")

    def test_has_infection_transition_flag(self):
        gen = Generator.from_yaml(SEIR_YAML)
        self.assertTrue(gen._config.has_infection_transition)


class TestInfectionStateTemplate(unittest.TestCase):

    def setUp(self):
        self.files = _render(SEIR_YAML)
        self.content = self.files["cpp/models/ode_seir/infection_state.h"]

    def test_include_guard(self):
        self.assertIn("#ifndef ODE_SEIR_INFECTIONSTATE_H", self.content)
        self.assertIn("#define ODE_SEIR_INFECTIONSTATE_H", self.content)
        self.assertIn("#endif // ODE_SEIR_INFECTIONSTATE_H", self.content)

    def test_namespace(self):
        self.assertIn("namespace oseir", self.content)

    def test_all_states_present(self):
        for state in [
                "Susceptible", "Exposed", "Infected", "Recovered", "Count"]:
            self.assertIn(state, self.content)

    def test_enum_class(self):
        self.assertIn("enum class InfectionState", self.content)


class TestParametersTemplate(unittest.TestCase):

    def setUp(self):
        self.files = _render(SEIR_YAML)
        self.content = self.files["cpp/models/ode_seir/parameters.h"]

    def test_include_guard(self):
        self.assertIn("#ifndef ODE_SEIR_PARAMETERS_H", self.content)

    def test_parameter_structs(self):
        for name in [
                "TransmissionProbabilityOnContact", "TimeExposed", "TimeInfected"]:
            self.assertIn(f"struct {name}", self.content)

    def test_contact_patterns_added(self):
        self.assertIn("struct ContactPatterns", self.content)

    def test_parameters_base(self):
        self.assertIn("using ParametersBase =", self.content)
        self.assertIn("ContactPatterns<FP>", self.content)

    def test_parameters_class(self):
        self.assertIn(
            "class Parameters : public ParametersBase<FP>", self.content)
        self.assertIn("apply_constraints", self.content)
        self.assertIn("check_constraints", self.content)

    def test_probability_constraint(self):
        self.assertIn(
            "TransmissionProbabilityOnContact<FP>>()[i] < 0.0", self.content)

    def test_time_constraint(self):
        self.assertIn("tol_times", self.content)


class TestModelTemplate(unittest.TestCase):

    def setUp(self):
        self.files = _render(SEIR_YAML)
        self.content = self.files["cpp/models/ode_seir/model.h"]

    def test_include_guard(self):
        self.assertIn("#ifndef ODE_SEIR_MODEL_H", self.content)

    def test_flows_typelist(self):
        self.assertIn("using Flows = TypeList<", self.content)
        self.assertIn(
            "Flow<InfectionState::Susceptible, InfectionState::Exposed>",
            self.content)
        self.assertIn(
            "Flow<InfectionState::Infected, InfectionState::Recovered>",
            self.content)

    def test_get_flows_method(self):
        self.assertIn("void get_flows(", self.content)

    def test_infection_flow_uses_contact_matrix(self):
        self.assertIn("ContactPatterns<FP>", self.content)
        self.assertIn("get_cont_freq_mat", self.content)

    def test_linear_flows(self):
        self.assertIn("TimeExposed<FP>>()[i]", self.content)
        self.assertIn("TimeInfected<FP>>()[i]", self.content)

    def test_serialize_deserialize(self):
        self.assertIn("void serialize(", self.content)
        self.assertIn("static IOResult<Model> deserialize(", self.content)

    def test_index_variables_for_all_states(self):
        for state in ["Susceptible", "Exposed", "Infected", "Recovered"]:
            self.assertIn(f"idx_{state}_i", self.content)


class TestPybindingsTemplate(unittest.TestCase):

    def setUp(self):
        self.files = _render(SEIR_YAML)
        key = "pycode/memilio-simulation/memilio/simulation/bindings/models/ode_seir.cpp"
        self.content = self.files[key]

    def test_module_name(self):
        self.assertIn("PYBIND11_MODULE(_simulation_oseir, m)", self.content)

    def test_enum_values(self):
        for state in ["Susceptible", "Exposed", "Infected", "Recovered"]:
            self.assertIn(f'.value("{state}"', self.content)

    def test_simulate_functions(self):
        self.assertIn('m.def("simulate"', self.content)
        self.assertIn('m.def("simulate_flows"', self.content)

    def test_model_init(self):
        self.assertIn("py::init<int>()", self.content)


class TestCMakeTemplate(unittest.TestCase):

    def setUp(self):
        self.files = _render(SEIR_YAML)
        self.content = self.files["cpp/models/ode_seir/CMakeLists.txt"]

    def test_add_library(self):
        self.assertIn("add_library(ode_seir", self.content)

    def test_source_files(self):
        for src in ["infection_state.h", "parameters.h", "model.h",
                    "model.cpp"]:
            self.assertIn(src, self.content)

    def test_link_libraries(self):
        self.assertIn(
            "target_link_libraries(ode_seir PUBLIC memilio)", self.content)


class TestValidation(unittest.TestCase):

    def _base(self):
        return {
            "model": {"name": "X", "namespace": "ox", "prefix": "ode_x"},
            "infection_states": ["S", "I"],
            "parameters": [
                {"name": "Rate", "description": "d",
                    "type": "time", "default": 1.0}
            ],
            "transitions": [
                {"from": "S", "to": "I", "type": "linear", "parameter": "Rate"}
            ],
        }

    def test_missing_model_section(self):
        d = self._base()
        del d["model"]
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_unknown_state_in_transition(self):
        d = self._base()
        d["transitions"][0]["from"] = "X_unknown"
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_unknown_parameter_in_transition(self):
        d = self._base()
        d["transitions"][0]["parameter"] = "NoSuchParam"
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_invalid_transition_type(self):
        d = self._base()
        d["transitions"][0]["type"] = "magic"
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_duplicate_states(self):
        d = self._base()
        d["infection_states"] = ["S", "S"]
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_self_loop_transition(self):
        d = self._base()
        d["transitions"][0]["to"] = "S"
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)


# CMakeLists patching
_CPP_CMAKE_STUB = """\
if(MEMILIO_BUILD_MODELS)
    add_subdirectory(models/ode_sir)
    add_subdirectory(models/ode_seir)
    add_subdirectory(models/ode_mseirs4)
endif()
"""

_SIM_CMAKE_STUB = """\
add_pymio_module(_simulation_oseir
    LINKED_LIBRARIES memilio ode_seir
    SOURCES memilio/simulation/bindings/models/oseir.cpp
)

# install all shared memilio libraries, which were given as "LINKED_LIBRARIES" to add_pymio_module
list(REMOVE_DUPLICATES PYMIO_MEMILIO_LIBS_LIST)
"""


class TestCMakePatching(unittest.TestCase):

    def _make_repo(self, cpp_cmake=_CPP_CMAKE_STUB, sim_cmake=_SIM_CMAKE_STUB):
        """Create a temporary directory tree that looks like a minimal MEmilio repo."""
        tmp = tempfile.mkdtemp()
        cpp_dir = os.path.join(tmp, "cpp")
        sim_dir = os.path.join(tmp, "pycode", "memilio-simulation")
        os.makedirs(cpp_dir)
        os.makedirs(sim_dir)
        with open(os.path.join(cpp_dir, "CMakeLists.txt"), "w") as f:
            f.write(cpp_cmake)
        with open(os.path.join(sim_dir, "CMakeLists.txt"), "w") as f:
            f.write(sim_cmake)
        return tmp

    def _gen(self):
        return Generator.from_yaml(SEIR_YAML)

    def test_cpp_cmake_gets_patched(self):
        # Use a model with a different prefix so it's not already present
        d = {
            "model": {"name": "SIR", "namespace": "osir_new", "prefix": "ode_sir_new"},
            "infection_states": ["Susceptible", "Infected", "Recovered"],
            "parameters": [
                {"name": "TransmissionRate", "description": "rate",
                    "type": "probability", "default": 0.3},
                {"name": "RecoveryTime", "description": "time",
                    "type": "time", "default": 7.0},
            ],
            "transitions": [
                {"from": "Susceptible", "to": "Infected", "type": "infection",
                 "parameter": "TransmissionRate", "infectious_state": "Infected"},
                {"from": "Infected", "to": "Recovered",
                    "type": "linear", "parameter": "RecoveryTime"},
            ],
        }
        gen = Generator.from_dict(d)
        tmp = self._make_repo()
        from pathlib import Path
        patches = gen.render_patches(Path(tmp))
        patched = patches[gen._CPP_CMAKE]
        self.assertIsNotNone(patched)
        assert patched is not None
        self.assertIn("add_subdirectory(models/ode_sir_new)", patched)
        # Original entries must still be there
        self.assertIn("add_subdirectory(models/ode_seir)", patched)

    def test_cpp_cmake_no_duplicate(self):
        """If the entry is already in the file, render_patches returns None."""
        gen = self._gen()  # prefix = ode_seir, already in stub
        tmp = self._make_repo()
        from pathlib import Path
        patches = gen.render_patches(Path(tmp))
        self.assertIsNone(patches[gen._CPP_CMAKE])

    def test_sim_cmake_gets_patched(self):
        d = {
            "model": {"name": "SIR", "namespace": "osir_new", "prefix": "ode_sir_new"},
            "infection_states": ["Susceptible", "Infected", "Recovered"],
            "parameters": [
                {"name": "TransmissionRate", "description": "rate",
                    "type": "probability", "default": 0.3},
                {"name": "RecoveryTime", "description": "time",
                    "type": "time", "default": 7.0},
            ],
            "transitions": [
                {"from": "Susceptible", "to": "Infected", "type": "infection",
                 "parameter": "TransmissionRate", "infectious_state": "Infected"},
                {"from": "Infected", "to": "Recovered",
                    "type": "linear", "parameter": "RecoveryTime"},
            ],
        }
        gen = Generator.from_dict(d)
        tmp = self._make_repo()
        from pathlib import Path
        patches = gen.render_patches(Path(tmp))
        patched = patches[gen._SIM_CMAKE]
        self.assertIsNotNone(patched)
        assert patched is not None
        self.assertIn("add_pymio_module(_simulation_osir_new", patched)
        self.assertIn("LINKED_LIBRARIES memilio ode_sir_new", patched)
        self.assertIn(
            "SOURCES memilio/simulation/bindings/models/ode_sir_new.cpp",
            patched)
        # Original content must still be there
        self.assertIn("_simulation_oseir", patched)
        self.assertIn("list(REMOVE_DUPLICATES", patched)

    def test_sim_cmake_no_duplicate(self):
        """If the module is already registered, render_patches returns None."""
        gen = self._gen()  # namespace = oseir -> _simulation_oseir, already in stub
        tmp = self._make_repo()
        from pathlib import Path
        patches = gen.render_patches(Path(tmp))
        self.assertIsNone(patches[gen._SIM_CMAKE])


if __name__ == "__main__":
    unittest.main()
