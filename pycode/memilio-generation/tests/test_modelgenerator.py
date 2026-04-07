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
import unittest

from memilio.modelgenerator import Generator
from memilio.modelgenerator.validator import ValidationError

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES_DIR = os.path.join(HERE, "..", "..", "examples", "modelgenerator")

SEIR_YAML = os.path.join(EXAMPLES_DIR, "seir.yaml")
SEIRD_YAML = os.path.join(EXAMPLES_DIR, "seird.yaml")
SEIR_TOML = os.path.join(EXAMPLES_DIR, "seir.toml")


def _render(yaml_path: str) -> dict:
    return Generator.from_yaml(yaml_path).render()


# Parsing
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

    def test_parameter_defaults(self):
        gen = Generator.from_yaml(SEIR_YAML)
        by_name = {p.name: p for p in gen._config.parameters}
        self.assertAlmostEqual(by_name["TimeExposed"].default, 5.2)
        self.assertAlmostEqual(by_name["TimeInfected"].default, 6.0)
        self.assertAlmostEqual(
            by_name["TransmissionProbabilityOnContact"].default, 1.0)

    def test_parameter_bounds(self):
        gen = Generator.from_yaml(SEIR_YAML)
        by_name = {p.name: p for p in gen._config.parameters}
        prob = by_name["TransmissionProbabilityOnContact"]
        self.assertEqual(prob.bounds, (0.0, 1.0))
        time_exp = by_name["TimeExposed"]
        self.assertAlmostEqual(time_exp.bounds[0], 0.1)
        self.assertIsNone(time_exp.bounds[1])


# TOML loading
class TestTomlLoading(unittest.TestCase):

    def test_toml_parses_same_meta_as_yaml(self):
        gen_yaml = Generator.from_yaml(SEIR_YAML)
        gen_toml = Generator.from_toml(SEIR_TOML)
        self.assertEqual(gen_toml._config.meta.name,
                         gen_yaml._config.meta.name)
        self.assertEqual(gen_toml._config.meta.namespace,
                         gen_yaml._config.meta.namespace)
        self.assertEqual(gen_toml._config.meta.prefix,
                         gen_yaml._config.meta.prefix)

    def test_toml_parses_same_states(self):
        gen_yaml = Generator.from_yaml(SEIR_YAML)
        gen_toml = Generator.from_toml(SEIR_TOML)
        self.assertEqual(gen_toml._config.infection_states,
                         gen_yaml._config.infection_states)

    def test_toml_parses_same_parameters(self):
        gen_yaml = Generator.from_yaml(SEIR_YAML)
        gen_toml = Generator.from_toml(SEIR_TOML)
        names_yaml = [p.name for p in gen_yaml._config.parameters]
        names_toml = [p.name for p in gen_toml._config.parameters]
        self.assertEqual(names_toml, names_yaml)

    def test_toml_parses_same_transitions(self):
        gen_yaml = Generator.from_yaml(SEIR_YAML)
        gen_toml = Generator.from_toml(SEIR_TOML)
        types_yaml = [t.type for t in gen_yaml._config.transitions]
        types_toml = [t.type for t in gen_toml._config.transitions]
        self.assertEqual(types_toml, types_yaml)

    def test_toml_renders_identical_model_h(self):
        files_yaml = Generator.from_yaml(SEIR_YAML).render()
        files_toml = Generator.from_toml(SEIR_TOML).render()
        self.assertEqual(
            files_toml["cpp/models/ode_seir/model.h"],
            files_yaml["cpp/models/ode_seir/model.h"])

    def test_toml_renders_identical_infection_state_h(self):
        files_yaml = Generator.from_yaml(SEIR_YAML).render()
        files_toml = Generator.from_toml(SEIR_TOML).render()
        self.assertEqual(
            files_toml["cpp/models/ode_seir/infection_state.h"],
            files_yaml["cpp/models/ode_seir/infection_state.h"])


# infection_state.h template
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


# parameters.h template
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

    def test_default_values_in_get_default(self):
        # TimeExposed default = 5.2, TimeInfected = 6.0
        self.assertIn("5.2", self.content)
        self.assertIn("6.0", self.content)

    def test_no_contact_patterns_without_infection(self):
        # A model with only linear transitions must not get ContactPatterns
        d = {
            "model": {"name": "SI", "namespace": "osi", "prefix": "ode_si"},
            "infection_states": ["S", "I"],
            "parameters": [
                {"name": "Rate", "description": "d",
                    "type": "time", "default": 5.0}
            ],
            "transitions": [
                {"from": "S", "to": "I", "type": "linear", "parameter": "Rate"}
            ],
        }
        content = Generator.from_dict(d).render()[
            "cpp/models/ode_si/parameters.h"]
        self.assertNotIn("ContactPatterns", content)


# model.h template
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

    def test_infection_flow_supports_multiple_infectious_states(self):
        d = {
            "model": {"name": "SEIIR", "namespace": "oseiir", "prefix": "ode_seiir"},
            "infection_states": ["S", "E", "I1", "I2", "R"],
            "parameters": [
                {"name": "Beta", "description": "b",
                    "type": "probability", "default": 0.5},
                {"name": "TimeExposed", "description": "t",
                    "type": "time", "default": 4.0},
            ],
            "transitions": [
                {"from": "S", "to": "E", "type": "infection", "parameter": "Beta",
                 "infectious_state": ["I1", "I2"]},
                {"from": "E", "to": "R", "type": "linear",
                    "parameter": "TimeExposed"},
            ],
        }
        content = Generator.from_dict(
            d).render()["cpp/models/ode_seiir/model.h"]
        self.assertIn("pop[idx_I1_j] +", content)
        self.assertIn("pop[idx_I2_j]", content)

    def test_infection_flow_supports_infectious_states_key(self):
        d = {
            "model": {"name": "SEIIR", "namespace": "oseiir", "prefix": "ode_seiir"},
            "infection_states": ["S", "E", "I1", "I2", "R"],
            "parameters": [
                {"name": "Beta", "description": "b",
                    "type": "probability", "default": 0.5},
                {"name": "TimeExposed", "description": "t",
                    "type": "time", "default": 4.0},
            ],
            "transitions": [
                {"from": "S", "to": "E", "type": "infection", "parameter": "Beta",
                 "infectious_states": ["I1", "I2"]},
                {"from": "E", "to": "R", "type": "linear",
                    "parameter": "TimeExposed"},
            ],
        }
        content = Generator.from_dict(
            d).render()["cpp/models/ode_seiir/model.h"]
        self.assertIn("pop[idx_I1_j] +", content)
        self.assertIn("pop[idx_I2_j]", content)

    def test_linear_flows(self):
        self.assertIn("TimeExposed<FP>>()[i]", self.content)
        self.assertIn("TimeInfected<FP>>()[i]", self.content)

    def test_serialize_deserialize(self):
        self.assertIn("void serialize(", self.content)
        self.assertIn("static IOResult<Model> deserialize(", self.content)

    def test_index_variables_for_all_states(self):
        for state in ["Susceptible", "Exposed", "Infected", "Recovered"]:
            self.assertIn(f"idx_{state}_i", self.content)

    def test_seird_custom_transition_todo(self):
        content = Generator.from_yaml(SEIRD_YAML).render()[
            "cpp/models/ode_seird/model.h"]
        self.assertIn("TODO", content)
        self.assertIn("YOUR EXPRESSION HERE", content)

    def test_seird_custom_formula_hint(self):
        content = Generator.from_yaml(SEIRD_YAML).render()[
            "cpp/models/ode_seird/model.h"]
        self.assertIn("DeathRate[i] * y[idx_Infected_i]", content)


# pybindings.cpp template
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


# CMakeLists.txt template
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


# Python example template
class TestExampleTemplate(unittest.TestCase):

    def setUp(self):
        self.files = _render(SEIR_YAML)
        self.key = "pycode/examples/simulation/ode_seir_simple.py"
        self.content = self.files[self.key]

    def test_example_key_in_render(self):
        self.assertIn(self.key, self.files)

    def test_imports_numpy(self):
        self.assertIn("import numpy as np", self.content)

    def test_imports_agegroup(self):
        self.assertIn("from memilio.simulation import AgeGroup", self.content)

    def test_imports_correct_module(self):
        self.assertIn(
            "from memilio.simulation.oseir import", self.content)

    def test_imports_damping_for_infection_model(self):
        self.assertIn("from memilio.simulation import Damping", self.content)

    def test_simulate_call(self):
        self.assertIn("simulate(t0, tmax, dt, model)", self.content)

    def test_interpolate_call(self):
        self.assertIn("interpolate_simulation_result(result)", self.content)

    def test_default_parameter_values(self):
        # TransmissionProbabilityOnContact default = 1.0
        self.assertIn(
            "model.parameters.TransmissionProbabilityOnContact[A0] = 1.0",
            self.content)
        # TimeExposed default = 5.2
        self.assertIn(
            "model.parameters.TimeExposed[A0] = 5.2", self.content)
        # TimeInfected default = 6.0
        self.assertIn(
            "model.parameters.TimeInfected[A0] = 6.0", self.content)

    def test_contact_patterns_setup(self):
        self.assertIn("cont_freq_mat[0].baseline", self.content)
        self.assertIn("cont_freq_mat[0].minimum", self.content)

    def test_initial_conditions(self):
        self.assertIn("State.Exposed", self.content)
        self.assertIn("set_difference_from_total", self.content)
        self.assertIn("State.Susceptible", self.content)

    def test_print_table(self):
        self.assertIn("get_num_time_points", self.content)
        self.assertIn("get_time", self.content)
        self.assertIn("get_value", self.content)

    def test_run_simulation_function(self):
        self.assertIn("def run_simulation(", self.content)

    def test_main_guard(self):
        self.assertIn('if __name__ == "__main__":', self.content)

    def test_tmax_10_days_default(self):
        self.assertIn("tmax=10.0", self.content)

    def test_no_damping_for_linear_only_model(self):
        d = {
            "model": {"name": "SI", "namespace": "osi", "prefix": "ode_si"},
            "infection_states": ["S", "I"],
            "parameters": [
                {"name": "Rate", "description": "d",
                    "type": "time", "default": 5.0}
            ],
            "transitions": [
                {"from": "S", "to": "I", "type": "linear", "parameter": "Rate"}
            ],
        }
        content = Generator.from_dict(d).render()[
            "pycode/examples/simulation/ode_si_simple.py"]
        self.assertNotIn("Damping", content)
        self.assertNotIn("ContactPatterns", content)

    def test_seird_example_uses_second_state(self):
        content = Generator.from_yaml(SEIRD_YAML).render()[
            "pycode/examples/simulation/ode_seird_simple.py"]
        # Second state is Exposed, initial conditions should seed it
        self.assertIn("State.Exposed", content)
        self.assertIn("State.Susceptible", content)


# simulation_py
class TestSimulationPyTemplate(unittest.TestCase):

    def setUp(self):
        self.files = _render(SEIR_YAML)
        self.key = "pycode/memilio-simulation/memilio/simulation/oseir.py"
        self.content = self.files[self.key]

    def test_key_in_render(self):
        self.assertIn(self.key, self.files)

    def test_imports_compiled_module(self):
        self.assertIn(
            "from memilio.simulation._simulation_oseir import *", self.content)

    def test_namespace_in_key(self):
        # namespace drives the filename
        files = Generator.from_yaml(SEIRD_YAML).render()
        self.assertIn(
            "pycode/memilio-simulation/memilio/simulation/oseird.py", files)

    def test_seird_imports_correct_module(self):
        files = Generator.from_yaml(SEIRD_YAML).render()
        content = files["pycode/memilio-simulation/memilio/simulation/oseird.py"]
        self.assertIn(
            "from memilio.simulation._simulation_oseird import *", content)


# Validation
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

    def test_too_few_states(self):
        d = self._base()
        d["infection_states"] = ["S"]
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_missing_infectious_state_for_infection_transition(self):
        d = self._base()
        d["transitions"] = [
            {"from": "S", "to": "I", "type": "infection",
             "parameter": "Rate", "infectious_state": "Unknown"}
        ]
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_validation_error_lists_all_errors(self):
        d = self._base()
        del d["model"]
        d["infection_states"] = ["S"]
        try:
            Generator.from_dict(d)
            self.fail("Expected ValidationError")
        except ValidationError as exc:
            self.assertGreater(len(exc.errors), 1)

    def test_description_must_be_string(self):
        d = self._base()
        d["parameters"][0]["description"] = 42
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_duplicate_parameter_names(self):
        d = self._base()
        d["parameters"].append(
            {"name": "Rate", "description": "dup", "type": "time", "default": 2.0}
        )
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_infectious_state_list_accepts_multiple_states(self):
        d = self._base()
        d["infection_states"] = ["S", "I1", "I2"]
        d["transitions"] = [
            {"from": "S", "to": "I1", "type": "infection",
             "parameter": "Rate", "infectious_state": ["I1", "I2"]}
        ]
        Generator.from_dict(d)

    def test_infectious_states_key_accepts_multiple_states(self):
        d = self._base()
        d["infection_states"] = ["S", "I1", "I2"]
        d["transitions"] = [
            {"from": "S", "to": "I1", "type": "infection",
             "parameter": "Rate", "infectious_states": ["I1", "I2"]}
        ]
        Generator.from_dict(d)

    def test_empty_infectious_state_list_rejected(self):
        d = self._base()
        d["transitions"] = [
            {"from": "S", "to": "I", "type": "infection",
             "parameter": "Rate", "infectious_state": []}
        ]
        with self.assertRaises(ValidationError):
            Generator.from_dict(d)

    def test_conflicting_infectious_state_keys_rejected(self):
        d = self._base()
        d["transitions"] = [
            {"from": "S", "to": "I", "type": "infection", "parameter": "Rate",
             "infectious_state": "I", "infectious_states": ["I"]}
        ]
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

_SIM_INIT_STUB = """\
from memilio.simulation._simulation import *


def __getattr__(attr):
    if attr == "oseir":
        import memilio.simulation.oseir as oseir
        return oseir
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, attr))
"""


class TestCMakePatching(unittest.TestCase):

    def _make_repo(self, cpp_cmake=_CPP_CMAKE_STUB, sim_cmake=_SIM_CMAKE_STUB,
                   sim_init=_SIM_INIT_STUB):
        from pathlib import Path
        tmp = tempfile.mkdtemp()
        cpp_dir = os.path.join(tmp, "cpp")
        sim_dir = os.path.join(tmp, "pycode", "memilio-simulation")
        sim_pkg_dir = os.path.join(sim_dir, "memilio", "simulation")
        os.makedirs(cpp_dir)
        os.makedirs(sim_pkg_dir)
        with open(os.path.join(cpp_dir, "CMakeLists.txt"), "w") as f:
            f.write(cpp_cmake)
        with open(os.path.join(sim_dir, "CMakeLists.txt"), "w") as f:
            f.write(sim_cmake)
        with open(os.path.join(sim_pkg_dir, "__init__.py"), "w") as f:
            f.write(sim_init)
        return tmp

    def _gen(self):
        return Generator.from_yaml(SEIR_YAML)

    def test_cpp_cmake_gets_patched(self):
        from pathlib import Path
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
        patches = gen.render_patches(Path(tmp))
        patched = patches[gen._CPP_CMAKE]
        self.assertIsNotNone(patched)
        self.assertIn("add_subdirectory(models/ode_sir_new)", patched)
        self.assertIn("add_subdirectory(models/ode_seir)", patched)

    def test_cpp_cmake_no_duplicate(self):
        from pathlib import Path
        gen = self._gen()
        tmp = self._make_repo()
        patches = gen.render_patches(Path(tmp))
        self.assertIsNone(patches[gen._CPP_CMAKE])

    def test_sim_cmake_gets_patched(self):
        from pathlib import Path
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
        patches = gen.render_patches(Path(tmp))
        patched = patches[gen._SIM_CMAKE]
        self.assertIsNotNone(patched)
        self.assertIn("add_pymio_module(_simulation_osir_new", patched)
        self.assertIn("LINKED_LIBRARIES memilio ode_sir_new", patched)
        self.assertIn(
            "SOURCES memilio/simulation/bindings/models/ode_sir_new.cpp",
            patched)
        self.assertIn("_simulation_oseir", patched)
        self.assertIn("list(REMOVE_DUPLICATES", patched)

    def test_sim_cmake_no_duplicate(self):
        from pathlib import Path
        gen = self._gen()
        tmp = self._make_repo()
        patches = gen.render_patches(Path(tmp))
        self.assertIsNone(patches[gen._SIM_CMAKE])

    def test_write_creates_all_files(self):
        from pathlib import Path
        gen = self._gen()
        tmp = self._make_repo()
        gen.write(tmp)
        prefix = gen._config.meta.prefix
        expected = [
            f"cpp/models/{prefix}/infection_state.h",
            f"cpp/models/{prefix}/parameters.h",
            f"cpp/models/{prefix}/model.h",
            f"cpp/models/{prefix}/model.cpp",
            f"cpp/models/{prefix}/CMakeLists.txt",
        ]
        for rel in expected:
            self.assertTrue(
                (Path(tmp) / rel).exists(), f"Missing: {rel}")

    def test_write_raises_if_model_dir_exists(self):
        from pathlib import Path
        gen = self._gen()
        tmp = self._make_repo()
        # First write succeeds
        gen.write(tmp)
        # Second write without overwrite=True must fail
        with self.assertRaises(FileExistsError):
            gen.write(tmp)

    def test_write_overwrite_flag_allows_second_write(self):
        from pathlib import Path
        gen = self._gen()
        tmp = self._make_repo()
        gen.write(tmp)
        # Should not raise
        gen.write(tmp, overwrite=True)

    def test_sim_init_gets_patched(self):
        from pathlib import Path
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
        patches = gen.render_patches(Path(tmp))
        patched = patches[gen._SIM_INIT]
        self.assertIsNotNone(patched)
        self.assertIn('attr == "osir_new"', patched)
        self.assertIn(
            "import memilio.simulation.osir_new as osir_new", patched)
        # original entry must still be there
        self.assertIn('attr == "oseir"', patched)
        # raise AttributeError must still be there
        self.assertIn("raise AttributeError", patched)

    def test_sim_init_no_duplicate(self):
        from pathlib import Path
        gen = self._gen()  # namespace = oseir, already in stub
        tmp = self._make_repo()
        patches = gen.render_patches(Path(tmp))
        self.assertIsNone(patches[gen._SIM_INIT])


if __name__ == "__main__":
    unittest.main()
