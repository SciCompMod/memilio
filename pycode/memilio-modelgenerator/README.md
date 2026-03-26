# MEmilio Model Generator

This package provides an automatic code generator for ODE compartment models in the [MEmilio](https://github.com/SciCompMod/memilio) C++ library. Given a YAML configuration file describing infection states, parameters, and transitions, it generates:

- `cpp/models/<prefix>/infection_state.h` – C++ enum for compartments
- `cpp/models/<prefix>/parameters.h` – C++ parameter structs and `Parameters` class
- `cpp/models/<prefix>/model.h` – C++ `FlowModel` with `get_flows()`
- `cpp/models/<prefix>/model.cpp` – C++ translation unit
- `cpp/models/<prefix>/CMakeLists.txt` – CMake target definition
- `pycode/memilio-simulation/memilio/simulation/bindings/models/<prefix>.cpp` – pybind11 module

## Installation

```bash
pip install -e .[dev]
```

## Usage

### Command Line

```bash
memilio-modelgenerator path/to/my_model.yaml
```

By default, the files are written to the MEmilio repository root (auto-detected from the package location). Use `--output-dir` to override:

```bash
memilio-modelgenerator path/to/my_model.yaml --output-dir /path/to/memilio
```

Use `--preview` to print generated files without writing them to disk:

```bash
memilio-modelgenerator path/to/my_model.yaml --preview
```

### Python API

```python
from memilio.modelgenerator import Generator

gen = Generator.from_yaml("examples/seir.yaml")
files = gen.render()          # dict: relative_path -> content
gen.write(output_dir=".")     # writes files to disk
```

## YAML Configuration Format

```yaml
model:
  name: "SEIR"          # human-readable name
  namespace: "oseir"    # C++ inner namespace  (mio::<namespace>)
  prefix: "ode_seir"    # folder and CMake target prefix

infection_states:
  - Susceptible
  - Exposed
  - Infected
  - Recovered

parameters:
  - name: TransmissionProbabilityOnContact
    description: "probability of getting infected from a contact"
    type: probability         # probability | time | custom
    default: 1.0
    per_age_group: true

  - name: TimeExposed
    description: "latent time in days"
    type: time
    default: 5.2
    per_age_group: true

  - name: TimeInfected
    description: "infectious time in days"
    type: time
    default: 6.0
    per_age_group: true

transitions:
  - from: Susceptible
    to: Exposed
    type: infection            # force-of-infection via contact matrix
    parameter: TransmissionProbabilityOnContact
    infectious_state: Infected

  - from: Exposed
    to: Infected
    type: linear               # rate = (1 / parameter) * source_compartment
    parameter: TimeExposed

  - from: Infected
    to: Recovered
    type: linear
    parameter: TimeInfected
```

### Transition types

| Type | Description |
|---|---|
| `infection` | Force-of-infection term using `ContactPatterns` and `S * I / N`. Adds `ContactPatterns` to the parameter set automatically. |
| `linear` | Simple outflow: `(1 / parameter) * source_compartment` |
| `custom` | Leaves a `TODO` comment with a placeholder for a user-defined expression. |

### Parameter types

| Type | C++ storage | Constraint check |
|---|---|---|
| `probability` | `CustomIndexArray<UncertainValue<FP>, AgeGroup>` | in `[0, 1]` |
| `time` | `CustomIndexArray<UncertainValue<FP>, AgeGroup>` | `>= 0.1` |
| `custom` | `CustomIndexArray<UncertainValue<FP>, AgeGroup>` | none |

## Examples

See the `pycode/examples/modelgenerator/` directory for `seir.yaml` and `seird.yaml`.
