import flags
import memilio.simulation
import numpy
from typing import ClassVar, Iterator, overload

__version__: str

class FlowSimulation:
    integrator: memilio.simulation.IntegratorCore
    def __init__(self, model: Model, t0: float = ..., dt: float = ...) -> None:
        """__init__(self: memilio.simulation.oseir.FlowSimulation, model: memilio.simulation.oseir.Model, t0: float = 0, dt: float = 0.1) -> None"""
    def advance(self, tmax: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """advance(self: memilio.simulation.oseir.FlowSimulation, tmax: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]"""
    @property
    def dt(self) -> float: ...
    @property
    def flows(self) -> memilio.simulation.TimeSeries: ...
    @property
    def model(self) -> Model: ...
    @property
    def result(self) -> memilio.simulation.TimeSeries: ...

class Index_InfectionState:
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.oseir.Index_InfectionState, value: int) -> None"""
    def __eq__(self, arg0: Index_InfectionState) -> bool:
        """__eq__(self: memilio.simulation.oseir.Index_InfectionState, arg0: memilio.simulation.oseir.Index_InfectionState) -> bool"""
    def __ne__(self, arg0: Index_InfectionState) -> bool:
        """__ne__(self: memilio.simulation.oseir.Index_InfectionState, arg0: memilio.simulation.oseir.Index_InfectionState) -> bool"""

class InfectionState:
    __members__: ClassVar[dict] = ...  # read-only
    Exposed: ClassVar[InfectionState] = ...
    Infected: ClassVar[InfectionState] = ...
    Recovered: ClassVar[InfectionState] = ...
    Susceptible: ClassVar[InfectionState] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.oseir.InfectionState, value: int) -> None"""
    @staticmethod
    def values() -> _InfectionStateValues:
        """values() -> memilio.simulation.oseir._InfectionStateValues"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object) -> int"""
    def __index__(self) -> int:
        """__index__(self: memilio.simulation.oseir.InfectionState) -> int"""
    def __int__(self) -> int:
        """__int__(self: memilio.simulation.oseir.InfectionState) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object) -> bool"""
    def __next__(self) -> InfectionState:
        """__next__(self: memilio.simulation.oseir.InfectionState) -> memilio.simulation.oseir.InfectionState"""
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Model(ModelBase):
    def __init__(self, num_agegroups: int) -> None:
        """__init__(self: memilio.simulation.oseir.Model, num_agegroups: int) -> None"""

class ModelBase:
    parameters: Parameters
    populations: Populations
    def __init__(self, arg0: Populations, arg1: Parameters) -> None:
        """__init__(self: memilio.simulation.oseir.ModelBase, arg0: memilio.simulation.oseir.Populations, arg1: memilio.simulation.oseir.Parameters) -> None"""
    def apply_constraints(self) -> bool:
        """apply_constraints(self: memilio.simulation.oseir.ModelBase) -> bool"""
    def check_constraints(self) -> bool:
        """check_constraints(self: memilio.simulation.oseir.ModelBase) -> bool"""
    def get_initial_values(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """get_initial_values(self: memilio.simulation.oseir.ModelBase) -> numpy.ndarray[numpy.float64[m, 1]]"""

class MultiIndex_PopulationsArray:
    @overload
    def __init__(self, arg0: memilio.simulation.Index_AgeGroup, arg1: Index_InfectionState) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.oseir.Index_InfectionState) -> None

        2. __init__(self: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg0: tuple) -> None
        """
    @overload
    def __init__(self, arg0: tuple) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.oseir.Index_InfectionState) -> None

        2. __init__(self: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg0: tuple) -> None
        """

class Parameters(ParametersBase):
    def __init__(self, arg0: memilio.simulation.AgeGroup) -> None:
        """__init__(self: memilio.simulation.oseir.Parameters, arg0: memilio.simulation.AgeGroup) -> None"""
    def check_constraints(self) -> bool:
        """check_constraints(self: memilio.simulation.oseir.Parameters) -> bool"""

class ParametersBase:
    ContactPatterns: memilio.simulation.UncertainContactMatrix
    TimeExposed: memilio.simulation.AgeGroupArray
    TimeInfected: memilio.simulation.AgeGroupArray
    TransmissionProbabilityOnContact: memilio.simulation.AgeGroupArray
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Populations(PopulationsArray):
    @overload
    def __init__(self, arg0: MultiIndex_PopulationsArray, arg1: float) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: float) -> None

        2. __init__(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray) -> None
        """
    @overload
    def __init__(self, arg0: MultiIndex_PopulationsArray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: float) -> None

        2. __init__(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray) -> None
        """
    def get_compartments(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """get_compartments(self: memilio.simulation.oseir.Populations) -> numpy.ndarray[numpy.float64[m, 1]]"""
    def get_group_total_AgeGroup(self, arg0: memilio.simulation.Index_AgeGroup) -> float:
        """get_group_total_AgeGroup(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.Index_AgeGroup) -> float"""
    def get_group_total_InfectionState(self, arg0: Index_InfectionState) -> float:
        """get_group_total_InfectionState(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.oseir.Index_InfectionState) -> float"""
    def get_num_compartments(self) -> int:
        """get_num_compartments(self: memilio.simulation.oseir.Populations) -> int"""
    def get_total(self) -> float:
        """get_total(self: memilio.simulation.oseir.Populations) -> float"""
    def set_difference_from_group_total_AgeGroup(self, arg0: MultiIndex_PopulationsArray, arg1: float) -> None:
        """set_difference_from_group_total_AgeGroup(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: float) -> None"""
    def set_difference_from_group_total_InfectionState(self, arg0: MultiIndex_PopulationsArray, arg1: float) -> None:
        """set_difference_from_group_total_InfectionState(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: float) -> None"""
    def set_difference_from_total(self, arg0: MultiIndex_PopulationsArray, arg1: float) -> None:
        """set_difference_from_total(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: float) -> None"""
    def set_group_total_AgeGroup(self, arg0: memilio.simulation.Index_AgeGroup, arg1: float) -> None:
        """set_group_total_AgeGroup(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.Index_AgeGroup, arg1: float) -> None"""
    def set_group_total_InfectionState(self, arg0: Index_InfectionState, arg1: float) -> None:
        """set_group_total_InfectionState(self: memilio.simulation.oseir.Populations, arg0: memilio.simulation.oseir.Index_InfectionState, arg1: float) -> None"""
    def set_total(self, arg0: float) -> None:
        """set_total(self: memilio.simulation.oseir.Populations, arg0: float) -> None"""

class PopulationsArray:
    @overload
    def __init__(self, arg0: MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __init__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray) -> None
        """
    @overload
    def __init__(self, arg0: MultiIndex_PopulationsArray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __init__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray) -> None
        """
    def get_flat_index(self, arg0: MultiIndex_PopulationsArray) -> int:
        """get_flat_index(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray) -> int"""
    def numel(self) -> int:
        """numel(self: memilio.simulation.oseir.PopulationsArray) -> int"""
    def resize(self, arg0: MultiIndex_PopulationsArray) -> None:
        """resize(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray) -> None"""
    def resize_AgeGroup(self, arg0: memilio.simulation.Index_AgeGroup) -> None:
        """resize_AgeGroup(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.Index_AgeGroup) -> None"""
    def resize_InfectionState(self, arg0: Index_InfectionState) -> None:
        """resize_InfectionState(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.Index_InfectionState) -> None"""
    def size(self) -> tuple[memilio.simulation.Index_AgeGroup, Index_InfectionState]:
        """size(self: memilio.simulation.oseir.PopulationsArray) -> Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.oseir.Index_InfectionState]"""
    def size_AgeGroup(self) -> memilio.simulation.Index_AgeGroup:
        """size_AgeGroup(self: memilio.simulation.oseir.PopulationsArray) -> memilio.simulation.Index_AgeGroup"""
    def size_InfectionState(self) -> Index_InfectionState:
        """size_InfectionState(self: memilio.simulation.oseir.PopulationsArray) -> memilio.simulation.oseir.Index_InfectionState"""
    @overload
    def __getitem__(self, arg0: MultiIndex_PopulationsArray) -> memilio.simulation.UncertainValue:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray) -> memilio.simulation.UncertainValue

        2. __getitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.oseir.Index_InfectionState]) -> memilio.simulation.UncertainValue
        """
    @overload
    def __getitem__(self, arg0: tuple[memilio.simulation.Index_AgeGroup, Index_InfectionState]) -> memilio.simulation.UncertainValue:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray) -> memilio.simulation.UncertainValue

        2. __getitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.oseir.Index_InfectionState]) -> memilio.simulation.UncertainValue
        """
    def __iter__(self) -> Iterator:
        """__iter__(self: memilio.simulation.oseir.PopulationsArray) -> Iterator"""
    @overload
    def __setitem__(self, arg0: MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.oseir.Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: tuple[memilio.simulation.Index_AgeGroup, Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.oseir.Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: memilio.simulation.UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.oseir.Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: float) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: memilio.simulation.oseir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.oseir.Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.oseir.PopulationsArray, arg0: object, arg1: float) -> None
        """

class Simulation:
    integrator: memilio.simulation.IntegratorCore
    def __init__(self, model: Model, t0: float = ..., dt: float = ...) -> None:
        """__init__(self: memilio.simulation.oseir.Simulation, model: memilio.simulation.oseir.Model, t0: float = 0, dt: float = 0.1) -> None"""
    def advance(self, tmax: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """advance(self: memilio.simulation.oseir.Simulation, tmax: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]"""
    @property
    def dt(self) -> float: ...
    @property
    def model(self) -> Model: ...
    @property
    def result(self) -> memilio.simulation.TimeSeries: ...

class _InfectionStateValues:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __iter__(self) -> InfectionState:
        """__iter__(self: memilio.simulation.oseir._InfectionStateValues) -> memilio.simulation.oseir.InfectionState"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.oseir._InfectionStateValues) -> int"""

def interpolate_ensemble_results(arg0: list[memilio.simulation.TimeSeries]) -> list[memilio.simulation.TimeSeries]:
    """interpolate_ensemble_results(arg0: List[memilio.simulation.TimeSeries]) -> List[memilio.simulation.TimeSeries]"""
@overload
def interpolate_simulation_result(ts: memilio.simulation.TimeSeries, abs_tol: float = ...) -> memilio.simulation.TimeSeries:
    """interpolate_simulation_result(*args, **kwargs)
    Overloaded function.

    1. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, abs_tol: float = 1e-14) -> memilio.simulation.TimeSeries

    2. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, interpolation_times: List[float]) -> memilio.simulation.TimeSeries
    """
@overload
def interpolate_simulation_result(ts: memilio.simulation.TimeSeries, interpolation_times: list[float]) -> memilio.simulation.TimeSeries:
    """interpolate_simulation_result(*args, **kwargs)
    Overloaded function.

    1. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, abs_tol: float = 1e-14) -> memilio.simulation.TimeSeries

    2. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, interpolation_times: List[float]) -> memilio.simulation.TimeSeries
    """
def simulate(t0: float, tmax: float, dt: float, model: Model, integrator: memilio.simulation.IntegratorCore = ...) -> memilio.simulation.TimeSeries:
    """simulate(t0: float, tmax: float, dt: float, model: memilio.simulation.oseir.Model, integrator: memilio.simulation.IntegratorCore = None) -> memilio.simulation.TimeSeries

    Simulates an ODE SEIR from t0 to tmax.
    """
def simulate_flows(t0: float, tmax: float, dt: float, model: Model, integrator: memilio.simulation.IntegratorCore = ...) -> list[memilio.simulation.TimeSeries]:
    """simulate_flows(t0: float, tmax: float, dt: float, model: memilio.simulation.oseir.Model, integrator: memilio.simulation.IntegratorCore = None) -> List[memilio.simulation.TimeSeries]

    Simulates an ODE SEIR with flows from t0 to tmax.
    """

