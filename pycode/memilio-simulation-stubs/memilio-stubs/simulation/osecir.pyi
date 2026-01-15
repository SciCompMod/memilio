import flags
import memilio.simulation
import numpy
from typing import Callable, ClassVar, Iterator, overload

__version__: str

class ContactLocation:
    __members__: ClassVar[dict] = ...  # read-only
    Home: ClassVar[ContactLocation] = ...
    Other: ClassVar[ContactLocation] = ...
    School: ClassVar[ContactLocation] = ...
    Work: ClassVar[ContactLocation] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.osecir.ContactLocation, value: int) -> None"""
    @staticmethod
    def values() -> _ContactLocationValues:
        """values() -> memilio.simulation.osecir._ContactLocationValues"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object) -> int"""
    def __index__(self) -> int:
        """__index__(self: memilio.simulation.osecir.ContactLocation) -> int"""
    def __int__(self) -> int:
        """__int__(self: memilio.simulation.osecir.ContactLocation) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object) -> bool"""
    def __next__(self) -> ContactLocation:
        """__next__(self: memilio.simulation.osecir.ContactLocation) -> memilio.simulation.osecir.ContactLocation"""
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class EnsembleGraphResults:
    def __init__(self) -> None:
        """__init__(self: memilio.simulation.osecir.EnsembleGraphResults) -> None"""
    def __bool__(self) -> bool:
        """__bool__(self: memilio.simulation.osecir.EnsembleGraphResults) -> bool

        Check whether the list is nonempty
        """
    def __getitem__(self, arg0: int) -> MobilityGraph:
        """__getitem__(self: memilio.simulation.osecir.EnsembleGraphResults, arg0: int) -> memilio.simulation.osecir.MobilityGraph"""
    def __iter__(self) -> Iterator:
        """__iter__(self: memilio.simulation.osecir.EnsembleGraphResults) -> Iterator"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.osecir.EnsembleGraphResults) -> int"""

class FlowSimulation:
    integrator: memilio.simulation.IntegratorCore
    def __init__(self, model: Model, t0: float = ..., dt: float = ...) -> None:
        """__init__(self: memilio.simulation.osecir.FlowSimulation, model: memilio.simulation.osecir.Model, t0: float = 0, dt: float = 0.1) -> None"""
    def advance(self, tmax: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """advance(self: memilio.simulation.osecir.FlowSimulation, tmax: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]"""
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
        """__init__(self: memilio.simulation.osecir.Index_InfectionState, value: int) -> None"""
    def __eq__(self, arg0: Index_InfectionState) -> bool:
        """__eq__(self: memilio.simulation.osecir.Index_InfectionState, arg0: memilio.simulation.osecir.Index_InfectionState) -> bool"""
    def __ne__(self, arg0: Index_InfectionState) -> bool:
        """__ne__(self: memilio.simulation.osecir.Index_InfectionState, arg0: memilio.simulation.osecir.Index_InfectionState) -> bool"""

class InfectionState:
    __members__: ClassVar[dict] = ...  # read-only
    Dead: ClassVar[InfectionState] = ...
    Exposed: ClassVar[InfectionState] = ...
    InfectedCritical: ClassVar[InfectionState] = ...
    InfectedNoSymptoms: ClassVar[InfectionState] = ...
    InfectedNoSymptomsConfirmed: ClassVar[InfectionState] = ...
    InfectedSevere: ClassVar[InfectionState] = ...
    InfectedSymptoms: ClassVar[InfectionState] = ...
    InfectedSymptomsConfirmed: ClassVar[InfectionState] = ...
    Recovered: ClassVar[InfectionState] = ...
    Susceptible: ClassVar[InfectionState] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.osecir.InfectionState, value: int) -> None"""
    @staticmethod
    def values() -> _InfectionStateValues:
        """values() -> memilio.simulation.osecir._InfectionStateValues"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object) -> int"""
    def __index__(self) -> int:
        """__index__(self: memilio.simulation.osecir.InfectionState) -> int"""
    def __int__(self) -> int:
        """__int__(self: memilio.simulation.osecir.InfectionState) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object) -> bool"""
    def __next__(self) -> InfectionState:
        """__next__(self: memilio.simulation.osecir.InfectionState) -> memilio.simulation.osecir.InfectionState"""
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class MobilityGraph:
    def __init__(self) -> None:
        """__init__(self: memilio.simulation.osecir.MobilityGraph) -> None"""
    @overload
    def add_edge(self, arg0: int, arg1: int, arg2: memilio.simulation.MobilityParameters) -> memilio.simulation.MobilityEdge:
        """add_edge(*args, **kwargs)
        Overloaded function.

        1. add_edge(self: memilio.simulation.osecir.MobilityGraph, arg0: int, arg1: int, arg2: memilio.simulation.MobilityParameters) -> memilio.simulation.MobilityEdge

        2. add_edge(self: memilio.simulation.osecir.MobilityGraph, arg0: int, arg1: int, arg2: numpy.ndarray[numpy.float64[m, 1]]) -> memilio.simulation.MobilityEdge
        """
    @overload
    def add_edge(self, arg0: int, arg1: int, arg2: numpy.ndarray[numpy.float64[m, 1]]) -> memilio.simulation.MobilityEdge:
        """add_edge(*args, **kwargs)
        Overloaded function.

        1. add_edge(self: memilio.simulation.osecir.MobilityGraph, arg0: int, arg1: int, arg2: memilio.simulation.MobilityParameters) -> memilio.simulation.MobilityEdge

        2. add_edge(self: memilio.simulation.osecir.MobilityGraph, arg0: int, arg1: int, arg2: numpy.ndarray[numpy.float64[m, 1]]) -> memilio.simulation.MobilityEdge
        """
    def add_node(self, id: int, model: Model, t0: float = ..., dt: float = ...) -> SimulationNode:
        """add_node(self: memilio.simulation.osecir.MobilityGraph, id: int, model: memilio.simulation.osecir.Model, t0: float = 0.0, dt: float = 0.1) -> memilio.simulation.osecir.SimulationNode"""
    def get_edge(self, arg0: int) -> memilio.simulation.MobilityEdge:
        """get_edge(self: memilio.simulation.osecir.MobilityGraph, arg0: int) -> memilio.simulation.MobilityEdge"""
    def get_node(self, arg0: int) -> SimulationNode:
        """get_node(self: memilio.simulation.osecir.MobilityGraph, arg0: int) -> memilio.simulation.osecir.SimulationNode"""
    def get_num_out_edges(self, arg0: int) -> int:
        """get_num_out_edges(self: memilio.simulation.osecir.MobilityGraph, arg0: int) -> int"""
    def get_out_edge(self, arg0: int, arg1: int) -> memilio.simulation.MobilityEdge:
        """get_out_edge(self: memilio.simulation.osecir.MobilityGraph, arg0: int, arg1: int) -> memilio.simulation.MobilityEdge"""
    @property
    def num_edges(self) -> int: ...
    @property
    def num_nodes(self) -> int: ...

class MobilitySimulation:
    def __init__(self, graph: MobilityGraph, t0: float = ..., dt: float = ...) -> None:
        """__init__(self: memilio.simulation.osecir.MobilitySimulation, graph: memilio.simulation.osecir.MobilityGraph, t0: float = 0.0, dt: float = 1.0) -> None"""
    def advance(self, tmax: float) -> None:
        """advance(self: memilio.simulation.osecir.MobilitySimulation, tmax: float) -> None"""
    @property
    def graph(self) -> MobilityGraph: ...
    @property
    def t(self) -> float: ...

class Model(ModelBase):
    def __init__(self, num_agegroups: int) -> None:
        """__init__(self: memilio.simulation.osecir.Model, num_agegroups: int) -> None"""

class ModelBase:
    parameters: Parameters
    populations: Populations
    def __init__(self, arg0: Populations, arg1: Parameters) -> None:
        """__init__(self: memilio.simulation.osecir.ModelBase, arg0: memilio.simulation.osecir.Populations, arg1: memilio.simulation.osecir.Parameters) -> None"""
    def apply_constraints(self) -> bool:
        """apply_constraints(self: memilio.simulation.osecir.ModelBase) -> bool"""
    def check_constraints(self) -> bool:
        """check_constraints(self: memilio.simulation.osecir.ModelBase) -> bool"""
    def get_initial_values(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """get_initial_values(self: memilio.simulation.osecir.ModelBase) -> numpy.ndarray[numpy.float64[m, 1]]"""

class ModelGraph:
    def __init__(self) -> None:
        """__init__(self: memilio.simulation.osecir.ModelGraph) -> None"""
    @overload
    def add_edge(self, start_node_idx: int, end_node_idx: int, mobility_parameters: memilio.simulation.MobilityParameters) -> memilio.simulation.MobilityParameterEdge:
        """add_edge(*args, **kwargs)
        Overloaded function.

        1. add_edge(self: memilio.simulation.osecir.ModelGraph, start_node_idx: int, end_node_idx: int, mobility_parameters: memilio.simulation.MobilityParameters) -> memilio.simulation.MobilityParameterEdge

        2. add_edge(self: memilio.simulation.osecir.ModelGraph, arg0: int, arg1: int, arg2: numpy.ndarray[numpy.float64[m, 1]]) -> memilio.simulation.MobilityParameterEdge
        """
    @overload
    def add_edge(self, arg0: int, arg1: int, arg2: numpy.ndarray[numpy.float64[m, 1]]) -> memilio.simulation.MobilityParameterEdge:
        """add_edge(*args, **kwargs)
        Overloaded function.

        1. add_edge(self: memilio.simulation.osecir.ModelGraph, start_node_idx: int, end_node_idx: int, mobility_parameters: memilio.simulation.MobilityParameters) -> memilio.simulation.MobilityParameterEdge

        2. add_edge(self: memilio.simulation.osecir.ModelGraph, arg0: int, arg1: int, arg2: numpy.ndarray[numpy.float64[m, 1]]) -> memilio.simulation.MobilityParameterEdge
        """
    def add_node(self, id: int, model: Model) -> ModelNode:
        """add_node(self: memilio.simulation.osecir.ModelGraph, id: int, model: memilio.simulation.osecir.Model) -> memilio.simulation.osecir.ModelNode"""
    def get_edge(self, arg0: int) -> memilio.simulation.MobilityParameterEdge:
        """get_edge(self: memilio.simulation.osecir.ModelGraph, arg0: int) -> memilio.simulation.MobilityParameterEdge"""
    def get_node(self, arg0: int) -> ModelNode:
        """get_node(self: memilio.simulation.osecir.ModelGraph, arg0: int) -> memilio.simulation.osecir.ModelNode"""
    def get_num_out_edges(self, arg0: int) -> int:
        """get_num_out_edges(self: memilio.simulation.osecir.ModelGraph, arg0: int) -> int"""
    def get_out_edge(self, arg0: int, arg1: int) -> memilio.simulation.MobilityParameterEdge:
        """get_out_edge(self: memilio.simulation.osecir.ModelGraph, arg0: int, arg1: int) -> memilio.simulation.MobilityParameterEdge"""
    @property
    def num_edges(self) -> int: ...
    @property
    def num_nodes(self) -> int: ...

class ModelNode:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    @property
    def id(self) -> int: ...
    @property
    def property(self) -> Model: ...

class MultiIndex_PopulationsArray:
    @overload
    def __init__(self, arg0: memilio.simulation.Index_AgeGroup, arg1: Index_InfectionState) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.osecir.Index_InfectionState) -> None

        2. __init__(self: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg0: tuple) -> None
        """
    @overload
    def __init__(self, arg0: tuple) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.osecir.Index_InfectionState) -> None

        2. __init__(self: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg0: tuple) -> None
        """

class ParameterStudy:
    num_runs: int
    t0: float
    tmax: float
    @overload
    def __init__(self, model: Model, t0: float, tmax: float, num_runs: int) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.osecir.ParameterStudy, model: memilio.simulation.osecir.Model, t0: float, tmax: float, num_runs: int) -> None

        2. __init__(self: memilio.simulation.osecir.ParameterStudy, model_graph: memilio.simulation.osecir.ModelGraph, t0: float, tmax: float, dt: float, num_runs: int) -> None
        """
    @overload
    def __init__(self, model_graph: ModelGraph, t0: float, tmax: float, dt: float, num_runs: int) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.osecir.ParameterStudy, model: memilio.simulation.osecir.Model, t0: float, tmax: float, num_runs: int) -> None

        2. __init__(self: memilio.simulation.osecir.ParameterStudy, model_graph: memilio.simulation.osecir.ModelGraph, t0: float, tmax: float, dt: float, num_runs: int) -> None
        """
    @overload
    def run(self, handle_result_func: Callable[[MobilityGraph, int], None]) -> None:
        """run(*args, **kwargs)
        Overloaded function.

        1. run(self: memilio.simulation.osecir.ParameterStudy, handle_result_func: Callable[[memilio.simulation.osecir.MobilityGraph, int], None]) -> None

        2. run(self: memilio.simulation.osecir.ParameterStudy) -> memilio.simulation.osecir.EnsembleGraphResults
        """
    @overload
    def run(self) -> EnsembleGraphResults:
        """run(*args, **kwargs)
        Overloaded function.

        1. run(self: memilio.simulation.osecir.ParameterStudy, handle_result_func: Callable[[memilio.simulation.osecir.MobilityGraph, int], None]) -> None

        2. run(self: memilio.simulation.osecir.ParameterStudy) -> memilio.simulation.osecir.EnsembleGraphResults
        """
    @overload
    def run_single(self, handle_result_func: Callable[[Simulation, int], None]) -> None:
        """run_single(*args, **kwargs)
        Overloaded function.

        1. run_single(self: memilio.simulation.osecir.ParameterStudy, handle_result_func: Callable[[memilio.simulation.osecir.Simulation, int], None]) -> None

        2. run_single(self: memilio.simulation.osecir.ParameterStudy) -> List[memilio.simulation.osecir.Simulation]
        """
    @overload
    def run_single(self) -> list[Simulation]:
        """run_single(*args, **kwargs)
        Overloaded function.

        1. run_single(self: memilio.simulation.osecir.ParameterStudy, handle_result_func: Callable[[memilio.simulation.osecir.Simulation, int], None]) -> None

        2. run_single(self: memilio.simulation.osecir.ParameterStudy) -> List[memilio.simulation.osecir.Simulation]
        """
    @property
    def model(self) -> Model: ...
    @property
    def model_graph(self) -> ModelGraph: ...

class Parameters(ParametersBase):
    def __init__(self, arg0: memilio.simulation.AgeGroup) -> None:
        """__init__(self: memilio.simulation.osecir.Parameters, arg0: memilio.simulation.AgeGroup) -> None"""
    def apply_constraints(self) -> bool:
        """apply_constraints(self: memilio.simulation.osecir.Parameters) -> bool"""
    def check_constraints(self) -> bool:
        """check_constraints(self: memilio.simulation.osecir.Parameters) -> bool"""

class ParametersBase:
    ContactPatterns: memilio.simulation.UncertainContactMatrix
    CriticalPerSevere: memilio.simulation.AgeGroupArray
    DeathsPerCritical: memilio.simulation.AgeGroupArray
    DynamicNPIsImplementationDelay: memilio.simulation.UncertainValue
    DynamicNPIsInfectedSymptoms: memilio.simulation.DynamicNPIs
    ICUCapacity: memilio.simulation.UncertainValue
    MaxRiskOfInfectionFromSymptomatic: memilio.simulation.AgeGroupArray
    RecoveredPerInfectedNoSymptoms: memilio.simulation.AgeGroupArray
    RelativeTransmissionNoSymptoms: memilio.simulation.AgeGroupArray
    RiskOfInfectionFromSymptomatic: memilio.simulation.AgeGroupArray
    Seasonality: memilio.simulation.UncertainValue
    SeverePerInfectedSymptoms: memilio.simulation.AgeGroupArray
    StartDay: float
    TestAndTraceCapacity: memilio.simulation.UncertainValue
    TestAndTraceCapacityMaxRisk: memilio.simulation.UncertainValue
    TimeExposed: memilio.simulation.AgeGroupArray
    TimeInfectedCritical: memilio.simulation.AgeGroupArray
    TimeInfectedNoSymptoms: memilio.simulation.AgeGroupArray
    TimeInfectedSevere: memilio.simulation.AgeGroupArray
    TimeInfectedSymptoms: memilio.simulation.AgeGroupArray
    TransmissionProbabilityOnContact: memilio.simulation.AgeGroupArray
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Populations(PopulationsArray):
    @overload
    def __init__(self, arg0: MultiIndex_PopulationsArray, arg1: float) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: float) -> None

        2. __init__(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray) -> None
        """
    @overload
    def __init__(self, arg0: MultiIndex_PopulationsArray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: float) -> None

        2. __init__(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray) -> None
        """
    def get_compartments(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """get_compartments(self: memilio.simulation.osecir.Populations) -> numpy.ndarray[numpy.float64[m, 1]]"""
    def get_group_total_AgeGroup(self, arg0: memilio.simulation.Index_AgeGroup) -> float:
        """get_group_total_AgeGroup(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.Index_AgeGroup) -> float"""
    def get_group_total_InfectionState(self, arg0: Index_InfectionState) -> float:
        """get_group_total_InfectionState(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.osecir.Index_InfectionState) -> float"""
    def get_num_compartments(self) -> int:
        """get_num_compartments(self: memilio.simulation.osecir.Populations) -> int"""
    def get_total(self) -> float:
        """get_total(self: memilio.simulation.osecir.Populations) -> float"""
    def set_difference_from_group_total_AgeGroup(self, arg0: MultiIndex_PopulationsArray, arg1: float) -> None:
        """set_difference_from_group_total_AgeGroup(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: float) -> None"""
    def set_difference_from_group_total_InfectionState(self, arg0: MultiIndex_PopulationsArray, arg1: float) -> None:
        """set_difference_from_group_total_InfectionState(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: float) -> None"""
    def set_difference_from_total(self, arg0: MultiIndex_PopulationsArray, arg1: float) -> None:
        """set_difference_from_total(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: float) -> None"""
    def set_group_total_AgeGroup(self, arg0: memilio.simulation.Index_AgeGroup, arg1: float) -> None:
        """set_group_total_AgeGroup(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.Index_AgeGroup, arg1: float) -> None"""
    def set_group_total_InfectionState(self, arg0: Index_InfectionState, arg1: float) -> None:
        """set_group_total_InfectionState(self: memilio.simulation.osecir.Populations, arg0: memilio.simulation.osecir.Index_InfectionState, arg1: float) -> None"""
    def set_total(self, arg0: float) -> None:
        """set_total(self: memilio.simulation.osecir.Populations, arg0: float) -> None"""

class PopulationsArray:
    @overload
    def __init__(self, arg0: MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __init__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray) -> None
        """
    @overload
    def __init__(self, arg0: MultiIndex_PopulationsArray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __init__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray) -> None
        """
    def get_flat_index(self, arg0: MultiIndex_PopulationsArray) -> int:
        """get_flat_index(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray) -> int"""
    def numel(self) -> int:
        """numel(self: memilio.simulation.osecir.PopulationsArray) -> int"""
    def resize(self, arg0: MultiIndex_PopulationsArray) -> None:
        """resize(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray) -> None"""
    def resize_AgeGroup(self, arg0: memilio.simulation.Index_AgeGroup) -> None:
        """resize_AgeGroup(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.Index_AgeGroup) -> None"""
    def resize_InfectionState(self, arg0: Index_InfectionState) -> None:
        """resize_InfectionState(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.Index_InfectionState) -> None"""
    def size(self) -> tuple[memilio.simulation.Index_AgeGroup, Index_InfectionState]:
        """size(self: memilio.simulation.osecir.PopulationsArray) -> Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.osecir.Index_InfectionState]"""
    def size_AgeGroup(self) -> memilio.simulation.Index_AgeGroup:
        """size_AgeGroup(self: memilio.simulation.osecir.PopulationsArray) -> memilio.simulation.Index_AgeGroup"""
    def size_InfectionState(self) -> Index_InfectionState:
        """size_InfectionState(self: memilio.simulation.osecir.PopulationsArray) -> memilio.simulation.osecir.Index_InfectionState"""
    @overload
    def __getitem__(self, arg0: MultiIndex_PopulationsArray) -> memilio.simulation.UncertainValue:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray) -> memilio.simulation.UncertainValue

        2. __getitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.osecir.Index_InfectionState]) -> memilio.simulation.UncertainValue
        """
    @overload
    def __getitem__(self, arg0: tuple[memilio.simulation.Index_AgeGroup, Index_InfectionState]) -> memilio.simulation.UncertainValue:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray) -> memilio.simulation.UncertainValue

        2. __getitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.osecir.Index_InfectionState]) -> memilio.simulation.UncertainValue
        """
    def __iter__(self) -> Iterator:
        """__iter__(self: memilio.simulation.osecir.PopulationsArray) -> Iterator"""
    @overload
    def __setitem__(self, arg0: MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.osecir.Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: tuple[memilio.simulation.Index_AgeGroup, Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.osecir.Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: memilio.simulation.UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.osecir.Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: float) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: memilio.simulation.osecir.MultiIndex_PopulationsArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.osecir.Index_InfectionState], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.osecir.PopulationsArray, arg0: object, arg1: float) -> None
        """

class Simulation:
    integrator: memilio.simulation.IntegratorCore
    def __init__(self, model: Model, t0: float = ..., dt: float = ...) -> None:
        """__init__(self: memilio.simulation.osecir.Simulation, model: memilio.simulation.osecir.Model, t0: float = 0, dt: float = 0.1) -> None"""
    def advance(self, tmax: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """advance(self: memilio.simulation.osecir.Simulation, tmax: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]"""
    @property
    def dt(self) -> float: ...
    @property
    def model(self) -> Model: ...
    @property
    def result(self) -> memilio.simulation.TimeSeries: ...

class SimulationNode:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    @property
    def id(self) -> int: ...
    @property
    def property(self) -> Simulation: ...

class _ContactLocationValues:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __iter__(self) -> ContactLocation:
        """__iter__(self: memilio.simulation.osecir._ContactLocationValues) -> memilio.simulation.osecir.ContactLocation"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.osecir._ContactLocationValues) -> int"""

class _InfectionStateValues:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __iter__(self) -> InfectionState:
        """__iter__(self: memilio.simulation.osecir._InfectionStateValues) -> memilio.simulation.osecir.InfectionState"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.osecir._InfectionStateValues) -> int"""

def draw_sample(model: Model) -> None:
    """draw_sample(model: memilio.simulation.osecir.Model) -> None"""
def ensemble_mean(arg0: list[list[memilio.simulation.TimeSeries]]) -> list[memilio.simulation.TimeSeries]:
    """ensemble_mean(arg0: List[List[memilio.simulation.TimeSeries]]) -> List[memilio.simulation.TimeSeries]"""
def ensemble_percentile(arg0: list[list[memilio.simulation.TimeSeries]], arg1: float) -> list[memilio.simulation.TimeSeries]:
    """ensemble_percentile(arg0: List[List[memilio.simulation.TimeSeries]], arg1: float) -> List[memilio.simulation.TimeSeries]"""
@overload
def interpolate_ensemble_results(arg0: list[memilio.simulation.TimeSeries]) -> list[memilio.simulation.TimeSeries]:
    """interpolate_ensemble_results(*args, **kwargs)
    Overloaded function.

    1. interpolate_ensemble_results(arg0: List[memilio.simulation.TimeSeries]) -> List[memilio.simulation.TimeSeries]

    2. interpolate_ensemble_results(arg0: memilio.simulation.osecir.EnsembleGraphResults) -> List[List[memilio.simulation.TimeSeries]]
    """
@overload
def interpolate_ensemble_results(arg0: EnsembleGraphResults) -> list[list[memilio.simulation.TimeSeries]]:
    """interpolate_ensemble_results(*args, **kwargs)
    Overloaded function.

    1. interpolate_ensemble_results(arg0: List[memilio.simulation.TimeSeries]) -> List[memilio.simulation.TimeSeries]

    2. interpolate_ensemble_results(arg0: memilio.simulation.osecir.EnsembleGraphResults) -> List[List[memilio.simulation.TimeSeries]]
    """
@overload
def interpolate_simulation_result(ts: memilio.simulation.TimeSeries, abs_tol: float = ...) -> memilio.simulation.TimeSeries:
    """interpolate_simulation_result(*args, **kwargs)
    Overloaded function.

    1. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, abs_tol: float = 1e-14) -> memilio.simulation.TimeSeries

    2. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, interpolation_times: List[float]) -> memilio.simulation.TimeSeries

    3. interpolate_simulation_result(arg0: memilio.simulation.osecir.MobilityGraph) -> List[memilio.simulation.TimeSeries]
    """
@overload
def interpolate_simulation_result(ts: memilio.simulation.TimeSeries, interpolation_times: list[float]) -> memilio.simulation.TimeSeries:
    """interpolate_simulation_result(*args, **kwargs)
    Overloaded function.

    1. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, abs_tol: float = 1e-14) -> memilio.simulation.TimeSeries

    2. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, interpolation_times: List[float]) -> memilio.simulation.TimeSeries

    3. interpolate_simulation_result(arg0: memilio.simulation.osecir.MobilityGraph) -> List[memilio.simulation.TimeSeries]
    """
@overload
def interpolate_simulation_result(arg0: MobilityGraph) -> list[memilio.simulation.TimeSeries]:
    """interpolate_simulation_result(*args, **kwargs)
    Overloaded function.

    1. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, abs_tol: float = 1e-14) -> memilio.simulation.TimeSeries

    2. interpolate_simulation_result(ts: memilio.simulation.TimeSeries, interpolation_times: List[float]) -> memilio.simulation.TimeSeries

    3. interpolate_simulation_result(arg0: memilio.simulation.osecir.MobilityGraph) -> List[memilio.simulation.TimeSeries]
    """
def read_input_data_county(arg0: list[Model], arg1: memilio.simulation.Date, arg2: list[int], arg3: list[float], arg4: float, arg5: str, arg6: int, arg7: bool) -> None:
    """read_input_data_county(arg0: List[memilio.simulation.osecir.Model], arg1: memilio.simulation.Date, arg2: List[int], arg3: List[float], arg4: float, arg5: str, arg6: int, arg7: bool) -> None"""
def set_edges(arg0: str, arg1: ModelGraph, arg2: int) -> None:
    """set_edges(arg0: str, arg1: memilio.simulation.osecir.ModelGraph, arg2: int) -> None"""
def set_nodes(arg0: Parameters, arg1: memilio.simulation.Date, arg2: memilio.simulation.Date, arg3: str, arg4: str, arg5: bool, arg6: ModelGraph, arg7: list[float], arg8: float, arg9: float, arg10: int, arg11: bool) -> None:
    """set_nodes(arg0: memilio.simulation.osecir.Parameters, arg1: memilio.simulation.Date, arg2: memilio.simulation.Date, arg3: str, arg4: str, arg5: bool, arg6: memilio.simulation.osecir.ModelGraph, arg7: List[float], arg8: float, arg9: float, arg10: int, arg11: bool) -> None"""
def set_params_distributions_normal(model: Model, t0: float, tmax: float, dev_rel: float) -> None:
    """set_params_distributions_normal(model: memilio.simulation.osecir.Model, t0: float, tmax: float, dev_rel: float) -> None"""
def simulate(t0: float, tmax: float, dt: float, model: Model, integrator: memilio.simulation.IntegratorCore = ...) -> memilio.simulation.TimeSeries:
    """simulate(t0: float, tmax: float, dt: float, model: memilio.simulation.osecir.Model, integrator: memilio.simulation.IntegratorCore = None) -> memilio.simulation.TimeSeries

    Simulates an ODE SECIHURD model from t0 to tmax.
    """
def simulate_flows(t0: float, tmax: float, dt: float, model: Model, integrator: memilio.simulation.IntegratorCore = ...) -> list[memilio.simulation.TimeSeries]:
    """simulate_flows(t0: float, tmax: float, dt: float, model: memilio.simulation.osecir.Model, integrator: memilio.simulation.IntegratorCore = None) -> List[memilio.simulation.TimeSeries]

    Simulates an ODE SECIHURD model with flows from t0 to tmax.
    """
def write_graph(arg0: ModelGraph, arg1: str) -> None:
    """write_graph(arg0: memilio.simulation.osecir.ModelGraph, arg1: str) -> None"""

