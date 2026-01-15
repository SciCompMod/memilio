from memilio.simulation import (
    osecir as osecir,
    abm as abm,
    osecirvvs as osecirvvs,
    oseir as oseir,
    osir as osir,
)
import flags
import numpy
import typing
from typing import ClassVar, Iterator, overload

__version__: str

class AgeGroup(Index_AgeGroup):
    def __init__(self, arg0: int) -> None:
        """__init__(self: memilio.simulation.AgeGroup, arg0: int) -> None"""

class AgeGroupArray:
    @overload
    def __init__(self, arg0: Index_AgeGroup, arg1: UncertainValue) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.UncertainValue) -> None

        2. __init__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup) -> None
        """
    @overload
    def __init__(self, arg0: Index_AgeGroup) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.UncertainValue) -> None

        2. __init__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup) -> None
        """
    def get_flat_index(self, arg0: Index_AgeGroup) -> int:
        """get_flat_index(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup) -> int"""
    def numel(self) -> int:
        """numel(self: memilio.simulation.AgeGroupArray) -> int"""
    def resize(self, arg0: Index_AgeGroup) -> None:
        """resize(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup) -> None"""
    def size(self) -> Index_AgeGroup:
        """size(self: memilio.simulation.AgeGroupArray) -> memilio.simulation.Index_AgeGroup"""
    def size_AgeGroup(self) -> Index_AgeGroup:
        """size_AgeGroup(self: memilio.simulation.AgeGroupArray) -> memilio.simulation.Index_AgeGroup"""
    @overload
    def __getitem__(self, arg0: Index_AgeGroup) -> UncertainValue:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup) -> memilio.simulation.UncertainValue

        2. __getitem__(self: memilio.simulation.AgeGroupArray, arg0: Tuple[memilio.simulation.Index_AgeGroup]) -> memilio.simulation.UncertainValue
        """
    @overload
    def __getitem__(self, arg0: tuple[Index_AgeGroup]) -> UncertainValue:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup) -> memilio.simulation.UncertainValue

        2. __getitem__(self: memilio.simulation.AgeGroupArray, arg0: Tuple[memilio.simulation.Index_AgeGroup]) -> memilio.simulation.UncertainValue
        """
    def __iter__(self) -> Iterator:
        """__iter__(self: memilio.simulation.AgeGroupArray) -> Iterator"""
    @overload
    def __setitem__(self, arg0: Index_AgeGroup, arg1: UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: Tuple[memilio.simulation.Index_AgeGroup], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: tuple[Index_AgeGroup], arg1: UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: Tuple[memilio.simulation.Index_AgeGroup], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: Tuple[memilio.simulation.Index_AgeGroup], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: float) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: Tuple[memilio.simulation.Index_AgeGroup], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.AgeGroupArray, arg0: object, arg1: float) -> None
        """

class AgeGroupSimulationDayArray:
    @overload
    def __init__(self, arg0: MultiIndex_AgeGroupSimulationDayArray, arg1: float) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg1: float) -> None

        2. __init__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray) -> None
        """
    @overload
    def __init__(self, arg0: MultiIndex_AgeGroupSimulationDayArray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg1: float) -> None

        2. __init__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray) -> None
        """
    def get_flat_index(self, arg0: MultiIndex_AgeGroupSimulationDayArray) -> int:
        """get_flat_index(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray) -> int"""
    def numel(self) -> int:
        """numel(self: memilio.simulation.AgeGroupSimulationDayArray) -> int"""
    def resize(self, arg0: MultiIndex_AgeGroupSimulationDayArray) -> None:
        """resize(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray) -> None"""
    def resize_AgeGroup(self, arg0: Index_AgeGroup) -> None:
        """resize_AgeGroup(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.Index_AgeGroup) -> None"""
    def resize_SimulationDay(self, arg0: Index_SimulationDay) -> None:
        """resize_SimulationDay(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.Index_SimulationDay) -> None"""
    def size(self) -> tuple[Index_AgeGroup, Index_SimulationDay]:
        """size(self: memilio.simulation.AgeGroupSimulationDayArray) -> Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.Index_SimulationDay]"""
    def size_AgeGroup(self) -> Index_AgeGroup:
        """size_AgeGroup(self: memilio.simulation.AgeGroupSimulationDayArray) -> memilio.simulation.Index_AgeGroup"""
    def size_SimulationDay(self) -> Index_SimulationDay:
        """size_SimulationDay(self: memilio.simulation.AgeGroupSimulationDayArray) -> memilio.simulation.Index_SimulationDay"""
    @overload
    def __getitem__(self, arg0: MultiIndex_AgeGroupSimulationDayArray) -> float:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray) -> float

        2. __getitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.Index_SimulationDay]) -> float
        """
    @overload
    def __getitem__(self, arg0: tuple[Index_AgeGroup, Index_SimulationDay]) -> float:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray) -> float

        2. __getitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.Index_SimulationDay]) -> float
        """
    def __iter__(self) -> Iterator:
        """__iter__(self: memilio.simulation.AgeGroupSimulationDayArray) -> Iterator"""
    @overload
    def __setitem__(self, arg0: MultiIndex_AgeGroupSimulationDayArray, arg1: float) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg1: float) -> None

        2. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.Index_SimulationDay], arg1: float) -> None

        3. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: object, arg1: float) -> None

        4. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: tuple[Index_AgeGroup, Index_SimulationDay], arg1: float) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg1: float) -> None

        2. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.Index_SimulationDay], arg1: float) -> None

        3. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: object, arg1: float) -> None

        4. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: float) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg1: float) -> None

        2. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.Index_SimulationDay], arg1: float) -> None

        3. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: object, arg1: float) -> None

        4. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: float) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg1: float) -> None

        2. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: Tuple[memilio.simulation.Index_AgeGroup, memilio.simulation.Index_SimulationDay], arg1: float) -> None

        3. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: object, arg1: float) -> None

        4. __setitem__(self: memilio.simulation.AgeGroupSimulationDayArray, arg0: object, arg1: float) -> None
        """

class ContactMatrix:
    baseline: numpy.ndarray[numpy.float64[m, n]]
    minimum: numpy.ndarray[numpy.float64[m, n]]
    @overload
    def __init__(self, baseline: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous], minimum: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous]) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.ContactMatrix, baseline: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous], minimum: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous]) -> None

        2. __init__(self: memilio.simulation.ContactMatrix, baseline: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous]) -> None

        3. __init__(self: memilio.simulation.ContactMatrix, size: int) -> None
        """
    @overload
    def __init__(self, baseline: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous]) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.ContactMatrix, baseline: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous], minimum: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous]) -> None

        2. __init__(self: memilio.simulation.ContactMatrix, baseline: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous]) -> None

        3. __init__(self: memilio.simulation.ContactMatrix, size: int) -> None
        """
    @overload
    def __init__(self, size: int) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.ContactMatrix, baseline: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous], minimum: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous]) -> None

        2. __init__(self: memilio.simulation.ContactMatrix, baseline: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous]) -> None

        3. __init__(self: memilio.simulation.ContactMatrix, size: int) -> None
        """
    def add_damping(self, arg0: Damping) -> None:
        """add_damping(self: memilio.simulation.ContactMatrix, arg0: memilio.simulation.Damping) -> None"""
    def get_dampings(self) -> list[Damping]:
        """get_dampings(self: memilio.simulation.ContactMatrix) -> List[memilio.simulation.Damping]"""
    def get_matrix_at(self, arg0: float) -> numpy.ndarray[numpy.float64[m, n]]:
        """get_matrix_at(self: memilio.simulation.ContactMatrix, arg0: float) -> numpy.ndarray[numpy.float64[m, n]]"""
    @property
    def num_groups(self) -> int: ...
    @property
    def shape(self) -> tuple: ...

class ContactMatrixGroup:
    def __init__(self, num_matrices: int, size: int) -> None:
        """__init__(self: memilio.simulation.ContactMatrixGroup, num_matrices: int, size: int) -> None"""
    def add_damping(self, arg0: Damping) -> None:
        """add_damping(self: memilio.simulation.ContactMatrixGroup, arg0: memilio.simulation.Damping) -> None"""
    def get_matrix_at(self, arg0: float) -> numpy.ndarray[numpy.float64[m, n]]:
        """get_matrix_at(self: memilio.simulation.ContactMatrixGroup, arg0: float) -> numpy.ndarray[numpy.float64[m, n]]"""
    def __getitem__(self, arg0: int) -> ContactMatrix:
        """__getitem__(self: memilio.simulation.ContactMatrixGroup, arg0: int) -> memilio.simulation.ContactMatrix"""
    def __iter__(self) -> typing.Iterator[ContactMatrix]:
        """def __iter__(self) -> typing.Iterator[ContactMatrix]"""
    def __setitem__(self, arg0: int, arg1: ContactMatrix) -> None:
        """__setitem__(self: memilio.simulation.ContactMatrixGroup, arg0: int, arg1: memilio.simulation.ContactMatrix) -> None"""
    @property
    def num_groups(self) -> int: ...
    @property
    def num_matrices(self) -> int: ...
    @property
    def shape(self) -> tuple: ...

class Damping:
    coeffs: numpy.ndarray[numpy.float64[m, n]]
    level: int
    time: float
    type: int
    @overload
    def __init__(self, size: int) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.Damping, size: int) -> None

        2. __init__(self: memilio.simulation.Damping, coeffs: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous], t: float, level: int = 0, type: int = 0) -> None
        """
    @overload
    def __init__(self, coeffs: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous], t: float, level: int = ..., type: int = ...) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.Damping, size: int) -> None

        2. __init__(self: memilio.simulation.Damping, coeffs: numpy.ndarray[numpy.float64[m, n], flags.f_contiguous], t: float, level: int = 0, type: int = 0) -> None
        """
    @property
    def shape(self) -> tuple: ...

class DampingSampling:
    group_weights: numpy.ndarray[numpy.float64[m, 1]]
    level: int
    matrix_indices: list[int]
    time: float
    type: int
    value: UncertainValue
    def __init__(self, value: UncertainValue, level: int, type: int, time: float, matrix_indices: list[int], group_weights: numpy.ndarray[numpy.float64[m, 1]]) -> None:
        """__init__(self: memilio.simulation.DampingSampling, value: memilio.simulation.UncertainValue, level: int, type: int, time: float, matrix_indices: List[int], group_weights: numpy.ndarray[numpy.float64[m, 1]]) -> None"""

class Dampings:
    def __init__(self, size: int) -> None:
        """__init__(self: memilio.simulation.Dampings, size: int) -> None"""
    def add(self, arg0: Damping) -> None:
        """add(self: memilio.simulation.Dampings, arg0: memilio.simulation.Damping) -> None"""
    def get_matrix_at(self, arg0: float) -> numpy.ndarray[numpy.float64[m, n]]:
        """get_matrix_at(self: memilio.simulation.Dampings, arg0: float) -> numpy.ndarray[numpy.float64[m, n]]"""
    @property
    def shape(self) -> tuple: ...

class Date:
    day: int
    month: int
    year: int
    def __init__(self, year: int, month: int, day: int) -> None:
        """__init__(self: memilio.simulation.Date, year: int, month: int, day: int) -> None"""
    def __add__(self, arg0: int) -> Date:
        """__add__(self: memilio.simulation.Date, arg0: int) -> memilio.simulation.Date"""
    def __eq__(self, arg0: Date) -> bool:
        """__eq__(self: memilio.simulation.Date, arg0: memilio.simulation.Date) -> bool"""
    def __ge__(self, arg0: Date) -> bool:
        """__ge__(self: memilio.simulation.Date, arg0: memilio.simulation.Date) -> bool"""
    def __gt__(self, arg0: Date) -> bool:
        """__gt__(self: memilio.simulation.Date, arg0: memilio.simulation.Date) -> bool"""
    def __iadd__(self, arg0: int) -> Date:
        """__iadd__(self: memilio.simulation.Date, arg0: int) -> memilio.simulation.Date"""
    def __isub__(self, arg0: int) -> Date:
        """__isub__(self: memilio.simulation.Date, arg0: int) -> memilio.simulation.Date"""
    def __le__(self, arg0: Date) -> bool:
        """__le__(self: memilio.simulation.Date, arg0: memilio.simulation.Date) -> bool"""
    def __lt__(self, arg0: Date) -> bool:
        """__lt__(self: memilio.simulation.Date, arg0: memilio.simulation.Date) -> bool"""
    def __ne__(self, arg0: Date) -> bool:
        """__ne__(self: memilio.simulation.Date, arg0: memilio.simulation.Date) -> bool"""
    @overload
    def __sub__(self, arg0: int) -> Date:
        """__sub__(*args, **kwargs)
        Overloaded function.

        1. __sub__(self: memilio.simulation.Date, arg0: int) -> memilio.simulation.Date

        2. __sub__(self: memilio.simulation.Date, arg0: memilio.simulation.Date) -> int
        """
    @overload
    def __sub__(self, arg0: Date) -> int:
        """__sub__(*args, **kwargs)
        Overloaded function.

        1. __sub__(self: memilio.simulation.Date, arg0: int) -> memilio.simulation.Date

        2. __sub__(self: memilio.simulation.Date, arg0: memilio.simulation.Date) -> int
        """
    @property
    def day_in_year(self) -> int: ...

class DynamicNPIs:
    base_value: float
    duration: float
    interval: float
    def __init__(self) -> None:
        """__init__(self: memilio.simulation.DynamicNPIs) -> None"""
    def set_threshold(self, arg0: float, arg1: list[DampingSampling]) -> None:
        """set_threshold(self: memilio.simulation.DynamicNPIs, arg0: float, arg1: List[memilio.simulation.DampingSampling]) -> None"""
    @property
    def threshold(self) -> _ThresholdRange: ...

class EulerIntegratorCore(IntegratorCore):
    def __init__(self) -> None:
        """__init__(self: memilio.simulation.EulerIntegratorCore) -> None"""
    def step(self, f: function, yt: numpy.ndarray[numpy.float64[m, 1]], t: float, dt: float, ytp1: numpy.ndarray[numpy.float64[m, 1], flags.writeable]) -> bool:
        """step(self: memilio.simulation.EulerIntegratorCore, f: function, yt: numpy.ndarray[numpy.float64[m, 1]], t: float, dt: float, ytp1: numpy.ndarray[numpy.float64[m, 1], flags.writeable]) -> bool"""

class Index_AgeGroup:
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.Index_AgeGroup, value: int) -> None"""
    def __eq__(self, arg0: Index_AgeGroup) -> bool:
        """__eq__(self: memilio.simulation.Index_AgeGroup, arg0: memilio.simulation.Index_AgeGroup) -> bool"""
    def __ne__(self, arg0: Index_AgeGroup) -> bool:
        """__ne__(self: memilio.simulation.Index_AgeGroup, arg0: memilio.simulation.Index_AgeGroup) -> bool"""

class Index_SimulationDay:
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.Index_SimulationDay, value: int) -> None"""
    def __eq__(self, arg0: Index_SimulationDay) -> bool:
        """__eq__(self: memilio.simulation.Index_SimulationDay, arg0: memilio.simulation.Index_SimulationDay) -> bool"""
    def __ne__(self, arg0: Index_SimulationDay) -> bool:
        """__ne__(self: memilio.simulation.Index_SimulationDay, arg0: memilio.simulation.Index_SimulationDay) -> bool"""

class IntegratorCore:
    dt_max: float
    dt_min: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class LogLevel:
    __members__: ClassVar[dict] = ...  # read-only
    Critical: ClassVar[LogLevel] = ...
    Debug: ClassVar[LogLevel] = ...
    Error: ClassVar[LogLevel] = ...
    Info: ClassVar[LogLevel] = ...
    Off: ClassVar[LogLevel] = ...
    Trace: ClassVar[LogLevel] = ...
    Warning: ClassVar[LogLevel] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.LogLevel, value: int) -> None"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object) -> int"""
    def __index__(self) -> int:
        """__index__(self: memilio.simulation.LogLevel) -> int"""
    def __int__(self) -> int:
        """__int__(self: memilio.simulation.LogLevel) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object) -> bool"""
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Mobility:
    @overload
    def __init__(self, coeffs: numpy.ndarray[numpy.float64[m, 1]]) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.Mobility, coeffs: numpy.ndarray[numpy.float64[m, 1]]) -> None

        2. __init__(self: memilio.simulation.Mobility, params: memilio.simulation.MobilityParameters) -> None
        """
    @overload
    def __init__(self, params: MobilityParameters) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.Mobility, coeffs: numpy.ndarray[numpy.float64[m, 1]]) -> None

        2. __init__(self: memilio.simulation.Mobility, params: memilio.simulation.MobilityParameters) -> None
        """
    @property
    def parameters(self) -> MobilityParameters: ...

class MobilityCoefficientGroup:
    def __init__(self, num_matrices: int, size: int) -> None:
        """__init__(self: memilio.simulation.MobilityCoefficientGroup, num_matrices: int, size: int) -> None"""
    def add_damping(self, arg0: MobilityDamping) -> None:
        """add_damping(self: memilio.simulation.MobilityCoefficientGroup, arg0: memilio.simulation.MobilityDamping) -> None"""
    def get_matrix_at(self, arg0: float) -> numpy.ndarray[numpy.float64[m, n]]:
        """get_matrix_at(self: memilio.simulation.MobilityCoefficientGroup, arg0: float) -> numpy.ndarray[numpy.float64[m, n]]"""
    def __getitem__(self, arg0: int) -> MobilityCoefficients:
        """__getitem__(self: memilio.simulation.MobilityCoefficientGroup, arg0: int) -> memilio.simulation.MobilityCoefficients"""
    def __iter__(self) -> typing.Iterator[MobilityCoefficients]:
        """def __iter__(self) -> typing.Iterator[MobilityCoefficients]"""
    def __setitem__(self, arg0: int, arg1: MobilityCoefficients) -> None:
        """__setitem__(self: memilio.simulation.MobilityCoefficientGroup, arg0: int, arg1: memilio.simulation.MobilityCoefficients) -> None"""
    @property
    def num_matrices(self) -> int: ...
    @property
    def shape(self) -> tuple: ...

class MobilityCoefficients:
    baseline: numpy.ndarray[numpy.float64[m, 1]]
    minimum: numpy.ndarray[numpy.float64[m, 1]]
    @overload
    def __init__(self, baseline: numpy.ndarray[numpy.float64[m, 1]], minimum: numpy.ndarray[numpy.float64[m, 1]]) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.MobilityCoefficients, baseline: numpy.ndarray[numpy.float64[m, 1]], minimum: numpy.ndarray[numpy.float64[m, 1]]) -> None

        2. __init__(self: memilio.simulation.MobilityCoefficients, baseline: numpy.ndarray[numpy.float64[m, 1]]) -> None

        3. __init__(self: memilio.simulation.MobilityCoefficients, size: int) -> None
        """
    @overload
    def __init__(self, baseline: numpy.ndarray[numpy.float64[m, 1]]) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.MobilityCoefficients, baseline: numpy.ndarray[numpy.float64[m, 1]], minimum: numpy.ndarray[numpy.float64[m, 1]]) -> None

        2. __init__(self: memilio.simulation.MobilityCoefficients, baseline: numpy.ndarray[numpy.float64[m, 1]]) -> None

        3. __init__(self: memilio.simulation.MobilityCoefficients, size: int) -> None
        """
    @overload
    def __init__(self, size: int) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.MobilityCoefficients, baseline: numpy.ndarray[numpy.float64[m, 1]], minimum: numpy.ndarray[numpy.float64[m, 1]]) -> None

        2. __init__(self: memilio.simulation.MobilityCoefficients, baseline: numpy.ndarray[numpy.float64[m, 1]]) -> None

        3. __init__(self: memilio.simulation.MobilityCoefficients, size: int) -> None
        """
    def add_damping(self, arg0: MobilityDamping) -> None:
        """add_damping(self: memilio.simulation.MobilityCoefficients, arg0: memilio.simulation.MobilityDamping) -> None"""
    def get_dampings(self) -> list[MobilityDamping]:
        """get_dampings(self: memilio.simulation.MobilityCoefficients) -> List[memilio.simulation.MobilityDamping]"""
    def get_matrix_at(self, arg0: float) -> numpy.ndarray[numpy.float64[m, 1]]:
        """get_matrix_at(self: memilio.simulation.MobilityCoefficients, arg0: float) -> numpy.ndarray[numpy.float64[m, 1]]"""
    @property
    def shape(self) -> tuple: ...

class MobilityDamping:
    coeffs: numpy.ndarray[numpy.float64[m, 1]]
    level: int
    time: float
    type: int
    @overload
    def __init__(self, size: int) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.MobilityDamping, size: int) -> None

        2. __init__(self: memilio.simulation.MobilityDamping, coeffs: numpy.ndarray[numpy.float64[m, 1]], t: float, level: int = 0, type: int = 0) -> None
        """
    @overload
    def __init__(self, coeffs: numpy.ndarray[numpy.float64[m, 1]], t: float, level: int = ..., type: int = ...) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.MobilityDamping, size: int) -> None

        2. __init__(self: memilio.simulation.MobilityDamping, coeffs: numpy.ndarray[numpy.float64[m, 1]], t: float, level: int = 0, type: int = 0) -> None
        """
    @property
    def shape(self) -> tuple: ...

class MobilityDampings:
    def __init__(self, size: int) -> None:
        """__init__(self: memilio.simulation.MobilityDampings, size: int) -> None"""
    def add(self, arg0: MobilityDamping) -> None:
        """add(self: memilio.simulation.MobilityDampings, arg0: memilio.simulation.MobilityDamping) -> None"""
    def get_matrix_at(self, arg0: float) -> numpy.ndarray[numpy.float64[m, 1]]:
        """get_matrix_at(self: memilio.simulation.MobilityDampings, arg0: float) -> numpy.ndarray[numpy.float64[m, 1]]"""
    @property
    def shape(self) -> tuple: ...

class MobilityEdge:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    @property
    def end_node_idx(self) -> int: ...
    @property
    def property(self) -> Mobility: ...
    @property
    def start_node_idx(self) -> int: ...

class MobilityParameterEdge:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    @property
    def end_node_idx(self) -> int: ...
    @property
    def property(self) -> Mobility: ...
    @property
    def start_node_idx(self) -> int: ...

class MobilityParameters:
    coefficients: MobilityCoefficientGroup
    @overload
    def __init__(self, coeffs: numpy.ndarray[numpy.float64[m, 1]]) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.MobilityParameters, coeffs: numpy.ndarray[numpy.float64[m, 1]]) -> None

        2. __init__(self: memilio.simulation.MobilityParameters, coeffs: memilio.simulation.MobilityCoefficientGroup) -> None
        """
    @overload
    def __init__(self, coeffs: MobilityCoefficientGroup) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.MobilityParameters, coeffs: numpy.ndarray[numpy.float64[m, 1]]) -> None

        2. __init__(self: memilio.simulation.MobilityParameters, coeffs: memilio.simulation.MobilityCoefficientGroup) -> None
        """

class MultiIndex_AgeGroupSimulationDayArray:
    @overload
    def __init__(self, arg0: Index_AgeGroup, arg1: Index_SimulationDay) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.Index_SimulationDay) -> None

        2. __init__(self: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg0: tuple) -> None
        """
    @overload
    def __init__(self, arg0: tuple) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg0: memilio.simulation.Index_AgeGroup, arg1: memilio.simulation.Index_SimulationDay) -> None

        2. __init__(self: memilio.simulation.MultiIndex_AgeGroupSimulationDayArray, arg0: tuple) -> None
        """

class ParameterDistribution:
    lower_bound: float
    upper_bound: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def add_predefined_sample(self, arg0: float) -> None:
        """add_predefined_sample(self: memilio.simulation.ParameterDistribution, arg0: float) -> None"""
    def get_sample(self) -> float:
        """get_sample(self: memilio.simulation.ParameterDistribution) -> float"""
    def remove_predefined_samples(self) -> None:
        """remove_predefined_samples(self: memilio.simulation.ParameterDistribution) -> None"""

class ParameterDistributionNormal(ParameterDistribution):
    mean: float
    standard_dev: float
    @overload
    def __init__(self, lb: float, ub: float, mean: float, std_dev: float) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.ParameterDistributionNormal, lb: float, ub: float, mean: float, std_dev: float) -> None

        2. __init__(self: memilio.simulation.ParameterDistributionNormal, lb: float, ub: float, mean: float) -> None
        """
    @overload
    def __init__(self, lb: float, ub: float, mean: float) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.ParameterDistributionNormal, lb: float, ub: float, mean: float, std_dev: float) -> None

        2. __init__(self: memilio.simulation.ParameterDistributionNormal, lb: float, ub: float, mean: float) -> None
        """

class ParameterDistributionUniform(ParameterDistribution):
    @overload
    def __init__(self) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.ParameterDistributionUniform) -> None

        2. __init__(self: memilio.simulation.ParameterDistributionUniform, lb: float, ub: float) -> None
        """
    @overload
    def __init__(self, lb: float, ub: float) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.ParameterDistributionUniform) -> None

        2. __init__(self: memilio.simulation.ParameterDistributionUniform, lb: float, ub: float) -> None
        """

class RKIntegratorCore(IntegratorCore):
    @overload
    def __init__(self) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.RKIntegratorCore) -> None

        2. __init__(self: memilio.simulation.RKIntegratorCore, abs_tol: float = 1e-10, rel_tol: float = 1e-05, dt_min: float = 2.2250738585072014e-308, dt_max: float = 1.7976931348623157e+308) -> None
        """
    @overload
    def __init__(self, abs_tol: float = ..., rel_tol: float = ..., dt_min: float = ..., dt_max: float = ...) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.RKIntegratorCore) -> None

        2. __init__(self: memilio.simulation.RKIntegratorCore, abs_tol: float = 1e-10, rel_tol: float = 1e-05, dt_min: float = 2.2250738585072014e-308, dt_max: float = 1.7976931348623157e+308) -> None
        """
    def set_abs_tolerance(self, tol: float) -> None:
        """set_abs_tolerance(self: memilio.simulation.RKIntegratorCore, tol: float) -> None"""
    def set_rel_tolerance(self, tol: float) -> None:
        """set_rel_tolerance(self: memilio.simulation.RKIntegratorCore, tol: float) -> None"""

class RungeKuttaCashKarp54IntegratorCore(IntegratorCore):
    @overload
    def __init__(self) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.RungeKuttaCashKarp54IntegratorCore) -> None

        2. __init__(self: memilio.simulation.RungeKuttaCashKarp54IntegratorCore, abs_tol: float, rel_tol: float, dt_min: float, dt_max: float) -> None
        """
    @overload
    def __init__(self, abs_tol: float, rel_tol: float, dt_min: float, dt_max: float) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.RungeKuttaCashKarp54IntegratorCore) -> None

        2. __init__(self: memilio.simulation.RungeKuttaCashKarp54IntegratorCore, abs_tol: float, rel_tol: float, dt_min: float, dt_max: float) -> None
        """
    def set_abs_tolerance(self, tol: float) -> None:
        """set_abs_tolerance(self: memilio.simulation.RungeKuttaCashKarp54IntegratorCore, tol: float) -> None"""
    def set_rel_tolerance(self, tol: float) -> None:
        """set_rel_tolerance(self: memilio.simulation.RungeKuttaCashKarp54IntegratorCore, tol: float) -> None"""

class SimulationDay(Index_SimulationDay):
    def __init__(self, arg0: int) -> None:
        """__init__(self: memilio.simulation.SimulationDay, arg0: int) -> None"""

class TimeSeries:
    def __init__(self, num_elements: int) -> None:
        """__init__(self: memilio.simulation.TimeSeries, num_elements: int) -> None"""
    @overload
    def add_time_point(self) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """add_time_point(*args, **kwargs)
        Overloaded function.

        1. add_time_point(self: memilio.simulation.TimeSeries) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]

        2. add_time_point(self: memilio.simulation.TimeSeries, arg0: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]

        3. add_time_point(self: memilio.simulation.TimeSeries, arg0: float, arg1: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]
        """
    @overload
    def add_time_point(self, arg0: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """add_time_point(*args, **kwargs)
        Overloaded function.

        1. add_time_point(self: memilio.simulation.TimeSeries) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]

        2. add_time_point(self: memilio.simulation.TimeSeries, arg0: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]

        3. add_time_point(self: memilio.simulation.TimeSeries, arg0: float, arg1: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]
        """
    @overload
    def add_time_point(self, arg0: float, arg1: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """add_time_point(*args, **kwargs)
        Overloaded function.

        1. add_time_point(self: memilio.simulation.TimeSeries) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]

        2. add_time_point(self: memilio.simulation.TimeSeries, arg0: float) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]

        3. add_time_point(self: memilio.simulation.TimeSeries, arg0: float, arg1: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]
        """
    def as_ndarray(self) -> numpy.ndarray[numpy.float64[m, n], flags.writeable, flags.f_contiguous]:
        """as_ndarray(self: memilio.simulation.TimeSeries) -> numpy.ndarray[numpy.float64[m, n], flags.writeable, flags.f_contiguous]"""
    def get_last_time(self) -> float:
        """get_last_time(self: memilio.simulation.TimeSeries) -> float"""
    def get_last_value(self) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """get_last_value(self: memilio.simulation.TimeSeries) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]"""
    def get_num_elements(self) -> int:
        """get_num_elements(self: memilio.simulation.TimeSeries) -> int"""
    def get_num_time_points(self) -> int:
        """get_num_time_points(self: memilio.simulation.TimeSeries) -> int"""
    def get_time(self, index: int) -> float:
        """get_time(self: memilio.simulation.TimeSeries, index: int) -> float"""
    def get_value(self, index: int) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """get_value(self: memilio.simulation.TimeSeries, index: int) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]"""
    def print_table(self, arg0: list[str], arg1: int, arg2: int) -> str:
        """print_table(self: memilio.simulation.TimeSeries, arg0: List[str], arg1: int, arg2: int) -> str"""
    def __getitem__(self, index: int) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]:
        """__getitem__(self: memilio.simulation.TimeSeries, index: int) -> numpy.ndarray[numpy.float64[m, 1], flags.writeable]"""
    def __iter__(self) -> typing.Iterator[numpy.ndarray[numpy.float64[m, 1], flags.writeable]]:
        """def __iter__(self) -> typing.Iterator[numpy.ndarray[numpy.float64[m, 1], flags.writeable]]"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.TimeSeries) -> int"""
    def __setitem__(self, index: int, v: numpy.ndarray[numpy.float64[m, 1]]) -> None:
        """__setitem__(self: memilio.simulation.TimeSeries, index: int, v: numpy.ndarray[numpy.float64[m, 1]]) -> None"""

class UncertainContactMatrix:
    cont_freq_mat: ContactMatrixGroup
    dampings: list[DampingSampling]
    school_holiday_damping: DampingSampling
    school_holidays: list[tuple[float, float]]
    @overload
    def __init__(self) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.UncertainContactMatrix) -> None

        2. __init__(self: memilio.simulation.UncertainContactMatrix, arg0: memilio.simulation.ContactMatrixGroup) -> None
        """
    @overload
    def __init__(self, arg0: ContactMatrixGroup) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.UncertainContactMatrix) -> None

        2. __init__(self: memilio.simulation.UncertainContactMatrix, arg0: memilio.simulation.ContactMatrixGroup) -> None
        """

class UncertainValue:
    value: float
    def __init__(self, value: float = ...) -> None:
        """__init__(self: memilio.simulation.UncertainValue, value: float = 0.0) -> None"""
    def draw_sample(self) -> float:
        """draw_sample(self: memilio.simulation.UncertainValue) -> float"""
    @overload
    def get_distribution(self) -> ParameterDistribution:
        """get_distribution(*args, **kwargs)
        Overloaded function.

        1. get_distribution(self: memilio.simulation.UncertainValue) -> memilio.simulation.ParameterDistribution

        2. get_distribution(self: memilio.simulation.UncertainValue) -> memilio.simulation.ParameterDistribution
        """
    @overload
    def get_distribution(self) -> ParameterDistribution:
        """get_distribution(*args, **kwargs)
        Overloaded function.

        1. get_distribution(self: memilio.simulation.UncertainValue) -> memilio.simulation.ParameterDistribution

        2. get_distribution(self: memilio.simulation.UncertainValue) -> memilio.simulation.ParameterDistribution
        """
    def set_distribution(self, arg0: ParameterDistribution) -> None:
        """set_distribution(self: memilio.simulation.UncertainValue, arg0: memilio.simulation.ParameterDistribution) -> None"""

class _Iter_ThresholdRange:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __next__(self) -> tuple[float, list[DampingSampling]]:
        """__next__(self: memilio.simulation._Iter_ThresholdRange) -> Tuple[float, List[memilio.simulation.DampingSampling]]"""

class _ThresholdRange:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __getitem__(self, arg0: int) -> tuple[float, list[DampingSampling]]:
        """__getitem__(self: memilio.simulation._ThresholdRange, arg0: int) -> Tuple[float, List[memilio.simulation.DampingSampling]]"""
    def __iter__(self) -> _Iter_ThresholdRange:
        """__iter__(self: memilio.simulation._ThresholdRange) -> memilio.simulation._Iter_ThresholdRange"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation._ThresholdRange) -> int"""

def get_holidays_de(state_id: int, start_date: Date = ..., end_date: Date = ...) -> list[tuple[Date, Date]]:
    """get_holidays_de(state_id: int, start_date: memilio.simulation.Date = <memilio.simulation.Date object at 0x7e0e4bbf5a30>, end_date: memilio.simulation.Date = <memilio.simulation.Date object at 0x7e0e4bbf5670>) -> List[Tuple[memilio.simulation.Date, memilio.simulation.Date]]"""
def get_node_ids(arg0: str, arg1: bool) -> list[int]:
    """get_node_ids(arg0: str, arg1: bool) -> List[int]"""
def get_state_id_de(county_id: int) -> int:
    """get_state_id_de(county_id: int) -> int"""
def read_mobility_plain(arg0: str) -> numpy.ndarray[numpy.float64[m, n]]:
    """read_mobility_plain(arg0: str) -> numpy.ndarray[numpy.float64[m, n]]"""
def seed_random_number_generator() -> None:
    """seed_random_number_generator() -> None"""
def set_log_level(arg0: LogLevel) -> None:
    """set_log_level(arg0: memilio.simulation.LogLevel) -> None"""


def __getattr__(attr):
    """ The __getattr__ function is used here to implement lazy loading for the submodules within the memilio.simulation package. 
    Submodules are only imported when they are first accessed, which can save memory and reduce startup time, if not all submodules 
    are needed for every execution.
    """
