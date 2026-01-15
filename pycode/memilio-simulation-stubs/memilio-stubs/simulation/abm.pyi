import memilio.simulation
from _typeshed import Incomplete
from typing import ClassVar, Iterator, overload

__version__: str

class Index_N3mio3abm12VirusVariantE:
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, value: int) -> None"""
    def __eq__(self, arg0: Index_N3mio3abm12VirusVariantE) -> bool:
        """__eq__(self: memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, arg0: memilio.simulation.abm.Index_N3mio3abm12VirusVariantE) -> bool"""
    def __ne__(self, arg0: Index_N3mio3abm12VirusVariantE) -> bool:
        """__ne__(self: memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, arg0: memilio.simulation.abm.Index_N3mio3abm12VirusVariantE) -> bool"""

class Index_N3mio3abm8TestTypeE:
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.abm.Index_N3mio3abm8TestTypeE, value: int) -> None"""
    def __eq__(self, arg0: Index_N3mio3abm8TestTypeE) -> bool:
        """__eq__(self: memilio.simulation.abm.Index_N3mio3abm8TestTypeE, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE) -> bool"""
    def __ne__(self, arg0: Index_N3mio3abm8TestTypeE) -> bool:
        """__ne__(self: memilio.simulation.abm.Index_N3mio3abm8TestTypeE, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE) -> bool"""

class InfectionState:
    __members__: ClassVar[dict] = ...  # read-only
    Dead: ClassVar[InfectionState] = ...
    Exposed: ClassVar[InfectionState] = ...
    InfectedCritical: ClassVar[InfectionState] = ...
    InfectedNoSymptoms: ClassVar[InfectionState] = ...
    InfectedSevere: ClassVar[InfectionState] = ...
    InfectedSymptoms: ClassVar[InfectionState] = ...
    Recovered: ClassVar[InfectionState] = ...
    Susceptible: ClassVar[InfectionState] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.abm.InfectionState, value: int) -> None"""
    @staticmethod
    def values() -> _InfectionStateValues:
        """values() -> memilio.simulation.abm._InfectionStateValues"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object) -> int"""
    def __index__(self) -> int:
        """__index__(self: memilio.simulation.abm.InfectionState) -> int"""
    def __int__(self) -> int:
        """__int__(self: memilio.simulation.abm.InfectionState) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object) -> bool"""
    def __next__(self) -> InfectionState:
        """__next__(self: memilio.simulation.abm.InfectionState) -> memilio.simulation.abm.InfectionState"""
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class LocalInfectionParameters:
    ContactRates: Incomplete
    MaximumContacts: float
    UseLocationCapacityForTransmissions: bool
    def __init__(self, arg0: int) -> None:
        """__init__(self: memilio.simulation.abm.LocalInfectionParameters, arg0: int) -> None"""

class Location:
    infection_parameters: LocalInfectionParameters
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    @property
    def id(self) -> LocationId: ...
    @property
    def type(self) -> LocationType: ...

class LocationId:
    def __init__(self, id: int) -> None:
        """__init__(self: memilio.simulation.abm.LocationId, id: int) -> None"""
    def index(self) -> int:
        """index(self: memilio.simulation.abm.LocationId) -> int"""

class LocationType:
    __members__: ClassVar[dict] = ...  # read-only
    BasicsShop: ClassVar[LocationType] = ...
    Car: ClassVar[LocationType] = ...
    Home: ClassVar[LocationType] = ...
    Hospital: ClassVar[LocationType] = ...
    ICU: ClassVar[LocationType] = ...
    PublicTransport: ClassVar[LocationType] = ...
    School: ClassVar[LocationType] = ...
    SocialEvent: ClassVar[LocationType] = ...
    TransportWithoutContact: ClassVar[LocationType] = ...
    Work: ClassVar[LocationType] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.abm.LocationType, value: int) -> None"""
    @staticmethod
    def values() -> _LocationTypeValues:
        """values() -> memilio.simulation.abm._LocationTypeValues"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object) -> int"""
    def __index__(self) -> int:
        """__index__(self: memilio.simulation.abm.LocationType) -> int"""
    def __int__(self) -> int:
        """__int__(self: memilio.simulation.abm.LocationType) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object) -> bool"""
    def __next__(self) -> LocationType:
        """__next__(self: memilio.simulation.abm.LocationType) -> memilio.simulation.abm.LocationType"""
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Model:
    parameters: Parameters
    testing_strategy: TestingStrategy
    trip_list: TripList
    use_mobility_rules: bool
    def __init__(self, arg0: int) -> None:
        """__init__(self: memilio.simulation.abm.Model, arg0: int) -> None"""
    def add_location(self, location_type: LocationType, num_cells: int = ...) -> LocationId:
        """add_location(self: memilio.simulation.abm.Model, location_type: memilio.simulation.abm.LocationType, num_cells: int = 1) -> memilio.simulation.abm.LocationId"""
    def add_person(self, location_id: LocationId, age_group: memilio.simulation.AgeGroup) -> PersonId:
        """add_person(self: memilio.simulation.abm.Model, location_id: memilio.simulation.abm.LocationId, age_group: memilio.simulation.AgeGroup) -> memilio.simulation.abm.PersonId"""
    def assign_location(self, person_id: PersonId, location_id: LocationId) -> None:
        """assign_location(self: memilio.simulation.abm.Model, person_id: memilio.simulation.abm.PersonId, location_id: memilio.simulation.abm.LocationId) -> None"""
    @property
    def locations(self) -> _ModelLocationsRange: ...
    @property
    def persons(self) -> _ModelPersonsRange: ...

class MultiIndex__AgeParameterArray:
    @overload
    def __init__(self, arg0: Index_N3mio3abm12VirusVariantE, arg1: memilio.simulation.Index_AgeGroup) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg0: memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, arg1: memilio.simulation.Index_AgeGroup) -> None

        2. __init__(self: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg0: tuple) -> None
        """
    @overload
    def __init__(self, arg0: tuple) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg0: memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, arg1: memilio.simulation.Index_AgeGroup) -> None

        2. __init__(self: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg0: tuple) -> None
        """

class Parameters(ParametersBase):
    def __init__(self, arg0: int) -> None:
        """__init__(self: memilio.simulation.abm.Parameters, arg0: int) -> None"""
    def check_constraints(self) -> bool:
        """check_constraints(self: memilio.simulation.abm.Parameters) -> bool"""

class ParametersBase:
    AerosolTransmissionRates: Incomplete
    AgeGroupGotoSchool: Incomplete
    AgeGroupGotoWork: Incomplete
    BasicShoppingRate: memilio.simulation.AgeGroupArray
    CriticalToDead: _AgeParameterArray
    CriticalToRecovered: _AgeParameterArray
    DetectInfection: _AgeParameterArray
    GotoSchoolTimeMaximum: Incomplete
    GotoSchoolTimeMinimum: Incomplete
    GotoWorkTimeMaximum: Incomplete
    GotoWorkTimeMinimum: Incomplete
    HighViralLoadProtectionFactor: Incomplete
    IncubationPeriod: _AgeParameterArray
    InfectedNoSymptomsToRecovered: _AgeParameterArray
    InfectedNoSymptomsToSymptoms: _AgeParameterArray
    InfectedSymptomsToRecovered: _AgeParameterArray
    InfectedSymptomsToSevere: _AgeParameterArray
    InfectionProtectionFactor: Incomplete
    InfectivityDistributions: Incomplete
    LockdownDate: TimePoint
    MaskProtection: Incomplete
    QuarantineDuration: TimeSpan
    RecoveredToSusceptible: _AgeParameterArray
    SchoolRatio: memilio.simulation.MobilityCoefficients
    SevereToCritical: _AgeParameterArray
    SevereToDead: _AgeParameterArray
    SevereToRecovered: _AgeParameterArray
    SeverityProtectionFactor: Incomplete
    SocialEventRate: memilio.simulation.MobilityCoefficients
    TestData: _TestData
    ViralLoadDistributions: Incomplete
    WorkRatio: memilio.simulation.MobilityCoefficients
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Person:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def set_assigned_location(self, arg0: LocationType, arg1: LocationId, arg2: int) -> None:
        """set_assigned_location(self: memilio.simulation.abm.Person, arg0: memilio.simulation.abm.LocationType, arg1: memilio.simulation.abm.LocationId, arg2: int) -> None"""
    @property
    def age(self) -> memilio.simulation.AgeGroup: ...
    @property
    def is_in_quarantine(self) -> bool: ...
    @property
    def location(self) -> LocationId: ...

class PersonId:
    def __init__(self, id: int) -> None:
        """__init__(self: memilio.simulation.abm.PersonId, id: int) -> None"""
    def index(self) -> int:
        """index(self: memilio.simulation.abm.PersonId) -> int"""

class ProtectionEvent:
    time: TimePoint
    type: ProtectionType
    def __init__(self, type: ProtectionType, time: TimePoint) -> None:
        """__init__(self: memilio.simulation.abm.ProtectionEvent, type: memilio.simulation.abm.ProtectionType, time: memilio.simulation.abm.TimePoint) -> None"""

class ProtectionType:
    __members__: ClassVar[dict] = ...  # read-only
    GenericVaccine: ClassVar[ProtectionType] = ...
    NaturalInfection: ClassVar[ProtectionType] = ...
    NoProtection: ClassVar[ProtectionType] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.abm.ProtectionType, value: int) -> None"""
    @staticmethod
    def values() -> _ProtectionTypeValues:
        """values() -> memilio.simulation.abm._ProtectionTypeValues"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object) -> int"""
    def __index__(self) -> int:
        """__index__(self: memilio.simulation.abm.ProtectionType) -> int"""
    def __int__(self) -> int:
        """__int__(self: memilio.simulation.abm.ProtectionType) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object) -> bool"""
    def __next__(self) -> ProtectionType:
        """__next__(self: memilio.simulation.abm.ProtectionType) -> memilio.simulation.abm.ProtectionType"""
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class ProtectionTypeIndex:
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.abm.ProtectionTypeIndex, value: int) -> None"""
    def __eq__(self, arg0: ProtectionTypeIndex) -> bool:
        """__eq__(self: memilio.simulation.abm.ProtectionTypeIndex, arg0: memilio.simulation.abm.ProtectionTypeIndex) -> bool"""
    def __ne__(self, arg0: ProtectionTypeIndex) -> bool:
        """__ne__(self: memilio.simulation.abm.ProtectionTypeIndex, arg0: memilio.simulation.abm.ProtectionTypeIndex) -> bool"""

class Simulation:
    def __init__(self, arg0: TimePoint, arg1: int) -> None:
        """__init__(self: memilio.simulation.abm.Simulation, arg0: memilio.simulation.abm.TimePoint, arg1: int) -> None"""
    def advance(self, tmax: TimePoint) -> None:
        """advance(self: memilio.simulation.abm.Simulation, tmax: memilio.simulation.abm.TimePoint) -> None"""
    @property
    def model(self) -> Model: ...

class TestParameters:
    required_time: TimeSpan
    sensitivity: memilio.simulation.UncertainValue
    specificity: memilio.simulation.UncertainValue
    type: TestType
    def __init__(self, arg0: float, arg1: float, arg2: TimeSpan, arg3: TestType) -> None:
        """__init__(self: memilio.simulation.abm.TestParameters, arg0: float, arg1: float, arg2: memilio.simulation.abm.TimeSpan, arg3: memilio.simulation.abm.TestType) -> None"""

class TestType:
    __members__: ClassVar[dict] = ...  # read-only
    Antigen: ClassVar[TestType] = ...
    Generic: ClassVar[TestType] = ...
    PCR: ClassVar[TestType] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.abm.TestType, value: int) -> None"""
    @staticmethod
    def values() -> _TestTypeValues:
        """values() -> memilio.simulation.abm._TestTypeValues"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object) -> int"""
    def __index__(self) -> int:
        """__index__(self: memilio.simulation.abm.TestType) -> int"""
    def __int__(self) -> int:
        """__int__(self: memilio.simulation.abm.TestType) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object) -> bool"""
    def __next__(self) -> TestType:
        """__next__(self: memilio.simulation.abm.TestType) -> memilio.simulation.abm.TestType"""
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class TestingCriteria:
    def __init__(self, age_groups: list[memilio.simulation.AgeGroup], infection_states: list[InfectionState]) -> None:
        """__init__(self: memilio.simulation.abm.TestingCriteria, age_groups: List[memilio.simulation.AgeGroup], infection_states: List[memilio.simulation.abm.InfectionState]) -> None"""

class TestingScheme:
    def __init__(self, testing_criteria: TestingCriteria, testing_validity_period: TimeSpan, start_date: TimePoint, end_date: TimePoint, test_parameters: TestParameters, probability: float) -> None:
        """__init__(self: memilio.simulation.abm.TestingScheme, testing_criteria: memilio.simulation.abm.TestingCriteria, testing_validity_period: memilio.simulation.abm.TimeSpan, start_date: memilio.simulation.abm.TimePoint, end_date: memilio.simulation.abm.TimePoint, test_parameters: memilio.simulation.abm.TestParameters, probability: float) -> None"""
    @property
    def active(self) -> bool: ...

class TestingStrategy:
    def __init__(self, arg0) -> None:
        """__init__(self: memilio.simulation.abm.TestingStrategy, arg0: List[mio::abm::TestingStrategy::LocalStrategy]) -> None"""

class TimePoint:
    def __init__(self, seconds: int = ...) -> None:
        """__init__(self: memilio.simulation.abm.TimePoint, seconds: int = 0) -> None"""
    def __add__(self, arg0: TimeSpan) -> TimePoint:
        """__add__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimeSpan) -> memilio.simulation.abm.TimePoint"""
    def __eq__(self, arg0: TimePoint) -> bool:
        """__eq__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimePoint) -> bool"""
    def __ge__(self, arg0: TimePoint) -> bool:
        """__ge__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimePoint) -> bool"""
    def __gt__(self, arg0: TimePoint) -> bool:
        """__gt__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimePoint) -> bool"""
    def __iadd__(self, arg0: TimeSpan) -> TimePoint:
        """__iadd__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimeSpan) -> memilio.simulation.abm.TimePoint"""
    def __isub__(self, arg0: TimeSpan) -> TimePoint:
        """__isub__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimeSpan) -> memilio.simulation.abm.TimePoint"""
    def __le__(self, arg0: TimePoint) -> bool:
        """__le__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimePoint) -> bool"""
    def __lt__(self, arg0: TimePoint) -> bool:
        """__lt__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimePoint) -> bool"""
    def __ne__(self, arg0: TimePoint) -> bool:
        """__ne__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimePoint) -> bool"""
    @overload
    def __sub__(self, arg0: TimePoint) -> TimeSpan:
        """__sub__(*args, **kwargs)
        Overloaded function.

        1. __sub__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimePoint) -> memilio.simulation.abm.TimeSpan

        2. __sub__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimeSpan) -> memilio.simulation.abm.TimePoint
        """
    @overload
    def __sub__(self, arg0: TimeSpan) -> TimePoint:
        """__sub__(*args, **kwargs)
        Overloaded function.

        1. __sub__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimePoint) -> memilio.simulation.abm.TimeSpan

        2. __sub__(self: memilio.simulation.abm.TimePoint, arg0: memilio.simulation.abm.TimeSpan) -> memilio.simulation.abm.TimePoint
        """
    @property
    def day_of_week(self) -> int: ...
    @property
    def days(self) -> float: ...
    @property
    def hour_of_day(self) -> int: ...
    @property
    def hours(self) -> float: ...
    @property
    def seconds(self) -> int: ...
    @property
    def time_since_midnight(self) -> TimeSpan: ...

class TimeSpan:
    def __init__(self, seconds: int = ...) -> None:
        """__init__(self: memilio.simulation.abm.TimeSpan, seconds: int = 0) -> None"""
    def __add__(self, arg0: TimeSpan) -> TimeSpan:
        """__add__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> memilio.simulation.abm.TimeSpan"""
    def __eq__(self, arg0: TimeSpan) -> bool:
        """__eq__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> bool"""
    def __gt__(self, arg0: TimeSpan) -> bool:
        """__gt__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> bool"""
    def __iadd__(self, arg0: TimeSpan) -> TimeSpan:
        """__iadd__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> memilio.simulation.abm.TimeSpan"""
    def __imul__(self, arg0: int) -> TimeSpan:
        """__imul__(self: memilio.simulation.abm.TimeSpan, arg0: int) -> memilio.simulation.abm.TimeSpan"""
    def __isub__(self, arg0: TimeSpan) -> TimeSpan:
        """__isub__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> memilio.simulation.abm.TimeSpan"""
    def __itruediv__(self, arg0: int) -> TimeSpan:
        """__itruediv__(self: memilio.simulation.abm.TimeSpan, arg0: int) -> memilio.simulation.abm.TimeSpan"""
    @overload
    def __le__(self, arg0: TimeSpan) -> bool:
        """__le__(*args, **kwargs)
        Overloaded function.

        1. __le__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> bool

        2. __le__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> bool
        """
    @overload
    def __le__(self, arg0: TimeSpan) -> bool:
        """__le__(*args, **kwargs)
        Overloaded function.

        1. __le__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> bool

        2. __le__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> bool
        """
    def __lt__(self, arg0: TimeSpan) -> bool:
        """__lt__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> bool"""
    def __mul__(self, arg0: int) -> TimeSpan:
        """__mul__(self: memilio.simulation.abm.TimeSpan, arg0: int) -> memilio.simulation.abm.TimeSpan"""
    def __ne__(self, arg0: TimeSpan) -> bool:
        """__ne__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> bool"""
    def __sub__(self, arg0: TimeSpan) -> TimeSpan:
        """__sub__(self: memilio.simulation.abm.TimeSpan, arg0: memilio.simulation.abm.TimeSpan) -> memilio.simulation.abm.TimeSpan"""
    def __truediv__(self, arg0: int) -> TimeSpan:
        """__truediv__(self: memilio.simulation.abm.TimeSpan, arg0: int) -> memilio.simulation.abm.TimeSpan"""
    @property
    def days(self) -> float: ...
    @property
    def hours(self) -> float: ...
    @property
    def seconds(self) -> int: ...

class Trip:
    cells: list[int]
    destination: LocationId
    destination_type: LocationType
    origin: LocationId
    person_id: PersonId
    time: TimePoint
    def __init__(self, person_id: int, time: TimePoint, destination: LocationId, origin: LocationId, type_of_activity: LocationType, cells: list[int] = ...) -> None:
        """__init__(self: memilio.simulation.abm.Trip, person_id: int, time: memilio.simulation.abm.TimePoint, destination: memilio.simulation.abm.LocationId, origin: memilio.simulation.abm.LocationId, type_of_activity: memilio.simulation.abm.LocationType, cells: List[int] = []) -> None"""

class TripList:
    def __init__(self) -> None:
        """__init__(self: memilio.simulation.abm.TripList) -> None"""
    def add_trip(self, trip: Trip, weekend: bool = ...) -> None:
        """add_trip(self: memilio.simulation.abm.TripList, trip: memilio.simulation.abm.Trip, weekend: bool = False) -> None"""
    def next_trip(self, weekend: bool = ...) -> Trip:
        """next_trip(self: memilio.simulation.abm.TripList, weekend: bool = False) -> memilio.simulation.abm.Trip"""
    def num_trips(self, weekend: bool = ...) -> int:
        """num_trips(self: memilio.simulation.abm.TripList, weekend: bool = False) -> int"""

class VirusVariant:
    __members__: ClassVar[dict] = ...  # read-only
    Wildtype: ClassVar[VirusVariant] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None:
        """__init__(self: memilio.simulation.abm.VirusVariant, value: int) -> None"""
    @staticmethod
    def values() -> _VirusVariantValues:
        """values() -> memilio.simulation.abm._VirusVariantValues"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object) -> int"""
    def __index__(self) -> int:
        """__index__(self: memilio.simulation.abm.VirusVariant) -> int"""
    def __int__(self) -> int:
        """__int__(self: memilio.simulation.abm.VirusVariant) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object) -> bool"""
    def __next__(self) -> VirusVariant:
        """__next__(self: memilio.simulation.abm.VirusVariant) -> memilio.simulation.abm.VirusVariant"""
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class _AgeParameterArray:
    @overload
    def __init__(self, arg0: MultiIndex__AgeParameterArray, arg1: memilio.simulation.UncertainValue) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __init__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray) -> None
        """
    @overload
    def __init__(self, arg0: MultiIndex__AgeParameterArray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __init__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray) -> None
        """
    def get_flat_index(self, arg0: MultiIndex__AgeParameterArray) -> int:
        """get_flat_index(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray) -> int"""
    def numel(self) -> int:
        """numel(self: memilio.simulation.abm._AgeParameterArray) -> int"""
    def resize(self, arg0: MultiIndex__AgeParameterArray) -> None:
        """resize(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray) -> None"""
    def resize_N3mio3abm12VirusVariantE(self, arg0: Index_N3mio3abm12VirusVariantE) -> None:
        """resize_N3mio3abm12VirusVariantE(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.Index_N3mio3abm12VirusVariantE) -> None"""
    def resize_N3mio8AgeGroupE(self, arg0: memilio.simulation.Index_AgeGroup) -> None:
        """resize_N3mio8AgeGroupE(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.Index_AgeGroup) -> None"""
    def size(self) -> tuple[Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup]:
        """size(self: memilio.simulation.abm._AgeParameterArray) -> Tuple[memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup]"""
    def size_N3mio3abm12VirusVariantE(self) -> Index_N3mio3abm12VirusVariantE:
        """size_N3mio3abm12VirusVariantE(self: memilio.simulation.abm._AgeParameterArray) -> memilio.simulation.abm.Index_N3mio3abm12VirusVariantE"""
    def size_N3mio8AgeGroupE(self) -> memilio.simulation.Index_AgeGroup:
        """size_N3mio8AgeGroupE(self: memilio.simulation.abm._AgeParameterArray) -> memilio.simulation.Index_AgeGroup"""
    @overload
    def __getitem__(self, arg0: MultiIndex__AgeParameterArray) -> memilio.simulation.UncertainValue:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray) -> memilio.simulation.UncertainValue

        2. __getitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup]) -> memilio.simulation.UncertainValue
        """
    @overload
    def __getitem__(self, arg0: tuple[Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup]) -> memilio.simulation.UncertainValue:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray) -> memilio.simulation.UncertainValue

        2. __getitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup]) -> memilio.simulation.UncertainValue
        """
    def __iter__(self) -> Iterator:
        """__iter__(self: memilio.simulation.abm._AgeParameterArray) -> Iterator"""
    @overload
    def __setitem__(self, arg0: MultiIndex__AgeParameterArray, arg1: memilio.simulation.UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: tuple[Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup], arg1: memilio.simulation.UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: memilio.simulation.UncertainValue) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: object, arg1: float) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: float) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: memilio.simulation.abm.MultiIndex__AgeParameterArray, arg1: memilio.simulation.UncertainValue) -> None

        2. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm12VirusVariantE, memilio.simulation.Index_AgeGroup], arg1: memilio.simulation.UncertainValue) -> None

        3. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: object, arg1: memilio.simulation.UncertainValue) -> None

        4. __setitem__(self: memilio.simulation.abm._AgeParameterArray, arg0: object, arg1: float) -> None
        """

class _InfectionStateValues:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __iter__(self) -> InfectionState:
        """__iter__(self: memilio.simulation.abm._InfectionStateValues) -> memilio.simulation.abm.InfectionState"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.abm._InfectionStateValues) -> int"""

class _Iter_ModelLocationsRange:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __next__(self) -> Location:
        """__next__(self: memilio.simulation.abm._Iter_ModelLocationsRange) -> memilio.simulation.abm.Location"""

class _Iter_ModelPersonsRange:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __next__(self) -> Person:
        """__next__(self: memilio.simulation.abm._Iter_ModelPersonsRange) -> memilio.simulation.abm.Person"""

class _LocationTypeValues:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __iter__(self) -> LocationType:
        """__iter__(self: memilio.simulation.abm._LocationTypeValues) -> memilio.simulation.abm.LocationType"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.abm._LocationTypeValues) -> int"""

class _ModelLocationsRange:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __getitem__(self, arg0: int) -> Location:
        """__getitem__(self: memilio.simulation.abm._ModelLocationsRange, arg0: int) -> memilio.simulation.abm.Location"""
    def __iter__(self) -> _Iter_ModelLocationsRange:
        """__iter__(self: memilio.simulation.abm._ModelLocationsRange) -> memilio.simulation.abm._Iter_ModelLocationsRange"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.abm._ModelLocationsRange) -> int"""

class _ModelPersonsRange:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __getitem__(self, arg0: int) -> Person:
        """__getitem__(self: memilio.simulation.abm._ModelPersonsRange, arg0: int) -> memilio.simulation.abm.Person"""
    def __iter__(self) -> _Iter_ModelPersonsRange:
        """__iter__(self: memilio.simulation.abm._ModelPersonsRange) -> memilio.simulation.abm._Iter_ModelPersonsRange"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.abm._ModelPersonsRange) -> int"""

class _ProtectionTypeValues:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __iter__(self) -> ProtectionType:
        """__iter__(self: memilio.simulation.abm._ProtectionTypeValues) -> memilio.simulation.abm.ProtectionType"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.abm._ProtectionTypeValues) -> int"""

class _TestData:
    @overload
    def __init__(self, arg0: Index_N3mio3abm8TestTypeE, arg1: TestParameters) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE, arg1: memilio.simulation.abm.TestParameters) -> None

        2. __init__(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE) -> None
        """
    @overload
    def __init__(self, arg0: Index_N3mio3abm8TestTypeE) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE, arg1: memilio.simulation.abm.TestParameters) -> None

        2. __init__(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE) -> None
        """
    def get_flat_index(self, arg0: Index_N3mio3abm8TestTypeE) -> int:
        """get_flat_index(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE) -> int"""
    def numel(self) -> int:
        """numel(self: memilio.simulation.abm._TestData) -> int"""
    def resize(self, arg0: Index_N3mio3abm8TestTypeE) -> None:
        """resize(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE) -> None"""
    def size(self) -> Index_N3mio3abm8TestTypeE:
        """size(self: memilio.simulation.abm._TestData) -> memilio.simulation.abm.Index_N3mio3abm8TestTypeE"""
    def size_N3mio3abm8TestTypeE(self) -> Index_N3mio3abm8TestTypeE:
        """size_N3mio3abm8TestTypeE(self: memilio.simulation.abm._TestData) -> memilio.simulation.abm.Index_N3mio3abm8TestTypeE"""
    @overload
    def __getitem__(self, arg0: Index_N3mio3abm8TestTypeE) -> TestParameters:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE) -> memilio.simulation.abm.TestParameters

        2. __getitem__(self: memilio.simulation.abm._TestData, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm8TestTypeE]) -> memilio.simulation.abm.TestParameters
        """
    @overload
    def __getitem__(self, arg0: tuple[Index_N3mio3abm8TestTypeE]) -> TestParameters:
        """__getitem__(*args, **kwargs)
        Overloaded function.

        1. __getitem__(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE) -> memilio.simulation.abm.TestParameters

        2. __getitem__(self: memilio.simulation.abm._TestData, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm8TestTypeE]) -> memilio.simulation.abm.TestParameters
        """
    def __iter__(self) -> Iterator:
        """__iter__(self: memilio.simulation.abm._TestData) -> Iterator"""
    @overload
    def __setitem__(self, arg0: Index_N3mio3abm8TestTypeE, arg1: TestParameters) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE, arg1: memilio.simulation.abm.TestParameters) -> None

        2. __setitem__(self: memilio.simulation.abm._TestData, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm8TestTypeE], arg1: memilio.simulation.abm.TestParameters) -> None

        3. __setitem__(self: memilio.simulation.abm._TestData, arg0: object, arg1: memilio.simulation.abm.TestParameters) -> None
        """
    @overload
    def __setitem__(self, arg0: tuple[Index_N3mio3abm8TestTypeE], arg1: TestParameters) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE, arg1: memilio.simulation.abm.TestParameters) -> None

        2. __setitem__(self: memilio.simulation.abm._TestData, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm8TestTypeE], arg1: memilio.simulation.abm.TestParameters) -> None

        3. __setitem__(self: memilio.simulation.abm._TestData, arg0: object, arg1: memilio.simulation.abm.TestParameters) -> None
        """
    @overload
    def __setitem__(self, arg0: object, arg1: TestParameters) -> None:
        """__setitem__(*args, **kwargs)
        Overloaded function.

        1. __setitem__(self: memilio.simulation.abm._TestData, arg0: memilio.simulation.abm.Index_N3mio3abm8TestTypeE, arg1: memilio.simulation.abm.TestParameters) -> None

        2. __setitem__(self: memilio.simulation.abm._TestData, arg0: Tuple[memilio.simulation.abm.Index_N3mio3abm8TestTypeE], arg1: memilio.simulation.abm.TestParameters) -> None

        3. __setitem__(self: memilio.simulation.abm._TestData, arg0: object, arg1: memilio.simulation.abm.TestParameters) -> None
        """

class _TestTypeValues:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __iter__(self) -> TestType:
        """__iter__(self: memilio.simulation.abm._TestTypeValues) -> memilio.simulation.abm.TestType"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.abm._TestTypeValues) -> int"""

class _VirusVariantValues:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __iter__(self) -> VirusVariant:
        """__iter__(self: memilio.simulation.abm._VirusVariantValues) -> memilio.simulation.abm.VirusVariant"""
    def __len__(self) -> int:
        """__len__(self: memilio.simulation.abm._VirusVariantValues) -> int"""

def days(arg0: int) -> TimeSpan:
    """days(arg0: int) -> memilio.simulation.abm.TimeSpan"""
def hours(arg0: int) -> TimeSpan:
    """hours(arg0: int) -> memilio.simulation.abm.TimeSpan"""
def minutes(arg0: int) -> TimeSpan:
    """minutes(arg0: int) -> memilio.simulation.abm.TimeSpan"""
def seconds(arg0: int) -> TimeSpan:
    """seconds(arg0: int) -> memilio.simulation.abm.TimeSpan"""

