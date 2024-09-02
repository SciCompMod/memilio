from dataclasses import dataclass, field


@dataclass
class InferenceConfig:
    T: int
    N: int
    intervention_model: bool = False
    observation_model: bool = False
    obs_data: list[float] = field(default_factory=list, init=False)

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, item, value):
        setattr(self, item, value)
