import os
import logging
from dataclasses import dataclass, field
from dataclasses_json import dataclass_json


@dataclass
class TrainerParameters:
    epochs: int
    summary_dim: int
    num_coupling_layers: int
    coupling_design: str = "affine"

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, item, value):
        setattr(self, item, value)


@dataclass_json
@dataclass
class InferenceConfig:
    T: int
    N: int
    trainer_parameters: TrainerParameters
    intervention_model: bool = False
    observation_model: bool = False
    obs_data: list[float] = field(default_factory=list, init=False)

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, item, value):
        setattr(self, item, value)

    def save(self, path, overwrite=True):
        # Logger init
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)

        filepath = os.path.join(path, "InferenceConfig.json")
        if not overwrite and os.path.exists(filepath):
            logger.info(
                f"Could not save InferenceConfig, because file {filepath} already exists")
            return
        with open(filepath, "w") as f:
            f.write(self.to_json())

    @ classmethod
    def load(cls, path):
        # Logger init
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)

        with open(os.path.join(path, "InferenceConfig.json")) as f:
            json_string = f.read()

        logger.info(f"Loaded loss history from {path}.")
        return cls.from_json(json_string)
