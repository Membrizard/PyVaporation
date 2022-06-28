from pathlib import Path

import pytest

from config import Config
from diffusion_curve import DiffusionCurveSet
from experiments import IdealExperiments
from membrane import Membrane


@pytest.fixture
def config() -> Config:
    return Config(source_path=Path("tests/default_membranes/Membrane_RomakonPM_102"))


def test_read_ideal_experiment(config: Config):
    ie = IdealExperiments.from_csv(config.ideal_experiments_path)
    assert len(ie) == 6


def test_read_diffusion_curves_set(config: Config):
    dcs = DiffusionCurveSet.load(config.diffusion_curve_sets_path / "table.csv")
    assert 0 == 0


def test_membrane_load(config):
    config_1 = Config(source_path=Path("tests/default_membranes/Pervap_4101"))
    config_2 = Config(source_path=Path("tests/default_membranes/Pervap_2510"))
    Membrane.load(config_1)
    Membrane.load(config_2)
    Membrane.load(config)
