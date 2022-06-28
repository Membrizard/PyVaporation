from pathlib import Path

from diffusion_curve import DiffusionCurveSet
from experiments import IdealExperiments
from membrane import Membrane


def test_read_ideal_experiment():
    membrane = Membrane.load(Path("tests/default_membranes/Membrane_RomakonPM_102"), create_results_dir=False)
    assert len(membrane.ideal_experiments) == 6


def test_read_diffusion_curves_set():
    membrane = Membrane.load(Path("tests/default_membranes/Membrane_RomakonPM_102"), create_results_dir=False)
    assert len(membrane.diffusion_curve_sets) > 0


def test_membrane_load():
    Membrane.load(Path("tests/default_membranes/Pervap_4101"), create_results_dir=False)
    Membrane.load(Path("tests/default_membranes/Pervap_2510"), create_results_dir=False)

