from pytest import fixture

from components import Components
from conditions import Conditions, TemperatureProgram
from diffusion_curve import DiffusionCurveSet
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixtures import Composition, CompositionType, Mixtures
from permeance import Permeance
from pervaporation import Pervaporation
from config import Config
from pathlib import Path
from optimizer import find_best_fit
from matplotlib.testing.decorators import image_comparison

@fixture
def membrane():
    config = Config(
        source_path=Path("tests/data/Pervap_2510"),
    )

    return Membrane.load(config)


@fixture
def diffusin_curve_set(membrane):
    return membrane.diffusion_curve_sets[0]


@fixture
def process_model(membrane):
    pv = Pervaporation(
        membrane=membrane,
        mixture=Mixtures.H2O_iPOH,
    )

    con = Conditions(
        membrane_area=0.017,
        initial_feed_temperature=368.15,
        initial_feed_amount=1.5,
        initial_feed_composition=Composition(p=0.15, type=CompositionType.weight),
        permeate_pressure=0,
    )

    return pv.non_ideal_isothermal_process(
        conditions=con,
        diffusion_curve_set=membrane.diffusion_curve_sets[0],
        number_of_steps=50,
        delta_hours=0.2,
    )


def test_plot_diffusion_curve():
    assert 0 == 0


def test_plot_pervaporation_function():
    assert 0 == 0


def test_plot_process_model():
    assert 0 == 0
