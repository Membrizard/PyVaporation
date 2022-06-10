import pytest
from pytest import fixture

from component import AllComponents
from conditions import CalculationType, Conditions, TemperatureProgram
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixture import AllMixtures, Composition, CompositionType
from permeance import Permeance
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def all_mixtures(all_components):
    return AllMixtures.load("mixtures.yml", all_components)


@fixture
def romakon_pm102_real(all_components, all_mixtures):
    experiment_h2o_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=all_components.h2o,
        permeance=Permeance(0.036091),
        activation_energy=19944,
    )

    experiment_etoh_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=all_components.etoh,
        permeance=Permeance(0.0000282),
        activation_energy=110806,
    )

    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_h2o_1,
            experiment_etoh_1,
        ]
    )

    membrane = Membrane(ideal_experiments=ideal_experiments, name="Romakon-PM102")
    diffusion_curve = Pervaporation(
        membrane=membrane,
        mixture=all_mixtures.h2o_etoh,
    ).ideal_diffusion_curve(
        compositions=[
            Composition(p=i / 10, type=CompositionType("weight")) for i in range(1, 10)
        ],
        feed_temperature=323.15,
    )

    return Membrane(
        ideal_experiments=ideal_experiments,
        diffusion_curve_sets=[DiffusionCurveSet("ideal curve", [diffusion_curve])],
        name="Romakon-PM102",
    )


@fixture
def test_conditions():
    return Conditions(
        membrane_area=0.04155,
        initial_feed_temperature=333.15,
        permeate_temperature=293.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType("weight")),
    )


@fixture()
def pervaporation(romakon_pm102_real, all_mixtures):
    return Pervaporation(
        membrane=romakon_pm102_real,
        mixture=all_mixtures.h2o_etoh,
    )


def test_validate_against_ideal_process():
    assert 0 == 0


def test_non_ideal_isothermal_process():
    assert 0 == 0
