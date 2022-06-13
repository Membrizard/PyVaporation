import pytest
from pytest import fixture

from components import AllComponents
from conditions import CalculationType, Conditions, TemperatureProgram
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixtures import AllMixtures, Composition, CompositionType
from permeance import Permeance
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet


import pytest
from pytest import fixture

from components import AllComponents
from conditions import Conditions
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixtures import AllMixtures, Composition, CompositionType
from permeance import Permeance
from pervaporation import Pervaporation


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def all_mixtures(all_components):
    return AllMixtures.load("mixtures.yml", all_components)


@fixture
def pervaporation(all_components, all_mixtures):
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

    membrane = Membrane(
        ideal_experiments=ideal_experiments,
        diffusion_curve_sets=[DiffusionCurveSet("ideal curve", [diffusion_curve])],
        name="Romakon-PM102",
    )

    return Pervaporation(membrane, all_mixtures.h2o_etoh)


@fixture()
def conditions():
    return Conditions(
        membrane_area=0.4155,
        initial_feed_temperature=333.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType("weight")),
    )


def test_validate_against_ideal_process(pervaporation, conditions):
    number_of_steps = 200
    delta_hours = 0.125
    ideal_model = pervaporation.ideal_isothermal_process(
        conditions=conditions, number_of_steps=number_of_steps, delta_hours=delta_hours
    )
    non_ideal_model = pervaporation.non_ideal_isothermal_process(
        conditions=conditions,
        diffusion_curves=pervaporation.membrane.diffusion_curve_sets[0],
        number_of_steps=number_of_steps,
        delta_hours=delta_hours,
        n_first=0,
        m_first=0,
        n_second=0,
        m_second=0,
    )

    for i in range(number_of_steps):
        assert (
            abs(non_ideal_model.partial_fluxes[i][0] - ideal_model.partial_fluxes[i][0])
            < 8e-3
        )
        assert (
            abs(non_ideal_model.partial_fluxes[i][1] - ideal_model.partial_fluxes[i][1])
            < 8e-3
        )
        assert (
            abs(
                non_ideal_model.permeances[i][0].value
                - ideal_model.permeances[i][0].value
            )
            < 5e-4
        )
        assert (
            abs(
                non_ideal_model.permeances[i][1].value
                - ideal_model.permeances[i][1].value
            )
            < 5e-4
        )
        assert abs(non_ideal_model.feed_mass[i] - ideal_model.feed_mass[i]) < 6.5e-2
        assert (
            abs(
                non_ideal_model.feed_evaporation_heat[i]
                - ideal_model.feed_evaporation_heat[i]
            )
            < 1
        )
        assert (
            non_ideal_model.permeate_condensation_heat[i]
            == ideal_model.permeate_condensation_heat[i]
        )
        assert (
            abs(
                non_ideal_model.feed_composition[i].first
                - ideal_model.feed_composition[i].first
            )
            < 1e-2
        )
