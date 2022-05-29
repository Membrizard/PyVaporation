from pytest import fixture

from component import AllComponents
from conditions import Conditions
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixture import AllMixtures, Composition, CompositionType
from pervaporation import Pervaporation


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def all_mixtures(all_components):
    return AllMixtures.load("mixtures.yml", all_components)


@fixture
def romakon_pm102_real(all_components):
    experiment_h2o_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=all_components.h2o,
        permeance=0.027,
        activation_energy=19944,
    )

    experiment_etoh_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=all_components.etoh,
        permeance=0.0000282,
        activation_energy=110806,
    )

    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_h2o_1,
            experiment_etoh_1,
        ]
    )

    return Membrane(ideal_experiments=ideal_experiments, name="Romakon-PM102")


@fixture
def test_conditions():
    return Conditions(
        membrane_area=0.04155,
        feed_temperature=333.15,
        permeate_temperature=293.15,
        feed_amount=72,
        initial_feed_composition=Composition(p=0.94, type=CompositionType("weight")),
    )


@fixture()
def pervaporation(romakon_pm102_real, all_mixtures):
    return Pervaporation(
        membrane=romakon_pm102_real,
        mixture=all_mixtures.h2o_etoh,
    )


def test_ideal_non_isothermal_process(pervaporation, test_conditions):
    model = pervaporation.model_ideal_non_isothermal_process(
        conditions=test_conditions, number_of_steps=8, delta_hours=0.125
    )

    validation_permeances_h2o = [0.034, 0.034, 0.034, 0.033, 0.033, 0.033, 0.032, 0.032]
    validation_fluxes_h2o = [
        1.480173991,
        1.429324803,
        1.381641193,
        1.336840606,
        1.294672953,
        1.254916091,
        1.217372039,
        1.181863786,
    ]

    validation_temperatures = [60, 59.5, 59.1, 58.7, 58.2, 57.8, 57.5, 57.1]
    validation_evaporation_heat = [
        40.68,
        40.68,
        40.68,
        40.69,
        40.69,
        40.69,
        40.70,
        40.70,
    ]

    for i in range(len(validation_permeances_h2o)):
        assert abs(validation_permeances_h2o[i] - model.permeances[i][0]) < 0
