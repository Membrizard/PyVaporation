from pytest import fixture

from component import AllComponents
from conditions import Conditions
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixture import AllMixtures, Composition, CompositionType
from pervaporation import Pervaporation
from process import ProcessModel


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
        temperature=313.15,
        component=all_components.h2o,
        permeance=0.036091,
        activation_energy=19944,
    )

    experiment_etoh_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=313.15,
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
def romakon_al2(all_components):
    experiment_h2o_1 = IdealExperiment(
        name="Romakon-Al2",
        temperature=319.65,
        component=all_components.h2o,
        permeance=0.016876,
    )
    experiment_etoh_1 = IdealExperiment(
        name="Romakon-Al2",
        temperature=319.65,
        component=all_components.etoh,
        permeance=0.000139,
    )

    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_h2o_1,
            experiment_etoh_1,
        ]
    )

    return Membrane(ideal_experiments=ideal_experiments, name="Romakon-Al2")


@fixture
def romakon_al2_experiment_conditions():
    return Conditions(
        membrane_area=0.0048,
        feed_temperature=319.65,
        feed_amount=0.047,
        initial_feed_composition=Composition(p=0.04, type=CompositionType("weight")),
        permeate_temperature=198,
    )


@fixture
def romakon_pm102_experiment_conditions():
    return Conditions(
        membrane_area=0.01,
        feed_temperature=320.4,
        feed_amount=0.14706,
        initial_feed_composition=Composition(p=0.949, type=CompositionType("weight")),
        permeate_temperature=200,
    )


@fixture
def romakon_al2_pervaporation(
    romakon_al2, romakon_al2_experiment_conditions, all_mixtures
):
    return Pervaporation(
        membrane=romakon_al2,
        mixture=all_mixtures.h2o_etoh,
    )


@fixture
def romakon_pm102_real_pervaporation(
    romakon_pm102_real, romakon_pm102_experiment_conditions, all_mixtures
):
    return Pervaporation(
        membrane=romakon_pm102_real,
        mixture=all_mixtures.h2o_etoh,
    )


def test_experimet_romakon_al2(
    romakon_al2_pervaporation, romakon_al2_experiment_conditions
):
    model = romakon_al2_pervaporation.model_ideal_isothermal_process(
        number_of_steps=5,
        delta_hours=1,
        conditions=romakon_al2_experiment_conditions,
        precision=5e-5,
    )
    experiment_partial_fluxes = [
        (0.0621, 0.0061),
        (0.0399, 0.0039),
        (0.0372, 0.0041),
        (0.0047, 0.0008),
        (0.0056, 0.0010),
    ]
    experiment_permeate_composition = [0.91, 0.91, 0.90, 0.86, 0.86]
    for i in range(5):

        assert (
            abs(model.partial_fluxes[i][0] - experiment_partial_fluxes[i][0]) < 2.5e-2
        )
        assert (
            abs(model.partial_fluxes[i][1] - experiment_partial_fluxes[i][1]) < 3.2e-3
        )
        assert (
            abs(
                model.permeate_composition[i].first - experiment_permeate_composition[i]
            )
            < experiment_permeate_composition[i] * 0.05
        )


def test_experimet_romakon_pm102(romakon_pm102_real_pervaporation, romakon_pm102_experiment_conditions):
    model = romakon_al2_pervaporation.model_ideal_isothermal_process(
        number_of_steps=14,
        d_time_step_hours=1,
        conditions=romakon_al2_experiment_conditions,
        precision=5e-5,
    )
    assert 0 == 0
