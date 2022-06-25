import pytest
from pytest import fixture

from components import Components
from conditions import Conditions
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixtures import Composition, CompositionType, Mixtures
from permeance import Permeance
from pervaporation import Pervaporation


@fixture
def romakon_al2():
    experiment_h2o_1 = IdealExperiment(
        name="Romakon-Al2",
        temperature=319.65,
        component=Components.H2O,
        permeance=Permeance(0.016876),
    )
    experiment_etoh_1 = IdealExperiment(
        name="Romakon-Al2",
        temperature=319.65,
        component=Components.EtOH,
        permeance=Permeance(0.000139),
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
        initial_feed_temperature=319.65,
        initial_feed_amount=0.047,
        initial_feed_composition=Composition(p=0.04, type=CompositionType.weight),
        permeate_temperature=198,
    )


@fixture
def romakon_al2_experiment_conditions_2():
    return Conditions(
        membrane_area=0.0048,
        initial_feed_temperature=319.65,
        initial_feed_amount=0.047,
        initial_feed_composition=Composition(p=0.04, type=CompositionType.weight),
    )


@fixture
def romakon_al2_experiment_conditions_3():
    return Conditions(
        membrane_area=0.0048,
        initial_feed_temperature=319.65,
        initial_feed_amount=0.047,
        initial_feed_composition=Composition(p=0.04, type=CompositionType.weight),
        permeate_pressure=0,
    )


@fixture
def romakon_al2_experiment_conditions_4():
    return Conditions(
        membrane_area=0.0048,
        initial_feed_temperature=319.65,
        initial_feed_amount=0.047,
        initial_feed_composition=Composition(p=0.04, type=CompositionType.weight),
        permeate_temperature=0,
        permeate_pressure=0,
    )


@fixture
def romakon_al2_pervaporation(
    romakon_al2,
    romakon_al2_experiment_conditions,
):
    return Pervaporation(
        membrane=romakon_al2,
        mixture=Mixtures.H2O_EtOH,
    )


def test_experimet_romakon_al2(
    romakon_al2_pervaporation, romakon_al2_experiment_conditions
):
    model = romakon_al2_pervaporation.ideal_isothermal_process(
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


def test_experimet_romakon_al2_no_permeate_params(
    romakon_al2_pervaporation, romakon_al2_experiment_conditions_2
):
    model = romakon_al2_pervaporation.ideal_isothermal_process(
        number_of_steps=5,
        delta_hours=1,
        conditions=romakon_al2_experiment_conditions_2,
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


def test_experimet_romakon_al2_permeate_pressure(
    romakon_al2_pervaporation, romakon_al2_experiment_conditions_3
):
    model = romakon_al2_pervaporation.ideal_isothermal_process(
        number_of_steps=5,
        delta_hours=1,
        conditions=romakon_al2_experiment_conditions_3,
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


def test_permeate_comp_and_temp(
    romakon_al2_pervaporation, romakon_al2_experiment_conditions_4
):
    with pytest.raises(ValueError):
        romakon_al2_pervaporation.ideal_isothermal_process(
            number_of_steps=5,
            delta_hours=1,
            conditions=romakon_al2_experiment_conditions_4,
            precision=5e-5,
        )
