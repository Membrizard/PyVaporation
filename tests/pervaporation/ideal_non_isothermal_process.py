from pytest import fixture

from component import AllComponents
from conditions import CalculationType, Conditions, TemperatureProgram
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixture import AllMixtures, Composition, CompositionType
from permeance import Permeance
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

    return Membrane(ideal_experiments=ideal_experiments, name="Romakon-PM102")


@fixture
def test_conditions():
    return Conditions(
        membrane_area=0.04155,
        initial_feed_temperature=333.15,
        permeate_temperature=293.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType("weight")),
    )


@fixture
def test_conditions_temp_program():
    return Conditions(
        membrane_area=0.04155,
        initial_feed_temperature=333.15,
        permeate_temperature=293.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType("weight")),
        temperature_program=TemperatureProgram(
            coefficients=[333.15, -1], type=CalculationType("polynomial")
        ),
    )


@fixture()
def pervaporation(romakon_pm102_real, all_mixtures):
    return Pervaporation(
        membrane=romakon_pm102_real,
        mixture=all_mixtures.h2o_etoh,
    )


def test_ideal_non_isothermal_process(pervaporation, test_conditions):
    model = pervaporation.ideal_non_isothermal_process(
        conditions=test_conditions, number_of_steps=8, delta_hours=0.125
    )

    validation_permeances_h2o = [
        0.0451084,
        0.0454345,
        0.0451811,
        0.0449334,
        0.0446911,
        0.0444539,
        0.0442218,
        0.0439945,
    ]
    validation_fluxes_h2o = [0.791, 0.775, 0.761, 0.747, 0.734, 0.721, 0.709, 0.697]

    validation_temperatures = [60, 59.5, 59.1, 58.7, 58.2, 57.8, 57.5, 57.1]

    for i in range(len(validation_permeances_h2o)):
        assert abs(validation_permeances_h2o[i] - model.permeances[i][0].value) < 2e-3
        assert (
            abs((validation_temperatures[i] + 273.15) - model.feed_temperature[i])
            < 3e-1
        )
        assert abs(validation_fluxes_h2o[i] - model.partial_fluxes[i][0]) < 6.3e-2
