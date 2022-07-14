import numpy
from pytest import fixture

from pyvaporation.components import Components
from pyvaporation.conditions import (CalculationType, Conditions,
                                     TemperatureProgram)
from pyvaporation.experiments import IdealExperiment, IdealExperiments
from pyvaporation.membrane import Membrane
from pyvaporation.mixtures import Composition, CompositionType, Mixtures
from pyvaporation.permeance import Permeance
from pyvaporation.pervaporation import Pervaporation


@fixture
def romakon_pm102_real():
    experiment_h2o_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=Components.H2O,
        permeance=Permeance(0.036091),
        activation_energy=19944,
    )

    experiment_etoh_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=Components.EtOH,
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
def test_conditions_isothermal():
    return Conditions(
        membrane_area=0.04155,
        initial_feed_temperature=333.15,
        permeate_temperature=293.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType.weight),
        temperature_program=TemperatureProgram(
            coefficients=[333.15], type=CalculationType.polynomial
        ),
    )


@fixture
def test_conditions_linear():
    return Conditions(
        membrane_area=0.04155,
        initial_feed_temperature=333.15,
        permeate_temperature=293.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType.weight),
        temperature_program=TemperatureProgram(
            coefficients=[333.15, -1], type=CalculationType.polynomial
        ),
    )


@fixture
def test_conditions_exp():
    return Conditions(
        membrane_area=0.04155,
        initial_feed_temperature=333.15,
        permeate_temperature=293.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType.weight),
        temperature_program=TemperatureProgram(
            coefficients=[333.15, 0, -1e-3], type=CalculationType.exponential
        ),
    )


@fixture
def test_conditions_ln():
    return Conditions(
        membrane_area=0.04155,
        initial_feed_temperature=333.15,
        permeate_temperature=293.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType.weight),
        temperature_program=TemperatureProgram(
            coefficients=[333.15, numpy.exp(1), 1e-2],
            type=CalculationType.logarithmic,
        ),
    )


@fixture()
def pervaporation(romakon_pm102_real):
    return Pervaporation(
        membrane=romakon_pm102_real,
        mixture=Mixtures.H2O_EtOH,
    )


def test_isothermal_temperature_program(pervaporation, test_conditions_isothermal):
    model = pervaporation.ideal_non_isothermal_process(
        conditions=test_conditions_isothermal, number_of_steps=8, delta_hours=0.125
    )

    validation_temperatures = [333.15] * 8
    validation_permeances = [0.0450987] * 8

    for i in range(len(validation_temperatures)):
        assert abs((validation_temperatures[i]) - model.feed_temperature[i]) == 0
        assert abs(validation_permeances[i] - model.permeances[i][0].value) < 1e-6


def test_linear_temperature_program(pervaporation, test_conditions_linear):
    model = pervaporation.ideal_non_isothermal_process(
        conditions=test_conditions_linear, number_of_steps=8, delta_hours=1
    )

    validation_temperatures = [
        333.15,
        332.15,
        331.15,
        330.15,
        329.15,
        328.15,
        327.15,
        326.15,
    ]
    validation_permeances = [
        0.0450987,
        0.0441316,
        0.0431796,
        0.04224253,
        0.041320303,
        0.04041277,
        0.0395198,
        0.038641277,
    ]

    for i in range(len(validation_temperatures)):
        assert abs((validation_temperatures[i]) - model.feed_temperature[i]) == 0
        assert abs(validation_permeances[i] - model.permeances[i][0].value) < 1e-6


def test_exp_temperature_program(pervaporation, test_conditions_exp):
    model = pervaporation.ideal_non_isothermal_process(
        conditions=test_conditions_exp, number_of_steps=8, delta_hours=1
    )

    validation_temperatures = [333.15 * numpy.exp(-1e-3 * i) for i in range(8)]
    validation_permeances = [
        0.0450987,
        0.0447750,
        0.0444533,
        0.0441336,
        0.0438158,
        0.0435000,
        0.0431863,
        0.0428744,
    ]

    for i in range(len(validation_temperatures)):
        assert abs((validation_temperatures[i]) - model.feed_temperature[i]) == 0
        assert abs(validation_permeances[i] - model.permeances[i][0].value) < 1e-6


def test_ln_temperature_program(pervaporation, test_conditions_ln):
    model = pervaporation.ideal_non_isothermal_process(
        conditions=test_conditions_ln, number_of_steps=8, delta_hours=1
    )

    validation_temperatures = [
        333.15 * numpy.log(1e-2 * i + numpy.exp(1)) for i in range(8)
    ]
    validation_permeances = [
        0.0450987,
        0.0463025,
        0.0475247,
        0.0487655,
        0.0500247,
        0.0513025,
        0.0525987,
        0.0539135,
    ]

    for i in range(len(validation_temperatures)):
        assert abs((validation_temperatures[i]) - model.feed_temperature[i]) < 1e-7
        assert abs(validation_permeances[i] - model.permeances[i][0].value) < 1e-7
