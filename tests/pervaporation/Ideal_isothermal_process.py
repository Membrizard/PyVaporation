from component import AllComponents
from experiments import IdealExperiments, IdealExperiment
from mixture import AllMixtures, Composition, CompositionType
from membrane import Membrane
from pervaporation import Pervaporation
from conditions import Conditions
from pytest import fixture


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def all_mixtures(all_components):
    return AllMixtures.load("mixtures.yml", all_components)


@fixture
def romakon_pm102_binary(all_components):
    experiment_h2o_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=313.15,
        component=all_components.h2o,
        permeance=0.05500,
        activation_energy=19944,
    )
    experiment_h2o_2 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=all_components.h2o,
        permeance=0.06713,
        activation_energy=19944,
    )
    experiment_h2o_3 = IdealExperiment(
        name="Romakon-PM102",
        temperature=333.15,
        component=all_components.h2o,
        permeance=0.08718,
        activation_energy=19944,
    )
    experiment_etoh_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=313.15,
        component=all_components.etoh,
        permeance=0.00002,
        activation_energy=110806,
    )
    experiment_etoh_2 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=all_components.etoh,
        permeance=0.00003,
        activation_energy=110806,
    )
    experiment_etoh_3 = IdealExperiment(
        name="Romakon-PM102",
        temperature=333.15,
        component=all_components.etoh,
        permeance=0.00027,
        activation_energy=110806,
    )

    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_h2o_1,
            experiment_h2o_2,
            experiment_h2o_3,
            experiment_etoh_1,
            experiment_etoh_2,
            experiment_etoh_3,
        ]
    )

    return Membrane(ideal_experiments=ideal_experiments, name="Romakon-PM102")


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
        temperature=0,
        component=all_components.h2o,
        permeance=0.036091,
        activation_energy=19944,
    )



    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_h2o_1,
            experiment_etoh_1,
        ]
    )

    return Membrane(ideal_experiments=ideal_experiments, name="Romakon-PM102")