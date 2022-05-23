from component import AllComponents
from experiments import IdealExperiments, IdealExperiment
from mixture import AllMixtures
from membrane import Membrane
from pytest import fixture


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def component_list(all_components):
    return [all_components.h2o, all_components.meoh, all_components.etoh]


@fixture
def romakon_pm102t(all_components):
    experiment_h2o_1 = IdealExperiment(
        name="Romakon PM102T",
        temperature=333,
        component=all_components.h2o,
        permeance=0.0449064,
        activation_energy=-23600,
    )
    experiment_h2o_2 = IdealExperiment(
        name="Romakon PM102T",
        temperature=343,
        component=all_components.h2o,
        permeance=0.0348624,
        activation_energy=-23600,
    )
    experiment_h2o_3 = IdealExperiment(
        name="Romakon PM102T",
        temperature=353,
        component=all_components.h2o,
        permeance=0.0282528,
        activation_energy=-23600,
    )
    experiment_etoh_1 = IdealExperiment(
        name="Romakon PM102T",
        temperature=333,
        component=all_components.etoh,
        permeance=0.0004743,
        activation_energy=-12600,
    )
    experiment_etoh_2 = IdealExperiment(
        name="Romakon PM102T",
        temperature=343,
        component=all_components.etoh,
        permeance=0.0004428,
        activation_energy=-12600,
    )
    experiment_etoh_3 = IdealExperiment(
        name="Romakon PM102T",
        temperature=353,
        component=all_components.etoh,
        permeance=0.0003698,
        activation_energy=-12600,
    )
    experiment_meoh_1 = IdealExperiment(
        name="Romakon PM102T",
        temperature=343,
        component=all_components.meoh,
        permeance=0.0018311,
    )
    experiment_meoh_2 = IdealExperiment(
        name="Romakon PM102T",
        temperature=353,
        component=all_components.meoh,
        permeance=0.0017302,
    )
    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_h2o_1,
            experiment_h2o_2,
            experiment_h2o_3,
            experiment_etoh_1,
            experiment_etoh_2,
            experiment_etoh_3,
            experiment_meoh_1,
            experiment_meoh_2,
        ]
    )

    return Membrane(ideal_experiments=ideal_experiments, name="Romakon-PM102T")


def test_get_penetrant_data(romakon_pm102t, component_list):
    h2o_experiments = romakon_pm102t.get_penetrant_data(component_list[0])
    meoh_experiments = romakon_pm102t.get_penetrant_data(component_list[1])
    etoh_experiments = romakon_pm102t.get_penetrant_data(component_list[2])

    for i in range(len(h2o_experiments)):
        assert h2o_experiments.experiments[i].component == component_list[0]

    for i in range(len(meoh_experiments)):
        assert meoh_experiments.experiments[i].component == component_list[1]

    for i in range(len(etoh_experiments)):
        assert etoh_experiments.experiments[i].component == component_list[2]


def test_calculate_activation_energy(romakon_pm102t, component_list):
    validation_energy_h2o = abs(
        romakon_pm102t.get_penetrant_data(component_list[0])
        .experiments[0]
        .activation_energy
    )
    validation_energy_etoh = abs(
        romakon_pm102t.get_penetrant_data(component_list[2])
        .experiments[0]
        .activation_energy
    )

    assert (
        abs(
            abs(romakon_pm102t.calculate_activation_energy(component_list[0]))
            - validation_energy_h2o
        )
        < validation_energy_h2o * 0.04
    )
    assert (
        abs(
            abs(romakon_pm102t.calculate_activation_energy(component_list[2]))
            - validation_energy_etoh
        )
        < validation_energy_etoh * 0.04
    )

    assert (
        abs(
            abs(romakon_pm102t.calculate_activation_energy(component_list[1]))
            - 0
        )
        < validation_energy_etoh * 0.04
    )
