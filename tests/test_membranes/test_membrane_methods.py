import pytest
from pytest import fixture

from pyvaporation.components import Components
from pyvaporation.experiments import IdealExperiment, IdealExperiments
from pyvaporation.membrane import Membrane
from pyvaporation.permeance import Permeance


@fixture
def component_list():
    return [Components.H2O, Components.MeOH, Components.EtOH]


@fixture
def romakon_pm102t():
    experiment_h2o_1 = IdealExperiment(
        name="Romakon-PM102T",
        temperature=333,
        component=Components.H2O,
        permeance=Permeance(value=0.0449064),
        activation_energy=-23600,
    )
    experiment_h2o_2 = IdealExperiment(
        name="Romakon-PM102T",
        temperature=343,
        component=Components.H2O,
        permeance=Permeance(value=0.0348624),
        activation_energy=-23600,
    )
    experiment_h2o_3 = IdealExperiment(
        name="Romakon-PM102T",
        temperature=353,
        component=Components.H2O,
        permeance=Permeance(value=0.0282528),
        activation_energy=-23600,
    )
    experiment_etoh_1 = IdealExperiment(
        name="Romakon-PM102T",
        temperature=333,
        component=Components.EtOH,
        permeance=Permeance(value=0.0004743),
        activation_energy=-12600,
    )
    experiment_etoh_2 = IdealExperiment(
        name="Romakon-PM102T",
        temperature=343,
        component=Components.EtOH,
        permeance=Permeance(value=0.0004428),
        activation_energy=-12600,
    )
    experiment_etoh_3 = IdealExperiment(
        name="Romakon-PM102T",
        temperature=353,
        component=Components.EtOH,
        permeance=Permeance(value=0.0003698),
        activation_energy=-12600,
    )
    experiment_meoh_1 = IdealExperiment(
        name="Romakon-PM102T",
        temperature=343,
        component=Components.MeOH,
        permeance=Permeance(value=0.0012226),
    )
    experiment_meoh_2 = IdealExperiment(
        name="Romakon-PM102T",
        temperature=353,
        component=Components.MeOH,
        permeance=Permeance(value=0.0014764),
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
    h2o_experiments = romakon_pm102t.get_penetrant_data(Components.H2O)
    meoh_experiments = romakon_pm102t.get_penetrant_data(Components.MeOH)
    etoh_experiments = romakon_pm102t.get_penetrant_data(Components.EtOH)

    for i in range(len(h2o_experiments)):
        assert h2o_experiments.experiments[i].component == component_list[0]

    for i in range(len(meoh_experiments)):
        assert meoh_experiments.experiments[i].component == component_list[1]

    for i in range(len(etoh_experiments)):
        assert etoh_experiments.experiments[i].component == component_list[2]


def test_calculate_activation_energy(romakon_pm102t, component_list):
    validation_energy_h2o = (
        romakon_pm102t.get_penetrant_data(component_list[0])
        .experiments[0]
        .activation_energy
    )
    validation_energy_etoh = (
        romakon_pm102t.get_penetrant_data(component_list[2])
        .experiments[0]
        .activation_energy
    )

    assert abs(
        romakon_pm102t.calculate_activation_energy(component_list[0])
        - validation_energy_h2o
    ) < abs(validation_energy_h2o * 0.04)
    assert abs(
        romakon_pm102t.calculate_activation_energy(component_list[2])
        - validation_energy_etoh
    ) < abs(validation_energy_etoh * 0.04)

    assert (
        abs(romakon_pm102t.calculate_activation_energy(component_list[1]) - 18985.56372)
        < 18985.56372 * 0.0005
    )


def test_get_permeance(romakon_pm102t):

    assert (
        abs(
            romakon_pm102t.get_permeance(343, Components.EtOH).value
            - romakon_pm102t.get_penetrant_data(Components.EtOH)
            .experiments[1]
            .permeance.value
        )
        == 0
    )

    assert (
        abs(romakon_pm102t.get_permeance(348, Components.EtOH).value - 0.0004155648872)
        < 1e-9
    )

    assert (
        abs(romakon_pm102t.get_permeance(349, Components.MeOH).value - 0.0013708803)
        < 1e-9
    )


def test_get_ideal_selectivity(romakon_pm102t):

    assert (
        abs(
            romakon_pm102t.get_ideal_selectivity(
                343, Components.H2O, Components.EtOH, "weight"
            )
            - 78.73170732
        )
        < 1e-5
    )
    assert (
        abs(
            romakon_pm102t.get_ideal_selectivity(
                348, Components.H2O, Components.EtOH, "weight"
            )
            - 74.48720837
        )
        < 1e-5
    )
    assert (
        abs(
            romakon_pm102t.get_ideal_selectivity(
                348, Components.H2O, Components.MeOH, "weight"
            )
            - 23.0084745
        )
        < 1e-5
    )

    assert (
        abs(
            romakon_pm102t.get_ideal_selectivity(343, Components.H2O, Components.EtOH)
            - 201.28578
        )
        < 1e-5
    )
    assert (
        abs(
            romakon_pm102t.get_ideal_selectivity(
                348,
                Components.H2O,
                Components.EtOH,
            )
            - 190.434278
        )
        < 1e-5
    )
    assert (
        abs(
            romakon_pm102t.get_ideal_selectivity(
                348,
                Components.H2O,
                Components.MeOH,
            )
            - 40.9096295
        )
        < 1e-5
    )


def test_get_pure_component_flux(romakon_pm102t):
    validation_fluxes_h2o = [0.726, 0.9923, 1.201]
    validation_temperatures = [333, 343, 353]
    for i in range(len(validation_temperatures)):
        assert (
            abs(
                romakon_pm102t.get_estimated_pure_component_flux(
                    validation_temperatures[i], Components.H2O, 298
                )
                - validation_fluxes_h2o[i]
            )
            < validation_fluxes_h2o[i] * 0.04
        )


def test_get_pure_component_flux_permeate_pressure(romakon_pm102t):
    validation_fluxes_h2o = [0.726, 0.9923, 1.201]
    validation_temperatures = [333, 343, 353]
    for i in range(len(validation_temperatures)):
        assert (
            abs(
                romakon_pm102t.get_estimated_pure_component_flux(
                    validation_temperatures[i],
                    Components.H2O,
                    permeate_pressure=3.165,
                )
                - validation_fluxes_h2o[i]
            )
            < validation_fluxes_h2o[i] * 0.04
        )

    with pytest.raises(ValueError):
        romakon_pm102t.get_estimated_pure_component_flux(
            333, Components.H2O, 298, permeate_pressure=0
        )
