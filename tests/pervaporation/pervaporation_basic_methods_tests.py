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
def pai_spi(all_components):
    experiment_meoh_1 = IdealExperiment(
        name="PAI:SPI 1wt%",
        temperature=293.15,
        component=all_components.meoh,
        permeance=0.02634,
        activation_energy=20250,
    )
    experiment_meoh_2 = IdealExperiment(
        name="PAI:SPI 1wt%",
        temperature=313.15,
        component=all_components.meoh,
        permeance=0.04479,
        activation_energy=-23600,
    )
    experiment_mtbe_1 = IdealExperiment(
        name="PAI:SPI 1wt%",
        temperature=293.15,
        component=all_components.mtbe,
        permeance=0.01102,
        activation_energy=-27980,
    )
    experiment_mtbe_2 = IdealExperiment(
        name="PAI:SPI 1wt%",
        temperature=313.15,
        component=all_components.mtbe,
        permeance=0.00529,
        activation_energy=-27980,
    )

    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_meoh_1,
            experiment_meoh_2,
            experiment_mtbe_1,
            experiment_mtbe_2,
        ]
    )

    return Membrane(ideal_experiments=ideal_experiments, name="PAI:SPI 1wt%")


@fixture
def test_conditions():
    return Conditions(
        membrane_area=0.05,
        feed_temperature=323.15,
        permeate_temperature=1,
        feed_amount=1,
        initial_feed_composition=Composition(p=0.15, type=CompositionType("weight")),
    )


@fixture()
def pervaporation(pai_spi, all_mixtures, test_conditions):
    return Pervaporation(
        membrane=pai_spi,
        mixture=all_mixtures.meoh_mtbe,
        conditions=test_conditions,
    )


def test_calculate_partial_fluxes(pervaporation):
    modelled_fluxes = pervaporation.calculate_partial_fluxes(
        feed_temperature=325.5,
        composition=Composition(p=0.03, type=CompositionType("weight")),
        precision=5e-5,
    )

    experimental_fluxes = (0.8282, 0.1818)

    assert abs(modelled_fluxes[0] - experimental_fluxes[0]) < 0
