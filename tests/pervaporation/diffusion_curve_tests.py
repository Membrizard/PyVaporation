from component import AllComponents
from conditions import Conditions
from experiments import IdealExperiments, IdealExperiment
from membrane import Membrane
from diffusion_curve import DiffusionCurve
from mixture import (
    AllMixtures,
    Composition,
    CompositionType,
    get_nrtl_partial_pressures,
)
from pervaporation import Pervaporation
from pytest import fixture


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def all_mixtures(all_components):
    return AllMixtures.load("mixtures.yml", all_components)


@fixture
def pai_spi(all_components, all_mixtures):
    experiment_meoh_1 = IdealExperiment(
        name="PAI_SPI_1(wt%)_asym",
        temperature=293.15,
        component=all_components.meoh,
        permeance=0.0267,
        activation_energy=20904.65064,
    )

    experiment_meoh_2 = IdealExperiment(
        name="PAI_SPI_1(wt%)_asym",
        temperature=313.15,
        component=all_components.meoh,
        permeance=0.04498,
        activation_energy=20904.65064,
    )

    experiment_meoh_3 = IdealExperiment(
        name="PAI_SPI_1(wt%)_asym",
        temperature=325.15,
        component=all_components.meoh,
        permeance=0.06246,
        activation_energy=20904.65064,
    )

    experiment_mtbe_1 = IdealExperiment(
        name="PAI_SPI_1(wt%)_asym",
        temperature=293.15,
        component=all_components.mtbe,
        permeance=0.01077,
        activation_energy=-39656.76576,
    )

    experiment_mtbe_2 = IdealExperiment(
        name="PAI_SPI_1(wt%)_asym",
        temperature=313.15,
        component=all_components.mtbe,
        permeance=0.00425,
        activation_energy=-39656.76576,
    )

    experiment_mtbe_3 = IdealExperiment(
        name="PAI_SPI_1(wt%)_asym",
        temperature=325.15,
        component=all_components.mtbe,
        permeance=0.00212,
        activation_energy=-39656.76576,
    )

    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_meoh_1,
            experiment_meoh_2,
            experiment_meoh_3,
            experiment_mtbe_1,
            experiment_mtbe_2,
            experiment_mtbe_3,
        ]
    )

    return Membrane(ideal_experiments=ideal_experiments, name="PAI_SPI_1(wt%)_asym")


@fixture
def meoh_mtbe_pervaporation(pai_spi, all_mixtures):
    test_conditions = Conditions(
        membrane_area=0.05,
        feed_temperature=323.15,
        permeate_temperature=1,
        feed_amount=1,
        initial_feed_composition=Composition(p=0.15, type=CompositionType("weight")),
    )
    return Pervaporation(membrane=pai_spi, mixture=all_mixtures.meoh_mtbe)


def test_get_basic_ideal_diffusion_curve(meoh_mtbe_pervaporation):
    feed_compositions = [
        Composition(p=(0.15 - i / 100), type=CompositionType("weight"))
        for i in range(15)
    ]
    modelled_curve = meoh_mtbe_pervaporation.get_ideal_diffusion_curve(
        feed_temperature=325.45,
        compositions=feed_compositions,
    )

    validation_fluxes = [
        (2.1455396637572135, 0.1498696525406852),
        (2.0877822978747695, 0.15178884869862294),
        (2.0251498762973883, 0.15378559236783537),
        (1.9568910640086834, 0.15587023268464897),
        (1.8821020223745228, 0.158054700141235),
        (1.7996878908215999, 0.16035281779686067),
        (1.7083126086484488, 0.16278068604796037),
        (1.6063329395891173, 0.16535716142772208),
        (1.4917108649143909, 0.16810445654515202),
        (1.3618960177123687, 0.17104889740066528),
        (1.2136661171946321, 0.17422188699827237),
        (1.042907747850469, 0.17766114200295094),
        (0.8443112067671265, 0.18141229456438168),
        (0.610939675677, 0.18553098801823115),
        (0.33361155071845805, 0.19008564869920502),
    ]
    for i in range(len(feed_compositions)):
        assert abs(modelled_curve.partial_fluxes[i][0] - validation_fluxes[i][0]) < 1e-5
        assert abs(modelled_curve.partial_fluxes[i][1] - validation_fluxes[i][1]) < 1e-5
        assert (
            abs(modelled_curve.feed_compositions[i].first - feed_compositions[i].first)
            < 1e-5
        )
