from pytest import fixture

from component import AllComponents
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixture import AllMixtures, Composition, CompositionType
from pervaporation import Pervaporation
from permeance import Permeance


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
        permeance=Permeance(0.0267),
        activation_energy=20904.65064,
    )

    experiment_meoh_2 = IdealExperiment(
        name="PAI_SPI_1(wt%)_asym",
        temperature=313.15,
        component=all_components.meoh,
        permeance=Permeance(0.04498),
        activation_energy=20904.65064,
    )

    experiment_mtbe_1 = IdealExperiment(
        name="PAI_SPI_1(wt%)_asym",
        temperature=293.15,
        component=all_components.mtbe,
        permeance=Permeance(0.01077),
        activation_energy=-39656.76576,
    )

    experiment_mtbe_2 = IdealExperiment(
        name="PAI_SPI_1(wt%)_asym",
        temperature=313.15,
        component=all_components.mtbe,
        permeance=Permeance(0.00425),
        activation_energy=-39656.76576,
    )

    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_meoh_1,
            experiment_meoh_2,
            experiment_mtbe_1,
            experiment_mtbe_2,
        ]
    )

    return Membrane(ideal_experiments=ideal_experiments, name="PAI_SPI_1(wt%)_asym")


@fixture
def meoh_mtbe_pervaporation(pai_spi, all_mixtures):
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
        (2.077983268028643, 0.171254353278905),
        (2.022044502627796, 0.1734473969756852),
        (1.9613841819297386, 0.17572905333464614),
        (1.8952746281296244, 0.17811114819656712),
        (1.8228404616712275, 0.18060731439962743),
        (1.7430212957481013, 0.18323334739697603),
        (1.6545231380702807, 0.18600764492913047),
        (1.5557544927900093, 0.18895175414283757),
        (1.444741512072742, 0.1920910571346426),
        (1.3190141321579323, 0.19545563632683113),
        (1.1754515318944574, 0.19908137557615233),
        (1.0100698144799902, 0.20301137328818128),
        (0.8177264630935732, 0.20729777280314604),
        (0.5917030783800551, 0.21200415713003018),
        (0.3231071567328018, 0.21720871626593188),
    ]

    assert modelled_curve.membrane_name == "PAI_SPI_1(wt%)_asym"

    for i in range(len(feed_compositions)):
        assert abs(modelled_curve.partial_fluxes[i][0] - validation_fluxes[i][0]) < 1e-1
        assert abs(modelled_curve.partial_fluxes[i][1] - validation_fluxes[i][1]) < 1e-1
        assert (
            abs(modelled_curve.feed_compositions[i].first - feed_compositions[i].first)
            < 1e-5
        )


def test_get_permeances(meoh_mtbe_pervaporation):
    feed_compositions = [
        Composition(p=(0.15 - i / 100), type=CompositionType("weight"))
        for i in range(15)
    ]
    modelled_curve = meoh_mtbe_pervaporation.get_ideal_diffusion_curve(
        feed_temperature=325.45,
        compositions=feed_compositions,
    )
    validation_permeances = (0.06093311, 0.002390)
    validation_selectivity = 25.4925
    modelled_permeances = modelled_curve.get_permeances
    calculated_selectivity = modelled_curve.get_selectivity
    for i in range(15):
        assert abs(modelled_permeances[i][0].value - validation_permeances[0]) < 1e-5
        assert abs(modelled_permeances[i][1].value - validation_permeances[1]) < 1e-5
        assert abs(calculated_selectivity[i] - validation_selectivity) < 1e-3
