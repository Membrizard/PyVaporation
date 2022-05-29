from pytest import fixture

from component import AllComponents
from conditions import Conditions
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixture import (
    AllMixtures,
    Composition,
    CompositionType,
    get_nrtl_partial_pressures,
)
from pervaporation import Pervaporation


@fixture
def all_components():
    return AllComponents.load("components.yml")


@fixture
def all_mixtures(all_components):
    return AllMixtures.load("mixtures.yml", all_components)


@fixture
def pervap_4101(all_components, all_mixtures):

    experiment_h2o_1 = IdealExperiment(
        name="Pervap 4101",
        temperature=333.15,
        component=all_components.h2o,
        permeance=0.01491,
    )

    experiment_h2o_2 = IdealExperiment(
        name="Pervap 4101",
        temperature=368.15,
        component=all_components.h2o,
        permeance=0.02460,
    )

    experiment_h2o_3 = IdealExperiment(
        name="Pervap 4101",
        temperature=378.15,
        component=all_components.h2o,
        permeance=0.03555,
    )

    experiment_etoh_1 = IdealExperiment(
        name="Pervap 4101",
        temperature=368.15,
        component=all_components.etoh,
        permeance=0.0000098,
    )

    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_h2o_1,
            experiment_h2o_2,
            experiment_h2o_3,
            experiment_etoh_1,
        ]
    )

    return Membrane(ideal_experiments=ideal_experiments, name="Pervap 4101")


@fixture
def h2o_etoh_pervaporation(pervap_4101, all_mixtures):
    return Pervaporation(membrane=pervap_4101, mixture=all_mixtures.h2o_etoh)


def test_validate_against_experimet(
    pervap_4101, meoh_mtbe_pervaporation, all_mixtures, all_components
):

    feed = [
        0.0080,
        0.0105,
        0.0193,
        0.0247,
        0.0315,
        0.0585,
        0.0728,
        0.0869,
        0.1000,
        0.1341,
        0.1536,
    ]

    feed_compositions = [
        Composition(p=feed[i], type=CompositionType("weight")) for i in range(len(feed))
    ]

    validation_fluxes = [
        (0.057, 0.010),
        (0.070, 0.018),
        (0.159, 0.016),
        (0.214, 0.018),
        (0.280, 0.018),
        (0.608, 0.018),
        (0.799, 0.027),
        (0.970, 0.027),
        (1.135, 0.025),
        (1.508, 0.033),
        (1.697, 0.033),
    ]

    for i in range(len(feed_compositions)):
        partial_fluxes = meoh_mtbe_pervaporation.calculate_partial_fluxes(
            feed_temperature=368.15, composition=feed_compositions[i]
        )
        assert abs(partial_fluxes[0] - validation_fluxes[i][0]) < 0.54
        assert abs(partial_fluxes[1] - validation_fluxes[i][1]) < 0.035
