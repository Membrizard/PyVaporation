from pytest import fixture
from optimizer import Measurements, fit, find_best_fit
import matplotlib.pyplot as plt
from components import Components
from conditions import Conditions
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixtures import Mixtures, Composition, CompositionType, get_nrtl_partial_pressures
from permeance import Permeance, Units
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet, DiffusionCurve


@fixture
def experimental_50():
    compositions_50 = [
        0.02408580572,
        0.2326338341,
        0.2532445626,
        0.2754824538,
        0.4110793513,
        0.4284357542,
        0.4517584206,
        0.4745386994,
    ][::-1]

    permeances_50_H2O = [
        2688.778689,
        5623.966085,
        5790.19941,
        6078.245449,
        7525.770499,
        7823.818024,
        7525.770499,
        7525.770499,
    ][::-1]
    permeances_50_EtOH = [
        23.76315076,
        107.036062,
        120.2633877,
        132.5265399,
        207.1487026,
        215.3525347,
        205.1470494,
        205.1470494,
    ][::-1]
    return DiffusionCurve(
        mixture=Mixtures.H2O_EtOH,
        membrane_name="Pervap 4100",
        feed_temperature=368.15,
        feed_compositions=[
            Composition(p=composition, type=CompositionType.weight)
            for composition in compositions_50
        ],
        permeances=[
            (
                Permeance(value=permeances_50_H2O[i], units=Units.GPU),
                Permeance(value=permeances_50_EtOH[i], units=Units.GPU),
            )
            for i in range(len(compositions_50))
        ],
    )


@fixture
def pervaporation(experimental_50):

    return Pervaporation(
        Membrane(
            name=experimental_50.membrane_name,
            diffusion_curve_sets=[
                DiffusionCurveSet(
                    name_of_the_set="Various Initial Feed",
                    diffusion_curves=[experimental_50],
                )
            ],
        ),
        Mixtures.H2O_EtOH,
    )


def test_pervap_4100_50(pervaporation, experimental_50):
    modelled_curve = pervaporation.non_ideal_diffusion_curve(
        diffusion_curves=pervaporation.membrane.diffusion_curve_sets[0],
        feed_temperature=368.15,
        initial_feed_composition=Composition(p=0.4745, type=CompositionType.weight),
        delta_composition=-0.005,
        number_of_steps=90
    )
    for i in range(len(experimental_50.feed_compositions)):
        d = 1
        index = 0
        for k in range(len(modelled_curve.feed_compositions)):
            if d > abs(modelled_curve.feed_compositions[k].first-experimental_50.feed_compositions[i].first):
                d = abs(modelled_curve.feed_compositions[k].first-experimental_50.feed_compositions[i].first)
                index = k

        assert abs(modelled_curve.permeances[index][0].value-experimental_50.permeances[i][0].value) < 6e-3
        assert abs(modelled_curve.permeances[index][1].value-experimental_50.permeances[i][1].value) < 4.2e-4
        assert abs(modelled_curve.partial_fluxes[index][0]-experimental_50.partial_fluxes[i][0]) < 4e-1
        assert abs(modelled_curve.partial_fluxes[index][1]-experimental_50.partial_fluxes[i][1]) < 7e-2


def test_pervap_4100_25(pervaporation):
    assert 0 == 0


def test_pervap_4100_15(pervaporation):
    assert 0 == 0


def test_pervap_4100_6(pervaporation):
    assert 0 == 0
