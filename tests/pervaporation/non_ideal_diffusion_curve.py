from pytest import fixture
from membrane import Membrane
from mixtures import Mixtures, Composition, CompositionType
from permeance import Permeance, Units
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet, DiffusionCurve


@fixture
def experimental_50():
    compositions_50 = [
        0.0240858057,
        0.2326338341,
        0.2532445626,
        0.2754824538,
        0.4110793513,
        0.4284357542,
        0.4517584206,
        0.4745386994,
    ]

    permeances_50_h2o = [
        2688.778689,
        5623.966085,
        5790.19941,
        6078.245449,
        7525.770499,
        7823.818024,
        7525.770499,
        7525.770499,
    ]
    permeances_50_etoh = [
        23.76315076,
        107.036062,
        120.2633877,
        132.5265399,
        207.1487026,
        215.3525347,
        205.1470494,
        205.1470494,
    ]
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
                Permeance(value=permeances_50_h2o[i], units=Units.GPU),
                Permeance(value=permeances_50_etoh[i], units=Units.GPU),
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
        number_of_steps=90,
    )
    for i in range(len(experimental_50.feed_compositions)):
        d = 1
        index = 0
        for k in range(len(modelled_curve.feed_compositions)):
            if d > abs(
                modelled_curve.feed_compositions[k].first
                - experimental_50.feed_compositions[i].first
            ):
                d = abs(
                    modelled_curve.feed_compositions[k].first
                    - experimental_50.feed_compositions[i].first
                )
                index = k

        assert (
            abs(
                modelled_curve.permeances[index][0].value
                - experimental_50.permeances[i][0].value
            )
            < 6e-3
        )
        assert (
            abs(
                modelled_curve.permeances[index][1].value
                - experimental_50.permeances[i][1].value
            )
            < 4.2e-4
        )
        assert (
            abs(
                modelled_curve.partial_fluxes[index][0]
                - experimental_50.partial_fluxes[i][0]
            )
            < 4e-1
        )
        assert (
            abs(
                modelled_curve.partial_fluxes[index][1]
                - experimental_50.partial_fluxes[i][1]
            )
            < 7e-2
        )


def test_pervap_4100_25(pervaporation):
    compositions_25 = [
        0.0053734337,
        0.0308656506,
        0.0376454955,
        0.1789374627,
        0.2239556327,
        0.2570412757,
    ]

    permeances_25_h2o = [
        1212.734594,
        1561.012113,
        1703.565186,
        3202.285314,
        3329.107307,
        3598.018023,
    ]
    permeances_25_etoh = [
        1.778198976,
        2.426177206,
        2.522262532,
        10.21010056,
        11.92617288,
        13.01528204,
    ]

    experimental_25 = DiffusionCurve(
        mixture=Mixtures.H2O_EtOH,
        membrane_name="Pervap 4100",
        feed_temperature=368.15,
        feed_compositions=[
            Composition(p=composition, type=CompositionType.weight)
            for composition in compositions_25
        ],
        permeances=[
            (
                Permeance(value=permeances_25_h2o[i], units=Units.GPU),
                Permeance(value=permeances_25_etoh[i], units=Units.GPU),
            )
            for i in range(len(compositions_25))
        ],
    )

    modelled_curve = pervaporation.non_ideal_diffusion_curve(
        diffusion_curves=pervaporation.membrane.diffusion_curve_sets[0],
        feed_temperature=368.15,
        initial_feed_composition=Composition(p=0.25705, type=CompositionType.weight),
        delta_composition=-0.005,
        number_of_steps=50,
        initial_permeances=(
            Permeance(value=3598.018023, units=Units.GPU),
            Permeance(value=13.01528204, units=Units.GPU),
        ),
    )

    for i in range(len(experimental_25.feed_compositions)):
        d = 1
        index = 0
        for k in range(len(modelled_curve.feed_compositions)):
            if d > abs(
                modelled_curve.feed_compositions[k].first
                - experimental_25.feed_compositions[i].first
            ):
                d = abs(
                    modelled_curve.feed_compositions[k].first
                    - experimental_25.feed_compositions[i].first
                )
                index = k

        assert (
            abs(
                modelled_curve.permeances[index][0].value
                - experimental_25.permeances[i][0].value
            )
            < 6.6e-3
        )
        assert (
            abs(
                modelled_curve.permeances[index][1].value
                - experimental_25.permeances[i][1].value
            )
            < 4.2e-4
        )
        assert (
            abs(
                modelled_curve.partial_fluxes[index][0]
                - experimental_25.partial_fluxes[i][0]
            )
            < 4e-1
        )
        assert (
            abs(
                modelled_curve.partial_fluxes[index][1]
                - experimental_25.partial_fluxes[i][1]
            )
            < 7e-2
        )


def test_pervap_4100_15(pervaporation):
    compositions_15 = [
        1.188208494,
        1.513641048,
        9.052828551,
        11.22237891,
        12.68682541
    ]

    permeances_15_h2o = [
        970.0097218,
        1018.264961,
        2171.612439,
        2171.612439,
        2586.349948
    ]
    permeances_15_etoh = [
        0.513111450,
        0.479394567,
        1.115750557,
        1.464338579,
        1.677562465
    ]

    experimental_15 = DiffusionCurve(
        mixture=Mixtures.H2O_EtOH,
        membrane_name="Pervap 4100",
        feed_temperature=368.15,
        feed_compositions=[
            Composition(p=composition/100, type=CompositionType.weight)
            for composition in compositions_15
        ],
        permeances=[
            (
                Permeance(value=permeances_15_h2o[i], units=Units.GPU),
                Permeance(value=permeances_15_etoh[i], units=Units.GPU),
            )
            for i in range(len(compositions_15))
        ],
    )

    modelled_curve = pervaporation.non_ideal_diffusion_curve(
        diffusion_curves=pervaporation.membrane.diffusion_curve_sets[0],
        feed_temperature=368.15,
        initial_feed_composition=Composition(p=0.12686, type=CompositionType.weight),
        delta_composition=-0.0025,
        number_of_steps=40,
        initial_permeances=(
            Permeance(value=2586.349948, units=Units.GPU),
            Permeance(value=1.677562465, units=Units.GPU),
        ),
    )

    for i in range(len(experimental_15.feed_compositions)):
        d = 1
        index = 0
        for k in range(len(modelled_curve.feed_compositions)):
            if d > abs(
                    modelled_curve.feed_compositions[k].first
                    - experimental_15.feed_compositions[i].first
            ):
                d = abs(
                    modelled_curve.feed_compositions[k].first
                    - experimental_15.feed_compositions[i].first
                )
                index = k

        assert (
                abs(
                    modelled_curve.permeances[index][0].value
                    - experimental_15.permeances[i][0].value
                )
                < 1.7e-2
        )
        assert (
                abs(
                    modelled_curve.permeances[index][1].value
                    - experimental_15.permeances[i][1].value
                )
                < 6.6e-6
        )
        assert (
                abs(
                    modelled_curve.partial_fluxes[index][0]
                    - experimental_15.partial_fluxes[i][0]
                )
                < 3.5e-1
        )
        assert (
                abs(
                    modelled_curve.partial_fluxes[index][1]
                    - experimental_15.partial_fluxes[i][1]
                )
                < 1.1e-3
        )


def test_pervap_4100_6(pervaporation):
    compositions_6 = [
        1.703476704,
        1.94755112,
        4.876444107,
        5.473070456,
        5.961219287
    ]

    permeances_6_h2o = [
        746.3089286,
        838.5364554,
        1038.232656,
        1122.096546,
        1133.045025

    ]
    permeances_6_etoh = [
        0.2100162934,
        0.2100162934,
        0.2141346138,
        0.2550303347,
        0.2100162934

    ]

    experimental_6 = DiffusionCurve(
        mixture=Mixtures.H2O_EtOH,
        membrane_name="Pervap 4100",
        feed_temperature=368.15,
        feed_compositions=[
            Composition(p=composition / 100, type=CompositionType.weight)
            for composition in compositions_6
        ],
        permeances=[
            (
                Permeance(value=permeances_6_h2o[i], units=Units.GPU),
                Permeance(value=permeances_6_etoh[i], units=Units.GPU),
            )
            for i in range(len(compositions_6))
        ],
    )

    modelled_curve = pervaporation.non_ideal_diffusion_curve(
        diffusion_curves=pervaporation.membrane.diffusion_curve_sets[0],
        feed_temperature=368.15,
        initial_feed_composition=Composition(p=0.05961, type=CompositionType.weight),
        delta_composition=-0.001,
        number_of_steps=40,
        initial_permeances=(
            Permeance(value=1133.045025, units=Units.GPU),
            Permeance(value=0.2100162934, units=Units.GPU),
        ),
    )

    for i in range(len(experimental_6.feed_compositions)):
        d = 1
        index = 0
        for k in range(len(modelled_curve.feed_compositions)):
            if d > abs(
                    modelled_curve.feed_compositions[k].first
                    - experimental_6.feed_compositions[i].first
            ):
                d = abs(
                    modelled_curve.feed_compositions[k].first
                    - experimental_6.feed_compositions[i].first
                )
                index = k

        assert (
                abs(
                    modelled_curve.permeances[index][0].value
                    - experimental_6.permeances[i][0].value
                )
                < 4.5e-3
        )
        assert (
                abs(
                    modelled_curve.permeances[index][1].value
                    - experimental_6.permeances[i][1].value
                )
                < 5e-6
        )
        assert (
                abs(
                    modelled_curve.partial_fluxes[index][0]
                    - experimental_6.partial_fluxes[i][0]
                )
                < 6e-2
        )
        assert (
                abs(
                    modelled_curve.partial_fluxes[index][1]
                    - experimental_6.partial_fluxes[i][1]
                )
                < 8e-4
        )

