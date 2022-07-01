import numpy
from pytest import fixture

from diffusion_curve import DiffusionCurve
from mixtures import Composition, CompositionType, Mixtures
from optimizer import Measurements, find_best_fit, fit
from permeance import Permeance, Units


@fixture
def romakon_102_diffusion_curve_set():
    compositions = [0.99, 0.89, 0.72]

    permeance_h2o_50 = [633e-9, 416e-9, 364e-9]
    permeance_h2o_60 = [608e-9, 463e-9, 388e-9]
    permeance_h2o_70 = [629e-9, 456e-9, 413e-9]
    permeance_h2o_80 = [630e-9, 480e-9, 426e-9]

    permeance_acetic_acid_50 = [58.4e-9, 3.75e-9, 1.75e-9]
    permeance_acetic_acid_60 = [55.4e-9, 3.09e-9, 3.63e-9]
    permeance_acetic_acid_70 = [64.7e-9, 5.3e-9, 3.98e-9]
    permeance_acetic_acid_80 = [37e-9, 10.3e-9, 4.3e-9]

    permeances_h2o = [
        permeance_h2o_50,
        permeance_h2o_60,
        permeance_h2o_70,
        permeance_h2o_80,
    ]
    permeances_acetic_acid = [
        permeance_acetic_acid_50,
        permeance_acetic_acid_60,
        permeance_acetic_acid_70,
        permeance_acetic_acid_80,
    ]

    temperatures = [323.15, 333.15, 343.15, 353.15]
    diffusion_curves = []

    for t in range(len(temperatures)):
        diffusion_curves.append(
            DiffusionCurve(
                mixture=Mixtures.H2O_AceticAcid,
                membrane_name="Romakon-PM102",
                feed_temperature=temperatures[t],
                feed_compositions=[
                    Composition(p=c, type=CompositionType.weight) for c in compositions
                ],
                permeances=[
                    (
                        Permeance(value=permeances_h2o[t][i], units=Units.SI),
                        Permeance(value=permeances_acetic_acid[t][i], units=Units.SI),
                    )
                    for i in range(len(permeances_h2o[t]))
                ],
            )
        )

    return diffusion_curves


@fixture
def spi_255_diffusion_curve():
    compositions = [
        0.9999,
        0.883,
        0.756,
        0.658,
        0.490,
        0.430,
        0.401,
        0.280,
        0.240,
        0.044,
    ]

    flux_h2o_40 = [
        0.7848,
        0.3480,
        0.2566,
        0.4406,
        0.3340,
        0.2773,
        0.2816,
        0.2292,
        0.1066,
        0.0247,
    ]

    flux_etoh_40 = [
        0,
        0.0017,
        0.0013,
        0.0136,
        0.0684,
        0.0568,
        0.0577,
        0.0573,
        0.0457,
        0.0165,
    ]

    return DiffusionCurve(
        mixture=Mixtures.H2O_EtOH,
        membrane_name="SPI 255 dense",
        feed_temperature=313.15,
        feed_compositions=[
            Composition(p=c, type=CompositionType.weight) for c in compositions
        ],
        partial_fluxes=[
            (flux_h2o_40[i], flux_etoh_40[i]) for i in range(len(compositions))
        ],
    )


def test_fit(romakon_102_diffusion_curve_set):
    measurements_h2o = Measurements.from_diffusion_curves_first(
        romakon_102_diffusion_curve_set
    )

    fit_1 = fit(measurements_h2o, n=1, m=0)
    fit_2 = fit(measurements_h2o, n=0, m=1)
    x = numpy.arange(0, 1, 0.05)

    for i in range(len(x)):
        assert abs(fit_1(x[i], 333.15) - fit_2(x[i], 333.15)) < 5e-4


def test_find_best_fit(romakon_102_diffusion_curve_set):
    measurements_h2o = Measurements.from_diffusion_curves_first(
        romakon_102_diffusion_curve_set
    )

    fit_h2o = find_best_fit(measurements_h2o)

    validation_b_h2o = [463.5322019863124, 24.18448749164457, -247.6594904722613]

    assert fit_h2o.n == 1
    assert fit_h2o.m == 2
    assert round(fit_h2o.alpha, 4) == round(0.04549912188685387, 4)
    assert round(fit_h2o.a[0], 4) == round(0.5774425746265356, 4)
    for i in range(len(validation_b_h2o)):
        assert abs(fit_h2o.b[i] - validation_b_h2o[i]) < 1e-5


def test_find_best_fit_spi(spi_255_diffusion_curve):
    measurements_h2o = Measurements.from_diffusion_curve_first(spi_255_diffusion_curve)
    measurements_etoh = Measurements.from_diffusion_curve_second(
        spi_255_diffusion_curve
    )

    fit_h2o = find_best_fit(measurements_h2o, n=9)
    fit_etoh = find_best_fit(measurements_etoh)

    for i in range(len(measurements_h2o)):
        assert (
            abs(fit_h2o(measurements_h2o[i].x, 313.15) - measurements_h2o[i].p) < 1.3e-2
        )
        assert (
            abs(fit_etoh(measurements_etoh[i].x, 313.15) - measurements_etoh[i].p)
            < 1.6e-3
        )
    assert fit_h2o.alpha == 0.34982213695528863

    validation_a = [
        5.363743190543987,
        -0.7698301270286048,
        -3.890383383394793,
        -2.1531077482555547,
        0.29257926263632544,
        0.7694100222167437,
        0.966517366730522,
        1.0129579184563662,
        0.9807091886465034,
    ]

    for i in range(len(fit_h2o.a)):
        assert fit_h2o.a[i] == validation_a[i]

    assert fit_h2o.b[0] == 1177.8598832639543
