import numpy
from pytest import fixture

from diffusion_curve import DiffusionCurve
from mixtures import Composition, CompositionType, Mixtures
from optimizer import Measurements, PervaporationFunction, find_best_fit


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


def test_find_best_fit_spi(spi_255_diffusion_curve):
    measurements_h2o = Measurements.from_diffusion_curve_first(spi_255_diffusion_curve)

    fit_h2o = find_best_fit(measurements_h2o, n=9)

    fit_h2o.save("tests/temp/test_pervaporation_function.pv")
    _fit_h2o = PervaporationFunction.load("tests/temp/test_pervaporation_function.pv")

    assert (numpy.array(fit_h2o.a) == numpy.array(_fit_h2o.a)).mean() == 1
    assert (numpy.array(fit_h2o.b) == numpy.array(_fit_h2o.b)).mean() == 1
    assert fit_h2o.alpha == _fit_h2o.alpha
