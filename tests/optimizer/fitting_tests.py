# plot SPI-255, Romakon-PM102, Pervap 4100
from pytest import fixture
from membrane import Membrane
from mixtures import Mixtures, Composition, CompositionType
from permeance import Permeance, Units
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet, DiffusionCurve
from experiments import IdealExperiment, IdealExperiments
from components import Components
from optimizer import Measurements, find_best_fit


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

    return Membrane(
        name="Romakon-PM102",
        diffusion_curve_sets=[
            DiffusionCurveSet(
                name_of_the_set="water/acetic acid", diffusion_curves=diffusion_curves
            )
        ],
    )


@fixture
def spi_255_diffusion_curve():
    compositions = [0.9999, 0.883, 0.756, 0.658, 0.490, 0.430, 0.401, 0.280, 0.240, 0.044]

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
    assert 0 == 0


def test_find_best_fit(romakon_102_diffusion_curve_set):
    # Best fit for H2O
    # PervaporationFunction(n=1, m=2, alpha=0.04549912188685387, a=array([0.57744257]), b=array([ 463.53220199,   24.18448749, -247.65949047]))
    # Best fit for Acetic Acid
    # PervaporationFunction(n=0, m=2, alpha=0.0005275351075504815, a=array([], dtype=float64), b=array([ 10097.8817611 , -15213.08201821,   3999.51281411]))

    assert 0 == 0


def test_find_best_fit_spi(spi_255_diffusion_curve):
    measurements_h2o = Measurements.from_diffusion_curve_first(spi_255_diffusion_curve)
    measurements_etoh = Measurements.from_diffusion_curve_second(spi_255_diffusion_curve)

    fit_h2o = find_best_fit(measurements_h2o, n=9)
    fit_etoh = find_best_fit(measurements_etoh)

    for i in range(len(measurements_h2o)):
        assert abs(fit_h2o(measurements_h2o[i].x, 313.15)-measurements_h2o[i].p) < 1.3e-2
        assert abs(fit_etoh(measurements_etoh[i].x, 313.15)-measurements_etoh[i].p) < 1.6e-3

