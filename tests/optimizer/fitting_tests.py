# plot SPI-255, Romakon-PM102, Pervap 4100
from pytest import fixture
from membrane import Membrane
from mixtures import Mixtures, Composition, CompositionType
from permeance import Permeance, Units
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurveSet, DiffusionCurve
from experiments import IdealExperiment, IdealExperiments
from components import Components


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


def test_fit(romakon_102_diffusion_curve_set):
    assert 0 == 0


def test_find_best_fit(romakon_102_diffusion_curve_set):
    assert 0 == 0
