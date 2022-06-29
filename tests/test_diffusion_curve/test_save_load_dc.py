from pytest import fixture
from pathlib import Path

from components import Components
from diffusion_curve import DiffusionCurve, DiffusionCurveSet
from experiments import IdealExperiment, IdealExperiments
from membrane import Membrane
from mixtures import Composition, CompositionType, Mixtures
from permeance import Permeance, Units
from pervaporation import Pervaporation


@fixture
def romakon_102_single_diffusion_curve():

    compositions = [0.99, 0.89, 0.72]
    permeances_h2o = [629e-9, 456e-9, 413e-9]
    permeances_acetic_acid = [64.7e-9, 5.3e-9, 3.98e-9]

    curve = DiffusionCurve(
        mixture=Mixtures.H2O_AceticAcid,
        membrane_name="Romakon-PM102",
        feed_temperature=343.15,
        feed_compositions=[
            Composition(p=c, type=CompositionType.weight) for c in compositions
        ],
        permeances=[
            (
                Permeance(value=permeances_h2o[i], units=Units.SI),
                Permeance(value=permeances_acetic_acid[i], units=Units.SI),
            )
            for i in range(len(permeances_h2o))
        ],
    )

    experiment_h2o_50 = IdealExperiment(
        name="H2O",
        temperature=323.15,
        component=Components.H2O,
        permeance=Permeance(value=416e-9, units=Units.SI),
    )

    experiment_h2o_60 = IdealExperiment(
        name="H2O",
        temperature=333.15,
        component=Components.H2O,
        permeance=Permeance(value=463e-9, units=Units.SI),
    )

    experiment_h2o_80 = IdealExperiment(
        name="H2O",
        temperature=353.15,
        component=Components.H2O,
        permeance=Permeance(value=480e-9, units=Units.SI),
    )

    experiment_acetic_acid_50 = IdealExperiment(
        name="Acetic acid",
        temperature=323.15,
        component=Components.AceticAcid,
        permeance=Permeance(value=3.75e-9, units=Units.SI),
    )

    experiment_acetic_acid_60 = IdealExperiment(
        name="Acetic acid",
        temperature=333.15,
        component=Components.AceticAcid,
        permeance=Permeance(value=3.09e-9, units=Units.SI),
    )

    experiment_acetic_acid_80 = IdealExperiment(
        name="Acetic acid",
        temperature=353.15,
        component=Components.AceticAcid,
        permeance=Permeance(value=10.3e-9, units=Units.SI),
    )

    experiments = IdealExperiments(
        experiments=[
            experiment_h2o_50,
            experiment_h2o_60,
            experiment_h2o_80,
            experiment_acetic_acid_50,
            experiment_acetic_acid_60,
            experiment_acetic_acid_80,
        ]
    )

    return Membrane(
        name="Romakon-PM102",
        ideal_experiments=experiments,
        diffusion_curve_sets=[
            DiffusionCurveSet(name="water/acetic acid", diffusion_curves=[curve])
        ],
    )


def test_semi_ideal_curve(romakon_102_single_diffusion_curve):
    pervaporation = Pervaporation(
        membrane=romakon_102_single_diffusion_curve, mixture=Mixtures.H2O_AceticAcid
    )

    compositions_100 = [0.89, 0.72]

    permeances_100_h2o = [491e-9, 510e-9]
    permeances_100_acetic_acid = [14.2e-9, 4.9e-9]

    experimental_100: DiffusionCurve = DiffusionCurve(
        mixture=Mixtures.H2O_AceticAcid,
        membrane_name="Romakon-PM102",
        feed_temperature=373.15,
        feed_compositions=[
            Composition(p=composition, type=CompositionType.weight)
            for composition in compositions_100
        ],
        permeances=[
            (
                Permeance(value=permeances_100_h2o[i], units=Units.SI),
                Permeance(value=permeances_100_acetic_acid[i], units=Units.SI),
            )
            for i in range(len(compositions_100))
        ],
    )

    modelled_curve = pervaporation.non_ideal_diffusion_curve(
        diffusion_curve_set=romakon_102_single_diffusion_curve.diffusion_curve_sets[0],
        feed_temperature=373.15,
        initial_feed_composition=Composition(p=0.99, type=CompositionType.weight),
        delta_composition=-0.0054,
        number_of_steps=50,
    )

    modelled_curve.to_csv("tests/temp/test_diffusion_curve.csv")
    modelled_curve = DiffusionCurveSet.load(Path("tests/temp/test_diffusion_curve.csv")).diffusion_curves[
        0
    ]

    for i in range(len(experimental_100.feed_compositions)):
        d = 1
        index = 0
        for k in range(len(modelled_curve.feed_compositions)):
            if d > abs(
                modelled_curve.feed_compositions[k].first
                - experimental_100.feed_compositions[i].first
            ):
                d = abs(
                    modelled_curve.feed_compositions[k].first
                    - experimental_100.feed_compositions[i].first
                )
                index = k

        assert (
            abs(
                modelled_curve.permeances[index][0].value
                - experimental_100.permeances[i][0].value
            )
            < 6.5e-3
        )
        assert (
            abs(
                modelled_curve.permeances[index][1].value
                - experimental_100.permeances[i][1].value
            )
            < 2.1e-3
        )
        assert (
            abs(
                modelled_curve.partial_fluxes[index][0]
                - experimental_100.partial_fluxes[i][0]
            )
            < 6e-1
        )
        assert (
            abs(
                modelled_curve.partial_fluxes[index][1]
                - experimental_100.partial_fluxes[i][1]
            )
            < 7e-3
        )
