from pathlib import Path

from pytest import fixture

from pyvaporation.components import Components
from pyvaporation.conditions import Conditions
from pyvaporation.diffusion_curve import DiffusionCurve, DiffusionCurveSet
from pyvaporation.experiments import IdealExperiment, IdealExperiments
from pyvaporation.membrane import Membrane
from pyvaporation.mixtures import Composition, CompositionType, Mixtures
from pyvaporation.permeance import Permeance
from pyvaporation.pervaporation import Pervaporation


@fixture
def ideal_pervaporation():
    experiment_h2o_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=Components.H2O,
        permeance=Permeance(0.036091),
        activation_energy=19944,
    )

    experiment_etoh_1 = IdealExperiment(
        name="Romakon-PM102",
        temperature=323.15,
        component=Components.EtOH,
        permeance=Permeance(0.0000282),
        activation_energy=110806,
    )

    ideal_experiments = IdealExperiments(
        experiments=[
            experiment_h2o_1,
            experiment_etoh_1,
        ]
    )

    membrane = Membrane(ideal_experiments=ideal_experiments, name="Romakon-PM102")
    diffusion_curve = Pervaporation(
        membrane=membrane,
        mixture=Mixtures.H2O_EtOH,
    ).ideal_diffusion_curve(
        compositions=[
            Composition(p=i / 10, type=CompositionType.weight) for i in range(1, 10)
        ],
        feed_temperature=323.15,
    )

    membrane = Membrane(
        ideal_experiments=ideal_experiments,
        diffusion_curve_sets=[DiffusionCurveSet("ideal curve", [diffusion_curve])],
        name="Romakon-PM102",
    )

    return Pervaporation(membrane, Mixtures.H2O_EtOH)


@fixture()
def ideal_conditions():
    return Conditions(
        membrane_area=0.4155,
        initial_feed_temperature=333.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType.weight),
    )


@fixture
def romakon_al2_pervaporation():
    compositions = [0.040, 0.037, 0.030, 0.026, 0.024, 0.021]

    flux_h2o = [
        0.062062,
        0.062062,
        0.03988709016,
        0.03720394737,
        0.004707604167,
        0.005613600498,
    ]

    flux_etoh = [
        0.006138,
        0.006138,
        0.003944877049,
        0.00413377193,
        0.0007663541667,
        0.0009520141196,
    ]

    curve = DiffusionCurve(
        mixture=Mixtures.H2O_EtOH,
        membrane_name="Romakon-Al2",
        feed_temperature=319.65,
        feed_compositions=[
            Composition(c, CompositionType.weight) for c in compositions
        ],
        partial_fluxes=[(flux_h2o[i], flux_etoh[i]) for i in range(len(flux_h2o))],
    )

    curve_set = DiffusionCurveSet(name=" ", diffusion_curves=[curve])

    romakon_al2 = Membrane(name="Romakon Al2", diffusion_curve_sets=[curve_set])

    return Pervaporation(romakon_al2, Mixtures.H2O_EtOH)


def test_validate_against_ideal_process(ideal_pervaporation, ideal_conditions):
    number_of_steps = 200
    delta_hours = 0.125
    ideal_model = ideal_pervaporation.ideal_isothermal_process(
        conditions=ideal_conditions,
        number_of_steps=number_of_steps,
        delta_hours=delta_hours,
    )
    non_ideal_model = ideal_pervaporation.non_ideal_isothermal_process(
        conditions=ideal_conditions,
        diffusion_curve_set=ideal_pervaporation.membrane.diffusion_curve_sets[0],
        number_of_steps=number_of_steps,
        delta_hours=delta_hours,
    )

    for i in range(number_of_steps):
        assert (
            abs(non_ideal_model.partial_fluxes[i][0] - ideal_model.partial_fluxes[i][0])
            < 9e-3
        )
        assert (
            abs(non_ideal_model.partial_fluxes[i][1] - ideal_model.partial_fluxes[i][1])
            < 9e-3
        )
        assert (
            abs(
                non_ideal_model.permeances[i][0].value
                - ideal_model.permeances[i][0].value
            )
            < 5e-4
        )
        assert (
            abs(
                non_ideal_model.permeances[i][1].value
                - ideal_model.permeances[i][1].value
            )
            < 5e-4
        )
        assert abs(non_ideal_model.feed_mass[i] - ideal_model.feed_mass[i]) < 8.5e-2
        assert (
            abs(
                non_ideal_model.feed_evaporation_heat[i]
                - ideal_model.feed_evaporation_heat[i]
            )
            < 1.1
        )
        assert (
            non_ideal_model.permeate_condensation_heat[i]
            == ideal_model.permeate_condensation_heat[i]
        )
        assert (
            abs(
                non_ideal_model.feed_compositions[i].first
                - ideal_model.feed_compositions[i].first
            )
            < 1e-2
        )


def test_real_process_modelling_romakon_al2(romakon_al2_pervaporation):
    experimental_time = [0, 1.03, 2.02, 2.97, 4.00, 5.02]

    conditions = Conditions(
        membrane_area=0.0048,
        initial_feed_composition=Composition(p=0.04, type=CompositionType.weight),
        initial_feed_amount=0.047,
        initial_feed_temperature=319.65,
    )

    modelled_experiment = romakon_al2_pervaporation.non_ideal_isothermal_process(
        conditions=conditions,
        diffusion_curve_set=romakon_al2_pervaporation.membrane.diffusion_curve_sets[0],
        number_of_steps=50,
        delta_hours=0.1,
    )

    curve = romakon_al2_pervaporation.membrane.diffusion_curve_sets[0].diffusion_curves[
        0
    ]

    for i in range(len(experimental_time)):
        d = 1
        index = 0
        for k in range(len(modelled_experiment.time)):
            if d > abs(modelled_experiment.time[k] - experimental_time[i]):
                d = abs(modelled_experiment.time[k] - experimental_time[i])
                index = k
        assert (
            abs(
                curve.feed_compositions[i].first
                - modelled_experiment.feed_compositions[index].first
            )
            < 3e-3
        )


def test_validate_against_literature_data_pervap_4101():
    """
    :return: Algorithms are validated against experimental default_membranes provided in doi:10.3390/membranes8010004
    """

    membrane = Membrane.load(Path("tests/default_membranes/Pervap_4101"))

    pv = Pervaporation(
        membrane=membrane,
        mixture=Mixtures.H2O_EtOH,
    )

    con = Conditions(
        membrane_area=0.017,
        initial_feed_temperature=368.15,
        initial_feed_amount=1.5,
        initial_feed_composition=Composition(p=0.1, type=CompositionType.weight),
        permeate_pressure=0,
    )

    model = pv.non_ideal_isothermal_process(
        conditions=con,
        diffusion_curve_set=membrane.diffusion_curve_sets[0],
        initial_permeances=(
            Permeance(0.0153),
            Permeance(0.00000632),
        ),
        number_of_steps=50,
        delta_hours=0.2,
    )

    experiment_time_hours = [
        0,
        0.5163,
        1.0148,
        1.5133,
        2.0118,
        3.0087,
        4.0057,
        5.0173,
        6.0143,
        7.0112,
        8.0082,
        9.0052,
        10,
    ]

    experiment_water_fraction = [
        0.1,
        0.0946,
        0.0938,
        0.0907,
        0.0887,
        0.0842,
        0.0791,
        0.0730,
        0.0702,
        0.0669,
        0.0625,
        0.0588,
        0.0573,
    ]

    experiment_water_flux = [
        0.4863,
        0.5209,
        0.5301,
        0.5865,
        0.4517,
        0.4736,
        0.4353,
        0.3789,
        0.3770,
        0.3370,
        0.2987,
        0.2823,
        0.2823,
    ]

    for i in range(len(experiment_time_hours)):
        d = 1
        index = 0
        for k in range(len(model.time)):
            if d > abs(model.time[k] - experiment_time_hours[i]):
                d = abs(model.time[k] - experiment_time_hours[i])
                index = k

        assert (
            abs(experiment_water_fraction[i] - model.feed_compositions[index].first)
            < 2e-3
        )

        assert abs(experiment_water_flux[i] - model.partial_fluxes[index][0]) < 1.1e-1
