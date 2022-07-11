from pytest import fixture

from pyvaporation.components import Components
from pyvaporation.conditions import Conditions, TemperatureProgram
from pyvaporation.diffusion_curve import DiffusionCurveSet
from pyvaporation.experiments import IdealExperiment, IdealExperiments
from pyvaporation.membrane import Membrane
from pyvaporation.mixtures import Composition, CompositionType, Mixtures
from pyvaporation.permeance import Permeance
from pyvaporation.pervaporation import Pervaporation
from pathlib import Path


@fixture
def pervaporation():
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
def basic_conditions():
    return Conditions(
        membrane_area=0.4155,
        initial_feed_temperature=333.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType.weight),
    )


def test_validate_against_ideal_process_self_cooling(pervaporation, basic_conditions):
    number_of_steps = 200
    delta_hours = 0.125
    ideal_model = pervaporation.ideal_non_isothermal_process(
        conditions=basic_conditions,
        number_of_steps=number_of_steps,
        delta_hours=delta_hours,
    )
    non_ideal_model = pervaporation.non_ideal_non_isothermal_process(
        conditions=basic_conditions,
        diffusion_curve_set=pervaporation.membrane.diffusion_curve_sets[0],
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
        assert abs(non_ideal_model.feed_mass[i] - ideal_model.feed_mass[i]) < 6.5e-2
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


def test_validate_against_ideal_process_temperature_program(pervaporation):
    temperature_program = Conditions(
        membrane_area=0.4155,
        initial_feed_temperature=333.15,
        initial_feed_amount=12,
        initial_feed_composition=Composition(p=0.94, type=CompositionType.weight),
        temperature_program=TemperatureProgram([333.15, -1]),
    )
    number_of_steps = 200
    delta_hours = 0.125
    ideal_model = pervaporation.ideal_non_isothermal_process(
        conditions=temperature_program,
        number_of_steps=number_of_steps,
        delta_hours=delta_hours,
    )
    non_ideal_model = pervaporation.non_ideal_non_isothermal_process(
        conditions=temperature_program,
        diffusion_curve_set=pervaporation.membrane.diffusion_curve_sets[0],
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
        assert abs(non_ideal_model.feed_mass[i] - ideal_model.feed_mass[i]) < 6.5e-2
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

def test_validate_against_experimental_data():
    """
    Data for validation is taken from: https://doi.org/10.1007/bf02705302.
    """
    membrane = Membrane.load(Path("tests/default_membranes/Chang_et_al_1998"))

    pv = Pervaporation(
        membrane=membrane,
        mixture=Mixtures.H2O_EtOH,
    )

    model_temp = []
    model_area = []
    model_fraction = [0.062]
    for i in range(4):
        conditions = Conditions(
            membrane_area=1,
            initial_feed_temperature=368.475,
            initial_feed_amount=12.106,
            initial_feed_composition=Composition(p=model_fraction[i], type=CompositionType.weight),
            permeate_pressure=1.3,
        )

        model = pv.non_ideal_non_isothermal_process(
            conditions=conditions,
            diffusion_curve_set=membrane.diffusion_curve_sets[0],
            number_of_steps=11,
            delta_hours=0.1,
        )

        model_temp.append(model.feed_temperature[-1]-273.15)
        model_area.append(i)
        model_fraction.append(model.feed_compositions[-1].first)

    exp_ret_temp = [75.59665724, 78.66904264, 84.62947032, 92.67912007]

    exp_fraction = [0.06162740899,
                    0.04331905782,
                    0.02886509636,
                    0.01794432548,
                    0.01141327623]

    for i in range(len(exp_ret_temp)):
        assert abs(exp_ret_temp[i]-model_temp[i]) < 5.5

    for i in range(len(exp_fraction)):
        assert abs(exp_fraction[i] - model_fraction[i]) < 1e-2