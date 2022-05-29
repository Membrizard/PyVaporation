import typing

import attr
import numpy

from conditions import Conditions
from pervaporation import Pervaporation
from diffusion_curve import DiffusionCurve
from membrane import Membrane
from mixture import Composition, CompositionType, Mixture, get_nrtl_partial_pressures
from process import ProcessModel


def model_ideal_isothermal_process(
    pervaporation: Pervaporation,
    number_of_steps: int,
    delta_hours: float,
    conditions: Conditions,
    precision: typing.Optional[float] = 5e-5,
) -> ProcessModel:
    time: typing.List[float] = [delta_hours * step for step in range(number_of_steps)]

    partial_fluxes: typing.List[typing.Tuple[float, float]] = [(0.0, 0.0)] * (
        number_of_steps + 1
    )

    first_component_permeance = pervaporation.membrane.get_permeance(
        conditions.feed_temperature, pervaporation.mixture.first_component
    )
    second_component_permeance = pervaporation.membrane.get_permeance(
        conditions.feed_temperature, pervaporation.mixture.second_component
    )
    permeances: typing.List[typing.Tuple[float, float]] = [
        (first_component_permeance, second_component_permeance)
    ] * (number_of_steps + 1)

    permeate_composition: typing.List[Composition] = []
    feed_composition: typing.List[Composition] = [conditions.initial_feed_composition]

    feed_evaporation_heat: typing.List[float] = []
    permeate_condensation_heat: typing.List[float] = []
    feed_mass: typing.List[float] = [conditions.feed_amount]

    evaporation_heat_1 = (
        pervaporation.mixture.first_component.get_vaporisation_heat(
            conditions.feed_temperature
        )
        / pervaporation.mixture.first_component.molecular_weight
        * 1000
    )
    evaporation_heat_2 = (
        pervaporation.mixture.second_component.get_vaporisation_heat(
            conditions.feed_temperature
        )
        / pervaporation.mixture.first_component.molecular_weight
        * 1000
    )
    condensation_heat_1 = (
        pervaporation.mixture.first_component.get_vaporisation_heat(
            conditions.permeate_temperature
        )
        / pervaporation.mixture.first_component.molecular_weight
        * 1000
    )
    condensation_heat_2 = (
        pervaporation.mixture.first_component.get_vaporisation_heat(
            conditions.permeate_temperature
        )
        / pervaporation.mixture.first_component.molecular_weight
        * 1000
    )

    cooling_heat_1 = pervaporation.mixture.first_component.get_cooling_heat(
        conditions.permeate_temperature, conditions.feed_temperature
    )
    cooling_heat_2 = pervaporation.mixture.second_component.get_cooling_heat(
        conditions.permeate_temperature, conditions.feed_temperature
    )

    for step in range(len(time)):
        partial_fluxes[step] = pervaporation.calculate_partial_fluxes(
            conditions.feed_temperature,
            feed_composition[step],
            precision,
            conditions.permeate_temperature,
            first_component_permeance,
            second_component_permeance,
        )

        permeate_composition.append(
            Composition(
                p=partial_fluxes[step][0] / (sum(partial_fluxes[step])),
                type=CompositionType("weight"),
            )
        )

        d_mass_1 = partial_fluxes[step][0] * conditions.membrane_area * delta_hours
        d_mass_2 = partial_fluxes[step][1] * conditions.membrane_area * delta_hours

        feed_evaporation_heat.append(
            evaporation_heat_1 * d_mass_1 + evaporation_heat_2 * d_mass_2
        )
        permeate_condensation_heat.append(
            condensation_heat_1 * d_mass_1
            + condensation_heat_2 * d_mass_2
            + (cooling_heat_1 * d_mass_1 + cooling_heat_2 * d_mass_2)
            * (conditions.feed_temperature - conditions.permeate_temperature)
        )

        feed_mass.append(feed_mass[step] - d_mass_1 - d_mass_2)

        feed_composition.append(
            Composition(
                p=(feed_composition[step].p * feed_mass[step] - d_mass_1)
                / feed_mass[step + 1],
                type=CompositionType("weight"),
            )
        )

        permeances[step] = (
            first_component_permeance,
            second_component_permeance,
        )

    return ProcessModel(
        mixture=pervaporation.mixture,
        membrane_name=pervaporation.membrane.name,
        isothermal=True,
        feed_temperature=[conditions.feed_temperature] * number_of_steps,
        feed_composition=feed_composition,
        permeate_composition=permeate_composition,
        permeate_temperature=[conditions.permeate_temperature] * number_of_steps,
        feed_mass=feed_mass,
        partial_fluxes=partial_fluxes,
        permeances=permeances,
        time=time,
        feed_evaporation_heat=feed_evaporation_heat,
        permeate_condensation_heat=permeate_condensation_heat,
        initial_conditions=conditions,
        IsTimeDefined=True,
    )
