import typing

import attr
import numpy

from diffusion import DiffusionCurve
from conditions import Conditions
from membrane import Membrane
from mixture import Composition, CompositionType, Mixture, get_nrtl_partial_pressures


def _get_permeate_composition_from_fluxes(
    fluxes: typing.Tuple[float, float],
) -> Composition:
    return Composition(
        p=fluxes[0] / sum(fluxes),
        type=CompositionType("weight"),
    )


@attr.s(auto_attribs=True)
class Pervaporation:
    membrane: Membrane
    mixture: Mixture
    conditions: Conditions
    ideal: bool = True

    def _get_partial_fluxes_from_permeate_composition(
        self,
        first_component_permeance: float,
        second_component_permeance: float,
        permeate_composition: Composition,
        feed_composition: Composition,
        feed_temperature: float,
        permeate_temperature: typing.Optional[float],
    ) -> typing.Tuple[float, float]:

        feed_nrtl_partial_pressures = get_nrtl_partial_pressures(
            feed_temperature, self.mixture, feed_composition
        )
        permeate_nrtl_partial_pressures = get_nrtl_partial_pressures(
            permeate_temperature, self.mixture, permeate_composition
        )

        return (
            first_component_permeance
            * (feed_nrtl_partial_pressures[0] - permeate_nrtl_partial_pressures[0]),
            second_component_permeance
            * (feed_nrtl_partial_pressures[1] - permeate_nrtl_partial_pressures[1]),
        )

    def calculate_partial_fluxes(
        self,
        feed_temperature: float,
        permeate_temperature: float,
        composition: Composition,
        precision: float = 5e-5,
    ) -> typing.Tuple[float, float]:
        first_component_permeance = self.membrane.get_permeance(
            feed_temperature, self.mixture.first_component
        )
        second_component_permeance = self.membrane.get_permeance(
            feed_temperature, self.mixture.second_component
        )

        initial_fluxes: typing.Tuple[float, float] = numpy.multiply(
            (first_component_permeance, second_component_permeance),
            get_nrtl_partial_pressures(feed_temperature, self.mixture, composition),
        )
        permeate_composition = _get_permeate_composition_from_fluxes(initial_fluxes)
        d = 1

        while d >= precision:
            permeate_composition_new = _get_permeate_composition_from_fluxes(
                self._get_partial_fluxes_from_permeate_composition(
                    first_component_permeance=first_component_permeance,
                    second_component_permeance=second_component_permeance,
                    permeate_composition=permeate_composition,
                    feed_composition=composition,
                    feed_temperature=feed_temperature,
                    permeate_temperature=permeate_temperature,
                )
            )
            d = max(
                numpy.absolute(
                    numpy.subtract(permeate_composition_new, permeate_composition)
                )
            )
            permeate_composition = permeate_composition_new
            # TODO: max iter and logs!!!
        return self._get_partial_fluxes_from_permeate_composition(
            first_component_permeance=first_component_permeance,
            second_component_permeance=second_component_permeance,
            permeate_composition=permeate_composition,
            feed_composition=composition,
            feed_temperature=feed_temperature,
            permeate_temperature=permeate_temperature,
        )

    # Calculate Permeate composition for at the given conditions
    def calculate_permeate_composition(
        self,
        feed_temperature: float,
        permeate_temperature: float,
        composition: Composition,
        precision: float,
    ) -> Composition:
        x = self.calculate_partial_fluxes(
            feed_temperature, permeate_temperature, composition, precision
        )
        return Composition(x[0] / numpy.sum(x), type=CompositionType("weight"))

    def calculate_separation_factor(
        self,
        feed_temperature: float,
        permeate_temperature: float,
        composition: Composition,
        precision: float,
    ) -> float:
        perm_comp = self.calculate_permeate_composition(
            feed_temperature, permeate_temperature, composition, precision
        )
        return (composition.second / (1 - composition.second)) / (
            perm_comp.second / (1 - perm_comp.second)
        )

    def get_ideal_diffusion_curve(
        self,
        feed_temperature: float,
        permeate_temperature: float,
        compositions: typing.List[Composition],
        precision,
    ) -> DiffusionCurve:
        return DiffusionCurve(
            mixture=self.mixture,
            membrane_name=self.membrane.name,
            feed_temperature=feed_temperature,
            permeate_temperature=permeate_temperature,
            compositions=compositions,
            partial_fluxes=[
                self.calculate_partial_fluxes(
                    feed_temperature,
                    permeate_temperature,
                    composition,
                    precision,
                )
                for composition in compositions
            ],
        )

    def model_ideal_isothermal_process(
        self,
        conditions: Conditions,
        number_of_steps: int,
        d_time_step_hours: float,
        precision: float = 5e-5,
    ) -> ProcessModel:
        isothermal = True
        time_steps = [d_time_step * step for step in number_of_steps]
        feed_temperature = conditions.feed_temperature
        permeate_temperature = conditions.permeate_temperature
        partial_fluxes = [tuple] * number_of_steps
        permeate_composition = [Composition] * number_of_steps
        feed_composition = [Composition] * number_of_steps
        feed_composition[0] = conditions.initial_feed_composition
        feed_evaporation_heat = [float] * number_of_steps
        permeate_condensation_heat = [float] * number_of_steps
        feed_mass = [float] * number_of_steps
        feed_mass[0] = conditions.feed_amount
        area = conditions.membrane_area

        first_component = self.mixture.first_component
        second_component = self.mixture.second_component
        evaporation_heat_1 = (
            first_component.get_vaporisation_heat(feed_temperature)
            / first_component.molecular_weight
            * 1000
        )
        evaporation_heat_2 = (
            second_component.get_vaporisation_heat(feed_temperature)
            / first_component.molecular_weight
            * 1000
        )
        condensation_heat_1 = (
            first_component.get_vaporisation_heat(permeate_temperature)
            / first_component.molecular_weight
            * 1000
        )
        condensation_heat_2 = (
            first_component.get_vaporisation_heat(permeate_temperature)
            / first_component.molecular_weight
            * 1000
        )

        heat_capacity_1 = first_component.get_specific_heat(feed_temperature)
        heat_capacity_2 = second_component.get_specific_heat(feed_temperature)

        for step in time_steps:
            partial_fluxes[step] = self.calculate_partial_fluxes(
                feed_temperature,
                permeate_temperature,
                feed_composition[step],
                precision,
            )

            permeate_composition[step] = Composition(
                partial_fluxes[step][0]
                / (sum(partial_fluxes[step]), CompositionType.weight)
            )

            d_mass_1 = partial_fluxes[step][0] * area * d_time_step_hours
            d_mass_2 = partial_fluxes[step][1] * area * d_time_step_hours

            feed_evaporation_heat[step] = (
                evaporation_heat_1 * d_mass_1 + evaporation_heat_2 * d_mass_2
            )
            permeate_condensation_heat[step] = (
                condensation_heat_1 * d_mass_1
                + condensation_heat_2 * d_mass_2
                + (heat_capacity_1 * d_mass_1 + heat_capacity_2 * d_mass_2)
                * (feed_temperature - permeate_temperature)
            )

            feed_mass[step + 1] = feed_mass[step] - d_mass_1 - d_mass_2

            feed_composition[step + 1] = Composition(
                (feed_composition[step].p * feed_mass[step] - d_mass_1)
                / feed_mass[step + 1],
                CompositionType.weight,
            )
        return ProcessModel()

    def model_ideal_non_isothermal_process(
        self,
        conditions: Conditions,
        number_of_steps: int,
        d_time_step_hours: float,
        precision: float = 5e-5,
    ) -> ProcessModel:
        return 0

    def model_non_ideal_process(
        self, conditions: Conditions, diffusion_curve: DiffusionCurve
    ):
        pass
