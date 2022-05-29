import typing

import attr
import numpy

from conditions import Conditions
from diffusion_curve import DiffusionCurve
from membrane import Membrane
from mixture import Composition, CompositionType, Mixture, get_nrtl_partial_pressures
from process import ProcessModel


def get_permeate_composition_from_fluxes(
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
    ideal: typing.Optional[bool] = True

    def get_partial_fluxes_from_permeate_composition(
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
        if permeate_temperature is None:
            permeate_nrtl_partial_pressures = (0, 0)

        else:
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
        composition: Composition,
        precision: float = 5e-5,
        permeate_temperature: typing.Optional[float] = None,
        first_component_permeance: typing.Optional[float] = None,
        second_component_permeance: typing.Optional[float] = None,
    ) -> typing.Tuple[float, float]:
        if second_component_permeance is None or first_component_permeance is None:
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
        permeate_composition = get_permeate_composition_from_fluxes(initial_fluxes)
        d = 1

        while d >= precision:
            try:
                permeate_composition_new = get_permeate_composition_from_fluxes(
                    self.get_partial_fluxes_from_permeate_composition(
                        first_component_permeance=first_component_permeance,
                        second_component_permeance=second_component_permeance,
                        permeate_composition=permeate_composition,
                        feed_composition=composition,
                        feed_temperature=feed_temperature,
                        permeate_temperature=permeate_temperature,
                    )
                )
                d = max(
                    abs(permeate_composition_new.first - permeate_composition.first),
                    abs(permeate_composition_new.second - permeate_composition.second),
                )
            except:
                raise ValueError(
                    "Partial fluxes are not defined in the stated conditions range"
                )
            else:
                permeate_composition = permeate_composition_new

            # TODO: max iter and logs!!!
        return self.get_partial_fluxes_from_permeate_composition(
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
        composition: Composition,
        precision: typing.Optional[float] = 5e-5,
        permeate_temperature: typing.Optional[float] = None,
    ) -> Composition:
        if permeate_temperature is None:
            x = self.calculate_partial_fluxes(feed_temperature, composition, precision)
        else:
            x = self.calculate_partial_fluxes(
                feed_temperature, composition, precision, permeate_temperature
            )
        return Composition(x[0] / numpy.sum(x), type=CompositionType("weight"))

    def calculate_separation_factor(
        self,
        feed_temperature: float,
        composition: Composition,
        permeate_temperature: typing.Optional[float] = None,
        precision: typing.Optional[float] = 5e-5,
    ) -> float:
        perm_comp = self.calculate_permeate_composition(
            feed_temperature, composition, precision, permeate_temperature
        )
        return (composition.second / (1 - composition.second)) / (
            perm_comp.second / (1 - perm_comp.second)
        )

    def get_ideal_diffusion_curve(
        self,
        feed_temperature: float,
        compositions: typing.List[Composition],
        permeate_temperature: typing.Optional[float] = None,
        precision: typing.Optional[float] = 5e-5,
    ) -> DiffusionCurve:
        return DiffusionCurve(
            mixture=self.mixture,
            membrane_name=self.membrane.name,
            feed_temperature=feed_temperature,
            permeate_temperature=permeate_temperature,
            feed_compositions=compositions,
            partial_fluxes=[
                self.calculate_partial_fluxes(
                    feed_temperature,
                    composition,
                    precision,
                    permeate_temperature,
                )
                for composition in compositions
            ],
        )

    def model_ideal_isothermal_process(
        self,
        number_of_steps: int,
        delta_hours: float,
        conditions: Conditions,
        precision: typing.Optional[float] = 5e-5,
    ) -> ProcessModel:

        time: typing.List[float] = [
            delta_hours * step for step in range(number_of_steps)
        ]

        partial_fluxes: typing.List[typing.Tuple[float, float]] = [(0.0, 0.0)] * (
            number_of_steps + 1
        )

        first_component_permeance = self.membrane.get_permeance(
            conditions.feed_temperature, self.mixture.first_component
        )
        second_component_permeance = self.membrane.get_permeance(
            conditions.feed_temperature, self.mixture.second_component
        )
        permeances: typing.List[typing.Tuple[float, float]] = [
            (first_component_permeance, second_component_permeance)
        ] * (number_of_steps + 1)

        permeate_composition: typing.List[Composition] = []
        feed_composition: typing.List[Composition] = [
            conditions.initial_feed_composition
        ]

        feed_evaporation_heat: typing.List[float] = []
        permeate_condensation_heat: typing.List[float] = []
        feed_mass: typing.List[float] = [conditions.feed_amount]

        evaporation_heat_1 = (
            self.mixture.first_component.get_vaporisation_heat(
                conditions.feed_temperature
            )
            / self.mixture.first_component.molecular_weight
            * 1000
        )
        evaporation_heat_2 = (
            self.mixture.second_component.get_vaporisation_heat(
                conditions.feed_temperature
            )
            / self.mixture.first_component.molecular_weight
            * 1000
        )
        condensation_heat_1 = (
            self.mixture.first_component.get_vaporisation_heat(
                conditions.permeate_temperature
            )
            / self.mixture.first_component.molecular_weight
            * 1000
        )
        condensation_heat_2 = (
            self.mixture.first_component.get_vaporisation_heat(
                conditions.permeate_temperature
            )
            / self.mixture.first_component.molecular_weight
            * 1000
        )

        cooling_heat_1 = self.mixture.first_component.get_cooling_heat(
            conditions.permeate_temperature, conditions.feed_temperature
        )
        cooling_heat_2 = self.mixture.second_component.get_cooling_heat(
            conditions.permeate_temperature, conditions.feed_temperature
        )

        for step in range(len(time)):
            partial_fluxes[step] = self.calculate_partial_fluxes(
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
            mixture=self.mixture,
            membrane_name=self.membrane.name,
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

    def model_ideal_non_isothermal_process(
        self,
        conditions: Conditions,
        number_of_steps: int,
        d_time_step_hours: float,
        precision: float = 5e-5,
    ) -> ProcessModel:
        time_steps = [d_time_step_hours * step for step in range(number_of_steps)]
        feed_temperature = [] * number_of_steps
        feed_temperature[0] = conditions.feed_temperature
        permeate_temperature = [float] * number_of_steps
        permeate_temperature[0] = conditions.permeate_temperature
        partial_fluxes = [] * number_of_steps
        permeances = [] * number_of_steps
        permeate_composition = [Composition] * number_of_steps
        feed_composition = [] * number_of_steps
        feed_composition[0] = conditions.initial_feed_composition.to_weight(
            self.mixture
        )
        feed_evaporation_heat = [] * number_of_steps
        permeate_condensation_heat = [] * number_of_steps
        feed_mass = [] * number_of_steps
        feed_mass[0] = conditions.feed_amount
        area = conditions.membrane_area

        first_component = self.mixture.first_component
        second_component = self.mixture.second_component

        for step in range(len(time_steps)):
            permeate_temperature[step] = permeate_temperature[0]

            evaporation_heat_1 = (
                first_component.get_vaporisation_heat(feed_temperature[step])
                / first_component.molecular_weight
                * 1000
            )
            evaporation_heat_2 = (
                second_component.get_vaporisation_heat(feed_temperature[step])
                / first_component.molecular_weight
                * 1000
            )
            condensation_heat_1 = (
                first_component.get_vaporisation_heat(permeate_temperature[step])
                / first_component.molecular_weight
                * 1000
            )
            condensation_heat_2 = (
                first_component.get_vaporisation_heat(permeate_temperature[step])
                / first_component.molecular_weight
                * 1000
            )

            specific_heat_1 = first_component.get_cooling_heat(
                feed_temperature[step], permeate_temperature[step]
            )
            specific_heat_2 = second_component.get_cooling_heat(
                feed_temperature[step], permeate_temperature[step]
            )
            heat_capacity_1 = (
                first_component.get_specific_heat(feed_temperature[step])
                / first_component.molecular_weight
            )
            heat_capacity_2 = (
                second_component.get_specific_heat(feed_temperature[step])
                / second_component.molecular_weight
            )
            feed_heat_capacity = (
                feed_composition[step].first * heat_capacity_1
                + feed_composition[step].second * heat_capacity_2
            )

            first_component_permeance = self.membrane.get_permeance(
                feed_temperature[step], self.mixture.first_component
            )
            second_component_permeance = self.membrane.get_permeance(
                feed_temperature[step], self.mixture.second_component
            )
            permeances[step] = (
                first_component_permeance,
                second_component_permeance,
            )

            partial_fluxes[step] = self.calculate_partial_fluxes(
                feed_temperature[step],
                feed_composition[step],
                precision,
                permeate_temperature[step],
            )

            permeate_composition[step] = Composition(
                p=partial_fluxes[step][0] / (sum(partial_fluxes[step])),
                type=CompositionType("weight"),
            )

            d_mass_1 = partial_fluxes[step][0] * area * d_time_step_hours
            d_mass_2 = partial_fluxes[step][1] * area * d_time_step_hours

            feed_evaporation_heat[step] = (
                evaporation_heat_1 * d_mass_1 + evaporation_heat_2 * d_mass_2
            )
            permeate_condensation_heat[step] = (
                condensation_heat_1 * d_mass_1
                + condensation_heat_2 * d_mass_2
                + (specific_heat_1 * d_mass_1 + specific_heat_2 * d_mass_2)
                * (feed_temperature[step] - permeate_temperature[step])
            )

            feed_mass[step + 1] = feed_mass[step] - d_mass_1 - d_mass_2

            feed_composition[step + 1] = Composition(
                p=(feed_composition[step].p * feed_mass[step] - d_mass_1)
                / feed_mass[step + 1],
                type=CompositionType("weight"),
            )

            feed_temperature[step + 1] = feed_temperature[step] - (
                feed_evaporation_heat[step] / (feed_heat_capacity * feed_mass[step])
            )

        return ProcessModel(
            mixture=self.mixture,
            membrane_name=self.membrane.name,
            isothermal=False,
            feed_temperature=feed_temperature,
            permeate_temperature=permeate_temperature,
            feed_composition=feed_composition,
            permeate_composition=permeate_composition,
            feed_mass=feed_mass,
            partial_fluxes=partial_fluxes,
            permeances=permeances,
            time=time_steps,
            feed_evaporation_heat=feed_evaporation_heat,
            permeate_condensation_heat=permeate_condensation_heat,
            initial_condtioins=conditions,
            IsTimeDefined=True,
            comments="",
        )

    def model_non_ideal_process(
        self,
        conditions: Conditions,
        diffusion_curve: DiffusionCurve,
    ):
        pass
