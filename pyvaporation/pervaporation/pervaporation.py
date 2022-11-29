import datetime
import typing
from datetime import datetime

import attr
import numpy

from ..conditions import Conditions
from ..diffusion_curve import DiffusionCurve, DiffusionCurveSet
from ..membrane import Membrane
from ..mixtures import Composition, CompositionType, Mixture, get_nrtl_partial_pressures
from ..optimizer import Measurements, find_best_fit
from ..permeance import Permeance, Units
from ..process import ProcessModel
from ..utils import R


def get_permeate_composition_from_fluxes(
    fluxes: typing.Tuple[float, float],
) -> Composition:
    return Composition(
        p=fluxes[0] / sum(fluxes),
        type=CompositionType.weight,
    )


@attr.s(auto_attribs=True)
class Pervaporation:
    membrane: Membrane
    mixture: Mixture

    def get_partial_fluxes_from_permeate_composition(
        self,
        first_component_permeance: Permeance,
        second_component_permeance: Permeance,
        permeate_composition: Composition,
        feed_composition: Composition,
        feed_temperature: float,
        permeate_temperature: typing.Optional[float] = None,
        permeate_pressure: typing.Optional[float] = None,
    ) -> typing.Tuple[float, float]:
        """
        Calculates partial fluxes at a given Permeate composition, accounting for the driving force change
        Either permeate temperature or permeate pressure could be stated
        :param first_component_permeance - Permeance of the first test_components
        :param second_component_permeance - Permeance of the second test_components
        :param permeate_composition - permeate composition
        :param feed_composition - feed composition
        :param feed_temperature - temperature, K
        :param permeate_temperature - permeate temperature, K , if not specified permeate pressure is considered 0 kPa
        :param permeate_pressure - permeate pressure, kPa , if not specified permeate pressure is considered 0 kPa
        """

        feed_nrtl_partial_pressures = get_nrtl_partial_pressures(
            feed_temperature, self.mixture, feed_composition
        )
        if permeate_temperature is None and permeate_pressure is None:
            permeate_nrtl_partial_pressures = (0, 0)

        elif permeate_temperature is not None and permeate_pressure is None:
            permeate_nrtl_partial_pressures = get_nrtl_partial_pressures(
                permeate_temperature, self.mixture, permeate_composition
            )
        elif permeate_pressure is not None and permeate_temperature is None:
            permeate_nrtl_partial_pressures = (
                permeate_pressure * permeate_composition.first,
                permeate_pressure * permeate_composition.second,
            )

        else:
            raise ValueError(
                "Either permeate temperature or permeate pressure could be stated not both"
            )

        return (
            first_component_permeance.value
            * (feed_nrtl_partial_pressures[0] - permeate_nrtl_partial_pressures[0]),
            second_component_permeance.value
            * (feed_nrtl_partial_pressures[1] - permeate_nrtl_partial_pressures[1]),
        )

    def calculate_partial_fluxes(
        self,
        feed_temperature: float,
        composition: Composition,
        precision: float = 3e-4,
        permeate_temperature: typing.Optional[float] = None,
        permeate_pressure: typing.Optional[float] = None,
        first_component_permeance: typing.Optional[Permeance] = None,
        second_component_permeance: typing.Optional[Permeance] = None,
    ) -> typing.Tuple[float, float]:
        """
        Calculates partial fluxes of the test_components at specified conditions.
        Either permeate temperature or permeate pressure could be stated
        :param feed_temperature: Feed temperature, K
        :param composition: Feed composition
        :param precision: Precision in obtained permeate composition, by default is 5e-5
        :param permeate_temperature: Permeate temperature, if not specified permeate pressure is set to 0 kPa
        :param permeate_pressure - permeate pressure, kPa , if not specified permeate pressure is considered 0 kPa
        :param first_component_permeance: Permeance of the first test_components, if not specified is calculated
        :param second_component_permeance: Permeance of the second test_components, if not specified is calculated
        :return: Partial fluxes of test_components as a tuple
        """
        if second_component_permeance is None or first_component_permeance is None:
            first_component_permeance = self.membrane.get_permeance(
                feed_temperature, self.mixture.first_component
            ).convert(
                to_units=Units().kg_m2_h_kPa, component=self.mixture.first_component
            )
            second_component_permeance = self.membrane.get_permeance(
                feed_temperature, self.mixture.second_component
            ).convert(
                to_units=Units().kg_m2_h_kPa, component=self.mixture.second_component
            )

        initial_fluxes: typing.Tuple[float, float] = numpy.multiply(
            (first_component_permeance.value, second_component_permeance.value),
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
                        permeate_pressure=permeate_pressure,
                    )
                )
                d = max(
                    abs(permeate_composition_new.first - permeate_composition.first),
                    abs(permeate_composition_new.second - permeate_composition.second),
                )
            except ():
                raise ValueError(
                    "Partial fluxes are not defined in the stated conditions range"
                )
            else:
                permeate_composition = permeate_composition_new

        return self.get_partial_fluxes_from_permeate_composition(
            first_component_permeance=first_component_permeance,
            second_component_permeance=second_component_permeance,
            permeate_composition=permeate_composition,
            feed_composition=composition,
            feed_temperature=feed_temperature,
            permeate_temperature=permeate_temperature,
            permeate_pressure=permeate_pressure,
        )

    def calculate_permeate_composition(
        self,
        feed_temperature: float,
        composition: Composition,
        precision: typing.Optional[float] = 5e-5,
        permeate_temperature: typing.Optional[float] = None,
        permeate_pressure: typing.Optional[float] = None,
    ) -> Composition:
        """
        Calculates permeate composition at given conditions
        Either permeate temperature or permeate pressure could be stated
        :param feed_temperature: Feed temperature, K
        :param composition: Feed Composition
        :param precision: Precision in obtained permeate composition, by default is 5e-5
        :param permeate_temperature: Permeate temperature, if not specified permeate pressure is set to 0 kPa
        :param permeate_pressure - permeate pressure, kPa , if not specified permeate pressure is considered 0 kPa
        :return: Permeate composition in weight %
        """

        x = self.calculate_partial_fluxes(
            feed_temperature,
            composition,
            precision,
            permeate_temperature,
            permeate_pressure,
        )
        return Composition(x[0] / numpy.sum(x), type=CompositionType.weight)

    def calculate_separation_factor(
        self,
        feed_temperature: float,
        composition: Composition,
        permeate_temperature: typing.Optional[float] = None,
        permeate_pressure: typing.Optional[float] = None,
        precision: typing.Optional[float] = 5e-5,
    ) -> float:
        """
        Calculates separation factor at given conditions
        """
        perm_comp = self.calculate_permeate_composition(
            feed_temperature,
            composition,
            precision,
            permeate_temperature,
            permeate_pressure,
        )
        return (composition.second / composition.first) / (
            perm_comp.second / perm_comp.first
        )

    def ideal_diffusion_curve(
        self,
        feed_temperature: float,
        compositions: typing.List[Composition],
        permeate_temperature: typing.Optional[float] = None,
        permeate_pressure: typing.Optional[float] = None,
        precision: typing.Optional[float] = 5e-5,
    ) -> DiffusionCurve:
        """
        Models Ideal Diffusion curve of a specified membrane, at a given temperature, for a given Mixture
        if Ideal experiments for both test_components are available.
        Either permeate temperature or permeate pressure could be stated;
        :param feed_temperature: Feed temperature, K
        :param compositions: List of compositions to model parameters at
        :param permeate_temperature: Permeate temperature, if not specified permeate pressure is set to 0 kPa
        :param permeate_pressure - Permeate pressure, kPa , if not specified permeate pressure is considered 0 kPa
        :param precision: Precision in obtained permeate composition, by default is 5e-5
        :return: A DiffusionCurve Object
        """
        return DiffusionCurve(
            mixture=self.mixture,
            membrane_name=self.membrane.name,
            feed_temperature=feed_temperature,
            permeate_temperature=permeate_temperature,
            permeate_pressure=permeate_pressure,
            feed_compositions=compositions,
            partial_fluxes=[
                self.calculate_partial_fluxes(
                    feed_temperature,
                    composition,
                    precision,
                    permeate_temperature,
                    permeate_pressure,
                )
                for composition in compositions
            ],
            comments=(
                str(self.membrane.name)
                + " "
                + str(self.mixture.first_component.name)
                + " / "
                + str(self.mixture.second_component.name)
                + " "
                + datetime.now().strftime("%m/%d/%Y, %H:%M")
            ),
        )

    def ideal_isothermal_process(
        self,
        number_of_steps: int,
        delta_hours: float,
        conditions: Conditions,
        precision: typing.Optional[float] = 5e-5,
    ) -> ProcessModel:
        """
        Models mass and heat balance of an Ideal (constant Permeance) Isothermal Pervaporation Process
        :param number_of_steps: Number of time steps to include in the model
        :param delta_hours: The duration of each step in hours
        :param conditions: Conditions object, where initial conditions are specified
        :param precision: Precision in obtained permeate composition, by default is 5e-5
        :return: A ProcessModel Object
        """

        time: typing.List[float] = [
            delta_hours * step for step in range(number_of_steps)
        ]

        partial_fluxes: typing.List[typing.Tuple[float, float]] = []

        first_component_permeance = self.membrane.get_permeance(
            conditions.initial_feed_temperature, self.mixture.first_component
        )
        second_component_permeance = self.membrane.get_permeance(
            conditions.initial_feed_temperature, self.mixture.second_component
        )
        permeances: typing.List[typing.Tuple[Permeance, Permeance]] = [
            (
                first_component_permeance,
                second_component_permeance,
            )
        ] * number_of_steps

        permeate_composition: typing.List[Composition] = []
        feed_composition: typing.List[Composition] = [
            conditions.initial_feed_composition.to_weight(self.mixture)
        ]

        feed_evaporation_heat: typing.List[float] = []
        permeate_condensation_heat: typing.List[typing.Optional[float]] = []
        feed_mass: typing.List[float] = [conditions.initial_feed_amount]

        evaporation_heat_1 = (
            self.mixture.first_component.get_vaporisation_heat(
                conditions.initial_feed_temperature
            )
            / self.mixture.first_component.molecular_weight
            * 1000
        )
        evaporation_heat_2 = (
            self.mixture.second_component.get_vaporisation_heat(
                conditions.initial_feed_temperature
            )
            / self.mixture.first_component.molecular_weight
            * 1000
        )
        if conditions.permeate_temperature is None:
            condensation_heat_1 = None
            condensation_heat_2 = None
            cooling_heat_1 = None
            cooling_heat_2 = None

        else:
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
                conditions.permeate_temperature, conditions.initial_feed_temperature
            )
            cooling_heat_2 = self.mixture.second_component.get_cooling_heat(
                conditions.permeate_temperature, conditions.initial_feed_temperature
            )

        for step in range(len(time)):
            partial_fluxes.append(
                self.calculate_partial_fluxes(
                    feed_temperature=conditions.initial_feed_temperature,
                    composition=feed_composition[step],
                    precision=precision,
                    permeate_temperature=conditions.permeate_temperature,
                    permeate_pressure=conditions.permeate_pressure,
                    first_component_permeance=first_component_permeance,
                    second_component_permeance=second_component_permeance,
                )
            )

            permeate_composition.append(
                Composition(
                    p=partial_fluxes[step][0] / (sum(partial_fluxes[step])),
                    type=CompositionType.weight,
                )
            )

            d_mass_1 = partial_fluxes[step][0] * conditions.membrane_area * delta_hours
            d_mass_2 = partial_fluxes[step][1] * conditions.membrane_area * delta_hours

            feed_evaporation_heat.append(
                evaporation_heat_1 * d_mass_1 + evaporation_heat_2 * d_mass_2
            )
            if conditions.permeate_temperature is None:
                permeate_condensation_heat.append(None)
            else:
                permeate_condensation_heat.append(
                    condensation_heat_1 * d_mass_1
                    + condensation_heat_2 * d_mass_2
                    + (cooling_heat_1 * d_mass_1 + cooling_heat_2 * d_mass_2)
                    * (
                        conditions.initial_feed_temperature
                        - conditions.permeate_temperature
                    )
                )

            feed_mass.append(feed_mass[step] - d_mass_1 - d_mass_2)

            feed_composition.append(
                Composition(
                    p=(feed_composition[step].p * feed_mass[step] - d_mass_1)
                    / feed_mass[step + 1],
                    type=CompositionType.weight,
                )
            )

        feed_composition.pop(-1)
        feed_mass.pop(-1)

        return ProcessModel(
            mixture=self.mixture,
            membrane_name=self.membrane.name,
            feed_temperature=[conditions.initial_feed_temperature] * number_of_steps,
            feed_compositions=feed_composition,
            permeate_composition=permeate_composition,
            permeate_temperature=[conditions.permeate_temperature] * number_of_steps,
            permeate_pressure=[conditions.permeate_pressure] * number_of_steps,
            feed_mass=feed_mass,
            partial_fluxes=partial_fluxes,
            permeances=permeances,
            time=time,
            feed_evaporation_heat=feed_evaporation_heat,
            permeate_condensation_heat=permeate_condensation_heat,
            initial_conditions=conditions,
            comments=(
                f"{self.membrane.name} {self.mixture.first_component.name} / {self.mixture.second_component.name}"
                f"Ideal Process Model" + datetime.now().strftime("%m/%d/%Y, %H:%M")
            ),
            membrane_path=self.membrane.path,
        )

    def ideal_non_isothermal_process(
        self,
        conditions: Conditions,
        number_of_steps: int,
        delta_hours: float,
        precision: typing.Optional[float] = 5e-5,
    ) -> ProcessModel:
        """
        Models mass and heat balance of an Ideal (constant Permeance) Non-Isothermal Pervaporation Process.
        The temperature program maybe specified in Conditions, by including a TemperatureProgram object.
        If temperature program is not specified, models self-cooling process;
        :param number_of_steps: Number of time steps to include in the model
        :param delta_hours: The duration of each step in hours
        :param conditions: Conditions object, where initial conditions are specified
        :param precision: Precision in obtained permeate composition, by default is 5e-5
        :return: A ProcessModel Object
        """

        time: typing.List[float] = [
            delta_hours * step for step in range(number_of_steps)
        ]

        feed_temperature: typing.List[float] = [conditions.initial_feed_temperature]

        partial_fluxes: typing.List[typing.Tuple[float, float]] = []

        permeances: typing.List[typing.Tuple[Permeance, Permeance]] = []

        permeate_composition: typing.List[Composition] = []

        feed_composition: typing.List[Composition] = [
            conditions.initial_feed_composition.to_weight(self.mixture)
        ]

        feed_evaporation_heat: typing.List[float] = []

        permeate_condensation_heat: typing.List[typing.Optional[float]] = []

        feed_mass: typing.List[float] = [conditions.initial_feed_amount]

        for step in range(len(time)):

            evaporation_heat_1 = (
                self.mixture.first_component.get_vaporisation_heat(
                    feed_temperature[step]
                )
                / self.mixture.first_component.molecular_weight
                * 1000
            )
            evaporation_heat_2 = (
                self.mixture.second_component.get_vaporisation_heat(
                    feed_temperature[step]
                )
                / self.mixture.second_component.molecular_weight
                * 1000
            )

            heat_capacity_1 = (
                self.mixture.first_component.get_specific_heat(feed_temperature[step])
                / self.mixture.first_component.molecular_weight
            )
            heat_capacity_2 = (
                self.mixture.second_component.get_specific_heat(feed_temperature[step])
                / self.mixture.second_component.molecular_weight
            )
            feed_heat_capacity = (
                feed_composition[step].first * heat_capacity_1
                + feed_composition[step].second * heat_capacity_2
            )

            permeances.append(
                (
                    self.membrane.get_permeance(
                        feed_temperature[step], self.mixture.first_component
                    ),
                    self.membrane.get_permeance(
                        feed_temperature[step], self.mixture.second_component
                    ),
                )
            )

            partial_fluxes.append(
                self.calculate_partial_fluxes(
                    feed_temperature=feed_temperature[step],
                    composition=feed_composition[step],
                    precision=precision,
                    permeate_temperature=conditions.permeate_temperature,
                    permeate_pressure=conditions.permeate_pressure,
                    first_component_permeance=permeances[step][0],
                    second_component_permeance=permeances[step][1],
                )
            )

            permeate_composition.append(
                Composition(
                    p=partial_fluxes[step][0] / (sum(partial_fluxes[step])),
                    type=CompositionType.weight,
                )
            )

            d_mass_1 = partial_fluxes[step][0] * conditions.membrane_area * delta_hours
            d_mass_2 = partial_fluxes[step][1] * conditions.membrane_area * delta_hours

            if conditions.permeate_temperature is None:
                permeate_condensation_heat.append(None)
            else:
                condensation_heat_1 = (
                    self.mixture.first_component.get_vaporisation_heat(
                        conditions.permeate_temperature
                    )
                    / self.mixture.first_component.molecular_weight
                    * 1000
                )
                condensation_heat_2 = (
                    self.mixture.second_component.get_vaporisation_heat(
                        conditions.permeate_temperature
                    )
                    / self.mixture.second_component.molecular_weight
                    * 1000
                )

                specific_heat_1 = self.mixture.first_component.get_cooling_heat(
                    feed_temperature[step], conditions.permeate_temperature
                )
                specific_heat_2 = self.mixture.second_component.get_cooling_heat(
                    feed_temperature[step], conditions.permeate_temperature
                )

                permeate_condensation_heat.append(
                    condensation_heat_1 * d_mass_1
                    + condensation_heat_2 * d_mass_2
                    + (specific_heat_1 * d_mass_1 + specific_heat_2 * d_mass_2)
                    * (feed_temperature[step] - conditions.permeate_temperature)
                )

            feed_evaporation_heat.append(
                evaporation_heat_1 * d_mass_1 + evaporation_heat_2 * d_mass_2
            )

            feed_mass.append(feed_mass[step] - d_mass_1 - d_mass_2)

            feed_composition.append(
                Composition(
                    p=(feed_composition[step].p * feed_mass[step] - d_mass_1)
                    / feed_mass[step + 1],
                    type=CompositionType.weight,
                )
            )

            if conditions.temperature_program is None:
                feed_temperature.append(
                    feed_temperature[step]
                    - (
                        feed_evaporation_heat[step]
                        / (feed_heat_capacity * feed_mass[step])
                    )
                )
            else:
                feed_temperature.append(
                    conditions.temperature_program.program(time[step] + delta_hours)
                )

        feed_mass.pop(-1)
        feed_composition.pop(-1)
        feed_temperature.pop(-1)

        return ProcessModel(
            mixture=self.mixture,
            membrane_name=self.membrane.name,
            feed_temperature=feed_temperature,
            permeate_temperature=[conditions.permeate_temperature] * number_of_steps,
            permeate_pressure=[conditions.permeate_pressure] * number_of_steps,
            feed_compositions=feed_composition,
            permeate_composition=permeate_composition,
            feed_mass=feed_mass,
            partial_fluxes=partial_fluxes,
            permeances=permeances,
            time=time,
            feed_evaporation_heat=feed_evaporation_heat,
            permeate_condensation_heat=permeate_condensation_heat,
            initial_conditions=conditions,
            comments=(
                f"{self.membrane.name} {self.mixture.first_component.name} / {self.mixture.second_component.name}"
                f"Ideal Process Model" + datetime.now().strftime("%m/%d/%Y, %H:%M")
            ),
            membrane_path=self.membrane.path,
        )

    def non_ideal_diffusion_curve(
        self,
        diffusion_curve_set: DiffusionCurveSet,
        feed_temperature: float,
        initial_feed_composition: Composition,
        delta_composition: float,
        number_of_steps: int,
        permeate_temperature: typing.Optional[float] = None,
        permeate_pressure: typing.Optional[float] = None,
        initial_permeances: typing.Optional[typing.Tuple[Permeance, Permeance]] = None,
        precision: typing.Optional[float] = 5e-5,
        n_first: typing.Optional[int] = None,
        n_second: typing.Optional[int] = None,
        m_first: typing.Optional[int] = None,
        m_second: typing.Optional[int] = None,
        include_zero: bool = False,
    ):
        """
        The Fucntion models Non-Ideal Diffusion curve
        Based on a set of Diffusion curves measured at different temperatures;
        The modelling could be also performed based on a single diffusion curve:
        In that case the apparent activation energy of transport is considered constant
        and is calculated for each component based on the IdealExperiments default_membranes provided for the Membrane.
        :param diffusion_curve_set: DiffusionCurveSet used as a basis for fitting of the PervaporationFunctions
        :param feed_temperature: temperature of the feed mixture
        :param initial_feed_composition: initial composition of the feed
        :param delta_composition: a composition step size
        :param number_of_steps: number of modelling steps
        :param permeate_temperature: temperature of permeate
        :param permeate_pressure: pressure of permeate
        :param initial_permeances: initial permeances of components, for calculation of an FR parameter using
        data from DiffusionCurveSet object
        :param precision: precision of permeate composition calculation
        :param n_first: n parameter of the PervaporationFunction of the first component
        :param n_second: n parameter of the PervaporationFunction of the second component
        :param m_first: m parameter of the PervaporationFunction of the first component
        :param m_second: m parameter of the PervaporationFunction of the second component
        :param include_zero: bool parameter to force default points while fitting the PervaporationFunction
        :return: non-ideal diffusion curve
        """

        measurements_first = Measurements.from_diffusion_curves_first(
            diffusion_curve_set
        )
        measurements_second = Measurements.from_diffusion_curves_second(
            diffusion_curve_set
        )

        if len(diffusion_curve_set.diffusion_curves) == 1:
            pervaporation_function_temperature = diffusion_curve_set.diffusion_curves[
                0
            ].feed_temperature

            pervaporation_function_first = find_best_fit(
                data=measurements_first,
                n=n_first,
                m=0,
                include_zero=include_zero,
                component_index=0,
            )

            pervaporation_function_second = find_best_fit(
                data=measurements_second,
                n=n_second,
                m=0,
                include_zero=include_zero,
                component_index=1,
            )

            if pervaporation_function_temperature == feed_temperature:
                pass
            else:
                activation_energy_first = self.membrane.calculate_activation_energy(
                    self.mixture.first_component
                )
                activation_energy_second = self.membrane.calculate_activation_energy(
                    self.mixture.second_component
                )

                pervaporation_function_first = pervaporation_function_first * numpy.exp(
                    -pervaporation_function_first.b[0]
                    / pervaporation_function_temperature
                    + activation_energy_first / (R * pervaporation_function_temperature)
                )

                pervaporation_function_second = (
                    pervaporation_function_second
                    * numpy.exp(
                        -pervaporation_function_second.b[0]
                        / pervaporation_function_temperature
                        + activation_energy_second
                        / (R * pervaporation_function_temperature)
                    )
                )

                pervaporation_function_first.b[0] = activation_energy_first / R
                pervaporation_function_second.b[0] = activation_energy_second / R

        else:
            pervaporation_function_first = find_best_fit(
                data=measurements_first,
                n=n_first,
                m=m_first,
                include_zero=include_zero,
                component_index=0,
            )
            pervaporation_function_second = find_best_fit(
                data=measurements_second,
                n=n_second,
                m=m_second,
                include_zero=include_zero,
                component_index=1,
            )

        if initial_permeances is None:
            first_component_permeance = Permeance(
                value=pervaporation_function_first(
                    initial_feed_composition.first, feed_temperature
                )
            )

            second_component_permeance = Permeance(
                value=pervaporation_function_second(
                    initial_feed_composition.first, feed_temperature
                )
            )

        else:
            first_component_permeance = initial_permeances[0].convert(
                to_units=Units.kg_m2_h_kPa, component=self.mixture.first_component
            )
            second_component_permeance = initial_permeances[1].convert(
                to_units=Units.kg_m2_h_kPa, component=self.mixture.second_component
            )

        permeances: typing.List[typing.Tuple[Permeance, Permeance]] = [
            (first_component_permeance, second_component_permeance)
        ]
        facilitation_rate_first = (
            first_component_permeance.value
            / pervaporation_function_first(
                x=initial_feed_composition.first, t=feed_temperature
            )
        )
        facilitation_rate_second = (
            second_component_permeance.value
            / pervaporation_function_second(
                x=initial_feed_composition.first, t=feed_temperature
            )
        )

        partial_fluxes: typing.List[typing.Tuple[float, float]] = []
        compositions: typing.List[Composition] = [
            initial_feed_composition.to_weight(self.mixture)
        ]

        for i in range(number_of_steps + 1):

            compositions.append(
                Composition(
                    p=compositions[i].first + delta_composition,
                    type=CompositionType.weight,
                )
            )

            partial_fluxes.append(
                self.calculate_partial_fluxes(
                    feed_temperature=feed_temperature,
                    composition=compositions[i],
                    precision=precision,
                    permeate_temperature=permeate_temperature,
                    permeate_pressure=permeate_pressure,
                    first_component_permeance=permeances[i][0],
                    second_component_permeance=permeances[i][1],
                )
            )

            permeances.append(
                (
                    Permeance(
                        value=(
                            pervaporation_function_first(
                                compositions[i + 1].first, feed_temperature
                            )
                            * facilitation_rate_first
                        )
                    ),
                    Permeance(
                        value=(
                            pervaporation_function_second(
                                compositions[i + 1].first, feed_temperature
                            )
                            * facilitation_rate_second
                        )
                    ),
                )
            )

        compositions.pop(-1)
        permeances.pop(-1)

        return DiffusionCurve(
            mixture=self.mixture,
            membrane_name=self.membrane.name,
            feed_temperature=feed_temperature,
            permeate_temperature=permeate_temperature,
            permeate_pressure=permeate_pressure,
            feed_compositions=compositions,
            partial_fluxes=partial_fluxes,
            permeances=permeances,
            comments=(
                f"{self.membrane.name} {self.mixture.first_component.name} "
                f"/ {self.mixture.second_component.name}"
                + datetime.now().strftime("%m/%d/%Y, %H:%M")
            ),
        )

    def non_ideal_isothermal_process(
        self,
        conditions: Conditions,
        diffusion_curve_set: DiffusionCurveSet,
        number_of_steps: int,
        delta_hours: float,
        precision: typing.Optional[float] = 5e-5,
        initial_permeances: typing.Optional[typing.Tuple[Permeance, Permeance]] = None,
        n_first: typing.Optional[int] = None,
        m_first: typing.Optional[int] = None,
        n_second: typing.Optional[int] = None,
        m_second: typing.Optional[int] = None,
        include_zero: bool = False,
    ):
        """
        The function models Non-Ideal Isothermal Process
        Based on a set of Diffusion curves measured at different temperatures;
        The modelling could be also performed based on a single diffusion curve:
        In that case the apparent activation energy of transport is considered constant
        and is calculated for each component based on the IdealExperiments default_membranes provided for the Membrane.
        :param conditions: Initial Conditions of the Process
        :param diffusion_curve_set: A set of Diffusion curves picked for the Modelling form the Membrane
        :param number_of_steps: Number of time steps for modelling
        :param delta_hours: Size of each step in hours
        :param precision: Precision in obtained permeate composition, by default is 5e-5
        :param initial_permeances: Initial Permeances, should be stated if the Membranes swelling history is significant
        :param n_first: Optional parameter,
        indicates the order of the polynomial of the composition part of the Permeance function for the
        first component
        :param m_first: Optional parameter,
        indicates the order of the polynomial of the temperature part of the Permeance function for
        the first component
        :param n_second: Optional parameter,
        indicates the order of the polynomial of the composition part of the Permeance function for
        the second component
        :param m_second: Optional parameter,
        indicates the order of the polynomial of the temperature part of the Permeance function for
        the second component
        :param include_zero: if True, points:
         first_component_fraction = 0 first_component_permeance=0 for the first test_components
         first_component_fraction = 1 second_component_permeance=0 for the second test_components
         for each temperature are added to the measurements in order to improve obtained fits
        :return: ProcessModel object
        """
        for curve in diffusion_curve_set.diffusion_curves:
            for c in curve.feed_compositions:
                if c.type == CompositionType.molar:
                    c.to_weight(curve.mixture)

        time: typing.List[float] = [
            delta_hours * step for step in range(number_of_steps)
        ]

        partial_fluxes: typing.List[typing.Tuple[float, float]] = []

        permeate_composition: typing.List[Composition] = []
        feed_composition: typing.List[Composition] = [
            conditions.initial_feed_composition.to_weight(self.mixture)
        ]

        feed_evaporation_heat: typing.List[float] = []
        permeate_condensation_heat: typing.List[typing.Optional[float]] = []
        feed_mass: typing.List[float] = [conditions.initial_feed_amount]

        measurements_first = Measurements.from_diffusion_curves_first(
            diffusion_curve_set
        )
        measurements_second = Measurements.from_diffusion_curves_second(
            diffusion_curve_set
        )

        if len(diffusion_curve_set.diffusion_curves) == 1:
            pervaporation_function_temperature = diffusion_curve_set.diffusion_curves[
                0
            ].feed_temperature

            pervaporation_function_first = find_best_fit(
                data=measurements_first,
                n=n_first,
                m=0,
                include_zero=False,
                component_index=0,
            )
            pervaporation_function_second = find_best_fit(
                data=measurements_second,
                n=n_second,
                m=0,
                include_zero=False,
                component_index=1,
            )

            if (
                pervaporation_function_temperature
                == conditions.initial_feed_temperature
            ):
                pass
            else:
                activation_energy_first = self.membrane.calculate_activation_energy(
                    self.mixture.first_component
                )
                activation_energy_second = self.membrane.calculate_activation_energy(
                    self.mixture.second_component
                )

                pervaporation_function_first = pervaporation_function_first * numpy.exp(
                    -pervaporation_function_first.b[0]
                    / pervaporation_function_temperature
                    + activation_energy_first / (R * pervaporation_function_temperature)
                )

                pervaporation_function_second = (
                    pervaporation_function_second
                    * numpy.exp(
                        -pervaporation_function_second.b[0]
                        / pervaporation_function_temperature
                        + activation_energy_second
                        / (R * pervaporation_function_temperature)
                    )
                )

                pervaporation_function_first.b[0] = activation_energy_first / R
                pervaporation_function_second.b[0] = activation_energy_second / R

        else:
            pervaporation_function_first = find_best_fit(
                data=measurements_first,
                n=n_first,
                m=m_first,
                include_zero=include_zero,
                component_index=0,
            )
            pervaporation_function_second = find_best_fit(
                data=measurements_second,
                n=n_second,
                m=m_second,
                include_zero=include_zero,
                component_index=1,
            )

        if initial_permeances is None:
            first_component_permeance = Permeance(
                value=pervaporation_function_first(
                    feed_composition[0].first, conditions.initial_feed_temperature
                )
            )

            second_component_permeance = Permeance(
                value=pervaporation_function_second(
                    feed_composition[0].first, conditions.initial_feed_temperature
                )
            )

        else:
            first_component_permeance = initial_permeances[0].convert(
                to_units=Units.kg_m2_h_kPa, component=self.mixture.first_component
            )
            second_component_permeance = initial_permeances[1].convert(
                to_units=Units.kg_m2_h_kPa, component=self.mixture.second_component
            )

        permeances: typing.List[typing.Tuple[Permeance, Permeance]] = [
            (first_component_permeance, second_component_permeance)
        ]

        facilitation_rate_first = (
            first_component_permeance.value
            / pervaporation_function_first(
                x=conditions.initial_feed_composition.first,
                t=conditions.initial_feed_temperature,
            )
        )
        facilitation_rate_second = (
            second_component_permeance.value
            / pervaporation_function_second(
                x=conditions.initial_feed_composition.first,
                t=conditions.initial_feed_temperature,
            )
        )

        evaporation_heat_1 = (
            self.mixture.first_component.get_vaporisation_heat(
                conditions.initial_feed_temperature
            )
            / self.mixture.first_component.molecular_weight
            * 1000
        )
        evaporation_heat_2 = (
            self.mixture.second_component.get_vaporisation_heat(
                conditions.initial_feed_temperature
            )
            / self.mixture.first_component.molecular_weight
            * 1000
        )
        if conditions.permeate_temperature is None:
            condensation_heat_1 = None
            condensation_heat_2 = None
            cooling_heat_1 = None
            cooling_heat_2 = None

        else:
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
                conditions.permeate_temperature, conditions.initial_feed_temperature
            )
            cooling_heat_2 = self.mixture.second_component.get_cooling_heat(
                conditions.permeate_temperature, conditions.initial_feed_temperature
            )

        for step in range(len(time)):

            partial_fluxes.append(
                self.calculate_partial_fluxes(
                    feed_temperature=conditions.initial_feed_temperature,
                    composition=feed_composition[step],
                    precision=precision,
                    permeate_temperature=conditions.permeate_temperature,
                    permeate_pressure=conditions.permeate_pressure,
                    first_component_permeance=permeances[step][0],
                    second_component_permeance=permeances[step][1],
                )
            )

            permeate_composition.append(
                Composition(
                    p=partial_fluxes[step][0] / (sum(partial_fluxes[step])),
                    type=CompositionType.weight,
                )
            )

            d_mass_1 = partial_fluxes[step][0] * conditions.membrane_area * delta_hours
            d_mass_2 = partial_fluxes[step][1] * conditions.membrane_area * delta_hours

            feed_evaporation_heat.append(
                evaporation_heat_1 * d_mass_1 + evaporation_heat_2 * d_mass_2
            )
            if conditions.permeate_temperature is None:
                permeate_condensation_heat.append(None)
            else:
                permeate_condensation_heat.append(
                    condensation_heat_1 * d_mass_1
                    + condensation_heat_2 * d_mass_2
                    + (cooling_heat_1 * d_mass_1 + cooling_heat_2 * d_mass_2)
                    * (
                        conditions.initial_feed_temperature
                        - conditions.permeate_temperature
                    )
                )

            feed_mass.append(feed_mass[step] - d_mass_1 - d_mass_2)

            feed_composition.append(
                Composition(
                    p=(feed_composition[step].p * feed_mass[step] - d_mass_1)
                    / feed_mass[step + 1],
                    type=CompositionType.weight,
                )
            )
            permeances.append(
                (
                    Permeance(
                        value=pervaporation_function_first(
                            x=feed_composition[step].first,
                            t=conditions.initial_feed_temperature,
                        )
                        * facilitation_rate_first
                    ),
                    Permeance(
                        value=pervaporation_function_second(
                            x=feed_composition[step].first,
                            t=conditions.initial_feed_temperature,
                        )
                        * facilitation_rate_second
                    ),
                ),
            )

        feed_composition.pop(-1)
        feed_mass.pop(-1)
        permeances.pop(-1)

        return ProcessModel(
            mixture=self.mixture,
            membrane_name=self.membrane.name,
            feed_temperature=[conditions.initial_feed_temperature] * number_of_steps,
            feed_compositions=feed_composition,
            permeate_composition=permeate_composition,
            permeate_temperature=[conditions.permeate_temperature] * number_of_steps,
            permeate_pressure=[conditions.permeate_pressure] * number_of_steps,
            feed_mass=feed_mass,
            partial_fluxes=partial_fluxes,
            permeances=permeances,
            time=time,
            feed_evaporation_heat=feed_evaporation_heat,
            permeate_condensation_heat=permeate_condensation_heat,
            initial_conditions=conditions,
            permeance_fits=(
                pervaporation_function_first,
                pervaporation_function_second,
            ),
            comments=(
                f"{self.membrane.name} {self.mixture.first_component.name} / {self.mixture.second_component.name}"
                f"Non-ideal Process Model" + datetime.now().strftime("%m/%d/%Y, %H:%M")
            ),
        )

    def non_ideal_non_isothermal_process(
        self,
        conditions: Conditions,
        diffusion_curve_set: DiffusionCurveSet,
        number_of_steps: int,
        delta_hours: float,
        precision: typing.Optional[float] = 5e-5,
        initial_permeances: typing.Optional[typing.Tuple[Permeance, Permeance]] = None,
        n_first: typing.Optional[int] = None,
        m_first: typing.Optional[int] = None,
        n_second: typing.Optional[int] = None,
        m_second: typing.Optional[int] = None,
        include_zero: bool = False,
    ) -> ProcessModel:
        """
        The function models Non-Ideal Non-Isothermal Process
        Based on a set of Diffusion curves measured at different temperatures;
        The modelling could be also performed based on a single diffusion curve:
        In that case the apparent activation energy of transport is considered constant
        and is calculated for each test_components based on the IdealExperiments default_membranes provided for the
        Membrane.
        :param conditions: Initial Conditions of the Process, a Temperature program may be added if necessary
        :param diffusion_curve_set: A set of Diffusion curves picked for the Modelling from the Membrane
        :param number_of_steps: Number of time steps for modelling
        :param delta_hours: Size of each step in hours
        :param precision: Precision in obtained permeate composition, by default is 5e-5
        :param initial_permeances: Initial Permeances, should be stated if the Membranes swelling history is significant
        :param n_first: Optional parameter,
        indicates the order of the polynomial of the composition part of the Permeance function for
         the first component
        :param m_first: Optional parameter,
        indicates the order of the polynomial of the temperature part of the Permeance function for
         the first component
        :param n_second: Optional parameter,
        indicates the order of the polynomial of the composition part of the Permeance function for
         the second component
        :param m_second: Optional parameter,
        indicates the order of the polynomial of the temperature part of the Permeance function for
         the second component
        :param include_zero: if True, points:
         first_component_fraction = 0 first_component_permeance=0 for the first test_components
         first_component_fraction = 1 second_component_permeance=0 for the second test_components
         for each temperature are added to the measurements in order to improve obtained fits
        :return: ProcessModel object
        """
        for curve in diffusion_curve_set.diffusion_curves:
            for c in curve.feed_compositions:
                if c.type == CompositionType.molar:
                    c.to_weight(curve.mixture)

        time: typing.List[float] = [
            delta_hours * step for step in range(number_of_steps)
        ]

        feed_temperature: typing.List[float] = [conditions.initial_feed_temperature]

        partial_fluxes: typing.List[typing.Tuple[float, float]] = []

        permeate_composition: typing.List[Composition] = []

        feed_composition: typing.List[Composition] = [
            conditions.initial_feed_composition.to_weight(self.mixture)
        ]

        feed_evaporation_heat: typing.List[float] = []

        permeate_condensation_heat: typing.List[typing.Optional[float]] = []

        feed_mass: typing.List[float] = [conditions.initial_feed_amount]

        measurements_first = Measurements.from_diffusion_curves_first(
            diffusion_curve_set
        )
        measurements_second = Measurements.from_diffusion_curves_second(
            diffusion_curve_set
        )

        if len(diffusion_curve_set.diffusion_curves) == 1:
            pervaporation_function_temperature = diffusion_curve_set.diffusion_curves[
                0
            ].feed_temperature
            activation_energy_first = self.membrane.calculate_activation_energy(
                self.mixture.first_component
            )
            activation_energy_second = self.membrane.calculate_activation_energy(
                self.mixture.second_component
            )

            pervaporation_function_first = find_best_fit(
                data=measurements_first,
                n=n_first,
                m=0,
                include_zero=False,
                component_index=0,
            )
            pervaporation_function_second = find_best_fit(
                data=measurements_second,
                n=n_second,
                m=0,
                include_zero=False,
                component_index=1,
            )

            pervaporation_function_first = pervaporation_function_first * numpy.exp(
                -pervaporation_function_first.b[0] / pervaporation_function_temperature
                + activation_energy_first / (R * pervaporation_function_temperature)
            )

            pervaporation_function_second = pervaporation_function_second * numpy.exp(
                -pervaporation_function_second.b[0] / pervaporation_function_temperature
                + activation_energy_second / (R * pervaporation_function_temperature)
            )

            pervaporation_function_first.b[0] = activation_energy_first / R
            pervaporation_function_second.b[0] = activation_energy_second / R

        else:
            pervaporation_function_first = find_best_fit(
                data=measurements_first,
                n=n_first,
                m=m_first,
                include_zero=include_zero,
                component_index=0,
            )
            pervaporation_function_second = find_best_fit(
                data=measurements_second,
                n=n_second,
                m=m_second,
                include_zero=include_zero,
                component_index=1,
            )

        if initial_permeances is None:
            first_component_permeance = Permeance(
                value=pervaporation_function_first(
                    feed_composition[0].p, conditions.initial_feed_temperature
                )
            )

            second_component_permeance = Permeance(
                value=pervaporation_function_second(
                    feed_composition[0].p, conditions.initial_feed_temperature
                )
            )

        else:
            first_component_permeance = initial_permeances[0].convert(
                to_units=Units.kg_m2_h_kPa, component=self.mixture.first_component
            )
            second_component_permeance = initial_permeances[1].convert(
                to_units=Units.kg_m2_h_kPa, component=self.mixture.second_component
            )

        permeances: typing.List[typing.Tuple[Permeance, Permeance]] = [
            (first_component_permeance, second_component_permeance)
        ]

        facilitation_rate_first = (
            first_component_permeance.value
            / pervaporation_function_first(
                x=conditions.initial_feed_composition.first,
                t=conditions.initial_feed_temperature,
            )
        )
        facilitation_rate_second = (
            second_component_permeance.value
            / pervaporation_function_second(
                x=conditions.initial_feed_composition.first,
                t=conditions.initial_feed_temperature,
            )
        )

        for step in range(len(time)):

            evaporation_heat_1 = (
                self.mixture.first_component.get_vaporisation_heat(
                    feed_temperature[step]
                )
                / self.mixture.first_component.molecular_weight
                * 1000
            )
            evaporation_heat_2 = (
                self.mixture.second_component.get_vaporisation_heat(
                    feed_temperature[step]
                )
                / self.mixture.second_component.molecular_weight
                * 1000
            )

            heat_capacity_1 = (
                self.mixture.first_component.get_specific_heat(feed_temperature[step])
                / self.mixture.first_component.molecular_weight
            )
            heat_capacity_2 = (
                self.mixture.second_component.get_specific_heat(feed_temperature[step])
                / self.mixture.second_component.molecular_weight
            )
            feed_heat_capacity = (
                feed_composition[step].first * heat_capacity_1
                + feed_composition[step].second * heat_capacity_2
            )

            partial_fluxes.append(
                self.calculate_partial_fluxes(
                    feed_temperature=feed_temperature[step],
                    composition=feed_composition[step],
                    precision=precision,
                    permeate_temperature=conditions.permeate_temperature,
                    permeate_pressure=conditions.permeate_pressure,
                    first_component_permeance=permeances[step][0],
                    second_component_permeance=permeances[step][1],
                )
            )

            permeate_composition.append(
                Composition(
                    p=partial_fluxes[step][0] / (sum(partial_fluxes[step])),
                    type=CompositionType.weight,
                )
            )

            d_mass_1 = partial_fluxes[step][0] * conditions.membrane_area * delta_hours
            d_mass_2 = partial_fluxes[step][1] * conditions.membrane_area * delta_hours

            if conditions.permeate_temperature is None:
                permeate_condensation_heat.append(None)
            else:
                condensation_heat_1 = (
                    self.mixture.first_component.get_vaporisation_heat(
                        conditions.permeate_temperature
                    )
                    / self.mixture.first_component.molecular_weight
                    * 1000
                )
                condensation_heat_2 = (
                    self.mixture.second_component.get_vaporisation_heat(
                        conditions.permeate_temperature
                    )
                    / self.mixture.second_component.molecular_weight
                    * 1000
                )

                specific_heat_1 = self.mixture.first_component.get_cooling_heat(
                    feed_temperature[step], conditions.permeate_temperature
                )
                specific_heat_2 = self.mixture.second_component.get_cooling_heat(
                    feed_temperature[step], conditions.permeate_temperature
                )

                permeate_condensation_heat.append(
                    condensation_heat_1 * d_mass_1
                    + condensation_heat_2 * d_mass_2
                    + (specific_heat_1 * d_mass_1 + specific_heat_2 * d_mass_2)
                    * (feed_temperature[step] - conditions.permeate_temperature)
                )

            feed_evaporation_heat.append(
                evaporation_heat_1 * d_mass_1 + evaporation_heat_2 * d_mass_2
            )

            feed_mass.append(feed_mass[step] - d_mass_1 - d_mass_2)

            feed_composition.append(
                Composition(
                    p=(feed_composition[step].p * feed_mass[step] - d_mass_1)
                    / feed_mass[step + 1],
                    type=CompositionType.weight,
                )
            )

            if conditions.temperature_program is None:
                feed_temperature.append(
                    feed_temperature[step]
                    - (
                        feed_evaporation_heat[step]
                        / (feed_heat_capacity * feed_mass[step])
                    )
                )
            else:
                feed_temperature.append(
                    conditions.temperature_program.program(time[step] + delta_hours)
                )

            permeances.append(
                (
                    Permeance(
                        value=pervaporation_function_first(
                            x=feed_composition[step + 1].first,
                            t=feed_temperature[step + 1],
                        )
                        * facilitation_rate_first
                    ),
                    Permeance(
                        value=pervaporation_function_second(
                            x=feed_composition[step + 1].first,
                            t=feed_temperature[step + 1],
                        )
                        * facilitation_rate_second
                    ),
                ),
            )

        feed_composition.pop(-1)
        feed_mass.pop(-1)
        permeances.pop(-1)
        feed_temperature.pop(-1)

        return ProcessModel(
            mixture=self.mixture,
            membrane_name=self.membrane.name,
            feed_temperature=feed_temperature,
            permeate_temperature=[conditions.permeate_temperature] * number_of_steps,
            permeate_pressure=[conditions.permeate_pressure] * number_of_steps,
            feed_compositions=feed_composition,
            permeate_composition=permeate_composition,
            feed_mass=feed_mass,
            partial_fluxes=partial_fluxes,
            permeances=permeances,
            time=time,
            feed_evaporation_heat=feed_evaporation_heat,
            permeate_condensation_heat=permeate_condensation_heat,
            initial_conditions=conditions,
            permeance_fits=(
                pervaporation_function_first,
                pervaporation_function_second,
            ),
            comments=(
                f"{self.membrane.name} {self.mixture.first_component.name} / {self.mixture.second_component.name}"
                f"Non-ideal Process Model" + datetime.now().strftime("%m/%d/%Y, %H:%M")
            ),
            membrane_path=self.membrane.path,
        )
