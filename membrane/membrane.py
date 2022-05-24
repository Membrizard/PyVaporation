import typing

import attr
import numpy

from component import Component
from experiments import IdealExperiments
from utils import R


@attr.s(auto_attribs=True)
class Membrane:
    ideal_experiments: IdealExperiments
    name: str

    # non_ideal_experiments: typing.List[NonIdealExperiment]  # TODO: non ideal experiments

    # @property
    # def ideal_experiments_names(self) -> typing.List[str]:
    #     return [ie.component.name for ie in self.ideal_experiments.experiments]

    # Get all the penetrants the membrane was tested for
    # def get_known_penetrants(self) -> typing.List[Component]:
    #     return numpy.unique(
    #         [ideal_experiment.component for ideal_experiment in self.ideal_experiments]
    #     )

    def get_penetrant_data(self, component: Component) -> IdealExperiments:
        return IdealExperiments(
            experiments=list(
                filter(
                    lambda x: x.component.name == component.name,  # TODO: fix names
                    self.ideal_experiments.experiments,
                )
            )
        )

    # Calculate an apparent activation energy of permeation
    def calculate_activation_energy(
        self,
        component: Component,
        # rewrite: bool = True,  # TODO: make rewrite great again
    ) -> float:

        component_experiments = self.get_penetrant_data(component)

        if len(component_experiments) < 2:
            raise ValueError(
                "At Least Two points are required for the calculation of Apparent Activation Energy"
            )

        # Get all the temperature values corresponding to the permeances for a given penetrant convert to 1/T
        x = numpy.divide(
            1,
            [
                ideal_experiment.temperature
                for ideal_experiment in component_experiments.experiments
            ],
        )
        # Get all the permeance values for a given Penetrant convert to Ln(permeance)
        y = numpy.log(
            [
                ideal_experiment.permeance
                for ideal_experiment in component_experiments.experiments
            ]
        )
        # Converting Ln(permeance) to the matrix equation form
        a = numpy.vstack([x, numpy.ones(len(x))]).T
        # Calculation of Least-squares Linear Fit of Ln(Permeance) versus 1/T
        activation_energy, c = numpy.linalg.lstsq(a, y, rcond=-1)[0] * R

        return -activation_energy

    def get_permeance(
        self,
        temperature: float,
        component: Component,
    ) -> float:

        component_experiments = self.get_penetrant_data(component)

        temperature_list = [
            ideal_experiment.temperature
            for ideal_experiment in component_experiments.experiments
        ]

        index = min(
            range(len(temperature_list)),
            key=lambda i: abs(temperature_list[i] - temperature),
        )

        given_permeance = component_experiments.experiments[index].permeance

        if component_experiments.experiments[index].activation_energy is None:
            # TODO: add logs
            activation_energy = self.calculate_activation_energy(component)
            return given_permeance * numpy.exp(
                -activation_energy / R * (1 / temperature - 1 / temperature_list[index])
            )
        else:
            return given_permeance * numpy.exp(
                -component_experiments.experiments[index].activation_energy
                / R
                * (1 / temperature - 1 / temperature_list[index])
            )

    # calculating the selectivity
    def get_ideal_selectivity(
        self,
        temperature: float,
        first_component: Component,
        second_component: Component,
        type: typing.Optional[str] = "molar",
    ) -> float:
        if type == "weight":
            return self.get_permeance(
                temperature, first_component
            ) / self.get_permeance(temperature, second_component)
        elif type == "molar":
            return (
                self.get_permeance(temperature, first_component)
                / self.get_permeance(temperature, second_component)
                * second_component.molecular_weight
                / first_component.molecular_weight
            )

    def get_pure_component_flux(
        self,
        temperature: float,
        component: Component,
        permeate_temperature: typing.Optional[float] = None,
    ) -> float:
        if permeate_temperature is None:
            return self.get_permeance(
                temperature, component
            ) * component.get_vapor_pressure(temperature)
        else:
            return self.get_permeance(temperature, component) * (
                component.get_vapor_pressure(temperature)
                - component.get_vapor_pressure(permeate_temperature)
            )
