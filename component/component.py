import typing
from pathlib import Path

import attr
import numpy
import yaml

from utils import AntoineConstants, HeatCapacityConstants, R


@attr.s(auto_attribs=True)
class Component:
    name: str
    molecular_weight: float = attr.ib(converter=lambda value: float(value))
    antoine_constants: AntoineConstants
    heat_capacity_constants: HeatCapacityConstants

    @classmethod
    def from_dict(cls, d: typing.Mapping) -> "Component":
        return Component(
            name=d["name"],
            molecular_weight=d["molecular_weight"],
            antoine_constants=AntoineConstants(**d["antoine_constants"]),
            heat_capacity_constants=HeatCapacityConstants(**d["heat_capacity_constants"])
        )

    def get_antoine_pressure(self, temperature: float) -> float:
        """
        Calculation of saturated pressure in kPa at a given temperature in K using Antoine equation (by the basis of 10)
        :param temperature: temperature in K
        :return: saturated pressure in kPa calculated with respect to constants and given temperature
        """
        return 10 ** (
            self.antoine_constants.a
            - self.antoine_constants.b / (temperature + self.antoine_constants.c)
        )

    def get_vaporisation_heat(self, temperature: float) -> float:
        """
        Calculation of Vaporisation heat in kJ/mol using Clapeyron-Clausius equation
        :param temperature: temperature in K
        :return: Vaporisation heat in kJ/mol
        """
        return (
            (temperature / (temperature + self.antoine_constants.c)) ** 2
            * R
            * self.antoine_constants.b
            * numpy.log(10)
        )

    def get_heat_capacity(self, temperature: float) -> float:
        """
        Calculation of Vaporisation heat in J/(mol*K) using polynomial isobaric heat capacity fit
        :param temperature: temperature in K
        :return: Isobaric Heat capacity in J/(mol*K)
        """
        return (
            self.heat_capacity_constants.a
            + self.heat_capacity_constants.b * temperature
            + self.heat_capacity_constants.c * temperature**2
            + self.heat_capacity_constants.d * temperature**3
        )


@attr.s(auto_attribs=True)
class AllComponents:
    components: typing.Mapping[str, Component]

    @classmethod
    def load(cls, path: typing.Union[str, Path]) -> "AllComponents":
        with open(path, "r") as handle:
            _components = yaml.load(handle, Loader=yaml.FullLoader)

        output = AllComponents(
            components={
                name: Component.from_dict(value) for name, value in _components.items()
            }
        )

        for name, component in output.components.items():
            setattr(output, name, component)

        return output
