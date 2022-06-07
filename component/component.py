import typing
from pathlib import Path

import attr
import numpy
import yaml

from utils import (HeatCapacityConstants, R, VaporPressureConstants,
                   VPConstantsType)


@attr.s(auto_attribs=True)
class Component:
    """
    A class to create Components,
    Defined with required Component's constants
    """

    name: str
    molecular_weight: float = attr.ib(converter=lambda value: float(value))
    vapour_pressure_constants: VaporPressureConstants
    heat_capacity_constants: HeatCapacityConstants

    @classmethod
    def from_dict(cls, d: typing.Mapping) -> "Component":
        return Component(
            name=d["name"],
            molecular_weight=d["molecular_weight"],
            vapour_pressure_constants=VaporPressureConstants(
                **d["vapor_pressure_constants"]
            ),
            heat_capacity_constants=HeatCapacityConstants(
                **d["heat_capacity_constants"]
            ),
        )

    def get_vapor_pressure(self, temperature: float) -> float:
        """
        Calculation of saturated pressure in kPa at a given temperature in K using Antoine equation (by the basis of 10)
        type of the calculation using Antoine (antoine) log10(P)=a+b/(T+C)
        or Frost equation (Frost) ln(P) = a+b/T+c/T^2
        :param temperature: temperature in K
        :return: saturated pressure in kPa calculated with respect to constants and given temperature
        """
        if self.vapour_pressure_constants.type == VPConstantsType.antoine:
            return 10 ** (
                self.vapour_pressure_constants.a
                + self.vapour_pressure_constants.b
                / (temperature + self.vapour_pressure_constants.c)
            )
        elif self.vapour_pressure_constants.type == VPConstantsType.frost:
            return numpy.exp(
                self.vapour_pressure_constants.a
                + self.vapour_pressure_constants.b / temperature
                + self.vapour_pressure_constants.c / temperature**2
            )
        else:
            raise ValueError("Type of calculation not supported")

    def get_vaporisation_heat(self, temperature: float) -> float:
        """
        Calculation of Vaporisation heat in kJ/mol using Clapeyron-Clausius equation
        for Antoine equation: H=b*R*ln(10)*(T/(T+C))^2
        for Frost Equation: H=R*(-b-2*c/T)
        :param temperature: temperature in K
        :return: Vaporisation heat in kJ/mol
        """
        if self.vapour_pressure_constants.type == VPConstantsType.antoine:
            return (
                -(
                    (temperature / (temperature + self.vapour_pressure_constants.c))
                    ** 2
                    * R
                    * self.vapour_pressure_constants.b
                    * numpy.log(10)
                )
                / 1000
            )
        elif self.vapour_pressure_constants.type == VPConstantsType.frost:
            return (
                -R
                * (
                    self.vapour_pressure_constants.b
                    + 2 * self.vapour_pressure_constants.c / temperature
                )
                / 1000
            )
        else:
            raise ValueError("Type of calculation not supported")

    def get_specific_heat(self, temperature: float) -> float:
        """
        Calculation of Heat Capacity in J/(mol*K) using polynomial isobaric heat capacity fit
        :param temperature: temperature in K
        :return: Isobaric Heat capacity in J/(mol*K)
        """
        return (
            self.heat_capacity_constants.a
            + self.heat_capacity_constants.b * temperature
            + self.heat_capacity_constants.c * temperature**2
            + self.heat_capacity_constants.d * temperature**3
        )

    def get_cooling_heat(self, temperature_1, temperature_2):
        """
        Calculation of Specific Heat in J/mol using Integral (T2-T1) (CpdT)
        :param temperature_1 and temperature_2: temperature in K (temperature_2 < temperature_1)
        :return: Specific Heat in J/mol
        """
        return (
            self.heat_capacity_constants.a * (temperature_1 - temperature_2)
            + self.heat_capacity_constants.b
            * (temperature_1**2 - temperature_2**2)
            / 2
            + self.heat_capacity_constants.c
            * (temperature_1**3 - temperature_2**3)
            / 3
            + self.heat_capacity_constants.d
            * (temperature_1**4 - temperature_2**4)
            / 4
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
